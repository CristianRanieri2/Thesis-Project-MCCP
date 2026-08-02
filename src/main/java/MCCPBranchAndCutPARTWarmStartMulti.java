import com.google.ortools.Loader;
import com.google.ortools.sat.BoolVar;
import com.google.ortools.sat.CpModel;
import com.google.ortools.sat.CpSolver;
import com.google.ortools.sat.CpSolverStatus;
import com.google.ortools.sat.LinearExpr;
import com.google.ortools.sat.LinearExprBuilder;

import java.util.*;
import java.util.concurrent.*;

/**
 * Solutore MCCP PART_{s-t} che sfrutta CP-SAT di Google OR-Tools.
 *
 * Parallelizzazione:
 * 1. Warm-Start VNS Parallelo (Multi-Thread)
 * 2. Solutore CP-SAT Parallelo (Worker su tutti i core CPU per il B&B)
 */
public class MCCPBranchAndCutPARTWarmStartMulti {

    static {
        Loader.loadNativeLibraries();
    }

    private final int numNodes;
    private final List<MCCPSolver.Edge> edges;
    private final int numColors;
    private final double[] colorCost;
    private final int s;
    private final int t;

    private long nodesExplored = 0;
    private boolean provenOptimal = false;
    private long vnsWarmStartCost = -1;
    private long vnsWarmStartTimeMillisUsed = 0;

    public MCCPBranchAndCutPARTWarmStartMulti(int numNodes, List<MCCPSolver.Edge> edges, int numColors,
                                              double[] colorCost, int sourceNode, int targetNode) {
        if (colorCost.length != numColors) {
            throw new IllegalArgumentException("colorCost deve avere esattamente numColors elementi");
        }
        if (sourceNode < 0 || sourceNode >= numNodes || targetNode < 0 || targetNode >= numNodes) {
            throw new IllegalArgumentException("sourceNode/targetNode fuori dal range dei nodi");
        }
        if (sourceNode == targetNode) {
            throw new IllegalArgumentException("sourceNode e targetNode devono essere due nodi distinti");
        }
        this.numNodes = numNodes;
        this.edges = edges;
        this.numColors = numColors;
        this.colorCost = colorCost.clone();
        this.s = sourceNode;
        this.t = targetNode;
    }

    public long getNodesExplored() { return nodesExplored; }
    public boolean isProvenOptimal() { return provenOptimal; }
    public long getVnsWarmStartCost() { return vnsWarmStartCost; }
    public long getVnsWarmStartTimeMillisUsed() { return vnsWarmStartTimeMillisUsed; }

    private static long vnsWarmStartTimeMillis(int numNodes) {
        if (numNodes <= 50) return 1_000L;
        if (numNodes <= 100) return 2_000L;
        if (numNodes <= 200) return 3_000L;
        if (numNodes <= 400) return 8_000L;
        if (numNodes <= 500) return 20_000L;
        if (numNodes <= 600) return 40_000L;
        if (numNodes <= 700) return 80_000L;
        if (numNodes <= 800) return 160_000L;
        if (numNodes <= 1000) return 280_000L;
        return 360_000L;
    }

    public MCCPSolver.MCCPResult solveExactWithVNSWarmStartParallel(MCCPSolver instance,
                                                                    long maxRunningTimeMillisForBC,
                                                                    int numParallelVnsRuns) {
        int defaultThreads = Math.max(1, Runtime.getRuntime().availableProcessors());
        return solveExactWithVNSWarmStartParallel(instance, maxRunningTimeMillisForBC, numParallelVnsRuns, defaultThreads);
    }

    public MCCPSolver.MCCPResult solveExactWithVNSWarmStartParallel(MCCPSolver instance,
                                                                    long maxRunningTimeMillisForBC,
                                                                    int numParallelVnsRuns,
                                                                    int numSolverThreads) {
        if (numParallelVnsRuns <= 0) {
            throw new IllegalArgumentException("numParallelVnsRuns deve essere positivo");
        }

        long vnsTimeMillis = vnsWarmStartTimeMillis(numNodes);
        System.out.println("[Warm start parallelo] " + numParallelVnsRuns
                + " esecuzioni VNS-Probabilistic (budget = " + vnsTimeMillis + " ms)");

        long tVnsStart = System.currentTimeMillis();
        ExecutorService executor = Executors.newFixedThreadPool(numParallelVnsRuns);
        List<Future<MCCPSolver.MCCPResult>> futures = new ArrayList<>();
        try {
            for (int i = 0; i < numParallelVnsRuns; i++) {
                futures.add(executor.submit(() -> instance.solveProbabilistic(vnsTimeMillis)));
            }

            MCCPSolver.MCCPResult best = null;
            for (Future<MCCPSolver.MCCPResult> f : futures) {
                MCCPSolver.MCCPResult r = f.get();
                if (best == null || r.cutCost < best.cutCost) {
                    best = r;
                }
            }
            long tVnsEnd = System.currentTimeMillis();

            this.vnsWarmStartCost = (long) best.cutCost;
            this.vnsWarmStartTimeMillisUsed = tVnsEnd - tVnsStart;

            System.out.println("[Warm start parallelo] Miglior costo VNS = " + best.cutCost
                    + " [" + vnsWarmStartTimeMillisUsed + " ms]");

            return solveExactCPSAT(maxRunningTimeMillisForBC, best.cutColors, numSolverThreads);
        } catch (InterruptedException | ExecutionException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Errore durante il warm start parallelo", e);
        } finally {
            executor.shutdown();
        }
    }

    // ========================================================================
    // COSTRUZIONE E RISOLUZIONE DEL MODELLO PART CON CP-SAT (PARALLELO NATIVO)
    // ========================================================================

    private MCCPSolver.MCCPResult solveExactCPSAT(long maxRunningTimeMillis, Set<Integer> warmStartCutColors, int numSolverThreads) {
        long searchStartMillis = System.currentTimeMillis();

        CpModel model = new CpModel();

        // 1. Variabili Colore: z_c in {0, 1}
        BoolVar[] z = new BoolVar[numColors];
        for (int c = 0; c < numColors; c++) {
            z[c] = model.newBoolVar("z_" + c);
        }

        // 2. Variabili Partizione Nodi: w_v in {0, 1}
        BoolVar[] w = new BoolVar[numNodes];
        for (int v = 0; v < numNodes; v++) {
            w[v] = model.newBoolVar("w_" + v);
        }

        // 3. Vincoli di Sorgente (s) e Pozzo (t)
        model.addEquality(w[s], 1);
        model.addEquality(w[t], 0);

        // 4. Vincoli d'Arco:
        // sum(z_c) + w_v >= w_u
        // sum(z_c) + w_u >= w_v
        for (MCCPSolver.Edge e : edges) {
            LinearExprBuilder expr1 = LinearExpr.newBuilder();
            for (int c : e.colors) {
                expr1.add(z[c]);
            }
            expr1.add(w[e.v]);
            model.addGreaterOrEqual(expr1.build(), w[e.u]);

            LinearExprBuilder expr2 = LinearExpr.newBuilder();
            for (int c : e.colors) {
                expr2.add(z[c]);
            }
            expr2.add(w[e.u]);
            model.addGreaterOrEqual(expr2.build(), w[e.v]);
        }

        // 5. Funzione Obiettivo: min sum(cost_c * z_c)
        LinearExprBuilder obj = LinearExpr.newBuilder();
        for (int c = 0; c < numColors; c++) {
            obj.addTerm(z[c], (long) Math.round(colorCost[c]));
        }
        model.minimize(obj.build());

        // 6. WARM START / HINT
        if (warmStartCutColors != null) {
            double[] wHintValues = computeWHint(warmStartCutColors);
            for (int c = 0; c < numColors; c++) {
                model.addHint(z[c], warmStartCutColors.contains(c) ? 1 : 0);
            }
            for (int v = 0; v < numNodes; v++) {
                model.addHint(w[v], wHintValues[v] > 0.5 ? 1 : 0);
            }
        }

        // 7. CONFIGURAZIONE SOLUTORE CP-SAT MULTI-THREAD
        CpSolver solver = new CpSolver();

        // Imposta direttamente il numero di worker/thread di CP-SAT
        solver.getParameters().setNumWorkers(numSolverThreads);
        System.out.println("[CP-SAT Solver] Allocati " + numSolverThreads + " worker paralleli su CPU.");

        if (maxRunningTimeMillis > 0) {
            solver.getParameters().setMaxTimeInSeconds(maxRunningTimeMillis / 1000.0);
        }

        // Risoluzione
        CpSolverStatus status = solver.solve(model);

        long totalTimeMs = System.currentTimeMillis() - searchStartMillis;
        this.nodesExplored = solver.numBranches();
        this.provenOptimal = (status == CpSolverStatus.OPTIMAL);

        Set<Integer> cutColors = new HashSet<>();
        Set<Integer> keptColors = new HashSet<>();
        double bestCost;

        if (status == CpSolverStatus.OPTIMAL || status == CpSolverStatus.FEASIBLE) {
            bestCost = solver.objectiveValue();
            for (int c = 0; c < numColors; c++) {
                if (solver.booleanValue(z[c])) {
                    cutColors.add(c);
                } else {
                    keptColors.add(c);
                }
            }
        } else {
            bestCost = Double.POSITIVE_INFINITY;
        }

        return new MCCPSolver.MCCPResult(cutColors, keptColors, bestCost, totalTimeMs, totalTimeMs);
    }

    private double[] computeWHint(Set<Integer> cutColors) {
        List<List<Integer>> adjacency = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) adjacency.add(new ArrayList<>());
        for (MCCPSolver.Edge e : edges) {
            if (Collections.disjoint(e.colors, cutColors)) {
                adjacency.get(e.u).add(e.v);
                adjacency.get(e.v).add(e.u);
            }
        }

        boolean[] visited = new boolean[numNodes];
        Deque<Integer> queue = new ArrayDeque<>();
        visited[s] = true;
        queue.add(s);
        while (!queue.isEmpty()) {
            int u = queue.poll();
            for (int v : adjacency.get(u)) {
                if (!visited[v]) {
                    visited[v] = true;
                    queue.add(v);
                }
            }
        }

        double[] wHint = new double[numNodes];
        for (int v = 0; v < numNodes; v++) wHint[v] = visited[v] ? 1.0 : 0.0;
        return wHint;
    }
}