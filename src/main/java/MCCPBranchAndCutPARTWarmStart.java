import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;
import java.util.concurrent.*;

public class MCCPBranchAndCutPARTWarmStart {

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
    private int rootEdgeConstraintCount = 0;
    private long vnsWarmStartCost = -1;
    private long vnsWarmStartTimeMillisUsed = 0;

    public MCCPBranchAndCutPARTWarmStart(int numNodes, List<MCCPSolver.Edge> edges, int numColors,
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
    public int getRootEdgeConstraintCount() { return rootEdgeConstraintCount; }
    public long getVnsWarmStartCost() { return vnsWarmStartCost; }
    public long getVnsWarmStartTimeMillisUsed() { return vnsWarmStartTimeMillisUsed; }

    private static long vnsWarmStartTimeMillis(int numNodes) {
        if (numNodes <= 50) return 1_000L;
        if (numNodes <= 100) return 20_000L;
        if (numNodes <= 200) return 30_000L;
        if (numNodes <= 400) return 80_000L;
        if (numNodes <= 500) return 200_000L;
        if (numNodes <= 600) return 400_000L;
        if (numNodes <= 700) return 800_000L;
        if (numNodes <= 800) return 1600_000L;
        if (numNodes <= 1000) return 2_800_000L;
        return 3_600_000L;
    }

    public MCCPSolver.MCCPResult solveExactWithVNSWarmStart(MCCPSolver instance, long maxRunningTimeMillisForBC) {
        long vnsTimeMillis = vnsWarmStartTimeMillis(numNodes);

        System.out.println("[Warm start] VNS-Probabilistic: budget = " + vnsTimeMillis
                + " ms (|V|=" + numNodes + ", secondo la Tabella 1 del paper)");
        long tVnsStart = System.currentTimeMillis();
        MCCPSolver.MCCPResult vnsResult = instance.solveProbabilistic(vnsTimeMillis);
        long tVnsEnd = System.currentTimeMillis();

        this.vnsWarmStartCost = (long) vnsResult.cutCost;
        this.vnsWarmStartTimeMillisUsed = tVnsEnd - tVnsStart;

        System.out.println("[Warm start] VNS-Probabilistic -> costo = " + vnsResult.cutCost
                + "  [" + vnsWarmStartTimeMillisUsed + " ms effettivi]");

        return solveExactInternal(maxRunningTimeMillisForBC, vnsResult.cutColors);
    }

    public MCCPSolver.MCCPResult solveExactWithVNSWarmStartParallel(MCCPSolver instance, long maxRunningTimeMillisForBC,
                                                                    int numParallelVnsRuns) {
        if (numParallelVnsRuns <= 0) {
            throw new IllegalArgumentException("numParallelVnsRuns deve essere positivo");
        }

        long vnsTimeMillis = vnsWarmStartTimeMillis(numNodes);
        System.out.println("[Warm start parallelo] " + numParallelVnsRuns
                + " esecuzioni concorrenti di VNS-Probabilistic, budget ciascuna = " + vnsTimeMillis
                + " ms (|V|=" + numNodes + ")");

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

            System.out.println("[Warm start parallelo] migliore fra " + numParallelVnsRuns
                    + " run -> costo = " + best.cutCost + "  [" + vnsWarmStartTimeMillisUsed + " ms di parete]");

            return solveExactInternal(maxRunningTimeMillisForBC, best.cutColors);
        } catch (InterruptedException | ExecutionException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Errore durante il warm start parallelo", e);
        } finally {
            executor.shutdown();
        }
    }

    public MCCPSolver.MCCPResult solveExactWithVNSWarmStart(MCCPSolver instance, long maxRunningTimeMillisForBC,
                                                            long vnsWarmStartTimeMillisOverride) {
        System.out.println("[Warm start] VNS-Probabilistic: budget = " + vnsWarmStartTimeMillisOverride
                + " ms (esplicito, |V|=" + numNodes + ")");
        long tVnsStart = System.currentTimeMillis();
        MCCPSolver.MCCPResult vnsResult = instance.solveProbabilistic(vnsWarmStartTimeMillisOverride);
        long tVnsEnd = System.currentTimeMillis();

        this.vnsWarmStartCost = (long) vnsResult.cutCost;
        this.vnsWarmStartTimeMillisUsed = tVnsEnd - tVnsStart;

        System.out.println("[Warm start] VNS-Probabilistic -> costo = " + vnsResult.cutCost
                + "  [" + vnsWarmStartTimeMillisUsed + " ms effettivi]");

        return solveExactInternal(maxRunningTimeMillisForBC, vnsResult.cutColors);
    }

    public MCCPSolver.MCCPResult solveExactWithVNSWarmStart(MCCPSolver2 instance2, long vnsMaxTimeMillis,
                                                            int vnsMaxStagnantIterations,
                                                            long maxRunningTimeMillisForBC) {
        System.out.println("[Warm start] VNS-Probabilistic (MCCPSolver2): budget tempo = " + vnsMaxTimeMillis
                + " ms, stagnazione max = " + vnsMaxStagnantIterations + " iterazioni consecutive");
        long tVnsStart = System.currentTimeMillis();
        MCCPSolver.MCCPResult vnsResult = instance2.solveProbabilistic(vnsMaxTimeMillis, vnsMaxStagnantIterations);
        long tVnsEnd = System.currentTimeMillis();

        this.vnsWarmStartCost = (long) vnsResult.cutCost;
        this.vnsWarmStartTimeMillisUsed = tVnsEnd - tVnsStart;

        System.out.println("[Warm start] VNS-Probabilistic (MCCPSolver2) -> costo = " + vnsResult.cutCost
                + "  [" + vnsWarmStartTimeMillisUsed + " ms effettivi, "
                + instance2.getLastRunTotalIterations() + " iterazioni, fermato per stagnazione="
                + instance2.isLastRunStoppedByStagnation() + "]");

        return solveExactInternal(maxRunningTimeMillisForBC, vnsResult.cutColors);
    }

    // ========================================================================
    // ALGORITMO APPROSSIMATO (Suddivisione budget con garanzia temporale)
    // ========================================================================

    public MCCPSolver.MCCPResult solveApproximate(MCCPSolver instance, long totalTimeMillis) {
        return solveApproximate(instance, totalTimeMillis, 0.5);
    }

    public MCCPSolver.MCCPResult solveApproximate(MCCPSolver instance, long totalTimeMillis, double vnsTimeFraction) {
        if (totalTimeMillis <= 0) {
            throw new IllegalArgumentException("totalTimeMillis deve essere positivo per l'algoritmo approssimato");
        }
        if (vnsTimeFraction <= 0.0 || vnsTimeFraction >= 1.0) {
            throw new IllegalArgumentException("vnsTimeFraction deve essere strettamente compreso tra 0 e 1");
        }

        long vnsTimeMillis = Math.round(totalTimeMillis * vnsTimeFraction);
        long partTimeMillis = totalTimeMillis - vnsTimeMillis;

        System.out.println("[Approssimato] budget totale = " + totalTimeMillis + " ms, suddiviso: warm start = "
                + vnsTimeMillis + " ms (" + Math.round(vnsTimeFraction * 100) + "%), PART = "
                + partTimeMillis + " ms (" + Math.round((1 - vnsTimeFraction) * 100) + "%)");

        long tVnsStart = System.currentTimeMillis();
        // Esegue il VNS con controlli temporali capillari attivi
        MCCPSolver.MCCPResult vnsResult = instance.solveProbabilistic(vnsTimeMillis);
        long tVnsEnd = System.currentTimeMillis();

        this.vnsWarmStartCost = (long) vnsResult.cutCost;
        this.vnsWarmStartTimeMillisUsed = tVnsEnd - tVnsStart;

        System.out.println("[Approssimato] Warm start -> costo = " + vnsResult.cutCost
                + "  [" + vnsWarmStartTimeMillisUsed + " ms effettivi]");

        MCCPSolver.MCCPResult finalResult = solveExactInternal(partTimeMillis, vnsResult.cutColors);

        if (provenOptimal) {
            System.out.println("[Approssimato] La fase PART ha comunque dimostrato l'ottimo entro il tempo concesso.");
        } else {
            System.out.println("[Approssimato] ATTENZIONE: la fase PART NON ha dimostrato l'ottimo entro i "
                    + partTimeMillis + " ms concessi -> il risultato restituito e' un'APPROSSIMAZIONE "
                    + "(upper bound), non l'ottimo certificato.");
        }

        return finalResult;
    }

    private MCCPSolver.MCCPResult solveExactInternal(long maxRunningTimeMillis, Set<Integer> warmStartCutColors) {
        long searchStartMillis = System.currentTimeMillis();

        MPSolver solver = MPSolver.createSolver("CBC");
        if (solver == null) {
            solver = MPSolver.createSolver("SCIP");
        }
        if (solver == null) {
            throw new RuntimeException("Impossibile caricare i solutori MILP di OR-Tools (CBC/SCIP).");
        }

        int numThreads = Runtime.getRuntime().availableProcessors();
        solver.setNumThreads(numThreads);

        if (maxRunningTimeMillis > 0) {
            solver.setTimeLimit(maxRunningTimeMillis);
        }

        MPVariable[] z = new MPVariable[numColors];
        for (int c = 0; c < numColors; c++) {
            z[c] = solver.makeBoolVar("z_" + c);
        }

        MPVariable[] w = new MPVariable[numNodes];
        for (int v = 0; v < numNodes; v++) {
            w[v] = solver.makeNumVar(0.0, 1.0, "w_" + v);
        }

        MPConstraint sConstraint = solver.makeConstraint(1.0, 1.0, "w_s_bound");
        sConstraint.setCoefficient(w[s], 1.0);

        MPConstraint tConstraint = solver.makeConstraint(0.0, 0.0, "w_t_bound");
        tConstraint.setCoefficient(w[t], 1.0);

        int totalConstraints = 2;

        for (MCCPSolver.Edge e : edges) {
            MPConstraint c1 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
            for (int c : e.colors) {
                c1.setCoefficient(z[c], 1.0);
            }
            c1.setCoefficient(w[e.u], -1.0);
            c1.setCoefficient(w[e.v], 1.0);
            totalConstraints++;

            MPConstraint c2 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
            for (int c : e.colors) {
                c2.setCoefficient(z[c], 1.0);
            }
            c2.setCoefficient(w[e.u], 1.0);
            c2.setCoefficient(w[e.v], -1.0);
            totalConstraints++;
        }

        this.rootEdgeConstraintCount = totalConstraints;

        MPObjective objective = solver.objective();
        for (int c = 0; c < numColors; c++) {
            objective.setCoefficient(z[c], colorCost[c]);
        }
        objective.setMinimization();

        if (warmStartCutColors != null) {
            double[] wHintValues = computeWHint(warmStartCutColors);

            MPVariable[] hintVars = new MPVariable[numColors + numNodes];
            double[] hintVals = new double[numColors + numNodes];

            for (int c = 0; c < numColors; c++) {
                hintVars[c] = z[c];
                hintVals[c] = warmStartCutColors.contains(c) ? 1.0 : 0.0;
            }
            for (int v = 0; v < numNodes; v++) {
                hintVars[numColors + v] = w[v];
                hintVals[numColors + v] = wHintValues[v];
            }

            solver.setHint(hintVars, hintVals);
        }

        MPSolver.ResultStatus status = solver.solve();

        long totalTimeMs = System.currentTimeMillis() - searchStartMillis;
        this.nodesExplored = solver.nodes();
        this.provenOptimal = (status == MPSolver.ResultStatus.OPTIMAL);

        Set<Integer> cutColors = new HashSet<>();
        Set<Integer> keptColors = new HashSet<>();
        double bestCost;

        if (status == MPSolver.ResultStatus.OPTIMAL || status == MPSolver.ResultStatus.FEASIBLE) {
            bestCost = objective.value();
            for (int c = 0; c < numColors; c++) {
                if (z[c].solutionValue() > 0.5) {
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