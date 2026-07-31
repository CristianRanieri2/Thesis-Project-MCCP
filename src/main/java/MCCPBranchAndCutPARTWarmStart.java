import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;

/**
 * Variante di {@link MCCPBranchAndCutPART} che usa come punto di partenza
 * (warm start / hint iniziale) la soluzione trovata da VNS-Probabilistic
 * (vedi {@link MCCPSolver#solveProbabilistic(long)}), invece di lasciare
 * che OR-Tools parta da zero.
 *
 * Risolve esattamente lo stesso modello PART_{s-t} di {@link MCCPBranchAndCutPART}
 * (stessa funzione obiettivo, stessi vincoli) — cambia solo la strategia di
 * ricerca: prima si genera una buona soluzione ammissibile con la metaeuristica
 * VNS, poi la si passa a CBC/SCIP come incumbent di partenza (MPSolver.setHint),
 * cosi' il branch-and-bound puo' fare pruning aggressivo fin dal primo nodo
 * invece di dover scoprire da solo una prima soluzione ammissibile.
 *
 * ============================================================================
 * BUDGET DI TEMPO PER IL WARM START
 * ============================================================================
 * Il tempo concesso a VNS-Probabilistic per il warm start e' determinato
 * automaticamente dal numero di nodi del grafo, secondo la Tabella 1 del
 * paper originale (Bordini & Protti, 2017, Sezione 2.1):
 *
 *   |V| <=   50  ->    1 s
 *   |V| <=  100  ->   20 s
 *   |V| <=  200  ->   30 s
 *   |V| <=  400  ->   80 s
 *   |V| <=  500  ->  200 s
 *   |V| <= 1000  -> 2800 s
 *   |V| >  1000  -> 3600 s
 *
 * Questo tempo e' AGGIUNTIVO rispetto al budget passato a
 * {@link #solveExactWithVNSWarmStart}, che viene concesso a OR-Tools DOPO
 * che il warm start e' gia' stato calcolato.
 *
 * ============================================================================
 * COME VIENE TRADOTTA LA SOLUZIONE DI VNS IN UN HINT
 * ============================================================================
 * VNS restituisce un insieme di colori tagliati (cutColors). Da questo si
 * ricava direttamente l'hint per le variabili z (1.0 se il colore e' nel
 * taglio, 0.0 altrimenti). L'hint per le variabili di partizione w si ricava
 * invece con una BFS sul sottografo che sopravvive alla rimozione di quei
 * colori: w_v = 1 se v e' ancora raggiungibile da s (coerente con w_s=1),
 * w_v = 0 altrimenti (coerente con w_t=0, dato che il taglio di VNS e', per
 * costruzione, sempre valido: t non e' mai raggiungibile da s).
 */
public class MCCPBranchAndCutPARTWarmStart {

    static {
        // Caricamento delle librerie native C++ di Google OR-Tools
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
    /** Costo della soluzione di warm start trovata da VNS-Probabilistic (-1 se non ancora eseguita). */
    public long getVnsWarmStartCost() { return vnsWarmStartCost; }
    /** Tempo (ms) effettivamente impiegato dal VNS-Probabilistic per il warm start. */
    public long getVnsWarmStartTimeMillisUsed() { return vnsWarmStartTimeMillisUsed; }

    // ========================================================================
    // Tabella 1 (Bordini & Protti, 2017): tempo massimo per VNS in funzione
    // del numero di nodi, riusata qui come budget per il warm start.
    // ========================================================================

    private static long vnsWarmStartTimeMillis(int numNodes) {
        if (numNodes <= 50) return 1_000L;
        if (numNodes <= 100) return 20_000L;
        if (numNodes <= 200) return 30_000L;
        if (numNodes <= 400) return 80_000L;
        if (numNodes <= 500) return 200_000L;
        if (numNodes <= 1000) return 2_800_000L;
        return 3_600_000L;
    }

    // ========================================================================
    // RISOLUZIONE TRAMITE OR-TOOLS, con warm start da VNS-Probabilistic
    // ========================================================================

    /**
     * Esegue prima VNS-Probabilistic su {@code instance} (con il budget di
     * tempo determinato automaticamente dal numero di nodi, Tabella 1 del
     * paper), poi risolve il modello PART_{s-t} con OR-Tools passando quella
     * soluzione come hint iniziale.
     *
     * @param instance                   l'istanza del problema (serve per
     *                                    poter chiamare solveProbabilistic);
     *                                    deve rappresentare LO STESSO
     *                                    grafo/costi/s/t con cui e' stato
     *                                    costruito questo solver
     * @param maxRunningTimeMillisForBC   tempo massimo (ms) concesso a
     *                                    OR-Tools DOPO il warm start (il
     *                                    tempo speso dal VNS e' aggiuntivo,
     *                                    non sottratto da questo budget)
     */
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

    /** Come sopra, ma con un budget esplicito per il warm start invece di quello della Tabella 1. */
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

    // ========================================================================
    // Costruzione del modello PART_{s-t} e risoluzione con hint
    // ========================================================================

    private MCCPSolver.MCCPResult solveExactInternal(long maxRunningTimeMillis, Set<Integer> warmStartCutColors) {
        long searchStartMillis = System.currentTimeMillis();

        MPSolver solver = MPSolver.createSolver("CBC");
        if (solver == null) {
            solver = MPSolver.createSolver("SCIP");
        }
        if (solver == null) {
            throw new RuntimeException("Impossibile caricare i solutori MILP di OR-Tools (CBC/SCIP).");
        }

        if (maxRunningTimeMillis > 0) {
            solver.setTimeLimit(maxRunningTimeMillis);
        }

        // 1. Variabili Colore: z_c in {0, 1}
        MPVariable[] z = new MPVariable[numColors];
        for (int c = 0; c < numColors; c++) {
            z[c] = solver.makeBoolVar("z_" + c);
        }

        // 2. Variabili Partizione Nodi: w_v in [0, 1]
        MPVariable[] w = new MPVariable[numNodes];
        for (int v = 0; v < numNodes; v++) {
            w[v] = solver.makeNumVar(0.0, 1.0, "w_" + v);
        }

        // 3. Vincoli di Confine per i nodi sorgente (s) e pozzo (t)
        MPConstraint sConstraint = solver.makeConstraint(1.0, 1.0, "w_s_bound");
        sConstraint.setCoefficient(w[s], 1.0);

        MPConstraint tConstraint = solver.makeConstraint(0.0, 0.0, "w_t_bound");
        tConstraint.setCoefficient(w[t], 1.0);

        int totalConstraints = 2;

        // 4. Vincoli d'Arco del Modello PART_{s-t}
        // Per ogni arco e = (u,v):
        //   (1) sum_{c in e.colors} z_c - w_u + w_v >= 0
        //   (2) sum_{c in e.colors} z_c + w_u - w_v >= 0
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

        // 5. Funzione Obiettivo
        MPObjective objective = solver.objective();
        for (int c = 0; c < numColors; c++) {
            objective.setCoefficient(z[c], colorCost[c]);
        }
        objective.setMinimization();

        // 6. WARM START: traduce la soluzione trovata da VNS in un hint per
        // z (1.0/0.0 sui colori) e per w (calcolato con una BFS sul
        // sottografo che sopravvive al taglio di VNS).
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

        // Esecuzione della risoluzione
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

    /**
     * Calcola il valore di hint per ciascuna variabile w_v a partire da un
     * insieme di colori tagliati: w_v = 1 se v e' raggiungibile da s nel
     * sottografo che sopravvive alla rimozione di quei colori, 0 altrimenti
     * (BFS pura, nessuna dipendenza da OR-Tools).
     */
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

    // ========================================================================
    // Main di Test e Benchmark: B&C-PART base vs B&C-PART con warm start VNS
    // ========================================================================

    public static void main(String[] args) {
        System.out.println("=== TEST: B&C-PART (base) vs B&C-PART con warm start VNS-Probabilistic ===");

        for (int trial = 0; trial < 3; trial++) {
            int numNodes = 200 + trial * 100;
            int numColors = 100 + 50 * trial;

            MCCPSolver instance = MCCPSolver.generateRandomInstance(numNodes, numColors, MCCPSolver.Density.MEDIUM, 0.10);

            // B&C-PART base (senza warm start), dal file MCCPBranchAndCutPART.java
            // Nessun limite di tempo: passando 0 non scatta mai setTimeLimit(...),
            // quindi OR-Tools gira finche' non dimostra l'ottimo.
            MCCPBranchAndCutPART partBase = new MCCPBranchAndCutPART(
                    instance.getNumNodes(), instance.getEdges(), instance.getNumColors(),
                    instance.getColorCost(), instance.getSourceNode(), instance.getTargetNode());
            long t0 = System.currentTimeMillis();
            MCCPSolver.MCCPResult baseResult = partBase.solveExact(0);
            long t1 = System.currentTimeMillis();

            // B&C-PART con warm start VNS-Probabilistic. Il tempo del warm start
            // resta quello della Tabella 1 (calcolato internamente in base a |V|);
            // la fase OR-Tools successiva gira anch'essa senza limite di tempo.
            MCCPBranchAndCutPARTWarmStart partWarm = new MCCPBranchAndCutPARTWarmStart(
                    instance.getNumNodes(), instance.getEdges(), instance.getNumColors(),
                    instance.getColorCost(), instance.getSourceNode(), instance.getTargetNode());
            long t2 = System.currentTimeMillis();
            MCCPSolver.MCCPResult warmResult = partWarm.solveExactWithVNSWarmStart(instance, 0);
            long t3 = System.currentTimeMillis();

            System.out.println("\n[Trial " + trial + "] " + numNodes + " nodi, " + numColors
                    + " colori, s=" + instance.getSourceNode() + ", t=" + instance.getTargetNode());

            System.out.println("B&C-PART (base)          -> costo=" + baseResult.cutCost
                    + " | nodi B&B=" + partBase.getNodesExplored()
                    + " | ottimo certificato=" + partBase.isProvenOptimal()
                    + " | tempo totale=" + (t1 - t0) + " ms");
            instance.verifyCutWithBFS(baseResult.cutColors, false);

            System.out.println("B&C-PART (warm start VNS) -> costo=" + warmResult.cutCost
                    + " | warm start VNS costo=" + partWarm.getVnsWarmStartCost()
                    + " (" + partWarm.getVnsWarmStartTimeMillisUsed() + " ms)"
                    + " | nodi B&B=" + partWarm.getNodesExplored()
                    + " | ottimo certificato=" + partWarm.isProvenOptimal()
                    + " | tempo totale (VNS incluso)=" + (t3 - t2) + " ms");
            instance.verifyCutWithBFS(warmResult.cutColors, false);

            boolean agree = Math.abs(baseResult.cutCost - warmResult.cutCost) < 1e-6;
            System.out.println(">>> I due risultati concordano sull'ottimo? " + (agree ? "SI (CORRETTO)" : "NO (ERRORE)"));
        }
    }
}