import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;

/**
 * Algoritmo ESATTO basato sul MODELLO COMPATTO PART_{s-t} (Silva et al., 2016).
 * Utilizza Google OR-Tools (Solutore CBC/SCIP) per una risoluzione robusta,
 * immune a errori numerici e velocissima.
 */
public class MCCPBranchAndCutPART {

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

    public MCCPBranchAndCutPART(int numNodes, List<MCCPSolver.Edge> edges, int numColors,
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

    // ========================================================================
    // RISOLUZIONE TRAMITE OR-TOOLS (CBC MILP Solver)
    // ========================================================================

    public MCCPSolver.MCCPResult solveExact(long maxRunningTimeMillis) {
        long searchStartMillis = System.currentTimeMillis();

        // Istanzia il solutore MILP CBC (Coin-or OR Branch and Cut)
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

        // w binaria (0/1): il nodo v sta dal lato di s (1) o dal lato di t (0)
        MPVariable[] w = new MPVariable[numNodes];
        for (int v = 0; v < numNodes; v++) {
            w[v] = solver.makeBoolVar("w_" + v);
        }

        // z continua, solo z_i >= 0 (nessun upper bound esplicito necessario)
        MPVariable[] z = new MPVariable[numColors];
        for (int c = 0; c < numColors; c++) {
            z[c] = solver.makeNumVar(0.0, Double.POSITIVE_INFINITY, "z_" + c);
        }

        // 3. Vincoli di Confine per i nodi sorgente (s) e pozzo (t)
        // w_s = 1 (s appartiene al lato S)
        MPConstraint sConstraint = solver.makeConstraint(1.0, 1.0, "w_s_bound");
        sConstraint.setCoefficient(w[s], 1.0);

        // w_t = 0 (t appartiene al lato T)
        MPConstraint tConstraint = solver.makeConstraint(0.0, 0.0, "w_t_bound");
        tConstraint.setCoefficient(w[t], 1.0);

        int totalConstraints = 2;

        // 4. Vincoli d'Arco del Modello PART_{s-t} (Silva et al., 2016)
        // Per ogni arco e = (u,v):
        //   (1) sum_{c in e.colors} z_c - w_u + w_v >= 0
        //   (2) sum_{c in e.colors} z_c + w_u - w_v >= 0
        for (MCCPSolver.Edge e : edges) {
            // Direzione u -> v
            MPConstraint c1 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
            for (int c : e.colors) {
                c1.setCoefficient(z[c], 1.0);
            }
            c1.setCoefficient(w[e.u], -1.0);
            c1.setCoefficient(w[e.v], 1.0);
            totalConstraints++;

            // Direzione v -> u
            MPConstraint c2 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
            for (int c : e.colors) {
                c2.setCoefficient(z[c], 1.0);
            }
            c2.setCoefficient(w[e.u], 1.0);
            c2.setCoefficient(w[e.v], -1.0);
            totalConstraints++;
        }

        this.rootEdgeConstraintCount = totalConstraints;

        // 5. Funzione Obiettivo: Minimizzare la somma dei costi dei colori tagliati
        MPObjective objective = solver.objective();
        for (int c = 0; c < numColors; c++) {
            objective.setCoefficient(z[c], colorCost[c]);
        }
        objective.setMinimization();

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

    // ========================================================================
    // Main di Test e Benchmark sullo stesso Grafo
    // ========================================================================

    public static void main(String[] args) {
        System.out.println("=== TEST BENCHMARK: B&C-Base vs B&C-PART (OR-Tools) ===");

        for (int trial = 0; trial < 5; trial++) {
            int numNodes = 80 + trial * 10;
            int numColors = 13 + 5 * trial;

            // Genera l'istanza del grafo condivisa (lo stesso grafo viene passato a entrambi)
            MCCPSolver instance = MCCPSolver.generateRandomInstance(numNodes, numColors, MCCPSolver.Density.MEDIUM, 0.10);

            // 1. Branch & Cut BASE
            MCCPBranchAndCut bc = new MCCPBranchAndCut(
                    instance.getNumNodes(),
                    instance.getEdges(),
                    instance.getNumColors(),
                    instance.getColorCost(),
                    instance.getSourceNode(),
                    instance.getTargetNode()
            );
            long t0 = System.currentTimeMillis();
            MCCPSolver.MCCPResult bcResult = bc.solveExact(200000);
            long t1 = System.currentTimeMillis();

            // 2. Modello PART_{s-t} tramite OR-Tools
            MCCPBranchAndCutPART part = new MCCPBranchAndCutPART(
                    instance.getNumNodes(),
                    instance.getEdges(),
                    instance.getNumColors(),
                    instance.getColorCost(),
                    instance.getSourceNode(),
                    instance.getTargetNode()
            );
            long t2 = System.currentTimeMillis();
            MCCPSolver.MCCPResult partResult = part.solveExact(200000);
            long t3 = System.currentTimeMillis();

            System.out.println("\n[Trial " + trial + "] Grafo: " + numNodes + " Nodi, " + numColors + " Colori | s=" + instance.getSourceNode() + ", t=" + instance.getTargetNode());

            System.out.println("B&C-Base      -> Costo: " + bcResult.cutCost + " | Nodi B&B: " + bc.getNodesExplored() + " | Tempo: " + (t1 - t0) + " ms");
            instance.verifyCutWithBFS(bcResult.cutColors, false);

            System.out.println("B&C-PART (OR) -> Costo: " + partResult.cutCost + " | Nodi B&B: " + part.getNodesExplored() + " | Tempo: " + (t3 - t2) + " ms");
            instance.verifyCutWithBFS(partResult.cutColors, false);

            boolean agree = Math.abs(bcResult.cutCost - partResult.cutCost) < 1e-6;
            System.out.println(">>> I due algoritmi concordano sull'ottimo? " + (agree ? "SI (CORRETTO)" : "NO (ERRORE)"));
        }
    }
}