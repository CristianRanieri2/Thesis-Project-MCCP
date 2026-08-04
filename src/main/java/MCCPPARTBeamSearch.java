import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;

/**
 * Algoritmo EURISTICO di tipo Beam Search per il problema MCstCP pesato,
 * costruito sopra il modello compatto PART_{s-t} con controlli temporali sia
 * sulla fase di Warm Start VNS che sulla fase di Beam Search.
 */
public class MCCPPARTBeamSearch {

    static {
        Loader.loadNativeLibraries();
    }

    private final int numNodes;
    private final List<MCCPSolver.Edge> edges;
    private final int numColors;
    private final double[] colorCost;
    private final int s;
    private final int t;

    private long vnsWarmStartCost = -1;
    private long vnsWarmStartTimeMillisUsed = 0;
    private long statesEvaluated = 0;
    private int levelsExplored = 0;

    public MCCPPARTBeamSearch(int numNodes, List<MCCPSolver.Edge> edges, int numColors,
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

    public long getVnsWarmStartCost() { return vnsWarmStartCost; }
    public long getVnsWarmStartTimeMillisUsed() { return vnsWarmStartTimeMillisUsed; }
    public long getStatesEvaluated() { return statesEvaluated; }
    public int getLevelsExplored() { return levelsExplored; }

    private boolean isFeasibleCut(Set<Integer> cutColors) {
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
        return !visited[t];
    }

    private double weight(Set<Integer> colors) {
        double total = 0.0;
        for (int c : colors) total += colorCost[c];
        return total;
    }

    private double lpBound(Set<Integer> fixedOnes) {
        MPSolver solver = MPSolver.createSolver("GLOP");
        if (solver == null) {
            throw new RuntimeException("Impossibile caricare il solutore LP GLOP di OR-Tools.");
        }

        MPVariable[] z = new MPVariable[numColors];
        for (int c = 0; c < numColors; c++) {
            if (fixedOnes.contains(c)) {
                z[c] = solver.makeNumVar(1.0, 1.0, "z_" + c);
            } else {
                z[c] = solver.makeNumVar(0.0, 1.0, "z_" + c);
            }
        }

        MPVariable[] w = new MPVariable[numNodes];
        for (int v = 0; v < numNodes; v++) {
            w[v] = solver.makeNumVar(0.0, 1.0, "w_" + v);
        }

        MPConstraint sConstraint = solver.makeConstraint(1.0, 1.0, "w_s_bound");
        sConstraint.setCoefficient(w[s], 1.0);
        MPConstraint tConstraint = solver.makeConstraint(0.0, 0.0, "w_t_bound");
        tConstraint.setCoefficient(w[t], 1.0);

        for (MCCPSolver.Edge e : edges) {
            MPConstraint c1 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
            for (int c : e.colors) c1.setCoefficient(z[c], 1.0);
            c1.setCoefficient(w[e.u], -1.0);
            c1.setCoefficient(w[e.v], 1.0);

            MPConstraint c2 = solver.makeConstraint(0.0, Double.POSITIVE_INFINITY);
            for (int c : e.colors) c2.setCoefficient(z[c], 1.0);
            c2.setCoefficient(w[e.u], 1.0);
            c2.setCoefficient(w[e.v], -1.0);
        }

        MPObjective objective = solver.objective();
        for (int c = 0; c < numColors; c++) {
            objective.setCoefficient(z[c], colorCost[c]);
        }
        objective.setMinimization();

        MPSolver.ResultStatus status = solver.solve();
        if (status != MPSolver.ResultStatus.OPTIMAL) {
            return Double.POSITIVE_INFINITY;
        }
        return objective.value();
    }

    private List<Integer> buildExpansionOrder(Set<Integer> vnsCutColors) {
        List<Integer> order = new ArrayList<>();
        List<Integer> vnsFirst = new ArrayList<>(vnsCutColors);
        vnsFirst.sort(Comparator.comparingDouble(c -> colorCost[c]));
        order.addAll(vnsFirst);

        List<Integer> rest = new ArrayList<>();
        for (int c = 0; c < numColors; c++) {
            if (!vnsCutColors.contains(c)) rest.add(c);
        }
        rest.sort(Comparator.comparingDouble(c -> colorCost[c]));
        order.addAll(rest);

        return order;
    }

    private static final class Candidate {
        final Set<Integer> colors;
        final double bound;
        Candidate(Set<Integer> colors, double bound) {
            this.colors = colors;
            this.bound = bound;
        }
    }

    /**
     * Esegue l'algoritmo Beam Search vincolato sia dal tempo per il Warm Start che dal timeout del Beam Search.
     *
     * @param instance             Istanza del problema per il Warm Start
     * @param vnsWarmStartTimeMs   Tempo max (ms) concesso a VNS-Probabilistic (100%)
     * @param beamSearchTimeoutMs  Tempo max (ms) concesso alla sola fase di Beam Search (100%)
     * @param beamWidth            Larghezza k del fascio
     */
    public MCCPSolver.MCCPResult solveBeamSearch(MCCPSolver instance, long vnsWarmStartTimeMs, long beamSearchTimeoutMs, int beamWidth) {
        if (beamWidth <= 0) {
            throw new IllegalArgumentException("beamWidth deve essere positivo");
        }

        long searchStart = System.currentTimeMillis();
        this.statesEvaluated = 0;
        this.levelsExplored = 0;

        // 1. Warm start VNS-Probabilistic
        System.out.println("[Beam Search] Warm start VNS-Probabilistic: budget = " + vnsWarmStartTimeMs + " ms");
        long tVnsStart = System.currentTimeMillis();
        MCCPSolver.MCCPResult vnsResult = instance.solveProbabilistic(vnsWarmStartTimeMs);
        long tVnsEnd = System.currentTimeMillis();

        this.vnsWarmStartCost = (long) vnsResult.cutCost;
        this.vnsWarmStartTimeMillisUsed = tVnsEnd - tVnsStart;

        System.out.println("[Beam Search] Warm start VNS-Probabilistic -> costo = " + vnsResult.cutCost
                + "  [" + vnsWarmStartTimeMillisUsed + " ms effettivi]");

        double bestCost = vnsResult.cutCost;
        Set<Integer> bestSolution = new HashSet<>(vnsResult.cutColors);

        List<Integer> expansionOrder = buildExpansionOrder(vnsResult.cutColors);

        // 2. Beam Search
        System.out.println("[Beam Search] Avvio fase Beam Search: budget = " + beamSearchTimeoutMs + " ms");
        long beamSearchStart = System.currentTimeMillis();

        List<Set<Integer>> beam = new ArrayList<>();
        beam.add(new HashSet<>());

        while (!beam.isEmpty() && levelsExplored < numColors) {
            // Controllo limite di tempo all'inizio del livello
            if (System.currentTimeMillis() - beamSearchStart >= beamSearchTimeoutMs) {
                System.out.println("[Beam Search] Timeout raggiunto (" + beamSearchTimeoutMs + " ms). Interruzione ricerca.");
                break;
            }

            levelsExplored++;
            List<Candidate> nextCandidates = new ArrayList<>();
            boolean timeoutReached = false;

            for (Set<Integer> state : beam) {
                if (timeoutReached) break;

                for (int c : expansionOrder) {
                    // Controllo capillare del timeout durante l'espansione dei candidati
                    if (System.currentTimeMillis() - beamSearchStart >= beamSearchTimeoutMs) {
                        timeoutReached = true;
                        break;
                    }

                    if (state.contains(c)) continue;

                    Set<Integer> child = new HashSet<>(state);
                    child.add(c);
                    statesEvaluated++;

                    if (isFeasibleCut(child)) {
                        double cost = weight(child);
                        if (cost < bestCost) {
                            bestCost = cost;
                            bestSolution = child;
                            System.out.println("[Beam Search] Nuovo incumbent al livello " + levelsExplored
                                    + ": costo = " + bestCost);
                        }
                        continue; // Stato terminale
                    }

                    double bound = lpBound(child);
                    if (bound >= bestCost - 1e-9) {
                        continue; // Pruning by bound
                    }
                    nextCandidates.add(new Candidate(child, bound));
                }
            }

            if (timeoutReached) {
                System.out.println("[Beam Search] Timeout raggiunto durante la valutazione del livello " + levelsExplored + ".");
                break;
            }

            nextCandidates.sort(Comparator.comparingDouble(cand -> cand.bound));
            beam = new ArrayList<>();
            for (int i = 0; i < Math.min(beamWidth, nextCandidates.size()); i++) {
                beam.add(nextCandidates.get(i).colors);
            }
        }

        long totalTimeMs = System.currentTimeMillis() - searchStart;

        Set<Integer> cutColors = new HashSet<>(bestSolution);
        Set<Integer> keptColors = new HashSet<>();
        for (int c = 0; c < numColors; c++) {
            if (!cutColors.contains(c)) keptColors.add(c);
        }

        System.out.println("[Beam Search] Terminata: " + levelsExplored + " livelli esplorati, "
                + statesEvaluated + " stati valutati, costo finale = " + bestCost);

        return new MCCPSolver.MCCPResult(cutColors, keptColors, bestCost, vnsWarmStartTimeMillisUsed, totalTimeMs);
    }
}