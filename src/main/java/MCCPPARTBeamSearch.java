import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;

/**
 * Algoritmo EURISTICO di tipo Beam Search per il problema MCstCP pesato,
 * costruito sopra il modello compatto PART_{s-t} (vedi {@link MCCPBranchAndCutPART}).
 *
 * ============================================================================
 * IDEA GENERALE
 * ============================================================================
 * A differenza del Branch-and-Cut (che esplora l'INTERO albero, garantendo
 * l'ottimo ma potenzialmente in tempo esponenziale), il Beam Search esplora
 * il problema LIVELLO PER LIVELLO (un colore aggiunto al taglio per volta),
 * mantenendo ad ogni livello SOLO i k stati parziali piu' promettenti (il
 * "beam", di larghezza k) e scartando tutti gli altri. Questo rende il
 * tempo di esecuzione limitato e prevedibile (O(livelli * k * |C|)
 * valutazioni), al prezzo di perdere la garanzia di ottimalita': e' un
 * algoritmo euristico, non esatto.
 *
 * Due componenti chiave, entrambe richieste esplicitamente:
 *
 * 1) WARM START DA VNS-PROBABILISTIC: prima di iniziare il beam search, si
 *    esegue VNS-Probabilistic per ottenere (a) un incumbent iniziale valido
 *    (mai peggiorato durante la ricerca: il risultato finale non e' MAI
 *    peggiore di quello di VNS da solo), e (b) un ordine di priorita' dei
 *    colori (quelli scelti da VNS vengono provati per primi in fase di
 *    espansione, una euristica di ordinamento che tende a produrre
 *    candidati piu' promettenti prima).
 *
 * 2) VALUTAZIONE VIA MODELLO PART: ogni stato parziale (un sottoinsieme di
 *    colori gia' scelti per il taglio) viene valutato risolvendo il
 *    RILASSAMENTO LINEARE del modello PART_{s-t} con quei colori fissati a
 *    z_c = 1 e tutti gli altri liberi in [0,1] (variabili w libere in
 *    [0,1] come nel modello originale). Il valore ottimo di questo
 *    rilassamento e' un lower bound sul costo di QUALSIASI completamento di
 *    quello stato parziale, e viene usato per classificare e potare i
 *    candidati. Si usa GLOP (il solver LP puro di OR-Tools, senza
 *    branch-and-bound) per queste valutazioni, molto piu' veloce di CBC/SCIP
 *    dato che ogni valutazione e' solo un rilassamento continuo.
 *
 * ============================================================================
 * ALGORITMO
 * ============================================================================
 * 1. Esegui VNS-Probabilistic -> incumbent iniziale (bestCost, bestSolution).
 * 2. beam <- { insieme vuoto di colori }
 * 3. finche' il beam non e' vuoto e non si supera il numero massimo di livelli:
 *    a. per ogni stato nel beam, genera un candidato per ciascun colore non
 *       ancora scelto (aggiungendolo allo stato);
 *    b. per ogni candidato:
 *       - se e' gia' un taglio FEASIBLE (verificato con una BFS indipendente,
 *         non con il modello PART): e' una soluzione completa. Se il suo
 *         costo migliora l'incumbent, aggiorna bestCost/bestSolution. Non
 *         viene espanso ulteriormente (aggiungere altri colori non puo' mai
 *         migliorare un taglio gia' valido, solo peggiorarne il costo);
 *       - altrimenti, calcola il bound LP di quel candidato; se il bound e'
 *         gia' >= bestCost, scartalo (pruning, non potra' mai migliorare
 *         l'incumbent anche completandolo nel modo migliore possibile);
 *         altrimenti tienilo come candidato per il prossimo livello;
 *    c. ordina i candidati superstiti per bound crescente e tieni solo i
 *       migliori k (larghezza del beam); questi diventano il beam del
 *       livello successivo.
 * 4. Restituisci il miglior taglio trovato (mai peggiore dell'incumbent di
 *    VNS, dato che quello resta sempre disponibile come rete di sicurezza).
 *
 * ============================================================================
 * NOTA SULLA SCALABILITA'
 * ============================================================================
 * Ogni livello valuta fino a k * (|C| - livello) candidati, ciascuno con una
 * risoluzione LP. Su istanze con centinaia di colori e beam largo, il numero
 * di risoluzioni LP puo' diventare significativo: usare un beamWidth
 * contenuto (es. 5-20) per istanze grandi, piu' largo (es. 50-100) solo per
 * istanze piccole/medie dove serve piu' qualita' e il tempo lo consente.
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

    // ========================================================================
    // Verifica di feasibility indipendente (BFS), identica nello spirito a
    // MCCPSolver.verifyCutWithBFS: un taglio e' valido se t non e' raggiungibile
    // da s nel sottografo che sopravvive alla rimozione dei colori in cutColors.
    // ========================================================================

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

    // ========================================================================
    // Bound LP: rilassamento continuo del modello PART_{s-t} con i colori in
    // fixedOnes fissati a z_c = 1 (in [1,1]) e gli altri liberi in [0,1].
    // Usa GLOP (LP puro, nessun branch-and-bound): valutazione rapida.
    // ========================================================================

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
            // non dovrebbe succedere (il rilassamento e' sempre feasible: basta
            // z_c=1 ovunque, w_v=1 solo per v=s), ma per sicurezza restituiamo
            // un bound "infinito" cosi' il candidato viene comunque scartato
            return Double.POSITIVE_INFINITY;
        }
        return objective.value();
    }

    // ========================================================================
    // Ordine di priorita' dei colori per l'espansione: prima quelli usati da
    // VNS (che ha gia' dimostrato sapessero produrre un taglio valido), poi
    // gli altri in ordine di costo crescente.
    // ========================================================================

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

    // ========================================================================
    // Beam Search
    // ========================================================================

    private static final class Candidate {
        final Set<Integer> colors;
        final double bound;
        Candidate(Set<Integer> colors, double bound) {
            this.colors = colors;
            this.bound = bound;
        }
    }

    /**
     * @param instance             istanza del problema, usata per il warm start
     *                             (deve rappresentare LO STESSO grafo/costi/s/t
     *                             con cui e' stato costruito questo solver)
     * @param vnsWarmStartTimeMillis tempo (ms) concesso a VNS-Probabilistic per il warm start
     * @param beamWidth            larghezza k del beam (numero di stati mantenuti ad ogni livello)
     */
    public MCCPSolver.MCCPResult solveBeamSearch(MCCPSolver instance, long vnsWarmStartTimeMillis, int beamWidth) {
        if (beamWidth <= 0) {
            throw new IllegalArgumentException("beamWidth deve essere positivo");
        }

        long searchStart = System.currentTimeMillis();
        this.statesEvaluated = 0;
        this.levelsExplored = 0;

        // 1. Warm start VNS-Probabilistic
        System.out.println("[Beam Search] warm start VNS-Probabilistic: budget = " + vnsWarmStartTimeMillis + " ms");
        long tVnsStart = System.currentTimeMillis();
        MCCPSolver.MCCPResult vnsResult = instance.solveProbabilistic(vnsWarmStartTimeMillis);
        long tVnsEnd = System.currentTimeMillis();

        this.vnsWarmStartCost = (long) vnsResult.cutCost;
        this.vnsWarmStartTimeMillisUsed = tVnsEnd - tVnsStart;

        System.out.println("[Beam Search] warm start VNS-Probabilistic -> costo = " + vnsResult.cutCost
                + "  [" + vnsWarmStartTimeMillisUsed + " ms effettivi]");

        double bestCost = vnsResult.cutCost;
        Set<Integer> bestSolution = new HashSet<>(vnsResult.cutColors);

        List<Integer> expansionOrder = buildExpansionOrder(vnsResult.cutColors);

        // 2. Beam Search: si parte da un unico stato, l'insieme vuoto di colori
        List<Set<Integer>> beam = new ArrayList<>();
        beam.add(new HashSet<>());

        while (!beam.isEmpty() && levelsExplored < numColors) {
            levelsExplored++;
            List<Candidate> nextCandidates = new ArrayList<>();

            for (Set<Integer> state : beam) {
                for (int c : expansionOrder) {
                    if (state.contains(c)) continue;

                    Set<Integer> child = new HashSet<>(state);
                    child.add(c);
                    statesEvaluated++;

                    if (isFeasibleCut(child)) {
                        double cost = weight(child);
                        if (cost < bestCost) {
                            bestCost = cost;
                            bestSolution = child;
                            System.out.println("[Beam Search] nuovo incumbent al livello " + levelsExplored
                                    + ": costo = " + bestCost);
                        }
                        // stato terminale: non va espanso ulteriormente
                        continue;
                    }

                    double bound = lpBound(child);
                    if (bound >= bestCost - 1e-9) {
                        continue; // pruning: non puo' mai migliorare l'incumbent
                    }
                    nextCandidates.add(new Candidate(child, bound));
                }
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

        System.out.println("[Beam Search] terminato: " + levelsExplored + " livelli esplorati, "
                + statesEvaluated + " stati valutati, costo finale = " + bestCost);

        return new MCCPSolver.MCCPResult(cutColors, keptColors, bestCost, vnsWarmStartTimeMillisUsed, totalTimeMs);
    }

    // ========================================================================
    // Esempio di utilizzo
    // ========================================================================

    public static void main(String[] args) {
        System.out.println("=== TEST: Beam Search (warm start VNS + bound PART) vs VNS-Greedy da solo ===");

        for (int trial = 0; trial < 4; trial++) {
            int numNodes = 60 + trial * 20;
            int numColors = 15 + trial * 5;

            MCCPSolver instance = MCCPSolver.generateRandomInstance(numNodes, numColors, MCCPSolver.Density.MEDIUM, 0.10);

            MCCPPARTBeamSearch beamSearch = new MCCPPARTBeamSearch(
                    instance.getNumNodes(), instance.getEdges(), instance.getNumColors(),
                    instance.getColorCost(), instance.getSourceNode(), instance.getTargetNode());

            long t0 = System.currentTimeMillis();
            MCCPSolver.MCCPResult beamResult = beamSearch.solveBeamSearch(instance, 5000, 20);
            long t1 = System.currentTimeMillis();

            MCCPSolver.MCCPResult vnsResult = instance.solve(5000);

            System.out.println("\n[Trial " + trial + "] " + numNodes + " nodi, " + numColors
                    + " colori, s=" + instance.getSourceNode() + ", t=" + instance.getTargetNode());

            System.out.println("Beam Search   -> costo=" + beamResult.cutCost
                    + " | warm start VNS costo=" + beamSearch.getVnsWarmStartCost()
                    + " | livelli=" + beamSearch.getLevelsExplored()
                    + " | stati valutati=" + beamSearch.getStatesEvaluated()
                    + " | tempo totale=" + (t1 - t0) + " ms");
            instance.verifyCutWithBFS(beamResult.cutColors, false);

            System.out.println("VNS-Greedy    -> costo=" + vnsResult.cutCost);
            instance.verifyCutWithBFS(vnsResult.cutColors, false);

            boolean beamBetterOrEqual = beamResult.cutCost <= vnsResult.cutCost + 1e-6;
            System.out.println(">>> Beam Search e' migliore o uguale a VNS-Greedy? " + beamBetterOrEqual);
        }
    }
}