import java.util.*;

/**
 * Implementazione basata su Variable Neighborhood Search (VNS) per una
 * variante PESATA del Minimum Color s-t Cut Problem (MCstCP), a partire
 * dall'Algoritmo 1 (algoritmo generale) descritto in:
 * "New algorithms for the Minimum Coloring Cut Problem" (Bordini & Protti, 2017).
 *
 * DIFFERENZE RISPETTO AL PROBLEMA ORIGINALE DEL PAPER (MCCP)
 * ------------------------------------------------------------------
 * 1) Nel paper l'obiettivo è disconnettere l'INTERO grafo con il minor
 *    numero di colori. Qui invece si vogliono separare due nodi specifici,
 *    sorgente s e destinazione t (variante MCstCP, già citata nel paper
 *    in Introduzione come "Minimum Color s-t Cut Problem", Coudert et al. 2007).
 *    "Feasible" quindi non significa più "il grafo è disconnesso", ma
 *    "s e t si trovano in componenti diverse".
 * 2) Ogni colore ha un COSTO di utilizzo (colorCost[c]). L'obiettivo non è
 *    più minimizzare il NUMERO di colori nel taglio, ma il COSTO TOTALE
 *    dei colori scelti per il taglio che separa s da t.
 *
 * Come nel paper, si segue la strategia "duale": invece di cercare
 * direttamente l'insieme di colori da rimuovere (il taglio), si cerca
 * l'insieme di colori da MANTENERE (BestS) che:
 *   - lascia s e t in componenti diverse (feasible),
 *   - massimizza il PESO totale (somma dei costi) dei colori mantenuti.
 * Il taglio finale è C \ BestS, e il suo costo è costoTotale(C) - peso(BestS)
 * (minimizzare il costo del taglio equivale quindi a massimizzare il peso
 * di BestS, esattamente come nel paper si massimizzava |BestS| invece di
 * minimizzare |C \ BestS|).
 *
 * ADATTAMENTO DELLE SUBROUTINE (Generate-Initial-Solution, New-Solution,
 * Local-Search)
 * ------------------------------------------------------------------
 * Nel paper, ad ogni passo si sceglie il colore c che MASSIMIZZA
 * Number-of-Components(BestS U {c}) fra TUTTI i candidati, e solo DOPO si
 * verifica se quella scelta è feasible (altrimenti ci si ferma). Questo
 * funziona nel paper perché "più componenti" è correlato empiricamente
 * con "restare disconnessi": il candidato migliore secondo quella metrica
 * è tipicamente anche quello più sicuro dal punto di vista della
 * feasibility.
 *
 * Con un costo arbitrario per colore questa correlazione NON esiste: il
 * colore più costoso può benissimo essere quello che, se aggiunto,
 * ricongiunge s e t, mentre un colore più economico potrebbe restare
 * feasible. Applicare la stessa struttura "scegli il migliore, poi
 * verifica" causerebbe quindi terminazioni premature e soluzioni molto
 * peggiori. Per questo, in questa implementazione, la selezione del
 * prossimo colore da aggiungere è ristretta FIN DA SUBITO ai soli
 * candidati che mantengono la feasibility, e fra questi si sceglie quello
 * di costo massimo (versione greedy) oppure lo si estrae con la funzione
 * di Boltzmann sui costi (versione probabilistica). La struttura generale
 * dell'algoritmo (Algoritmo 1) e delle subroutine Shake/Fix resta
 * invariata rispetto al paper.
 *
 * ESTENSIONE: archi con uno o più colori (multi-color / "label set")
 * ------------------------------------------------------------------
 * Ogni arco può avere un INSIEME di colori invece di uno solo (es. un
 * collegamento condiviso da più operatori o soggetto a più "risk group",
 * cfr. applicazione "Shared Risk Link Group" citata nell'Introduzione del
 * paper). Un arco e è presente nel sottografo indotto da un insieme
 * mantenuto C' se e solo se TUTTI i suoi colori appartengono a C'
 * (Colors(e) ⊆ C'): basta che UN SOLO colore assegnato all'arco venga
 * rimosso perché l'arco cada (semantica SRLG-like). Per un arco
 * monocolore (caso originale del paper) questo equivale esattamente a
 * colorSet.contains(e.color).
 */
public class MCCPSolver {

    // ---------- Rappresentazione del grafo ----------

    static class Edge {
        final int u, v;
        final Set<Integer> colors;

        Edge(int u, int v, int... colors) {
            this.u = u;
            this.v = v;
            this.colors = new HashSet<>();
            for (int c : colors) this.colors.add(c);
        }

        Edge(int u, int v, Collection<Integer> colors) {
            this.u = u;
            this.v = v;
            this.colors = new HashSet<>(colors);
        }
    }

    private final int numNodes;
    private final List<Edge> edges;
    private final int numColors;      // colori numerati da 0 a numColors-1
    private final double[] colorCost; // costo di utilizzo di ciascun colore
    private final int s;              // nodo sorgente
    private final int t;              // nodo destinazione

    /**
     * @param numNodes   numero di nodi del grafo (indicizzati da 0 a numNodes-1)
     * @param edges      lista degli archi (ognuno con uno o più colori)
     * @param numColors  numero totale di colori (indicizzati da 0 a numColors-1)
     * @param colorCost  costo di utilizzo di ciascun colore, colorCost[c] per il colore c
     * @param sourceNode nodo sorgente (s) da separare dal nodo destinazione
     * @param targetNode nodo destinazione (t) da separare dal nodo sorgente
     */
    public MCCPSolver(int numNodes, List<Edge> edges, int numColors,
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

    // ---------- Union-Find (disjoint-set) con path compression e union by rank ----------

    static class UnionFind {
        final int[] parent;
        final int[] rank;

        UnionFind(int n) {
            parent = new int[n];
            rank = new int[n];
            for (int i = 0; i < n; i++) parent[i] = i;
        }

        int find(int x) {
            while (parent[x] != x) {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        }

        void union(int a, int b) {
            int ra = find(a), rb = find(b);
            if (ra == rb) return;
            if (rank[ra] < rank[rb]) { int tmp = ra; ra = rb; rb = tmp; }
            parent[rb] = ra;
            if (rank[ra] == rank[rb]) rank[ra]++;
        }
    }

    /**
     * Restituisce true se i nodi a e b appartengono alla stessa componente
     * connessa del sottografo indotto dai colori in colorSet.
     */
    public boolean connected(int a, int b, Set<Integer> colorSet) {
        UnionFind uf = new UnionFind(numNodes);
        for (Edge e : edges) {
            if (colorSet.containsAll(e.colors)) {
                uf.union(e.u, e.v);
            }
        }
        return uf.find(a) == uf.find(b);
    }

    /**
     * Una soluzione C' (insieme di colori mantenuti) è "feasible" se
     * mantenendo solo quei colori, s e t risultano in componenti diverse.
     */
    public boolean isFeasible(Set<Integer> colorSet) {
        return !connected(s, t, colorSet);
    }

    /** Peso (somma dei costi) di un insieme di colori. */
    public double weight(Set<Integer> colorSet) {
        double total = 0.0;
        for (int c : colorSet) total += colorCost[c];
        return total;
    }

    // ---------- Selezione greedy del miglior candidato feasible ----------

    private boolean tryAddBestFeasibleColor(Set<Integer> base, Set<Integer> candidatePool) {
        int bestColor = -1;
        double bestCost = Double.NEGATIVE_INFINITY;

        for (int c : candidatePool) {
            Set<Integer> candidate = new HashSet<>(base);
            candidate.add(c);
            if (isFeasible(candidate) && colorCost[c] > bestCost) {
                bestCost = colorCost[c];
                bestColor = c;
            }
        }

        if (bestColor == -1) return false;
        base.add(bestColor);
        return true;
    }

    // ---------- 2.1 Generate-Initial-Solution (versione greedy, adattata dall'Algoritmo 2) ----------

    public Set<Integer> generateInitialSolutionGreedy() {
        Set<Integer> bestS = new HashSet<>();

        boolean progressed = true;
        while (progressed) {
            Set<Integer> remaining = diff(allColors(), bestS);
            progressed = tryAddBestFeasibleColor(bestS, remaining);
        }
        return bestS;
    }

    // ---------- 2.2 New-Solution (versione greedy, adattata dall'Algoritmo 4) ----------

    public Set<Integer> newSolutionGreedy(Set<Integer> bestS) {
        Set<Integer> s = new HashSet<>();

        boolean progressed = true;
        while (isFeasible(s) && !diff(bestS, s).isEmpty() && progressed) {
            progressed = tryAddBestFeasibleColor(s, diff(bestS, s));
        }

        progressed = true;
        while (isFeasible(s) && progressed) {
            progressed = tryAddBestFeasibleColor(s, diff(bestS, s));
        }

        return s;
    }

    // ---------- 2.3 Shake (Algoritmo 6, invariato) ----------

    public void shake(Set<Integer> sPrime, int k, Set<Integer> currentS, Random rnd) {
        for (int i = 0; i < k; i++) {
            double delta = rnd.nextDouble();
            if (delta < 0.5 && !sPrime.isEmpty()) {
                List<Integer> intersection = new ArrayList<>(sPrime);
                intersection.retainAll(currentS);
                if (!intersection.isEmpty()) {
                    int toRemove = intersection.get(rnd.nextInt(intersection.size()));
                    sPrime.remove(toRemove);
                }
            } else {
                List<Integer> candidates = new ArrayList<>();
                for (int c = 0; c < numColors; c++) {
                    if (!currentS.contains(c) && !sPrime.contains(c)) {
                        candidates.add(c);
                    }
                }
                if (!candidates.isEmpty()) {
                    int toAdd = candidates.get(rnd.nextInt(candidates.size()));
                    sPrime.add(toAdd);
                }
            }
        }
    }

    // ---------- 2.4 Fix (Algoritmo 7, adattato) ----------

    public void fix(Set<Integer> sPrime, Random rnd) {
        while (!isFeasible(sPrime)) {
            if (sPrime.isEmpty()) break;
            List<Integer> list = new ArrayList<>(sPrime);
            int toRemove = list.get(rnd.nextInt(list.size()));
            sPrime.remove(toRemove);
        }
    }

    // ---------- 2.5 Local-Search (versione greedy, adattata dall'Algoritmo 8) ----------

    public void localSearchGreedy(Set<Integer> sPrime) {
        while (isFeasible(sPrime)) {
            Set<Integer> complement = diff(allColors(), sPrime);
            if (!tryAddBestFeasibleColor(sPrime, complement)) break;
        }
    }

    // ---------- Funzione di Boltzmann (usata dalle versioni probabilistiche) ----------
    private static final double TEMPERATURE = 1.0;

    private Integer pickProbabilistic(Set<Integer> base, Set<Integer> candidatePool, Random rnd) {
        if (candidatePool.isEmpty()) return null;

        double maxCost = Double.NEGATIVE_INFINITY;
        for (int c : candidatePool) {
            if (colorCost[c] > maxCost) maxCost = colorCost[c];
        }

        List<Integer> feasibleCandidates = new ArrayList<>();
        List<Double> weights = new ArrayList<>();
        double totalWeight = 0.0;

        for (int c : candidatePool) {
            Set<Integer> candidate = new HashSet<>(base);
            candidate.add(c);
            if (isFeasible(candidate)) {
                double delta = colorCost[c] - maxCost;
                double w = Math.exp(delta / TEMPERATURE);
                feasibleCandidates.add(c);
                weights.add(w);
                totalWeight += w;
            }
        }

        if (feasibleCandidates.isEmpty()) return null;

        double r = rnd.nextDouble() * totalWeight;
        double cumulative = 0.0;
        for (int i = 0; i < feasibleCandidates.size(); i++) {
            cumulative += weights.get(i);
            if (r <= cumulative) {
                return feasibleCandidates.get(i);
            }
        }
        return feasibleCandidates.get(feasibleCandidates.size() - 1);
    }

    // ---------- 2.1 Generate-Initial-Solution (versione probabilistica, adattata dall'Algoritmo 3) ----------

    public Set<Integer> generateInitialSolutionProbabilistic(Random rnd) {
        Set<Integer> bestS = new HashSet<>();
        Set<Integer> remaining = allColors();

        boolean endLoop = false;
        while (!endLoop) {
            Integer chosen = pickProbabilistic(bestS, remaining, rnd);
            if (chosen != null) {
                bestS.add(chosen);
                remaining.remove(chosen);
            } else {
                endLoop = true;
            }
        }
        return bestS;
    }

    // ---------- 2.2 New-Solution (versione probabilistica, adattata dall'Algoritmo 5) ----------

    public Set<Integer> newSolutionProbabilistic(Set<Integer> bestS, Random rnd) {
        Set<Integer> s = new HashSet<>();

        boolean progressed = true;
        while (isFeasible(s) && !diff(bestS, s).isEmpty() && progressed) {
            Integer chosen = pickProbabilistic(s, diff(bestS, s), rnd);
            if (chosen != null) {
                s.add(chosen);
            } else {
                progressed = false;
            }
        }

        progressed = true;
        while (isFeasible(s) && progressed) {
            Integer chosen = pickProbabilistic(s, diff(bestS, s), rnd);
            if (chosen != null) {
                s.add(chosen);
            } else {
                progressed = false;
            }
        }

        return s;
    }

    // ---------- 2.5 Local-Search (versione probabilistica, adattata dall'Algoritmo 9) ----------

    public void localSearchProbabilistic(Set<Integer> sPrime, Random rnd) {
        while (isFeasible(sPrime)) {
            Set<Integer> complement = diff(allColors(), sPrime);
            Integer chosen = pickProbabilistic(sPrime, complement, rnd);
            if (chosen != null) {
                sPrime.add(chosen);
            } else {
                break;
            }
        }
    }

    // ---------- Risultato restituito da solve() / solveProbabilistic() ----------

    public static class MCCPResult {
        public final double cutCost;
        public final int cutSize;
        public final Set<Integer> cutColors;
        public final Set<Integer> keptColors;

        MCCPResult(Set<Integer> cutColors, Set<Integer> keptColors, double cutCost) {
            this.cutColors = Collections.unmodifiableSet(new TreeSet<>(cutColors));
            this.keptColors = Collections.unmodifiableSet(new TreeSet<>(keptColors));
            this.cutSize = cutColors.size();
            this.cutCost = cutCost;
        }

        @Override
        public String toString() {
            return "costo del taglio = " + cutCost + ", colori nel taglio = " + cutColors
                    + " (" + cutSize + " colori; colori mantenuti: " + keptColors + ")";
        }
    }

    private MCCPResult buildResult(Set<Integer> bestS) {
        Set<Integer> cutColors = diff(allColors(), bestS);
        double cutCost = weight(cutColors);
        return new MCCPResult(cutColors, bestS, cutCost);
    }

    // ---------- Algoritmo 1: General algorithm ----------

    public MCCPResult solve(long maxRunningTimeMillis) {
        Random rnd = new Random();
        long startTime = System.currentTimeMillis();

        Set<Integer> bestS = generateInitialSolutionGreedy();
        int maxNeighborhood = numColors - bestS.size();

        do {
            Set<Integer> s = newSolutionGreedy(bestS);

            while (weight(s) > weight(bestS)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
                s = newSolutionGreedy(bestS);
            }

            int k = 1;

            while (k < maxNeighborhood) {
                Set<Integer> sPrime = new HashSet<>(s);
                shake(sPrime, k, s, rnd);

                if (!isFeasible(sPrime)) {
                    fix(sPrime, rnd);
                }

                localSearchGreedy(sPrime);

                if (weight(sPrime) > weight(s)) {
                    s = sPrime;
                    k = 1;
                } else {
                    k++;
                }

                if (System.currentTimeMillis() - startTime > maxRunningTimeMillis) break;
            }

            if (weight(s) > weight(bestS)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
            }

        } while (System.currentTimeMillis() - startTime <= maxRunningTimeMillis);

        return buildResult(bestS);
    }

    public MCCPResult solveProbabilistic(long maxRunningTimeMillis) {
        Random rnd = new Random();
        long startTime = System.currentTimeMillis();

        Set<Integer> bestS = generateInitialSolutionProbabilistic(rnd);
        int maxNeighborhood = numColors - bestS.size();

        do {
            Set<Integer> s = newSolutionProbabilistic(bestS, rnd);

            while (weight(s) > weight(bestS)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
                s = newSolutionProbabilistic(bestS, rnd);
            }

            int k = 1;

            while (k < maxNeighborhood) {
                Set<Integer> sPrime = new HashSet<>(s);
                shake(sPrime, k, s, rnd);

                if (!isFeasible(sPrime)) {
                    fix(sPrime, rnd);
                }

                localSearchProbabilistic(sPrime, rnd);

                if (weight(sPrime) > weight(s)) {
                    s = sPrime;
                    k = 1;
                } else {
                    k++;
                }

                if (System.currentTimeMillis() - startTime > maxRunningTimeMillis) break;
            }

            if (weight(s) > weight(bestS)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
            }

        } while (System.currentTimeMillis() - startTime <= maxRunningTimeMillis);

        return buildResult(bestS);
    }

    // ---------- Validatore esatto (solo per un numero contenuto di colori) ----------

    /**
     * Calcola la soluzione ESATTA per enumerazione di tutti i 2^numColors
     * sottoinsiemi di colori. Utile per validare il risultato della VNS su
     * istanze di test con pochi colori (indicativamente numColors <= 22-24).
     * Non pensato per l'uso in produzione su istanze reali (crescita esponenziale).
     */
    public MCCPResult bruteForceOptimal() {
        if (numColors > 24) {
            throw new IllegalStateException(
                    "bruteForceOptimal e' pensato solo per istanze di test con pochi colori (<= 24)");
        }
        Set<Integer> bestKept = new HashSet<>();
        double bestWeight = -1.0;
        int totalMasks = 1 << numColors;
        for (int mask = 0; mask < totalMasks; mask++) {
            Set<Integer> candidate = new HashSet<>();
            for (int c = 0; c < numColors; c++) {
                if ((mask & (1 << c)) != 0) candidate.add(c);
            }
            if (isFeasible(candidate)) {
                double w = weight(candidate);
                if (w > bestWeight) {
                    bestWeight = w;
                    bestKept = candidate;
                }
            }
        }
        return buildResult(bestKept);
    }

    // ---------- utilità ----------

    private Set<Integer> allColors() {
        Set<Integer> all = new HashSet<>();
        for (int c = 0; c < numColors; c++) all.add(c);
        return all;
    }

    private static Set<Integer> diff(Set<Integer> a, Set<Integer> b) {
        Set<Integer> result = new HashSet<>(a);
        result.removeAll(b);
        return result;
    }

    // ---------- esempio di utilizzo ----------

    public static void main(String[] args) {
        runWeightedStExample();
        System.out.println();
        runLargeRandomExample(
                "Esempio grande #1 (50 nodi)",
                50,      // numNodes
                10,      // numColors
                90,      // numExtraEdges (oltre allo spanning tree, che garantisce la connettivita')
                0.10,    // probabilita' che un arco extra sia multi-colore
                0, 49,   // sourceNode, targetNode
                12345L   // seed (per riproducibilita')
        );
        System.out.println();
        runLargeRandomExample(
                "Esempio grande #2 (50 nodi)",
                50,
                12,
                110,
                0.15,
                5, 44,
                67890L
        );
    }

    private static void runWeightedStExample() {
        List<Edge> edges = new ArrayList<>();

        edges.add(new Edge(0, 1, 0));
        edges.add(new Edge(0, 2, 1));
        edges.add(new Edge(1, 2, 2));
        edges.add(new Edge(1, 3, 0));
        edges.add(new Edge(3, 2, 1));

        edges.add(new Edge(4, 5, 0));
        edges.add(new Edge(4, 6, 1));
        edges.add(new Edge(5, 6, 2));
        edges.add(new Edge(5, 7, 0));
        edges.add(new Edge(7, 6, 1));

        edges.add(new Edge(0, 4, 2));
        edges.add(new Edge(2, 6, 1));
        edges.add(new Edge(3, 7, 2));

        double[] colorCost = { 5.0, 2.0, 8.0 }; // colore1=5, colore2=2, colore3=8

        int sourceNode = 0;
        int targetNode = 6;

        MCCPSolver solver = new MCCPSolver(8, edges, 3, colorCost, sourceNode, targetNode);

        MCCPResult resultGreedy = solver.solve(500);
        MCCPResult resultProbabilistic = solver.solveProbabilistic(500);

        System.out.println("[Minimum Color s-t Cut pesato: separare il nodo " + sourceNode
                + " dal nodo " + targetNode + "]");
        System.out.println("VNS-Greedy        -> " + resultGreedy);
        System.out.println("VNS-Probabilistic -> " + resultProbabilistic);
        System.out.println("Atteso (verificato per enumerazione esaustiva su questo grafo/costi): "
                + "costo = 7.0, colori nel taglio = [0, 1] (colore1 + colore2)");
    }

    /**
     * Genera un grafo casuale ma CONNESSO: prima uno spanning tree casuale
     * (che garantisce la connettivita' complessiva), poi numExtraEdges archi
     * aggiuntivi tra coppie di nodi casuali, ciascuno con probabilita'
     * multiColorProb di avere un secondo colore (arco multi-colore).
     */
    private static List<Edge> generateRandomGraph(int numNodes, int numColors, int numExtraEdges,
                                                  double multiColorProb, Random rnd) {
        List<Edge> edges = new ArrayList<>();

        // 1) spanning tree casuale: connette il nodo i-esimo (in un ordine casuale)
        // a un nodo scelto a caso fra quelli gia' inseriti, cosi' il grafo finale
        // e' garantito connesso.
        List<Integer> order = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) order.add(i);
        Collections.shuffle(order, rnd);

        for (int i = 1; i < numNodes; i++) {
            int u = order.get(i);
            int v = order.get(rnd.nextInt(i));
            int color = rnd.nextInt(numColors);
            edges.add(new Edge(u, v, color));
        }

        // 2) archi extra casuali, per rendere il grafo piu' denso e interessante
        for (int i = 0; i < numExtraEdges; i++) {
            int u = rnd.nextInt(numNodes);
            int v = rnd.nextInt(numNodes);
            if (u == v) continue;
            int color1 = rnd.nextInt(numColors);
            if (rnd.nextDouble() < multiColorProb && numColors > 1) {
                int color2;
                do { color2 = rnd.nextInt(numColors); } while (color2 == color1);
                edges.add(new Edge(u, v, color1, color2));
            } else {
                edges.add(new Edge(u, v, color1));
            }
        }
        return edges;
    }

    /**
     * Genera un'istanza casuale di dimensioni maggiori (es. 50 nodi), la
     * risolve con VNS-Greedy e VNS-Probabilistic, e — dato che il numero di
     * colori resta piccolo (10-12) — calcola anche la soluzione ESATTA per
     * confronto tramite bruteForceOptimal().
     */
    private static void runLargeRandomExample(String title, int numNodes, int numColors,
                                              int numExtraEdges, double multiColorProb,
                                              int sourceNode, int targetNode, long seed) {
        Random rnd = new Random(seed);

        List<Edge> edges = generateRandomGraph(numNodes, numColors, numExtraEdges, multiColorProb, rnd);

        double[] colorCost = new double[numColors];
        for (int c = 0; c < numColors; c++) {
            colorCost[c] = 1 + rnd.nextInt(20); // costo intero fra 1 e 20
        }

        MCCPSolver solver = new MCCPSolver(numNodes, edges, numColors, colorCost, sourceNode, targetNode);

        long t0 = System.currentTimeMillis();
        MCCPResult resultGreedy = solver.solve(1000);
        long t1 = System.currentTimeMillis();
        MCCPResult resultProbabilistic = solver.solveProbabilistic(1000);
        long t2 = System.currentTimeMillis();
        MCCPResult resultExact = solver.bruteForceOptimal();
        long t3 = System.currentTimeMillis();

        System.out.println("[" + title + "] " + numNodes + " nodi, " + edges.size() + " archi, "
                + numColors + " colori, s=" + sourceNode + ", t=" + targetNode);
        System.out.println("Costi dei colori: " + Arrays.toString(colorCost));
        System.out.println("VNS-Greedy        -> " + resultGreedy + "  [" + (t1 - t0) + " ms]");
        System.out.println("VNS-Probabilistic -> " + resultProbabilistic + "  [" + (t2 - t1) + " ms]");
        System.out.println("Ottimo esatto     -> " + resultExact + "  [" + (t3 - t2) + " ms, forza bruta]");
    }
}