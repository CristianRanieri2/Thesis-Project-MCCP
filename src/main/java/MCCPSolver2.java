import java.util.*;

/**
 * Variante di {@link MCCPSolver} che implementa SOLO l'algoritmo
 * VNS-Probabilistic (nessuna versione greedy), con un criterio di
 * terminazione DOPPIO:
 *
 *   1) TEMPORALE: come in MCCPSolver, l'esecuzione non supera mai
 *      maxRunningTimeMillis.
 *   2) STAGNAZIONE: l'esecuzione si ferma anche prima del tempo massimo se
 *      per maxStagnantIterations iterazioni CONSECUTIVE la soluzione
 *      migliore (BestS) non migliora.
 *
 * L'algoritmo termina al verificarsi del PRIMO dei due eventi. Questo evita
 * di continuare a far girare la ricerca quando e' evidentemente bloccata su
 * un ottimo locale (o sull'ottimo globale) ben prima della scadenza del
 * tempo, risparmiando tempo di calcolo nei test sperimentali senza
 * sacrificare la qualita' della soluzione.
 *
 * Cosa si intende per "iterazione": un giro completo del ciclo principale
 * dell'Algoritmo 1 (righe 3-28 del paper) -- generazione di una nuova
 * soluzione candidata (New-Solution) seguita dal nucleo VNS (Shake + Fix +
 * Local-Search per tutti i valori di k). Il contatore di stagnazione si
 * azzera ogni volta che, in una qualunque fase di quell'iterazione, BestS
 * viene effettivamente migliorata; altrimenti si incrementa di 1.
 */
public class MCCPSolver2 {

    // ---------- Rappresentazione del grafo (riusa i tipi di MCCPSolver) ----------

    private final int numNodes;
    private final List<MCCPSolver.Edge> edges;
    private final int numColors;
    private final double[] colorCost;
    private final int s;
    private final int t;

    // diagnostica sull'ultima esecuzione di solveProbabilistic
    private boolean lastRunStoppedByStagnation = false;
    private int lastRunStagnantIterationsAtStop = 0;
    private long lastRunTotalIterations = 0;

    public MCCPSolver2(int numNodes, List<MCCPSolver.Edge> edges, int numColors,
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

    // ---------- getter ----------

    public int getNumNodes() { return numNodes; }
    public int getNumEdges() { return edges.size(); }
    public int getNumColors() { return numColors; }
    public double[] getColorCost() { return colorCost.clone(); }
    public int getSourceNode() { return s; }
    public int getTargetNode() { return t; }
    public List<MCCPSolver.Edge> getEdges() { return new ArrayList<>(edges); }

    /** True se l'ultima chiamata a solveProbabilistic si e' fermata per stagnazione (non per il tempo). */
    public boolean isLastRunStoppedByStagnation() { return lastRunStoppedByStagnation; }
    /** Quante iterazioni consecutive senza miglioramento sono state contate quando si e' fermata l'ultima esecuzione. */
    public int getLastRunStagnantIterationsAtStop() { return lastRunStagnantIterationsAtStop; }
    /** Numero totale di iterazioni del ciclo principale eseguite nell'ultima chiamata a solveProbabilistic. */
    public long getLastRunTotalIterations() { return lastRunTotalIterations; }

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

    public boolean connected(int a, int b, Set<Integer> colorSet) {
        UnionFind uf = new UnionFind(numNodes);
        for (MCCPSolver.Edge e : edges) {
            if (colorSet.containsAll(e.colors)) {
                uf.union(e.u, e.v);
            }
        }
        return uf.find(a) == uf.find(b);
    }

    public boolean isFeasible(Set<Integer> colorSet) {
        return !connected(s, t, colorSet);
    }

    public double weight(Set<Integer> colorSet) {
        double total = 0.0;
        for (int c : colorSet) total += colorCost[c];
        return total;
    }

    // ---------- Funzione di Boltzmann (versione probabilistica) ----------

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

    // ---------- 2.1 Generate-Initial-Solution (versione probabilistica) ----------

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

    // ---------- 2.2 New-Solution (versione probabilistica) ----------

    public Set<Integer> newSolutionProbabilistic(Set<Integer> bestS, Random rnd) {
        Set<Integer> sSol = new HashSet<>();

        boolean progressed = true;
        while (isFeasible(sSol) && !diff(bestS, sSol).isEmpty() && progressed) {
            Integer chosen = pickProbabilistic(sSol, diff(bestS, sSol), rnd);
            if (chosen != null) {
                sSol.add(chosen);
            } else {
                progressed = false;
            }
        }

        progressed = true;
        while (isFeasible(sSol) && progressed) {
            Integer chosen = pickProbabilistic(sSol, diff(bestS, sSol), rnd);
            if (chosen != null) {
                sSol.add(chosen);
            } else {
                progressed = false;
            }
        }

        return sSol;
    }

    // ---------- 2.3 Shake (comune) ----------

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

    // ---------- 2.4 Fix (comune) ----------

    public void fix(Set<Integer> sPrime, Random rnd) {
        while (!isFeasible(sPrime)) {
            if (sPrime.isEmpty()) break;
            List<Integer> list = new ArrayList<>(sPrime);
            int toRemove = list.get(rnd.nextInt(list.size()));
            sPrime.remove(toRemove);
        }
    }

    // ---------- 2.5 Local-Search (versione probabilistica) ----------

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

    private MCCPSolver.MCCPResult buildResult(Set<Integer> bestS, long timeToBestMs, long totalTimeMs) {
        Set<Integer> cutColors = diff(allColors(), bestS);
        double cutCost = weight(cutColors);
        return new MCCPSolver.MCCPResult(cutColors, bestS, cutCost, timeToBestMs, totalTimeMs);
    }

    // ========================================================================
    // Algoritmo 1 (VNS-Probabilistic) con doppio criterio di terminazione
    // ========================================================================

    /**
     * Esegue VNS-Probabilistic fino al PRIMO fra questi due eventi:
     *  - il tempo trascorso supera maxRunningTimeMillis;
     *  - per maxStagnantIterations iterazioni consecutive del ciclo
     *    principale, BestS non migliora.
     *
     * @param maxRunningTimeMillis   tempo massimo di esecuzione, in millisecondi
     * @param maxStagnantIterations  numero massimo di iterazioni consecutive
     *                               senza miglioramento prima di fermarsi
     *                               (deve essere > 0)
     */
    public MCCPSolver.MCCPResult solveProbabilistic(long maxRunningTimeMillis, int maxStagnantIterations) {
        if (maxStagnantIterations <= 0) {
            throw new IllegalArgumentException("maxStagnantIterations deve essere un intero positivo");
        }

        Random rnd = new Random();
        long startTime = System.currentTimeMillis();

        Set<Integer> bestS = generateInitialSolutionProbabilistic(rnd);
        int maxNeighborhood = numColors - bestS.size();
        long timeToBestMs = System.currentTimeMillis() - startTime;

        int stagnantIterations = 0;
        long totalIterations = 0;
        boolean stoppedByStagnation = false;

        do {
            totalIterations++;
            boolean improvedThisIteration = false;

            Set<Integer> sSol = newSolutionProbabilistic(bestS, rnd);

            while (weight(sSol) > weight(bestS)) {
                bestS = new HashSet<>(sSol);
                maxNeighborhood = numColors - bestS.size();
                timeToBestMs = System.currentTimeMillis() - startTime;
                improvedThisIteration = true;
                sSol = newSolutionProbabilistic(bestS, rnd);
            }

            int k = 1;
            while (k < maxNeighborhood) {
                Set<Integer> sPrime = new HashSet<>(sSol);
                shake(sPrime, k, sSol, rnd);

                if (!isFeasible(sPrime)) {
                    fix(sPrime, rnd);
                }

                localSearchProbabilistic(sPrime, rnd);

                if (weight(sPrime) > weight(sSol)) {
                    sSol = sPrime;
                    k = 1;
                } else {
                    k++;
                }

                if (System.currentTimeMillis() - startTime > maxRunningTimeMillis) break;
            }

            if (weight(sSol) > weight(bestS)) {
                bestS = new HashSet<>(sSol);
                maxNeighborhood = numColors - bestS.size();
                timeToBestMs = System.currentTimeMillis() - startTime;
                improvedThisIteration = true;
            }

            if (improvedThisIteration) {
                stagnantIterations = 0;
            } else {
                stagnantIterations++;
            }

            if (stagnantIterations >= maxStagnantIterations) {
                stoppedByStagnation = true;
                break;
            }

        } while (System.currentTimeMillis() - startTime <= maxRunningTimeMillis);

        this.lastRunStoppedByStagnation = stoppedByStagnation;
        this.lastRunStagnantIterationsAtStop = stagnantIterations;
        this.lastRunTotalIterations = totalIterations;

        long totalTimeMs = System.currentTimeMillis() - startTime;
        return buildResult(bestS, timeToBestMs, totalTimeMs);
    }

    // ---------- Verificatore indipendente della soluzione (rimozione archi + BFS) ----------

    public boolean verifyCutWithBFS(Set<Integer> cutColors, boolean verbose) {
        List<List<Integer>> adjacency = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) adjacency.add(new ArrayList<>());

        int removedEdges = 0;
        for (MCCPSolver.Edge e : edges) {
            if (Collections.disjoint(e.colors, cutColors)) {
                adjacency.get(e.u).add(e.v);
                adjacency.get(e.v).add(e.u);
            } else {
                removedEdges++;
            }
        }

        boolean[] visited = new boolean[numNodes];
        Deque<Integer> queue = new ArrayDeque<>();
        List<Integer> visitOrder = new ArrayList<>();

        visited[s] = true;
        queue.add(s);
        while (!queue.isEmpty()) {
            int u = queue.poll();
            visitOrder.add(u);
            for (int v : adjacency.get(u)) {
                if (!visited[v]) {
                    visited[v] = true;
                    queue.add(v);
                }
            }
        }

        boolean tRaggiunto = visited[t];
        boolean valida = !tRaggiunto;

        if (verbose) {
            if (numNodes <= 60) {
                System.out.println("  Nodi raggiunti dalla BFS a partire da s=" + s + ": "
                        + visitOrder.size() + "/" + numNodes + " -> " + visitOrder);
            } else {
                System.out.println("  Nodi raggiunti dalla BFS a partire da s=" + s + ": "
                        + visitOrder.size() + "/" + numNodes);
            }
            System.out.println("Esito della soluzione: " + (valida ? "VALIDA" : "NON VALIDA"));
        }

        return valida;
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

    // ---------- generazione di istanze casuali (identica a MCCPSolver) ----------

    private static List<MCCPSolver.Edge> generateRandomGraph(int numNodes, int numColors, double density,
                                                             double multiColorProb, Random rnd) {
        List<MCCPSolver.Edge> edges = new ArrayList<>();

        List<Integer> order = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) order.add(i);
        Collections.shuffle(order, rnd);

        for (int i = 1; i < numNodes; i++) {
            int u = order.get(i);
            int v = order.get(rnd.nextInt(i));
            int color = rnd.nextInt(numColors);
            edges.add(new MCCPSolver.Edge(u, v, color));
        }

        long targetTotalEdges = Math.round(density * numNodes * (numNodes - 1) / 2.0);
        int numExtraEdges = (int) Math.max(0, targetTotalEdges - (numNodes - 1));

        for (int i = 0; i < numExtraEdges; i++) {
            int u = rnd.nextInt(numNodes);
            int v = rnd.nextInt(numNodes);
            if (u == v) continue;
            int color1 = rnd.nextInt(numColors);
            if (rnd.nextDouble() < multiColorProb && numColors > 1) {
                int color2;
                do { color2 = rnd.nextInt(numColors); } while (color2 == color1);
                edges.add(new MCCPSolver.Edge(u, v, color1, color2));
            } else {
                edges.add(new MCCPSolver.Edge(u, v, color1));
            }
        }
        return edges;
    }

    public static MCCPSolver2 generateRandomInstance(int numNodes, int numColors,
                                                     MCCPSolver.Density density, double multiColorProb) {
        return generateRandomInstance(numNodes, numColors, density, multiColorProb, System.nanoTime());
    }

    public static MCCPSolver2 generateRandomInstance(int numNodes, int numColors, MCCPSolver.Density density,
                                                     double multiColorProb, long seed) {
        if (numNodes < 2) {
            throw new IllegalArgumentException("numNodes deve essere almeno 2 per avere s e t distinti");
        }

        Random rnd = new Random(seed);

        List<MCCPSolver.Edge> edges = generateRandomGraph(numNodes, numColors, density.value, multiColorProb, rnd);

        double[] colorCost = new double[numColors];
        for (int c = 0; c < numColors; c++) {
            colorCost[c] = 1 + rnd.nextInt(50);
        }

        int sourceNode = rnd.nextInt(numNodes);
        int targetNode;
        do {
            targetNode = rnd.nextInt(numNodes);
        } while (targetNode == sourceNode);

        System.out.println("(istanza generata con seed=" + seed + ", densita'=" + density
                + " (" + density.value + "), per riprodurla esattamente passare questo seed "
                + "a generateRandomInstance(..., seed))");

        return new MCCPSolver2(numNodes, edges, numColors, colorCost, sourceNode, targetNode);
    }

    // ---------- esempio di utilizzo ----------

    public static void main(String[] args) {
        MCCPSolver2 solver = generateRandomInstance(200, 40, MCCPSolver.Density.MEDIUM, 0.10);

        // tempo massimo generoso (60s) per far emergere chiaramente l'effetto
        // del criterio di stagnazione: se la ricerca si blocca su un ottimo
        // (locale o globale) prima della scadenza, si fermera' prima grazie a
        // maxStagnantIterations, senza dover aspettare i 60s pieni.
        long maxTimeMillis = 60_000;
        int maxStagnantIterations = 300;

        MCCPSolver.MCCPResult result = solver.solveProbabilistic(maxTimeMillis, maxStagnantIterations);

        System.out.println("[MCCPSolver2 - VNS-Probabilistic con doppio criterio di terminazione]");
        System.out.println(result);
        System.out.println("Iterazioni totali eseguite: " + solver.getLastRunTotalIterations());
        System.out.println("Fermato per stagnazione? " + solver.isLastRunStoppedByStagnation()
                + " (soglia=" + maxStagnantIterations + ", raggiunte=" + solver.getLastRunStagnantIterationsAtStop() + ")");
        solver.verifyCutWithBFS(result.cutColors, true);
    }
}