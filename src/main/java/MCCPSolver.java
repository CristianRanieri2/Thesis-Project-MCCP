import java.util.*;

/**
 * Implementazione basata su Variable Neighborhood Search (VNS) per una
 * variante PESATA del Minimum Color s-t Cut Problem (MCstCP).
 */
public class MCCPSolver {

    // ---------- Rappresentazione del grafo ----------

    public static class Edge {
        final int u, v;
        final Set<Integer> colors;

        public Edge(int u, int v, int... colors) {
            this.u = u;
            this.v = v;
            this.colors = new HashSet<>();
            for (int c : colors) this.colors.add(c);
        }

        public Edge(int u, int v, Collection<Integer> colors) {
            this.u = u;
            this.v = v;
            this.colors = new HashSet<>(colors);
        }
    }

    private final int numNodes;
    private final List<Edge> edges;
    private final int numColors;
    private final double[] colorCost;
    private final int s;
    private final int t;

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

    public int getNumNodes() { return numNodes; }
    public int getNumEdges() { return edges.size(); }
    public int getNumColors() { return numColors; }
    public double[] getColorCost() { return colorCost.clone(); }
    public int getSourceNode() { return s; }
    public int getTargetNode() { return t; }
    public List<Edge> getEdges() { return new ArrayList<>(edges); }

    // ---------- Helper temporale ----------
    private boolean isTimeExpired(long startTime, long maxRunningTimeMillis) {
        return (System.currentTimeMillis() - startTime) >= maxRunningTimeMillis;
    }

    // ---------- Union-Find ----------

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
        for (Edge e : edges) {
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

    // ---------- Selezione greedy con check temporale ----------

    private boolean tryAddBestFeasibleColor(Set<Integer> base, Set<Integer> candidatePool, long startTime, long maxRunningTimeMillis) {
        int bestColor = -1;
        double bestCost = Double.NEGATIVE_INFINITY;

        for (int c : candidatePool) {
            if (isTimeExpired(startTime, maxRunningTimeMillis)) break;
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

    public Set<Integer> generateInitialSolutionGreedy(long startTime, long maxRunningTimeMillis) {
        Set<Integer> bestS = new HashSet<>();
        boolean progressed = true;
        while (progressed && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            Set<Integer> remaining = diff(allColors(), bestS);
            progressed = tryAddBestFeasibleColor(bestS, remaining, startTime, maxRunningTimeMillis);
        }
        return bestS;
    }

    public Set<Integer> newSolutionGreedy(Set<Integer> bestS, long startTime, long maxRunningTimeMillis) {
        Set<Integer> s = new HashSet<>();
        boolean progressed = true;
        while (isFeasible(s) && !diff(bestS, s).isEmpty() && progressed && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            progressed = tryAddBestFeasibleColor(s, diff(bestS, s), startTime, maxRunningTimeMillis);
        }

        progressed = true;
        while (isFeasible(s) && progressed && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            progressed = tryAddBestFeasibleColor(s, diff(bestS, s), startTime, maxRunningTimeMillis);
        }

        return s;
    }

    public void localSearchGreedy(Set<Integer> sPrime, long startTime, long maxRunningTimeMillis) {
        while (isFeasible(sPrime) && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            Set<Integer> complement = diff(allColors(), sPrime);
            if (!tryAddBestFeasibleColor(sPrime, complement, startTime, maxRunningTimeMillis)) break;
        }
    }

    // ---------- Subroutine Probabilistiche con check temporale ----------

    private static final double TEMPERATURE = 1.0;

    private Integer pickProbabilistic(Set<Integer> base, Set<Integer> candidatePool, Random rnd, long startTime, long maxRunningTimeMillis) {
        if (candidatePool.isEmpty() || isTimeExpired(startTime, maxRunningTimeMillis)) return null;

        double maxCost = Double.NEGATIVE_INFINITY;
        for (int c : candidatePool) {
            if (colorCost[c] > maxCost) maxCost = colorCost[c];
        }

        List<Integer> feasibleCandidates = new ArrayList<>();
        List<Double> weights = new ArrayList<>();
        double totalWeight = 0.0;

        for (int c : candidatePool) {
            if (isTimeExpired(startTime, maxRunningTimeMillis)) break;
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

    public Set<Integer> generateInitialSolutionProbabilistic(Random rnd, long startTime, long maxRunningTimeMillis) {
        Set<Integer> bestS = new HashSet<>();
        Set<Integer> remaining = allColors();

        boolean endLoop = false;
        while (!endLoop && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            Integer chosen = pickProbabilistic(bestS, remaining, rnd, startTime, maxRunningTimeMillis);
            if (chosen != null) {
                bestS.add(chosen);
                remaining.remove(chosen);
            } else {
                endLoop = true;
            }
        }
        return bestS;
    }

    public Set<Integer> newSolutionProbabilistic(Set<Integer> bestS, Random rnd, long startTime, long maxRunningTimeMillis) {
        Set<Integer> s = new HashSet<>();

        boolean progressed = true;
        while (isFeasible(s) && !diff(bestS, s).isEmpty() && progressed && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            Integer chosen = pickProbabilistic(s, diff(bestS, s), rnd, startTime, maxRunningTimeMillis);
            if (chosen != null) {
                s.add(chosen);
            } else {
                progressed = false;
            }
        }

        progressed = true;
        while (isFeasible(s) && progressed && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            Integer chosen = pickProbabilistic(s, diff(bestS, s), rnd, startTime, maxRunningTimeMillis);
            if (chosen != null) {
                s.add(chosen);
            } else {
                progressed = false;
            }
        }

        return s;
    }

    public void localSearchProbabilistic(Set<Integer> sPrime, Random rnd, long startTime, long maxRunningTimeMillis) {
        while (isFeasible(sPrime) && !isTimeExpired(startTime, maxRunningTimeMillis)) {
            Set<Integer> complement = diff(allColors(), sPrime);
            Integer chosen = pickProbabilistic(sPrime, complement, rnd, startTime, maxRunningTimeMillis);
            if (chosen != null) {
                sPrime.add(chosen);
            } else {
                break;
            }
        }
    }

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

    public void fix(Set<Integer> sPrime, Random rnd) {
        while (!isFeasible(sPrime)) {
            if (sPrime.isEmpty()) break;
            List<Integer> list = new ArrayList<>(sPrime);
            int toRemove = list.get(rnd.nextInt(list.size()));
            sPrime.remove(toRemove);
        }
    }

    // ---------- Risultati ----------

    public static class MCCPResult {
        public final double cutCost;
        public final int cutSize;
        public final Set<Integer> cutColors;
        public final Set<Integer> keptColors;
        public final long timeToBestMs;
        public final long totalTimeMs;

        public MCCPResult(Set<Integer> cutColors, Set<Integer> keptColors, double cutCost,
                          long timeToBestMs, long totalTimeMs) {
            this.cutColors = Collections.unmodifiableSet(new TreeSet<>(cutColors));
            this.keptColors = Collections.unmodifiableSet(new TreeSet<>(keptColors));
            this.cutSize = cutColors.size();
            this.cutCost = cutCost;
            this.timeToBestMs = timeToBestMs;
            this.totalTimeMs = totalTimeMs;
        }

        @Override
        public String toString() {
            return "costo del taglio = " + cutCost + ", colori nel taglio = " + cutColors
                    + " (" + cutSize + " colori; colori mantenuti: " + keptColors + ")"
                    + "\n [time-to-best = " + timeToBestMs + " ms, tempo totale = " + totalTimeMs + " ms]";
        }
    }

    private MCCPResult buildResult(Set<Integer> bestS, long timeToBestMs, long totalTimeMs) {
        Set<Integer> cutColors = diff(allColors(), bestS);
        double cutCost = weight(cutColors);
        return new MCCPResult(cutColors, bestS, cutCost, timeToBestMs, totalTimeMs);
    }

    // ---------- Algoritmi con controllo temporale ad ogni passo ----------

    public MCCPResult solve(long maxRunningTimeMillis) {
        Random rnd = new Random();
        long startTime = System.currentTimeMillis();

        Set<Integer> bestS = generateInitialSolutionGreedy(startTime, maxRunningTimeMillis);
        int maxNeighborhood = numColors - bestS.size();
        long timeToBestMs = System.currentTimeMillis() - startTime;

        do {
            if (isTimeExpired(startTime, maxRunningTimeMillis)) break;

            Set<Integer> s = newSolutionGreedy(bestS, startTime, maxRunningTimeMillis);

            while (weight(s) > weight(bestS) && !isTimeExpired(startTime, maxRunningTimeMillis)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
                timeToBestMs = System.currentTimeMillis() - startTime;
                s = newSolutionGreedy(bestS, startTime, maxRunningTimeMillis);
            }

            int k = 1;

            while (k < maxNeighborhood && !isTimeExpired(startTime, maxRunningTimeMillis)) {
                Set<Integer> sPrime = new HashSet<>(s);
                shake(sPrime, k, s, rnd);

                if (!isFeasible(sPrime)) {
                    fix(sPrime, rnd);
                }

                localSearchGreedy(sPrime, startTime, maxRunningTimeMillis);

                if (weight(sPrime) > weight(s)) {
                    s = sPrime;
                    k = 1;
                } else {
                    k++;
                }
            }

            if (weight(s) > weight(bestS)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
                timeToBestMs = System.currentTimeMillis() - startTime;
            }

        } while (!isTimeExpired(startTime, maxRunningTimeMillis));

        long totalTimeMs = System.currentTimeMillis() - startTime;
        return buildResult(bestS, timeToBestMs, totalTimeMs);
    }

    public MCCPResult solveProbabilistic(long maxRunningTimeMillis) {
        Random rnd = new Random();
        long startTime = System.currentTimeMillis();

        Set<Integer> bestS = generateInitialSolutionProbabilistic(rnd, startTime, maxRunningTimeMillis);
        int maxNeighborhood = numColors - bestS.size();
        long timeToBestMs = System.currentTimeMillis() - startTime;

        do {
            if (isTimeExpired(startTime, maxRunningTimeMillis)) break;

            Set<Integer> s = newSolutionProbabilistic(bestS, rnd, startTime, maxRunningTimeMillis);

            while (weight(s) > weight(bestS) && !isTimeExpired(startTime, maxRunningTimeMillis)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
                timeToBestMs = System.currentTimeMillis() - startTime;
                s = newSolutionProbabilistic(bestS, rnd, startTime, maxRunningTimeMillis);
            }

            int k = 1;

            while (k < maxNeighborhood && !isTimeExpired(startTime, maxRunningTimeMillis)) {
                Set<Integer> sPrime = new HashSet<>(s);
                shake(sPrime, k, s, rnd);

                if (!isFeasible(sPrime)) {
                    fix(sPrime, rnd);
                }

                localSearchProbabilistic(sPrime, rnd, startTime, maxRunningTimeMillis);

                if (weight(sPrime) > weight(s)) {
                    s = sPrime;
                    k = 1;
                } else {
                    k++;
                }
            }

            if (weight(s) > weight(bestS)) {
                bestS = new HashSet<>(s);
                maxNeighborhood = numColors - bestS.size();
                timeToBestMs = System.currentTimeMillis() - startTime;
            }

        } while (!isTimeExpired(startTime, maxRunningTimeMillis));

        long totalTimeMs = System.currentTimeMillis() - startTime;
        return buildResult(bestS, timeToBestMs, totalTimeMs);
    }

    public boolean verifyCutWithBFS(Set<Integer> cutColors, boolean verbose) {
        List<List<Integer>> adjacency = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) adjacency.add(new ArrayList<>());

        for (Edge e : edges) {
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

        boolean valida = !visited[t];
        if (verbose) {
            System.out.println("Esito della soluzione: " + (valida ? "VALIDA" : "NON VALIDA"));
        }
        return valida;
    }

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

    public enum Density {
        LOW(0.2), MEDIUM(0.5), HIGH(0.8);
        public final double value;
        Density(double value) { this.value = value; }
    }

    private static List<Edge> generateRandomGraph(int numNodes, int numColors, double density, double multiColorProb, Random rnd) {
        List<Edge> edges = new ArrayList<>();
        List<Integer> order = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) order.add(i);
        Collections.shuffle(order, rnd);

        for (int i = 1; i < numNodes; i++) {
            int u = order.get(i);
            int v = order.get(rnd.nextInt(i));
            int color = rnd.nextInt(numColors);
            edges.add(new Edge(u, v, color));
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
                edges.add(new Edge(u, v, color1, color2));
            } else {
                edges.add(new Edge(u, v, color1));
            }
        }
        return edges;
    }

    public static MCCPSolver generateRandomInstance(int numNodes, int numColors, Density density, double multiColorProb) {
        return generateRandomInstance(numNodes, numColors, density, multiColorProb, System.nanoTime());
    }

    public static MCCPSolver generateRandomInstance(int numNodes, int numColors, Density density, double multiColorProb, long seed) {
        if (numNodes < 2) {
            throw new IllegalArgumentException("numNodes deve essere almeno 2 per avere s e t distinti");
        }
        Random rnd = new Random(seed);
        List<Edge> edges = generateRandomGraph(numNodes, numColors, density.value, multiColorProb, rnd);
        double[] colorCost = new double[numColors];
        for (int c = 0; c < numColors; c++) {
            colorCost[c] = 1 + rnd.nextInt(50);
        }
        int sourceNode = rnd.nextInt(numNodes);
        int targetNode;
        do {
            targetNode = rnd.nextInt(numNodes);
        } while (targetNode == sourceNode);

        return new MCCPSolver(numNodes, edges, numColors, colorCost, sourceNode, targetNode);
    }
}