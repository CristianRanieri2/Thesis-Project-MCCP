import java.util.*;

/**
 * Algoritmo ESATTO di tipo Branch-and-Cut basato sul MODELLO A POTENZIALI
 * descritto nella tesi (capitolo "Il Modello Matematico di SRTIP"), nella sua
 * versione corretta per grafi non orientati, rinforzato con gli stessi tagli
 * di cammino lazy usati in {@link MCCPBranchAndCut} (che invece parte da un
 * rilassamento VUOTO fatto solo di tagli generati lazy, senza vincoli di
 * potenziale). Questo file esiste per poter confrontare sperimentalmente le
 * due formulazioni.
 *
 * ============================================================================
 * MODELLO MATEMATICO DI BASE (dalla tesi, con la correzione bidirezionale)
 * ============================================================================
 * Variabili:
 *   z_i in {0,1}   per ogni colore c_i: 1 se il colore viene interdetto (scelto per il taglio)
 *   x_v in [0,1]   per ogni nodo v: potenziale, determina da che lato del taglio sta v
 *
 * Obiettivo:
 *   min  sum_i w_i * z_i
 *
 * Vincoli:
 *   x_s = 0,  x_t = 1
 *   x_v - x_u <= sum_{i : c_i in colori(e)} z_i     per OGNI arco NON ORIENTATO
 *   x_u - x_v <= sum_{i : c_i in colori(e)} z_i     e = {u,v}, ENTRAMBE le direzioni
 *
 * (la doppia direzione e' la correzione discussa: senza di essa il modello
 * sarebbe valido solo per grafi orientati).
 *
 * ============================================================================
 * PERCHE' SERVONO COMUNQUE DEI TAGLI (anche con un modello gia' "esatto")
 * ============================================================================
 * Il sistema di vincoli sopra e' una formulazione ESATTA del problema quando
 * z e' vincolato a {0,1}: per ogni z intero esiste una x ammissibile se e solo
 * se z e' un taglio valido (dimostrazione per assurdo via telescoping lungo un
 * cammino, vedi tesi). Il RILASSAMENTO LINEARE (z continuo in [0,1]), pero',
 * NON e' in generale "tight": quando un colore e' condiviso da piu' archi, il
 * poliedro del rilassamento puo' avere vertici frazionari che non
 * corrispondono a nessun taglio reale, indebolendo il bound e quindi il
 * pruning del branch-and-bound.
 *
 * Per rinforzare il rilassamento si aggiungono, in modo lazy, gli stessi
 * VINCOLI DI CAMMINO gia' usati in {@link MCCPBranchAndCut}:
 *
 *     sum_{c in Colori(P)} z_c >= 1     per ogni cammino s-t P
 *
 * Questi vincoli sono impliciti nel modello a potenziali SOLO in forma
 * indebolita (sommando i vincoli d'arco lungo un cammino si ottiene una
 * somma CON RIPETIZIONE dei colori, non l'unione): il vincolo di cammino con
 * l'unione dei colori (senza ripetizioni) e' invece STRETTAMENTE PIU' FORTE
 * ogni volta che un colore ricorre su piu' archi dello stesso cammino, e va
 * quindi aggiunto esplicitamente come taglio.
 *
 * ============================================================================
 * SEPARAZIONE
 * ============================================================================
 * Identica a {@link MCCPBranchAndCut}: dato un punto (anche frazionario) del
 * rilassamento, si pesano gli archi con la somma degli z correnti e si cerca
 * il cammino minimo s-t con Dijkstra. Se il peso trovato e' < 1, l'unione dei
 * colori lungo quel cammino da' un vincolo violato, che viene aggiunto e la
 * LP viene riottimizzata.
 *
 * ============================================================================
 * UN'INFEASIBILITA' "IN PIU'" RISPETTO AL MODELLO SOLO-TAGLI
 * ============================================================================
 * A differenza di {@link MCCPBranchAndCut} (dove l'unico modo per un nodo di
 * essere infeasibile e' che un taglio noto non abbia piu' colori liberi), qui
 * il rilassamento stesso puo' risultare INFEASIBILE anche PRIMA di generare
 * alcun taglio di cammino: se durante il branching si fissano a 0 tutti i
 * colori di un intero "ponte" che collegherebbe s a t, i soli vincoli
 * d'arco (potenziali) sono gia' in contraddizione (forzano x_s = x_t). Il
 * Simplex deve quindi saper segnalare l'infeasibilita' della LP stessa (non
 * solo la violazione di un taglio), gestita qui restituendo un risultato
 * "infeasibile" invece di lanciare un'eccezione, e trattata come pruning del
 * nodo esattamente come le altre infeasibilita'.
 *
 * ============================================================================
 * NOTA SULLA SCALABILITA'
 * ============================================================================
 * Il rilassamento di base ha O(|V|) variabili ma O(|E|) VINCOLI (uno per ogni
 * arco, in entrambe le direzioni) fin dal nodo radice, a differenza di
 * {@link MCCPBranchAndCut} che parte da ZERO vincoli e li aggiunge solo se
 * violati. Su grafi densi con molti nodi questo rende ogni singola
 * risoluzione della LP molto piu' pesante (il Simplex qui implementato e' un
 * tableau denso scritto da zero, non un solver professionale). Questo file va
 * quindi usato per il confronto su istanze piccole/medie, non per le stesse
 * dimensioni con cui si stressa VNS.
 */
public class MCCPBranchAndCutPotential {

    private final int numNodes;
    private final List<MCCPSolver.Edge> edges;
    private final int numColors;
    private final double[] colorCost;
    private final int s;
    private final int t;

    private final List<Set<Integer>> globalCuts = new ArrayList<>();
    private double bestCost = Double.POSITIVE_INFINITY;
    private int[] bestSolution;
    private long nodesExplored = 0;
    private boolean provenOptimal = true;

    private long searchStartMillis;
    private long maxRunningTimeMillis;
    private long bestFoundAtMs = 0;

    /** Numero di vincoli d'arco (potenziale) presenti nella LP fin dal nodo radice: O(|E|). */
    private int rootEdgeConstraintCount = -1;

    public MCCPBranchAndCutPotential(int numNodes, List<MCCPSolver.Edge> edges, int numColors,
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
    // Dijkstra di separazione (identico, per costruzione, a MCCPBranchAndCut)
    // ========================================================================

    private static final class SeparationResult {
        final double pathWeight;
        final Set<Integer> colorsOnPath;
        SeparationResult(double pathWeight, Set<Integer> colorsOnPath) {
            this.pathWeight = pathWeight;
            this.colorsOnPath = colorsOnPath;
        }
    }

    private SeparationResult dijkstraSeparation(double[] yFull) {
        List<List<int[]>> adjNeighborEdge = new ArrayList<>();
        for (int i = 0; i < numNodes; i++) adjNeighborEdge.add(new ArrayList<>());
        double[] edgeWeight = new double[edges.size()];
        for (int idx = 0; idx < edges.size(); idx++) {
            MCCPSolver.Edge e = edges.get(idx);
            double w = 0.0;
            for (int c : e.colors) w += yFull[c];
            edgeWeight[idx] = w;
            adjNeighborEdge.get(e.u).add(new int[]{e.v, idx});
            adjNeighborEdge.get(e.v).add(new int[]{e.u, idx});
        }

        double[] dist = new double[numNodes];
        Arrays.fill(dist, Double.POSITIVE_INFINITY);
        int[] parentNode = new int[numNodes];
        int[] parentEdge = new int[numNodes];
        Arrays.fill(parentNode, -1);
        Arrays.fill(parentEdge, -1);

        dist[s] = 0.0;
        PriorityQueue<double[]> pq = new PriorityQueue<>(Comparator.comparingDouble(a -> a[0]));
        pq.add(new double[]{0.0, s});
        boolean[] settled = new boolean[numNodes];

        while (!pq.isEmpty()) {
            double[] top = pq.poll();
            int u = (int) top[1];
            if (settled[u]) continue;
            settled[u] = true;
            if (u == t) break;
            for (int[] ne : adjNeighborEdge.get(u)) {
                int v = ne[0];
                int edgeIdx = ne[1];
                double nd = dist[u] + edgeWeight[edgeIdx];
                if (nd < dist[v] - 1e-12) {
                    dist[v] = nd;
                    parentNode[v] = u;
                    parentEdge[v] = edgeIdx;
                    pq.add(new double[]{nd, v});
                }
            }
        }

        if (dist[t] == Double.POSITIVE_INFINITY) {
            return new SeparationResult(Double.POSITIVE_INFINITY, null);
        }

        Set<Integer> colorsOnPath = new HashSet<>();
        int cur = t;
        while (cur != s) {
            int edgeIdx = parentEdge[cur];
            colorsOnPath.addAll(edges.get(edgeIdx).colors);
            cur = parentNode[cur];
        }
        return new SeparationResult(dist[t], colorsOnPath);
    }

    // ========================================================================
    // Simplex Big-M GENERALIZZATO: righe sia <= (con slack) sia >= (con
    // surplus + artificiale). Serve rispetto a MCCPBranchAndCut perche' i
    // vincoli d'arco del modello a potenziali, dopo la sostituzione dei nodi
    // fissati (s,t) e dei colori fissati, possono avere RHS negativo (vanno
    // quindi capovolti in vincoli >=).
    // ========================================================================

    private enum RowType { LE, GE }

    private static final class LPRow {
        final Map<Integer, Double> coeffs;
        final RowType type;
        final double rhs;
        LPRow(Map<Integer, Double> coeffs, RowType type, double rhs) {
            this.coeffs = coeffs;
            this.type = type;
            this.rhs = rhs;
        }
    }

    private static final class LPResult {
        final boolean feasible;
        final double objective;
        final double[] x;
        LPResult(boolean feasible, double objective, double[] x) {
            this.feasible = feasible;
            this.objective = objective;
            this.x = x;
        }
    }

    private static LPResult solveLP(int n, List<LPRow> rows, double[] cost) {
        int m = rows.size();
        if (m == 0) {
            return new LPResult(true, 0.0, new double[n]);
        }

        final double BIG_M = 1.0e7;
        int ncols = n + m + m; // variabili reali + surplus/slack + artificiali
        int rhsCol = ncols;
        int rowsCount = m + 1;

        double[][] tab = new double[rowsCount][ncols + 1];
        int[] basis = new int[m];

        for (int i = 0; i < m; i++) {
            LPRow row = rows.get(i);
            for (Map.Entry<Integer, Double> entry : row.coeffs.entrySet()) {
                tab[i][entry.getKey()] += entry.getValue();
            }
            tab[i][rhsCol] = row.rhs;
            if (row.type == RowType.LE) {
                tab[i][n + i] = 1.0; // slack
                basis[i] = n + i;
            } else {
                tab[i][n + i] = -1.0;      // surplus
                tab[i][n + m + i] = 1.0;   // artificiale
                basis[i] = n + m + i;
            }
        }

        double[] obj = new double[ncols + 1];
        for (int j = 0; j < n; j++) obj[j] = cost[j];
        for (int i = 0; i < m; i++) {
            if (rows.get(i).type == RowType.GE) obj[n + m + i] = BIG_M;
        }
        tab[m] = obj;

        for (int i = 0; i < m; i++) {
            if (basis[i] >= n + m) {
                double factor = tab[m][basis[i]];
                if (factor != 0.0) {
                    for (int k = 0; k <= ncols; k++) tab[m][k] -= factor * tab[i][k];
                }
            }
        }

        final double EPS = 1e-9;
        final int MAX_ITER = 20000;
        int iter = 0;
        while (iter++ < MAX_ITER) {
            int entering = -1;
            for (int j = 0; j < ncols; j++) {
                if (tab[m][j] < -EPS) { entering = j; break; }
            }
            if (entering == -1) break;

            int leaving = -1;
            double bestRatio = Double.POSITIVE_INFINITY;
            for (int i = 0; i < m; i++) {
                if (tab[i][entering] > EPS) {
                    double ratio = tab[i][rhsCol] / tab[i][entering];
                    if (ratio < bestRatio - 1e-12 ||
                            (Math.abs(ratio - bestRatio) < 1e-12 && (leaving == -1 || basis[i] < basis[leaving]))) {
                        bestRatio = ratio;
                        leaving = i;
                    }
                }
            }
            if (leaving == -1) {
                // illimitato: non dovrebbe accadere in questo modello (tutte le
                // variabili hanno un upper bound esplicito o un costo positivo)
                return new LPResult(false, Double.POSITIVE_INFINITY, null);
            }

            double pivotVal = tab[leaving][entering];
            for (int k = 0; k <= ncols; k++) tab[leaving][k] /= pivotVal;
            for (int i = 0; i < rowsCount; i++) {
                if (i == leaving) continue;
                double factor = tab[i][entering];
                if (Math.abs(factor) > EPS) {
                    for (int k = 0; k <= ncols; k++) tab[i][k] -= factor * tab[leaving][k];
                }
            }
            basis[leaving] = entering;
        }

        // se un'artificiale resta basica con valore positivo, la LP e' infeasibile
        for (int i = 0; i < m; i++) {
            if (basis[i] >= n + m && tab[i][rhsCol] > 1e-6) {
                return new LPResult(false, Double.POSITIVE_INFINITY, null);
            }
        }

        double[] x = new double[n];
        for (int i = 0; i < m; i++) {
            if (basis[i] < n) x[basis[i]] = tab[i][rhsCol];
        }
        double objective = 0.0;
        for (int j = 0; j < n; j++) objective += cost[j] * x[j];
        return new LPResult(true, objective, x);
    }

    // ========================================================================
    // Costruzione della LP (modello a potenziali) per un nodo dell'albero
    // ========================================================================

    private static final class BuildResult {
        final List<Integer> freeColors;   // indici ORIGINALI dei colori liberi
        final int[] freeColorLocalIndex;  // colore originale -> indice locale (-1 se non libero)
        final int nVars;
        final List<LPRow> rows;           // mutabile: la separazione ne aggiunge altre
        final double[] costVec;

        BuildResult(List<Integer> freeColors, int[] freeColorLocalIndex, int nVars,
                    List<LPRow> rows, double[] costVec) {
            this.freeColors = freeColors;
            this.freeColorLocalIndex = freeColorLocalIndex;
            this.nVars = nVars;
            this.rows = rows;
            this.costVec = costVec;
        }
    }

    /** @return null se il nodo e' gia' infeasibile per costruzione (un taglio noto non ha piu' colori liberi) */
    private BuildResult buildLP(int[] fixed) {
        List<Integer> freeColors = new ArrayList<>();
        int[] freeColorLocalIndex = new int[numColors];
        Arrays.fill(freeColorLocalIndex, -1);
        for (int c = 0; c < numColors; c++) {
            if (fixed[c] == -1) {
                freeColorLocalIndex[c] = freeColors.size();
                freeColors.add(c);
            }
        }
        int nColorsFree = freeColors.size();

        List<Integer> freeNodes = new ArrayList<>();
        int[] nodeLocalIndex = new int[numNodes];
        Arrays.fill(nodeLocalIndex, -1);
        for (int v = 0; v < numNodes; v++) {
            if (v != s && v != t) {
                nodeLocalIndex[v] = nColorsFree + freeNodes.size();
                freeNodes.add(v);
            }
        }

        int nVars = nColorsFree + freeNodes.size();
        List<LPRow> rows = new ArrayList<>();

        // upper bound: z_c <= 1 per ogni colore libero
        for (int i = 0; i < nColorsFree; i++) {
            Map<Integer, Double> coeffs = new HashMap<>();
            coeffs.put(i, 1.0);
            rows.add(new LPRow(coeffs, RowType.LE, 1.0));
        }
        // upper bound: x_v <= 1 per ogni nodo libero
        for (int v : freeNodes) {
            Map<Integer, Double> coeffs = new HashMap<>();
            coeffs.put(nodeLocalIndex[v], 1.0);
            rows.add(new LPRow(coeffs, RowType.LE, 1.0));
        }

        // vincoli di taglio sugli archi, ENTRAMBE le direzioni (modello a potenziali corretto)
        for (MCCPSolver.Edge e : edges) {
            List<Integer> freeInEdge = new ArrayList<>();
            int fixed1InEdge = 0;
            for (int c : e.colors) {
                if (fixed[c] == -1) freeInEdge.add(freeColorLocalIndex[c]);
                else if (fixed[c] == 1) fixed1InEdge++;
            }

            int[][] directions = { {e.u, e.v}, {e.v, e.u} };
            for (int[] dir : directions) {
                int a = dir[0], b = dir[1]; // vincolo: x_b - x_a <= sum_free z_c(e) + fixed1InEdge
                Map<Integer, Double> coeffs = new HashMap<>();
                double constB = 0.0, constA = 0.0;
                if (b == s) { constB = 0.0; }
                else if (b == t) { constB = 1.0; }
                else { coeffs.merge(nodeLocalIndex[b], 1.0, Double::sum); }

                if (a == s) { constA = 0.0; }
                else if (a == t) { constA = 1.0; }
                else { coeffs.merge(nodeLocalIndex[a], -1.0, Double::sum); }

                for (int localColor : freeInEdge) {
                    coeffs.merge(localColor, -1.0, Double::sum);
                }

                double rhs = fixed1InEdge - constB + constA;
                if (rhs >= 0) {
                    rows.add(new LPRow(coeffs, RowType.LE, rhs));
                } else {
                    Map<Integer, Double> flipped = new HashMap<>();
                    for (Map.Entry<Integer, Double> en : coeffs.entrySet()) flipped.put(en.getKey(), -en.getValue());
                    rows.add(new LPRow(flipped, RowType.GE, -rhs));
                }
            }
        }

        // tagli di cammino lazy noti finora (globali), ristretti FRESCHI ai colori liberi di questo nodo
        for (Set<Integer> cut : globalCuts) {
            boolean satisfiedByFixed1 = false;
            for (int c : cut) if (fixed[c] == 1) { satisfiedByFixed1 = true; break; }
            if (satisfiedByFixed1) continue;

            Map<Integer, Double> coeffs = new HashMap<>();
            for (int c : cut) if (fixed[c] == -1) coeffs.put(freeColorLocalIndex[c], 1.0);
            if (coeffs.isEmpty()) return null; // nodo infeasibile
            rows.add(new LPRow(coeffs, RowType.GE, 1.0));
        }

        double[] costVec = new double[nVars];
        for (int i = 0; i < nColorsFree; i++) costVec[i] = colorCost[freeColors.get(i)];
        // i costi delle variabili di potenziale x_v restano 0 (non compaiono nell'obiettivo)

        return new BuildResult(freeColors, freeColorLocalIndex, nVars, rows, costVec);
    }

    // ========================================================================
    // Warm start: identico a MCCPBranchAndCut (euristica greedy sui cammini scoperti)
    // ========================================================================

    private void warmStart() {
        double[] y = new double[numColors];
        while (true) {
            SeparationResult sep = dijkstraSeparation(y);
            if (sep.colorsOnPath == null || sep.pathWeight >= 1.0 - 1e-9) break;
            int cheapest = -1;
            double cheapestCost = Double.POSITIVE_INFINITY;
            for (int c : sep.colorsOnPath) {
                if (colorCost[c] < cheapestCost) {
                    cheapestCost = colorCost[c];
                    cheapest = c;
                }
            }
            y[cheapest] = 1.0;
        }
        double cost = 0.0;
        int[] sol = new int[numColors];
        for (int c = 0; c < numColors; c++) {
            sol[c] = y[c] >= 0.5 ? 1 : 0;
            if (sol[c] == 1) cost += colorCost[c];
        }
        bestCost = cost;
        bestSolution = sol;
        bestFoundAtMs = System.currentTimeMillis() - searchStartMillis;
    }

    // ========================================================================
    // Branch-and-Cut
    // ========================================================================

    public MCCPSolver.MCCPResult solveExact(long maxRunningTimeMillis) {
        this.searchStartMillis = System.currentTimeMillis();
        this.maxRunningTimeMillis = maxRunningTimeMillis;
        this.nodesExplored = 0;
        this.provenOptimal = true;
        this.globalCuts.clear();
        this.rootEdgeConstraintCount = -1;

        warmStart();

        int[] fixed = new int[numColors];
        Arrays.fill(fixed, -1);
        branch(fixed);

        long totalTimeMs = System.currentTimeMillis() - searchStartMillis;

        Set<Integer> keptColors = new HashSet<>();
        Set<Integer> cutColors = new HashSet<>();
        for (int c = 0; c < numColors; c++) {
            if (bestSolution[c] == 1) cutColors.add(c);
            else keptColors.add(c);
        }
        return new MCCPSolver.MCCPResult(cutColors, keptColors, bestCost, bestFoundAtMs, totalTimeMs);
    }

    private void branch(int[] fixed) {
        nodesExplored++;

        if (!provenOptimal) return;
        if (maxRunningTimeMillis > 0 && System.currentTimeMillis() - searchStartMillis > maxRunningTimeMillis) {
            provenOptimal = false;
            return;
        }

        BuildResult built = buildLP(fixed);
        if (built == null) return; // nodo infeasibile: pruning

        if (rootEdgeConstraintCount == -1) {
            // registra, solo la prima volta (nodo radice), quante righe ha la LP di base
            rootEdgeConstraintCount = built.rows.size();
        }

        List<LPRow> rows = built.rows; // mutabile: la separazione ne aggiunge altre

        final int MAX_SEPARATION_ITER = 200;
        for (int iter = 0; iter < MAX_SEPARATION_ITER; iter++) {
            if (maxRunningTimeMillis > 0 && System.currentTimeMillis() - searchStartMillis > maxRunningTimeMillis) {
                provenOptimal = false;
                return;
            }

            LPResult lp = built.nVars == 0 ? new LPResult(true, 0.0, new double[0])
                    : solveLP(built.nVars, rows, built.costVec);
            if (!lp.feasible) return; // nodo infeasibile: pruning (puo' accadere gia' dai soli vincoli d'arco)

            double fixed1Cost = 0.0;
            for (int c = 0; c < numColors; c++) if (fixed[c] == 1) fixed1Cost += colorCost[c];
            double totalBound = fixed1Cost + lp.objective;

            if (totalBound >= bestCost - 1e-7) return; // pruning

            double[] yFull = new double[numColors];
            for (int c = 0; c < numColors; c++) if (fixed[c] == 1) yFull[c] = 1.0;
            for (int i = 0; i < built.freeColors.size(); i++) yFull[built.freeColors.get(i)] = lp.x[i];

            SeparationResult sep = dijkstraSeparation(yFull);

            if (sep.colorsOnPath == null || sep.pathWeight >= 1.0 - 1e-7) {
                // nessun cammino violato: controlla integralita' delle SOLE variabili z (le x sono ausiliarie)
                boolean integral = true;
                for (int i = 0; i < built.freeColors.size(); i++) {
                    double v = lp.x[i];
                    if (v > 1e-7 && v < 1 - 1e-7) { integral = false; break; }
                }

                if (integral) {
                    if (totalBound < bestCost - 1e-9) {
                        bestCost = totalBound;
                        int[] sol = fixed.clone();
                        for (int i = 0; i < built.freeColors.size(); i++) {
                            sol[built.freeColors.get(i)] = lp.x[i] > 0.5 ? 1 : 0;
                        }
                        bestSolution = sol;
                        bestFoundAtMs = System.currentTimeMillis() - searchStartMillis;
                    }
                    return;
                }

                int branchLocalIdx = -1;
                double branchValue = -1.0;
                for (int i = 0; i < built.freeColors.size(); i++) {
                    double v = lp.x[i];
                    if (v > 1e-7 && v < 1 - 1e-7 && v > branchValue) {
                        branchValue = v;
                        branchLocalIdx = i;
                    }
                }
                int branchColor = built.freeColors.get(branchLocalIdx);

                fixed[branchColor] = 1;
                branch(fixed);
                fixed[branchColor] = 0;
                branch(fixed);
                fixed[branchColor] = -1;
                return;

            } else {
                // vincolo violato: il taglio VALIDO GLOBALMENTE usa l'insieme COMPLETO
                // di colori del cammino (vedi nota nel bugfix di MCCPBranchAndCut: un
                // colore fissato a 0 qui potrebbe essere libero in un altro ramo).
                if (!globalCuts.contains(sep.colorsOnPath)) {
                    globalCuts.add(new HashSet<>(sep.colorsOnPath));
                }

                Map<Integer, Double> coeffs = new HashMap<>();
                for (int c : sep.colorsOnPath) {
                    if (fixed[c] == -1) coeffs.put(built.freeColorLocalIndex[c], 1.0);
                }
                if (coeffs.isEmpty()) return; // nodo infeasibile: pruning
                rows.add(new LPRow(coeffs, RowType.GE, 1.0));
                // ripete il ciclo: risolve di nuovo la LP con il nuovo taglio
            }
        }
    }

    // ========================================================================
    // Esempio di utilizzo: confronto a tre vie su istanze piccole/medie
    // ========================================================================

    public static void main(String[] args) {
        System.out.println("=== Confronto a tre vie: B&C-Potenziali vs B&C-base (solo tagli) vs VNS-Greedy ===");
        for (int trial = 0; trial < 5; trial++) {
            int numNodes = 25 + trial * 4;
            int numColors = 6 + trial;
            MCCPSolver instance = MCCPSolver.generateRandomInstance(numNodes, numColors, MCCPSolver.Density.MEDIUM, 0.10);

            MCCPBranchAndCutPotential bcp = new MCCPBranchAndCutPotential(instance.getNumNodes(), instance.getEdges(),
                    instance.getNumColors(), instance.getColorCost(), instance.getSourceNode(), instance.getTargetNode());
            long t0 = System.currentTimeMillis();
            MCCPSolver.MCCPResult bcpResult = bcp.solveExact(10000);
            long t1 = System.currentTimeMillis();

            MCCPBranchAndCut bc = new MCCPBranchAndCut(instance.getNumNodes(), instance.getEdges(),
                    instance.getNumColors(), instance.getColorCost(), instance.getSourceNode(), instance.getTargetNode());
            long t2 = System.currentTimeMillis();
            MCCPSolver.MCCPResult bcResult = bc.solveExact(10000);
            long t3 = System.currentTimeMillis();

            MCCPSolver.MCCPResult vnsResult = instance.solve(1000);

            System.out.println("\n[trial " + trial + "] " + numNodes + " nodi, " + numColors
                    + " colori, s=" + instance.getSourceNode() + ", t=" + instance.getTargetNode());

            System.out.println("B&C-Potenziali -> " + bcpResult.toString()
                    + "  [nodi albero=" + bcp.getNodesExplored() + ", ottimo certificato=" + bcp.isProvenOptimal()
                    + ", vincoli d'arco alla radice=" + bcp.getRootEdgeConstraintCount()
                    + ", tempo wall-clock=" + (t1 - t0) + " ms]");
            instance.verifyCutWithBFS(bcpResult.cutColors, false);

            System.out.println("B&C-base       -> " + bcResult.toString()
                    + "  [nodi albero=" + bc.getNodesExplored() + ", ottimo certificato=" + bc.isProvenOptimal()
                    + ", tempo wall-clock=" + (t3 - t2) + " ms]");
            instance.verifyCutWithBFS(bcResult.cutColors, false);

            System.out.println("VNS-Greedy     -> " + vnsResult.toString());
            instance.verifyCutWithBFS(vnsResult.cutColors, false);

            boolean agree = Math.abs(bcpResult.cutCost - bcResult.cutCost) < 1e-6;
            System.out.println("I due Branch-and-Cut concordano sull'ottimo? " + agree);
        }
    }
}