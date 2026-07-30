import java.util.*;

/**
 * Algoritmo ESATTO di tipo Branch-and-Cut per la variante pesata del
 * Minimum Color s-t Cut Problem (MCstCP) risolta euristicamente da
 * {@link MCCPSolver} (VNS-Greedy / VNS-Probabilistic).
 *
 * ============================================================================
 * MODELLO MATEMATICO
 * ============================================================================
 * Variabili: y_c in {0,1} per ogni colore c = 1 se e solo se il colore c
 * viene scelto per il taglio (rimosso).
 *
 * Obiettivo:      minimizzare  sum_c costo(c) * y_c
 *
 * Vincoli:        per OGNI cammino P da s a t nel grafo:
 *
 *                     sum_{c in Colori(P)} y_c  >=  1
 *
 *                 dove Colori(P) e' l'unione (senza ripetizioni) dei colori
 *                 di tutti gli archi del cammino P. Il vincolo dice: "lungo
 *                 ogni cammino che collega s a t deve essere rimosso almeno
 *                 un colore", cioe' e' l'esatta traduzione della condizione
 *                 "dopo aver rimosso gli archi dei colori scelti, s e t non
 *                 sono piu' connessi".
 *
 * Questi sono i classici vincoli di tipo "cutset"/"path" usati nei modelli
 * esatti per problemi di taglio (analoghi ai vincoli subtour-elimination
 * del TSP o ai vincoli di connettivita' nei problemi di rete): sono in
 * numero ESPONENZIALE (uno per ogni cammino s-t), quindi non si possono
 * elencare tutti a priori. Si usa percio' la tecnica standard della
 * GENERAZIONE LAZY DI VINCOLI (row generation / cutting planes): si parte
 * senza vincoli di cammino, si risolve il rilassamento, e ogni volta che la
 * soluzione corrente lascia ancora s e t connessi si individua il cammino
 * "colpevole" e si aggiunge il relativo vincolo (taglio), per poi
 * riottimizzare. Questo e' esattamente cio' che rende l'algoritmo un
 * Branch-and-CUT (branch-and-bound + generazione di tagli specifici sul
 * modello) e non un semplice branch-and-bound.
 *
 * ============================================================================
 * SEPARAZIONE (come si trova un vincolo violato)
 * ============================================================================
 * Dato un vettore y (anche frazionario, cioe' preso dal rilassamento
 * lineare), si costruisce un peso per ogni arco e:
 *
 *     peso(e) = somma di y_c per ogni colore c assegnato all'arco e
 *
 * e si calcola il cammino minimo da s a t con Dijkstra (i pesi sono sempre
 * >= 0). Se il cammino minimo trovato ha peso < 1, allora l'insieme dei
 * colori usati lungo quel cammino (senza ripetizioni) soddisfa
 * sum_{c in Colori(P)} y_c <= peso(cammino) < 1, cioe' il vincolo
 * corrispondente a quel cammino e' VIOLATO dalla y corrente: lo si aggiunge
 * al pool di vincoli e si riottimizza. Se il cammino minimo ha peso >= 1,
 * nessun vincolo di questo tipo e' violato per la y corrente (fino a dove
 * la separazione riesce a vedere).
 *
 * Caso particolare importante: se y e' un vettore 0/1 (intero), un arco ha
 * peso 0 se e solo se NESSUNO dei suoi colori e' stato scelto per il
 * taglio, cioe' se e solo se l'arco sopravvive alla rimozione. Quindi per y
 * interi, "il cammino minimo pesato ha peso 0" equivale ESATTAMENTE a "s e
 * t sono ancora connessi nel sottografo che sopravvive al taglio y" (e se
 * il minimo e' > 0 vuol dire che ogni cammino usa almeno un colore
 * rimosso). La stessa routine di separazione via Dijkstra e' quindi, per
 * y interi, un controllo di fattibilita' ESATTO (equivalente alla BFS di
 * {@link MCCPSolver#verifyCutWithBFS}), non solo un'euristica.
 *
 * ============================================================================
 * RILASSAMENTO LINEARE E BOUNDING
 * ============================================================================
 * Ad ogni nodo dell'albero di branch-and-bound alcuni colori sono FISSATI
 * (a 0 = "mai rimosso in questo ramo", o a 1 = "sempre rimosso in questo
 * ramo") e gli altri sono LIBERI. Si risolve il rilassamento lineare
 *
 *     minimizzare  sum_{c libero} costo(c) * y_c  [+ costo dei colori fissati a 1]
 *     soggetto a   sum_{c in taglio_i, c libero} y_c >= 1   per ogni taglio noto
 *                  y_c >= 0
 *
 * con un Simplex Big-M scritto ad hoc (nessuna libreria esterna). Il valore
 * ottimo di questo rilassamento e' un lower bound valido per il nodo (dato
 * che si tratta di un sottoinsieme dei vincoli del problema vero, il suo
 * ottimo e' sempre <= all'ottimo vero, quindi utilizzabile per il pruning).
 * Se il bound non batte la miglior soluzione intera trovata finora,
 * il nodo viene tagliato (pruning). Se la soluzione del rilassamento e'
 * intera e la separazione non trova piu' cammini violati, e' la soluzione
 * ottima per quel nodo. Altrimenti si fa BRANCHING su un colore frazionario
 * (si sceglie quello con valore piu' vicino a 1, creando due figli:
 * fissato a 1 e fissato a 0).
 *
 * ============================================================================
 * WARM START
 * ============================================================================
 * Prima di iniziare il branch-and-cut si calcola una soluzione iniziale
 * ammissibile con un semplice euristico greedy (varianti dell'idea usata da
 * VNS-Greedy in {@link MCCPSolver}): finche' esiste un cammino s-t
 * "scoperto" (nessun suo colore ancora scelto), si sceglie il colore piu'
 * economico fra quelli del cammino e lo si aggiunge al taglio. Questo da'
 * un limite superiore iniziale che rende il pruning efficace fin dal primo
 * nodo.
 */
public class MCCPBranchAndCut {

    private final int numNodes;
    private final List<MCCPSolver.Edge> edges;
    private final int numColors;
    private final double[] colorCost;
    private final int s;
    private final int t;

    // stato della ricerca
    private final List<Set<Integer>> globalCuts = new ArrayList<>(); // vincoli di cammino scoperti finora (validi globalmente)
    private double bestCost = Double.POSITIVE_INFINITY;
    private int[] bestSolution; // 0/1 per ciascun colore (1 = rimosso/nel taglio)
    private long nodesExplored = 0;
    private boolean provenOptimal = true; // diventa false se si esce per timeout prima di aver esplorato tutto

    private long searchStartMillis;
    private long maxRunningTimeMillis;

    public MCCPBranchAndCut(int numNodes, List<MCCPSolver.Edge> edges, int numColors,
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

    /** Numero di nodi dell'albero di branch-and-bound effettivamente esplorati nell'ultima chiamata a solveExact. */
    public long getNodesExplored() { return nodesExplored; }
    /** True se l'ultima chiamata a solveExact ha esplorato l'intero albero (ottimalita' certificata). */
    public boolean isProvenOptimal() { return provenOptimal; }

    // ========================================================================
    // Dijkstra "pesato sui colori": usato sia per la separazione (su y
    // frazionario) sia, quando y e' 0/1, come controllo di fattibilita' ESATTO.
    // ========================================================================

    private static final class SeparationResult {
        final double pathWeight;          // Double.POSITIVE_INFINITY se t irraggiungibile
        final Set<Integer> colorsOnPath;  // unione dei colori lungo il cammino trovato (null se irraggiungibile)

        SeparationResult(double pathWeight, Set<Integer> colorsOnPath) {
            this.pathWeight = pathWeight;
            this.colorsOnPath = colorsOnPath;
        }
    }

    private SeparationResult dijkstraSeparation(double[] yFull) {
        // adiacenza con peso dell'arco = somma di yFull[c] per c nei colori dell'arco
        // rappresentazione: per ogni nodo, lista di (vicino, indiceArco)
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
    // Simplex Big-M (metodo delle penalita') per il rilassamento lineare:
    //   minimizzare  sum_j cost[j] * x_j
    //   soggetto a   sum_{j in cut} x_j >= 1   per ogni cut in "cuts"
    //                x_j >= 0
    // Ogni "cut" e' un insieme di INDICI LOCALI (0..n-1) fra le variabili
    // libere di questo nodo (gia' filtrati/rimappati dal chiamante).
    // ========================================================================

    private static final class LPResult {
        final double objective;
        final double[] x;
        LPResult(double objective, double[] x) { this.objective = objective; this.x = x; }
    }

    private static LPResult solveCoveringLP(int n, List<Set<Integer>> cuts, double[] cost) {
        int m = cuts.size();
        if (m == 0) {
            return new LPResult(0.0, new double[n]);
        }

        final double BIG_M = 1.0e7;
        int ncols = n + m + m; // variabili reali + surplus + artificiali
        int rhsCol = ncols;
        int rows = m + 1; // + riga obiettivo

        double[][] tab = new double[rows][ncols + 1];

        for (int i = 0; i < m; i++) {
            for (int j : cuts.get(i)) {
                tab[i][j] = 1.0;
            }
            tab[i][n + i] = -1.0;      // surplus
            tab[i][n + m + i] = 1.0;   // artificiale
            tab[i][rhsCol] = 1.0;      // rhs
        }

        double[] obj = new double[ncols + 1];
        for (int j = 0; j < n; j++) obj[j] = cost[j];
        for (int i = 0; i < m; i++) obj[n + m + i] = BIG_M;
        tab[m] = obj;

        int[] basis = new int[m];
        for (int i = 0; i < m; i++) basis[i] = n + m + i; // artificiale i basica nella riga i

        // azzera i coefficienti delle variabili basiche (artificiali) nella riga obiettivo
        for (int i = 0; i < m; i++) {
            double factor = tab[m][n + m + i];
            if (factor != 0.0) {
                for (int k = 0; k <= ncols; k++) {
                    tab[m][k] -= factor * tab[i][k];
                }
            }
        }

        final double EPS = 1e-9;
        final int MAX_ITER = 20000;
        int iter = 0;
        while (iter++ < MAX_ITER) {
            // regola di Bland: variabile entrante = indice piu' piccolo con costo ridotto negativo
            int entering = -1;
            for (int j = 0; j < ncols; j++) {
                if (tab[m][j] < -EPS) { entering = j; break; }
            }
            if (entering == -1) break; // ottimo raggiunto

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
                throw new IllegalStateException("LP illimitato: non dovrebbe accadere in questo modello");
            }

            // pivot
            double pivotVal = tab[leaving][entering];
            for (int k = 0; k <= ncols; k++) tab[leaving][k] /= pivotVal;
            for (int i = 0; i < rows; i++) {
                if (i == leaving) continue;
                double factor = tab[i][entering];
                if (Math.abs(factor) > EPS) {
                    for (int k = 0; k <= ncols; k++) {
                        tab[i][k] -= factor * tab[leaving][k];
                    }
                }
            }
            basis[leaving] = entering;
        }

        double[] x = new double[n];
        for (int i = 0; i < m; i++) {
            if (basis[i] < n) {
                x[basis[i]] = tab[i][rhsCol];
            } else if (basis[i] >= n + m && tab[i][rhsCol] > 1e-6) {
                throw new IllegalStateException("LP infeasibile: non dovrebbe accadere in questo modello");
            }
        }

        double objective = 0.0;
        for (int j = 0; j < n; j++) objective += cost[j] * x[j];
        return new LPResult(objective, x);
    }

    // ========================================================================
    // Warm start: soluzione ammissibile iniziale via euristica greedy
    // (copre iterativamente il cammino s-t "scoperto" piu' economico da chiudere).
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
    // Branch-and-Cut: ricerca ricorsiva in profondita' con backtracking su
    // un array "fixed" (-1 = libero, 0 = fissato a 0, 1 = fissato a 1).
    // ========================================================================

    /**
     * Risolve il problema in modo ESATTO.
     * @param maxRunningTimeMillis tempo massimo di ricerca (safety valve: se
     *        superato, la ricerca si interrompe e restituisce la miglior
     *        soluzione trovata finora, MA isProvenOptimal() tornera' false
     *        perche' l'albero non e' stato completamente esplorato).
     */
    public MCCPSolver.MCCPResult solveExact(long maxRunningTimeMillis) {
        this.searchStartMillis = System.currentTimeMillis();
        this.maxRunningTimeMillis = maxRunningTimeMillis;
        this.nodesExplored = 0;
        this.provenOptimal = true;
        this.globalCuts.clear();

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
        // NOTA: timeToBestMs qui e' approssimato al tempo del warm start se la
        // ricerca esatta non lo ha piu' migliorato; per una stima precisa del
        // vero "time-to-best" sull'incumbent definitivo servirebbe tracciarlo
        // dentro branch(): lo facciamo con bestFoundAtMs qui sotto.
        return new MCCPSolver.MCCPResult(cutColors, keptColors, bestCost, bestFoundAtMs, totalTimeMs);
    }

    private long bestFoundAtMs = 0;

    private void branch(int[] fixed) {
        nodesExplored++;

        if (!provenOptimal) return; // gia' interrotti per timeout, non esplorare oltre
        if (maxRunningTimeMillis > 0 && System.currentTimeMillis() - searchStartMillis > maxRunningTimeMillis) {
            provenOptimal = false;
            return;
        }

        // colori liberi in questo nodo e mappa colore-originale -> indice locale
        List<Integer> free = new ArrayList<>();
        for (int c = 0; c < numColors; c++) if (fixed[c] == -1) free.add(c);
        int[] freeIndex = new int[numColors];
        Arrays.fill(freeIndex, -1);
        for (int i = 0; i < free.size(); i++) freeIndex[free.get(i)] = i;

        double fixed1Cost = 0.0;
        for (int c = 0; c < numColors; c++) if (fixed[c] == 1) fixed1Cost += colorCost[c];

        // filtra i tagli globali noti: rimuove quelli gia' soddisfatti da un
        // colore fissato a 1; se un taglio non ha piu' alcun colore libero
        // (tutti i suoi colori sono fissati a 0), il nodo e' infeasibile.
        List<Set<Integer>> localCuts = new ArrayList<>();
        for (Set<Integer> cut : globalCuts) {
            boolean satisfiedByFixed1 = false;
            for (int c : cut) if (fixed[c] == 1) { satisfiedByFixed1 = true; break; }
            if (satisfiedByFixed1) continue;

            Set<Integer> restricted = new HashSet<>();
            for (int c : cut) if (fixed[c] == -1) restricted.add(freeIndex[c]);
            if (restricted.isEmpty()) return; // nodo infeasibile: pruning
            localCuts.add(restricted);
        }

        double[] localCost = new double[free.size()];
        for (int i = 0; i < free.size(); i++) localCost[i] = colorCost[free.get(i)];

        final int MAX_SEPARATION_ITER = 200;
        for (int iter = 0; iter < MAX_SEPARATION_ITER; iter++) {
            // controllo del tempo GRANULARE: un singolo nodo puo' risolvere la LP
            // fino a MAX_SEPARATION_ITER volte, quindi non basta controllare il
            // tempo solo all'ingresso di branch() (vedi sopra) -- va controllato
            // anche qui dentro, altrimenti su istanze grandi un solo nodo puo'
            // restare "intrappolato" nel ciclo di separazione ben oltre il limite.
            if (maxRunningTimeMillis > 0 && System.currentTimeMillis() - searchStartMillis > maxRunningTimeMillis) {
                provenOptimal = false;
                return;
            }

            LPResult lp = free.isEmpty() ? new LPResult(0.0, new double[0])
                    : solveCoveringLP(free.size(), localCuts, localCost);
            double totalBound = fixed1Cost + lp.objective;

            if (totalBound >= bestCost - 1e-7) return; // pruning

            // costruisce y completo (fissati + LP) per la separazione
            double[] yFull = new double[numColors];
            for (int c = 0; c < numColors; c++) if (fixed[c] == 1) yFull[c] = 1.0;
            for (int i = 0; i < free.size(); i++) yFull[free.get(i)] = lp.x[i];

            SeparationResult sep = dijkstraSeparation(yFull);

            if (sep.colorsOnPath == null || sep.pathWeight >= 1.0 - 1e-7) {
                // nessun cammino violato trovato: la soluzione del rilassamento e' "feasible"
                boolean integral = true;
                for (double v : lp.x) {
                    if (v > 1e-7 && v < 1 - 1e-7) { integral = false; break; }
                }

                if (integral) {
                    if (totalBound < bestCost - 1e-9) {
                        bestCost = totalBound;
                        int[] sol = fixed.clone();
                        for (int i = 0; i < free.size(); i++) {
                            sol[free.get(i)] = lp.x[i] > 0.5 ? 1 : 0;
                        }
                        bestSolution = sol;
                        bestFoundAtMs = System.currentTimeMillis() - searchStartMillis;
                    }
                    return;
                }

                // branching: scegli il colore libero frazionario con valore piu' vicino a 1
                int branchLocalIdx = -1;
                double branchValue = -1.0;
                for (int i = 0; i < lp.x.length; i++) {
                    double v = lp.x[i];
                    if (v > 1e-7 && v < 1 - 1e-7 && v > branchValue) {
                        branchValue = v;
                        branchLocalIdx = i;
                    }
                }
                int branchColor = free.get(branchLocalIdx);

                fixed[branchColor] = 1;
                branch(fixed);
                fixed[branchColor] = 0;
                branch(fixed);
                fixed[branchColor] = -1;
                return;

            } else {
                // VINCOLO VIOLATO. Il taglio VALIDO GLOBALMENTE e' sull'INTERO
                // insieme di colori del cammino trovato (sep.colorsOnPath), non
                // solo sui colori liberi in QUESTO nodo: un colore qui fissato a
                // 0 potrebbe essere libero in un ramo diverso dell'albero. Se
                // memorizzassimo nel pool GLOBALE solo il sottoinsieme ristretto
                // ai colori liberi di questo nodo, il taglio diventerebbe una
                // disuguaglianza troppo forte (quindi NON valida) quando riusato
                // in un altro ramo dove quei colori sono liberi — causando
                // pruning scorretto e soluzioni sub-ottimali dichiarate ottime.
                if (!globalCuts.contains(sep.colorsOnPath)) {
                    globalCuts.add(new HashSet<>(sep.colorsOnPath));
                }

                // Per l'uso IMMEDIATO nella LP di QUESTO nodo, invece, e'
                // corretto restringere ai soli colori liberi: i colori fissati a
                // 0 non possono comunque aiutare in questo sottoalbero (e nessun
                // colore fissato a 1 puo' comparire qui, dato che peso(sep) < 1).
                Set<Integer> restrictedLocal = new HashSet<>();
                for (int c : sep.colorsOnPath) if (fixed[c] == -1) restrictedLocal.add(freeIndex[c]);
                if (restrictedLocal.isEmpty()) return; // nodo infeasibile: pruning
                if (!localCuts.contains(restrictedLocal)) {
                    localCuts.add(restrictedLocal);
                }
                // ripete il ciclo: risolve di nuovo la LP con il nuovo taglio
            }
        }
        // se non si converge entro MAX_SEPARATION_ITER (non dovrebbe succedere), ci si ferma qui per sicurezza
    }

    // ========================================================================
    // Esempio di utilizzo
    // ========================================================================

    public static void main(String[] args) {
        System.out.println("=== Confronto Branch-and-Cut esatto vs VNS-Greedy su istanze casuali piccole/medie ===");
        for (int trial = 0; trial < 5; trial++) {
            int numNodes = 50 + trial * 5;
            int numColors = 8 + trial;
            MCCPSolver instance = MCCPSolver.generateRandomInstance(numNodes, numColors, MCCPSolver.Density.MEDIUM, 0.10);

            MCCPBranchAndCut bc = new MCCPBranchAndCut(instance.getNumNodes(), instance.getEdges(),
                    instance.getNumColors(), instance.getColorCost(), instance.getSourceNode(), instance.getTargetNode());

            long t0 = System.currentTimeMillis();
            MCCPSolver.MCCPResult bcResult = bc.solveExact(10000);
            long t1 = System.currentTimeMillis();

            System.out.println("[trial " + trial + "] " + numNodes + " nodi, " + numColors + " colori");
            System.out.println("Branch-and-Cut (esatto) -> " + bcResult.toString()
                    + "  [nodi albero=" + bc.getNodesExplored() + ", ottimo certificato=" + bc.isProvenOptimal()
                    + ", tempo wall-clock=" + (t1 - t0) + " ms]");
            instance.verifyCutWithBFS(bcResult.cutColors, true);

            MCCPSolver.MCCPResult vnsRes = instance.solve(1000);
            System.out.println("VNS-Greedy (euristico) -> " + vnsRes.toString());
            instance.verifyCutWithBFS(vnsRes.cutColors, true);

            boolean match = Math.abs(bcResult.cutCost - vnsRes.cutCost) < 1e-6;
            System.out.println("VNS-Greedy ha trovato l'ottimo? " + match);
            System.out.println();
        }

        System.out.println("=== Istanza piu' grande (Branch-and-Cut vs VNS-Greedy euristico) ===");
        MCCPSolver grande = MCCPSolver.generateRandomInstance(200, 50, MCCPSolver.Density.MEDIUM, 0.10);
        MCCPBranchAndCut bcGrande = new MCCPBranchAndCut(grande.getNumNodes(), grande.getEdges(),
                grande.getNumColors(), grande.getColorCost(), grande.getSourceNode(), grande.getTargetNode());

        MCCPSolver.MCCPResult bcRes = bcGrande.solveExact(15000);
        System.out.println("Branch-and-Cut (esatto) -> " + bcRes.toString()
                + "  [nodi albero=" + bcGrande.getNodesExplored() + ", ottimo certificato=" + bcGrande.isProvenOptimal() + "]");
        grande.verifyCutWithBFS(bcRes.cutColors, true);

        MCCPSolver.MCCPResult vnsRes = grande.solve(3000);
        System.out.println("VNS-Greedy (euristico) -> " + vnsRes.toString());
        grande.verifyCutWithBFS(vnsRes.cutColors, true);
    }
}