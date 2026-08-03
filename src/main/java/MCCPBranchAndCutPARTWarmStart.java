import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import java.util.*;
import java.util.concurrent.*;

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
 *   |V| <=  600  ->  400 s
 *   |V| <=  700  ->  800 s
 *   |V| <=  800  -> 1600 s
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
 *
 * ============================================================================
 * PARALLELIZZAZIONE
 * ============================================================================
 * Due livelli di parallelismo, indipendenti fra loro:
 *
 * 1) MULTITHREADING DEL SOLUTORE (fase esatta): in solveExactInternal, subito
 *    dopo la creazione del solver, si chiama solver.setNumThreads(...) per
 *    abilitare la ricerca B&B multi-thread nativa di CBC/SCIP. Se la tua
 *    versione di OR-Tools non espone questo metodo, sostituiscilo con
 *    solver.setSolverSpecificParametersAsString("parallel/maxnthreads = N")
 *    (sintassi nativa di SCIP), lasciato commentato subito sotto come
 *    alternativa pronta all'uso.
 *
 * 2) WARM START PARALLELO (fase euristica): {@link #solveExactWithVNSWarmStartParallel}
 *    lancia N esecuzioni INDIPENDENTI e CONCORRENTI di VNS-Probabilistic
 *    (una per thread) e tiene la migliore fra tutte, ottenendo un incumbent
 *    iniziale di qualita' piu' alta nello stesso tempo di parete invece che
 *    in un tempo N volte piu' lungo. E' sicuro perche' MCCPSolver non ha
 *    stato mutabile condiviso fra le chiamate (i campi dell'istanza sono
 *    impostati una sola volta nel costruttore; ogni chiamata a
 *    solveProbabilistic crea il proprio Random locale).
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
        if (numNodes <= 600) return 400_000L;
        if (numNodes <= 700) return 800_000L;
        if (numNodes <= 800) return 1600_000L;
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

    /**
     * Come {@link #solveExactWithVNSWarmStart(MCCPSolver, long)}, ma il warm
     * start viene calcolato lanciando {@code numParallelVnsRuns} esecuzioni
     * INDIPENDENTI e CONCORRENTI di VNS-Probabilistic (ciascuna con lo stesso
     * budget di tempo della Tabella 1), tenendo il risultato migliore fra
     * tutte. Essendo VNS-Probabilistic stocastico, esecuzioni indipendenti
     * convergono tipicamente a soluzioni leggermente diverse: eseguirne
     * diverse in parallelo (una per thread) permette di ottenere un
     * incumbent iniziale migliore nello STESSO tempo di parete, invece che
     * un tempo N volte piu' lungo come accadrebbe eseguendole in sequenza.
     *
     * @param numParallelVnsRuns numero di esecuzioni VNS-Probabilistic da
     *                            lanciare in parallelo (tipicamente pari al
     *                            numero di core disponibili)
     */
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

    /**
     * Come {@link #solveExactWithVNSWarmStart(MCCPSolver, long, long)}, ma il
     * warm start viene calcolato con {@link MCCPSolver2#solveProbabilistic(long, int)},
     * cioe' con il doppio criterio di terminazione (tempo E stagnazione)
     * invece del solo criterio temporale di {@link MCCPSolver#solveProbabilistic(long)}.
     *
     * @param instance2                  istanza MCCPSolver2 che rappresenta
     *                                    LO STESSO grafo/costi/s/t con cui e'
     *                                    stato costruito questo solver
     * @param vnsMaxTimeMillis            tempo massimo (ms) concesso al warm start
     * @param vnsMaxStagnantIterations    numero massimo di iterazioni consecutive
     *                                    senza miglioramento prima che il warm
     *                                    start si fermi (deve essere > 0)
     * @param maxRunningTimeMillisForBC   tempo massimo (ms) concesso a OR-Tools
     *                                    DOPO il warm start
     */
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
    // ALGORITMO APPROSSIMATO: tempo totale suddiviso tra warm start e PART
    // ========================================================================

    /**
     * Versione APPROSSIMATA del solutore: a differenza di
     * {@link #solveExactWithVNSWarmStart(MCCPSolver, long)} — dove il tempo
     * del warm start e' fisso (Tabella 1) e la fase PART gira SENZA limite
     * fino a dimostrare l'ottimo — qui si parte da un unico budget di tempo
     * TOTALE che viene suddiviso ESATTAMENTE A META':
     *
     *   - meta' del tempo per VNS-Probabilistic (warm start),
     *   - meta' del tempo per la fase PART con OR-Tools, questa volta con un
     *     LIMITE DI TEMPO effettivo (non piu' illimitato).
     *
     * Poiche' la fase PART ha ora un tempo finito, non e' garantito che
     * arrivi a dimostrare l'ottimo: il risultato restituito puo' quindi
     * essere solo un upper bound (una buona approssimazione), non l'ottimo
     * certificato. Controlla sempre {@link #isProvenOptimal()} dopo la
     * chiamata per sapere se, nonostante il tempo ridotto, l'ottimo e' stato
     * comunque dimostrato (puo' succedere su istanze facili).
     *
     * @param instance         l'istanza del problema (stesso grafo/costi/s/t
     *                          con cui e' stato costruito questo solver)
     * @param totalTimeMillis  budget di tempo TOTALE (ms), suddiviso a meta'
     *                          fra warm start e fase PART
     */
    public MCCPSolver.MCCPResult solveApproximate(MCCPSolver instance, long totalTimeMillis) {
        return solveApproximate(instance, totalTimeMillis, 0.5);
    }

    /**
     * Come {@link #solveApproximate(MCCPSolver, long)}, ma con una frazione
     * di tempo per il warm start configurabile invece del 50% fisso.
     *
     * @param vnsTimeFraction frazione (0,1) di totalTimeMillis da concedere
     *                        al warm start; il resto va alla fase PART
     */
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

    // ========================================================================
    // Costruzione del modello PART_{s-t} e risoluzione con hint
    // ========================================================================

    private MCCPSolver.MCCPResult solveExactInternal(long maxRunningTimeMillis, Set<Integer> warmStartCutColors) {
        long searchStartMillis = System.currentTimeMillis();

        MPSolver solver = MPSolver.createSolver("CBC");
        if (solver == null) {
            solver = MPSolver.createSolver("SCIP");
            System.out.println("Solver:" + solver);
        }
        if (solver == null) {
            throw new RuntimeException("Impossibile caricare i solutori MILP di OR-Tools (CBC/SCIP).");
        }

        // Abilita la ricerca B&B multi-thread nativa del solutore.
        // Prova prima l'opzione A (se la tua versione di OR-Tools la espone
        // direttamente); se NON COMPILA, commentala e usa l'opzione B.
        int numThreads = Runtime.getRuntime().availableProcessors();
        // --- Opzione A ---
        solver.setNumThreads(numThreads);
        // --- Opzione B (alternativa, se la A non compila) ---
        // solver.setSolverSpecificParametersAsString("parallel/maxnthreads = " + numThreads);
        System.out.println("[Solver] multithreading richiesto: " + numThreads + " thread");

        //solver.enableOutput(); // Stampa a schermo i log nativi C++ di CBC/SCIP

        if (maxRunningTimeMillis > 0) {
            solver.setTimeLimit(maxRunningTimeMillis);
        }


        //VERSIONE VELOCE
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

/*
        //VERSIONE LENTA
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
 */

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
        System.out.println("=== TEST: B&C-PART (base) vs B&C-PART con warm start VNS-Probabilistic (parallelo) ===");

        for (int trial = 0; trial < 5; trial++) {
            int numNodes = 80 + trial * 10;
            int numColors = 13 + 5 * trial;

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

            // B&C-PART con warm start VNS-Probabilistic PARALLELO (N run concorrenti).
            // Il tempo del warm start resta quello della Tabella 1 (calcolato
            // internamente in base a |V|); la fase OR-Tools successiva gira
            // anch'essa senza limite di tempo.
            int numParallelVnsRuns = Runtime.getRuntime().availableProcessors();
            MCCPBranchAndCutPARTWarmStart partWarm = new MCCPBranchAndCutPARTWarmStart(
                    instance.getNumNodes(), instance.getEdges(), instance.getNumColors(),
                    instance.getColorCost(), instance.getSourceNode(), instance.getTargetNode());
            long t2 = System.currentTimeMillis();
            MCCPSolver.MCCPResult warmResult = partWarm.solveExactWithVNSWarmStartParallel(instance, 0, numParallelVnsRuns);
            long t3 = System.currentTimeMillis();

            System.out.println("\n[Trial " + trial + "] " + numNodes + " nodi, " + numColors
                    + " colori, s=" + instance.getSourceNode() + ", t=" + instance.getTargetNode());

            System.out.println("B&C-PART (base)          -> costo=" + baseResult.cutCost
                    + " | nodi B&B=" + partBase.getNodesExplored()
                    + " | ottimo certificato=" + partBase.isProvenOptimal()
                    + " | tempo totale=" + (t1 - t0) + " ms");
            instance.verifyCutWithBFS(baseResult.cutColors, false);

            System.out.println("B&C-PART (warm start VNS parallelo, " + numParallelVnsRuns + " thread) -> costo=" + warmResult.cutCost
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