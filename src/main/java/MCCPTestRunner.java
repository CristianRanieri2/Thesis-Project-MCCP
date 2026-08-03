import java.util.*;

/**
 * Test Runner isolato per eseguire e validare i singoli algoritmi MCCP
 * utilizzando la struttura ufficiale di MCCPSolver.java.
 */
public class MCCPTestRunner {

    // Seleziona l'algoritmo da testare singolarmente:
    public enum SolverType {
        BASE,               // Branch & Cut Base (Pure Java)
        PART,               // Modello PART_{s-t} con Google OR-Tools
        PART_WARM,          // Modello PART_{s-t} con Google OR-Tools + warm start da VNS-Probabilistic (MCCPSolver)
        PART_WARM_PARALLEL, // Come PART_WARM, ma con N esecuzioni di VNS-Probabilistic in parallelo per il warm start
        PART_WARM2,         // Modello PART_{s-t} con Google OR-Tools + warm start da VNS-Probabilistic (MCCPSolver2, doppio criterio di terminazione)
        PART_APPROX,        // Versione APPROSSIMATA: budget di tempo totale suddiviso a meta' tra warm start e PART
        VNS                 // Metaeuristica VNS-Greedy nativa di MCCPSolver
    }

/*
    public static void main(String[] args) {

        // ========================================================================
        // 1. CONFIGURAZIONE DEL TEST
        // ========================================================================

        // MODIFICA QUI per scegliere quale solutore eseguire:
        SolverType selectedSolver = SolverType.PART;

        int numNodes = 100;
        int numColors = 25;
        long timeoutMs = vnsWarmStartTimeMillisFromTable(numNodes);     // Timeout in ms (5 secondi per VNS, o più per B&C)
        long seed = 29575L;                                             // Seed per la perfetta riproducibilità del grafo

        // Usato solo da PART_WARM2: numero massimo di iterazioni consecutive
        // senza miglioramento prima che il warm start (MCCPSolver2) si fermi.
        int vnsMaxStagnantIterations = 100;

        // Usato solo da PART_WARM_PARALLEL: numero di esecuzioni concorrenti
        // di VNS-Probabilistic per il warm start (tipicamente = numero di core).
        int numParallelVnsRuns = Runtime.getRuntime().availableProcessors();

        // Usato solo da PART_APPROX: budget di tempo TOTALE, suddiviso a
        // meta' tra warm start e fase PART (limitata nel tempo).
        long approxTotalTimeMillis = 60_000L;

        System.out.println("====== ESECUZIONE ISOLATA ALGORITMO: " + selectedSolver + " ======");
        System.out.println("Nodi: " + numNodes + " | Colori: " + numColors + " | Seed: " + seed);
        System.out.println("------------------------------------------------------------");

        // ========================================================================
        // 2. GENERAZIONE GRAFO USANDO MCCPSolver NATIVO
        // ========================================================================

        MCCPSolver instance = MCCPSolver.generateRandomInstance(
                numNodes,
                numColors,
                MCCPSolver.Density.MEDIUM,
                0.10,
                seed
        );

        System.out.println("Sorgente (s): " + instance.getSourceNode() + " | Pozzo (t): " + instance.getTargetNode());
        System.out.println("Archi totali: " + instance.getNumEdges());
        System.out.println("------------------------------------------------------------");

        // ========================================================================
        // 3. ESECUZIONE DEL SOLUTORE SELEZIONATO
        // ========================================================================

        switch (selectedSolver) {
            case BASE:
                runBaseSolver(instance, timeoutMs);
                break;

            case PART:
                runPartSolver(instance, timeoutMs);
                break;

            case PART_WARM:
                runPartWarmStartSolver(instance, timeoutMs);
                break;

            case PART_WARM_PARALLEL:
                runPartWarmStartParallelSolver(instance, timeoutMs, numParallelVnsRuns);
                break;

            case PART_WARM2:
                runPartWarmStart2Solver(instance, timeoutMs, vnsMaxStagnantIterations);
                break;

            case PART_APPROX:
                runPartApproxSolver(instance, approxTotalTimeMillis);
                break;

            case VNS:
                runVnsSolver(instance, timeoutMs);
                break;
        }
    }

 */




    public static void main(String[] args) {

        // ========================================================================
        // 1. CONFIGURAZIONE GENERALE DEL TEST
        // ========================================================================

        // Modifica qui per scegliere il solutore da testare
        SolverType selectedSolver = SolverType.PART_APPROX;

        int numNodes = 500;                                             // Numero di nodi del grafo |V|
        long timeoutMs = vnsWarmStartTimeMillisFromTable(numNodes);     // Timeout globale in ms
        long seed = 167023L;                                             // Seed per la riproducibilità
        int vnsMaxStagnantIterations = 100;

        // Usato solo da PART_WARM_PARALLEL: numero di esecuzioni concorrenti
        // di VNS-Probabilistic per il warm start (tipicamente = numero di core).
        int numParallelVnsRuns = Runtime.getRuntime().availableProcessors();

        // Usato solo da PART_APPROX: budget di tempo TOTALE, suddiviso a
        // meta' tra warm start e fase PART (limitata nel tempo).
        long approxTotalTimeMillis = timeoutMs;

        // Fattori per i colori:
        double[] colorRatios = {0.5, 1, 1.25 };

        System.out.println("============================================================");
        System.out.println("====== BATTERIA DI TEST SPERIMENTALI: " + selectedSolver + " ======");
        System.out.println("Nodi base (|V|): " + numNodes + " | Seed: " + seed);
        System.out.println("============================================================\n");

        int testCounter = 1;

        // ========================================================================
        // 2. CICLI ANNIDATI SPERIMENTALI (8 TEST TOTALI)
        // ========================================================================

        // 1° FOR: Variazione del numero di colori in funzione del numero di nodi
        for (double ratio : colorRatios) {

            // Calcolo dinamico del numero di colori (arrotondato all'intero più vicino)
            int numColors = (int) Math.round(numNodes * ratio);



            System.out.println("------------------------------------------------------------");
            System.out.printf("TEST #%d | Nodi: %d | Colori: %d (%.0f%% |V|) | Densità: %s%n",
                    testCounter++, numNodes, numColors, ratio * 100, MCCPSolver.Density.MEDIUM);
            System.out.println("------------------------------------------------------------");

            // Generazione dell'istanza del grafo per la combinazione corrente
            MCCPSolver instance = MCCPSolver.generateRandomInstance(
                    numNodes,
                    numColors,
                    MCCPSolver.Density.LOW,
                    0.35,
                    seed
            );

            System.out.println("Sorgente (s): " + instance.getSourceNode() + " | Pozzo (t): " + instance.getTargetNode());
            System.out.println("Archi totali: " + instance.getNumEdges());
            System.out.println("------------------------------------------------------------");

            // Esecuzione del solutore selezionato
            switch (selectedSolver) {
                case BASE:
                    runBaseSolver(instance, timeoutMs);
                    break;

                case PART:
                    runPartSolver(instance, timeoutMs);
                    break;

                case PART_WARM:
                    runPartWarmStartSolver(instance, timeoutMs);
                    break;

                case PART_WARM_PARALLEL:
                    runPartWarmStartParallelSolver(instance, timeoutMs, numParallelVnsRuns);
                    break;

                case PART_WARM2:
                    runPartWarmStart2Solver(instance, timeoutMs, vnsMaxStagnantIterations);
                    break;

                case PART_APPROX:
                    runPartApproxSolver(instance, approxTotalTimeMillis);
                    break;

                case VNS:
                    runVnsSolver(instance, timeoutMs);
                    break;
            }

            System.out.println("\n");
        }

        System.out.println("============================================================");
        System.out.println("====== BATTERIA DI TEST COMPLETATA CON SUCCESSO ======");
        System.out.println("============================================================");
    }





    // ========================================================================
    // ESECUZIONE DEI SINGOLI SOLUTORI
    // ========================================================================

    private static void runBaseSolver(MCCPSolver instance, long timeoutMs) {
        System.out.println("Avvio Branch & Cut BASE...");

        MCCPBranchAndCut bc = new MCCPBranchAndCut(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = bc.solveExact(timeoutMs);
        long totalTime = System.currentTimeMillis() - startTime;

        printReport("B&C-Base", result, bc.getNodesExplored(), bc.isProvenOptimal(), totalTime, instance);
    }

    private static void runPartSolver(MCCPSolver instance, long timeoutMs) {
        System.out.println("Avvio Branch & Cut PART (Google OR-Tools)...");

        MCCPBranchAndCutPART part = new MCCPBranchAndCutPART(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        // Nessun limite di tempo: passando 0, MCCPBranchAndCutPART.solveExact
        // non chiama solver.setTimeLimit(...) (vedi il controllo "if (maxRunningTimeMillis > 0)"
        // al suo interno), quindi OR-Tools gira finche' non dimostra l'ottimo,
        // senza fermarsi per timeout. timeoutMs (parametro del test runner)
        // viene volutamente ignorato qui.
        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = part.solveExact(0);
        long totalTime = System.currentTimeMillis() - startTime;

        printReport("B&C-PART (OR-Tools)", result, part.getNodesExplored(), part.isProvenOptimal(), totalTime, instance);
    }

    private static void runPartWarmStartSolver(MCCPSolver instance, long timeoutMs) {
        System.out.println("Avvio Branch & Cut PART (Google OR-Tools) con warm start da VNS-Probabilistic...");

        MCCPBranchAndCutPARTWarmStart partWarm = new MCCPBranchAndCutPARTWarmStart(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        // Nota: il tempo del warm start (VNS-Probabilistic) e' determinato
        // internamente in base al numero di nodi (Tabella 1 del paper) e resta
        // invariato. La parte OR-Tools, invece, viene lanciata SENZA limite di
        // tempo (0 = nessun setTimeLimit, vedi MCCPBranchAndCutPARTWarmStart):
        // gira finche' non dimostra l'ottimo. timeoutMs viene qui ignorato.
        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = partWarm.solveExactWithVNSWarmStart(instance, 0);
        long totalTime = System.currentTimeMillis() - startTime;

        System.out.println("------------------------------------------------------------");
        System.out.println("Costo Warm Start (VNS):     " + partWarm.getVnsWarmStartCost()
                + "  [" + partWarm.getVnsWarmStartTimeMillisUsed() + " ms]");

        printReport("B&C-PART + Warm Start VNS", result, partWarm.getNodesExplored(),
                partWarm.isProvenOptimal(), totalTime, instance);
    }

    private static void runPartWarmStartParallelSolver(MCCPSolver instance, long timeoutMs, int numParallelVnsRuns) {
        System.out.println("Avvio Branch & Cut PART (Google OR-Tools) con warm start PARALLELO ("
                + numParallelVnsRuns + " esecuzioni concorrenti di VNS-Probabilistic)...");

        MCCPBranchAndCutPARTWarmStart partWarm = new MCCPBranchAndCutPARTWarmStart(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        // Come PART_WARM: fase OR-Tools senza limite di tempo (0); timeoutMs ignorato.
        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = partWarm.solveExactWithVNSWarmStartParallel(instance, 0, numParallelVnsRuns);
        long totalTime = System.currentTimeMillis() - startTime;

        System.out.println("------------------------------------------------------------");
        System.out.println("Costo Warm Start (migliore fra " + numParallelVnsRuns + " run paralleli): "
                + partWarm.getVnsWarmStartCost() + "  [" + partWarm.getVnsWarmStartTimeMillisUsed() + " ms di parete]");

        printReport("B&C-PART + Warm Start VNS Parallelo", result, partWarm.getNodesExplored(),
                partWarm.isProvenOptimal(), totalTime, instance);
    }

    private static void runPartApproxSolver(MCCPSolver instance, long totalTimeMillis) {
        System.out.println("Avvio versione APPROSSIMATA (budget totale = " + totalTimeMillis
                + " ms, suddiviso a meta' tra warm start e fase PART)...");

        MCCPBranchAndCutPARTWarmStart partApprox = new MCCPBranchAndCutPARTWarmStart(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = partApprox.solveApproximate(instance, totalTimeMillis);
        long totalTime = System.currentTimeMillis() - startTime;

        System.out.println("------------------------------------------------------------");
        System.out.println("Costo Warm Start (meta' del budget):     " + partApprox.getVnsWarmStartCost()
                + "  [" + partApprox.getVnsWarmStartTimeMillisUsed() + " ms]");
        System.out.println("Ottimo dimostrato nonostante il tempo limitato? " + partApprox.isProvenOptimal()
                + " (se NO, il risultato e' un'approssimazione, non l'ottimo certificato)");

        printReport("B&C-PART Approssimato (50% warm start + 50% PART)", result, partApprox.getNodesExplored(),
                partApprox.isProvenOptimal(), totalTime, instance);
    }

    private static void runPartWarmStart2Solver(MCCPSolver instance, long timeoutMs, int vnsMaxStagnantIterations) {
        System.out.println("Avvio Branch & Cut PART (Google OR-Tools) con warm start da MCCPSolver2 "
                + "(VNS-Probabilistic con doppio criterio di terminazione: tempo + stagnazione)...");

        // MCCPSolver2 sullo STESSO grafo/costi/s/t dell'istanza gia' generata,
        // cosi' il confronto con gli altri solutori resta sulla stessa istanza.
        MCCPSolver2 instance2 = new MCCPSolver2(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        // Stesso budget di tempo per il warm start della Tabella 1 del paper,
        // usato anche da MCCPBranchAndCutPARTWarmStart.solveExactWithVNSWarmStart(MCCPSolver,...).
        long vnsMaxTimeMillis = vnsWarmStartTimeMillisFromTable(instance.getNumNodes());

        MCCPBranchAndCutPARTWarmStart partWarm2 = new MCCPBranchAndCutPARTWarmStart(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        // La fase OR-Tools gira senza limite di tempo (0), come per PART_WARM;
        // timeoutMs viene qui ignorato.
        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = partWarm2.solveExactWithVNSWarmStart(
                instance2, vnsMaxTimeMillis, vnsMaxStagnantIterations, 0);
        long totalTime = System.currentTimeMillis() - startTime;

        System.out.println("------------------------------------------------------------");
        System.out.println("Costo Warm Start (VNS via MCCPSolver2): " + partWarm2.getVnsWarmStartCost()
                + "  [" + partWarm2.getVnsWarmStartTimeMillisUsed() + " ms]");
        System.out.println("Iterazioni totali (warm start):        " + instance2.getLastRunTotalIterations());
        System.out.println("Warm start fermato per stagnazione?    " + instance2.isLastRunStoppedByStagnation()
                + " (soglia=" + vnsMaxStagnantIterations
                + ", raggiunte=" + instance2.getLastRunStagnantIterationsAtStop() + ")");

        printReport("B&C-PART + Warm Start VNS (MCCPSolver2)", result, partWarm2.getNodesExplored(),
                partWarm2.isProvenOptimal(), totalTime, instance);
    }

    /** Tabella 1 (Bordini & Protti, 2017): tempo massimo per VNS in funzione del numero di nodi. */
    private static long vnsWarmStartTimeMillisFromTable(int numNodes) {
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

    private static void runVnsSolver(MCCPSolver instance, long timeoutMs) {
        System.out.println("Avvio Metaeuristica VNS-Greedy...");

        long startTime = System.currentTimeMillis();
        // Esegue il metodo solve() nativo della tua classe MCCPSolver
        MCCPSolver.MCCPResult result = instance.solve(timeoutMs);
        long totalTime = System.currentTimeMillis() - startTime;

        printReport("VNS-Greedy", result, -1, false, totalTime, instance);
    }

    // ========================================================================
    // REPORT E VERIFICA INDIPENDENTE BFS
    // ========================================================================

    private static void printReport(String name, MCCPSolver.MCCPResult result, long nodes, boolean optimal, long totalTimeMs, MCCPSolver instance) {
        System.out.println("\n================ RISULTATI: " + name + " ================");
        System.out.println("Costo Taglio Trovato:       " + result.cutCost);
        System.out.println("Colori Tagliati:           " + result.cutColors);
        System.out.println("Colori Mantenuti:          " + result.keptColors);
        System.out.println("------------------------------------------------------------");
        System.out.println("Tempo Ultimo Miglioramento: " + result.timeToBestMs + " ms  <-- (timeToBestMs)");
        System.out.println("Tempo Totale Esecuzione:    " + totalTimeMs + " ms");
        System.out.println("------------------------------------------------------------");

        if (nodes >= 0) {
            System.out.println("Nodi B&B Esplorati:        " + nodes);
            System.out.println("Ottimo Dimostrato:         " + (optimal ? "SI" : "NO (Timeout)"));
        }

        // Utilizza direttamente la funzione di verifica BFS nativa della tua classe MCCPSolver
        System.out.print("Verifica BFS Grafo:        ");
        instance.verifyCutWithBFS(result.cutColors, false);
        System.out.println("============================================================\n");
    }
}