import java.util.*;

/**
 * Test Runner isolato per eseguire e validare i singoli algoritmi MCCP
 * utilizzando la struttura ufficiale di MCCPSolver.java.
 */
public class MCCPTestRunner {

    public enum SolverType {
        BASE,               // Branch & Cut Base (Pure Java)
        PART,               // Modello PART_{s-t} con Google OR-Tools
        PART_WARM,          // Modello PART_{s-t} con Google OR-Tools + warm start da VNS-Probabilistic
        PART_WARM_PARALLEL, // Come PART_WARM, ma con N esecuzioni di VNS-Probabilistic in parallelo
        PART_WARM2,         // Modello PART_{s-t} + warm start da MCCPSolver2
        PART_APPROX,        // Versione APPROSSIMATA con controlli temporali capillari
        VNS,                // Metaeuristica VNS-Greedy nativa
        VNS_PROBABILISTIC   // Metaeuristica VNS-Probabilistica nativa
    }

    public static void main(String[] args) {

        // Seleziona qui il solutore da testare (es. SolverType.VNS_PROBABILISTIC)
        SolverType selectedSolver = SolverType.VNS_PROBABILISTIC;

        int numNodes = 700;
        long timeoutMs = vnsWarmStartTimeMillisFromTable(numNodes);
        long seed = 29575L;
        int vnsMaxStagnantIterations = 100;

        int numParallelVnsRuns = Runtime.getRuntime().availableProcessors();
        long approxTotalTimeMillis = timeoutMs;

        double[] colorRatios = { 1 };
        MCCPSolver.Density[] densitys = {MCCPSolver.Density.LOW, MCCPSolver.Density.MEDIUM, MCCPSolver.Density.HIGH};

        System.out.println("============================================================");
        System.out.println("====== BATTERIA DI TEST SPERIMENTALI: " + selectedSolver + " ======");
        System.out.println("Nodi base (|V|): " + numNodes + " | Seed: " + seed);
        System.out.println("============================================================\n");

        int testCounter = 1;

        for (double ratio : colorRatios) {

            for(MCCPSolver.Density density : densitys){

                int numColors = (int) Math.round(numNodes * ratio);

                System.out.println("------------------------------------------------------------");
                System.out.printf("TEST #%d | Nodi: %d | Colori: %d (%.0f%% |V|) | Densità: %s%n",
                        testCounter++, numNodes, numColors, ratio * 100, density);
                System.out.println("------------------------------------------------------------");

                MCCPSolver instance = MCCPSolver.generateRandomInstance(
                        numNodes,
                        numColors,
                        density,
                        0.2,
                        seed
                );

                System.out.println("Sorgente (s): " + instance.getSourceNode() + " | Pozzo (t): " + instance.getTargetNode());
                System.out.println("Archi totali: " + instance.getNumEdges());
                System.out.println("------------------------------------------------------------");

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

                    case VNS_PROBABILISTIC:
                        runVnsProbabilisticSolver(instance, timeoutMs);
                        break;
                }

                System.out.println("\n");
            }
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
        System.out.println("Avvio Branch & Cut PART con warm start PARALLELO (" + numParallelVnsRuns + " thread)...");

        MCCPBranchAndCutPARTWarmStart partWarm = new MCCPBranchAndCutPARTWarmStart(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

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
        System.out.println("Avvio versione APPROSSIMATA con controlli temporali rigidi (budget totale = " + totalTimeMillis
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
        System.out.println("Avvio Branch & Cut PART con warm start da MCCPSolver2...");

        MCCPSolver2 instance2 = new MCCPSolver2(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        long vnsMaxTimeMillis = vnsWarmStartTimeMillisFromTable(instance.getNumNodes());

        MCCPBranchAndCutPARTWarmStart partWarm2 = new MCCPBranchAndCutPARTWarmStart(
                instance.getNumNodes(),
                instance.getEdges(),
                instance.getNumColors(),
                instance.getColorCost(),
                instance.getSourceNode(),
                instance.getTargetNode()
        );

        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = partWarm2.solveExactWithVNSWarmStart(
                instance2, vnsMaxTimeMillis, vnsMaxStagnantIterations, 0);
        long totalTime = System.currentTimeMillis() - startTime;

        System.out.println("------------------------------------------------------------");
        System.out.println("Costo Warm Start (VNS via MCCPSolver2): " + partWarm2.getVnsWarmStartCost()
                + "  [" + partWarm2.getVnsWarmStartTimeMillisUsed() + " ms]");

        printReport("B&C-PART + Warm Start VNS (MCCPSolver2)", result, partWarm2.getNodesExplored(),
                partWarm2.isProvenOptimal(), totalTime, instance);
    }

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
        MCCPSolver.MCCPResult result = instance.solve(timeoutMs);
        long totalTime = System.currentTimeMillis() - startTime;

        printReport("VNS-Greedy", result, -1, false, totalTime, instance);
    }

    private static void runVnsProbabilisticSolver(MCCPSolver instance, long timeoutMs) {
        System.out.println("Avvio Metaeuristica VNS-Probabilistica...");

        long startTime = System.currentTimeMillis();
        MCCPSolver.MCCPResult result = instance.solveProbabilistic(timeoutMs);
        long totalTime = System.currentTimeMillis() - startTime;

        printReport("VNS-Probabilistic", result, -1, false, totalTime, instance);
    }

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

        System.out.print("Verifica BFS Grafo:        ");
        instance.verifyCutWithBFS(result.cutColors, false);
        System.out.println("============================================================\n");
    }
}