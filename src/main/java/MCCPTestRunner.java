import java.util.*;

/**
 * Test Runner isolato per eseguire e validare i singoli algoritmi MCCP
 * utilizzando la struttura ufficiale di MCCPSolver.java.
 */
public class MCCPTestRunner {

    // Seleziona l'algoritmo da testare singolarmente:
    public enum SolverType {
        BASE,        // Branch & Cut Base (Pure Java)
        PART,        // Modello PART_{s-t} con Google OR-Tools
        PART_WARM,   // Modello PART_{s-t} con Google OR-Tools + warm start da VNS-Probabilistic
        VNS          // Metaeuristica VNS-Greedy nativa di MCCPSolver
    }

    public static void main(String[] args) {

        // ========================================================================
        // 1. CONFIGURAZIONE DEL TEST
        // ========================================================================

        // MODIFICA QUI per scegliere quale solutore eseguire:
        SolverType selectedSolver = SolverType.PART_WARM;

        int numNodes = 1000;
        int numColors = 500;
        long timeoutMs = 500000000;   // Timeout in ms (5 secondi per VNS, o più per B&C)
        long seed = 12345L;      // Seed per la perfetta riproducibilità del grafo

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

            case VNS:
                runVnsSolver(instance, timeoutMs);
                break;
        }
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