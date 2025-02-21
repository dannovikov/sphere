import java.io.*;
import java.util.*;

import javafx.util.Pair;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import picocli.CommandLine;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

/**
 * Created by Sergey Knyazev, continued by Daniel Novikov on June 16th, 2020
 */
@CommandLine.Command(name = "sphere", mixinStandardHelpOptions = true, version = "0.0")
public class Main implements Runnable {
    @CommandLine.Option(names = {"-i", "--inFile"}, description = "input fasta file",
            paramLabel = "FILE", required = true)
    private File inputFile;
    @CommandLine.Option(names = {"-r", "--reference"}, description = "reference fasta file containing root sequence",
            paramLabel = "FILE", required = true)
    private File referenceFile;
    @CommandLine.Option(names = {"-e", "--edgeOutFile"},
            description = "output file for edges of phylogeny as an edge list",
            paramLabel = "FILE", required = true)
    private File eOutputFile;
    @CommandLine.Option(names = {"-v", "--vertexOutFile"},
            description = "output path for vertices and their strain IDs",
            paramLabel = "FILE", required = true)
    private File vOutputFile;
//     @CommandLine.Option(names = {"-mp", "--multipleParents"},
//             description = "Use this flag to allow multiple parents, otherwise don't include it", defaultValue="false")
//     private Boolean multiple_parents;
    @CommandLine.Option(names = {"-s", "-seqsOutputFile"},
            description = "output fasta path for gap filled sequences",
            paramLabel = "FILE", required = true)
    private File seqsOutputFile;

    public void run() {
        try {
            // Scanner wait = new Scanner(System.in);
            // System.out.println("Hello.");
            // System.out.println(wait.nextLine());


            final long startTime = System.currentTimeMillis();

            System.out.println("Reading seqs...");

            LinkedHashMap<String, DNASequence> seqs = FastaReaderHelper.readFastaDNASequence(inputFile);

            System.out.println("Reading reference...");

            LinkedHashMap<String, DNASequence> ref_map = FastaReaderHelper.readFastaDNASequence(referenceFile);
            Map.Entry<String, DNASequence> ref = ref_map.entrySet().iterator().next(); 

            System.out.println("Getting SNP positions...");

            LinkedList<Integer> variable_positions = SeqPreproc.getVariablePositions(seqs, ref.getValue());

            System.out.println("Getting distance from reference...");

            Pair<Map<Integer, Set<String>>, Map<String, HashSet<Integer>>> p = SeqPreproc.getHamDistToRef(
                    ref_map.entrySet().iterator().next().getValue(),
                    seqs, variable_positions);

            Map<Integer,Set<String>> h_dist_map = p.getKey();
            Map<String, HashSet<Integer>> ref_diff_positions = p.getValue(); 

            System.out.println("Building tree...");

            PhyloTree g = new PhyloTree(variable_positions);

            g.buildPhylo(h_dist_map, ref_diff_positions, seqs, ref.getKey(), ref.getValue(), false);
            
            g.exportEdgesCsv(eOutputFile);
            g.exportNodesCsv(vOutputFile);
            g.exportSeqsCsv(seqsOutputFile);
            // g.exportDistanceMatrix(dmOutputFile);

            final long endTime = System.currentTimeMillis();
            System.out.println("\n\nComplete. " + (endTime - startTime) + " milliseconds.");
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
    }

}
