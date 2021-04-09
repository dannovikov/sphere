import javafx.util.Pair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.traverse.BreadthFirstIterator;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

import java.io.File;
import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class PhyloTree {
    Graph<PhyloNode, PhyloEdge> g;
    public static LinkedList<Integer> var_pos;

    LinkedList<PhyloNode> nodes = new LinkedList<PhyloNode>();
    public static PriorityQueue<PhyloNode> pq = new PriorityQueue<PhyloNode>();
    static LinkedList<PhyloNode> collapsed_nodes = new LinkedList<PhyloNode>();
    static PhyloNode last;

    public PhyloTree(LinkedList<Integer> var_pos) {
        g = new DefaultDirectedGraph<>(PhyloEdge.class);
        this.var_pos = var_pos;
    }

    public void buildPhylo(Map<Integer, Set<String>> h_dist_map, 
                           Map<String, HashSet<Integer>> ref_diff_positions,
                           LinkedHashMap<String, DNASequence> seqs,
                           String ref_id, DNASequence ref_seq, boolean multiple_parents) {

        //Create root node from reference sequence
        PhyloNode root = new PhyloNode(ref_id, ref_seq, 0, new HashSet<Integer>());
        nodes.add(root);

        //Create nodes for all sequences in order of distance to root
        Object[] distance_keys = h_dist_map.keySet().toArray();
        Arrays.sort(distance_keys);

        int i = 1;
        for (Object distance_key : distance_keys) 
        {
            Set<String> seq_ids = h_dist_map.get(distance_key);
            for (String id : seq_ids) {
                HashSet<Integer> ref_diffs = ref_diff_positions.get(id);
                nodes.add(makeNode(id, 
                                   seqs.get(id),
                                   Integer.parseInt(distance_key.toString()),
                                   ref_diffs));
            }
            //System.out.println(String.format("Creating nodes with dist=%d to root. There are %d such sequences.", distance_key, h_dist_map.get(distance_key).size()));
        }
        //last = root;

        //Add all nodes to priority queue. 
        for(PhyloNode n : nodes)
        {
            n.parents.add(root);
            pq.add(n);
        }
        
        //Build Tree
        PhyloNode v;
        while (!pq.isEmpty()) {
            v = pq.poll(); 
            if(v.asInt() == root.asInt()) {
                g.addVertex(v);
                continue;
            } 
            // if chooseParents finds a parent with hamming distance 0, then it will collapse the node to that parent
            // otherwise it will find the best parent, set the parent of v, and return false.
            boolean node_is_collapsed = chooseParents(v, g, multiple_parents);
            if(!node_is_collapsed) 
            { 
                // last = v;
                // System.out.println("Adding " + v.parents.size() + "edges.");

                g.addVertex(v);
                // add edges from all parents to v
                boolean fill = true;  //fill gaps only from first parent
                for (PhyloNode parent: v.parents) 
                {
                    addEdgeAndFillGaps(parent, v, fill);
                    fill = false;
                }
            }
            //Limiting output to once every thousand seqs
            if (i < 1000 || i % 1000 == 0 || i >= (seqs.size()+1) - 1000) 
                System.out.print(String.format("Inserting sequence #%d/%d into graph with %d nodes\r",i,seqs.size()+1, g.vertexSet().size()));
            ++i;
        }
    }



    private PhyloNode makeNode(String seq_name, DNASequence seq, int h_dist_ref, HashSet<Integer> ref_diffs) 
    {
        PhyloNode n = new PhyloNode(seq_name, seq, h_dist_ref, ref_diffs);
        return n;
    }


    private static boolean chooseParents(PhyloNode x,
                                         Graph<PhyloNode, PhyloEdge> g,
                                         boolean multiple_parents )
    {
        //System.out.println("Finding parents for node: " + x);
        Object[] tree_array = g.vertexSet().toArray();
        //int count = 0;
        //System.out.println();
        int max_dist_ref = -99;
        for(int i = g.vertexSet().size()-1; i >= 0; i--) {
            //count += 1;

            PhyloNode candidate = (PhyloNode) tree_array[i];
            if (candidate.asInt() == x.asInt()) continue;

            int distance_to_candidate = SeqAlgs.hamDist(x, candidate);
            int distance_to_current_parent = SeqAlgs.hamDist(x.parents.getFirst(), x);

            if(distance_to_candidate == 0) {
                collapseNodeToParent(x, candidate);
                return true;
            }

            //Constrains multi-parent search to same distance from ref as first found parent
            if (max_dist_ref != -99 && candidate.dist_ref != max_dist_ref) {
                return false;
            }

            if (candidate.dist_ref + distance_to_candidate == x.dist_ref && candidate.asInt() != 0) {
                if (distance_to_candidate < distance_to_current_parent)
                    x.parents.clear();
                x.parents.add(candidate); 
                
                if (!multiple_parents) {
                    //System.out.println("found in:" + count + " steps");
                    return false;
                }
                else {
                    max_dist_ref = candidate.dist_ref;
                }
            }
            //System.out.print(String.format("Searching for parent #%d\r", count));
        }
        //System.out.println("No parent candidate found");
        //System.out.println("found in:" + count + " steps");
        return false;
    }





    private void addEdgeAndFillGaps(PhyloNode source, PhyloNode target, boolean fill) 
    {
        if(fill) {
            //System.out.println("Filling gaps in node: " + target + " from " + source +  "...");
            LinkedList<Integer> changedPositions = fillGapsFromParent(target, source);
        }
        Set<Pair<Integer, Character>> mutations = SeqAlgs.findMutations(var_pos, source.seq, target.seq);
        //System.out.println("Adding edge from " + source + " to " + target);
        // System.out.println(mutations.size() + " mutations.");
        g.addEdge(source, target, new PhyloEdge(mutations));
    }


    private static void collapseNodeToParent(PhyloNode v, PhyloNode p) {
        // Instead of adding v to tree, we add its sequence ID to its parent p
        //System.out.println("Collapsing node: " + v + " to " + p);
        String seq_id = v.seq_ids.stream().findFirst().get();
        DNASequence seq = v.getSeq();
        p.updateSeq(seq_id, seq);
        v.parents.clear();
        v.parents.add(p); 
        collapsed_nodes.add(v);
    }


    private LinkedList<Integer> fillGapsFromParent(PhyloNode node, PhyloNode parent)
    {
        char[] node_seq = node.seq.getSequenceAsString().toUpperCase().toCharArray();
        LinkedList<Integer> changedPositions = new LinkedList<Integer>();
        for(int i = 0; i < node.seq.getLength(); i++)
        {
            if (node.seq.getCompoundAt(i+1).getShortName().equals("N"))
            {
                if (!parent.seq.getCompoundAt(i+1).getShortName().equals("N"))
                {
                    // Fill position i from parent
                    
                    node_seq[i] = parent.seq.getCompoundAt(i+1).getShortName().charAt(0);
                    changedPositions.add(i);
                }
            }
        }
        try {
            node.seq = new DNASequence(new String(node_seq));
        } catch (CompoundNotFoundException e) {e.printStackTrace();}
        
        //System.out.println(changedPositions.size() + " gap positions filled.");
        return changedPositions;
    }


    public void exportEdgesCsv(File out_file) throws FileNotFoundException {
        PrintWriter f = new PrintWriter(out_file);
        f.println("Source,Target,Number_mutations");
        for (PhyloEdge e: g.edgeSet()) {
            PhyloNode source = g.getEdgeSource(e);
            PhyloNode target = g.getEdgeTarget(e);
            f.println(String.format("%s,%s,%d", source, target, e.mutations.size()));
        }
        f.close();
    }
    public void exportNodesCsv(File out_file) throws FileNotFoundException {
        PrintWriter f = new PrintWriter(out_file);
        f.println("Strain,Vertex");
        for (PhyloNode v: g.vertexSet())
            for(String s: v.seq_ids) f.println(String.format("%s,%s", s, v));
        f.close();
    }
    public void exportSeqsCsv(File out_file) throws FileNotFoundException {
        //Save node_seqs
        PrintWriter f = new PrintWriter(out_file);
        f.println("Node,Seq");
        for(PhyloNode node : g.vertexSet())
            f.println( node + "," + node.seq);
        f.close();
    }  
}

