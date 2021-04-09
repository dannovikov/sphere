import javafx.util.Pair;
import org.biojava.nbio.core.sequence.DNASequence;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

public class SeqAlgs {


    public static int hamDist(LinkedList<Integer> var_pos, DNASequence s1, DNASequence s2) {
        int h_dist = 0;
        for (Integer i : var_pos) {
            if (s1.getCompoundAt(i+1).getShortName().equals("N") || s2.getCompoundAt(i+1).getShortName().equals("N")) {
                continue;
            }
            if (!s1.getCompoundAt(i + 1).equals(s2.getCompoundAt(i + 1))) ++h_dist;
        }
        return h_dist;
    }

    public static int hamDist(HashSet<Integer> var_pos, DNASequence s1, DNASequence s2) {
        int h_dist = 0;
        for (Integer i : var_pos) {
            if (s1.getCompoundAt(i+1).getShortName().equals("N") || s2.getCompoundAt(i+1).getShortName().equals("N")) {
                continue;
            }
            if (!s1.getCompoundAt(i + 1).equals(s2.getCompoundAt(i + 1))) ++h_dist;
        }
        return h_dist;
    }

    public static int hamDist(PhyloNode a, PhyloNode b) {
        int h_dist = 0;
        for (Integer i : a.ref_diffs) { 
            if (b.ref_diffs.contains(i)) {
                if (!(a.seq.getCompoundAt(i+1).equals(b.seq.getCompoundAt(i+1)))) {
                    h_dist+=1;
                }
            } 
            else {
                h_dist +=1;
            }
        }
        for (Integer i : b.ref_diffs) 
            if (!(a.ref_diffs.contains(i)))
                h_dist+=1;
        return h_dist;
    }

    public static int hamDistOnIntersection(PhyloNode a, PhyloNode b) {
        //A, B are the set of positions where sequences a, b differ from reference sequence.
        //Hamming distance = |A| + |B| - 2*|A∩B| + d(A∩B)
        HashSet<Integer> intersection = new HashSet<>(a.ref_diffs);
        intersection.retainAll(b.ref_diffs);
        int h_dist = 0;
        h_dist += a.ref_diffs.size() + b.ref_diffs.size();
        h_dist -= 2 * intersection.size();
        h_dist += hamDist(intersection, a.seq, b.seq);
        return h_dist;
    }


    public static HashSet<Integer> hamDistWithPositions(DNASequence s1, DNASequence s2) {
        int h_dist = 0;
        HashSet<Integer> diff_indices = new HashSet<>();
        for (int i = 0; i < s1.getLength(); i++) {
            if (s1.getCompoundAt(i+1).getShortName().equals("N") || s2.getCompoundAt(i+1).getShortName().equals("N")) {
                continue;
            }
            if (!s1.getCompoundAt(i + 1).equals(s2.getCompoundAt(i + 1))) {
                ++h_dist;
                diff_indices.add(i);
            }
        }
        return diff_indices;
    }



    public static Set<Pair<Integer, Character>> findMutations(LinkedList<Integer> var_pos,
                                                              DNASequence parent, DNASequence child) {
        Set<Pair<Integer, Character>> mutations = new HashSet<>();
        for (Integer i : var_pos) {
            if (parent.getCompoundAt(i+1).getShortName().equals("N") || child.getCompoundAt(i+1).getShortName().equals("N")) continue;
            if (!parent.getCompoundAt(i+1).equals(child.getCompoundAt(i+1))) {
                mutations.add(new Pair<>(i, child.getCompoundAt(i+1).getShortName().charAt(0)));
            }
        }
        return mutations;
    }
}
