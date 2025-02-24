#sphere

A scalable phylogenetic tree inference algorithm for SARS-CoV-2 sequences. 

```
java -jar sphere.jar 
```

Run SPHERE with the following arguments:  

```
  -i:  input path to aligned sequences (fasta file) from nucleotide alphabet {A,C,T,G}, or "N" for missing positions/ambiguities.  
  -r:  input path to reference sequence with no missing positions/ambiguities ({A,C,T,G} only, no N's)  
  -e:  output path for edges (.csv)  
  -v:  output path for vertices (.csv)  
  -s:  output path for gap-filled sequences (.fasta)  
```

One way to run the sample files would be:

```
java -jar sphere/sphere.jar -i sample_inputs/jan31.fasta -r sample_inputs/ref.fas -v vertices_out.csv -e edges_out.csv -s sequences_out.csv 
```

You can run SPHERE with the precompiled jar in the "sphere" folder.  

Or, if you have Maven installed on your system, run 

```mvn clean install```

to recompile the source code into a new jar. 


Daniel Novikov, Sergey Knyazev, Mark Grinshpon, Pelin Icer, Pavel Skums, and Alex Zelikovsky. Journal of Computational Biology. Nov 2021.1130-1141. http://doi.org/10.1089/cmb.2021.0306
