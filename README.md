#spy

A phylogenetic tree inference algorithm for SARS-CoV-2 Genomes.

java -jar sphere.jar 

This will show the help menu with all the necessary options and their descriptions.

One way to run the sample files would be:

```
java -jar sphere/sphere.jar -i sample_inputs/c2c.fasta -r sample_inputs/ref.fas -v vertices_out.csv -e edges_out.csv -s sequences_out.fasta 
```


You can use the included precompiled jar in the sphere folder, or if you have Maven installed, run 

```mvn clean install```

to recompile the source code into a new jar. 
