# SpotONE

The following repository pertains the code used to develop the back-end construction of the results present in http://www.moreiralab.com/resources/spotone/. If you use this code please cite:

- **Preto, A.J. and Moreira I. S.**, *SPOTONE: hot SPOTs ON protein complexes with Extremely randomized trees via sequence-only features*, 2020, submitted

**Abstract**: *Protein Hot-Spots (HS) are experimentally determined amino-acids, key to small ligand binding and that tend to be structural landmarks on protein-protein interactions. As such, they were extensively approached by structure-based Machine Learning (ML) prediction methods. However, the availability of a much larger array of protein sequences in comparison to determined tree-dimensional structures, indicates that a sequence-based HS predictor has the potential to be more useful for the scientific community. Herein, we present SPOTONE, a new ML predictor, able to accurately classify protein HS via sequence-only features. This algorithm shows an accuracy, AUROC, precision, recall and F1-score of 0.82, 0.83, 0.91, 0.82 and 0.85, respectively, in an independent testing set. The algorithm is deployed within a free-to-use webserver at http://moreiralab.com/resources/spotone, and only requiring the user to submit a FASTA file with one or more protein sequences.*

Provided you have the **psiblast** installation on your own computer, you will need to split the input fasta file in individual fasta files with the individual chains and then run, for each of them, the following, in the terminal:
```
psiblast -query input_fasta.fasta -evalue 0.001 -num_iterations 3 -db path_to_nr_database -outfmt 5 -out input_fasta.txt -out_ascii_pssm input_fasta.pssm -num_threads 32  
```
