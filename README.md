# SpotONE

## 1. Introduction
The following repository pertains the code used to develop the back-end construction SPOTONE - http://www.moreiralab.com/resources/spotone/. If you use this code please cite:

- **Preto, A.J. and Moreira I. S.**, *SPOTONE: hot SPOTs ON protein complexes with Extremely randomized trees via sequence-only features*, IJMS, 2020, submitted

**Abstract**: *Protein Hot-Spots (HS) are experimentally determined amino-acids, key to small ligand binding and that tend to be structural landmarks on protein-protein interactions. As such, they were extensively approached by structure-based Machine Learning (ML) prediction methods. However, the availability of a much larger array of protein sequences in comparison to determined tree-dimensional structures, indicates that a sequence-based HS predictor has the potential to be more useful for the scientific community. Herein, we present SPOTONE, a new ML predictor, able to accurately classify protein HS via sequence-only features. This algorithm shows an accuracy, AUROC, precision, recall and F1-score of 0.82, 0.83, 0.91, 0.82 and 0.85, respectively, in an independent testing set. The algorithm is deployed within a free-to-use webserver at http://moreiralab.com/resources/spotone, and only requiring the user to submit a FASTA file with one or more protein sequences.*


## 2. Requirements
To follow this pipeline you need to install:
	- **psiblast** executable, along with the non-redundant (**nr**) dataset of proteins sequences, both made available by NCBI.
	- **iFeature** package, available on *https://github.com/Superzchen/iFeature*
	- Python packages **os**, **sys**, **pandas**, **string**, **sklearn**, **pickle** and **numpy**, most of which are available through standard Python installation with Anaconda. It is advised to create a Anaconda environment for this project.

The folder structure should include, inside the **spotone** folder:
	- a **datasets** folder, to store your final datasets combinations.
	- a **resources** folder, location of files **amino_properties.csv** and **encoding.csv**.
	- a **results** folder, to output the iFeature results.

## 3. Setup

Please make sure you adapt the **variables.py** script to your case, namely the **DEFAULT_LOCATION** variable, which will be used throughout the whole pipeline. This variable should be the folder where calculations will be run. No further changes should be needed.

## 4. Feature Extraction

### 4.1. PSSM features
Provided you have the **psiblast** installation on your own computer, you will need to split the input fasta file in individual fasta files with the individual chains and then run, for each of them, the following, in the terminal:

```
psiblast -query input_fasta.fasta -evalue 0.001 -num_iterations 3 -db path_to_nr_database -outfmt 5 -out input_fasta.txt -out_ascii_pssm input_fasta.pssm -num_threads 32  
```

Each of your unique **.fasta** files will generate an output **.pssm** file that you can then process into tables.

### 4.2. In-house features
You need to run an input **.fasta** file with all the sequences, in this case **spotone.fasta** (in the *input_files* folder) with the **feature_extraction.py** script, which will output a **.csv** file with all the features for all the protein chains.

### 4.3. iFeature features
After you have download the iFeature, and put the folder in the same folder where you run the reamining scripts, you call the **call_ifeature.py** script with the same fasta file you used to run the In-house features. This will output to the **results** folder a **.csv** file for each feature.

### 4.4. Join features
In the **input_files** folder, there is the **spoton.csv** file, which you can use to get the class, as described in Moreira, IS, *et al.*, 2017, pertaining the original "SpotOn" algorithm and dataset. After fetching the files, you need to stitch this table together to include the columns regarding the **PSSM**, **In-house** and **iFeature** features. The files used for the protocol steps are on the **datasets** folder.

## 5. Run all Machine Learning Models
Use **spotone_all_models.py** to run default version of Extreme Randomized Trees, Multi-layer Perceptron (Neural Network), AdaBoosting and Support Vector Machine classifiers. Pick the best model to follow up.

## 6. Optimize the best model
Use **spotone_best_model.py** to optmize, with Grid Search, the model you previously selected, as well as detect the best fine tuned amino acid type parameters.
