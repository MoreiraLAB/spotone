#!/usr/bin/env python

"""
Call the iFeature commands
The iFeature module must be in the same folder as this scripts.
Download it using: git clone https://github.com/Superzchen/iFeature
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "SPOTONE"

from variables import DEFAULT_LOCATION, SYSTEM_SEP, SEP, IFEATURE_FOLDER, \
						RESULTS_FOLDER, TXT_TERMINATION, IFEATURE_FEATURES, \
						IFEATURE_TERMINATION, PSEUDO_IFEATURES, IFEATURES_PCA_DICT
import os
import sys
import pandas as pd

def run_ifeature(input_file, output_prefix = ""):

	output_loc = DEFAULT_LOCATION + SYSTEM_SEP + RESULTS_FOLDER + SYSTEM_SEP
	#excluded: EAAC, BINARY, EGAAC, AAINDEX, ZSCALE, BLOSUM62

	for feature in IFEATURE_FEATURES:
		run_command = "python iFeature.py --file " + input_file + \
						" --type " + feature + " --out " + output_loc + \
						output_prefix + feature + IFEATURE_TERMINATION
		os.system(run_command)
		if IFEATURES_PCA_DICT[feature] != None:
			cluster_command = "python scripts/pcaAnalysis.py --file " + output_loc + \
							output_prefix + feature + IFEATURE_TERMINATION + \
							" --ncomponents " + str(IFEATURES_PCA_DICT[feature]) + " --out " + output_loc + \
							output_prefix + feature + "_reduced" + IFEATURE_TERMINATION
			os.system(cluster_command)
	#excluded: type5, type13
	for pseudo_feature in PSEUDO_IFEATURES:
		run_second_command = "python iFeaturePseKRAAC.py --file " + input_file + \
						" --type " + pseudo_feature + " --gap_lambda 5 --raactype 5 --out " + output_loc + \
						output_prefix + pseudo_feature + IFEATURE_TERMINATION
		os.system(run_second_command)
		if IFEATURES_PCA_DICT[pseudo_feature] != None:
			second_cluster_command = "python scripts/pcaAnalysis.py --file " + output_loc + \
							output_prefix + pseudo_feature + IFEATURE_TERMINATION + \
							" --ncomponents " + str(IFEATURES_PCA_DICT[pseudo_feature]) + " --out " + output_loc + \
							output_prefix + pseudo_feature + "_reduced" + IFEATURE_TERMINATION
			os.system(second_cluster_command)

ifeature_path = DEFAULT_LOCATION + SYSTEM_SEP + IFEATURE_FOLDER
fasta_file = DEFAULT_LOCATION + INPUT_FILES_FOLDER + SYSTEM_SEP + "spotone.fasta"
os.chdir(ifeature_path)
run_ifeature(fasta_file)