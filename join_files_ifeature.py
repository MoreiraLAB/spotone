#!/usr/bin/env python

"""
Extract sequence based features on a protocol including window size
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "SPOTONN"

import pandas as pd
import os
import sys
from variables import DEFAULT_LOCATION, SYSTEM_SEP, SEP, \
						COMPLEX_NAME_COL, RES_NAME_COL, \
                        RES_NUMBER_COL, CHAIN_NAME_COL, \
                        INPUT_FEATURES_SEQ_NOPYDPI_FILE, \
                        INPUT_FEATURES_SEQ_FILE, RESULTS_FOLDER, \
                        IFEATURE_FEATURES, CLASSES_DICT, \
                        BASE_HEADER, INTERMEDIATE_SEP, \
                        BINARY_CLASS_NAME, IFEATURE_TERMINATION, \
                        IFEATURE_SEP, PDB_FILE_FORMAT, \
                        IFEATURE_ID_COL, AMINO_PROPERTIES_FILE, \
                        AMINO_PROPERTIES_THREE_CODE_COL, \
                        PSEUDO_IFEATURES, IFEATURE_REDUCED_ID_COL, \
                        IFEATURES_PCA_DICT

def locate_ifeatures(input_complex, retrieved_header = False):

	file_location = DEFAULT_LOCATION + SYSTEM_SEP + RESULTS_FOLDER + SYSTEM_SEP
	complex_name_ifeature = PDB_FILE_FORMAT + input_complex
	header, output_row = [], []
	for feature in IFEATURE_FEATURES:
		if IFEATURES_PCA_DICT[feature] != None:
			proper_file_name = file_location + feature + "_reduced" + IFEATURE_TERMINATION
			id_column = IFEATURE_REDUCED_ID_COL
		else:
			proper_file_name = file_location + feature + IFEATURE_TERMINATION
			id_column = IFEATURE_ID_COL
		opened_file = pd.read_csv(proper_file_name, sep = IFEATURE_SEP, header = 0)
		results_row = opened_file.loc[opened_file[id_column] == complex_name_ifeature].drop([id_column], axis = 1)
		if retrieved_header == False: header += list(results_row)
		output_row += list(results_row.values[0])
	for pseudo_feature in PSEUDO_IFEATURES:
		if IFEATURES_PCA_DICT[pseudo_feature] != None:
			proper_file_name = file_location + pseudo_feature + "_reduced" + IFEATURE_TERMINATION
			id_column = IFEATURE_REDUCED_ID_COL
		else:
			proper_file_name = file_location + pseudo_feature + IFEATURE_TERMINATION
			id_column = IFEATURE_ID_COL
		opened_file = pd.read_csv(proper_file_name, sep = IFEATURE_SEP, header = 0)
		results_row = opened_file.loc[opened_file[id_column] == complex_name_ifeature].drop([id_column], axis = 1)
		if retrieved_header == False: header += list(results_row)
		output_row += list(results_row.values[0])
	if retrieved_header == False: return header, output_row
	elif retrieved_header == True: return output_row 

def extract_amino_properties(input_res, retrieved_header = False):

	opened_holder = pd.read_csv(AMINO_PROPERTIES_FILE, sep = SEP, header = 0)
	res_row = opened_holder.loc[opened_holder[AMINO_PROPERTIES_THREE_CODE_COL] == input_res]
	if retrieved_header == False: return list(res_row)[3:], list(res_row.values[0])[3:]
	elif retrieved_header == True: return list(res_row.values[0])[3:]

output_name = "spoton_ifeature_22_01_2020_reduced.csv"
opened_base_file = pd.read_csv("spoton_seq_clean.csv", sep = SEP, header = 0)
output_table = []
written_header = False
count = 1
for index, row in opened_base_file.iterrows():
	identifier = row["0"].lower() + INTERMEDIATE_SEP + row["1"]
	residue_number, residue_type, current_class = row["2"], row["3"], row["4"]
	proper_class = CLASSES_DICT[current_class]
	if written_header == True:
		calculated_features = locate_ifeatures(identifier, retrieved_header = True)
		amino_features = extract_amino_properties(residue_type, retrieved_header = True)
	if written_header == False:
		header_ifeature, calculated_features = locate_ifeatures(identifier, retrieved_header = False)
		header_amino, amino_features = extract_amino_properties(residue_type, retrieved_header = False)
		header = BASE_HEADER  + ["binary_" + str(i) for i in range(1, 21)] + ["rel_dist_1", "rel_dist2"] + header_amino + header_ifeature + [BINARY_CLASS_NAME]
		written_header = True
	output_row = list(row.values) + amino_features + calculated_features + [proper_class]
	output_table.append(output_row)
	print("Current row:", count, " Complex:", identifier)
	count += 1


output_df = pd.DataFrame(output_table, columns = header)
output_df.to_csv(output_name, sep = SEP, index = False)