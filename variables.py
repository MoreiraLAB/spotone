#!/usr/bin/env python

"""
Store overall usage variables for the Hot-spot determination with sequence only features
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "SPOTONE"

"""
Path variables
"""

TXT_TERMINATION = ".txt"
IFEATURE_TERMINATION = ".tsv"
IFEATURE_SEP = "\t"
SYSTEM_SEP = "/"
DEFAULT_LOCATION = "/home/imoreira/Desktop/DL_chapter"
RESOURCES_FOLDER = "resources"
RESULTS_FOLDER = "results"
SUPPORT_FOLDER = "support"
EVALUATION_FOLDER = "evaluation"
IFEATURE_FOLDER = "iFeature"
AMINO_PROPERTIES_FILE = DEFAULT_LOCATION + SYSTEM_SEP + RESOURCES_FOLDER + SYSTEM_SEP + "amino_properties.csv"
ENCODING_FILE = DEFAULT_LOCATION + SYSTEM_SEP + RESOURCES_FOLDER + SYSTEM_SEP + "encoding.csv"
RAW_FILE = DEFAULT_LOCATION + SYSTEM_SEP + "spoton.csv"
PDB_FOLDER_NAME = "PDB"
EVALUATION_TRAIN_FILE = EVALUATION_FOLDER + SYSTEM_SEP + "train"
EVALUATION_TEST_FILE = EVALUATION_FOLDER + SYSTEM_SEP + "test"
INPUT_FEATURES_FILE = DEFAULT_LOCATION + SYSTEM_SEP + "spoton_clean.csv"
INPUT_FEATURES_SEQ_FILE = DEFAULT_LOCATION + SYSTEM_SEP + "spoton_seq_clean.csv"
INPUT_FEATURES_SEQ_NOPYDPI_FILE = DEFAULT_LOCATION + SYSTEM_SEP + "spoton_seq_no_pydpi.csv"
WINDOWS_FEATURES_FILE = DEFAULT_LOCATION + SYSTEM_SEP + "spoton_windows.csv"
INPUT_CLASS_FILE = DEFAULT_LOCATION + SYSTEM_SEP + "class_clean.csv"
CLASS_DISTRIBUTION_FILE = EVALUATION_FOLDER + SYSTEM_SEP + "class_distribution.csv"
FINAL_MODEL = "final_model"
PKL_TERMINATION = ".pkl"
MODELS_LOCATION = "models_folder"
ENCODING = "utf-8"
SEP = ","
PDB_FILE_FORMAT = "pdb"
PDB_FILE_TERMINATION = ".ent"
INTERMEDIATE_SEP = "_"

"""
Table associated variables
"""
ID_COLS = ["0","1","2","3","4","5"]
CENTRED_FEATURES_WINDOWS = [2, 5, 7, 10, 25, 50, 75]
SECTIONS_FEATURES_SPLITS = [100]
FEATURES_HOLDER_COLS = ['Helix_propensity', 'Sheet_propensity', 'Helix_propensity_values', 'Sheet_propensity_values', 'MW', 'pka_carboxylate', 'pka_amine', 'pka_side_chain', 'C_number', 'H_number', 'N_number', 'O_number', 'S_number', 'Standard_area_free', 'Standard_area_protein', 'folded_buried_area', 'mean_fractional_area_loss', 'residue_mass', 'monoisotopic_mass']
IFEATURE_FEATURES = ['AAC', 'CKSAAP', 'DPC', 'DDE', 'TPC',
				'GAAC', 'CKSAAGP', 'GDPC', 'GTPC',
				'NMBroto', 'Moran', 'Geary',
				'CTDC', 'CTDT', 'CTDD', 'CTriad', 'KSCTriad', 'SOCNumber',
				'QSOrder', 'PAAC', 'APAAC']
PSEUDO_IFEATURES = ['type1', 'type2', 'type3A', 'type3B', 'type4', 'type6A', 'type6B', 'type6C',
								 'type7', 'type8', 'type9', 'type10', 'type11', 'type12', 'type14', 'type15',
								 'type16']
IFEATURES_PCA_DICT = {'AAC': None, 'CKSAAP': 50, 'DPC': 40, 'DDE': 60, 'TPC': 100,
				'GAAC': None, 'CKSAAGP' : 40, 'GDPC': None, 'GTPC': 40,
				'NMBroto': 50, 'Moran': 50, 'Geary': 40,
				'CTDC': 20, 'CTDT': 20, 'CTDD': 40, 'CTriad': 40, 'KSCTriad': 40, 'SOCNumber': 30,
				'QSOrder': 30, 'PAAC': 30, 'APAAC': 30, 'type1': None, 'type2': None, 'type3A': None,
				'type3B': 20, 'type4': None, 'type6A': None, 'type6B': None, 'type6C': None, 'type7': None,
				'type8': None, 'type9': None, 'type10': None, 'type11': None, 'type12': None, 'type14': None,
				'type15': None, 'type16': None}
CLASSES_DICT = {"NS": 0, "HS": 1}
BASE_HEADER = ["CPX", "PDBChain", "PDBResNo", "PDBResName", "Classe"]
AMINO_LIST = ['ASP', 'SER', 'GLN', 'LYS',
                 'ILE', 'PRO', 'THR', 'CYS', 'PHE', 'ASN', 
                 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 
                 'ALA', 'VAL', 'GLU', 'MET', 'TYR']
BINARY_CLASS_NAME = "binary_class"
RAW_CLASS_NAME = "Classe"
RES_NAME_COL = "PDBResName"
RES_NUMBER_COL = "PDBResNo"
CHAIN_NAME_COL = "PDBChain"
COMPLEX_NAME_COL = "CPX"
ALPHA_CARBON = "CA"
IFEATURE_ID_COL = "#"
IFEATURE_REDUCED_ID_COL = "protein"
AMINO_PROPERTIES_THREE_CODE_COL = "three_code"
FEAT_PYDPI_PROTEIN_TERMINATION = "_protein"
