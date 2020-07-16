#!/usr/bin/env python

"""
Deploys a machine learning model in a classification pipeline
In this case, prediction probability thresholds are updated depending on a specific column (amino acid type)
Several models and datasets can be tested in this pipeline
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__group_leader__email = "irina.moreira@cnc.uc.pt"
__project__ = "SPOTONE"

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, recall_score, precision_score, f1_score, make_scorer
from sklearn.metrics import confusion_matrix
import os
import sys
import numpy as np
from numpy.random import seed
import matplotlib
import random
import pickle
from variables import ID_COLS, RANDOM_SEED, EVALUATION_TRAIN_FILE, \
    EVALUATION_TEST_FILE, CLASS_DISTRIBUTION_FILE, \
    SEP, BINARY_CLASS_NAME, BASE_HEADER, ENSEMBLE_ML_MODEL, \
    AMINO_LIST, RES_NAME_COL, BINARY_CLASS_NAME, \
    RAW_CLASS_NAME, INTERMEDIATE_SEP, PKL_TERMINATION, \
    MODELS_LOCATION, SYSTEM_SEP

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, precision_recall_curve, \
    auc, make_scorer, recall_score, accuracy_score, \
    precision_score, confusion_matrix

from sklearn.ensemble import ExtraTreesClassifier, AdaBoostClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier


os.environ['PYTHONHASHSEED'] = str(RANDOM_SEED)
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)


def prepare_file_split(input_name, amino_list=AMINO_LIST, \
                       split_threshold=0.4):
    """
    Split the dataset into three: identifiers, features and class
    The split is performed according to the unique proteins provided
    """
    opened_table = pd.read_csv(input_name, header=0, sep=SEP)
    HS_table = opened_table.loc[opened_table[RAW_CLASS_NAME] == "HS"]
    NS_table = opened_table.loc[opened_table[RAW_CLASS_NAME] == "NS"]
    HS_dict, NS_dict = {}, {}
    start_merge = False
    amino_list = AMINO_LIST
    for entry in amino_list:
        current_table_HS = HS_table.loc[HS_table[RES_NAME_COL] == entry]
        current_table_NS = NS_table.loc[NS_table[RES_NAME_COL] == entry]
        try:
            HS_train, HS_test = train_test_split(current_table_HS, test_size=split_threshold, \
                                                 random_state=RANDOM_SEED)
            NS_train, NS_test = train_test_split(current_table_NS, test_size=split_threshold, \
                                                 random_state=RANDOM_SEED)
        except:
            train = pd.concat([train, current_table_HS, current_table_NS], axis=0)
            continue
        if start_merge is True:
            train = pd.concat([train, HS_train, NS_train], axis=0)
            test = pd.concat([test, HS_test, NS_test], axis=0)
        if start_merge is False:
            train = pd.concat([HS_train, NS_train], axis=0)
            test = pd.concat([HS_test, NS_test], axis=0)
            start_merge = True
    train_ids = train[BASE_HEADER]
    test_ids = test[BASE_HEADER]
    train_col = train[[BINARY_CLASS_NAME]]
    test_col = test[[BINARY_CLASS_NAME]]
    dropable_columns = BASE_HEADER + [BINARY_CLASS_NAME]
    train_features = train.drop(dropable_columns, axis=1)
    test_features = test.drop(dropable_columns, axis=1)
    return train_ids, test_ids, train_features, test_features, train_col, test_col


class evaluation:
    """
    Input a vector with actual values and predictions to have its evaluation:
    - Accuracy
    - AUC
    - Precision
    - Recall
    - F1-score
    """

    def __init__(self, real_vector, prediction_vector):
        try:
            self.accuracy = round(accuracy_score(real_vector, prediction_vector), 2)
        except:
            self.accuracy = 0.00
        try:
            self.AUC = round(roc_auc_score(real_vector, prediction_vector, average="weighted"), 2)
        except:
            self.AUC = 0.00
        try:
            self.precision = round(precision_score(real_vector, prediction_vector, average="weighted"), 2)
        except:
            self.precision = 0.00
        try:
            self.recall = round(recall_score(real_vector, prediction_vector, average="weighted"), 2)
        except:
            self.recall = 0.00
        try:
            self.f_score = round(f1_score(real_vector, prediction_vector, average="weighted"), 2)
        except:
            self.f_score = 0.00

        self.all_scores = [self.accuracy, self.AUC, self.precision, self.recall, self.f_score]


def prepare_for_tuning(features, ids, class_column, classifier, new_column="raw_pred"):
    """
    Prepare simple operations to output the table ready to be tuned according to amino acid type
    """
    raw_pred = pd.DataFrame(classifier.predict(features))
    pred_proba = classifier.predict_proba(features)
    full_dataset = pd.concat([pd.DataFrame(ids), pd.DataFrame(class_column)], axis=1).reset_index()
    full_dataset[new_column] = raw_pred
    full_dataset = pd.concat([full_dataset, pd.DataFrame(pred_proba)], axis=1)
    return full_dataset


def split_by_amino(input_table, amino_list=AMINO_LIST, res_name_column=RES_NAME_COL):
    """
    Split the table by amino acid type
    """
    output_tables = {}
    for amino_acid in amino_list:
        table_subset = input_table.loc[input_table[res_name_column] == amino_acid]
        output_tables[amino_acid] = table_subset
    return output_tables


def generate_update_indexes(input_amino_table, amino_list=AMINO_LIST, \
                            raw_prediction_column="raw_pred", target_tune_column=1,
                            positive_value=1, negative_value=0, original_class=BINARY_CLASS_NAME,
                            unwanted_prediction=0, wanted_prediction=1):
    """
    Detect new threshold according to amino acid type and maximum 
    """
    tuned_parameters = {}
    for amino_acid in amino_list:
        current_table = input_amino_table[amino_acid]
        positive_subset = current_table.loc[current_table[original_class] == wanted_prediction]
        negative_subset = current_table.loc[current_table[original_class] == unwanted_prediction]
        false_negatives = positive_subset.loc[positive_subset[raw_prediction_column] == negative_value]
        true_negatives = negative_subset.loc[negative_subset[raw_prediction_column] == negative_value]
        false_negatives_max = false_negatives[target_tune_column].max()
        true_negatives_max = true_negatives[target_tune_column].max()
        if false_negatives_max > true_negatives_max:
            tuned_parameters[amino_acid] = 0.5 - false_negatives_max
        else:
            tuned_parameters[amino_acid] = 0.00
    return tuned_parameters


def update_tuned_class(input_table, predictions_column, tuning_parameters, \
                       amino_list=AMINO_LIST, updated_class_name="adapted_pred", \
                       res_name_column=RES_NAME_COL):
    """
    Update the model prediction threshold 
    """
    adapted_predictions = []
    for index, row in input_table.iterrows():
        adapted_value = predictions_column[index][1] + tuning_parameters[row[res_name_column]]
        if adapted_value >= 0.5:
            final_value = 1
        else:
            final_value = 0
        adapted_predictions.append(final_value)
    input_table[updated_class_name] = adapted_predictions
    return input_table


def apply_ML(input_datasets, input_methods, \
             output_name="full_evaluation.csv",
             output_header=["Dataset", "Classifier name", "subset", "Accuracy", "AUC", "Precision", "Recall",
                            "F1-score"], \
             adapted_class_name="adapted_pred", verbose=False, split_threshold=0.5):
    """
    Pipeline to deploy the machine learning models on the datasets provided
    and store the models in the "models_folder" folder using pickle
    """
    output_results = []
    for dataset in list(input_datasets.keys()):
        ids_train, ids_test, X_train, X_test, y_train, y_test = prepare_file_split(input_datasets[dataset])
        ids_test, ids_test_independent, X_test, X_test_independent, y_test, y_test_independent = train_test_split(
            ids_test, X_test, y_test, \
            test_size=split_threshold, random_state=RANDOM_SEED)
        for classifier_name in list(input_methods.keys()):
            classifier_pkl_name = MODELS_LOCATION + SYSTEM_SEP + \
                                  INTERMEDIATE_SEP.join(
                                      (dataset + INTERMEDIATE_SEP + classifier_name).split(" ")) + PKL_TERMINATION
            classifier = input_methods[classifier_name]
            classifier.fit(X_train, y_train.values.ravel())
            train_pred = classifier.predict(X_train)
            test_pred = classifier.predict(X_test)
            train_scores = evaluation(y_train, train_pred).all_scores
            test_scores = evaluation(y_test, test_pred).all_scores
            train_vector = [dataset, classifier_name, "Train"] + train_scores
            test_vector = [dataset, classifier_name, "Test"] + test_scores
            output_results.append(train_vector)
            output_results.append(test_vector)
            try:
                """
                Use the classifier to predict the 
                probability of each sample belonging to each class
                """
                pred_proba_test = classifier.predict_proba(X_test)
                pred_proba_train = classifier.predict_proba(X_train)
                pred_proba_independent = classifier.predict_proba(X_test_independent)

                """
                Prepare all three subsets for probability
                tuning according to amino acid type
                """
                tunable_train = prepare_for_tuning(X_train, ids_train, y_train, classifier)
                tunable_independent = prepare_for_tuning(X_test_independent, ids_test_independent, y_test_independent,
                                                         classifier)
                tunable_test = prepare_for_tuning(X_test, ids_test, y_test, classifier)

                """
                Split the portion of the data to be used for
                probability tuning (test set)
                """
                split_test_aa = split_by_amino(tunable_test)
                tuned_probabilities = generate_update_indexes(split_test_aa)

                """
                Update all prediction probabilities according to new amino acid
                type based thresholds
                """
                updated_test = update_tuned_class(tunable_test, pred_proba_test, tuned_probabilities)
                updated_train = update_tuned_class(tunable_train, pred_proba_train, tuned_probabilities)
                updated_independent = update_tuned_class(tunable_independent, pred_proba_independent,
                                                         tuned_probabilities)

                """
                Evaluate all the datasets with the new thresholds
                """
                updated_train_scores = evaluation(y_train, updated_train[adapted_class_name]).all_scores
                updated_test_scores = evaluation(y_test, updated_test[adapted_class_name]).all_scores
                updated_independent_scores = evaluation(y_test_independent,
                                                        updated_independent[adapted_class_name]).all_scores
                updated_train_vector = [dataset, classifier_name, "Train after tuning"] + updated_train_scores
                updated_test_vector = [dataset, classifier_name, "Test after tuning"] + updated_test_scores
                updated_independent_vector = [dataset, classifier_name,
                                              "Independent Test after tuning"] + updated_independent_scores

                """
                Store the updated vector
                """
                output_results.append(updated_train_vector)
                output_results.append(updated_test_vector)
                output_results.append(updated_independent_vector)
            except:
                if verbose is True:
                    print("Dataset:", dataset, "\nClassifier:", classifier_name, "Was not possible to write")
                continue
            if verbose is True:
                print(train_vector)
                print(test_vector)
                print(updated_train_vector)
                print(updated_test_vector)
                print(updated_independent_vector)
                print(tuned_probabilities)

            """
            Store the models in the "models_folder"
            """
            with open(classifier_pkl_name, 'wb') as f:
                pickle.dump(classifier, f)

    """
    Store a summary of the results in a csv file
    """
    classifier_dataframe = pd.DataFrame(output_results, columns=output_header)
    classifier_dataframe.to_csv(output_name, index=False)


"""
Define the datasets to be used
"""
datasets_dictionary = {"In-house features": "datasets/1_homemade.csv",
                       "PSSM features": "datasets/2_pssm.csv",
                       "In-house + PSSM": "datasets/3_home_pssm.csv",
                       "In-house + iFeatures": "datasets/4_home_ifeature.csv",
                       "In-house + PSSM + iFeatures": "datasets/7_all_features.csv",
                       "iFeatures": "datasets/8_ifeatures.csv"}

"""
Initialize the classification methods to be tested
"""
ML_methods = {"Extreme Randomized Trees": ExtraTreesClassifier(random_state=RANDOM_SEED),
              "Neural Network": MLPClassifier(random_state=RANDOM_SEED),
              "Ada Boosting": AdaBoostClassifier(random_state=RANDOM_SEED),
              "Support Vector Machine": SVC(random_state=RANDOM_SEED)}

apply_ML(datasets_dictionary, ML_methods, verbose=False)
