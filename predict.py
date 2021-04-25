"""
The entry point for your prediction algorithm.
"""

from __future__ import annotations
from Bio.PDB.Polypeptide import PPBuilder
import enum
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from enum import unique
import os
import argparse
import csv
import itertools
from pathlib import Path
import pprint
from typing import Any
import zipfile
import re
import freesasa
from DSSPparser import parseDSSP
from joblib import dump, load
from collections import Counter,defaultdict, OrderedDict
import sklearn
import numpy as np
from feature_computation_functions import *
#from keras import model, load_model

#MODEL_PATH = "data/first_model.bin"  # under /home/biolib
MODEL_PATH = "data/lasso_model_v1.bin"  # under /home/biolib
KERAS_MODEL_PATH = "data/dnn_model.h5"
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.Structure import Structure
from Bio import SeqIO
import temppathlib


test_mode  = True
# SASA_sol
map_surface = {'A': 118.1, 'R': 256.0, 'N': 165.5, 'D': 158.7, 'C': 146.1, 'Q': 193.2,
               'E': 186.2, 'G': 88.1, 'H': 202.5, 'I': 181.0, 'L': 193.1, 'K': 225.8,
               'M': 203.4, 'F': 222.8, 'P': 146.8, 'S': 129.8, 'T': 152.5, 'W': 266.3,
               'Y': 236.8, 'V': 164.5, 'X': 88.1}

DSSP_codes = {"H" # = α-helix
              "B", #= residue in isolated β-bridge
              "E",  # = extended strand, participates in β ladder
              "G",  # = 3-helix(310 helix)
              "I",  # = 5 helix(π-helix)
              "T",  # = hydrogen bonded turn
              "S",  # = bend
}

def predict(pdb_file: Path) -> float:
    """
    The function that puts it all together: parsing the PDB file, generating
    features from it and performing inference with the ML model.
    """

    # parse PDB
    parser = PDBParser()

    init_feat_vec = pdb_to_feat_vec(pdb_file)
    structure = parser.get_structure(pdb_file.stem, pdb_file)

    # featurize + perform inference
    all_features = featurize(structure,init_feat_vec)

    predicted_solubility = ml_inference(all_features)

    return predicted_solubility


def pdb_to_feat_vec(pdb_path):
        #Extracting sequence
        pdb_name= pdb_path.parts[-1]
        structure = PDBParser().get_structure(pdb_name, pdb_path)
        ppb = PPBuilder()
        seq = ""  # Re-initializing to prevent any redundancy
        for pp in ppb.build_peptides(structure):
                seq += str(pp.get_sequence())

        #Adding features
        individual_protein_dict = {"PDB File": pdb_name,
                                    "Sequence": seq, "Length": len(seq)}

        calculate_SAS(individual_protein_dict, pdb_path, len(seq))
        calculate_aa_combos(individual_protein_dict, seq)
        calculate_physiochemical_features(individual_protein_dict, seq)
        calculate_residue_features(individual_protein_dict, seq)
        compute_DSSP(individual_protein_dict, str(pdb_path))




def pdb_to_feat_vec (pdb_path):
    pdb_feat_dict = OrderedDict()
    try:
        if pdb_path.is_file():
            structure = freesasa.Structure(str(pdb_path))
            result = freesasa.calc(structure=structure)
            area_classes = freesasa.classifyResults(result, structure)
        else:
            print(str(pdb_path)+' does not exist. failed to extract freeSASA features')
            area_classes = {'Polar': 0.0, "Apolar": 0.0}

    except Exception as e:
        print(str(pdb_path)+' failed to extract freeSASA features')
        print(e)
        area_classes = {'Polar' : 0.0, "Apolar":0.0}
    try:
        sec_str_based_features=compute_dssp_based(str(pdb_path))
    except Exception as exc_obj:
        print(str(pdb_path)+' failed to extract secondary structure features')
        print(exc_obj)
        sec_str_based_features = {}
    with open(str(pdb_path), 'r') as pdb_file:

        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            seq = str(record.seq)
            LysArg, AspGlu, AspGluLysArg, PheTyrTrp = calculate_aa_combos(seq)
            pdb_name = str(pdb_path.parts[-1])

            pdb_feat_dict = OrderedDict([("PDB File", pdb_name),  ("Sequence", seq), ("Length",  len(seq)),
                                   ("Lys+Arg/Len",  LysArg), ("Asp+Glu/Len",  AspGlu), ("Asp+Glu+Lys+Arg/Len",  AspGluLysArg),
                                   ("Phe+Tyr+Trp/Len",  PheTyrTrp),
                                   ("Polar",  area_classes['Polar']),
                                   ("Apolar",  area_classes['Apolar']) ])
            pdb_feat_dict.update(sec_str_based_features)
            # if  not test_mode:
            #     pdb_feat_dict["Solubility Score"] = solubility_map[pdb_name]

    return(pdb_feat_dict)


def calculate_aa_combos (sequence):
    seq_length = len(sequence)

    lys_arg = round((sequence.count("K") + sequence.count("R")) /seq_length, 3)
    asp_glu = round((sequence.count("D") + sequence.count("E")) /seq_length, 3)
    asp_glu_lys_arg = round(lys_arg + asp_glu, 3)
    phe_tyr_trp = round((sequence.count("F")+ sequence.count("Y") + sequence.count("W")) /seq_length, 3)

    return lys_arg, asp_glu, asp_glu_lys_arg, phe_tyr_trp


def featurize(structure: Structure, nonstruct_feats: OrderedDict) -> list[Any]:
    """
    Calculates 3D ML features from the `structure`.
    """


    # calculate some random 3D features (you should be smarter here!)


    # create the feature vector
    remove_columns= ["PDB File", "Sequence"]
    features = [v for (k, v) in nonstruct_feats.items() if k not in remove_columns]
    #print(features)


    return (features)


def ml_inference(features: list[Any]) -> float:
    """
    This would be a function where you normalize/standardize your features and
    then feed them to your trained ML model (which you would load from a file).
    """

    loaded_model = load(MODEL_PATH)

    feature_vec = np.zeros(shape = (1,len(features)) )
    for feat_i, feat_val in enumerate(features):
        feature_vec[0,feat_i] = float(feat_val)
    pred_single = loaded_model.predict(feature_vec)

    return (pred_single)



if __name__ == "__main__":

    # set up argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, default="data/test.zip")
    args = parser.parse_args()

    predictions = []
    # use a temporary directory so we don't pollute our repo
    with temppathlib.TemporaryDirectory() as tmpdir:
        # unzip the file with all the test PDBs
        with zipfile.ZipFile(args.infile, "r") as zip_:
            zip_.extractall(tmpdir.path)

        # iterate over all test PDBs and generate predictions
        for test_pdb in tmpdir.path.glob("*.pdb"):
            predictions.append({"protein": test_pdb.stem, "solubility": predict(test_pdb)})

    # save to csv file, this will be used for benchmarking
    outpath = "predictions.csv"
    with open(outpath, "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=["protein", "solubility"])
        writer.writeheader()
        writer.writerows(predictions)

    # print predictions to screen
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(predictions)
