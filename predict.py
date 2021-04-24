"""
The entry point for your prediction algorithm.
"""

from __future__ import annotations
import os
import argparse
import csv
import itertools
from pathlib import Path
import pprint
from typing import Any
import zipfile
import re
import DSSPparser
from joblib import dump, load
from collections import Counter

MODEL_PATH = ""
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.Structure import Structure
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
	#generate_pdb_list([pdb_file])
	init_feat_vec = pdb_to_feat_vec([pdb_file])
	structure = parser.get_structure(pdb_file.stem, pdb_file)

	# featurize + perform inference
	all_features = featurize(structure,init_feat_vec)

	predicted_solubility = ml_inference(all_features)

	return predicted_solubility

def pdb_to_feat_vec (pdb_path):
	pdb_feat_dict = defaultdict()
	with open(pdb_path, 'r') as pdb_file:
		structure = freesasa.Structure(pdb_file)
		result = freesasa.calc(structure)
		area_classes = freesasa.classifyResults(result, structure)

		compute_dssp(pdb_path)
		for record in SeqIO.parse(pdb_file, 'pdb-atom'):
			seq = str(record.seq)
			LysArg, AspGlu, AspGluLysArg, PheTyrTrp = calculate_aa_combos(seq)
			pdb_name = pdb.rstrip(".pdb")
			if test_mode:
				pdb_feat_dict.update({"PDB File": pdb_name,  "Sequence": seq, "Length": len(seq),
								   "Lys+Arg/Len": LysArg, "Asp+Glu/Len": AspGlu, "Asp+Glu+Lys+Arg/Len": AspGluLysArg,
								   "Phe+Tyr+Trp/Len": PheTyrTrp,
								   "Polar": area_classes['Polar'],
								   "Apolar": area_classes['Apolar']})
			else:
				pdb_feat_dict.update({"PDB File": pdb_name, "Solubility Score": solubility_map[pdb_name], "Sequence": seq, "Length": len(seq),
								   "Lys+Arg/Len": LysArg, "Asp+Glu/Len": AspGlu, "Asp+Glu+Lys+Arg/Len": AspGluLysArg,
								   "Phe+Tyr+Trp/Len": PheTyrTrp,
								   "Polar": area_classes['Polar'],
								   "Apolar": area_classes['Apolar']})
	return(pdb_feat_dict)



def structure_based_feats(pdb_path):
	sec_str_feats = compute_dssp(pdb_path)


def compute_dssp(pdb_path, pathmkdssp="/usr/bin/mkdssp"):
	pdb_file = os.basename(pdb_path)
	pdb_name = os.path.splitext(pdb_file,".pdb")[0]


	allclass = []
	for i in map_surface.keys():
		allclass = allclass+['buried_'+i, 'mod_buried_'+i, 'exposed_'+i]
	pathoutput ="/tmp/{}".format(pdb_file)
	try:
		os.system('%s %s %s/%s.dssp' %
				(pathmkdssp, pdb_file, pathoutput, pdb_name))
	except Exception as e:
		print (e)
		print("DSSP problem with pdb_file %s".format(pdb_path))
	dssp_file = "%s/%s.dssp"
	#fw.write('PDBid\t%s\n' %('\t'.join(['buried_L','buried_charged_sum','buried_aromatic'])))
	parse_dssp = parseDSSP(dssp_file)
	parse_dssp.parse()


	pddict = parser.dictTodataframe()

	pdb = pdb_name
	chains = 1
	daac = defaultdict(list)

	# Denominator: len_protein
	count_pro = set()
	with open(pdb_path, 'r') as f:
		for row in f:
			count_pro.add((row[72:-1], row[22:27].strip()))
		pro_len = len(count_pro) #make sure this really computes len_protein

	# wt

#parse_.parse() pddict = parse_.dictTodataframe()

	str_res_dict = defaultdict()
	for chain in chains:

		for row_dict in pddict.iterrows():
			if row_dict['chain']==chain:
				resnum = row_dict['resnum']
				res = row_dict['aa']
				reschain = row_dict['chain']
				resacc = row_dict['acc']
		#with open(pathoutput+'/'+pdb+'.dssp') as fdssp:
		#    for ldssp in fdssp:

				# if re.match(r'^\s+\d+\s+\d+ {}'.format(chain), ldssp):
							# res = ldssp[13]
							# resnum = ldssp[5:10].strip()
							# reschain = ldssp[11]

#                            resacc = ldssp[35:39].strip()
				resacc_tri = float(resacc)/map_surface[res]
				str_code = row_dict['struct']
				if str_code =='' or str_code =='C':
					str_code == 'C'
				cur_str_val = str_res_dict.get(str_code, [])

				if resacc_tri <= 0.2:
					resloc = 'buried'
					cur_str_val.append(resloc)
					str_res_dict.update({str_code: cur_str_val})
				elif resacc_tri < 0.5:
					resloc = 'mod_buried'
					cur_str_val.append(resloc)
					str_res_dict.update({str_code: cur_str_val})
				else:
					resloc = 'exposed'
					cur_str_val.append(resloc)
					str_res_dict.update({str_code: cur_str_val})
				daac[resloc+'_'+res].append(reschain+'_'+resnum)
	#for H, B and E get fracs
	h_cts = Counter(str_res_dict['H'])
	b_cts = Counter(str_res_dict['B'])
	e_cts = Counter(str_res_dict['C'])
	h_fracs = [h_cts['buried']/pro_len,
			   h_cts['mod_buried']/pro_len, h_cts['exposed']/pro_len]
	bury_locs = ['buried','mod_buried', 'exposed']

	aacdic = defaultdict(float)
	str_sec = 'helix'
	for loc_i,loc in enumerate(bury_locs):
		aacdic['_'.join(str_sec,loc)] = h_fracs[loc_i]
	b_fracs = [b_cts['buried']/pro_len,
			   b_cts['mod_buried']/pro_len, b_cts['exposed']]
	str_sec = 'beta'
	for loc_i, loc in enumerate(bury_locs):
		aacdic['_'.join(str_sec, loc)] = b_fracs[loc_i]
	e_fracs = [e_cts['buried']/pro_len,
			   e_cts['mod_buried']/pro_len, e_cts['exposed']/pro_len]
	str_sec = 'coil'
	for loc_i, loc in enumerate(bury_locs):
		aacdic['_'.join(str_sec, loc)] = e_fracs[loc_i]


	# Amino acid composition

	for aac in allclass:
		if aac in daac.keys():
			aacdic[aac] = float(len(daac[aac]))/pro_len
		else:
			aacdic[aac] = 0

	# classify
	#aacresult = [aacdic['buried_L'],
	#				aacdic['buried_K']+aacdic['buried_R']+aacdic['buried_D']+aacdic['buried_E'],  # charged sum
	#				aacdic['buried_F']+aacdic['buried_W']+aacdic['buried_Y']]  # aromatic
	# just return all aacdic for now
	return(aacdic)

def calculate_aa_combos (sequence):
	seq_length = len(sequence)

	lys_arg = round((sequence.count("K") + sequence.count("R")) /seq_length, 3)
	asp_glu = round((sequence.count("D") + sequence.count("E")) /seq_length, 3)
	asp_glu_lys_arg = round(lys_arg + asp_glu, 3)
	phe_tyr_trp = round((sequence.count("F")+ sequence.count("Y") + sequence.count("W")) /seq_length, 3)

	return lys_arg, asp_glu, asp_glu_lys_arg, phe_tyr_trp


def featurize(structure: Structure,nonstruct_feats: list[Any]) -> list[Any]:
	"""
	Calculates 3D ML features from the `structure`.
	"""

	# get all the residues
	residues = [res for res in structure.get_residues()]

	# calculate some random 3D features (you should be smarter here!)
	protein_length = residues[1]["CA"] - residues[-2]["CA"]
	angle = calc_dihedral(
		residues[1]["CA"].get_vector(),
		residues[2]["CA"].get_vector(),
		residues[-3]["CA"].get_vector(),
		residues[-2]["CA"].get_vector(),
	)

	# create the feature vector


	features = [protein_length, angle]


	features.extend(nonstruct_feats)

	return (features)


def ml_inference(features: list[Any]) -> float:
	"""
	This would be a function where you normalize/standardize your features and
	then feed them to your trained ML model (which you would load from a file).
	"""

	loaded_model = load(MODEL_PATH)
	pred_single = loaded_model.predict(features)

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
