from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import freesasa
import os
from os.path import join as join_p
#import predict
from collections import Counter, defaultdict, OrderedDict
from DSSPparser import parseDSSP

ROOT_DIR = os.getcwd()
hydrophobic = ['V', 'I', 'L', 'F', 'M', 'W', 'Y', 'C']
ambigous_aa = ['B', 'X', 'J', 'Z']
negatively_charged = ['D', 'E']
positively_charged = ['R', 'K']
normal_format_pro = ['CYS', 'GLN', 'ILE', 'SER', 'VAL', 'MET', 'ASN', 'PRO', 'LYS', 'THR', 'PHE', 'ALA', 'HIS', 'GLY', 'ASP', 'LEU',
                     'ARG', 'TRP', 'GLU', 'TYR']

# map residue name three letters to one
map_three_one = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
                 "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
                 "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
                 "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
                 "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
                 "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
                 "PTR": "X", "XLE": "X", "XAA": "X", "HSD": "H", "HID": "H",
                 "HSE": "H"}

# map residue name one letter to three
map_one_three = {"G": "GLY", "A": "ALA", "S": "SER", "T": "THR", "C": "CYS",
                 "V": "VAL", "L": "LEU", "I": "ILE", "M": "MET", "P": "PRO",
                 "F": "PHE", "Y": "TYR", "W": "TRP", "D": "ASP", "E": "GLU",
                 "N": "ASN", "Q": "GLN", "H": "HIS", "K": "LYS", "R": "ARG"}

# SASA_sol
map_surface = {'A': 118.1, 'R': 256.0, 'N': 165.5, 'D': 158.7, 'C': 146.1, 'Q': 193.2,
               'E': 186.2, 'G': 88.1, 'H': 202.5, 'I': 181.0, 'L': 193.1, 'K': 225.8,
               'M': 203.4, 'F': 222.8, 'P': 146.8, 'S': 129.8, 'T': 152.5, 'W': 266.3,
               'Y': 236.8, 'V': 164.5, 'X': 88.1}

DSSP_codes = {"H"  # = α-helix
              "B",  # = residue in isolated β-bridge
              "E",  # = extended strand, participates in β ladder
              "G",  # = 3-helix(310 helix)
              "I",  # = 5 helix(π-helix)
              "T",  # = hydrogen bonded turn
              "S",  # = bend
              }

# Kyte & Doolittle index of hydrophobicity
# J. Mol. Biol. 157:105-132(1982).
kd = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
      "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
      "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

# Hydrophilicity
# 1 Hopp & Wood
# Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).
hw = {"A": -0.5, "R": 3.0, "N": 0.2, "D": 3.0, "C": -1.0,
      "Q": 0.2, "E": 3.0, "G": 0.0, "H": -0.5, "I": -1.8,
      "L": -1.8, "K": 3.0, "M": -1.3, "F": -2.5, "P": 0.0,
      "S": 0.3, "T": -0.4, "W": -3.4, "Y": -2.3, "V": -1.5}

# Surface accessibility
# Vergoten G & Theophanides T, Biomolecular Structure and Dynamics,
# pg.138 (1997).
# 1 Emini Surface fractional probability
em = {"A": 0.815, "R": 1.475, "N": 1.296, "D": 1.283, "C": 0.394,
      "Q": 1.348, "E": 1.445, "G": 0.714, "H": 1.180, "I": 0.603,
      "L": 0.603, "K": 1.545, "M": 0.714, "F": 0.695, "P": 1.236,
      "S": 1.115, "T": 1.184, "W": 0.808, "Y": 1.089, "V": 0.606}

# 2 Janin Interior to surface transfer energy scale
ja = {"A": 0.28, "R": -1.14, "N": -0.55, "D": -0.52, "C": 0.97,
      "Q": -0.69, "E": -1.01, "G": 0.43, "H": -0.31, "I": 0.60,
      "L": 0.60, "K": -1.62, "M": 0.43, "F": 0.46, "P": -0.42,
      "S": -0.19, "T": -0.32, "W": 0.29, "Y": -0.15, "V": 0.60}

def calculate_SAS(temp_dict, pdb_path,seq_len):
    struct = freesasa.Structure(pdb_path)
    result = freesasa.calc(struct)
    area_classes = freesasa.classifyResults(result, struct)
    polar = area_classes['Polar']
    apolar = area_classes['Apolar']
    sasa_fraction = (polar+apolar)/ seq_len
    temp_dict.update({"Polar": polar, "Apolar": apolar, "SASA Fraction": sasa_fraction})

def calculate_aa_combos (temp_dict, sequence):
    seq_length = len(sequence)

    lys_arg = (sequence.count("K") + sequence.count("R")) /seq_length
    lys_minus_arg = (sequence.count("K") - sequence.count("R")) /seq_length

    asp_glu = (sequence.count("D") + sequence.count("E")) /seq_length
    asp_minus_glu = (sequence.count("D") - sequence.count("E")) /seq_length

    lys_arg_asp_glu = lys_arg + asp_glu
    lys_arg_asp_minus_glu = lys_arg + asp_minus_glu

    phe_tyr_trp = (sequence.count("F")+ sequence.count("Y") + sequence.count("W")) /seq_length

    aa_combo_dict = {"Lys+Arg":lys_arg, "Asp+Glu":asp_glu, "Lys+Arg+Asp+Glu":lys_arg_asp_glu,"Phe+Tyr+Trp": phe_tyr_trp,
                     "Lys-Arg":lys_minus_arg, "Asp-Glu":asp_minus_glu, "Lys+Arg+Asp-Glu":lys_arg_asp_minus_glu }
    temp_dict.update(aa_combo_dict)

def calculate_physiochemical_features(temp_dict,sequence):
    analyzed_seq = ProteinAnalysis(sequence)

    charge_at_pH7 = analyzed_seq.charge_at_pH(7)
    instability_index = analyzed_seq.instability_index()
    molecular_weight = analyzed_seq.molecular_weight()
    aromaticity = analyzed_seq.aromaticity()
    molar_extinction_coefficient = analyzed_seq.molar_extinction_coefficient()
    gravy = analyzed_seq.gravy() #Grand Average Hyrdopathy - Higher value = More Hydrophobic
    isoelectric_point = analyzed_seq.isoelectric_point()
    helix_fraction, turn_fraction, sheet_fraction = analyzed_seq.secondary_structure_fraction()

    physiochem_dict = {"Charge at pH7": charge_at_pH7, "Instability Index": instability_index,"Molecular Wt": molecular_weight, "Aromaticity": aromaticity, "Molar Extinction Coeff": molar_extinction_coefficient,
    "Gravy":gravy, "Isoelectric pt": isoelectric_point, "Helix Fraction": helix_fraction, "Turn Fraction": turn_fraction,
    "Sheet Fraction": sheet_fraction }
    temp_dict.update(physiochem_dict)

    #Adding separately because get_amino_acids_percent() generates a dictionary on its own
    aa_percent = analyzed_seq.get_amino_acids_percent()
    temp_dict.update(aa_percent)

def calculate_residue_features(temp_dict, sequence):
    analyzed_seq = ProteinAnalysis(sequence)
    aa_percent = analyzed_seq.get_amino_acids_percent()

    hydrophobicity = 0
    hydrophilicity = 0
    interior__surface_transfer_energy_scale = 0
    surface_fractional_probability = 0

    for key in aa_percent.keys():
        hydrophobicity += aa_percent[key]*kd[key]
        hydrophilicity += aa_percent[key]*hw[key]
        surface_fractional_probability += aa_percent[key]*em[key]
        interior__surface_transfer_energy_scale += aa_percent[key]*ja[key]

    temp_dict.update({"Hydrophobicity": hydrophobicity,"Hydrophilicity": hydrophilicity,
                      "Surface Fractional Probability": surface_fractional_probability,
                     "I2S Transfer Energy Scale": interior__surface_transfer_energy_scale})
    temp_dict.update(aa_percent)

def compute_DSSP(temp_dict, pdb_path):
    """compute DSSP and related features"""
    with open(pdb_path, 'r') as pdb_file:
        try:
            sec_str_based_features = compute_dssp_based(pdb_path)
            temp_dict.update(sec_str_based_features)
        except Exception as e:
            print(pdb_path+' failed to extract secondary structure features')
            print(e)


def compute_dssp_based(pdb_path, pathmkdssp="/usr/bin/mkdssp"):
    pdb_file = os.path.basename(pdb_path)
    pdb_name = os.path.splitext(pdb_file)[0]

    allclass = []
    for i in map_surface.keys():
        allclass = allclass+['buried_'+i, 'mod_buried_'+i, 'exposed_'+i]
    pathoutput = "/tmp/{}".format(pdb_file)
    if not os.path.exists(pathoutput):
        os.makedirs(pathoutput)
    try:
        cmd_str = '%s -i %s -o %s/%s.dssp' % (pathmkdssp,
                                              pdb_path, pathoutput, pdb_name)
        ret_code = os.system(cmd_str)
        #print("DSSP run:\n", cmd_str)
        #print ("DSSP ret code", ret_code)
    except Exception as e:
        print(e)
        print("DSSP problem with pdb_file %s".format(pdb_path))
    dssp_file = "%s/%s.dssp" % (pathoutput, pdb_name)
    #fw.write('PDBid\t%s\n' %('\t'.join(['buried_L','buried_charged_sum','buried_aromatic'])))
    parse_dssp = parseDSSP(dssp_file)
    parse_dssp.parse()

    pddict = parse_dssp.dictTodataframe()

    pdb = pdb_name
    chains = set(pddict['chain'].to_list())
    #print(chains)
    daac = defaultdict(list)

    # Denominator: len_protein
    count_pro = set()
    with open(pdb_path, 'r') as f:
        for row in f:
            count_pro.add((row[72:-1], row[22:27].strip()))
        pro_len = len(count_pro)  # make sure this really computes len_protein

    # wt

#parse_.parse() pddict = parse_.dictTodataframe()

    str_res_dict = defaultdict()
    for chain in chains:
        for pd_row in pddict.iterrows():
            row_dict = pd_row[1]
            if row_dict['chain'] == chain:
                resnum = row_dict['resnum']
                res = row_dict['aa']
                if res == 'a':  # correct for lower case cysteines
                    res = 'C'
                if res == '!':
                   continue  # chain break. skip this line
                reschain = row_dict['chain']
                resacc = row_dict['acc']
        #with open(pathoutput+'/'+pdb+'.dssp') as fdssp:
        #    for ldssp in fdssp:

                # if re.match(r'^\s+\d+\s+\d+ {}'.format(chain), ldssp):
                # res = ldssp[13]
                # resnum = ldssp[5:10].strip()
                # reschain = ldssp[11]

#                            resacc = ldssp[35:39].strip()
                if res not in map_surface.keys():
                   print("residue not valid in DSSP file:" + res)
                   continue

                resacc_tri = float(resacc)/map_surface[res]
                str_code = row_dict['struct'].strip()

                if str_code == '' or str_code == 'C':
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
    #print(str_res_dict)
    if 'H' in str_res_dict:
        h_cts = Counter(str_res_dict['H'])
        h_fracs = [h_cts['buried']/pro_len,
                   h_cts['mod_buried']/pro_len, h_cts['exposed']/pro_len]
    else:
        h_fracs = [0.0, 0.0, 0.0]
    if 'B' in str_res_dict:
        b_cts = Counter(str_res_dict['B'])
        b_fracs = [b_cts['buried']/pro_len,
                   b_cts['mod_buried']/pro_len, b_cts['exposed']]
    else:
        b_fracs = [0.0, 0.0, 0.0]
    if 'C' in str_res_dict:
        e_cts = Counter(str_res_dict['C'])
        e_fracs = [e_cts['buried']/pro_len,
                   e_cts['mod_buried']/pro_len, e_cts['exposed']/pro_len]
    else:
        e_fracs = [0.0, 0.0, 0.0]

    bury_locs = ['buried', 'mod_buried', 'exposed']

    aacdic = OrderedDict()
    str_sec = 'helix'
    for loc_i, loc in enumerate(bury_locs):
        aacdic['_'.join([str_sec, loc])] = h_fracs[loc_i]

    str_sec = 'beta'
    for loc_i, loc in enumerate(bury_locs):
        aacdic['_'.join([str_sec, loc])] = b_fracs[loc_i]

    str_sec = 'coil'
    for loc_i, loc in enumerate(bury_locs):
        aacdic['_'.join([str_sec, loc])] = e_fracs[loc_i]

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
