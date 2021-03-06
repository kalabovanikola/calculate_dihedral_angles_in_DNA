#!/usr/bin/python3
import sys
import math
from math import degrees
import gemmi
from gemmi import cif
import numpy as np
import csv
import locale
from copy import *
from locale import atof

"""
must_have = ["O3'", "P", "O5'", "C5'", "C4'", "C3'", "C1'", "N9" or "N1", "O4'", "C2'"]
torsions_def = {
    "alpha": ["O3'", "P", "O5'", "C5'"],
    "beta": ["P", "O5'", "C5'", "C4'"],
    "gamma": ["O5'", "C5'", "C4'", "C3'"],
    "delta": ["C5'", "C4'", "C3'", "O3'"],
    "epsilon": ["C4'", "C3'", "O3'", "P"],
    "zeta": ["C3'", "O3'", "P", "O5'"],
    "chi_pur": ["O4'", "C1'", "N9", "C4"],
    "chi_pyr": ["O4'", "C1'", "N1", "C2"],
    "tau0": ["C4'", "O4'", "C1'", "C2'"],
    "tau1": ["O4'", "C1'", "C2'", "C3'"],
    "tau2": ["C1'", "C2'", "C3'", "C4'"],
    "tau3": ["C2'", "C3'", "C4'", "O4'"],
    "tau4": ["C3'", "C4'", "O4'", "C1'"],
}
############################PB############################
NNdist_def = {
    "NN_RY": ["N9", "N1"],
    "NN_RR": ["N9", "N9"],
    "NN_YR": ["N1", "N9"],
    "NN_YY": ["N1", "N1"],
}

NCCN_torsion_def = {
    "tors_RY": ["N9", "C1'", "C1'", "N1"],
    "tors_RR": ["N9", "C1'", "C1'", "N9"],
    "tors_YR": ["N1", "C1'", "C1'", "N9"],
    "tors_YY": ["N1", "C1'", "C1'", "N1"],
}
"""
# purine and pyrimidine arrays, generated by backassign/standard/sort_RY
purines = ["0AD", "0AV", "0SP", "0UH", "12A", "18M", "1AP", "1MA", "1MG", "2BU", "2EG", "2MA", "2MG", "2PR", "6FK",
           "6HA", "6HB", "6HG", "6IA", "6MA", "6MC", "6MT", "6OG", "6PO", "7AT", "7BG", "7DA", "7GU", "7MG", "8AA",
           "8BA", "8FG", "8MG", "8OG", "8PY", "A", "A23", "A2L", "A2M", "A44", "A5L", "A5O", "A6A", "A6G", "A7E", "A9Z",
           "ABR", "ABS", "AET", "AF2", "AFG", "AP7", "ATP", "AVC", "B8K", "B8W", "B9B", "BGH", "BGM", "C6G", "CG1",
           "DA", "DCG", "DFG", "DG", "DZM", "E", "E1X", "E6G", "E7G", "EDA", "EHG", "F74", "FDG", "FMG", "G25", "G2L",
           "G48", "G7M", "GAO", "GF2", "GMS", "GOM", "GRB", "GSR", "GSS", "GX1", "HN0", "HN1", "IG", "IGU", "KAG",
           "LCG", "LG", "LGP", "M1G", "M2G", "M7A", "MA6", "MA7", "MAD", "MFO", "MG1", "MGV", "MHG", "MIA", "MRG",
           "MTU", "N6G", "O2G", "OGX", "OMG", "P5P", "P7G", "PG7", "PGN", "PGP", "PPW", "PR5", "PRN", "QUO", "R", "RDG",
           "RIA", "S6G", "SDH", "T6A", "V3L", "X", "XPB", "XUA", "XUG", "YG", "YYG", "AD2", "A3P", "DDG", "GDP", "GFL",
           "G", "GNE", "GTP", "TGP", "2FI", "2DA"]
pyrimidines = ["0AP", "0AU", "125", "126", "127", "18Q", "1CC", "1FC", "1RN", "2AT", "2AU", "2BT", "2GT", "2MU", "2NT",
               "2OM", "2OT", "2ST", "3AU", "3ME", "3MU", "47C", "4OC", "4PC", "4PD", "4PE", "4SU", "4U3", "5BU", "5CM",
               "5FC", "5HC", "5HM", "5HU", "5IC", "5IT", "5IU", "5MC", "5MU", "5NC", "5PC", "5PY", "5SE", "64T", "6HC",
               "6HT", "70U", "75B", "77Y", "85Y", "8DT", "8RO", "94O", "9QV", "A5M", "A6C", "A6U", "ATD", "ATL", "B7C",
               "B8Q", "B8T", "B9H", "BOE", "BRU", "C25", "C2L", "C34", "C43", "C5L", "CAR", "CB2", "CBR", "CCC", "CDW",
               "CH", "CTG", "D00", "DC", "DDN", "DFC", "DHU", "DNR", "DRT", "DT", "DU", "DUZ", "E3C", "EIT", "F3H",
               "GMU", "H2U", "HEU", "I4U", "I5C", "IC", "IMC", "IU", "JDT", "JMH", "LC", "LHU", "LSH", "LST", "M5M",
               "MNU", "N5M", "NCU", "NDU", "NMS", "NMT", "NTT", "OHU", "OMC", "OMU", "ONE", "P2T", "P4U", "PDU", "PST",
               "PYO", "RPC", "RSQ", "S2M", "S4U", "SMT", "SUR", "T39", "T3P", "T4S", "T5O", "T5S", "TA3", "TAF", "TC1",
               "TDY", "TED", "TFE", "TLC", "TLN", "TTM", "U25", "U2L", "U2N", "U36", "U8U", "UAR", "UBI", "UBR", "UD5",
               "UFT", "UMS", "UMX", "UPE", "UPV", "UR3", "URX", "US3", "USM", "UVX", "XCY", "YCO", "Z", "ZDU", "CFL",
               "C", "C38", "DOC", "F2T", "ME6", "SPT", "T", "TPC", "TCP", "U", "UMP", "2DT", "CFZ"]


class Coords:
    def __init__(self, row):
        self.x = row[4]
        self.y = row[5]
        self.z = row[6]


class Atom:
    def __init__(self, row):
        self.atom_name = remove_quotation_marks(row[3])
        self.coords = Coords(row)


class Model:
    def __init__(self, number):
        self.entities: [Entity] = []
        self.number = number


class Entity:
    def __init__(self):
        self.sequences: [Sequence] = []
        self.steps: [Step] = []


class Step:
    def __init__(self, model_number, is_the_only_one, seq1=None, seq2=None):
        self.name = get_step_name(seq1, seq2, model_number, is_the_only_one)
        self.angles = Angles(seq1, seq2)


class Angles:
    def __init__(self, seq1=None, seq2=None):
        """
        function computes all the angles for one step if there are enough data provided
        :param seq1: dict structure of the first sequence
        :param seq2: dict structure of the second sequence
        """
        if seq1 is None:
            self.has_angles = False
        else:
            self.has_angles = True
            self.alpha = calculate_dihedral_angles(seq1, seq2, ["O3'"], ["P", "O5'", "C5'"])
            self.beta = calculate_dihedral_angles(seq1, seq2, [], ["P", "O5'", "C5'", "C4'"])
            self.delta = calculate_dihedral_angles(seq1, seq2, ["C5'", "C4'", "C3'", "O3'"], [])
            self.epsilon = calculate_dihedral_angles(seq1, seq2, ["C4'", "C3'", "O3'"], ["P"])
            self.zeta = calculate_dihedral_angles(seq1, seq2, ["C3'", "O3'"], ["P", "O5'"])
            self.gamma = calculate_dihedral_angles(seq1, seq2, [], ["O5'", "C5'", "C4'", "C3'"])
            self.delta2 = calculate_dihedral_angles(seq1, seq2, [], ["C5'", "C4'", "C3'", "O3'"])
            self.chi = calculate_dihedral_angles(seq1, seq2, ["O4'", "C1'", "N9", "C4"], []) \
                if seq1.base_type == "purine" \
                else calculate_dihedral_angles(seq1, seq2, ["O4'", "C1'", "N1", "C2"], [])
            self.chi2 = calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "N9", "C4"]) \
                if seq2.base_type == "purine" \
                else calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "N1", "C2"])
            self.NCCN_tors = calculate_NCCN_torsion(seq1, seq2)
            self.taus = Taus(seq1, seq2)
            self.NNdist = get_NNdist(seq1, seq2)
            self.CCdist = get_distance(seq1.atoms["C1'"].coords, seq2.atoms["C1'"].coords)


class Taus:
    def __init__(self, seq1, seq2):
        self.tau0_1 = calculate_dihedral_angles(seq1, seq2, ["C4'", "O4'", "C1'", "C2'"], [])
        self.tau1_1 = calculate_dihedral_angles(seq1, seq2, ["O4'", "C1'", "C2'", "C3'"], [])
        self.tau2_1 = calculate_dihedral_angles(seq1, seq2, ["C1'", "C2'", "C3'", "C4'"], [])
        self.tau3_1 = calculate_dihedral_angles(seq1, seq2, ["C2'", "C3'", "C4'", "O4'"], [])
        self.tau4_1 = calculate_dihedral_angles(seq1, seq2, ["C3'", "C4'", "O4'", "C1'"], [])
        self.tau0_2 = calculate_dihedral_angles(seq1, seq2, [], ["C4'", "O4'", "C1'", "C2'"])
        self.tau1_2 = calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "C2'", "C3'"])
        self.tau2_2 = calculate_dihedral_angles(seq1, seq2, [], ["C1'", "C2'", "C3'", "C4'"])
        self.tau3_2 = calculate_dihedral_angles(seq1, seq2, [], ["C2'", "C3'", "C4'", "O4'"])
        self.tau4_2 = calculate_dihedral_angles(seq1, seq2, [], ["C3'", "C4'", "O4'", "C1'"])
        self.tau1, self.t1_max = get_tau(self.tau0_1, self.tau1_1, self.tau2_1, self.tau3_1, self.tau4_1)
        self.tau1_type = get_pseudorotation_name(self.tau1)
        self.tau2, self.t2_max = get_tau(self.tau0_2, self.tau1_2, self.tau2_2, self.tau3_2, self.tau4_2)
        self.tau2_type = get_pseudorotation_name(self.tau2)


class Sequence:
    def __init__(self, asym_id, comp_id, seq_number, inscode, is_hetatom):
        self.SequenceVariants: {str: SequenceVariant} = {
            ".": SequenceVariant(asym_id, comp_id, ".", seq_number, inscode, is_hetatom)}


class Name:
    def __init__(self, asym_id, comp_id, alt_id, seq_number, inscode):
        self.seq_number = seq_number
        self.auth_asym_id = asym_id
        self.auth_comp_id = comp_id
        self.label_alt_id = alt_id
        self.inscode = inscode

    # redefining "==" to equality of properties
    def __eq__(self, other):
        if other is None:
            return False
        return (int(self.seq_number), self.auth_asym_id, self.auth_comp_id, self.label_alt_id, self.inscode) == \
               (int(other.seq_number), other.auth_asym_id, other.auth_comp_id, other.label_alt_id, other.inscode)


class SequenceVariant:
    def __init__(self, asym_id, comp_id, alt_id, seq_number, inscode, is_hetatom):
        self.atoms: {str: Atom} = {}
        self.is_valid = True
        self.base_type = None
        self.name = Name(asym_id, comp_id, alt_id, seq_number, inscode)
        self.is_hetatom = is_hetatom


class WrongNumberOfArgumentsException(Exception):
    def __init__(self, amount, message="Not the right amount of arguments"):
        self.amount = amount
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.amount} -> {self.message}'


def get_NNdist(seq1, seq2):
    """
    get the distance between two N on bases
    :param seq1:
    :param seq2:
    :return: euclidean distance
    """
    if seq1.base_type == "purine" and seq2.base_type == "purine":
        return get_distance(seq1.atoms["N9"].coords, seq2.atoms["N9"].coords)
    if seq1.base_type == "pyrimidine" and seq2.base_type == "purine":
        return get_distance(seq1.atoms["N1"].coords, seq2.atoms["N9"].coords)
    if seq1.base_type == "purine" and seq2.base_type == "pyrimidine":
        return get_distance(seq1.atoms["N9"].coords, seq2.atoms["N1"].coords)
    if seq1.base_type == "pyrimidine" and seq2.base_type == "pyrimidine":
        return get_distance(seq1.atoms["N1"].coords, seq2.atoms["N1"].coords)


def get_distance(atom1: Coords, atom2: Coords):
    """

    :param atom1: coords of the first atom
    :param atom2: coords of the second one
    :return: euclidean distance between two atoms
    """
    a = np.array((float(atom1.x), float(atom1.y), float(atom1.z)))
    b = np.array((float(atom2.x), float(atom2.y), float(atom2.z)))
    return round(np.linalg.norm(a - b), 2)


def get_tau(tau0, tau1, tau2, tau3, tau4):
    tau2 = 0.1 if tau2 == 0 else tau2
    tau0, tau1, tau2, tau3, tau4 = normalize_taus([tau0, tau1, tau2, tau3, tau4])
    tan_p = ((tau4 + tau1) - (tau3 + tau0)) / (
                2 * tau2 * (math.sin(36.0 * math.pi / 180.0) + math.sin(72.0 * math.pi / 180.0)))
    p = math.atan(tan_p) * 180.0 / math.pi
    if tau2 < 0:
        p = p + 180.0
    elif tan_p < 0:
        p = p + 360.0
    tmax = abs(tau2 / math.cos(p * math.pi / 180.0))
    return "%.1f" % p, "%.1f" % tmax


def normalize_taus(taus):
    taus = [tau - 360.0 if tau > 180.0 else tau for tau in taus]
    return taus


def get_pseudorotation_name(P):
    """
    For the given pseudorotation angle return its name. Based on Altona, Sundaralingam, JACS 94:23, 1972, 8205-8212
    :param P: pseudorotation
    :return: type of pseudorotation
    """
    P = atof(P)
    if ((P >= 0.0) and (P < 36.0)) or (P == 360.0):
        return "C3end"
    if (P >= 36.0) and (P < 72.0):
        return "C4exo"
    if (P >= 72.0) and (P < 108.0):
        return "O4end"
    if (P >= 108.0) and (P < 144.0):
        return "C1exo"
    if (P >= 144.0) and (P < 180.0):
        return "C2end"
    if (P >= 180.0) and (P < 216.0):
        return "C3exo"
    if (P >= 216.0) and (P < 252.0):
        return "C4end"
    if (P >= 252.0) and (P < 288.0):
        return "O4exo"
    if (P >= 288.0) and (P < 324.0):
        return "C1end"
    if (P >= 324.0) and (P < 360.0):
        return "C2exo"
    return ""


def get_step_name(seq1: SequenceVariant, seq2: SequenceVariant, model_number, is_the_only_one):
    seq1_alt = "." + seq1.name.label_alt_id if seq1.name.label_alt_id != "." else ""
    seq2_alt = "." + seq2.name.label_alt_id if seq2.name.label_alt_id != "." else ""
    model_num = "-m" + str(model_number) if not is_the_only_one else ""
    inscode1 = "." + seq1.name.inscode if seq1.name.inscode != "?" else ""
    inscode2 = "." + seq2.name.inscode if seq2.name.inscode != "?" else ""
    return (
            get_cif_name() + model_num + "_" + seq1.name.auth_asym_id + "_" + seq1.name.auth_comp_id + seq1_alt + "_" +
            seq1.name.seq_number + inscode1 + "_" + seq2.name.auth_comp_id + seq2_alt + "_" + seq2.name.seq_number
            + inscode2)


def calculate_NCCN_torsion(seq1, seq2):
    """
    calculates the NCCN torsion angle, which depends on purine or pyrimidine structure of sequences
    :param seq1: dict structure of the first sequence
    :param seq2: dict structure of the second sequence
    :return: a dihedral angle
    """
    if seq1.base_type == "purine" and seq2.base_type == "purine":
        return calculate_dihedral_angles(seq1, seq2, ["N9", "C1'"], ["C1'", "N9"])
    elif seq1.base_type == "pyrimidine" and seq2.base_type == "purine":
        return calculate_dihedral_angles(seq1, seq2, ["N1", "C1'"], ["C1'", "N9"])
    elif seq1.base_type == "purine" and seq2.base_type == "pyrimidine":
        return calculate_dihedral_angles(seq1, seq2, ["N9", "C1'"], ["C1'", "N1"])
    elif seq1.base_type == "pyrimidine" and seq2.base_type == "pyrimidine":
        return calculate_dihedral_angles(seq1, seq2, ["N1", "C1'"], ["C1'", "N1"])


def calculate_dihedral_angles(seq1: SequenceVariant, seq2: SequenceVariant, atom_names_seq1: [str],
                              atom_names_seq2: [str]):
    """
    calculates a dihedral angles of given points from two following sequences.
    The arrays may be empty, but it should be 4 points all together
    :param seq1: dict structure of the first sequence
    :param seq2: dict structure of the first sequence
    :param atom_names_seq1: the atoms that are needed from the first sequence
    :param atom_names_seq2: the atoms that are needed from the second sequence
    :return: an dihedral angle
    """

    points = [make_gemmi_position_format_from_coords(seq1.atoms[atom_name].coords) for atom_name in atom_names_seq1] + \
             [make_gemmi_position_format_from_coords(seq2.atoms[atom_name].coords) for atom_name in atom_names_seq2]
    if len(points) == 4:
        result = degrees(gemmi.calculate_dihedral(points[0], points[1], points[2], points[3]))
        return round(result, 1) if result >= 0 else round(360 + result, 1)  # angles should be > 0
    else:
        raise WrongNumberOfArgumentsException(len(points))


def remove_quotation_marks(atom):
    """
    solves the problem that cif stores the atoms with ' (ex. O4') with quotation marks, which needs to be removed
    :param atom: the name of the atom
    :return: name without quotation marks if there were any
    """
    return (atom[1:-1]) if '"' in atom else atom


def calculate_angles(entity: Entity, model_number, is_just_one_model):
    """
    Calculates angles if all needed atoms are present in the sequence, for all the sequence variants - for the same alt
    positions and for an alt position with canonical sequence.
    :param entity: one Entity
    :param model_number:
    :param is_just_one_model: for names sake (step ID)
    :return: an array of steps for the entity.
    """
    for i in range(len(entity.sequences) - 1):
        for sequence_var, sequence_var_object in entity.sequences[i].SequenceVariants.items():  # key and value
            if sequence_var_object.is_valid:
                if sequence_var == ".":  # if it is a canonical sequence
                    # calculates angles for all valid seq. variants of + 1 sequence.
                    for sequence_var2 in entity.sequences[i + 1].SequenceVariants.values():
                        if sequence_var2.is_valid:
                            entity.steps.append(
                                Step(model_number, is_just_one_model,
                                     sequence_var_object, sequence_var2))
                else:
                    # calculates angles for the matching alt position (if there is any valid one)
                    if sequence_var in entity.sequences[i + 1].SequenceVariants and \
                            entity.sequences[i + 1].SequenceVariants[sequence_var].is_valid:
                        entity.steps.append(Step(model_number, is_just_one_model, sequence_var_object,
                                                 entity.sequences[i + 1].SequenceVariants[sequence_var]))
                    # calculates angles for the canonical basis (if there is any valid one)
                    if "." in entity.sequences[i + 1].SequenceVariants and \
                            entity.sequences[i + 1].SequenceVariants["."].is_valid:
                        entity.steps.append(Step(model_number, is_just_one_model, sequence_var_object,
                                                 entity.sequences[i + 1].SequenceVariants["."]))


previous_comp_id = ""


def parse_to_entities(table):
    """
    takes the whole structure and parses it to entities
    :param table: gemmi_table with all necessary information
    :return: array of entities
    """
    global previous_comp_id
    array_of_models = []
    numpy_table = np.array(table)
    m, n, k = "", "", ""
    for row in numpy_table:
        if row[10] != k:  # models
            array_of_models.append(Model(row[10]))
            k = row[10]
        if row[1] != n or (row[2] != m and row[2] == "1"):  # entities
            array_of_models[-1].entities.append(Entity())
            n = row[1]
            # append the first sequence of an entity
            array_of_models[-1].entities[-1].sequences.append(
                Sequence(row[8], row[9], row[11], row[12], row[0] == "HETATM"))
            m = row[2]
        if row[2] != m:  # sequences
            array_of_models[-1].entities[-1].sequences.append(
                Sequence(row[8], row[9], row[11], row[12], row[0] == "HETATM"))
            m = row[2]
        if row[7] not in array_of_models[-1].entities[-1].sequences[-1].SequenceVariants:  # new sequence variant
            array_of_models[-1].entities[-1].sequences[-1].SequenceVariants[row[7]] = \
                SequenceVariant(row[8], row[9], row[7], row[11], row[12], row[0] == "HETATM")
        # adding an atom into a sequence variant
        array_of_models[-1].entities[-1].sequences[-1].SequenceVariants[row[7]].atoms.update(
            {remove_quotation_marks(row[3]): Atom(row)})
    handle_altpos(array_of_models)
    return array_of_models


def handle_altpos(array_of_models):
    """
    removes void "." altpositions, adds the atoms from "." to the alt positions
    :param array_of_models: the whole structure
    :return:
    """
    for model in array_of_models:
        for entity in model.entities:
            for sequence in entity.sequences:
                if len(sequence.SequenceVariants) > 1:
                    without_altpos = copy(sequence.SequenceVariants["."])
                    del sequence.SequenceVariants["."]
                    for variant in sequence.SequenceVariants.values():
                        variant.atoms.update(without_altpos.atoms)


def make_gemmi_position_format_from_coords(coords):
    return gemmi.Position(float(coords.x), float(coords.y), float(coords.z))


def operate_cif_file():
    doc = cif.read(sys.argv[1])
    block = doc[0]
    return block.find(
        ["_atom_site.group_PDB", "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.label_atom_id",
         "_atom_site.Cartn_x",
         "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.label_alt_id", "_atom_site.auth_asym_id",
         "_atom_site.auth_comp_id", "_atom_site.pdbx_PDB_model_num", "_atom_site.auth_seq_id",
         "_atom_site.pdbx_PDB_ins_code"])


def check_if_sequence_is_valid_and_add_basic_information(array_of_models):
    """
    the sequence is valid only if it contains all of the atoms given in variable "must_have". The function goes through
    all sequences in all entities and checks if sequences are valid. In addition to that it classifies the base as
    a purine or a pyrimidine.
    :param array_of_models: all entities in given structure
    :return: is_valid and base_type property of a sequence
    """
    must_have = ["O3'", "P", "O5'", "C5'", "C4'", "C3'", "C1'", "O4'", "C2'"]
    for model in array_of_models:
        for entity in model.entities:
            # the first sequence in entity does not have "P", that is why it is handled separately.
            for sequence_var in entity.sequences[0].SequenceVariants.values():
                all_atoms = sequence_var.atoms.keys()
                for atom in must_have:
                    if atom not in all_atoms and atom != "P":
                        sequence_var.is_valid = False
                if "N9" not in all_atoms and "N1" not in all_atoms:  # purine/pyrimidine validity
                    sequence_var.is_valid = False
                sequence_var.base_type = purine_or_pyrimidine(all_atoms)
                if sequence_var.is_hetatom:
                    # for hetatoms there is a special table with vald ones
                    if sequence_var.name.auth_comp_id not in pyrimidines and \
                            sequence_var.name.auth_comp_id not in purines:
                        sequence_var.is_valid = False

            for sequence in entity.sequences[1:]:
                for sequence_var in sequence.SequenceVariants.values():
                    all_atoms = [x for x in sequence_var.atoms.keys()]
                    for atom in must_have:
                        if atom not in all_atoms:
                            sequence_var.is_valid = False
                    if "N9" not in all_atoms and "N1" not in all_atoms:  # purine/pyrimidine validity
                        sequence_var.is_valid = False
                    sequence_var.base_type = purine_or_pyrimidine(all_atoms)
                    if sequence_var.is_hetatom:
                        if sequence_var.name.auth_comp_id not in pyrimidines and \
                                sequence_var.name.auth_comp_id not in purines:
                            sequence_var.is_valid = False


def purine_or_pyrimidine(all_atoms):
    if "N9" in all_atoms and "N1" in all_atoms:
        base_type = "purine"
    elif "N1" in all_atoms:
        base_type = "pyrimidine"
    else:
        base_type = "incomplete sequence"
    return base_type


def split_into_entities_and_calculate_parameters():
    models_array = parse_to_entities(operate_cif_file())
    check_if_sequence_is_valid_and_add_basic_information(models_array)
    for model in models_array:
        for entity in model.entities:
            calculate_angles(entity, models_array.index(model) + 1, len(models_array) == 1)
    return models_array


def get_cif_name():
    # the cif can have ".cif" or ".gz" ending
    if sys.argv[1][-2:] != "gz":
        return sys.argv[1][:-4]
    else:
        return sys.argv[1][:-3]


def write_to_csv():
    with open(get_cif_name() + '.csv', 'w', newline='') as file:
        writer = csv.writer(file, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["step_ID", "d1", "e1", "z1", "a2", "b2", "g2", "d2", "ch1", "ch2", "NN", "CC", "mu", "P1",
                         "t1", "Pn1", "P2", "t2", "Pn2", "nu11", "nu12", "nu13", "nu14", "nu15", "nu21", "nu22", "nu23",
                         "nu24", "nu25"])
        for _model in _models_array:
            for _entity in _model.entities:
                for _atom in _entity.steps:
                    writer.writerow(
                        [_atom.name, _atom.angles.delta, _atom.angles.epsilon, _atom.angles.zeta, _atom.angles.alpha,
                         _atom.angles.beta,
                         _atom.angles.gamma, _atom.angles.delta2, _atom.angles.chi, _atom.angles.chi2,
                         _atom.angles.NNdist, _atom.angles.CCdist, _atom.angles.NCCN_tors, _atom.angles.taus.tau1,
                         _atom.angles.taus.t1_max, _atom.angles.taus.tau1_type, _atom.angles.taus.tau2,
                         _atom.angles.taus.t2_max, _atom.angles.taus.tau2_type, _atom.angles.taus.tau0_1,
                         _atom.angles.taus.tau1_1, _atom.angles.taus.tau2_1, _atom.angles.taus.tau3_1,
                         _atom.angles.taus.tau4_1, _atom.angles.taus.tau0_2,
                         _atom.angles.taus.tau1_2, _atom.angles.taus.tau2_2, _atom.angles.taus.tau3_2,
                         _atom.angles.taus.tau4_2])


_models_array = split_into_entities_and_calculate_parameters()
write_to_csv()
print("done")
