#!/usr/bin/python3
import sys
import math
from math import degrees
import gemmi
from gemmi import cif
import numpy as np
import csv
print(gemmi.__version__)
"""
must_have = ["O3'", "P", "O5'", "C5'", "C4'", "C3'", "C1'", "N9" or "N1", "O4'", "C2'"]
# which atoms take place in the torsion angle
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


# TODO: neměly by být steps u sekvencí?

class Coords:
    def __init__(self, row):
        self.x = row[4]
        self.y = row[5]
        self.z = row[6]


class Atom:
    def __init__(self, row):
        self.atom_name = remove_quotation_marks(row[3])
        self.coords = Coords(row)


class Entity:
    def __init__(self):
        self.sequences: [Sequence] = []
        self.steps: [Angles] = []


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


class Sequence:
    def __init__(self):
        self.sequence_variants: {str: Sequence_variant} = {}


class Sequence_variant:
    def __init__(self):
        self.atoms: {str: Atom} = {}
        self.is_valid = True
        self.base_type = None


class WrongNumberOfArgumentsException(Exception):
    def __init__(self, amount, message="Not the right amount of arguments"):
        self.amount = amount
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.amount} -> {self.message}'


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


def make_array_form_coords(coords):
    return [float(coords.x), float(coords.y), float(coords.z)]


def calculate_dihedral_angles(seq1: Sequence_variant, seq2: Sequence_variant, atom_names_seq1: [str],
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

        return round(result, 1) if result >= 0 else round(360 + result, 1)
    else:
        raise WrongNumberOfArgumentsException(len(points))


def remove_quotation_marks(atom):
    """
    solves the problem that cif stores the atoms with ' (ex. O4') with quotation marks, which needs to be removed
    :param atom: the name of the atom
    :return: name without quotation marks if there were any
    """
    return (atom[1:-1]) if '"' in atom else atom


def calculate_angles(entity: Entity):
    """
    calculates angles if all needed atoms are present in the sequence, else adds an empty step
    :param entity: one Entity
    :return: an array of steps for the entity.
    """
    for i in range(len(entity.sequences) - 1):
        for sequence_var, sequence_var_object in entity.sequences[i].sequence_variants.items():
            if sequence_var_object.is_valid:
                if sequence_var == ".":
                    for sequence_var2 in entity.sequences[i + 1].sequence_variants.values():
                        if sequence_var2.is_valid:
                            entity.steps.append(
                                Angles(entity.sequences[i].sequence_variants[sequence_var], sequence_var2))
                else:
                    if sequence_var in entity.sequences[i + 1].sequence_variants and \
                            entity.sequences[i + 1].sequence_variants[sequence_var].is_valid:
                        entity.steps.append(Angles(sequence_var_object,
                                                   entity.sequences[i + 1].sequence_variants[sequence_var]))
                    if "." in entity.sequences[i + 1].sequence_variants and entity.sequences[i + 1].sequence_variants[
                        "."].is_valid:
                        entity.steps.append(Angles(sequence_var_object,entity.sequences[i + 1].sequence_variants["."]))


def parse_to_entities(table):
    """
    takes the whole structure and parses it to entities according to cif distribution
    :param table: gemmi_table with all necessary information
    :return: array of entities
    """
    array_of_entities = []
    numpy_table = np.array(table)
    print(numpy_table)
    m, n = 0, 0
    for row in numpy_table:
        if row[0] == "ATOM":
            if int(row[1]) != n:
                array_of_entities.append(Entity())
                n = int(row[1])
            if int(row[2]) != m:
                array_of_entities[-1].sequences.append(Sequence())
                m = int(row[2])
            if row[7] not in array_of_entities[-1].sequences[-1].sequence_variants:
                array_of_entities[-1].sequences[-1].sequence_variants[row[7]] = Sequence_variant()
            array_of_entities[-1].sequences[-1].sequence_variants[row[7]].atoms.update(
            {remove_quotation_marks(row[3]): Atom(row)})
    return array_of_entities


def make_gemmi_position_format_from_coords(coords):
    return gemmi.Position(float(coords.x), float(coords.y), float(coords.z))


def operate_cif_file():
    doc = cif.read(sys.argv[1])
    block = doc[0]
    return block.find(
        ["_atom_site.group_PDB", "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.label_atom_id",
         "_atom_site.Cartn_x",
         "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.label_alt_id"])


def check_if_sequence_is_valid_and_add_basic_information(array_of_entities):
    """
    the sequence is valid only if it contains all of the atoms given in variable "must_have". The function goes through
    all sequences in all entities and checks if sequences are valid. In addition to that it classifies the base as
    a purine or a pyrimidine.
    :param array_of_entities: all entities in given structure
    :return: is_valid and base_type property of a sequence
    """
    must_have = ["O3'", "P", "O5'", "C5'", "C4'", "C3'", "C1'", "O4'", "C2'"]

    for entity in array_of_entities:

        for sequence_var in entity.sequences[0].sequence_variants.values():
            all_atoms = sequence_var.atoms.keys()
            for atom in must_have:
                if atom not in all_atoms and atom != "P":
                    sequence_var.is_valid = False
        for sequence in entity.sequences[1:]:
            for sequence_var in sequence.sequence_variants.values():
                all_atoms = [x for x in sequence_var.atoms.keys()]
                for atom in must_have:
                    if atom not in all_atoms:
                        sequence_var.is_valid = False
                if "N9" not in all_atoms and "N1" not in all_atoms:
                    sequence_var.is_valid = False
                if "N9" in all_atoms and "N1" in all_atoms:
                    sequence_var.base_type = "purine"
                elif "N1" in all_atoms:
                    sequence_var.base_type = "pyrimidine"
                else:
                    sequence_var.base_type = "incomplete sequence"


def split_into_entities_and_calculate_parameters():
    entities_array = parse_to_entities(operate_cif_file())
    check_if_sequence_is_valid_and_add_basic_information(entities_array)
    for entity in entities_array:
        calculate_angles(entity)
    return entities_array


def write_it_out(array_of_entities):
    """
    just a testing function
    :param array_of_entities:
    :return:
    """
    a = 0
    for entity in array_of_entities:
        print("new entity")
        b = 1
        a += 1
        for sequence in entity.sequences:
            print("new_sequence")
            print(f'{a}.{b}')
            print(sequence.base_type)
            print(sequence.is_valid)
            b += 1
            for atom in sequence.atoms.values():
                print(atom.atom_name, end=" ")
                print(atom.coords.x)


_entities_array = split_into_entities_and_calculate_parameters()
# write_it_out(_entities_array)
with open(sys.argv[1][:-3] + 'csv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["d1", "e1", "z1", "a2", "b2", "g2", "d2", "ch1", "ch2", "NCCN_tors"])
    for _entity in _entities_array:
        for _atom in _entity.steps:
            writer.writerow(
                [_atom.delta, _atom.epsilon, _atom.zeta, _atom.alpha, _atom.beta, _atom.gamma, _atom.delta2, _atom.chi,
                 _atom.chi2, _atom.NCCN_tors])
