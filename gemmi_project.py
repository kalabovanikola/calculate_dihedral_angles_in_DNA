#!/bin/python3
import math
from math import degrees
import gemmi
from gemmi import cif
import numpy as np

"""
must_have = ["O3'", "P", "O5'", "C5'", "C4'", "C3'", "C1'", "N9" nebo "N1", "O4'", "C2'"]
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


# TODO: možná smazat Atom a nechat jen Coords - zatím ne, mohlo by se hodit
# TODO: neměly by být steps u sekvencí?

class Coords:
    def __init__(self, row):
        self.x = row[3]
        self.y = row[4]
        self.z = row[5]


class Atom:
    def __init__(self, row):
        self.atom_name = remove_quotation_marks(row[2])
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
                if seq1.pur_or_pyr == "purine" \
                else calculate_dihedral_angles(seq1, seq2, ["O4'", "C1'", "N1", "C2"], [])
            self.chi2 = calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "N9", "C4"]) \
                if seq2.pur_or_pyr == "purine" \
                else calculate_dihedral_angles(seq1, seq2, [], ["O4'", "C1'", "N1", "C2"])
            self.NCCN_tors = calculate_NCCN_torsion(seq1, seq2)


class Sequence:
    def __init__(self):
        self.atoms: {str: Atom} = {}
        self.is_valid = True
        self.pur_or_pyr = None


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
    if seq1.pur_or_pyr == "purine" and seq2.pur_or_pyr == "purine":
        return calculate_dihedral_angles(seq1, seq2, ["N9", "C1'"], ["C1'", "N9"])
    elif seq1.pur_or_pyr == "pyrimidine" and seq2.pur_or_pyr == "purine":
        return calculate_dihedral_angles(seq1, seq2, ["N1", "C1'"], ["C1'", "N9"])
    elif seq1.pur_or_pyr == "purine" and seq2.pur_or_pyr == "pyrimidine":
        return calculate_dihedral_angles(seq1, seq2, ["N9", "C1'"], ["C1'", "N1"])
    elif seq1.pur_or_pyr == "pyrimidine" and seq2.pur_or_pyr == "pyrimidine":
        return calculate_dihedral_angles(seq1, seq2, ["N1", "C1'"], ["C1'", "N1"])


def make_array_form_coords(coords):
    return [float(coords.x), float(coords.y), float(coords.z)]


def calculate_dihedral_angles(seq1: Sequence, seq2: Sequence, atom_names_seq1: [str], atom_names_seq2: [str]):
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
        return result
    else:
        raise WrongNumberOfArgumentsException(len(points))


def remove_quotation_marks(atom):
    """
    solves the problem that cif stores the atoms with ' (ex. O4') with quotation marks, which needs to be removed
    :param atom: the name of the atom
    :return: name without quotation marks if there were any
    """
    return (atom[1:-1]) if '"' in atom else atom


# TODO: Možná předělat prázdné Angles na None... a nebo to vyřešit nějak chytřeji
def calculate_angles(entity: Entity):
    """
    calculates angles if all needed atoms are present in the sequence, else adds an empty step
    :param entity: one Entity
    :return: an array of steps for the entity.
    """
    for i in range(len(entity.sequences) - 1):
        if entity.sequences[i - 1].is_valid and entity.sequences[i].is_valid:
            entity.steps.append(Angles(entity.sequences[i - 1], entity.sequences[i]))
        else:
            entity.steps.append(Angles())


def parse_to_entities(table):
    table_of_entities = []
    numpy_table = np.array(table)
    m,n = 0,0
    for row in numpy_table:
        if int(row[0]) != n:
            table_of_entities.append(Entity())
            n = int(row[0])
        if int(row[1]) != m:
            table_of_entities[-1].sequences.append(Sequence())
            m = int(row[1])
        table_of_entities[-1].sequences[-1].atoms.update({remove_quotation_marks(row[2]): Atom(row)})
    return table_of_entities


def make_gemmi_position_format_from_coords(coords):
    return gemmi.Position(float(coords.x), float(coords.y), float(coords.z))


def is_valid(array_of_entities):
    must_have = ["O3'", "P", "O5'", "C5'", "C4'", "C3'", "C1'", "O4'", "C2'"]
    for entity in array_of_entities:
        for sequence in entity.sequences:
            allatoms = [x for x in sequence.atoms.keys()]
            for atom in must_have:
                if atom not in allatoms:
                    sequence.is_valid = False
            if "N9" not in allatoms and "N1" not in allatoms:
                sequence.is_valid = False
            if "N9" in allatoms and "N1" in allatoms:
                sequence.pur_or_pyr = "purine"
            elif "N1" in allatoms:
                sequence.pur_or_pyr = "pyrimidine"
            else:
                sequence.pur_or_pyr = "incomplete sequence"


def write_it_out(array_of_entities):
    a = 0
    for entity in array_of_entities:
        print("new entity")
        b = 1
        a += 1
        for sequence in entity.sequences:
            print("new_sequence")
            print(f'{a}.{b}')
            print(sequence.pur_or_pyr)
            print(sequence.is_valid)
            b += 1
            for atom in sequence.atoms.values():
                print(atom.atom_name, end=" ")
                print(atom.coords.x)


doc = cif.read_file("cif_example.cif")
block = doc[0]
gemmi_table = block.find(
    ["atom_site.labelentity_id", "atom_site.label_seq_id", "atom_site.labelatom_id", "atom_site.Cartn_x",
     "atom_site.Cartn_y", "atom_site.Cartn_z"])

print(np.array(gemmi_table))
entities_array = parse_to_entities(gemmi_table)
is_valid(entities_array)
write_it_out(entities_array)
print(type(make_gemmi_position_format_from_coords(entities_array[0].sequences[1].atoms["P"].coords)))
for _entity in entities_array:
    calculate_angles(_entity)
for _entity in entities_array:
    for _atom in _entity.steps:
        if _atom.has_angles:
            print(_atom.chi)
        else:
            pass
