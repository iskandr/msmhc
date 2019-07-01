from math import log, ceil

classical_amino_acids = "RHKDESTNQCGPAVILMFYW"
assert len(classical_amino_acids) == 20

amino_acids = classical_amino_acids + "U"

amino_acid_order = {aa: i for (i, aa) in enumerate(amino_acids)}

n_aa = len(amino_acids)
base = n_aa + 1

def peptide_to_int(sequence):
    result = 0
    for i, aa in enumerate(sequence):
        idx = amino_acid_order[aa]
        result += (base ** i) * (idx + 1)
    return result

def int_to_peptide(numerical_representation):
    letters = []
    while numerical_representation >= n_aa:
        current_idx = numerical_representation % base - 1
        numerical_representation = numerical_representation // base
        letters.append(amino_acids[current_idx])
    letters.append(amino_acids[numerical_representation - 1])
    return "".join(letters)


def drop_first_letter(numerical_representation):
    return numerical_representation // base

def peptide_length(numerical_representation):
    if numerical_representation < base:
        return 1
    exponent = log(numerical_representation, base)
    return int(ceil(exponent))

def drop_last_letter(numerical_representation):
    return numerical_representation % (base ** (peptide_length(numerical_representation) - 1))
