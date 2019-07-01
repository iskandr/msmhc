# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from collections import Counter, defaultdict

from progressbar import progressbar

from .sequence import Sequence
from .mutant_sequence import MutantSequence
from .reference_sequence import ReferenceSequence
from .alt_orf import AltORF, UpstreamORF, DownstreamORF

class_priority_list = [
    ReferenceSequence,
    MutantSequence,
    UpstreamORF,
    DownstreamORF,
    AltORF,
    Sequence,
]

class_priority_dict = {t: i for (i, t) in enumerate(class_priority_list)}

def sequence_type_priority_key_function(sequence_type):
    """
    Returns higher number for higher priority sequence type,
    so that e.g. ReferenceSequence should have a higher number
    than AltORF.

    Parameters
    ----------
    sequence_type : type

    Returns
    -------
    int
    """
    if sequence_type in class_priority_dict:
        i = class_priority_dict[sequence_type]
        # negate the priority since 0 is the "highest", so we want
        # all lower priority types to have lower numerical values
        return -i
    raise ValueError("Unrecognzied type %s" % (sequence_type,))

def keep_max_priority_sequences(sequences):
    """
    Determine the highest priority type in the input list of sequences
    and keep only sequences which have the same type.
    Parameters
    ----------
    sequences : list of objects derived from Sequence

    Returns
    -------
    list of objects derived from Sequence, all of same type
    """
    sequences_types = map(type, sequences)
    highest_priority_type = max(
        sequences_types,
        key=sequence_type_priority_key_function)

    return [
        s for s in sequences
        if type(s) is highest_priority_type
    ]

def collapse_peptide_sources(peptide_dict):
    """
    Given a dictionary mapping from peptide sequences to all of the
    different protein sequences which generated that peptide, collapse
    the protein sequences into an aggregate Sequence object whose
    attribute dictionary maps field to sets of values.

    Parameters
    ----------
    peptide_dict : dict
        Dictionary from peptide to list of objects derived from Sequence

    Returns
    -------
    List of Sequence corresponding to unique peptides
    """
    name_group_counts = Counter()
    peptide_sequences = []
    for peptide, sources in peptide_dict.items():
        assert len(sources) > 0
        filtered_sources = keep_max_priority_sequences(sources)
        assert len(filtered_sources) > 0
        attribute_dicts = [
            s.attributes for s in filtered_sources
        ]
        combined_attributes = defaultdict(set)
        for d in attribute_dicts:
            for k, v in d.items():
                combined_attributes[k].add(v)
        type_name = filtered_sources[0].__class__.__name__

        if "gene_name" in combined_attributes:
            gene_names = combined_attributes["gene_name"]
            # names will look like ReferenceSequence-TP53
            name_group = "%s-%s" % (
                    type_name,
                    "-".join(sorted(gene_names)))
        else:
            name_group = type_name
        name_group_counts[name_group] += 1
        # add a numerical ID to make the sequence names unique
        name = "%s-%d" % (
            name_group,
            name_group_counts[name_group])
        peptide_sequences.append(Sequence(
            name=name,
            amino_acids=peptide,
            attributes=combined_attributes))
    return peptide_sequences

def extract_peptides(sequences, min_length=7, max_length=20):
    """
    Extract subsequences from full protein sequences

    Parameters
    ----------
    sequences : list of Sequence
        All generated protein sequences

    min_length : int
        Smallest peptide length to include

    max_length : int
        Largest peptide length to include

    Returns
    -------
    Dictionary from str to list of Sequence objects which contained that peptide
    """
    peptide_dict = defaultdict(list)
    for sequence_obj in progressbar(sequences):
        amino_acids = sequence_obj.amino_acids
        n_aa = len(amino_acids)
        for k in range(min_length, max_length + 1):
            if n_aa >= k:
                for i in range(n_aa - k + 1):
                    peptide_dict[amino_acids[i:i + k]].append(sequence_obj)
    return peptide_dict