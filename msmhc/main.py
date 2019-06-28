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

from collections import defaultdict

import progressbar

from .reference_sequence import ReferenceSequence
from .mutant_sequence import MutantSequence

def generate_reference_sequences(genome):
    """
    Generate list of ReferenceSequence objects which may
    repeat the same protein sequence.

    Parameters
    ----------
    genome : pyensembl.Genome

    Returns
    -------
    list of ReferenceTranscript
    """
    sequences = []
    print("Gathering transcripts...")
    transcripts = genome.transcripts()
    for t in progressbar.progressbar(transcripts):
        if t.biotype == "protein_coding" and t.complete and t.protein_sequence is not None:
            sequences.append(ReferenceSequence(t))
    return sequences

def generate_mutant_sequences(variants):
    """
    Parameters
    ----------
    variants : varcode.VariantCollection

    Returns
    -------
    list of MutantSequence
    """
    sequences = []
    effects = variants.effects()
    for effect in effects.top_priority_effect_per_variant().values():
        if effect.modifies_protein_sequence:
            if effect.mutant_protein_sequence is not None:
                sequences.append(MutantSequence(effect))
    return sequences

def generate_upstream_reading_frames(sequences):
    """

    Parameters
    ----------
    sequences : list of ReferenceSequence

    Returns
    -------
    list of AltORF
    """
    pass


def generate_downstream_reading_frames(sequences):
    """
    Parameters
    ----------
    sequences : list of ReferenceSequence

    Returns
    -------
    list of AltORF
    """
    pass

def generate_skipped_exon_sequences(sequences):
    """
    Parameters
    ----------
    sequences : list of ReferenceSequence

    Returns
    -------
    list of ExonSkipSequence
    """
    pass

def extract_peptides(sequences, min_length=7, max_length=20):
    """
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
    for sequence_obj in sequences:
        amino_acids = sequence_obj.sequence
        n_aa = len(amino_acids)
        for k in range(min_length, max_length + 1):
            if n_aa >= k:
                for i in range(n_aa - k + 1):
                    peptide = amino_acids[i:i + k]
                    peptide_dict[peptide].append(sequence_obj)
    return peptide_dict


def generate_sequences(
        genome,
        variants=[],
        upstream_reading_frames=False,
        downstream_reading_frames=False,
        skip_exons=False):
    """

    Parameters
    ----------
    genome : pyensembl.Genome

    variants : varcode.VariantCollection

    upstream_reading_frames : bool

    downstream_reading_frames : bool

    skip_exons : bool

    Returns list of msmhc.Sequence
    """
    print("Generating sequences from reference transcripts")
    reference_sequences = generate_reference_sequences(genome)
    sequences = reference_sequences.copy()
    if upstream_reading_frames:
        print("Generating sequences from upstream reading frames")
        sequences.extend(generate_upstream_reading_frames(reference_sequences))
    if downstream_reading_frames:
        print("Generating sequences from downstream reading frames")
        sequences.extend(generate_downstream_reading_frames(reference_sequences))
    if skip_exons:
        print("Generating sequences from skipped exons")
        sequences.extedn(generate_skipped_exon_sequences(reference_sequences))
    if variants:
        print("Generating sequences from %d variants" % len(variants))
        sequences.extend(generate_mutant_sequences(variants))
    return sequences