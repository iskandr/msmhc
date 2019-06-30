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
from .peptides import collapse_peptide_sources, extract_peptides

def generate_reference_sequences(
        genome,
        restrict_sources_to_gene_name=None):
    """
    Generate list of ReferenceSequence objects which may
    repeat the same protein sequence.

    Parameters
    ----------
    genome : pyensembl.Genome

    restrict_sources_to_gene_name : str or None

    Returns
    -------
    list of ReferenceTranscript
    """
    sequences = []
    print("Gathering transcripts...")
    if restrict_sources_to_gene_name:
        genes = genome.genes_by_name(restrict_sources_to_gene_name)
        transcripts = []
        for g in genes:
            transcripts.extend(g.transcripts)
    else:
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


def generate_protein_sequences(
        genome,
        variants=[],
        upstream_reading_frames=False,
        downstream_reading_frames=False,
        skip_exons=False,
        restrict_sources_to_gene_name=None):
    """

    Parameters
    ----------
    genome : pyensembl.Genome

    variants : varcode.VariantCollection

    upstream_reading_frames : bool

    downstream_reading_frames : bool

    skip_exons : bool

    restrict_sources_to_gene_name : str or None

    Returns list of msmhc.Sequence
    """
    print("Generating sequences from reference transcripts")
    reference_sequences = generate_reference_sequences(
        genome,
        restrict_sources_to_gene_name=restrict_sources_to_gene_name)
    sequences = reference_sequences.copy()
    if upstream_reading_frames:
        print("Generating sequences from upstream reading frames")
        sequences.extend(generate_upstream_reading_frames(reference_sequences))
    if downstream_reading_frames:
        print("Generating sequences from downstream reading frames")
        sequences.extend(generate_downstream_reading_frames(reference_sequences))
    if skip_exons:
        print("Generating sequences from skipped exons")
        sequences.extend(generate_skipped_exon_sequences(reference_sequences))
    if variants:
        print("Generating sequences from %d variants" % len(variants))
        sequences.extend(generate_mutant_sequences(variants))
    return sequences


def generate_peptide_sequences(
        genome,
        variants=[],
        upstream_reading_frames=False,
        downstream_reading_frames=False,
        skip_exons=False,
        restrict_sources_to_gene_name=None,
        min_peptide_length=7,
        max_peptide_length=20):
    """
    Parameters
    ----------
    genome : pyensembl.Genome

    variants : varcode.VariantCollection

    upstream_reading_frames : bool

    downstream_reading_frames : bool

    skip_exons : bool

    restrict_sources_to_gene_name : str or None

    min_peptide_length : int

    max_peptide_length : int

    Returns list of msmhc.Sequence containing peptide sequences
    """
    protein_sequences = generate_protein_sequences(
        genome=genome,
        variants=variants,
        upstream_reading_frames=upstream_reading_frames,
        downstream_reading_frames=downstream_reading_frames,
        skip_exons=skip_exons,
        restrict_sources_to_gene_name=restrict_sources_to_gene_name)
    peptide_sequence_dict = extract_peptides(
        protein_sequences,
        min_length=min_peptide_length,
        max_length=max_peptide_length)
    peptide_sequences = collapse_peptide_sources(peptide_sequence_dict)
    return peptide_sequences