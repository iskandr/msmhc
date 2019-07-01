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

from progressbar import progressbar

from .alt_orf import generate_alt_reading_frames
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
    print("Gathering reference genes...")
    if restrict_sources_to_gene_name:
        genes = genome.genes_by_name(restrict_sources_to_gene_name)
    else:
        genes = genome.genes()

    for g in progressbar(genes):
        if g.is_protein_coding:
            for t in g.transcripts:
                if t.is_protein_coding and t.complete and t.protein_sequence is not None:
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


def generate_upstream_reading_frames(
        sequences,
        min_peptide_length=7):
    """

    Parameters
    ----------
    sequences : list of ReferenceSequence

    min_peptide_length : int

    Returns
    -------
    list of AltORF
    """
    results = []
    for sequence in progressbar(sequences):
        results.extend(
            generate_alt_reading_frames(
                sequence,
                min_peptide_length=min_peptide_length,
                search_start_offset=None,
                search_end_offset=0))
    print("Generated %d upstream reading frames" % len(results))
    return results


def generate_downstream_reading_frames(
        sequences,
        min_peptide_length=7,
        search_end_offset=500):
    """
    Parameters
    ----------
    sequences : list of ReferenceSequence

    min_peptide_length : int

    search_end_offset : int
        Maximum nucleotides past the original start codon to look for new ORFs.
        Setting the default to 500bp since the paper
        "miniMAVS, you complete me!" managed to find an alternative start
        site 400bp downstream from the annotated start of MAVS.

    Returns
    -------
    list of DownstreamORF
    """
    results = []
    for sequence in progressbar(sequences):
        results.extend(
            generate_alt_reading_frames(
                sequence,
                min_peptide_length=min_peptide_length,
                search_start_offset=3,
                search_end_offset=search_end_offset))
    print("Generated %d downstream reading frames" % len(results))
    return results

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
        restrict_sources_to_gene_name=None,
        min_peptide_length=7):
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

    Returns list of msmhc.Sequence
    """
    print("Generating sequences from reference transcripts")
    reference_sequences = generate_reference_sequences(
        genome,
        restrict_sources_to_gene_name=restrict_sources_to_gene_name)
    sequences = reference_sequences.copy()
    if upstream_reading_frames:
        print("Generating sequences from upstream reading frames")
        sequences.extend(
            generate_upstream_reading_frames(
                reference_sequences,
                min_peptide_length=min_peptide_length))
    if downstream_reading_frames:
        print("Generating sequences from downstream reading frames")
        sequences.extend(
            generate_downstream_reading_frames(
                reference_sequences,
                min_peptide_length=min_peptide_length))
    if skip_exons:
        print("Generating sequences from skipped exons")
        sequences.extend(generate_skipped_exon_sequences(reference_sequences))
    if variants:
        print("Generating sequences from %d variants" % len(variants))
        sequences.extend(generate_mutant_sequences(variants))
    return sequences


