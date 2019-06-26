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

import progressbar

from .reference_sequence import ReferenceSequence
from .mutant_sequence import MutantSequence

def generate_reference_sequences(genome):
    """
    Generate list of ReferenceSequence objects which may
    repeat the same protein sequence.
    """
    sequences = []
    for t in progressbar.progressbar(genome.transcripts()):
        if t.biotype == "protein_coding" and t.complete and t.protein_sequence is not None:
            sequences.append(ReferenceSequence(t))
    return sequences

def generate_mutant_sequences(variants):
    sequences = []
    effects = variants.effects()
    for effect in effects.top_priority_effect_per_variant().values():
        if effect.modifies_protein_sequence:
            if effect.mutant_protein_sequence is not None:
                sequences.append(MutantSequence(effect))
    return sequences

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
    sequences = generate_reference_sequences(genome)
    if variants:
        sequences.extend(generate_mutant_sequences(variants))
    if upstream_reading_frames:
        pass
    if downstream_reading_frames:
        pass
    if skip_exons:
        pass
    return sequences