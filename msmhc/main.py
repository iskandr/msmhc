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

from collections import OrderedDict
from .sequence_from_reference_transcript import SequenceFromReferenceTranscript

def generate_redundant_reference_sequences(genome):
    """
    Generate list of SequenceFromReferenceTranscript objects which may
    repeat the same protein sequence.
    """
    sequences = []
    for t in genome.transcripts():
        if t.is_protein_coding and t.complete and t.protein_sequence is not None:
            sequences(SequenceFromReferenceTranscript(t))
    return sequences

def collpase_redundant_sequences(sequences):
    str_to_objects = OrderedDict()
    for obj in sequences:
        if obj in str_to_objects

def generate_sequences(
        genome,
        variants=[],
        upstream_reading_frames=False,
        downstream_reading_frames=False,
        skip_exons=False):
    sequences = generate_redundant_reference_sequences(genome)
    if variants:
        pass
    if upstream_reading_frames:
        pass
    if downstream_reading_frames:
        pass
    if skip_exons:
        pass
    return sequences