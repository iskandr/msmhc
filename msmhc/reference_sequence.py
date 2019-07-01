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

from .sequence import Sequence

class ReferenceSequence(Sequence):
    """
    Representation of protein sequences from reference transcripts
    """
    __slots__ = [
        "transcript",
    ]

    def __init__(self, transcript):
        self.transcript = transcript
        Sequence.__init__(
            self,
            name="ref-%s" % transcript.transcript_id,
            amino_acids=transcript.protein_sequence,
            attributes={
                "source": "reference",
                "transcript_id": transcript.transcript_id,
                "transcript_name": transcript.transcript_name,
                "gene_id": transcript.gene_id,
                "gene_name": transcript.gene_name,
            })
