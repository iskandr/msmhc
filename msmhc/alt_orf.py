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

class AltORF(Sequence):
    """
    Representation of upstream and downstream reading frames relative
    to annotated start codons
    """
    def __init__(self, transcript, relative_start, amino_acids, coding_sequence):
        self.transcript = transcript
        self.relative_start = relative_start
        self.downstream = relative_start > 0
        self.upstream = relative_start < 0
        self.coding_sequence = coding_sequence
        if self.upstream:
            name = "Upstream-%dnt-%s" % (
                -self.relative_start,
                self.transcript.id,
            )
        else:
            name = "Downstream-%dnt-%s" % (
                self.relative_start,
                self.transcript.id
            )
        Sequence.__init__(
            self,
            name=name,
            amino_acids=amino_acids,
            attributes={
                "upstream-orf": "1" if self.upstream else "0",
                "downstream-orf": "0" if self.upstream else "1",
                "relative_start": str(self.relative_start),
                "transcript_name": self.transcript.transcript_name,
                "transcript_id": self.transcript.transcript_id,
                "gene_name": self.transcript.gene_name,
                "gene_id": self.transcript.gene_id,
                "coding_sequence": self.coding_sequence})