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

class MutantSequence(Sequence):
    """
    Representation of protein sequences from reference transcripts
    """

    __slots__ = [
        "effect",
        "_sanitized_effect_description",
    ]

    def __init__(self, effect):
        self.effect = effect
        self._sanitized_effect_description = \
            self.effect.short_description.replace(" ", "-")
        Sequence.__init__(
            self,
            name="mut-%s-%s" % (
                self._sanitized_effect_description,
                effect.transcript.transcript_id),
            amino_acids=self.effect.mutant_protein_sequence,
            attributes={
                "source": "mutation",
                "genomic_variant": effect.variant.short_description,
                "protein_effect": self._sanitized_effect_description,
                "transcript_id": effect.transcript.transcript_id,
                "transcript_name": effect.transcript.transcript_name,
                "gene_id": effect.transcript.gene.gene_id,
                "gene_name": effect.transcript.gene.gene_name,
            })
