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

from isovar.genetic_code import translate_cdna

from .sequence import Sequence
from .reference_sequence import ReferenceSequence
from .kozak import translation_initiation_score, ALT_START_CODONS

class AltORF(Sequence):
    """
    Representation of upstream and downstream reading frames relative
    to annotated start codons
    """
    __slots__ = [
        "original_sequence_name",
        "transcript_id",
        "transcript_name",
        "gene_id",
        "gene_name",
        "relative_start",
        "sequence_before_start",
        "start_codon",
        "sequence_after_start",
        "ends_with_stop_codon",
        "translation_initiation_score",
    ]
    def __init__(
            self,
            original_sequence_name,
            transcript_id,
            transcript_name,
            gene_id,
            gene_name,
            relative_start,
            amino_acids,
            sequence_before_start,
            start_codon,
            sequence_after_start,
            ends_with_stop_codon):
        assert len(amino_acids) > 0
        assert len(start_codon) == 3
        assert len(sequence_before_start) <= 6
        assert len(sequence_after_start) <= 3

        self.original_sequence_name = original_sequence_name
        self.transcript_id = transcript_id
        self.transcript_name = transcript_name
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.relative_start = relative_start


        self.sequence_before_start = sequence_before_start
        self.start_codon = start_codon
        self.sequence_after_start = sequence_after_start
        self.ends_with_stop_codon = ends_with_stop_codon
        self.translation_initiation_score = translation_initiation_score(
            self.sequence_before_start,
            self.start_codon,
            self.sequence_after_start)
        is_upstream = relative_start < 0

        if is_upstream:
            name = "Upstream-%s-%dnt" % (
                self.original_sequence_name,
                -self.relative_start,
            )
        else:
            name = "Downstream-%s-%dnt" % (
                self.original_sequence_name,
                self.relative_start,
            )
        Sequence.__init__(
            self,
            name=name,
            amino_acids=amino_acids,
            attributes={
                "source": "upstream-orf" if is_upstream else "downstream-orf",
                "relative_start": str(self.relative_start),
                "original_sequence_name": self.original_sequence_name,
                "gene_name": self.gene_name,
                "gene_id": self.gene_id,
                "transcript_name": self.transcript_name,
                "transcript_id": self.transcript_id,
                "sequence_before_start": self.sequence_before_start,
                "sequence_after_start": self.sequence_after_start,
                "start_codon": self.start_codon,
                "translation_initiation_score": self.translation_initiation_score,
                "ends_with_stop_codon": self.ends_with_stop_codon,
            })


class DownstreamORF(AltORF):
    pass

class UpstreamORF(AltORF):
    pass


START_CODONS = {"ATG"}.union(ALT_START_CODONS)

def generate_alt_reading_frames(
        sequence_obj,
        min_peptide_length=7,
        search_start_offset=None,
        search_end_offset=None,
        kozak_length_before_start_codon=6,
        kozak_length_after_start_codon=3):
    """
    Parameters
    ----------
    sequence_obj : ReferenceSequence

    min_peptide_length : int

    search_start_offset : int or None
        Offset relative to annotated start codon to start search,
        if None then start at beginning of sequence

    search_end_offset : int
        Offset relative to annotated start codon to end search, if None
        then go to end of original coding sequence

    kozak_length_before_start_codon : int

    kozak_length_after_start_codon : int

    Returns
    -------
    list of AltORF
    """
    if type(sequence_obj) != ReferenceSequence:
        return []
    results = []
    transcript = sequence_obj.transcript
    original_sequence_name = sequence_obj.name
    original_gene_name = transcript.gene_name
    original_gene_id = transcript.gene_id
    original_transcript_name = transcript.transcript_name
    original_transcript_id = transcript.transcript_id
    cdna_sequence = transcript.sequence
    start_codon_offset = min(transcript.start_codon_spliced_offsets)

    if search_start_offset is None:
        start = 0
    else:
        start = start_codon_offset + search_start_offset

    if search_end_offset is None:
        end = start_codon_offset + len(transcript.coding_sequence)
    else:
        end = start_codon_offset + search_end_offset

    for i in range(start, end):
        if i == 0:
            # skip the original reading frame
            continue
        first_codon = cdna_sequence[i:i + 3]
        relative_start = i - start_codon_offset
        if relative_start < 0:
            result_class = UpstreamORF
        else:
            result_class = DownstreamORF

        sequence_before_start_codon = cdna_sequence[
            max(0, i - kozak_length_before_start_codon):i]
        second_codon = cdna_sequence[i + 3:i + 3 + kozak_length_after_start_codon]
        if first_codon == "ATG":
            cds = cdna_sequence[i:]
            amino_acids, ends_with_stop_codon = translate_cdna(
                cdna_sequence=cds,
                first_codon_is_start=True)
            n_aa = len(amino_acids)
            if n_aa >= min_peptide_length:
                results.append(result_class(
                    relative_start=relative_start,
                    amino_acids=amino_acids,
                    original_sequence_name=original_sequence_name,
                    transcript_id=original_transcript_id,
                    transcript_name=original_transcript_name,
                    gene_id=original_gene_id,
                    gene_name=original_gene_name,
                    sequence_before_start=sequence_before_start_codon,
                    start_codon=first_codon,
                    sequence_after_start=second_codon,
                    ends_with_stop_codon=ends_with_stop_codon))
    return results