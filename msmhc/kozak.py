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

KOZAK_CONSENSUS_SEQUENCE_BEFORE_START = "gccRcc"
ALT_START_CODONS = {"CTG", "GTG", "TTG"}
KOZAK_CONSENSEUS_SEQUENCE_AFTER_START = "Gcg"


def translation_initiation_score(sequence_before_start, start_codon, sequence_after_start):
    if start_codon == "ATG":
        score = 100
    elif start_codon in ALT_START_CODONS:
        score = 10
    else:
        score = 0

    for pair in [
            zip(KOZAK_CONSENSUS_SEQUENCE_BEFORE_START[-len(sequence_before_start):], sequence_before_start),
            zip(KOZAK_CONSENSEUS_SEQUENCE_AFTER_START, sequence_after_start)]:
        for (consensus, observed) in pair:
            consensus_upper = consensus.upper()
            if consensus_upper == "R":
                valid = {"A", "G"}
            else:
                valid = {consensus_upper}
            match = observed in valid
            if consensus == consensus_upper:
                score += 10 * match
            else:
                score += match
    max_score = 100 + 5 + 10 + 10 + 2
    return score / max_score
