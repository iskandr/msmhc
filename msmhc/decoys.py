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
from random import shuffle, seed

from .sequence import Sequence

class Decoy(Sequence):
    decoy_counter = 0

    def __init__(self, amino_acids):
        Decoy.decoy_counter += 1
        name = "Decoy-%d" % (Decoy.decoy_counter,)
        Sequence.__init__(
            self,
            name=name,
            amino_acids=amino_acids,
            attributes={"source": "decoy"})


def generate_decoys(
        sequences,
        n_decoys=None,
        max_scrambling_attempts_per_decoy=3,
        random_seed=0):
    """
    Parameters
    ----------
    sequences : list of objects derived from Sequence

    n_decoys : int
        Number of decoys to generate, if None then make one for each hit.

    max_scrambling_attempts_per_decoy : int
        Number of times to try to scramble a sequence before giving up if all
        scrambled results were already in the hit sequences.

    random_seed : int
        Initial random seed to make scrambling determinstic


    Returns list of Decoy
    -------

    """
    if n_decoys <= 0:
        return []

    seed(random_seed)
    real_peptide_set = {s.amino_acids for s in sequences}
    real_peptide_list = list(real_peptide_set)
    shuffle(real_peptide_list)
    n_hits = len(real_peptide_list)
    decoys = []
    print("Generating decoy sequences by scrambling...")
    if n_decoys is None:
        n_decoys = len(real_peptide_list)

    outer_iter = 0
    while len(decoys) < n_decoys:
        outer_iter += 1
        expected_iters = n_decoys / n_hits

        if outer_iter > 100 * expected_iters:
            print("Warning: failed to generate sufficient decoys")

        for peptide in progressbar(real_peptide_list):
            if len(decoys) >= n_decoys:
                break
            list_of_amino_acids = list(peptide)
            for attempt in range(max_scrambling_attempts_per_decoy):
                shuffle(list_of_amino_acids)
                scrambled = "".join(list_of_amino_acids)
                if scrambled not in real_peptide_set:
                    decoys.append(Decoy(amino_acids=scrambled))
                    break

    print("Generated %d decoy sequences" % len(decoys))

    return decoys