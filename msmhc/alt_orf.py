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
    def __init__(self, transcript):
        self.transcript = transcript
        Sequence.__init__(
            self,
            name="alt-orf-?",
            amino_acids="?",
            attributes={})