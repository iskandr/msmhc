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

from .mass import mass_of_peptide
from .common import convert_to_string

class Sequence(object):
    """
    Base class used for reference and modified protein sequences
    """
    __slots__ = [
        "name",
        "amino_acids",
        "attributes",
    ]

    def __init__(self, name, amino_acids, attributes={}):
        self.name = name
        self.amino_acids = amino_acids

        attributes = attributes.copy()
        attributes["length"] = len(amino_acids)
        attributes["mass"] = mass_of_peptide(amino_acids)
        self.attributes = attributes

    def __str__(self):
        return "%s(name='%s', amino_acids='%s', attributes=%s)" % (
            self.__class__.__name__,
            self.name,
            self.amino_acids,
            self.attributes)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if self.__class__ is not other.__class__:
            return False
        if self.name != other.name:
            return False
        if self.amino_acids != other.amino_acids:
            return False
        for k, v in self.attributes.items():
            if other.attributes.get(k) != v:
                return False
        return True

    def __hash__(self):
        return hash(self.name)

    def sorted_attribute_list(self):
        """
        Returns list of attribute key/value pairs in alphabetical order
        of the keys.
        """
        return sorted(self.attributes.items(), key=lambda x: x[0])

    def attribute_string(self):
        return " ".join([
            "%s=%s" % (k, convert_to_string(v))
            for (k, v)
            in self.sorted_attribute_list()])

    def sequence_split_into_lines(self, maxwidth=80):
        lines = []
        for i in range(len(self.amino_acids) // maxwidth + 1):
            subseq = self.amino_acids[i * maxwidth:(i + 1) * maxwidth]
            if len(subseq) > 0:
                lines.append(subseq)
        return lines

    def fasta_string(self):
        return ">%s %s\n%s\n" % (
            self.name,
            self.attribute_string(),
            "\n".join(self.sequence_split_into_lines()))

    def write_to_fasta_file(self, file_handle):
        file_handle.write(self.fasta_string())