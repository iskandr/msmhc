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

class Sequence(object):
    """
    Base class used for reference and modified protein sequences
    """
    def __init__(self, name, amino_acids, attributes={}):
        self.name = name
        self.amino_acids = amino_acids

        attributes = attributes.copy()
        attributes["length"] = len(amino_acids)
        self.attributes = attributes

    def fasta_key(self):
        return self.name

    def sorted_attribute_list(self):
        """
        Returns list of attribute key/value pairs in alphabetical order
        of the keys.
        """
        return sorted(self.attributes.items(), key=lambda x: x[0])

    def attribute_string(self):
        return " ".join([
            "%s=%s" % (k, v)
            for (k, v)
            in self.sorted_attribute_list()])

    def sequence_split_into_lines(self, maxwidth=60):
        lines = []
        for i in range(len(self.amino_acids) // maxwidth + 1):
            subseq = self.amino_acids[i * maxwidth:(i + 1) * maxwidth]
            if len(subseq) > 0:
                lines.append(subseq)
        return lines

    def fasta_string(self):
        return ">%s %s\n%s\n" % (
            self.fasta_key(),
            self.fasta_attribute_string(),
            self.sequence_split_into_lines())

    def write_to_fasta_file(self, file_handle):
        file_handle.write(self.fasta_string())