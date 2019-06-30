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

mass_table = """
# letter mono-isotropic-mass average-mass three-letters
A	71.037114	71.0779	Ala
R	156.101111	156.1857	Arg
N	114.042927	114.1026	Asn
D	115.026943	115.0874	Asp
C	103.009185	103.1429    Cys
E	129.042593	129.114	Glu
Q	128.058578	128.1292	Gln
G	57.021464	57.0513	Gly
H	137.058912	137.1393	His
I	113.084064	113.1576	Ile
L	113.084064	113.1576	Leu
K	128.094963	128.1723	Lys
M	131.040485	131.1961	Met
F	147.068414	147.1739	Phe
P	97.052764	97.1152	Pro
S	87.032028	87.0773	Ser
T	101.047679	101.1039	Thr
U	150.95363	150.0379	Sec
W	186.079313	186.2099	Trp
Y	163.06332	163.1733	Tyr
V	99.068414	99.1311	Val
"""

mass_dict = {}

for line in mass_table.split("\n"):
    if line and not line.startswith("#"):
        parts = line.split()
        letter = parts[0]
        average_mass = float(parts[2])
        mass_dict[letter] = average_mass

def mass_of_peptide(peptide):
    """
    Add up average monomeric masses for each residue in peptide

    Parameters
    ----------
    peptide : str

    Returns
    -------
    float
    """
    return sum([mass_dict[amino_acid] for amino_acid in peptide])



