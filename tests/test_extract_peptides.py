from msmhc.sequence import Sequence
from msmhc.peptides import extract_peptides
from nose.tools import eq_

def test_extract_peptides():
    seq = Sequence(name="test-seq", amino_acids="SIINFEKL")
    peptide_dict = extract_peptides([seq], min_length=7, max_length=8)
    eq_(set(peptide_dict.keys()), {
        "SIINFEKL",
        "SIINFEK",
        "IINFEKL"
    })