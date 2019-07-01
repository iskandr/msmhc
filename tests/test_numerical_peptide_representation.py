from msmhc.numerical_peptide_representation import (
    peptide_to_int,
    int_to_peptide,
    drop_first_letter,
    drop_last_letter
)
from nose.tools import eq_

def test_SIINFEKL_round_trip():
    eq_(int_to_peptide(peptide_to_int("SIINFEKL")), "SIINFEKL")

def test_SIINFEKL_drop_first_letter():
    eq_(int_to_peptide(
            drop_first_letter(peptide_to_int("SIINFEKL"))), "IINFEKL")

def test_SIINFEKL_drop_last_letter():
    eq_(int_to_peptide(
            drop_last_letter(peptide_to_int("SIINFEKL"))), "SIINFEK")

def test_SIINFEKL_drop_multiple_letters():
    eq_(int_to_peptide(
            drop_last_letter(
                drop_first_letter(
                    drop_last_letter(peptide_to_int("SIINFEKL"))))), "IINFE")

def test_RRR_round_trip():
    eq_(int_to_peptide(peptide_to_int("RRR")), "RRR")

