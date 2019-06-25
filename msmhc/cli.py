from . import __version__ 

from argparse import ArgumentParser
from sys import argv

from .fasta import write_fasta
from collections import OrderedDict

from varcode.reference import genome_for_reference_name
parser = ArgumentParser("MS-MHC")

parser.add_argument("--output", required=True, help="Name of output FASTA file")
parser.add_argument("--reference", default="grch37", help="Name of reference genome (e.g. 'grch37')")

parser.add_argument(
    "--upstream-reading-frames",
    default=False,
    action="store_true",
    help="Include upstream reading frames from 5' UTR start codons")

parser.add_argument(
    "--downstream-reading-frames",
    default=False,
    action="store_true",
    help="Include downstream reading frames from use of alternative start codons")

parser.add_argument(
    "--skip-exons",
    default=False,
    action="store_true",
    help="Include sequences by skipping exons in a coding sequence")

def generate_sequences(genome):
    sequence_to_names = OrderedDict()

    for t in genome.transcripts():
        if t.is_protein_coding and t.complete and t.protein_sequence is not None:
            if t.protein_sequence in sequence_to_names:
                sequence_to_names[t.protein_sequence].add(t.transcript_name)
            else:
                sequence_to_names[t.protein_sequence] = {t.transcript_name}
    names_to_sequences = OrderedDict()
    for sequence, transcript_names in sequence_to_names.items():
        key = ";".join(transcript_names)
        names_to_sequences[key] = sequence
    return names_to_sequences

def run(args_list=None):
    if args_list is None:
        args_list = argv[1:]
    args = parser.parse_args(args_list)
    print("MS-MHC version %s" % __version__)
    reference_genome = genome_for_reference_name(args.reference)
    name_to_seq_dict = generate_sequences(reference_genome)
    print("Writing %d records" % len(name_to_seq_dict))
    write_fasta(name_to_seq_dict)
    print("Done.")

