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


from . import __version__

from argparse import ArgumentParser
from sys import argv

from .main import generate_sequences

from varcode.reference import genome_for_reference_name

from varcode.cli import variant_collection_from_args, add_variant_args

parser = ArgumentParser("MS-MHC")

parser.add_argument("--output", required=True, help="Name of output FASTA file")

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

add_variant_args(parser)

def run(args_list=None):
    if args_list is None:
        args_list = argv[1:]
    args = parser.parse_args(args_list)
    print("MS-MHC version %s" % __version__)
    reference_genome = genome_for_reference_name(args.genome)
    if args.vcf or args.maf or args.variant or args.json_variants:
        variants = variant_collection_from_args(args)
    else:
        variants = []
    sequences = generate_sequences(
        genome=reference_genome,
        variants=variants,
        upstream_reading_frames=args.upstream_reading_frames,
        downstream_reading_frames=args.downstream_reading_frames,
        skip_exons=args.skip_exons)
    print("Writing %d records" % len(sequences))
    with open(args.outputs, "w") as f:
        for seq in sequences:
            seq.write_to_fasta_file(f)
    print("Done.")

