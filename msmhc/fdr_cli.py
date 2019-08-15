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

import argparse
import sys

from .fdr import load_byonic_output, fdr_curve

parser = argparse.ArgumentParser(
    "msmhc-fdr",
    description="Extract peptide hits at specified FDR cutoff from Byonic output")

input_group = parser.add_argument_group("Input")
input_group.add_argument(
    "--excel",
    required=True)

fdr_group = parser.add_argument_group("FDR")
fdr_group.add_argument(
    "--fdr-cutoff",
    default=0.01,
    type=float,
    help="Maximum desired false discovery rate")

fdr_group.add_argument(
    "--fdr-score-column",
    default="|Log Prob|",
    help="Which column of Byonic PSM data should be used as score for FDR determination")

output_group = parser.add_argument_group("--output")
output_group.add_argument("--output-csv", required=True)

def run(args_list=None):
    if args_list is None:
        args_list = sys.argv[1:]
    args = parser.parse_args(args_list)
    print("Loading %s..." % args.excel)
    df = load_byonic_output(args.excel)

    print("Peptide source counts before filtering:")
    print(df["Source"].value_counts())

    print("Computing FDR curve...")
    unique_scores, counts, fdrs = fdr_curve(
        df,
        score_column=args.fdr_score_column)
    mask = fdrs < args.fdr_cutoff
    fdr_estimate = fdrs[mask][-1]
    best_score = unique_scores[mask][-1]
    df_kept = df[df[args.fdr_score_column] >= best_score]
    print("Keeping %d peptide entries at %s=%0.2f (FDR ~= %0.4f)" % (
        len(df_kept),
        args.fdr_score_column,
        best_score,
        fdr_estimate))
    print("Peptide source counts before filtering:")
    print(df_kept["Source"].value_counts())
    df_kept.to_csv(args.output_csv, index=False)


