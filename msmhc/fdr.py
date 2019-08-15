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

"""
Module for analyzing output from Byonic mass spec software and selecting
peptide hits up to a specified FDR rate threshold.
"""

import numpy as np
import pandas as pd
import progressbar

def keep_unique_peptides(df):
    """
    Group DataFrame by peptide sequence, keep peptide with highest
    log p-value and if there are ties keep one with highest score.

    Parameters
    ----------
    df : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame
    """
    rows = []
    for peptide, df_group in df.groupby("Peptide"):
        if len(df_group) > 1:
            probs = df_group["|Log Prob|"]
            max_prob = probs.max()
            df_group = df_group[probs == max_prob]
            scores = df_group["Score"]
            max_score = scores.max()
            df_group = df_group[scores == max_score]
        rows.append(df_group.iloc[0])
    return pd.DataFrame(rows)


def load_byonic_output(
        filename="011419_5e7ceq_Cervical_15501T2_Rep3_GRCH37upstreamdownstream.xlsx",
        sheet_name="Spectra",
        raw_peptide_column="Peptide\n< ProteinMetrics Confidential >",
        modification_column="Modification Type(s)"):
    """
    Load Byonic Excel file as Pandas DataFrame, and group by peptide sequence.
    Also creates the following extra columns:
        - Peptide : peptide sequence without protein context
        - Before : amino acid before peptide in full protein
        - After : amino acid after peptide in full protein
        - Source : what kind of sequence did each peptide come from?

    Parameters
    ----------
    filename : str

    sheet_name : str

    raw_peptide_column : str

    modification_column : str

    Returns
    -------
    pandas.DataFrame
    """
    df = pd.read_excel(
        filename,
        sheet_name=sheet_name)
    df["Peptide"] = df[raw_peptide_column].map(lambda x: x[2:-2])
    df["Before"] = df[raw_peptide_column].map(lambda x: x[0])
    df["After"] = df[raw_peptide_column].map(lambda x: x[-1])
    df["Source"] = df["Protein Name"].map(
        lambda x: x[1:].split(" ")[0].split("-")[0])
    df["Is_Modified"] = ~(df[modification_column].isnull())
    return keep_unique_peptides(df)

def fdr_curve(df, score_column="|Log Prob|"):
    """
    Estimate FDR for every unique score value.

    Parameters
    ----------
    df : pandas.DataFrame

    score_column : str

    Returns
    -------
    Three arrays:
        - score thresholds
        - number of samples kept at each threshold
        - FDR estimate for each threshold
    """
    unique_scores = np.array(list(reversed(sorted(set(df[score_column])))))
    counts = []
    fdrs = []
    for score in progressbar.progressbar(unique_scores):
        df_above = df[df[score_column] >= score]
        source_counts = df_above["Source"].value_counts()
        n_decoy = source_counts.get("Decoy", 0)
        n_other = source_counts.sum() - n_decoy
        fdr = n_decoy / n_other
        counts.append(n_other)
        fdrs.append(fdr)
    counts = np.array(counts)
    fdrs = np.array(fdrs)
    return unique_scores, counts, fdrs

