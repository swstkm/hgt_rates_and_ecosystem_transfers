# This file contains functions that are used by the notebooks in the notebooks folder.

import dis
import ete3
import numpy as np
import pandas as pd
import os
from IPython.display import display
import traceback


def map_output_to_input_nodes(input_tree, output_tree, ale=False):
    """
    Given two ete3 trees, this function returns a dictionary that maps the nodes of the output tree to the nodes of the input tree.
    """
    # first we create a dictionary for each tree. The keys are the descendant leaves of each internal node, and the values are the nodes themselves.
    # The keys are strings of the form "leaf1,leaf2,leaf3" where the leaf names are sorted lexicographically.
    input_tree_dict = {}
    for input_node in input_tree.traverse(strategy="postorder"):
        if not input_node.is_leaf():  # internal node
            descendants = input_node.get_leaf_names()
            descendants.sort()
            input_tree_dict[",".join(descendants)] = input_node.name
            if input_node.name == "":
                print(
                    "Empty name found. Number of descendants: {}.".format(
                        len(descendants)
                    )
                )
        else:
            input_tree_dict[input_node.name] = input_node.name
    output_tree_dict = {}
    for output_node in output_tree.traverse(
        strategy="postorder"
    ):  # we traverse the tree in postorder
        if not output_node.is_leaf():  # internal node
            descendants = output_node.get_leaf_names()
            descendants.sort()
            if ale:
                output_tree_dict[",".join(descendants)] = f"N{output_node.name}"
            else:
                output_tree_dict[",".join(descendants)] = output_node.name
        else:
            output_tree_dict[output_node.name] = output_node.name

    # now we iterate over the keys of the input tree dictionary, and for each key we check if it is also a key in the output tree dictionary
    # if it is, we add an entry to the mapping dictionary (output_node -> input_node)
    # if not, we add an entry to the mapping dictionary (input_node -> None), and raise an error
    node_mapping = {}
    node_not_found_list = []
    for output_key in output_tree_dict.keys():
        if (
            output_key in input_tree_dict
        ):  # if the output key is found in the input tree, add the mapping
            node_mapping[output_tree_dict[output_key]] = input_tree_dict[output_key]
        else:
            node_mapping[output_tree_dict[output_key]] = None
            node_not_found_list.append(output_key)
    # if node_not_found_list is not empty, print length and contents of it
    if node_not_found_list:
        print(f"Nodes not found in input tree: {len(node_not_found_list)}")
        print(node_not_found_list)

    return node_mapping


def compile_ale_outputs(
    output_dir: str, input_tree: str, compiled_results_dir: str, var_str: str
):
    """
    Given the output_dir of an ALE run, and the input_tree, this function compiles the outputs of the ALE run into a set of TSV files.
    """

    # extract the specific ALE dir name from the output_dir
    var_str = f"{var_str}.ale"

    # read in the first uml file you can find in output_dir.
    uml_file = [f for f in os.listdir(output_dir) if f.endswith(".uml_rec")][0]
    # Split the third line and take the last element as the newick string
    with open(os.path.join(output_dir, uml_file), "r") as uml_fo:
        ale_tree_string = uml_fo.readlines()[2].split()[-1].strip()

    ale_tree = ete3.Tree(ale_tree_string, format=1)
    ale_node_mapping = map_output_to_input_nodes(input_tree, ale_tree, ale=True)

    # now read in all the *uTs files in the output_dir, and process them.
    # each file contains the columns: 'source_branch', 'recipient_branch', and 'freq' (of transfer).
    # for each tuple of branches (from, to),  take the sum of the freqs across all the files

    # first, get the list of all the .uTs files (not .uml_rec)
    all_uts_files = [f for f in os.listdir(output_dir) if f.endswith(".uTs")]

    # first we create a dict of dfs, mapping nog to file
    varwise_hgt_dfs_dict = {}
    for uts_file in all_uts_files:
        nog_id = uts_file.split("_")[-1].split(".")[0]
        # read in the file
        try:
            with open(os.path.join(output_dir, uts_file), "r") as uts_fo:
                # split each line and take last 3 elements, join the first 2 of them for the from-to tuple, and create a dict between that and the last element
                uts_lines = [
                    tuple(l.split()[-3:])
                    for l in uts_fo.readlines()
                    if not l.startswith("#")
                ]
                new_uts_lines = []
                for i, uts_line in enumerate(uts_lines):
                    from_branch, to_branch, freq = uts_line
                    if "(" in from_branch:
                        new_from_branch = from_branch.split("(")[0]
                    else:
                        new_from_branch = f"N{from_branch}"
                    if "(" in to_branch:
                        new_to_branch = to_branch.split("(")[0]
                    else:
                        new_to_branch = f"N{to_branch}"
                    new_uts_lines.append((new_from_branch, new_to_branch, freq))
                uts_records = [
                    (
                        ale_node_mapping[from_branch],
                        ale_node_mapping[to_branch],
                        float(freq),
                    )
                    for from_branch, to_branch, freq in new_uts_lines
                ]
                # create a df from this list of records, with columns 'source_branch'(str), 'recipient_branch'(str), 'transfers'(float)
                uts_df = pd.DataFrame.from_records(
                    uts_records,
                    columns=["source_branch", "recipient_branch", "transfers"],
                )
                # make sure the from and to columns are string type only and transfers is float type
                uts_df["source_branch"] = uts_df["source_branch"].astype(str)
                uts_df["recipient_branch"] = uts_df["recipient_branch"].astype(str)
                uts_df["transfers"] = uts_df["transfers"].astype(float)
                # add nog_id column
                uts_df["nog_id"] = nog_id
                # add this df to the dict
                varwise_hgt_dfs_dict[nog_id] = uts_df
        except Exception as e:
            print(f"Error reading in {uts_file}")
            raise e

    # then we can create a varwise.branchwise file by concatenating the dfs.
    varwise_branchwise_hgt_df = pd.concat(
        varwise_hgt_dfs_dict.values(), ignore_index=True
    )

    # columns should be 'nog_id', 'source_branch', 'recipient_branch', 'transfers' in that order
    varwise_branchwise_hgt_df = varwise_branchwise_hgt_df[
        ["nog_id", "source_branch", "recipient_branch", "transfers"]
    ]
    # show it
    print("varwise, branchwise transfers:")
    print(varwise_branchwise_hgt_df)

    # # debug: display original ALE node IDs using ale_node_mapping for branch in rows of this df where source==recipient
    # reverse_ale_node_mapping = {v: k for k, v in ale_node_mapping.items()}
    # print("Original ALE node IDs for branches where source==recipient:")
    # debug_varwise_branchwise_hgt_df = varwise_branchwise_hgt_df.copy()
    # debug_varwise_branchwise_hgt_df["ale_branch_id"] = debug_varwise_branchwise_hgt_df[
    #     "source_branch"
    # ].map(reverse_ale_node_mapping)
    # debug_varwise_branchwise_hgt_df = debug_varwise_branchwise_hgt_df[
    #     debug_varwise_branchwise_hgt_df["source_branch"]
    #     == debug_varwise_branchwise_hgt_df["recipient_branch"]
    # ]
    # display(debug_varwise_branchwise_hgt_df)

    # write it out
    varwise_branchwise_hgt_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.varwise.branchpairwise.{var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # create a varwise transfers file by summing the transfers for each nog_id for a df with columns 'nog_id', 'transfers'. Drop other columns
    varwise_hgt_df = (
        varwise_branchwise_hgt_df.groupby("nog_id")["transfers"].sum().reset_index()
    )
    # show it
    print("varwise transfers:")
    print(varwise_hgt_df)
    # write it out
    varwise_hgt_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.varwise.{var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # create a branch-pair-wise transfers file by taking mean of the transfers for each branch pair
    branchpairwise_hgt_df = (
        varwise_branchwise_hgt_df.groupby(["source_branch", "recipient_branch"])[
            "transfers"
        ]
        .mean()
        .reset_index()
    )
    # show it
    print("Branchwise transfers for each pair of branches:")
    print(branchpairwise_hgt_df)
    # write it out
    branchpairwise_hgt_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.branchpairwise.{var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # create a branchwise transfers file by taking the sum of every recipient branch across all source branches
    rec_branchwise_hgt_df = (
        branchpairwise_hgt_df.groupby("recipient_branch")["transfers"].sum().reset_index()
    )
    # rename this recipient branch column to 'branch'
    rec_branchwise_hgt_df.rename(columns={"recipient_branch": "branch"}, inplace=True)
    # show it and write it out
    print("Recipient branchwise transfers:")
    print(rec_branchwise_hgt_df)
    rec_branchwise_hgt_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.branchwise.{
                                 var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )
    # similarly for varwise.branchwise df, create a recipient branchwise df by summing the transfers for each recipient branch
    rec_varwise_branchwise_hgt_df = (
        varwise_branchwise_hgt_df.groupby(["nog_id", "recipient_branch"])["transfers"]
        .sum()
        .reset_index()
    )
    # rename this recipient branch column to 'branch'
    rec_varwise_branchwise_hgt_df.rename(
        columns={"recipient_branch": "branch"}, inplace=True
    )
    # show it and write it out
    print("Recipient varwise branchwise transfers:")
    print(rec_varwise_branchwise_hgt_df)
    rec_varwise_branchwise_hgt_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.varwise.branchwise.{
                                         var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    ##########################################################################################################

    # now we compile the losses inferred by ALE. This needs to be read in from the .uml_rec files
    all_uml_rec_files = [f for f in os.listdir(output_dir) if f.endswith(".uml_rec")]

    varwise_losses_dfs_dict = {}
    for uml_rec_file in all_uml_rec_files:
        nog_id = uml_rec_file.split("_")[2].split(".")[0]
        with open(os.path.join(output_dir, uml_rec_file), "r") as uml_rec_fo:
            # read in the table after the line that contains "Originations"
            uml_rec_lines = uml_rec_fo.readlines()
            header_line = [i for i, l in enumerate(uml_rec_lines) if "Originations" in l][0]
            uml_rec_lines = uml_rec_lines[header_line + 1 :]
            # loop over lines and extract the branch name (with ale_node_mapping), and the number of losses
            uml_rec_records = []
            for l in uml_rec_lines:
                if l.strip() == "":
                    break
                l_split = l.split()
                branch_type, branch, losses = l_split[0], l_split[1], l_split[4]
                if branch_type == "S_terminal_branch":
                    new_branch = branch.split("(")[0]
                else:
                    new_branch = f"N{branch}"
                uml_rec_records.append((new_branch, float(losses)))
            # create a df from this list of records, with columns 'branch'(str), 'losses'(float)
            uml_rec_df = pd.DataFrame.from_records(
                uml_rec_records, columns=["branch", "losses"]
            )
            # make sure the branch column is string type and losses is float type
            uml_rec_df["branch"] = uml_rec_df["branch"].astype(str)
            uml_rec_df["losses"] = uml_rec_df["losses"].astype(float)
            # add nog_id column
            uml_rec_df["nog_id"] = nog_id
            # add this df to the dict
            varwise_losses_dfs_dict[nog_id] = uml_rec_df

    # then we can create a varwise.branchwise file by concatenating the dfs.
    # create a varwise branchwise losses df by concatenating the dfs
    varwise_branchwise_losses_df = pd.concat(
        varwise_losses_dfs_dict.values(), ignore_index=True
    )
    # show it
    print("varwise, branchwise losses:")
    print(varwise_branchwise_losses_df)
    # write it out
    varwise_branchwise_losses_df.to_csv(
        f"{compiled_results_dir}/compiled_losses.varwise.branchwise.{var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # group by branch and sum the losses to get branchwise losses
    branchwise_losses_df = (
        varwise_branchwise_losses_df.groupby("branch")["losses"].sum().reset_index()
    )
    # show it
    print("branchwise losses:")
    print(branchwise_losses_df)
    # write it out
    branchwise_losses_df.to_csv(
        f"{compiled_results_dir}/compiled_losses.branchwise.{var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # combine the branchwise transfers and losses to get branchwise dynamics
    branchwise_dynamics_df = pd.merge(
        rec_branchwise_hgt_df, branchwise_losses_df, on="branch", how="outer"
    )
    # fill NaNs with 0
    branchwise_dynamics_df = branchwise_dynamics_df.fillna(0)
    # show it
    print("branchwise dynamics:")
    print(branchwise_dynamics_df)
    # write it out
    branchwise_dynamics_df.to_csv(
        f"{compiled_results_dir}/compiled_dynamics.branchwise.{var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # similarly for varwise.branchwise dynamics, group by var_id and branch, and sum the transfers and losses
    varwise_branchwise_dynamics_df = (
        pd.merge(
            rec_varwise_branchwise_hgt_df,
            varwise_branchwise_losses_df,
            on=["nog_id", "branch"],
            how="outer",
        )
        .groupby(["nog_id", "branch"])[["transfers", "losses"]]
        .sum()
        .reset_index()
    )
    # fill NaNs on transfers and losses with 0
    varwise_branchwise_dynamics_df[["transfers", "losses"]] = varwise_branchwise_dynamics_df[
        ["transfers", "losses"]
    ].fillna(0)
    # show it
    print("varwise, branchwise dynamics:")
    print(varwise_branchwise_dynamics_df)
    # show any rows with nan values
    nan_rows = varwise_branchwise_dynamics_df[
        varwise_branchwise_dynamics_df.isnull().any(axis=1)
    ]
    if not nan_rows.empty:
        print("Rows with NaN values:")
        print(nan_rows)
    else:
        print("No NaN values found in varwise, branchwise dynamics")
    # write it out
    varwise_branchwise_dynamics_df.to_csv(
        f"{compiled_results_dir}/compiled_dynamics.varwise.branchwise.{var_str}.tsv",
        index=False,
        header=True,
        sep="\t",
    )


def compile_gloome_results(
    expectations_df: pd.DataFrame,
    gainloss_df: pd.DataFrame,
    gloome_node_mapping: dict,
    ml_mp: str,
    var_str: str,
    var_name_str: str,
    pa_matrix_tsv_filepath: str,
    compiled_results_dir: str,
    prob_cutoff: float
):

    print(f"Compiling GLOOME results for {ml_mp}.{var_str} for {var_name_str} values")

    if ml_mp not in ["ml", "mp"]:
        raise ValueError(
            "ml_mp must be either 'ml' or 'mp', for maximum-likelihood or maximum-parsimony"
        )

    # prepare varwise file
    # Replace the branch names with the input tree branch names.
    gainloss_df["gloome_branch_name"] = gainloss_df["branch"]
    gainloss_df["recipient_branch"] = gainloss_df["branch"].map(gloome_node_mapping)

    # we need to replace the POS column data with var-IDs.
    # First read in the tsv file of the PA matrix we created for Count (since the matrix was the same for GLOOME also)

    pa_matrix_df = pd.read_csv(
        os.path.join(pa_matrix_tsv_filepath), sep="\t", header=0, index_col=0
    )
    # the first column here contains the vars. The row number of the var corresponds to the POS IDs in the gainloss_df
    # create a dict of row number to var IDs in pa_matrix_df
    print(
        f"Var IDs in the PA matrix file for {ml_mp}.{
          var_str} looks like: {pa_matrix_df.index}"
    )
    pos_var_dict = {i + 1: var for i, var in enumerate(pa_matrix_df.index)}
    print(
        f"POS to var ID mapping for {ml_mp}.{
          var_str} looks like: {pos_var_dict}"
    )
    # now use the dict to replace POS column with var IDs
    gainloss_df["POS"] = gainloss_df["POS"].map(pos_var_dict)
    # rename the POS column to var_name_str and 'expectation' column to 'transfers'
    gainloss_df = gainloss_df.rename(
        columns={"POS": var_name_str, "expectation": "transfers"}
    )

    # if ml_mp is 'ml', filter the gains_df by the prob_cutoff on probability column
    if ml_mp == "ml" and prob_cutoff is not None:
        gainloss_df = gainloss_df[gainloss_df["probability"] >= prob_cutoff]

    # retain only the columns 'recipient_branch', 'transfers', var_name_str, 'gloome_branch_name' and 'G/L' in that order
    column_names = [
        var_name_str,
        "recipient_branch",
        "transfers",
        "gloome_branch_name",
        "G/L",
    ]
    gainloss_df = gainloss_df[column_names]
    gains_df = gainloss_df[gainloss_df["G/L"] == "gain"][column_names[:-1]]

    # groupby var_name_str and sum the transfers to get varwise transfers
    varwise_gains_df = gains_df.groupby(var_name_str).sum().reset_index()
    varwise_gains_df = varwise_gains_df[[var_name_str, "transfers"]]
    # show it
    print(f"varwise transfers for {ml_mp}.{var_str} looks like:")
    print(varwise_gains_df)

    # write it out
    varwise_gains_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.varwise.{var_str}.{ml_mp}.gloome.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # write out varwise.branchwise file
    # write out gainloss_df where G/L is gain.
    gains_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.varwise.branchwise.{var_str}.{ml_mp}.gloome.tsv",
        index=False,
        header=True,
        sep="\t",
    )
    print(
        f"varwise, branchwise transfers for {
        ml_mp}.{var_str} looks like:"
    )
    print(gains_df)
    # do the same for a losses_df
    losses_df = gainloss_df[gainloss_df["G/L"] == "loss"][column_names[:-1]]
    # rename transfers to losses
    losses_df = losses_df.rename(
        columns={"recipient_branch": "branch", "transfers": "losses"}
    )
    # show it
    print(f"varwise, branchwise losses for {ml_mp}.{var_str} looks like:")
    print(losses_df)
    # write it out
    losses_df.to_csv(
        f"{compiled_results_dir}/compiled_losses.varwise.branchwise.{var_str}.{ml_mp}.gloome.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # prepare branchwise file
    # rename the columns 'branch' and 'exp01' to 'gloome_branch_name' and 'transfers' respectively, and drop the other columns
    # use the gloome_node_mapping dict to create an additional column 'branch' with the internal node names from the input tree
    expectations_df = expectations_df.rename(
        columns={"branch": "gloome_branch_name", "exp01": "transfers"}
    )
    expectations_df = expectations_df[["gloome_branch_name", "transfers"]]
    expectations_df["branch"] = expectations_df["gloome_branch_name"].map(
        gloome_node_mapping
    )
    # reorder columns
    expectations_df = expectations_df[["branch", "transfers", "gloome_branch_name"]]
    print(
        f"branchwise transfers df for {ml_mp}.{
        var_str} looks like:"
    )
    print(expectations_df)
    # write to file. Write 'NA' for missing data
    expectations_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.branchwise.{var_str}.{ml_mp}.gloome.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # write a branchwise file with gains and losses in separate columns, using the gainloss_df
    # the column 'G/L' has the values 'gain' and 'loss'. Some of the rows have 'loss' in 'G/L' column, and some have 'gain', and they can overlap in terms of 'branch' value
    # we need to create a df with columns 'recipient_branch', 'transfers', 'losses'
    # rename back to 'branch' from 'recipient_branch' for gains_df
    gains_df = gains_df.rename(columns={"recipient_branch": "branch"})
    # create a column combining 'branch' and 'var_name_str' for both gains and losses
    gains_df[f"{var_name_str}_B-branch"] = gains_df.apply(
        lambda x: f"{x[var_name_str]}_B-{x['branch']}", axis=1
    )
    losses_df[f"{var_name_str}_B-branch"] = losses_df.apply(
        lambda x: f"{x[var_name_str]}_B-{x['branch']}", axis=1
    )
    varwise_branchwise_changes_df = pd.merge(
        gains_df,
        losses_df,
        on=[f"{var_name_str}_B-branch", var_name_str],
        how="outer",
        suffixes=("_gain", "_loss"),
    )
    # fill NaNs with 0 but only for columns transfers and losses
    varwise_branchwise_changes_df = varwise_branchwise_changes_df.fillna(
        {col: 0 for col in varwise_branchwise_changes_df.columns if ("transfers" in col) or ("losses" in col)}
    )
    # define branch column using the _B-branch column
    varwise_branchwise_changes_df["branch"] = varwise_branchwise_changes_df[
        f"{var_name_str}_B-branch"
    ].apply(lambda x: x.split("_B-")[1])

    #  drop the column 'var_name_str_B-branch'
    varwise_branchwise_changes_df = varwise_branchwise_changes_df.drop(
        columns=[f"{var_name_str}_B-branch", "branch_loss", "branch_gain"]
    )
    # show it
    print(f"varwise, branchwise changes for {ml_mp}.{var_str} looks like:")
    print(varwise_branchwise_changes_df)
    # write to file
    varwise_branchwise_changes_df.to_csv(
        f"{compiled_results_dir}/compiled_dynamics.varwise.branchwise.{
                                         var_str}.{ml_mp}.gloome.tsv",
        index=False,
        header=True,
        sep="\t",
    )
    # group it by branch and sum across var_name_str to get branchwise changes
    branchwise_changes_df = (
        varwise_branchwise_changes_df.groupby("branch")[["transfers", "losses"]]
        .sum()
        .reset_index()
    )
    print(f"branchwise changes for {ml_mp}.{var_str} looks like:")
    print(branchwise_changes_df)
    # write to file
    branchwise_changes_df.to_csv(
        f"{compiled_results_dir}/compiled_dynamics.branchwise.{var_str}.{ml_mp}.gloome.tsv",
        index=False,
        header=True,
        sep="\t",
    )


def read_and_compile_gloome_results(
    gloome_output_dir: str,
    input_tree: str,
    var_str,
    var_name_str,
    pa_matrix_tsv_filepath,
    compiled_results_dir: str,
    prob_cutoff: float = 0.5,
):
    """
    Given the output directory of a GLOOME run, this function reads in the output files and compiles the results into a set of TSV files.
    Args:
    - gloome_output_dir: the output directory of the GLOOME run
    - input_tree: the input genome tree file with a newick format tree
    - var_str: a string that represents the variable type (e.g. 'gene', or 'ecosystem'). This will be used in the output file names.
    - var_name_str: a string that represents the variable name (e.g. 'nog_id', or 'ecosystem_id'). This will be used in the column names of the output files.
    - pa_matrix_tsv_filepath: the path to the PA matrix TSV file. Typically I use the same PA matrix file for Count and GLOOME, numerical for former and binary for the latter.

    """
    # if species_tree_bool is True, there was a tree given to GLOOME as input,
    # so we can find branchwise transfers and compare it to other programs like ALE, Ranger, etc.
    # But if it's not available, branchwise transfers are not comparable, so we don't need to write out those files

    gloome_tree = ete3.Tree(
        os.path.join(gloome_output_dir, "TheTree.INodes.ph"), format=1
    )
    gloome_node_mapping = map_output_to_input_nodes(input_tree, gloome_tree)

    # first we process the ML results, and then the MP results
    # for each of them, we read and prepare expectations_df and probabilties_df
    # case ML: branchwise expectations are from "ExpectationPerBranch.txt" and probabilities are from "gainLossProbExpPerPosPerBranch.txt"
    # case MP: branchwise expectations are from "gainLossMP.1.PerBranch.txt" and probabilities are from "gainLossMP.1.PerPosPerBranch.txt"

    # case ML: branchwise expectations
    ml_expectations_file_path = os.path.join(
        gloome_output_dir, "ExpectationPerBranch.txt"
    )
    ml_expectations_df = pd.read_csv(ml_expectations_file_path, skiprows=[0], sep="\t")
    # print(f"original ML expectations file looks like:")
    # print(ml_expectations_df)
    # case ML: varwise and varwise.branchwise expectations
    ml_gainloss_file_path = os.path.join(
        gloome_output_dir, "gainLossProbExpPerPosPerBranch.txt"
    )
    ml_gainloss_df = pd.read_csv(ml_gainloss_file_path, skiprows=[0], sep="\t")
    # print("original ML gainloss file looks like:")
    # print(gainloss_df)

    # prepare all compiled files for ML
    compile_gloome_results(
        ml_expectations_df,
        ml_gainloss_df,
        gloome_node_mapping,
        "ml",
        var_str,
        var_name_str,
        pa_matrix_tsv_filepath,
        compiled_results_dir,
        prob_cutoff
    )

    # case MP: branchwise expectations
    mp_expectations_file_path = [
        os.path.join(gloome_output_dir, f)
        for f in os.listdir(gloome_output_dir)
        if f.startswith("gainLossMP.") and f.endswith(".PerBranch.txt")
    ][0]
    # this file has 6 rows to skip instead of just one
    mp_expectations_df = pd.read_csv(
        mp_expectations_file_path, skiprows=list(range(6)), sep="\t"
    )

    # case MP: varwise and varwise.branchwise expectations
    mp_gainloss_file_path = [
        os.path.join(gloome_output_dir, f)
        for f in os.listdir(gloome_output_dir)
        if f.startswith("gainLossMP.") and f.endswith(".PerPosPerBranch.txt")
    ][0]
    # skip first 5 rows and read in the file
    mp_gainloss_df = pd.read_csv(
        mp_gainloss_file_path, skiprows=list(range(5)), sep="\t"
    )

    # prepare all compiled files for MP
    compile_gloome_results(
        mp_expectations_df,
        mp_gainloss_df,
        gloome_node_mapping,
        "mp",
        var_str,
        var_name_str,
        pa_matrix_tsv_filepath,
        compiled_results_dir,
        None,
    )


def find_changes_varwise_rowwise(row: pd.Series, input_tree: ete3.Tree):
    # Initialize the list to store the results
    dynamics_list = []

    # Iterate over the branches
    for branch in input_tree.get_descendants(strategy="postorder"):
        # Get the number of genes in the branch and the ancestor
        genes_in_branch = row[branch.name]
        genes_in_ancestor = row[branch.up.name]

        # Calculate the transfers, losses, expansions, and reductions
        transfers, losses, expansions, reductions = 0, 0, 0, 0
        if genes_in_ancestor == 0 and genes_in_branch > 0:
            transfers = genes_in_branch
        elif genes_in_ancestor > 0 and genes_in_branch == 0:
            losses = genes_in_ancestor
        elif genes_in_ancestor > 0 and genes_in_branch > genes_in_ancestor:
            expansions = genes_in_branch - genes_in_ancestor
        elif genes_in_ancestor > 0 and genes_in_branch < genes_in_ancestor:
            reductions = genes_in_ancestor - genes_in_branch

        # Add the results to the list
        dynamics_list.append(
            (row["name"], branch.name, transfers, losses, expansions, reductions)
        )

    return dynamics_list


def compile_count_changes(
    count_dir: str,
    input_tree_filepath: str,
    var_str: str,
    var_name_str: str,
    compiled_results_dir: str,
):

    # read in the count file
    count_branchwise_transfers_filepath = os.path.join(count_dir, "Count_changes.tsv")

    count_branchwise_transfers_df = pd.read_csv(
        count_branchwise_transfers_filepath,
        sep="\t",
        usecols=[1, 2],
        names=["branch", "transfers"],
        header=0,
    )
    # remove the node that says 'total'
    count_branchwise_transfers_df = count_branchwise_transfers_df[
        count_branchwise_transfers_df["branch"] != "total"
    ]
    # print and write to csv
    print("Count branchwise transfers:")
    print(count_branchwise_transfers_df)
    count_branchwise_transfers_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.branchwise.{var_str}.count.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # family/varwise transfers
    count_varwise_transfers_filepath = os.path.join(count_dir, "Count_families.tsv")
    # this file has headers like name, branch1, branch2, ..., Gains, Losses, Expansions, Reductions
    count_full_varwise_changes_df = pd.read_csv(
        count_varwise_transfers_filepath, sep="\t"
    )
    count_full_varwise_changes_df = count_full_varwise_changes_df.sort_values("name")
    # keep only the columns 'name', 'Gains'
    count_varwise_transfers_df = count_full_varwise_changes_df[
        ["name", "Gains"]
    ].rename(columns={"name": var_name_str, "Gains": "transfers"})
    # print and write to csv
    print("Count varwise transfers:")
    print(count_varwise_transfers_df)
    count_varwise_transfers_df.to_csv(
        f"{compiled_results_dir}/compiled_transfers.varwise.{var_str}.count.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # compiling varwise.branchwise file is somewhat tricky, since we don't directly get info of varwise transfers per branch
    # but in the varwise transfers output file (_families.tsv), we have the number of genes in each var in each branch, including internal ones
    # we can use this info to calculate the transfers per branch for each var
    # basically, we traverse the input tree, and for each var, for each branch, the total number of gains is the genes in the branch of this var, minus the genes in the ancestor.

    # first read in input tree
    input_tree = ete3.Tree(input_tree_filepath, format=1)

    # apply find_changes_varwise_rowwise to each row of the count_full_varwise_changes_df
    # to get the list of dynamics for each var
    varwise_branchwise_dynamics_list = count_full_varwise_changes_df.apply(
        find_changes_varwise_rowwise, axis=1, input_tree=input_tree
    )
    # flatten the list of lists
    varwise_branchwise_dynamics_list = [
        item for sublist in varwise_branchwise_dynamics_list for item in sublist
    ]

    # create a df from this list
    varwise_branchwise_changes_df = pd.DataFrame(
        varwise_branchwise_dynamics_list,
        columns=[
            var_name_str,
            "branch",
            "transfers",
            "losses",
            "expansions",
            "reductions",
        ],
    )
    # drop rows where all the dynamics are zero
    varwise_branchwise_changes_df = varwise_branchwise_changes_df[
        (
            varwise_branchwise_changes_df[
                ["transfers", "losses", "expansions", "reductions"]
            ]
            != 0
        ).any(axis=1)
    ]

    # print and write to csv
    print("Count varwise, branchwise transfers:")
    print(varwise_branchwise_changes_df)

    varwise_branchwise_changes_df.to_csv(
        f"{compiled_results_dir}/compiled_dynamics.varwise.branchwise.{var_str}.count.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # another file for just the transfers and not the losses, expansions, reductions
    varwise_branchwise_changes_df[[var_name_str, "branch", "transfers"]].to_csv(
        f"{compiled_results_dir}/compiled_transfers.varwise.branchwise.{var_str}.count.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # create a varwise df by summing the dynamics for each var
    varwise_changes_df = (
        varwise_branchwise_changes_df.groupby(var_name_str)[
            ["transfers", "losses", "expansions", "reductions"]
        ]
        .sum()
        .reset_index()
    )
    # sort this by var_id
    varwise_changes_df = varwise_changes_df.sort_values(var_name_str)
    # print and write to csv
    print("Count varwise dynamics:")
    print(varwise_changes_df)
    varwise_changes_df.to_csv(
        f"{compiled_results_dir}/compiled_dynamics.varwise.{var_str}.count.tsv",
        index=False,
        header=True,
        sep="\t",
    )

    # check if the full branchwise transfers matches what we have in the branchwise transfers df
    print(varwise_changes_df.head())
    print(
        count_full_varwise_changes_df[
            ["name", "Gains", "Losses", "Expansions", "Reductions"]
        ].head()
    )
    # show all the rows of varwise_changes_df where any of the changes don't match with what is there in count_full_varwise_changes_df
    # changes here refer to the columns 'transfers', 'losses', 'expansions', 'reductions'
    print(
        "Rows where the full branchwise transfers don't match with the compiled branchwise transfers:"
    )
    print(
        varwise_changes_df[
            ~varwise_changes_df.apply(
                lambda row: row["transfers"]
                == count_full_varwise_changes_df.loc[
                    count_full_varwise_changes_df["name"] == row[var_name_str], "Gains"
                ].values[0]
                or row["losses"]
                == count_full_varwise_changes_df.loc[
                    count_full_varwise_changes_df["name"] == row[var_name_str], "Losses"
                ].values[0]
                or row["expansions"]
                == count_full_varwise_changes_df.loc[
                    count_full_varwise_changes_df["name"] == row[var_name_str],
                    "Expansions",
                ].values[0]
                or row["reductions"]
                == count_full_varwise_changes_df.loc[
                    count_full_varwise_changes_df["name"] == row[var_name_str],
                    "Reductions",
                ].values[0],
                axis=1,
            )
        ]
    )
    print(
        "\033[91mIf the above df is not empty, then the calculation of changes that we make for varwise.branchwise is only approximate for the above rows.\033[0m"
    )

    branchwise_changes_df = (
        varwise_branchwise_changes_df.groupby("branch")[
            ["transfers", "losses", "expansions", "reductions"]
        ]
        .sum()
        .reset_index()
    )
    # print and write to csv
    print("Count branchwise dynamics:")
    print(branchwise_changes_df)
    branchwise_changes_df.to_csv(
        f"{compiled_results_dir}/compiled_dynamics.branchwise.{var_str}.count.tsv",
        index=False,
        header=True,
        sep="\t",
    )
