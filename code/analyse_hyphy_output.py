import os
import json
import pandas as pd
import multiprocessing as mp
from tqdm import tqdm

# to suppress warning from ete3 because it's not up to date with python3.12
import warnings

warnings.filterwarnings("ignore", category=SyntaxWarning)
from ete3 import Tree

# Description: Helper functions to analyse the output of HyPhy analyses.

"""
Helper functions to analyse the output of HyPhy analyses.
"""


def extract_orthogroup_id(filepath_str):
    # Extract the orthogroup ID from the filename
    return os.path.basename(filepath_str).split(".")[0]


def read_json_file(args_tuple):
    """
    Read a JSON file and extract the orthogroup ID and specified keys from the JSON data.
    The arguments tuple should contain the following elements:
    - json_filepath_str: The path to the JSON file to read.
    - keys_to_extract_list: A list of keys to extract from the JSON data.
    - subkeys_to_exclude_dict: A dictionary with a list of subkeys to exclude for each key in `keys_to_extract_list`.
    """
    json_filepath_str, keys_to_extract_list, subkeys_to_exclude_dict = args_tuple
    try:
        # Open and read the JSON file
        with open(json_filepath_str, "r") as json_file:
            json_data_dict = json.load(json_file)
            # Extract the orthogroup ID and the specified keys from the JSON data
            extracted_data_dict = {}
            if keys_to_extract_list == "all":
                keys_to_extract_list = list(json_data_dict.keys())
            for key in keys_to_extract_list:
                if key in json_data_dict:
                    if key in subkeys_to_exclude_dict:
                        extracted_data_dict[key] = {
                            subkey: value
                            for subkey, value in json_data_dict[key].items()
                            if subkey not in subkeys_to_exclude_dict[key]
                        }
                    else:
                        extracted_data_dict[key] = json_data_dict[key]
                else:
                    extracted_data_dict[key] = None
            return 0, extract_orthogroup_id(json_filepath_str), extracted_data_dict
    except json.JSONDecodeError:
        # print(f"Error reading {json_filepath_str}")
        return 1, extract_orthogroup_id(json_filepath_str), "Error reading"
    except Exception as e:
        # print(f"Error reading {json_filepath_str}: {e}")
        return 1, extract_orthogroup_id(json_filepath_str), f"Exception: {e}"


def extract_hyphy_output_relax(
    output_directory_str,
    keys_to_extract_list="all",
    subkeys_to_exclude_dict={},
    max_threads=None,
    debug=False,
):
    # List all JSON files in the output directory
    json_filepaths_list = [
        os.path.join(output_directory_str, filename)
        for filename in os.listdir(output_directory_str)
        if filename.endswith(".json")
    ]
    if debug:
        print("DEBUG: Extracting only first 10 files")
        json_filepaths_list = json_filepaths_list[:10]
    hyphy_results_json_list = []

    if max_threads is None:
        max_threads = mp.cpu_count() - 2

    # Use a multiprocessing pool to read JSON files in parallel
    OGs_Error_reading = []
    OGs_Exception = {}
    with mp.Pool(max_threads) as pool:
        for pcode, orthogroup_id_str, extracted_data_dict in tqdm(
            pool.imap_unordered(
                read_json_file,
                [
                    (json_filepath_str, keys_to_extract_list, subkeys_to_exclude_dict)
                    for json_filepath_str in json_filepaths_list
                ],
            ),
            total=len(json_filepaths_list),
            desc="Extracting data",
        ):
            # Store the extracted data in the results dictionary
            if pcode == 0:
                hyphy_results_json_list.append(
                    {
                        "OG": orthogroup_id_str,
                        **extracted_data_dict,
                    }
                )
            elif pcode == 1:
                if extracted_data_dict == "Error reading":
                    OGs_Error_reading.append(orthogroup_id_str)
                elif extracted_data_dict.startswith("Exception"):
                    OGs_Exception[orthogroup_id_str] = extracted_data_dict
            else:
                raise Exception(f"Unknown error code: {pcode}")

    # for each key (OG), we want the following subkeys:
    # "fits": {"MG94xREV with separate rates for branch sets":{"AIC-c", "Log Likelihood", "Rate Distributions", "estimated parameters"},
    # "Unconstrained model": {"AIC-c", "Log Likelihood", "estimated parameters",
    # "Rate Distributions": {"Background": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}},
    #                       "Test": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}}}}
    # "Constrained model": {"AIC-c", "Log Likelihood", "estimated parameters",
    # "Rate Distributions": {"Background": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}},
    #                       "Test": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}}}}
    # }
    # "test results": {"p-value", "LRT"}
    # (Note that Constrained model may not exist for all OGs in which case we will have to skip it)
    hyphy_results_json_list = [
        {k: v for k, v in d.items() if k in ["OG", "fits", "test results"]}
        for d in hyphy_results_json_list
    ]

    # Take all of these keys and flatten them into single keys for each OG. E.g. "Unconstrained model|AIC-c", "Unconstrained model|Rate Distributions|Background|0|omega", etc.

    # Create a DataFrame from the extracted data and flatten the DataFrame
    hyphy_results_flat_df = pd.json_normalize(hyphy_results_json_list, sep="|")

    # choose only the columns that we want
    columns_to_keep = [
        "OG",
        # "fits|MG94xREV with separate rates for branch sets|AIC-c",
        # "fits|MG94xREV with separate rates for branch sets|Log Likelihood",
        # "fits|MG94xREV with separate rates for branch sets|Rate Distributions",
        # "fits|MG94xREV with separate rates for branch sets|estimated parameters",
        "fits|Unconstrained model|AIC-c",
        "fits|Unconstrained model|Log Likelihood",
        "fits|Unconstrained model|estimated parameters",
        "fits|Unconstrained model|Rate Distributions|Background|0|omega",
        "fits|Unconstrained model|Rate Distributions|Background|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Background|1|omega",
        "fits|Unconstrained model|Rate Distributions|Background|1|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|0|omega",
        "fits|Unconstrained model|Rate Distributions|Test|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|1|omega",
        "fits|Unconstrained model|Rate Distributions|Test|1|proportion",
        "fits|Constrained model|AIC-c",
        "fits|Constrained model|Log Likelihood",
        "fits|Constrained model|estimated parameters",
        "fits|Constrained model|Rate Distributions|Background|0|omega",
        "fits|Constrained model|Rate Distributions|Background|0|proportion",
        "fits|Constrained model|Rate Distributions|Background|1|omega",
        "fits|Constrained model|Rate Distributions|Background|1|proportion",
        "fits|Constrained model|Rate Distributions|Test|0|omega",
        "fits|Constrained model|Rate Distributions|Test|0|proportion",
        "fits|Constrained model|Rate Distributions|Test|1|omega",
        "fits|Constrained model|Rate Distributions|Test|1|proportion",
        "test results|p-value",
        "test results|LRT",
        "test results|relaxation or intensification parameter",
    ]
    # remove columns that are not in the list above. Some elements of the list may not be present in the columns
    columns_to_keep = [
        col for col in columns_to_keep if col in hyphy_results_flat_df.columns
    ]
    hyphy_results_flat_df = hyphy_results_flat_df[columns_to_keep]

    rename_col_dict = {
        "fits|MG94xREV with separate rates for branch sets|AIC-c": "MG94xREV|AIC-c",
        "fits|MG94xREV with separate rates for branch sets|Log Likelihood": "MG94xREV|Log Likelihood",
        "fits|MG94xREV with separate rates for branch sets|Rate Distributions": "MG94xREV|Rate Distributions",
        "fits|MG94xREV with separate rates for branch sets|estimated parameters": "MG94xREV|estimated parameters",
        "fits|Unconstrained model|AIC-c": "Unconstrained model|AIC-c",
        "fits|Unconstrained model|Log Likelihood": "Unconstrained model|Log Likelihood",
        "fits|Unconstrained model|estimated parameters": "Unconstrained model|estimated parameters",
        "fits|Unconstrained model|Rate Distributions|Background|0|omega": "Unconstrained model|Background|0|omega",
        "fits|Unconstrained model|Rate Distributions|Background|0|proportion": "Unconstrained model|Background|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Background|1|omega": "Unconstrained model|Background|1|omega",
        "fits|Unconstrained model|Rate Distributions|Background|1|proportion": "Unconstrained model|Background|1|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|0|omega": "Unconstrained model|Test|0|omega",
        "fits|Unconstrained model|Rate Distributions|Test|0|proportion": "Unconstrained model|Test|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|1|omega": "Unconstrained model|Test|1|omega",
        "fits|Unconstrained model|Rate Distributions|Test|1|proportion": "Unconstrained model|Test|1|proportion",
        "fits|Constrained model|AIC-c": "Constrained model|AIC-c",
        "fits|Constrained model|Log Likelihood": "Constrained model|Log Likelihood",
        "fits|Constrained model|estimated parameters": "Constrained model|estimated parameters",
        "fits|Constrained model|Rate Distributions|Background|0|omega": "Constrained model|Background|0|omega",
        "fits|Constrained model|Rate Distributions|Background|0|proportion": "Constrained model|Background|0|proportion",
        "fits|Constrained model|Rate Distributions|Background|1|omega": "Constrained model|Background|1|omega",
        "fits|Constrained model|Rate Distributions|Background|1|proportion": "Constrained model|Background|1|proportion",
        "fits|Constrained model|Rate Distributions|Test|0|omega": "Constrained model|Test|0|omega",
        "fits|Constrained model|Rate Distributions|Test|0|proportion": "Constrained model|Test|0|proportion",
        "fits|Constrained model|Rate Distributions|Test|1|omega": "Constrained model|Test|1|omega",
        "fits|Constrained model|Rate Distributions|Test|1|proportion": "Constrained model|Test|1|proportion",
        "test results|p-value": "RELAX test|p-value",
        "test results|LRT": "RELAX test|LRT",
        "test results|relaxation or intensification parameter": "RELAX test|K",
    }
    hyphy_results_flat_df = hyphy_results_flat_df.rename(columns=rename_col_dict)

    print(
        f"Error while reading JSON for {len(OGs_Error_reading)} OGs: {OGs_Error_reading}"
    )
    print(
        f"Exceptions while reading JSON for {len(OGs_Exception)} OGs: {OGs_Exception}"
    )

    return hyphy_results_flat_df


def extract_hyphy_output(
    output_directory_str,
    keys_to_extract_list="all",
    subkeys_to_exclude_dict={},
    max_threads=None,
    debug=False,
):
    # List all JSON files in the output directory
    json_filepaths_list = [
        os.path.join(output_directory_str, filename)
        for filename in os.listdir(output_directory_str)
        if filename.endswith(".json")
    ]
    if debug:
        print("DEBUG: Extracting only first 10 files")
        json_filepaths_list = json_filepaths_list[:10]
    hyphy_results_json_list = []

    if max_threads is None:
        max_threads = mp.cpu_count() - 2

    # Use a multiprocessing pool to read JSON files in parallel
    OGs_Error_reading = []
    OGs_Exception = {}
    with mp.Pool(max_threads) as pool:
        for pcode, orthogroup_id_str, extracted_data_dict in tqdm(
            pool.imap_unordered(
                read_json_file,
                [
                    (json_filepath_str, keys_to_extract_list, subkeys_to_exclude_dict)
                    for json_filepath_str in json_filepaths_list
                ],
            ),
            total=len(json_filepaths_list),
            desc="Extracting data",
        ):
            # Store the extracted data in the results dictionary
            if pcode == 0:
                hyphy_results_json_list.append(
                    {
                        "OG": orthogroup_id_str,
                        **extracted_data_dict,
                    }
                )
            elif pcode == 1:
                if extracted_data_dict == "Error reading":
                    OGs_Error_reading.append(orthogroup_id_str)
                elif extracted_data_dict.startswith("Exception"):
                    OGs_Exception[orthogroup_id_str] = extracted_data_dict
            else:
                raise Exception(f"Unknown error code: {pcode}")

    # for each key (OG), we want the following subkeys:
    # "fits": {"MG94xREV with separate rates for branch sets":{"AIC-c", "Log Likelihood", "Rate Distributions", "estimated parameters"},
    # "Unconstrained model": {"AIC-c", "Log Likelihood", "estimated parameters",
    # "Rate Distributions": {"Background": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}},
    #                       "Test": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}}}}
    # "Constrained model": {"AIC-c", "Log Likelihood", "estimated parameters",
    # "Rate Distributions": {"Background": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}},
    #                       "Test": {"0": {"omega", "proportion"}, "1": {"omega", "proportion"}}}}
    # }
    # "test results": {"p-value", "LRT"}
    # (Note that Constrained model may not exist for all OGs in which case we will have to skip it)
    hyphy_results_json_list = [
        {k: v for k, v in d.items() if k in ["OG", "fits", "test results"]}
        for d in hyphy_results_json_list
    ]

    # Take all of these keys and flatten them into single keys for each OG. E.g. "Unconstrained model|AIC-c", "Unconstrained model|Rate Distributions|Background|0|omega", etc.

    # Create a DataFrame from the extracted data and flatten the DataFrame
    hyphy_results_flat_df = pd.json_normalize(hyphy_results_json_list, sep="|")

    # choose only the columns that we want
    columns_to_keep = [
        "OG",
        # "fits|MG94xREV with separate rates for branch sets|AIC-c",
        # "fits|MG94xREV with separate rates for branch sets|Log Likelihood",
        # "fits|MG94xREV with separate rates for branch sets|Rate Distributions",
        # "fits|MG94xREV with separate rates for branch sets|estimated parameters",
        "fits|Unconstrained model|AIC-c",
        "fits|Unconstrained model|Log Likelihood",
        "fits|Unconstrained model|estimated parameters",
        "fits|Unconstrained model|Rate Distributions|Background|0|omega",
        "fits|Unconstrained model|Rate Distributions|Background|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Background|1|omega",
        "fits|Unconstrained model|Rate Distributions|Background|1|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|0|omega",
        "fits|Unconstrained model|Rate Distributions|Test|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|1|omega",
        "fits|Unconstrained model|Rate Distributions|Test|1|proportion",
        "fits|Constrained model|AIC-c",
        "fits|Constrained model|Log Likelihood",
        "fits|Constrained model|estimated parameters",
        "fits|Constrained model|Rate Distributions|Background|0|omega",
        "fits|Constrained model|Rate Distributions|Background|0|proportion",
        "fits|Constrained model|Rate Distributions|Background|1|omega",
        "fits|Constrained model|Rate Distributions|Background|1|proportion",
        "fits|Constrained model|Rate Distributions|Test|0|omega",
        "fits|Constrained model|Rate Distributions|Test|0|proportion",
        "fits|Constrained model|Rate Distributions|Test|1|omega",
        "fits|Constrained model|Rate Distributions|Test|1|proportion",
        "test results|p-value",
        "test results|LRT",
    ]
    # remove columns that are not in the list above. Some elements of the list may not be present in the columns
    columns_to_keep = [
        col for col in columns_to_keep if col in hyphy_results_flat_df.columns
    ]
    hyphy_results_flat_df = hyphy_results_flat_df[columns_to_keep]

    rename_col_dict = {
        "fits|MG94xREV with separate rates for branch sets|AIC-c": "MG94xREV|AIC-c",
        "fits|MG94xREV with separate rates for branch sets|Log Likelihood": "MG94xREV|Log Likelihood",
        "fits|MG94xREV with separate rates for branch sets|Rate Distributions": "MG94xREV|Rate Distributions",
        "fits|MG94xREV with separate rates for branch sets|estimated parameters": "MG94xREV|estimated parameters",
        "fits|Unconstrained model|AIC-c": "Unconstrained model|AIC-c",
        "fits|Unconstrained model|Log Likelihood": "Unconstrained model|Log Likelihood",
        "fits|Unconstrained model|estimated parameters": "Unconstrained model|estimated parameters",
        "fits|Unconstrained model|Rate Distributions|Background|0|omega": "Unconstrained model|Background|0|omega",
        "fits|Unconstrained model|Rate Distributions|Background|0|proportion": "Unconstrained model|Background|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Background|1|omega": "Unconstrained model|Background|1|omega",
        "fits|Unconstrained model|Rate Distributions|Background|1|proportion": "Unconstrained model|Background|1|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|0|omega": "Unconstrained model|Test|0|omega",
        "fits|Unconstrained model|Rate Distributions|Test|0|proportion": "Unconstrained model|Test|0|proportion",
        "fits|Unconstrained model|Rate Distributions|Test|1|omega": "Unconstrained model|Test|1|omega",
        "fits|Unconstrained model|Rate Distributions|Test|1|proportion": "Unconstrained model|Test|1|proportion",
        "fits|Constrained model|AIC-c": "Constrained model|AIC-c",
        "fits|Constrained model|Log Likelihood": "Constrained model|Log Likelihood",
        "fits|Constrained model|estimated parameters": "Constrained model|estimated parameters",
        "fits|Constrained model|Rate Distributions|Background|0|omega": "Constrained model|Background|0|omega",
        "fits|Constrained model|Rate Distributions|Background|0|proportion": "Constrained model|Background|0|proportion",
        "fits|Constrained model|Rate Distributions|Background|1|omega": "Constrained model|Background|1|omega",
        "fits|Constrained model|Rate Distributions|Background|1|proportion": "Constrained model|Background|1|proportion",
        "fits|Constrained model|Rate Distributions|Test|0|omega": "Constrained model|Test|0|omega",
        "fits|Constrained model|Rate Distributions|Test|0|proportion": "Constrained model|Test|0|proportion",
        "fits|Constrained model|Rate Distributions|Test|1|omega": "Constrained model|Test|1|omega",
        "fits|Constrained model|Rate Distributions|Test|1|proportion": "Constrained model|Test|1|proportion",
        "test results|p-value": "BUSTED-SR test|p-value",
        "test results|LRT": "BUSTED-SR test|LRT",
    }
    hyphy_results_flat_df = hyphy_results_flat_df.rename(columns=rename_col_dict)

    print(
        f"Error while reading JSON for {len(OGs_Error_reading)} OGs: {OGs_Error_reading}"
    )
    print(
        f"Exceptions while reading JSON for {len(OGs_Exception)} OGs: {OGs_Exception}"
    )

    return hyphy_results_flat_df


def extract_ogs_w_wo_eg(
    ecotype_varwise_branchwise_file,
    gene_varwise_branchwise_file,
    input_genome_tree_file,
):

    # first find the recipient_branch values that have ecotype gains and the ones without
    ecotype_varwise_branchwise_df = pd.read_csv(
        ecotype_varwise_branchwise_file, header=0, sep="\t"
    )
    gene_varwise_branchwise_df = pd.read_csv(
        gene_varwise_branchwise_file, header=0, sep="\t"
    )

    # # if it's not recipient branch, rename col with gene_ or eco_ prefix
    # ecotype_varwise_branchwise_df = ecotype_varwise_branchwise_df.rename(
    #     {
    #         k: f"eco_{k}"
    #         for k in ecotype_varwise_branchwise_df.columns
    #         if k != "recipient_branch"
    #     },
    #     axis=1,
    # )
    # gene_varwise_branchwise_df = gene_varwise_branchwise_df.rename(
    #     {
    #         k: f"gene_{k}"
    #         for k in gene_varwise_branchwise_df.columns
    #         if k != "recipient_branch"
    #     },
    #     axis=1,
    # )

    # groupby recipient_branch and list all nog_id values for each
    gene_varwise_branchwise_df = (
        gene_varwise_branchwise_df.groupby("recipient_branch")["nog_id"]
        .apply(set)
        .reset_index()
    )
    # same for ecotype but sum the transfers
    ecotype_varwise_branchwise_df = (
        ecotype_varwise_branchwise_df.groupby("recipient_branch")
        .sum()["transfers"]
        .reset_index()
    )
    # merge the two on the recipient branch column, and fill NaN with 0
    combined_gene_eco_df = pd.merge(
        ecotype_varwise_branchwise_df,
        gene_varwise_branchwise_df,
        on="recipient_branch",
        how="outer",
    ).fillna(0)

    # # # debug
    # print(combined_gene_eco_df)

    # read in all branch names from the input genome tree using ete3
    tree = Tree(input_genome_tree_file, format=1)
    branches_list = [node.name for node in tree.traverse()]

    # get the recipient_branch values that have ecotype gains
    w_eg_branches = combined_gene_eco_df[combined_gene_eco_df["transfers"] > 0][
        "recipient_branch"
    ]
    wo_eg_branches = combined_gene_eco_df[combined_gene_eco_df["transfers"] == 0][
        "recipient_branch"
    ]

    w_eg_nog_ids_dict_set = (
        combined_gene_eco_df[combined_gene_eco_df["transfers"] > 0]
        .set_index("recipient_branch")["nog_id"]
        .to_dict()
    )
    wo_eg_nog_ids_dict_set = (
        combined_gene_eco_df[combined_gene_eco_df["transfers"] == 0]
        .set_index("recipient_branch")["nog_id"]
        .to_dict()
    )

    # combine the sets

    # set of nog_ids for branches with ecotype gains
    w_eg_nog_ids_set = set(
        [
            nog_id
            for nog_id_set in w_eg_nog_ids_dict_set.values()
            for nog_id in nog_id_set
        ]
    )

    # set of nog_ids for branches without ecotype gains
    wo_eg_nog_ids_set = set(
        [
            nog_id
            for nog_id_set in wo_eg_nog_ids_dict_set.values()
            for nog_id in nog_id_set
        ]
    )

    return w_eg_nog_ids_set, wo_eg_nog_ids_set
