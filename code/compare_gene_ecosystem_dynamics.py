from re import M
import pandas as pd
import numpy as np
import scipy.stats as st
import statsmodels.api as sm
import ete3

import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import iplot
import pandas as pd
import numpy as np
import ete3
import ipywidgets as widgets
from IPython.display import display

# this is to make sure that the plots are rendered in the notebook properly
import plotly.io as pio

pio.renderers.default = "notebook"


def linear_modelling(linearmodelling_combined_df: pd.DataFrame, formula: str) -> None:
    model = sm.formula.ols(formula=formula, data=linearmodelling_combined_df)
    res = model.fit()
    res_sum = res.summary()

    display(res_sum)
    return None


def combine_gene_env_inferences(
    gene_dynamics_input,
    env_dynamics_input,
    input_tree,
    single_copy_nogs_list_filepath=None,
):
    """
    Combine gene and environment inferences for each branch.
    """
    # Check if inputs are file paths or dataframes
    if isinstance(gene_dynamics_input, str):
        gene_dynamics_df = pd.read_csv(gene_dynamics_input, sep="\t", header=0)
    else:
        gene_dynamics_df = gene_dynamics_input

    if isinstance(env_dynamics_input, str):
        env_dynamics_df = pd.read_csv(env_dynamics_input, sep="\t", header=0)
    else:
        env_dynamics_df = env_dynamics_input

    # If single copy NOGs list is provided, filter the gene dynamics dataframe
    if single_copy_nogs_list_filepath:
        # this file may be multiple columns, but we only need the first column
        with open(single_copy_nogs_list_filepath) as f:
            single_copy_nogs_list = [line.strip().split()[0] for line in f.readlines()]
        gene_dynamics_df = gene_dynamics_df[
            gene_dynamics_df["nog_id"].isin(single_copy_nogs_list)
        ]

        # since single copy NOGs list was provided, the dynamics files should be varwise.branchwise.
        # convert them to branchwise by grouping by branch and summing the dynamics. Keep only 'branch', and 'transfers' and/or 'losses' columns
        gene_dynamics_df = gene_dynamics_df.groupby("branch").sum().reset_index()
        gene_dynamics_df = gene_dynamics_df[
            ["branch"]
            + [
                col
                for col in gene_dynamics_df.columns
                if col == "transfers" or col == "losses"
            ]
        ]
        env_dynamics_df = env_dynamics_df.groupby("branch").sum().reset_index()
        env_dynamics_df = env_dynamics_df[
            ["branch"]
            + [
                col
                for col in env_dynamics_df.columns
                if col == "transfers" or col == "losses"
            ]
        ]

    # Rename columns to reflect whether they are gene or env dynamics
    gene_dynamics_df = gene_dynamics_df.rename(
        {x: f"gene_{x}" for x in gene_dynamics_df.columns[1:]}, axis=1
    )
    env_dynamics_df = env_dynamics_df.rename(
        {x: f"ecosystem_{x}" for x in env_dynamics_df.columns[1:]}, axis=1
    )

    # Merge gene and env dynamics
    combined_df = pd.merge(gene_dynamics_df, env_dynamics_df, on="branch", how="outer")
    combined_df = combined_df.fillna(0)
    combined_df.columns = [col.replace("_", " ").title() for col in combined_df.columns]

    # Read in input tree and get branch lengths
    if isinstance(input_tree, str):
        input_tree = ete3.Tree(input_tree, format=1)
    branch_lengths = {
        node.name: node.dist for node in input_tree.traverse() if node.name
    }
    branch_names = [node.name for node in input_tree.traverse() if node.name]

    # Add branches that are not present in the combined dataframe to the Branch column
    missing_branches = set(branch_names) - set(combined_df["Branch"])
    missing_branches_df = pd.DataFrame(missing_branches, columns=["Branch"])
    combined_df = pd.concat([combined_df, missing_branches_df])

    # fill the NaN values in the other gene_ or ecosystem_ columns with 0, but not the Branch column
    cols_to_fill = [col for col in combined_df.columns if col != "Branch"]
    combined_df[cols_to_fill] = combined_df[cols_to_fill].fillna(0)

    # Add branch lengths to the combined dataframe
    combined_df["Branch Length"] = combined_df["Branch"].map(branch_lengths)

    combined_df = combined_df.dropna()

    return combined_df


def ols_y_pred(x_vals, y_vals):
    x_vals = sm.add_constant(x_vals)
    model = sm.OLS(y_vals, x_vals)
    res = model.fit()
    y_pred = res.predict(x_vals)
    return y_pred


# Ensure that the plots are rendered in the notebook properly
pio.renderers.default = "notebook"


def plot_gene_env_comparison(
    gene_dynamics_filepath: str,
    env_dynamics_filepath: str,
    input_tree_filepath: str,
    output_dir: str,
    plots_dir: str,
    terminal_or_internal: str = "all",
    display_stats_plots: bool = True,
    single_copy_nogs_list_filepath: str = None,
):
    """
    Plot comparison between gene and environment dynamics.

    Parameters:
    - gene_dynamics_filepath: str, path to the gene dynamics file: either varwise.branchwise or branchwise
    - env_dynamics_filepath: str, path to the environment dynamics file: either varwise.branchwise or branchwise
    - input_tree_filepath: str, path to the input tree file
    - output_dir: str, path to the output directory where the stats files will be saved
    - plots_dir : str, path to the output directory where the plots will be saved
    - terminal_or_internal: str, filter for terminal or internal branches ('all', 'internal', 'terminal')
    - display_stats_plots: bool, whether to display stats and plots
    - single_copy_nogs_list_filepath: str, path to the single copy NOGs list. If provided, the dynamics files should be varwise.branchwise

    """
    gene_method = gene_dynamics_filepath.split(".branchwise.")[1].split(".tsv")[0]
    env_method = env_dynamics_filepath.split(".branchwise.")[1].split(".tsv")[0]

    # find out if the gene and env dynamics files are varwise.branchwise or branchwise
    if "varwise.branchwise" in gene_dynamics_filepath:
        if single_copy_nogs_list_filepath is None:
            raise ValueError(
                "The gene dynamics file is varwise.branchwise. Please provide the single copy NOGs list."
            )
        if "varwise.branchwise" not in env_dynamics_filepath:
            raise ValueError(
                "The gene dynamics file is varwise.branchwise. Please provide the environment dynamics file as varwise.branchwise."
            )
    if "varwise.branchwise" in env_dynamics_filepath:
        if single_copy_nogs_list_filepath is None:
            raise ValueError(
                "The environment dynamics file is varwise.branchwise. Please provide the single copy NOGs list."
            )
        if "varwise.branchwise" not in gene_dynamics_filepath:
            raise ValueError(
                "The environment dynamics file is varwise.branchwise. Please provide the gene dynamics file as varwise.branchwise."
            )

    if single_copy_nogs_list_filepath is not None:
        # prepare combined_df for single copy genes
        combined_df = combine_gene_env_inferences(
            gene_dynamics_filepath,
            env_dynamics_filepath,
            input_tree_filepath,
            single_copy_nogs_list_filepath,
        )
    else:
        # combine gene and env inferences for each branch and get branch lengths from input tree
        combined_df = combine_gene_env_inferences(
            gene_dynamics_filepath, env_dynamics_filepath, input_tree_filepath
        )

    # get the variable columns: figure out whether 'transfers' and/or 'losses' are present in the columns
    var_columns = []
    # if any('Transfers' in col for col in combined_df.columns):
    if "Gene Transfers" in combined_df.columns:
        var_columns.append("Transfers")
    # if any('Losses' in col for col in combined_df.columns):
    if "Gene Losses" in combined_df.columns:
        var_columns.append("Losses")
    if not var_columns:
        raise ValueError(
            'The columns in the input file do not contain either "transfers" or "losses".'
        )

    # filter for terminal or internal branches if specified
    if terminal_or_internal in ["internal", "terminal"]:
        combined_df = (
            combined_df[combined_df["Branch"].str.startswith("N")]
            if terminal_or_internal == "internal"
            else combined_df[~combined_df["Branch"].str.startswith("N")]
        )

    # model gene and env dynamics with branch length to correct for the effect of branch length by taking the residuals
    for y in [f"Gene {x}" for x in var_columns] + [
        f"Ecosystem {x}" for x in var_columns
    ]:
        formula = f'Q("{y}") ~ Q("Branch Length")'
        model = sm.OLS.from_formula(formula, data=combined_df).fit()
        combined_df[f"Corrected {y}"] = combined_df[y] - model.predict(combined_df)
        # new variable denoting the ratio of the gene/env dynamics to the branch length
        combined_df[f"{y}/Branch Length"] = (
            combined_df[y] / combined_df["Branch Length"]
        )

    # save the combined dataframe to a file
    combined_df.to_csv(
        f"{output_dir}/combined_df-{terminal_or_internal}-{gene_method}-{env_method}.tsv",
        sep="\t",
        index=False,
    )

    # prepare combinations of variables (out of the specified var_columns) to compare
    # first the gene and env dynamics, corrected and uncorrected, against each other
    combinations_list = [
        (f"{prefix} {x}", f"{suffix} {x}")
        for prefix, suffix in [
            ("Ecosystem", "Gene"),
            ("Corrected Gene", "Corrected Ecosystem"),
        ]
        for x in var_columns
    ]
    # then the gene and env dynamics, corrected and uncorrected, against branch length
    combinations_list += [
        (f"Branch Length", f"{prefix} {x}")
        for prefix in ["Gene", "Ecosystem", "Corrected Gene", "Corrected Ecosystem"]
        for x in var_columns
    ]

    if display_stats_plots:
        # Create tabs for displaying data and plots
        data_tabs = widgets.Accordion(
            children=[
                widgets.Output(),  # Combined DataFrame
                widgets.Output(),  # Spearman Correlations
                widgets.Output(),  # Linear Modelling
                widgets.Output(),  # Mann-Whitney U Test
                widgets.Output(),  # Bootstrapped Mean Differences
            ]
        )
        data_tabs.set_title(0, "Dataframe")
        data_tabs.set_title(
            1,
            "Spearman Correlations among Gene Dynamics, Environment Dynamics, and Branch Length",
        )
        data_tabs.set_title(
            2,
            "Linear Modelling of Gene Dynamics as a Function of Environment Dynamics and Branch Length",
        )
        data_tabs.set_title(
            3,
            "Mann-Whitney U Test for whether Gene Dynamics differ between branch sets with Zero and Nonzero Ecosystem Transfers",
        )
        data_tabs.set_title(
            4,
            "Bootstrapped Mean Differences Test for whether Gene Dynamics differ between branch sets with Zero and Nonzero Ecosystem Transfers",
        )

        plot_tabs = widgets.Accordion(
            children=[
                widgets.Output(),  # Gene and Environment Dynamics Comparison Plot
                widgets.Output(),  # Branch Sets Plot
                widgets.Output(),  # Histogram Plot
            ]
        )
        plot_tabs.set_title(0, "Plot:Gene and Environment Dynamics Comparison")
        plot_tabs.set_title(
            1, "Plot: Gene Dynamics vs Branch Length across Different Branch Sets"
        )
        plot_tabs.set_title(
            2, "Plot: Histogram of Gene Dynamics across Different Branch Sets"
        )

        # Set the size of the plots
        plot_size_side = 400

        num_cols = 4
        num_rows = int(np.ceil(len(combinations_list) / num_cols))
        fig_gd_ed = make_subplots(
            rows=num_rows, cols=num_cols, vertical_spacing=0.05, horizontal_spacing=0.1
        )
        fig_gd_ed.update_layout(
            width=plot_size_side * num_cols, height=plot_size_side * num_rows
        )

    # find spearman correlations between the variables and plot the comparisons
    spearman_correlations = {}

    for i, (x, y) in enumerate(combinations_list):
        if combined_df[x].nunique() == 1 or combined_df[y].nunique() == 1:
            print(
                f"Skipping Spearman correlation for {
                  x} and {y} since one the variables has a constant value."
            )
            continue
        spearman_correlations[(x, y)] = st.spearmanr(combined_df[x], combined_df[y])

        if display:
            fig_gd_ed.add_trace(
                go.Scatter(
                    x=combined_df[x], y=combined_df[y], mode="markers", showlegend=False
                ),
                row=i // num_cols + 1,
                col=i % num_cols + 1,
            )
            model = sm.OLS(combined_df[y], sm.add_constant(combined_df[x])).fit()
            fig_gd_ed.add_trace(
                go.Scatter(
                    x=combined_df[x],
                    y=model.predict(sm.add_constant(combined_df[x])),
                    mode="lines",
                    line=dict(color="black"),
                    showlegend=False,
                ),
                row=i // num_cols + 1,
                col=i % num_cols + 1,
            )
            fig_gd_ed.update_xaxes(
                title_text=x, row=i // num_cols + 1, col=i % num_cols + 1
            )
            fig_gd_ed.update_yaxes(
                title_text=y, row=i // num_cols + 1, col=i % num_cols + 1
            )

    # Save spearman correlations to a file and display
    spearman_correlations_df = pd.DataFrame.from_dict(
        spearman_correlations, orient="index", columns=["Rho", "p-value"]
    )
    spearman_correlations_df["significant"] = spearman_correlations_df["p-value"].apply(
        lambda x: "*" if x < 0.05 else ""
    )
    spearman_correlations_df.index.name = "Variable-Pair"
    spearman_correlations_df["Terminal or Internal"] = terminal_or_internal
    spearman_correlations_df.to_csv(
        f"{output_dir}/spearman_correlations-{terminal_or_internal}-{gene_method}-{env_method}.tsv",
        sep="\t",
    )

    if display_stats_plots:
        fig_gd_ed.update_layout(title_text="Gene and Environment Dynamics Comparison")

        with plot_tabs.children[0]:
            fig_gd_ed.show()
            # save the plot to svg
            fig_gd_ed.write_image(
                f"{plots_dir}/gene_env_dynamics_comparison-{terminal_or_internal}-{gene_method}-{env_method}.png"
            )

        # Display combined DataFrame
        with data_tabs.children[0]:
            display(combined_df)
        with data_tabs.children[1]:
            display(spearman_correlations_df)

        ################################################# MODELLING GD ~ ED + BL ########################################################

        # Linear modelling of gene dynamics as a function of env dynamics and branch length
        with data_tabs.children[2]:
            for x in var_columns:
                display(
                    linear_modelling(
                        combined_df,
                        f'Q("Gene {
                        x}") ~ Q("Ecosystem {x}") + Q("Branch Length")',
                    )
                )

        # Plot gene dynamics as a function of branch length for different branch sets
        # Each branch set is defined by the value of the 'Ecosystem Transfers' or 'Ecosystem Losses' column
        def plot_branch_sets(branch_sets, var_column, col, fig_obj):
            for branch_set_df, name, color in branch_sets:
                if branch_set_df.empty:
                    print(
                        f"Skipping {name} branch set for {
                        var_column} due to empty set."
                    )
                    continue

                fig_obj.add_trace(
                    go.Scatter(
                        x=branch_set_df["Branch Length"],
                        y=branch_set_df[
                            f"Gene {
                    var_column}"
                        ],
                        mode="markers",
                        name=name,
                        opacity=0.8,
                        showlegend=True,
                    ),
                    row=i + 1,
                    col=col,
                )
                fig_obj.add_trace(
                    go.Scatter(
                        x=branch_set_df["Branch Length"],
                        y=ols_y_pred(
                            branch_set_df["Branch Length"],
                            branch_set_df[f"Gene {var_column}"],
                        ),
                        name=f"OLS: {name}",
                        mode="lines",
                        line=dict(color=color),
                        showlegend=True,
                    ),
                    row=i + 1,
                    col=col,
                )

        if "Transfers" in var_columns:
            branch_sets_transfers = [
                (
                    combined_df[combined_df["Ecosystem Transfers"] == 0],
                    "Zero Ecosystem Transfers",
                    "black",
                ),
                (
                    combined_df[combined_df["Ecosystem Transfers"] != 0],
                    "Nonzero Ecosystem Transfers",
                    "red",
                ),
            ]
        if "Losses" in var_columns:
            branch_sets_losses = [
                (
                    combined_df[combined_df["Ecosystem Losses"] == 0],
                    "Zero Ecosystem Losses",
                    "black",
                ),
                (
                    combined_df[combined_df["Ecosystem Losses"] != 0],
                    "Nonzero Ecosystem Losses",
                    "red",
                ),
            ]

        # print('Number of branches with:\nZero Ecosystem Transfers:',
        #       branch_sets_transfers[0][0].shape[0], '\n',
        #       'Nonzero Ecosystem Transfers:', branch_sets_transfers[1][0].shape[0], '\n',
        #       'Zero Ecosystem Losses:', branch_sets_losses[0][0].shape[0], '\n',
        #       'Nonzero Ecosystem Losses:', branch_sets_losses[1][0].shape[0])

        fig_branch_sets = make_subplots(
            rows=len(var_columns), cols=2, vertical_spacing=0.05, horizontal_spacing=0.1
        )
        fig_branch_sets.update_layout(
            width=plot_size_side * len(var_columns) * 4,
            height=plot_size_side
        )
        for i, var_column in enumerate(var_columns):
            if "Transfers" in var_columns:
                plot_branch_sets(branch_sets_transfers, var_column, 1, fig_branch_sets)
            if "Losses" in var_columns:
                plot_branch_sets(branch_sets_losses, var_column, 2, fig_branch_sets)
            fig_branch_sets.update_xaxes(title_text="Branch Length", row=i + 1, col=1)
            fig_branch_sets.update_yaxes(
                title_text=f"Gene {var_column}", row=i + 1, col=1
            )
            fig_branch_sets.update_xaxes(title_text="Branch Length", row=i + 1, col=2)
            fig_branch_sets.update_yaxes(
                title_text=f"Gene {var_column}", row=i + 1, col=2
            )

        with plot_tabs.children[1]:
            fig_branch_sets.show()
            # save the plot to svg
            fig_branch_sets.write_image(
                f"{plots_dir}/gene_dynamics_vs_branch_length-{terminal_or_internal}-{gene_method}-{env_method}.png"
            )

        fig_histogram = make_subplots(
            rows=len(var_columns), cols=2, vertical_spacing=0.05, horizontal_spacing=0.1
        )
        fig_histogram.update_layout(
            width=plot_size_side * len(var_columns) * 4, 
            height=plot_size_side
        )
        for i, y in enumerate(
            [f"Gene {x}" for x in var_columns]
            + [f"Gene {x}/Branch Length" for x in var_columns]
        ):
            for branch_set, name, color in branch_sets_transfers:
                fig_histogram.add_trace(
                    go.Histogram(
                        x=branch_set[y], name=name, opacity=0.8, showlegend=True
                    ),
                    row=i // 2 + 1,
                    col=i % 2 + 1,
                )
            if np.abs(branch_sets_transfers[0][0][y].skew()) > 1:
                fig_histogram.update_yaxes(
                    type="log",
                    row=i // 2 + 1,
                    col=i % 2 + 1,
                    title_text="Frequency (log scale)",
                )
            else:
                fig_histogram.update_yaxes(
                    title_text="Frequency", row=i // 2 + 1, col=i % 2 + 1
                )
            fig_histogram.update_xaxes(title_text=y, row=i // 2 + 1, col=i % 2 + 1)
            fig_histogram.update_layout(barmode="overlay")

        with plot_tabs.children[2]:
            fig_histogram.show()
            # save the plot to svg
            fig_histogram.write_image(
                f"{plots_dir}/histogram_gene_dynamics-{terminal_or_internal}-{gene_method}-{env_method}.png"
            )

    ################################################# Compare Branch Sets ########################################################

    def mann_whitney_u(branch_sets, y):
        # remove NaN and infinite values
        branch_set1 = branch_sets[0][0][y].replace([np.inf, -np.inf], np.nan).dropna()
        branch_set2 = branch_sets[1][0][y].replace([np.inf, -np.inf], np.nan).dropna()
        return [
            {
                "Variable": y,
                "Mann-Whitney U": st.mannwhitneyu(
                    branch_set1, branch_set2, alternative="two-sided"
                )[0],
                "p-value": st.mannwhitneyu(
                    branch_set1, branch_set2, alternative="two-sided"
                )[1],
                "CLES": st.mannwhitneyu(
                    branch_set1, branch_set2, alternative="two-sided"
                )[0] / (branch_set1.shape[0] * branch_set2.shape[0]),
                f"Mean HGT with {branch_sets[0][1]}": branch_set1.mean(),
                f"Mean HGT with {branch_sets[1][1]}": branch_set2.mean(),
                f"Number of branches with {branch_sets[0][1]}": branch_set1.shape[0],
                f"Number of branches with {branch_sets[1][1]}": branch_set2.shape[0],
            }
        ]

    mwu_data = []
    for var in var_columns:
        mwu_data.append(mann_whitney_u(branch_sets_transfers, f"Gene {var}"))
        mwu_data.append(
            mann_whitney_u(branch_sets_transfers, f"Gene {var}/Branch Length")
        )
    mwu_df = pd.DataFrame([item for sublist in mwu_data for item in sublist])
    mwu_df["Terminal or Internal"] = terminal_or_internal
    mwu_df["significant"] = mwu_df["p-value"].apply(lambda x: "*" if x < 0.05 else "")
    # write this dataframe to a file
    mwu_df.to_csv(
        f"{output_dir}/mwu-{terminal_or_internal}-{gene_method}-{env_method}.tsv",
        sep="\t",
        index=False,
        header=True,
    )

    combined_set = pd.concat([branch_sets_transfers[0][0], branch_sets_transfers[1][0]])
    randomized_mean_diffs = []
    randomize_iter = 1e4
    for y in [f"Gene {x}" for x in var_columns] + [
        f"Gene {x}/Branch Length" for x in var_columns
    ]:
        # if y has a constant value, skip
        if combined_set[y].nunique() == 1:
            print(
                f"Skipping randomized mean difference for {y} due to constant value."
            )
            continue
        observed_mean_diff = (
            branch_sets_transfers[0][0][y].mean()
            - branch_sets_transfers[1][0][y].mean()
        )
        # sample from combined_set two sets of the same size as zero and nonzero transfers branch sets, then take the difference in means
        # do this for a large number of iterations to get the distribution of mean differences at random
        randomized_means = [
            combined_set.sample(branch_sets_transfers[0][0].shape[0], replace=True)[
                y
            ].mean()
            - combined_set.sample(branch_sets_transfers[1][0].shape[0], replace=True)[
                y
            ].mean()
            for _ in range(int(randomize_iter))
        ]
        randomized_mean_diffs.append(
            {
                "Variable": y,
                "Num branches w Zero Ecosystem Transfers": branch_sets_transfers[0][
                    0
                ].shape[0],
                "Num branches w Non Zero Ecosystem Transfers": branch_sets_transfers[1][
                    0
                ].shape[0],
                "Mean HGT with Zero Ecosystem Dynamics": branch_sets_transfers[0][0][
                    y
                ].mean(),
                "Mean HGT with Non Zero Ecosystem Dynamics": branch_sets_transfers[1][
                    0
                ][y].mean(),
                "Observed Mean Difference": observed_mean_diff,
                "p-value": (
                    np.abs(randomized_means) > np.abs(observed_mean_diff)
                ).mean(),
            }
        )
    randomized_mean_diffs_df = pd.DataFrame(randomized_mean_diffs)
    randomized_mean_diffs_df["significant"] = randomized_mean_diffs_df[
        "p-value"
    ].apply(lambda x: "*" if x < 0.05 else ".")
    randomized_mean_diffs_df["Terminal or Internal"] = terminal_or_internal

    # remove rows where variable doesn't have 'Branch Length' in the name
    randomized_mean_diffs_df2 = randomized_mean_diffs_df[
        randomized_mean_diffs_df["Variable"].str.contains("Branch Length")
    ]

    randomized_mean_diffs_df2.to_csv(
        f"{output_dir}/randomized_mean_diffs-{terminal_or_internal}-{
        gene_method}-{env_method}.tsv",
        sep="\t",
        index=False,
        header=True,
    )

    if display_stats_plots:
        with data_tabs.children[3]:
            print(
                f"Mann-Whitney U test for Gene Dynamics between branch sets with Zero and Nonzero Ecosystem Transfers:"
            )
            display(mwu_df)
        with data_tabs.children[4]:
            print(
                f"Bootstrapped mean differences comparing Gene Dynamics between branch sets with Zero and Nonzero Ecosystem Transfers:"
            )
            display(randomized_mean_diffs_df)

            # Display the tabs with data and plots
    display(data_tabs)
    display(plot_tabs)

    return None


def plot_box_plots_comparing_gene_env(combined_df_filepaths_list, plots_dir):
    """
    Plot box plots comparing gene and environment dynamics.
    Requires a list of filepaths to the `combined_df` files generated by the `plot_gene_env_comparison` function.
    For each file, this function extracts the gene transfers per unit branch length, and environment transfers, for each branch.
    It then plots box plots using pyplot with split boxs, such that each box is for a combination of gene transfer inference method and environment transfer inference method.
    Each box is split into two parts: one for branches with zero environment transfers and one for branches with non-zero environment transfers.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Read in the combined_df files
    all_combined_dfs_list = []
    for combined_df_filepath in combined_df_filepaths_list:
        gene_method = ' '.join(combined_df_filepath.split("-")[2].split(".")[1:]).upper()

        combined_df = pd.read_csv(combined_df_filepath, sep="\t")
        # keep only the columns of interest
        combined_df = combined_df[["Gene Transfers/Branch Length", "Ecosystem Transfers"]]
        combined_df['Ecosystem Transfers'] = combined_df["Ecosystem Transfers"].apply(lambda x: "Ecosystem transfers" if x > 0 else "No ecosystem transfers")
        combined_df['Gene Method'] = gene_method

        all_combined_dfs_list.append(combined_df)

    all_combined_dfs = pd.concat(all_combined_dfs_list).reset_index(drop=True)
    all_combined_dfs = all_combined_dfs.dropna()

    # Check for NaN and infinite values
    if all_combined_dfs.isnull().values.any():
        print("NaN values present in the data")
        display(all_combined_dfs[all_combined_dfs.isnull().any(axis=1)])
    if all_combined_dfs[all_combined_dfs.isin([np.nan, np.inf, -np.inf]).any(axis=1)].shape[0] > 0:
        print("Infinite values present in the data")
        display(all_combined_dfs[all_combined_dfs.isin([np.nan, np.inf, -np.inf]).any(axis=1)])

    # Plot the box plots. One box for each gene method, split by branches with zero and non-zero ecosystem transfers
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(
        data=all_combined_dfs,
        x="Gene Method",
        y="Gene Transfers/Branch Length",
        hue="Ecosystem Transfers",
        # data=all_combined_dfs,
        # split=True,
        # # gap=0.1,
        # inner="quartile",
        log_scale=True,
    )
    # plt.title("Violin plots comparing gene and environment dynamics")
    plt.ylabel("Rate of HGT")
    plt.xlabel("HGT Inference Method")
    plt.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()

    # save as png
    plt.savefig(f"{plots_dir}/box_plots_compare_gd_ed.png")

    plt.show()
