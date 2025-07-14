from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import scipy as sp
import scipy.stats as st
from IPython.display import display


def linear_regression_logy_x(x, y):
    log_y = np.log(y)
    x = sm.add_constant(x)
    model = sm.OLS(log_y, x).fit()
    a, b = np.exp(model.params[0]), model.params[1]
    y_hat = a * np.exp(b * x[:, 1])
    return y_hat, a, b


# 2 rows, 2 columns
# 1,1: histplot of #genes in branches with/without ecosystem gains
# 2,1: scatterplot of #genes vsf {rate}
# 2,2: vertical histplot off {rate} in branches with/without ecosystem gains
# columns share x axis
# rows share y axis

def analyze_genome_size_vs_hgt(
    meth,
    compiled_gene_dynamics_df,
    branch_lengths_dict,
    ecotype_count_df,
    plots_dir,
    genome_size_col="genes",
    env="pathogenicity",
):
    """
    Analyze the relationship between genome size and HGT rates.
    """
    # dicts to be returned
    # Spearman correlations
    rates_genomesize_spearmans_list = []
    rates_genomesize_linearcoeffs_list = []
    corrected_rates_genomesize_spearmans_records_list = []
    zero_rate_branchcounts_records_list = []
    # Mann-Whitney U test results
    rates_mwu_records_list = []
    corrected_rates_mwu_records_list = []
    genomesize_mwu_records_list = []

    # opacity of scatterplot and histogram
    scatter_alpha = 0.8
    hist_alpha = 0.4

    # take the first letter of `env`, capitalise it, and add a G to it
    env_label1 = env[0].upper() + "G"
    env_label0 = "N" + env_label1

    compiled_gene_dynamics_df_copy = compiled_gene_dynamics_df.copy()

    compiled_gene_dynamics_df_copy["branch_length"] = (
        compiled_gene_dynamics_df_copy["branch"].map(branch_lengths_dict.get).fillna(0)
    )
    compiled_gene_dynamics_df_copy["branch_type"] = compiled_gene_dynamics_df_copy["branch"].map(
        lambda x: (
            env_label1
            if x
            in ecotype_count_df.loc[
                ecotype_count_df["branch_type"] == env_label1, "branch"
            ].values
            else env_label0
        )
    )
    compiled_gene_dynamics_df_copy["gain_rate"] = (
        compiled_gene_dynamics_df_copy["transfers"]
        / compiled_gene_dynamics_df_copy["branch_length"]
    )
    compiled_gene_dynamics_df_copy["loss_rate"] = (
        compiled_gene_dynamics_df_copy["losses"] / compiled_gene_dynamics_df_copy["branch_length"]
    )
    compiled_gene_dynamics_gt0_df = compiled_gene_dynamics_df_copy[
        compiled_gene_dynamics_df_copy[genome_size_col] > 0
    ].copy()

    for rate in ["gain_rate", "loss_rate"]:
        fig = plt.figure(figsize=(12, 8))
        gs = GridSpec(2, 2, figure=fig, width_ratios=[4, 1], height_ratios=[1, 4])

        ax_marginal_x = fig.add_subplot(gs[0])
        ax_legend = fig.add_subplot(gs[1])
        ax_joint = fig.add_subplot(gs[2])
        ax_marginal_y = fig.add_subplot(gs[3])

        # increase gap between the joint and the top marginal plot
        gs.update(hspace=0.1, wspace=0.1)

        zero_rate_branchcounts_records_list.append(
            {
                "rate": rate,
                "method": meth,
                "zero_rate_branch_count": compiled_gene_dynamics_gt0_df[
                    compiled_gene_dynamics_gt0_df[rate] == 0
                ].shape[0],
            }
        )

        compiled_gene_dynamics_gt0_df = compiled_gene_dynamics_gt0_df[
            compiled_gene_dynamics_gt0_df[rate] > 0
        ].copy()

        # plot the scatterplot of #genes vs {rate}
        sns.scatterplot(
            data=compiled_gene_dynamics_gt0_df,
            x="genes",
            y=rate,
            hue="branch_type",
            ax=ax_joint,
            alpha=scatter_alpha,
        )
        # sns.despine(ax=ax_joint, offset=2)

        x, y = compiled_gene_dynamics_gt0_df["genes"].values, compiled_gene_dynamics_gt0_df[rate]
        y_hat, a, b = linear_regression_logy_x(x, y)
        sorted_indices = np.argsort(x)
        x_sorted, y_hat_sorted = x[sorted_indices], y_hat[sorted_indices]

        ax_joint.plot(x_sorted, y_hat_sorted, "r-")
        # Spearman correlation
        rates_genomesize_spearmans_list.append(
            {
                "rate": rate,
                "method": meth,
                "branch_set": "all_branches",
                "spearman_corr": sp.stats.spearmanr(
                    compiled_gene_dynamics_df_copy["genes"], compiled_gene_dynamics_df_copy[rate]
                )[0],
                "p_value": sp.stats.spearmanr(
                    compiled_gene_dynamics_df_copy["genes"], compiled_gene_dynamics_df_copy[rate]
                )[1],
            }
        )
        rates_genomesize_spearmans_list.append(
            {
                "rate": rate,
                "method": meth,
                "branch_set": "branches_with_gains",
                "spearman_corr": sp.stats.spearmanr(
                    compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, "genes"
                    ],
                    compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, rate
                    ],
                )[0],
                "p_value": sp.stats.spearmanr(
                    compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, "genes"
                    ],
                    compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, rate
                    ],
                )[1],
            }
        )
        rates_genomesize_linearcoeffs_list.append(
            {
                "rate": rate,
                "method": meth,
                "a": a,
                "b": b,
                "num_branches": compiled_gene_dynamics_gt0_df.shape[0],
            }
        )
        # Print the linear regression coefficients and R2
        print(f"{meth}, {rate}: log(y) = {a:.2f} + {b:.2e} * x",
                f"\n (Pearson r = {sp.stats.pearsonr(compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, "genes"
                    ],
                    compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, rate
                    ])[0]:.2f}, p-value = {sp.stats.pearsonr(compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, "genes"
                    ],
                    compiled_gene_dynamics_df_copy.loc[
                        compiled_gene_dynamics_df_copy[rate] > 0, rate
                    ])[1]:.2e})")

        ax_joint.get_legend().set_visible(False)
        ax_joint.set_xlabel("Number of genes")
        ax_joint.set_ylabel("HGT rate" if rate == "gain_rate" else "Gene loss rate")
        ax_joint.set_yscale("log")
        ax_joint.set_ylim(-1, ax_joint.get_ylim()[1] * 1.2)
        # ax_joint.grid(True, which="both", linestyle="--", linewidth=0.5)

        rate_residuals = compiled_gene_dynamics_gt0_df[rate] - y_hat
        compiled_gene_dynamics_gt0_df.loc[
            compiled_gene_dynamics_gt0_df[rate] > 0, f"{rate}_residuals"
        ] = rate_residuals

        ax_legend.axis("off")
        handles, labels = ax_joint.get_legend_handles_labels()
        ax_legend.legend(
            handles,
            labels,
            loc="upper right",
            fontsize=ax_joint.xaxis.label.get_size() * 0.8,
            facecolor="white",
            edgecolor="black",
            frameon=False,
            markerscale=2,
        )

        ############## MARGINAL PLOTS ##############

        # plot the marginal histogram of the x axis
        sns.histplot(
            data=compiled_gene_dynamics_gt0_df,
            x="genes",
            hue="branch_type",
            stat="percent",
            common_norm=False, # normalize each histogram separately
            ax=ax_marginal_x,
            bins=30,
            kde=True,
            alpha=hist_alpha,
        )
        sns.despine(ax=ax_marginal_x, offset={"bottom": 4})
        for line in ax_marginal_x.get_lines():
            line.set_markersize(0)
        ax_marginal_x.set_xlabel("")
        ax_marginal_x.set_ylabel("")
        medians = []
        for branch_type in ecotype_count_df["branch_type"].unique():
            median = compiled_gene_dynamics_df_copy[
                compiled_gene_dynamics_df_copy["branch_type"] == branch_type
            ]["genes"].median()
            ax_marginal_x.axvline(
                median, color="black", linestyle="--", label=f"{branch_type}"
            )
            ax_marginal_x.text(
                median,
                ax_marginal_x.get_ylim()[1] * 1.2,
                branch_type,
                color="black",
                ha="right",
                va="top",
                rotation=90,
                fontsize=13,
            )
            medians.append(median)
        ax_marginal_x.get_legend().set_visible(False)
        ax_marginal_x.set_xticklabels([])

        # Print Mann-Whitney U test results and medians
        mwu, pval = st.mannwhitneyu(
            compiled_gene_dynamics_df_copy[
                (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
            ][rate],
            compiled_gene_dynamics_df_copy[
                (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
            ][rate],
            alternative="two-sided",
        )
        rates_mwu_records_list.append(
            {
                "rate": rate,
                "method": meth,
                "MWU statistic": mwu,
                "p-value": pval,
                f"{env_label1} Median": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
                    & (compiled_gene_dynamics_df_copy[rate] > 0)
                ][rate].median(),
                f"{env_label0} Median": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
                    & (compiled_gene_dynamics_df_copy[rate] > 0)
                ][rate].median(),
                f"{env_label1} Mean": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
                ][rate].mean(),
                f"{env_label0} Mean": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
                ][rate].mean(),
                "Effect size (CLES)": mwu
                / (
                    compiled_gene_dynamics_df_copy[
                        (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
                        # & (compiled_gene_dynamics_df_copy[rate] > 0)
                    ][rate].shape[0]
                    * compiled_gene_dynamics_df_copy[
                        (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
                        # & (compiled_gene_dynamics_df_copy[rate] > 0)
                    ][rate].shape[0]
                ),
                f"{env_label1} Size": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
                ][rate].shape[0],
                f"{env_label0} Size": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
                ][rate].shape[0],
            }
        )
        # display rows with EG branches
        print(
            f"\n{meth}, {rate}:\n")
        display(
            compiled_gene_dynamics_df_copy[
                (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
            ]
            .sort_values(by=rate, ascending=True)
        )

        # find mwu for genome size
        mwu_genomesize, pval_mwu_genomesize = st.mannwhitneyu(
            compiled_gene_dynamics_df_copy[
                (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
            ]["genes"],
            compiled_gene_dynamics_df_copy[
                (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
            ]["genes"],
            alternative="two-sided",
        )
        genomesize_mwu_records_list.append(
            {
                "rate": rate,
                "method": meth,
                "MWU statistic": mwu_genomesize,
                "p-value": pval_mwu_genomesize,
                f"{env_label1} Median": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
                ]["genes"].median(),
                f"{env_label0} Median": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
                ]["genes"].median(),
                f"{env_label1} Mean": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
                ]["genes"].mean(),
                f"{env_label0} Mean": compiled_gene_dynamics_df_copy[
                    (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
                ]["genes"].mean(),
                "Effect size (CLES)": mwu_genomesize
                / (
                    compiled_gene_dynamics_df_copy[
                        (compiled_gene_dynamics_df_copy["branch_type"] == env_label1)
                    ]["genes"].shape[0]
                    * compiled_gene_dynamics_df_copy[
                        (compiled_gene_dynamics_df_copy["branch_type"] == env_label0)
                    ]["genes"].shape[0]
                ),
            }
        )

        # plot the marginal histogram of the y axis
        sns.histplot(
            data=compiled_gene_dynamics_gt0_df,
            y=rate,
            hue="branch_type",
            stat="percent",
            common_norm=False, # normalize each histogram separately
            ax=ax_marginal_y,
            bins=30,
            kde=True,
            log_scale=True,
            alpha=hist_alpha,
        )
        sns.despine(ax=ax_marginal_y, offset={"left": 4})
        for line in ax_marginal_y.get_lines():
            line.set_markersize(0)
        ax_marginal_y.set_ylabel("")
        ax_marginal_y.set_xlabel("")
        for branch_type in ecotype_count_df["branch_type"].unique():
            median = compiled_gene_dynamics_gt0_df[
            compiled_gene_dynamics_gt0_df["branch_type"] == branch_type
            ][rate].median()
            ax_marginal_y.axhline(
            median, color="black", linestyle="--", label=f"{branch_type}"
            )
            ax_marginal_y.text(
            ax_marginal_y.get_xlim()[1],
            median,
            branch_type,
            color="black",
            ha="right",
            va="bottom",
            rotation=0,
            fontsize=13,
            )
        ax_marginal_y.get_legend().set_visible(False)
        ax_marginal_y.set_ylim(ax_joint.get_ylim())
        ax_marginal_y.set_yscale("log")
        ax_marginal_y.set_yticklabels([])  # Ensure y-axis tick labels are not shown
        plt.subplots_adjust(hspace=0.05, wspace=0.05)
        # save the figure for this method
        plt.savefig(
            f"{plots_dir}/genome_size_analysis_combinedplot_{meth}_{rate}_{env}.png",
            dpi=300,
            bbox_inches="tight"
        )
        plt.show()

        return {
            "Spearman correlation between rates and genome size": rates_genomesize_spearmans_list,
            "Linear regression coefficients for rates vs genome size": rates_genomesize_linearcoeffs_list,
            "Spearman correlation for corrected rates vs genome size": corrected_rates_genomesize_spearmans_records_list,
            "Count of branches with zero rates": zero_rate_branchcounts_records_list,
            "Mann-Whitney U test results for rates between EG and NEG branches": rates_mwu_records_list,
            "Mann-Whitney U test results for corrected rates between EG and NEG branches": corrected_rates_mwu_records_list,
            "Mann-Whitney U test results for genome size between EG and NEG branches": genomesize_mwu_records_list,
        }


def combine_all_stats_records(rates_n_stats_all_dict):
    """
    Combine all records for each method into a list of record lists and create dataframes.

    Parameters:
    rates_n_stats_all_dict (dict): Dictionary containing rates and stats for each method.

    Returns:
    dict: A dictionary of dataframes for each combined stat.
    """
    # Combine all records for each method into a list of record lists
    all_stats_dict = {}
    for meth, rates_n_stats_meth_dict in rates_n_stats_all_dict.items():
        for key, records_list in rates_n_stats_meth_dict.items():
            if key not in all_stats_dict:
                all_stats_dict[key] = []
            all_stats_dict[key].extend(records_list)

    # Make dataframes for each of these stats
    all_stats_dfs_dict = {}
    for key, records_list in all_stats_dict.items():
        all_stats_dfs_dict[key] = pd.DataFrame.from_records(records_list)

    # Display the dataframes
    for key, df in all_stats_dfs_dict.items():
        print(f"\n{key}:\n")
        display(df)

    return all_stats_dfs_dict


def extract_branchwise_subset_dfs(
    compiled_varwise_branchwise_gene_dynamics_dfs_dict: dict,
    subset_var_list: list,
    compiled_res_dir: str,
    subset_name: str,
    var_name: str = "nog_id",
) -> dict:
    """
    Given a dict of method->varwise_branchwise_gene_dynamics_df values,
    prepare a dict of method->subset_branchwise_gene_dynamics_df values
    for the given subset_var_list.
    Each of these new dfs can then be fed into the analyze_genome_size_vs_hgt function.
    """
    subset_branchwise_gene_dynamics_dfs_dict = {}
    for method, df in compiled_varwise_branchwise_gene_dynamics_dfs_dict.items():
        try:
            subset_branchwise_gene_dynamics_dfs_dict[method] = df[
                df[var_name].isin(subset_var_list)
            ].copy()
            # groupby branch and sum the values
            subset_branchwise_gene_dynamics_dfs_dict[method] = subset_branchwise_gene_dynamics_dfs_dict[
                method
            ].groupby("branch").sum().reset_index()
            # keep only any of these columns if they exist
            cols_to_keep = ['branch', 'transfers', 'losses', 'expansions', 'reductions']
            cols_to_keep = [col for col in cols_to_keep if col in subset_branchwise_gene_dynamics_dfs_dict[method].columns]
            subset_branchwise_gene_dynamics_dfs_dict[method] = subset_branchwise_gene_dynamics_dfs_dict[
                method
            ][cols_to_keep]

            # write this subset branchwise df to a file in compiled_res_dir
            subset_branchwise_gene_dynamics_filepath = f"{compiled_res_dir}/compiled_transfers.branchwise.gene.{subset_name}.{method}.tsv"
            subset_branchwise_gene_dynamics_dfs_dict[method].to_csv(
                subset_branchwise_gene_dynamics_filepath,
                sep="\t",
                index=False,
            )
        except Exception as e:
            print(
                f"Error processing method {method} for subset {subset_name}: {e}"
            )
            raise

    return subset_branchwise_gene_dynamics_dfs_dict