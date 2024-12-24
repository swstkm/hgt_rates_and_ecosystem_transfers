import os
import argparse
import subprocess
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
import sys
import traceback

from loguru import logger


def run_hyphy_busted_single(args):
    """
    Run HyPhy BUSTED on a single alignment and treefile
    """
    alignment, treefile, output_dir, path_to_busted_script, max_retries, two_rate_classes = args

    # Create output file
    output_file = os.path.join(output_dir, os.path.basename(
        alignment).replace(".aln.fna", ".busted.json"))

    # Define command
    command = [
        "hyphy", path_to_busted_script,
        "--tree", treefile,
        "--alignment", alignment,
        "--output", output_file,
        "--kill-zero-lengths", "No",  # Do not remove zero-length branches
        "--srv", "Yes",  # Synonymous rate variation across sites is allowed
    ]

    if two_rate_classes:
        command += [
            # Branches marked with suffix {Test} will be tested as foreground branches, for selection
            "--branches", "Test",
            "--rates", "2",  # Two rate classes
            "--syn-rates", "2"  # Two synonymous rate classes
        ]

    for i in range(max_retries):
        try:
            # Run command and capture stdout and stderr
            logger.info(f"Running command: {' '.join(command)}")

            process = subprocess.Popen(
                command,
                # shell=True,  # if shell is True, hyphy seems to get into interactive mode instead of running the command
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # Log stdout and stderr to a file specific to this run
            log_file = output_file.replace(".json", ".log")
            with open(log_file, "w") as f:
                f.write(f"{'--' * 30}\nDumping stdout:\n{'--' * 30}\n{stdout.decode()
                                                                      }\n\n\n{
                    '--' * 30}\nDumping stderr\n{
                    '--' * 30}\n{stderr.decode()}")

            # If process returns non-zero exit code, log error and retry if max_retries not reached
            if process.returncode != 0:
                logger.error(
                    f"Error running HyPhy BUSTED for {alignment}: {stderr.decode()}")
                # retry
                if i < max_retries - 1:
                    logger.info(
                        f"Retrying {i+1}/{max_retries} for {alignment}")
                    continue
                elif i == max_retries - 1:
                    logger.error(
                        f"Failed to run HyPhy BUSTED for {alignment} after {max_retries} retries")
                    return None
            else:
                logger.info(f"HyPhy BUSTED run for {
                            alignment} completed successfully and output written to {output_file}")
                return output_file
        except Exception as e:
            logger.error(f"Error running HyPhy BUSTED for {alignment}: {e}")
            logger.error(traceback.format_exc())
            return None

    return output_file


def run_hyphy_busted_parallel(args):
    """
    Run HyPhy BUSTED in parallel on a list of alignments and treefiles
    """

    input_file = args.input_file
    output_dir = args.output_dir
    path_to_busted_script = args.path_to_busted_script
    debug = args.debug
    max_retries = args.retries
    two_rate_classes = args.two_rate_classes
    max_threads = args.threads

    # Read input file
    input_df = pd.read_csv(input_file, sep="\t", header=None,
                           names=["alignment", "treefile"])
    if debug:
        input_df = input_df.head(5)
        max_threads = min(5, max_threads)
        logger.info(
            "Debug mode enabled. Running on first 5 rows of the input file")

    # Create output directory
    if output_dir is None:
        output_dir = os.path.dirname(input_file)
    os.makedirs(output_dir, exist_ok=True)
    # get absolute path
    output_dir = os.path.abspath(output_dir)

    # Run HyPhy BUSTED
    args = [(row.alignment, row.treefile, output_dir, path_to_busted_script, max_retries, two_rate_classes)
            for row in input_df.itertuples()]
    with mp.Pool(max_threads) as pool:
        outfiles = list(tqdm(pool.imap(run_hyphy_busted_single, args), total=len(
            args), desc="Running HyPhy BUSTED in parallel"))

    return outfiles


def run_hyphy_relax_single(args):
    """
    Run HyPhy RELAX on a single alignment and treefile
    """
    alignment, treefile, output_dir, path_to_relax_script, max_retries = args

    # Create output file
    output_file = os.path.join(output_dir, os.path.basename(
        alignment).replace(".aln.fna", ".relax.json"))
    
    # Define command
    command = [
        "hyphy", path_to_relax_script,
        "--tree", treefile,
        "--alignment", alignment,
        "--output", output_file,
        "--test", "Test", # Branches marked with suffix {Test} will be tested as foreground branches, for relaxation/intensification of selection
        "--kill-zero-lengths", "No",  # Do not remove zero-length branches
    ]

    for i in range(max_retries):
        try:
            # Run command and capture stdout and stderr
            logger.info(f"Running command: {' '.join(command)}")

            process = subprocess.Popen(
                command,
                # shell=True,  # if shell is True, hyphy seems to get into interactive mode instead of running the command
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # Log stdout and stderr to a file specific to this run
            log_file = output_file.replace(".json", ".log")
            with open(log_file, "w") as f:
                f.write(f"{'--' * 30}\nDumping stdout:\n{'--' * 30}\n{stdout.decode()
                                                                      }\n\n\n{
                    '--' * 30}\nDumping stderr\n{
                    '--' * 30}\n{stderr.decode()}")

            # If process returns non-zero exit code, log error and retry if max_retries not reached
            if process.returncode != 0:
                logger.error(
                    f"Error running HyPhy RELAX for {alignment}: {stderr.decode()}")
                # retry
                if i < max_retries - 1:
                    logger.info(
                        f"Retrying {i+1}/{max_retries} for {alignment}")
                    continue
                elif i == max_retries - 1:
                    logger.error(
                        f"Failed to run HyPhy RELAX for {alignment} after {max_retries} retries")
                    return None
            else:
                logger.info(f"HyPhy RELAX run for {
                            alignment} completed successfully and output written to {output_file}")
                return output_file
        except Exception as e:
            logger.error(f"Error running HyPhy RELAX for {alignment}: {e}")
            logger.error(traceback.format_exc())
            return None
    
    return output_file


def run_hyphy_relax_parallel(args):
    """
    Run HyPhy RELAX in parallel on a list of alignments and treefiles
    """

    input_file = args.input_file
    output_dir = args.output_dir
    path_to_relax_script = args.path_to_relax_script
    debug = args.debug
    max_retries = args.retries
    max_threads = args.threads

    # Read input file
    input_df = pd.read_csv(input_file, sep="\t", header=None,
                           names=["alignment", "treefile"])
    if debug:
        input_df = input_df.head(5)
        max_threads = min(5, max_threads)
        logger.info(
            "Debug mode enabled. Running on first 5 rows of the input file")

    # Create output directory
    if output_dir is None:
        output_dir = os.path.dirname(input_file)
    os.makedirs(output_dir, exist_ok=True)
    # get absolute path
    output_dir = os.path.abspath(output_dir)

    # Run HyPhy RELAX
    args = [(row.alignment, row.treefile, output_dir, path_to_relax_script, max_retries)
            for row in input_df.itertuples()]
    with mp.Pool(max_threads) as pool:
        outfiles = list(tqdm(pool.imap(run_hyphy_relax_single, args), total=len(
            args), desc="Running HyPhy RELAX in parallel"))

    return outfiles


def run_hyphy_general_single(args):
    """
    Run a custom HyPhy script on a single alignment and treefile
    """
    alignment, treefile, output_dir, path_to_hyphy_script, general_arguments, max_retries = args

    # Create output file
    output_file = os.path.join(output_dir, os.path.basename(
        alignment).replace(".aln.fna", ".hyphy.json"))
    
    # Define command
    command = [
        "hyphy", path_to_hyphy_script,
        "--tree", treefile,
        "--alignment", alignment,
        "--output", output_file,
    ]
    if general_arguments:
        command += general_arguments.split(" ")

    for i in range(max_retries):
        try:
            # Run command and capture stdout and stderr
            logger.info(f"Running command: {' '.join(command)}")

            process = subprocess.Popen(
                command,
                # shell=True,  # if shell is True, hyphy seems to get into interactive mode instead of running the command
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # Log stdout and stderr to a file specific to this run
            log_file = output_file.replace(".json", ".log")
            with open(log_file, "w") as f:
                f.write(f"{'--' * 30}\nDumping stdout:\n{'--' * 30}\n{stdout.decode()
                                                                      }\n\n\n{
                    '--' * 30}\nDumping stderr\n{
                    '--' * 30}\n{stderr.decode()}")

            # If process returns non-zero exit code, log error and retry if max_retries not reached
            if process.returncode != 0:
                logger.error(
                    f"Error running HyPhy general script for {alignment}: {stderr.decode()}")
                # retry
                if i < max_retries - 1:
                    logger.info(
                        f"Retrying {i+1}/{max_retries} for {alignment}")
                    continue
                elif i == max_retries - 1:
                    logger.error(
                        f"Failed to run HyPhy general script for {alignment} after {max_retries} retries")
                    return None
            else:
                logger.info(f"HyPhy general script run for {
                            alignment} completed successfully and output written to {output_file}")
                return output_file
        except Exception as e:
            logger.error(f"Error running HyPhy general script for {alignment}: {e}")
            logger.error(traceback.format_exc())
            return None
        
    return output_file


def run_hyphy_general_parallel(args):
    """
    Run a custom HyPhy script in parallel on a list of alignments and treefiles
    """

    input_file = args.input_file
    output_dir = args.output_dir
    path_to_hyphy_script = args.path_to_hyphy_script
    general_arguments = args.general_arguments
    debug = args.debug
    max_retries = args.retries
    max_threads = args.threads

    # Read input file
    input_df = pd.read_csv(input_file, sep="\t", header=None,
                           names=["alignment", "treefile"])
    if debug:
        input_df = input_df.head(5)
        max_threads = min(5, max_threads)
        logger.info(
            "Debug mode enabled. Running on first 5 rows of the input file")

    # Create output directory
    if output_dir is None:
        output_dir = os.path.dirname(input_file)
    os.makedirs(output_dir, exist_ok=True)
    # get absolute path
    output_dir = os.path.abspath(output_dir)

    # Run HyPhy general script
    args = [(row.alignment, row.treefile, output_dir, path_to_hyphy_script, general_arguments, max_retries)
            for row in input_df.itertuples()]
    with mp.Pool(max_threads) as pool:
        outfiles = list(tqdm(pool.imap(run_hyphy_general_single, args), total=len(
            args), desc="Running HyPhy general script in parallel"))

    return outfiles


if __name__ == "__main__":
    import time
    start = time.time()

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Run HyPhy programs on a list of alignments and corresponding pruned genome trees')

    parser.add_argument('--input-file', '-i',
                        required=True, type=str,
                        help='[Required] File containing filepaths to codon alignment files and treefiles as a tab-separated list of two columns. \
                        The first column should contain the path to the codon alignment file and the second column should contain the path to the treefile.\
                        Example: `og_id1.aln.fna\tog_id1.treefile`')
    parser.add_argument('--threads', '-t',
                        type=int, default=None,
                        help='[Optional, Default=90 perc of available threads] Number of threads to use for parallel processing')
    parser.add_argument('--output-dir', '-o',
                        required=False, type=str, default=None,
                        help='[Optional, Default=Same directory as input file] Directory where output files will be written')
    parser.add_argument('--hyphy-algorithm', '-a',
                        required=True, type=str,
                        help='[Required] HyPhy algorithm to run. Options: `busted`, `relax`, `general`')
    parser.add_argument('--path-to-busted-script', '-pb',
                        required=False, type=str, default="/root/bin/hyphy-analyses/BUSTED-SR/BUSTED-SR.bf",
                        help='[Conditionally required, Default="/root/bin/hyphy-analyses/BUSTED-SR/BUSTED-SR.bf"] Path to the HyPhy BUSTED script if running BUSTED')
    parser.add_argument('--path-to-relax-script', '-pr',
                        required=False, type=str, default="relax",
                        help='[Conditionally required, Default="relax"] Path to the HyPhy RELAX script if running RELAX')
    parser.add_argument('--path-to-hyphy-scipt', '-ph',
                        required=False, type=str, 
                        help='[Conditionally required] Path to the HyPhy script if running a HyPhy algorithm in `general`')
    parser.add_argument('--general-arguments', '-ga',
                        required=False, type=str, 
                        help='[Optional] General arguments to pass to the HyPhy script')
    parser.add_argument('--debug', '-d',
                        action='store_true',
                        help='[Optional] Enable debug mode: runs the script only on first 5 rows of the input file')
    parser.add_argument('--retries', '-r',
                        type=int, default=5,
                        help='[Optional, Default=5] Number of times to retry running a failed HyPhy BUSTED run')
    parser.add_argument('--two-rate-classes', '-r2',
                        action='store_true',
                        help='[Optional] Use two rate classes each for synonymous and non-synonymous sites in HyPhy BUSTED')

    args = parser.parse_args()

    # if threads is not provided, use 90% of available threads
    if args.threads is None:
        max_threads = int(mp.cpu_count() * 0.9)
    else:
        max_threads = args.threads

    if args.debug:
        # set log level to debug
        logger.remove()
        # now add a new handler with the new configuration for stdout and stderr
        logger.add(sys.stdout, level="DEBUG",
                   format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | <level>{message}</level>")

    # log args
    logger.info(f"Arguments: {args}")

    # # Run HyPhy BUSTED
    # all_outfiles = run_hyphy_busted_parallel(
    #     args, max_threads)
    # logger.info(f"Results written to {args.output_dir}")

    if args.hyphy_algorithm == "busted":
        all_outfiles = run_hyphy_busted_parallel(
            args)
    elif args.hyphy_algorithm == "relax":
        all_outfiles = run_hyphy_relax_parallel(
            args)
    elif args.hyphy_algorithm == "general":
        raise NotImplementedError(
            "Running HyPhy with a custom script is not yet implemented")
        # all_outfiles = run_hyphy_general_parallel(
        #     args, max_threads)
    else:
        raise ValueError(
            f"Invalid HyPhy algorithm: {args.hyphy_algorithm}. Options: `busted`, `relax`")
    logger.info(f"Results for {args.hyphy_algorithm.upper()} written to {args.output_dir}")

    # log how many successful runs
    successful_runs = [f for f in all_outfiles if f is not None]
    logger.info(
        f"Successful runs: {len(successful_runs)} out of {len(all_outfiles)}")

    end = time.time()
    logger.info(f"Time taken: {end - start} seconds")
