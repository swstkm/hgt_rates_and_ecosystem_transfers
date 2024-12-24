# run IQtree on the MSA files output by MAFFT (via OrthoFinder, in previous step)
# we run IQtree in parallel on all the MSA files
import os
import multiprocessing
import subprocess
from loguru import logger
import random
import sys


def run_iqtree(msa_filepath, iqtree_bin, output_trees_dir, substitution_model, threads_max, threads_per_process=1):
    output_tree_prefix = f'{
        output_trees_dir}/{os.path.basename(msa_filepath).replace(".fa", "")}'
    # run iqtree
    cmd = f'{
        iqtree_bin} -s {msa_filepath} -m {substitution_model} -pre {output_tree_prefix}, -quiet -T {threads_max} -nt {threads_per_process}'
    logger.info(f'Running IQtree on {
                msa_filepath}, to write to {output_tree_prefix} as output prefix.')
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        return f'Error running IQtree on {msa_filepath}: {e}'
    logger.info(f'Finished running IQtree on {msa_filepath}')
    return None


def run_iqtree_wrapper(args):
    """
    Wrapper function to run IQtree with multiple arguments, for the multiprocessing.Pool.imap_unordered function.
    """
    return run_iqtree(*args)


def get_seq_length(msa_filepath):
    """
    Function to get the number of characters in the first sequence of the MSA file.
    """
    try:
        with open(msa_filepath, 'r') as f:
            # first check if the file is empty
            if os.stat(msa_filepath).st_size == 0:
                return 0
            # split file by >, and take the second element, which is the first sequence
            seq = f.read().split('>')[1]
            return len(seq)
    except Exception as e:
        logger.error(f'Error reading MSA file {msa_filepath}: {e}')
        raise e


if __name__ == '__main__':
    import time
    import argparse

    start_time = time.time()

    parser = argparse.ArgumentParser(
        description='Run IQtree on all MSA files in a directory. For large MSA files, run with multiple threads per process')

    parser.add_argument('-m', '--msa_dir', type=str, required=True,
                        help='Directory containing MSA files')
    parser.add_argument('-o', '--output_dir', type=str, default='./gene_trees/',
                        help='Output directory for gene trees')
    parser.add_argument('-T', '--threads-max', type=int,
                        default=100,
                        help='Maximum number of threads to use')
    parser.add_argument('-t', '--threads-per-process', type=int, default=5,
                        help='Number of threads per process to use, for large MSA files')
    parser.add_argument('-b', '--iqtree_bin', type=str, default='/root/bin/iqtree-2.2.2.6-Linux/bin/iqtree2',
                        help='Path to IQtree binary')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Enable debug mode (run on first 5 MSA files only)')
    parser.add_argument('-s', '--substitution_model', type=str, default='JTT',
                        help='Substitution model for IQtree')
    parser.add_argument('--filesize_cutoff', type=int, default=100,
                        help='Cutoff for number of sequences in MSA file to determine if it is a large or small file.\n \
                        Files with more than this cutoff will be run with multiple threads per process.')
    parser.add_argument('--min_seqs', type=int, default=20,
                        help='Minimum number of sequences in MSA file to run IQtree. Files with fewer sequences will be skipped.')
    parser.add_argument('--random', type=int, nargs='?', const=1000, default=None,
                        help='If given, pick this many MSA at random. \
                            If no number is provided, default to 1000. \
                            If this argument is not given, take all MSA above min_seqs in size.')
    parser.add_argument('--min_len', type=int, default=150,
                        help='Minimum length of longest sequence of MSA file (including gaps). Files with smaller genes will be skipped.')
    args = parser.parse_args()

    # remove output directory if it exists
    if os.path.exists(args.output_dir):
        os.system(f'rm -r {args.output_dir}')

    # create output directory new
    os.makedirs(args.output_dir, exist_ok=True)

    # set up logging
    logger.add(
        f'{args.output_dir}/iqtree_run_{time.strftime("%Y%m%d_%H%M%S")}.log')

    # check if IQtree binary exists and MSA directory exists
    if not os.path.exists(args.iqtree_bin):
        raise FileNotFoundError(f'IQtree binary not found: {args.iqtree_bin}')
    if not os.path.exists(args.msa_dir):
        raise FileNotFoundError(f'MSA directory not found: {args.msa_dir}')

    # get all MSA files (in FASTA format, ending with ".fa") in the directory specified by `args.msa_dir`
    msa_filepaths = [f for f in os.listdir(
        args.msa_dir) if f.endswith('.fa')]
    msa_filepaths = [os.path.join(args.msa_dir, f) for f in msa_filepaths]

    # first get the sizes (number of characters in the second line), of all the MSA files
    msa_lengths = [(f, get_seq_length(f)) for f in msa_filepaths]

    # remove all MSA files that have fewer than min_len characters in the first sequence
    short_seqs_to_remove = []
    for msa_filepath, length in msa_lengths:
        if length < args.min_len:
            short_seqs_to_remove.append(msa_filepath)
    msa_filepaths = [f for f in msa_filepaths if f not in short_seqs_to_remove]
    logger.info(
        f'Number of MSA files with more than {args.min_len} characters in seq: {len(msa_filepaths)}')

    # write this list of MSA files and their sequence counts to a file
    with open(f'../scrap/msa_files_sequence_counts.txt', 'w') as f:
        for msa_filepath, num_seqs in msa_lengths:
            f.write(f'{num_seqs}\t{msa_filepath}\n')

    # for all the MSA files, find the number of sequences they have by opening and reading the number of lines starting with ">".
    # Retain only MSA files that have more than `args.min_seqs` sequences
    small_files_to_remove = []
    msa_filepath_sizes = []
    for msa_filepath in msa_filepaths:
        with open(msa_filepath, 'r') as f:
            num_seqs = sum(1 for l in f if l.startswith('>'))
        if num_seqs < args.min_seqs:
            small_files_to_remove.append(msa_filepath)
        else:
            msa_filepath_sizes.append((msa_filepath, num_seqs))

    msa_filepaths = [
        f for f in msa_filepaths if f not in small_files_to_remove]
    logger.info(
        f'Number of MSA files with more than {args.min_seqs} sequences: {len(msa_filepaths)}')

    if args.random is not None:
        msa_filepath_sizes = random.sample(
            msa_filepath_sizes, min(args.random, len(msa_filepath_sizes)))
        logger.info(f'Picked {len(msa_filepath_sizes)} MSA files at random')

    # run IQtree on all MSA files
    if args.debug:
        msa_filepath_sizes = msa_filepath_sizes[:5]

    # we will use a hybrid approach, wherein large MSA files will have multiple threads per process
    # and small MSA files will have a single thread per process
    # large MSA files are those with more than 100 sequences
    small_msa_filepaths = [
        f for f in msa_filepath_sizes if f[1] <= args.filesize_cutoff]
    large_msa_filepaths = [
        f for f in msa_filepath_sizes if f[1] > args.filesize_cutoff]

    logger.info(
        f'Number of MSA files with <= {args.filesize_cutoff} sequences: {len(small_msa_filepaths)}, and > {args.filesize_cutoff} sequences: {len(large_msa_filepaths)}')
    # sys.exit(0)

    if small_msa_filepaths:
        logger.info(
            f'Running IQtree on {len(small_msa_filepaths)} MSA files each with <= 100 sequences')
        # run small MSA files first
        with multiprocessing.Pool(processes=args.threads_max) as pool:
            for pool_return in pool.imap_unordered(run_iqtree_wrapper, [(msa_filepath, args.iqtree_bin, args.output_dir, args.substitution_model, args.threads_max, 1)
                                                                        for msa_filepath, _ in small_msa_filepaths]):
                if pool_return is not None:
                    logger.error(f'Error running IQtree: {pool_return}')
    if large_msa_filepaths:
        logger.info(
            f'Running IQtree on {len(large_msa_filepaths)} MSA files each with > 100 sequences')
        # run large MSA files next
        with multiprocessing.Pool(processes=args.threads_max//args.threads_per_process) as pool:
            for pool_return in pool.imap_unordered(run_iqtree_wrapper, [(msa_filepath, args.iqtree_bin, args.output_dir, args.substitution_model, args.threads_max, args.threads_per_process)
                                                                        for msa_filepath, _ in large_msa_filepaths]):
                if pool_return is not None:
                    logger.error(f'Error running IQtree: {pool_return}')

    # log time taken
    time_taken = time.time() - start_time
    logger.info(f'Time taken to finish: {
                time.strftime("%H:%M:%S", time.gmtime(time_taken))}')
