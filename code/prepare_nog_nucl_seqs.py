from email.policy import default
import os
import json
from Bio import SeqIO
from multiprocessing import Pool
from loguru import logger
from collections import defaultdict
import pandas as pd
import traceback
from tqdm import tqdm


# Extract the nucleotide sequences for each accession_id, in parallel, from the genome files
def extract_nucl_seq(args):
    accession_id, locus_tag_features, genome_dir, output_dir = args
    try:
        # Create a mapping from seqid to locus_tag and features
        seqid_to_locus_tag = defaultdict(list)
        for locus_tag, features in locus_tag_features.items():
            seqid_to_locus_tag[features['seqid']].append((locus_tag, features))

        genome_filepath = os.path.join(genome_dir, f"{accession_id}.fna")
        if not os.path.exists(genome_filepath):
            # raise error
            raise FileNotFoundError(
                f"Genome file not found for {accession_id}. Make sure it exists as {genome_filepath}")
        # Read the genome file and extract sequences
        with open(genome_filepath, 'r') as f:
            genome_records = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))

        # Extract the sequences for each locus_tag
        for seqid, locus_tags in seqid_to_locus_tag.items():
            if seqid in genome_records:
                seqid_record = genome_records[seqid]
                for locus_tag, features in locus_tags:
                    output_filepath = os.path.join(
                        output_dir, f"{features['taxon_id']}.{locus_tag}.fna")
                    start, end, strand = features['start'], features['end'], features['strand']
                    gene_seq = seqid_record.seq[start-1:end]
                    if strand == '-':
                        gene_seq = gene_seq.reverse_complement()
                    with open(output_filepath, 'w') as out_f:
                        out_f.write(f">{features['taxon_id']}\n{gene_seq}\n")
            else:
                raise KeyError(
                    f"SEQID {seqid} not found in the genome file for {accession_id}")

        return None
    except Exception as e:
        logger.error(f"Error processing {accession_id}: {e}")
        traceback.print_exc()
        return e

        return None
    except Exception as e:
        logger.error(f"Error processing {accession_id}: {e}")
        traceback.print_exc()
        return e


def extract_nog_nucl_seqs(gene_family_file, taxon_locus_json, genome_dir, output_dir, threads):
    # Load the gene family mapping
    og_members_df = pd.read_csv(gene_family_file, sep='\t', header=None,
                                names=['nog', 'members'],
                                index_col='nog')
    og_members_dict = og_members_df['members'].apply(
        lambda x: x.split(',')).to_dict()

    # Load the taxon-locus_tag mapping
    with open(taxon_locus_json, 'r') as f:
        taxon_locus_tag_dict = json.load(f)
    logger.info(f"Loaded files {gene_family_file} and {taxon_locus_json}")

    # Create the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Output directory created at {output_dir}")

    # Flatten the taxon-locus_tag mapping for faster processing
    flattened_mapping = {
        (features['accession_id'], locus_tag): {**features, 'taxon_id': taxon_id}
        for taxon_id, locus_tag_features in taxon_locus_tag_dict.items()
        for locus_tag, features in locus_tag_features.items()
    }
    logger.info("Flattened the taxon-locus_tag mapping")

    # Group by accession_id
    accession_id_to_locus_tag_features = defaultdict(dict)
    for (accession_id, locus_tag), features in flattened_mapping.items():
        accession_id_to_locus_tag_features[accession_id][locus_tag] = features
    logger.info("Grouped the locus tags by accession_id")
    logger.info(f"Number of accession_ids: {len(accession_id_to_locus_tag_features)}. \
                 These are: {', '.join(sorted(list(accession_id_to_locus_tag_features.keys())))}")

    # Create the temporary directory to store the gene sequences
    tmp_genes_dir = os.path.join(output_dir, 'genes')
    if not os.path.exists(tmp_genes_dir):
        os.makedirs(tmp_genes_dir)
        logger.info(f"Gene output directory created at {
                    tmp_genes_dir}")

    # Create the arguments for the parallel processing
    args = [(accession_id, locus_tag_features, genome_dir, tmp_genes_dir)
            for accession_id, locus_tag_features in accession_id_to_locus_tag_features.items()]

    # Run the parallel processing
    with Pool(threads) as p:
        for pool_return in tqdm(p.imap_unordered(extract_nucl_seq, args), total=len(args), desc="Extracting nucleotide sequences"):
            if pool_return is not None:
                logger.error(f"Error in processing: {pool_return}")

    # # now all the taxon_id.locus_tag values are extracted and written to the output directory
    # now for each NOG, combine the sequences and write to a file
    for nog, members in tqdm(og_members_dict.items(), desc="Writing NOG FASTN files"):
        nog_seq_records = []
        for member in members:
            try:
                taxon_id, locus_tag = member.split('.', 1)
                seq_file = os.path.join(
                    output_dir, 'genes', f"{taxon_id}.{locus_tag}.fna")
                logger.info(f"Processing {taxon_id}.{locus_tag}")
                if os.path.exists(seq_file):
                    with open(seq_file, 'r') as f:
                        record = SeqIO.read(f, 'fasta')
                        nog_seq_records.append(record)
                else:
                    raise FileNotFoundError(
                        f"File not found for {taxon_id}.{locus_tag}")
            except Exception as e:
                logger.error(f"Error processing {member}: {e}")
                traceback.print_exc()

        # write the sequences to a file
        nog_output_file = os.path.join(output_dir, f"{nog}.fna")
        with open(nog_output_file, 'w') as f:
            SeqIO.write(nog_seq_records, f, 'fasta')


if __name__ == "__main__":
    import argparse
    import time
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--og-members-map", "-m",
        type=str, required=True,
        help=f"Path to the gene family mapping file. \
            Should be TSV file with one column as the gene family ID and the other is a comma-separated list of gene members, each in format taxon_id.locus_tag")
    parser.add_argument(
        "--gene-features-json", "-j",
        type=str, required=True,
        help=f"Path to the JSON file containing the taxon-locus_tag-features mapping")
    parser.add_argument(
        "--genome-dir", "-g",
        type=str, required=True,
        help="Path to the directory containing the genome sequences in FASTA format")
    parser.add_argument(
        "--output-dir", "-o",
        type=str, required=True,
        help="Path to the directory to save the extracted sequences")
    parser.add_argument(
        "--threads", "-t",
        type=int, default=100,
        help="Number of threads to use for parallel processing")
    args = parser.parse_args()

    # log arguments
    logger.info(f"Arguments: {args}")

    extract_nog_nucl_seqs(args.og_members_map, args.gene_features_json,
                          args.genome_dir, args.output_dir, args.threads)

    logger.info(f"Time taken: {time.time() - start_time} seconds")
