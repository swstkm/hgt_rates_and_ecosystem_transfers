import os
import re
import pandas as pd
from Bio import Entrez
import gffutils
from tqdm import tqdm
import argparse
from loguru import logger
import time
import datetime
import subprocess
import zipfile
from collections import defaultdict
import traceback

import xml.etree.ElementTree as ET

# This script downloads genome assemblies from NCBI datasets program.
# It takes a file with gene IDs in the format `taxon_id.locus_tag` and downloads the GFF files for the genome assemblies.
# This is performed in two rounds for each taxon ID:
# - In the first round, we only process the first locus tag of each locus tag prefix (e.g. `BACINT_00001 for BACINT_00001, BACINT_00002, etc.)
# - In the second round, we process the rest of the locus tags for each prefix
# At every round, we check if the locus tag is already found in the GFF files downloaded previously for the same taxon ID.
# If the locus tag is found, we skip downloading the GFF file for that locus tag. If not found, we try downloading the GFF file for the locus tag.
# The script writes the number of locus tags found in the GFF files for each taxon ID to a file `found_locus_tags_count.tsv`.


def get_genome_assembly(taxon_id, locus_tag, n_retries, api_key, bin_dir):
    '''
    Get the genome assembly accession for a given locus tag.
    This function uses the NCBI datasets program to fetch the genome assembly accession for a given locus tag.
    It first finds the protein ID from the locus tag, then finds the bioproject accession from the protein record.
    Then it finds the biosample ID from the bioproject record, and finally fetches the assembly accession from the biosample record using the `datasets` program.
    '''
    time_delay = 0.5  # seconds
    overall_retries = n_retries
    while overall_retries > 0:
        try:
            # find the protein ID from the locus tag
            retries = n_retries
            while retries > 0:
                handle = Entrez.esearch(db='protein', term=f'{
                                        locus_tag}')
                record = Entrez.read(handle)
                handle.close()
                if not record['IdList']:
                    # retry
                    retries -= 1
                    if retries == 0:
                        return None
                    logger.info(f"Retrying search for {
                                locus_tag}, retries left: {retries}")
                    time.sleep(time_delay)
                else:
                    break
            protein_id = record['IdList'][0]

            # find the bioproject accession from the protein record
            retries = n_retries
            while retries > 0:
                handle = Entrez.efetch(db='protein', id=protein_id,
                                       rettype='gb', retmode='xml')
                record = Entrez.read(handle)
                handle.close()

                bioproject_acc = record[0].get('GBSeq_project', None)
                if not bioproject_acc:
                    # retry
                    retries -= 1
                    if retries == 0:
                        return None
                    logger.info(f"Retrying fetch for protein ID {protein_id} with locus tag {
                                locus_tag}, retries left: {retries}")
                    time.sleep(time_delay)
                else:
                    break

            # find the biosample ID from the efetch bioproject xml records
            retries = n_retries
            while retries > 0:
                handle = Entrez.efetch(db='bioproject', id=f'{
                    bioproject_acc}', retmax=10)
                bioproject_record = handle.read()
                time.sleep(time_delay)
                handle.close()

                root = ET.fromstring(bioproject_record)
                biosample_id = None
                our_locus_tag_prefix = locus_tag.split('_')[0]
                # list how many LocusTagPrefix elements are there in total
                all_locus_tag_prefixes = [
                    elem.text for elem in root.iter('LocusTagPrefix')]
                logger.info(f"Bioproject accession {bioproject_acc} has {
                            len(all_locus_tag_prefixes)} LocusTagPrefix elements: {all_locus_tag_prefixes}")

                for elem in root.iter('LocusTagPrefix'):
                    if elem.text == our_locus_tag_prefix:
                        biosample_id = elem.attrib.get('biosample_id')
                        logger.info(f"Found biosample ID {biosample_id} for bioproject accession {
                            bioproject_acc} with locus tag {
                            locus_tag}")
                        # if assembly accession exists as attribute to this tag, return it directly
                        assembly_accession = elem.attrib.get('assembly_id')
                        if assembly_accession:
                            return assembly_accession
                        break

                if not biosample_id:
                    # retry
                    retries -= 1
                    if retries == 0:
                        return None
                    logger.info(f"Retrying fetch for bioproject accession {
                        bioproject_acc} with locus tag {
                        locus_tag}, retries left: {retries}. Bioproject: {bioproject_acc}")
                    time.sleep(time_delay)
                else:
                    break

            # find the assembly accession from the biosample record using `datasets`
            # basically datasets can be queried for summary of bioproject ID, and that can give us the assembly accession for the biosample ID
            dataformat_fields = "--fields accession,assminfo-name,source_database,assminfo-refseq-category,annotinfo-featcount-gene-total,assminfo-biosample-accession"
            postprocess_cmd = f"grep -w {biosample_id} | cut -f 1 | tail -n 1"
            as_json_lines_or_not = "--as-json-lines"
            api_key_or_not = f"--api-key {api_key}" if api_key else ""
            retries = n_retries
            assembly_accession = None
            while retries > 0:
                try:
                    # datasets summary call
                    datasets_call_output = subprocess.run(
                        f"{bin_dir}/datasets summary genome accession {bioproject_acc} {as_json_lines_or_not} {
                            api_key_or_not} | {bin_dir}/dataformat tsv genome {dataformat_fields} | {postprocess_cmd}",
                        shell=True,
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE
                    )
                    logger.info(f"datasets call output: {
                                datasets_call_output.stdout.decode()}")
                    assembly_accession = datasets_call_output.stdout.decode().strip()
                    if not assembly_accession:
                        # retry
                        retries -= 1
                        if retries == 0:
                            return None
                        logger.info(f"Retrying fetch for bioproject accession {
                            bioproject_acc} with locus tag {
                            locus_tag}, retries left: {retries}. Bioproject: {bioproject_acc}")
                        time.sleep(time_delay)
                    else:
                        return assembly_accession
                except subprocess.CalledProcessError as e:
                    logger.error(f"Error fetching genome assembly for {
                        # add info and traceback
                        taxon_id}.{locus_tag}: {e}")
                    logger.error(traceback.format_exc())
                    return None

            return None
        except Exception as e:
            logger.error(f"Error fetching genome assembly for {
                taxon_id}.{locus_tag}: {e}")
            # logger.error(traceback.format_exc())
            overall_retries -= 1
            if overall_retries == 0:
                logger.error(f"Failed to fetch genome assembly for {taxon_id}.{
                             locus_tag} after multiple retries")
                return None
            logger.info(f"Retrying fetch for {taxon_id}.{
                        locus_tag} after raised exception, retries left: {overall_retries}")


def download_gff_file(taxon_id, assembly_acc, output_dir, n_retries, api_key, bin_dir):
    # use datasets program to download the GFF file
    output_path = os.path.join(output_dir, f"{taxon_id}_{
                               assembly_acc}_genomic.gff")
    api_key_or_not = f"--api-key {api_key}" if api_key else ""
    progress_bar_or_not = "--no-progressbar"
    while n_retries > 0:
        try:
            subprocess.run(
                f"{bin_dir}/datasets download genome accession {
                    assembly_acc} --include gff3 --filename {output_path}.zip {progress_bar_or_not} {api_key_or_not}",
                shell=True,
                check=True
            )
            with zipfile.ZipFile(f"{output_path}.zip", 'r') as zip_ref:
                if len(zip_ref.namelist()) == 0:
                    logger.error(f"no files found in zip archive for assembly accession: {
                        assembly_acc}, taxon ID: {taxon_id}")
                    break
                for file_i in zip_ref.namelist():
                    if file_i.endswith(".json") or file_i.endswith(".jsonl") or file_i.endswith("README.md") or file_i.endswith("md5sum.txt"):
                        continue

                    # extract other files
                    logger.info(f"extracting {file_i} for assembly accession: {
                        assembly_acc}, taxon ID: {taxon_id} from zip archive {output_path}.zip")
                    file_i_outpath = None
                    with zip_ref.open(file_i) as f:
                        file_i_basename = os.path.basename(file_i)
                        # check if file_i_basename starts with assembly_acc
                        if file_i_basename.startswith(assembly_acc):
                            file_i_outpath = f"{
                                output_dir}/{taxon_id}_{file_i_basename}"
                        else:
                            file_i_outpath = f"{
                                output_dir}/{taxon_id}_{assembly_acc}_{file_i_basename}"
                        with open(f"{output_dir}/{taxon_id}_{assembly_acc}_{file_i_basename}", 'wb') as f_out:
                            f_out.write(f.read())
                    logger.info(f"extracted file {file_i} for assembly accession: {
                        assembly_acc}, taxon ID: {taxon_id} to {file_i_outpath}")
            os.remove(f"{output_path}.zip")

            break
        except subprocess.CalledProcessError as e:
            logger.error(f"Error downloading GFF file for {
                         taxon_id}.{assembly_acc}: {e}")
            logger.error(traceback.format_exc())
            n_retries -= 1
            if n_retries == 0:
                logger.error(f"Failed to download GFF file for {taxon_id}.{
                             assembly_acc} after multiple retries")
                return None
            logger.info(f"Retrying download for {taxon_id}.{
                        assembly_acc}, retries left: {n_retries}")

    return output_path


def extract_locus_tags_from_gff(gff_path):
    try:
        gff_db = gffutils.create_db(
            gff_path, ':memory:', merge_strategy='merge')
        all_found_tags = set()
        for feature in gff_db.features_of_type('gene'):
            if 'locus_tag' in feature.attributes:
                all_found_tags.add(feature.attributes['locus_tag'][0])
            if 'old_locus_tag' in feature.attributes:
                all_found_tags.add(feature.attributes['old_locus_tag'][0])
        return all_found_tags
    except Exception as e:
        logger.error(f"Error checking locus tags in GFF file {gff_path}: {e}")
        logger.error(traceback.format_exc())
        return set()


def rearrange_locus_tags(locus_tags):
    prefix2tags = defaultdict(list)

    for locus_tag in locus_tags:
        prefix = locus_tag.split('_')[0]
        prefix2tags[prefix].append(locus_tag)

    # we make two lists, one with the first locus tag of each prefix and
    # the other with the rest of the locus tags
    round1_locus_tags = [tags[0] for tags in prefix2tags.values()]
    round2_locus_tags = [tag for tags in prefix2tags.values()
                         for tag in tags[1:]]

    return round1_locus_tags, round2_locus_tags


def download_gff_of_locus_tags(gene_ids_file, output_dir, debug,
                               n_retries, bin_dir, api_key):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    gene_ids_df = pd.read_csv(gene_ids_file, header=None, names=['gene_id'])
    # sort df
    gene_ids_df = gene_ids_df.sort_values('gene_id')
    # split gene_id into taxon_id and locus_tag
    gene_ids_df['taxon_id'] = gene_ids_df['gene_id'].apply(
        lambda x: x.split('.')[0])
    gene_ids_df['locus_tag'] = gene_ids_df['gene_id'].apply(
        lambda x: x.split('.')[1])

    if debug:  # run only on a random sample of 100 gene IDs
        gene_ids_df = gene_ids_df.sample(1000)

    taxa2genes = gene_ids_df.groupby('taxon_id')

    def process_locus_tags(taxon_id, locus_tags, output_dir, n_retries, api_key, bin_dir,
                           all_gff_locus_tags=set()):
        '''
        Process locus tags for a given taxon_id. Download the GFF file for the assembly
        '''

        now_found_tags = set()
        already_found_tags = set()
        # make a copy of locus_tags
        searched_found_locus_tags = locus_tags.copy()

        for locus_tag in locus_tags:
            # check if locus tag is already found in the GFF files and if not, download the GFF file
            if locus_tag in all_gff_locus_tags or locus_tag in now_found_tags:
                already_found_tags.add(locus_tag)
                continue
            else:  # locus tag not found in any of the GFF files
                assembly_acc = get_genome_assembly(
                    taxon_id, locus_tag, n_retries, api_key, bin_dir)
                if not assembly_acc:  # assembly not found
                    logger.error(f"Assembly not found for {
                                 taxon_id}.{locus_tag}")
                    # remove it from searched_found_locus_tags
                    searched_found_locus_tags.remove(locus_tag)
                    continue
                gff_path = download_gff_file(
                    taxon_id, assembly_acc, output_dir, n_retries, api_key, bin_dir)
                if not gff_path:
                    # remove it from searched_found_locus_tags
                    searched_found_locus_tags.remove(locus_tag)
                    continue
                # extract locus tags from the GFF file
                now_found_tags.update(
                    # update with tags found in the GFF file
                    extract_locus_tags_from_gff(gff_path))
        logger.info(f"This round: Taxon {taxon_id}: Searched for {
                    len(locus_tags)} locus tags, found {
                        len(searched_found_locus_tags)} of them. Found {len(
                            now_found_tags)} locus tags additionally, in GFF files. {
                                len(already_found_tags)} locus tags already found in previous GFF files")

        return now_found_tags, searched_found_locus_tags

    # set of all locus tags found out of the searched locus tags
    searched_found_all_tags = set()

    # we make a set of taxon_ids for which we found at least 4 searched locus tags in the GFF files
    taxa_with_found_locus_tags = set()
    with open(f"{output_dir}/found_locus_tags_count.tsv", 'w') as f:
        f.write("taxon_id\tlocus_tags_found_total\tlocus_tags_searched_and_found\n")

    # iterate over each taxon_id and process locus tags
    for taxon_id, group in tqdm(taxa2genes):
        locus_tags = set(group['locus_tag'])
        if len(locus_tags) < 4:  # need a minimum of 4 locus tags
            logger.info(f"Skip processing {taxon_id} with only {
                        len(locus_tags)} locus tags")
            continue

        # we split the locus tags into two rounds
        # in the first round, we only process the first locus tag of each locus tag prefix
        # this is because of a single locus tag exists with the prefix, it is likely that others do too, in the same gff file
        # in the second round, we process the rest of the locus tags for each prefix
        round1_locus_tags, round2_locus_tags = rearrange_locus_tags(locus_tags)
        logger.info(f"Taxon {taxon_id} with {len(locus_tags)} locus tags. Round 1: {
                    len(round1_locus_tags)}, Round 2: {len(round2_locus_tags)}")

        # Process round1_locus_tags
        round1_n_retries = 3  # lower number of retries for round 1
        all_found_tags_round1, subset_found_tags_round1 = process_locus_tags(
            taxon_id, round1_locus_tags, output_dir, round1_n_retries, api_key, bin_dir)
        logger.info(f"Round 1: Taxon {taxon_id}: Found {
                    len(all_found_tags_round1)} locus tags across GFF files while searching for {len(round1_locus_tags)}")

        # if round1 returns 0 locus tags, we skip round2
        if not all_found_tags_round1:
            logger.info(f"Skip round 2 for Taxon {
                        taxon_id}. Here are examples of locus tags that were not found in round 1: {round1_locus_tags[:5]}")
            continue

        # Process round2_locus_tags
        all_found_tags_round2, subset_found_tags_round2 = process_locus_tags(
            taxon_id, round2_locus_tags, output_dir, n_retries, api_key, bin_dir,
            all_gff_locus_tags=all_found_tags_round1)
        logger.info(f"Found {len(all_found_tags_round2)} locus tags in round 2 for {
                    taxon_id}, out of {len(locus_tags)}")

        all_found_tags = all_found_tags_round1.union(all_found_tags_round2)
        searched_found_tags = set(subset_found_tags_round1).union(
            set(subset_found_tags_round2))
        logger.info(f"All rounds : Taxon {taxon_id}: Found {len(all_found_tags)} locus tags in total. Searched {
                    len(locus_tags)}, found {len(searched_found_tags)}")

        if not all_found_tags:
            logger.info(
                f"Could not find locus tags in GFF files for {taxon_id}")
        else:
            logger.info(f"Found {len(all_found_tags)} locus tags in GFF files for {
                        taxon_id}, out of {len(locus_tags)}")
            # write to a file, the number of locus tags found for each taxon_id
            with open(f"{output_dir}/found_locus_tags_count.tsv", 'a') as f:
                f.write(f"{taxon_id}\t{len(all_found_tags)}\t{
                        len(searched_found_tags)}\n")

        searched_found_all_tags.update(searched_found_tags)
        if len(all_found_tags) >= 4:
            taxa_with_found_locus_tags.add(taxon_id)

    logger.info(f"Searched for {gene_ids_df['locus_tag'].nunique(
    )} locus tags in total, found {len(searched_found_all_tags)}. Downloaded all associated GFF files")
    # write this list of searched and found locus tags to a file
    with open(f"{output_dir}/found_locus_tags.txt", 'w') as f:
        f.write("\n".join(searched_found_all_tags) + "\n")

    logger.info(f"Done processing {gene_ids_file}")


if __name__ == "__main__":

    # time
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description='Download genome assemblies from NCBI and check if locus tags are present in the GFF files')
    parser.add_argument('--gene_ids_file', '-i', type=str,
                        help='File with a gene ID per line in the format `taxon_id.locus_tag`')
    parser.add_argument('--output_dir', '-o', type=str,
                        help='Output directory to save the GFF files')
    parser.add_argument('--debug', '-d', action='store_true',
                        help='Debug mode: run only on a random sample of 1000 gene IDs')
    parser.add_argument('--n-retries', '-r', type=int, default=5,
                        help='Number of retries to download the GFF file')
    parser.add_argument('--credentials_file', '-c', type=str, required=True,
                        help='Path to a file containing the email and API key for NCBI datasets program')
    parser.add_argument('--bin-dir', '-b', type=str, default='/root/mambaforge/envs/hgt_analyses/bin/',
                        help='Path to the directory where the `datasets` and `dataformat` binaries are located')
    args = parser.parse_args()

    if os.path.exists(args.credentials_file):
        with open(args.credentials_file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                raise ValueError(
                    "Credentials file must contain at least two lines: email on the first line and API key on the second line")
            Entrez.email = lines[0].strip()
            api_key = lines[1].strip()
    else:
        raise FileNotFoundError(
            f"Credentials file {args.credentials_file} not found")

    download_gff_of_locus_tags(
        args.gene_ids_file, args.output_dir, args.debug, args.n_retries, args.bin_dir, api_key)

    # log time taken using datetime
    logger.info(f"Time taken: {datetime.timedelta(
        seconds=(time.time() - start_time))}")
