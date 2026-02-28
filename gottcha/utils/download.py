#!/usr/bin/env python3

#This script allows the user to pull the existing gottcha database
import requests
import sys
import os
import tarfile
import hashlib
from tqdm import *
import argparse
from gottcha import GOTTCHA_DB_LATEST
import logging

def calculate_md5(file_path, chunk_size=8192):
    md5 = hashlib.md5()
    with open(file_path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b''):
            md5.update(chunk)
    return md5.hexdigest()


def parse_params(args):
    parser = argparse.ArgumentParser(prog='download.py', description="""This script will pull the latest version of the Gottcha2 database.""")

    parser.add_argument('-u', '--url', default=GOTTCHA_DB_LATEST,
                    help='specify a URL to pull from (will override the default)')
    parser.add_argument('-r', '--rank', default='species',
                    help='taxonomic rank of the database (superkingdom, phylum, class, order, famiily, genus, species)')
    return parser.parse_args(args)


def download_db(argvs):
    if os.path.isdir('database'):
        sys.exit('Please make sure a database directory does not exist.')

    os.mkdir('database')

    download_url = argvs.url
    archive_name = os.path.basename(download_url)

    logging.info(f"Downloading GOTTCHA2 database from {download_url}...")

    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        with open(archive_name, 'wb') as f:
            total_size = int(r.headers.get('Content-Length', 0))
            pbar = tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024)
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
            pbar.close()

    logging.info(f"Downloading GOTTCHA2 database MD5 checksum from {download_url}.md5...")

    with requests.get(f"{download_url}.md5", stream=True) as r:
        r.raise_for_status()
        with open(f"{archive_name}.md5", 'wb') as f:
            total_size = int(r.headers.get('Content-Length', 0))
            pbar = tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024)
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
            pbar.close()

    # check md5
    logging.info("Verifying MD5 checksum...")

    with open(f"{archive_name}.md5", 'r') as f:
        expected_md5 = f.read().strip().split()[0]
    actual_md5 = calculate_md5(archive_name)
    if actual_md5 != expected_md5:
        os.remove(archive_name)
        os.remove(f"{archive_name}.md5")
        sys.exit("MD5 checksum does not match expected value. Download may be corrupted. Please try again.")

    logging.info("MD5 checksum verified successfully.")

    logging.info(f"Extracting GOTTCHA2 database from {archive_name}...")
    with tarfile.open(archive_name) as tar:
        tar.extractall('database')
    logging.info("Database extraction completed.")
    os.remove(archive_name)
    os.remove(f"{archive_name}.md5")
    logging.info("Temporary files removed. Database is ready to use.")

def main(args):
    argvs = parse_params(args)
    download_db(argvs)


if __name__ == '__main__':
    main(sys.argv[1:])
