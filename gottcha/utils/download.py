#!/usr/bin/env python3

#This script allows the user to pull the existing gottcha database
import requests
import sys
import os
import tarfile
import hashlib
import shutil
import subprocess
from urllib.parse import urlparse
from tqdm import *
import argparse
from gottcha import GOTTCHA_DB_FAST_LATEST, GOTTCHA_DB_STD_LATEST
import logging
from requests.exceptions import RequestException, SSLError


CHUNK_SIZE = 8192
REQUEST_TIMEOUT = (10, 60)


class DownloadError(Exception):
    pass

def calculate_sha256(file_path, chunk_size=8192):
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b''):
            sha256.update(chunk)
    return sha256.hexdigest()


def parse_params(args):
    parser = argparse.ArgumentParser(prog='download.py', description="""This script will pull the latest version of the Gottcha2 database.""")

    parser.add_argument('-u', '--url', required=False,
                    help='specify a URL to pull from (will override the default)')
    parser.add_argument('-d', '--database', default='fast', choices=['standard', 'fast'],
                    help='specify the type of database to download (standard or fast)')
    parser.add_argument('--ca-bundle', required=False,
                    help='path to a PEM certificate bundle to use for HTTPS verification')
    parser.add_argument('--no-verify-ssl', action='store_true',
                    help='disable HTTPS certificate verification (not recommended)')
    return parser.parse_args(args)


def get_archive_name(download_url):
    archive_name = os.path.basename(urlparse(download_url).path)
    if not archive_name:
        raise DownloadError(f"Unable to determine archive filename from URL: {download_url}")
    return archive_name


def get_verify_setting(argvs):
    if argvs.no_verify_ssl:
        logging.warning("HTTPS certificate verification is disabled. Use only on trusted networks.")
        requests.packages.urllib3.disable_warnings()
        return False

    if argvs.ca_bundle:
        if not os.path.isfile(argvs.ca_bundle):
            raise DownloadError(f"CA bundle does not exist: {argvs.ca_bundle}")
        return argvs.ca_bundle

    return True


def remove_if_exists(path):
    if os.path.exists(path):
        os.remove(path)


def response_content_length(response):
    try:
        return int(response.headers.get('Content-Length', 0))
    except (TypeError, ValueError):
        return 0


def download_file_with_requests(download_url, output_path, verify):
    with requests.get(
        download_url,
        stream=True,
        timeout=REQUEST_TIMEOUT,
        verify=verify,
    ) as r:
        r.raise_for_status()
        total_size = response_content_length(r)
        with open(output_path, 'wb') as f:
            pbar = tqdm(total=total_size, unit='B', unit_scale=True, unit_divisor=1024)
            for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
            pbar.close()


def download_file_with_curl(download_url, output_path, ca_bundle=None, insecure=False):
    curl_path = shutil.which('curl')
    if not curl_path:
        raise DownloadError("curl is not available for HTTPS certificate fallback")

    command = [
        curl_path,
        '--fail',
        '--location',
        '--show-error',
        '--progress-bar',
        '--output',
        output_path,
    ]

    if ca_bundle:
        command.extend(['--cacert', ca_bundle])

    if insecure:
        command.append('--insecure')

    command.append(download_url)

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        raise DownloadError(f"curl failed while downloading {download_url}") from exc


def download_file(download_url, output_path, verify, ca_bundle=None, insecure=False):
    try:
        download_file_with_requests(download_url, output_path, verify)
    except SSLError as exc:
        remove_if_exists(output_path)

        if verify is False:
            raise DownloadError(f"HTTPS download failed for {download_url}: {exc}") from exc

        logging.warning(
            "Python HTTPS certificate verification failed. Retrying with curl so "
            "the operating-system certificate store can be used."
        )
        try:
            download_file_with_curl(download_url, output_path, ca_bundle=ca_bundle)
        except DownloadError as curl_exc:
            raise DownloadError(
                "HTTPS certificate verification failed. Install the required CA "
                "certificate for Python, set REQUESTS_CA_BUNDLE/CURL_CA_BUNDLE, "
                "pass --ca-bundle, or retry with --no-verify-ssl if you accept "
                "the security risk."
            ) from curl_exc
    except RequestException as exc:
        remove_if_exists(output_path)
        raise DownloadError(f"Failed to download {download_url}: {exc}") from exc


def download_db(argvs):
    if os.path.isdir('database'):
        sys.exit('Please make sure a database directory does not exist.')

    if not argvs.url:
        if argvs.database == 'standard':
            download_url = GOTTCHA_DB_STD_LATEST
        elif argvs.database == 'fast':
            download_url = GOTTCHA_DB_FAST_LATEST
        else:
            sys.exit('Invalid database type specified. Please choose "standard" or "fast".')

    if argvs.url:
        download_url = argvs.url

    try:
        archive_name = get_archive_name(download_url)
        verify = get_verify_setting(argvs)
    except DownloadError as exc:
        sys.exit(str(exc))

    checksum_name = f"{archive_name}.sha256"

    try:
        os.mkdir('database')
    except OSError as exc:
        sys.exit(f"Failed to create database directory: {exc}")

    try:
        print(f"Downloading GOTTCHA2 database from {download_url}...")
        download_file(
            download_url,
            archive_name,
            verify,
            ca_bundle=argvs.ca_bundle,
            insecure=argvs.no_verify_ssl,
        )

        print(f"Downloading GOTTCHA2 database SHA256 checksum from {download_url}.sha256...")
        download_file(
            f"{download_url}.sha256",
            checksum_name,
            verify,
            ca_bundle=argvs.ca_bundle,
            insecure=argvs.no_verify_ssl,
        )

        print("Verifying SHA256 checksum...")

        with open(checksum_name, 'r') as f:
            expected_sha256 = f.read().strip().split()[0]
        actual_sha256 = calculate_sha256(archive_name)
        if actual_sha256 != expected_sha256:
            raise DownloadError("SHA256 checksum does not match expected value. Download may be corrupted. Please try again.")

        print("SHA256 checksum verified successfully.")

        print(f"Extracting GOTTCHA2 database from {archive_name}...")
        with tarfile.open(archive_name) as tar:
            tar.extractall('database')
        print("Database extraction completed.")
    except (DownloadError, OSError, tarfile.TarError, IndexError) as exc:
        remove_if_exists(archive_name)
        remove_if_exists(checksum_name)
        shutil.rmtree('database', ignore_errors=True)
        sys.exit(str(exc))

    remove_if_exists(archive_name)
    remove_if_exists(checksum_name)
    print("Temporary files removed. Database is ready to use.")

def main(args):
    argvs = parse_params(args)
    download_db(argvs)

if __name__ == '__main__':
    main(sys.argv[1:])
