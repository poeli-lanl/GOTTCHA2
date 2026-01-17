#!/usr/bin/env python3
import argparse
import gzip
import sys
from typing import Optional

def open_in(path: str):
    if path == "-":
        return sys.stdin.buffer
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb", buffering=5 * 1024 * 1024)

def open_out(path: str, gz_level: int):
    if path == "-":
        return sys.stdout.buffer
    if path.endswith(".gz"):
        # compresslevel=1 is much faster than default(9) for big outputs
        return gzip.open(path, "wb", compresslevel=gz_level)
    return open(path, "wb", buffering=5 * 1024 * 1024)

def detect_format(first_line: bytes) -> str:
    if not first_line:
        raise ValueError("Empty input")
    b0 = first_line[:1]
    if b0 == b">":
        return "fasta"
    if b0 == b"@":
        return "fastq"
    raise ValueError("Unknown format: expected FASTA (>) or FASTQ (@)")

def write_chunks_fasta(out, rid: bytes, seq: bytes, L: int, step: int, drop_tail: bool, min_tail: int, prefix: bytes):
    n = len(seq)
    if n == 0:
        return

    out_lines = []
    chunk_id = 0
    chunks_written = 0
    # full chunks
    last_full_start = n - L
    if last_full_start >= 0:
        for s in range(0, last_full_start + 1, step):
            e = s + L
            # >id|chunk=k
            header = b">" + prefix + rid + b"|chunk=" + str(chunk_id).encode() + b"\n"
            out_lines.append(header)
            out_lines.append(seq[s:e] + b"\n")
            chunk_id += 1
            chunks_written += 1

        # tail
        tail_start = ((last_full_start // step) + 1) * step
        if not drop_tail and tail_start < n:
            tail_len = n - tail_start
            if tail_len >= min_tail:
                header = b">" + prefix + rid + b"|chunk=" + str(chunk_id).encode() + b"\n"
                out_lines.append(header)
                out_lines.append(seq[tail_start:] + b"\n")
                chunk_id += 1
                chunks_written += 1
    else:
        # read shorter than L
        if not drop_tail and n >= min_tail:
            header = b">" + prefix + rid + b"|chunk=" + str(chunk_id).encode() + b"\n"
            out_lines.append(header)
            out_lines.append(seq + b"\n")
            chunk_id += 1
            chunks_written += 1

    out.writelines(out_lines)
    return chunks_written

def split_to_fasta(input_path: str, 
                   output_path: str, 
                   split_length: int, 
                   step_length: int,
                   drop_tail: bool, 
                   min_tail: int = 1,
                   prefix: str = "", 
                   gzip_level: int = 1) -> int:
    """
    Split long reads from input_path into fixed-length chunks and write to output_path.

    Returns the total number of chunks written.
    """
    prefix_b = prefix.encode()
    chunk_count = 0

    with open_in(input_path) as hin, open_out(output_path, gzip_level) as hout:
        first = hin.readline()
        if not first:
            return 0
        fmt = detect_format(first)

        if fmt == "fastq":
            # FASTQ assumed standard 4-line records (seq on one line).
            header = first
            while header:
                if not header.startswith(b"@"):
                    raise ValueError("Malformed FASTQ: expected '@' header")
                seq = hin.readline()
                plus = hin.readline()
                qual = hin.readline()
                if not seq or not plus or not qual:
                    break
                # read id up to first whitespace
                rid = header[1:].strip().split(None, 1)[0]
                seq = seq.strip()
                chunk_count += write_chunks_fasta(hout, rid, seq, split_length, step_length, drop_tail, min_tail, prefix_b)
                header = hin.readline()

        else:
            # FASTA: collect per record (still fast; avoids per-base shifting).
            header = first
            rid = None
            seq_parts = []
            while True:
                if header.startswith(b">"):
                    if rid is not None:
                        seq = b"".join(seq_parts)
                        chunk_count += write_chunks_fasta(hout, rid, seq, split_length, step_length, drop_tail, min_tail, prefix_b)
                    rid = header[1:].strip().split(None, 1)[0]
                    seq_parts = []
                else:
                    seq_parts.append(header.strip())

                header = hin.readline()
                if not header:
                    break

            if rid is not None:
                seq = b"".join(seq_parts)
                chunk_count += write_chunks_fasta(hout, rid, seq, split_length, step_length, drop_tail, min_tail, prefix_b)

    return chunk_count

def main():
    ap = argparse.ArgumentParser(description="Split long reads into fixed-length chunks; output FASTA only (input FASTA/FASTQ, .gz ok).")
    ap.add_argument("-i", "--input", required=True, help="Input FASTA/FASTQ (.gz ok) or '-'")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA (.gz ok) or '-'")
    ap.add_argument("-l", "--length", type=int, default=150, help="Chunk length (bp), default 150")
    ap.add_argument("--step", type=int, help="Step between chunk starts (default = length, no overlap)")
    ap.add_argument("--drop-tail", action="store_true", help="Drop final shorter tail chunk")
    ap.add_argument("--min-tail", type=int, default=1, help="If keeping tail, minimum tail length to emit (default 1)")
    ap.add_argument("--prefix", default="", help="Prefix added to read IDs (optional)")
    ap.add_argument("--gzip-level", type=int, default=1, help="Gzip compression level for .gz output (1 fastest, 9 smallest). Default 1")
    args = ap.parse_args()

    L = args.length
    step_length = args.step if args.step is not None else L
    if L <= 0 or step_length <= 0:
        raise ValueError("--length and --step must be positive")
    if not (1 <= args.gzip_level <= 9):
        raise ValueError("--gzip-level must be 1..9")

    split_to_fasta(
        input_path=args.input,
        output_path=args.output,
        split_length=L,
        step_length=step_length,
        drop_tail=args.drop_tail,
        min_tail=args.min_tail,
        prefix=args.prefix,
        gzip_level=args.gzip_level,
    )

if __name__ == "__main__":
    main()
