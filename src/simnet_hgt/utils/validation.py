from pathlib import Path

import pandas as pd
from Bio import SeqIO


def iter_fasta_paths(folder: Path):
    for p in sorted(folder.glob("*")):
        if p.is_file() and str(p.name).lower().endswith(".faa"):
            yield p


def count_sequences_in_file(p: Path) -> int:
    try:
        with open(p) as handle:
            # allow SeqIO to detect format via parser 'fasta'
            return sum(1 for _ in SeqIO.parse(handle, "fasta"))
    except Exception as e:
        # return -1 to flag failure; caller can inspect
        print(f"WARNING: failed to parse {p}: {e}")
        return -1


def summarize_folder(folder: Path) -> tuple[int, int, list]:
    """Return (num_files, total_sequences, list_of_file_records)"""
    file_records = []
    total = 0
    nfiles = 0
    for p in iter_fasta_paths(folder):
        nfiles += 1
        nseq = count_sequences_in_file(p)
        file_records.append({"folder": str(folder), "file": str(p), "sequences": nseq})
        if nseq > 0:
            total += nseq
    return nfiles, total, file_records


def summarize_protein_data(root: Path) -> None:
    folders = [root / "archaea", root / "bacteria", root / "plasmid", root / "virus"]
    labels = ["archaea", "bacteria", "plasmids", "viruses"]

    all_file_records = []
    summary_rows = []
    for lab, folder in zip(labels, folders, strict=False):
        if not folder.exists():
            print(f"WARNING: folder {folder} does not exist, skipping.")
            summary_rows.append({"category": lab, "n_files": 0, "total_proteins": 0})
            continue
        nfiles, total, file_records = summarize_folder(folder)
        summary_rows.append({"category": lab, "n_files": nfiles, "total_proteins": total})
        # add category field to file records
        for rec in file_records:
            rec["category"] = lab
        all_file_records.extend(file_records)

    df_summary = pd.DataFrame(summary_rows)
    print("Summary of protein counts by category:")
    print(df_summary.to_string(index=False))
