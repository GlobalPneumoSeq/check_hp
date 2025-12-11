#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import re
import natsort
import gffutils
import pyfaidx


def main(gff_path, promoter_length, repeats):
    fasta_path = Path(gff_path).with_suffix(".fasta")
    db_path = f"{gff_path}_db"

    try:
        fasta = extract_fasta(gff_path, fasta_path)
        db = gffutils.create_db(gff_path, db_path, force=True)

        bases = ("A", "T", "C", "G")
        counter = {"promoter": 0, **{f"poly{base}": 0 for base in bases}}

        for contig in natsort.natsorted(db.seqids()):
            contig_fasta = fasta[contig]
            contig_start, contig_end, cds_features = get_contig_details(db, contig, len(contig_fasta))

            for i, feature in enumerate(cds_features):
                upstream_feature, downstream_feature = get_adjacent_features(i, cds_features)
                upstream_end = upstream_feature.end if upstream_feature else contig_start
                downstream_start = downstream_feature.start if downstream_feature else contig_end

                if feature.strand == "+":
                    gap = feature.start - upstream_end
                elif feature.strand == "-":
                    gap = downstream_start - feature.end
                else:
                    raise ValueError("Malformatted GFF file with strand value neither + or -")

                if gap > promoter_length:
                    counter["promoter"] += 1
                    for base in bases:
                        counter[f"poly{base}"] += get_homopolymer_count(feature, contig_fasta, feature.strand, promoter_length, base, repeats)
       
        print(counter)

    finally:
        for file in (db_path, fasta_path, f"{fasta_path}.fai"):
            if os.path.isfile(file):
                os.remove(file)


# Extract FASTA from GFF3
def extract_fasta(gff_path, fasta_path):
    with open(gff_path, "r") as gff_file, open(fasta_path, "w") as fasta_file:
        reached_fasta = False
        for line in gff_file:
            if reached_fasta:
                fasta_file.write(line)
            if line.strip() == "##FASTA":
                reached_fasta = True
        if not reached_fasta:
            raise ValueError("No ##FASTA section found in the provided GFF3 file")
    return pyfaidx.Fasta(fasta_path)


# Extract details regarding the contig
def get_contig_details(db, contig, contig_fasta_length):
    try:
        contig_start = db[contig].start
        contig_end = db[contig].end
    except gffutils.exceptions.FeatureNotFoundError:
        contig_start = 1
        contig_end = contig_fasta_length
    
    cds_features = [feature for feature in db.region(seqid=contig, featuretype="CDS")]
    
    return contig_start, contig_end, cds_features


# Get upstream and downstream features if possible else returns None
def get_adjacent_features(i, features):
    features_count = len(features)
    upstream_feature = features[i - 1] if i > 0 else None
    downstream_feature = features[i + 1] if i + 1 < features_count else None

    return upstream_feature, downstream_feature


def get_homopolymer_count(feature, contig_fasta, strand, promoter_length, base, repeats):
    match_pattern = re.compile(rf"{base}{{{repeats},}}")
    if strand == "+":
        promoter_sequence = contig_fasta[(feature.start - promoter_length - 1) : (feature.start - 1)]
    elif strand == "-":
        promoter_sequence = contig_fasta[feature.end : (feature.end + promoter_length)].reverse.complement
    
    return len(match_pattern.findall(str(promoter_sequence)))


# Parse arguments
def parse_arguments():
    def check_file_exist(file_path):
        if os.path.isfile(file_path):
            return file_path
        else:
            raise argparse.ArgumentTypeError(f'{file_path} does not exist.')

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        'gff',
        type=check_file_exist,
        help='path to a GFF3 file',
    )

    parser.add_argument(
        '-l', '--length',
        default=100,
        type=int,
        help='length of promoter region',
    )

    parser.add_argument(
        '-r', '--repeats',
        default=8,
        type=int,
        help='minimum repeats of the homopolymer',
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    main(args.gff, args.length, args.repeats)
