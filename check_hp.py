#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import re
import natsort
import gffutils
import pyfaidx
import tempfile
from intervaltree import IntervalTree


def main(gff_path, promoter_length, repeats):
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = str(Path(tmpdir) / "tmp.fasta")
        fasta = extract_fasta(gff_path, fasta_path)
        db = gffutils.create_db(data=gff_path, dbfn=":memory:")

        bases = ("A", "T", "C", "G")
        counter = {"promoter": 0, **{f"poly{base}": 0 for base in bases}}

        for contig in natsort.natsorted(db.seqids()):
            contig_details = get_contig_details(db, contig, fasta, featuretype="CDS")

            # Create an interval tree for the features in the contig
            # feature.end + 1 as IntervalTree is non-inclusive of the upper limit
            features_intervals = IntervalTree.from_tuples([(feature.start, feature.end + 1) for feature in contig_details["features"]])

            for feature in contig_details["features"]:
                validate_feature_strand(feature, contig)

                if promoter_feature_overlapped(feature, features_intervals, promoter_length):
                    continue

                if not promoter_within_contig(feature, contig_details, promoter_length):
                    continue

                counter["promoter"] += 1
                for base in bases:
                    counter[f"poly{base}"] += get_homopolymer_count(feature, contig_details["fasta"], promoter_length, base, repeats)
       
        print(counter)


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


# Extract contig details
def get_contig_details(db, contig, fasta, featuretype):
    contig_fasta = fasta[contig]

    try:
        contig_start = db[contig].start
        contig_end = db[contig].end
    except gffutils.exceptions.FeatureNotFoundError:
        contig_start = 1
        contig_end = len(contig_fasta)
    
    features = list(db.region(seqid=contig, featuretype=featuretype))
    
    return {"fasta": contig_fasta, "start": contig_start, "end": contig_end, "features": features}


# Validate value of strand of a feature
def validate_feature_strand(feature, contig):
    if feature.strand not in ("+", "-"):
        raise ValueError(f"Malformed GFF3: expected strand value of + or -, got {feature.strand} for feature {feature.id} on contig {contig}")


# Check if the promoter region overlaps with any features
def promoter_feature_overlapped(feature, features_intervals, promoter_length):
    if feature.strand == "+":
        promoter_start = feature.start - promoter_length
        promoter_end = feature.start - 1
    elif feature.strand == "-":
        promoter_start = feature.end + 1
        promoter_end = feature.end + promoter_length
    # promoter_end + 1 as IntervalTree is non-inclusive of the upper limit
    return len(features_intervals[promoter_start:promoter_end + 1]) > 0


# Check if the promoter region within the contig boundaries
def promoter_within_contig(feature, contig_details, promoter_length):
    if feature.strand == "+":
        return feature.start - promoter_length >= contig_details["start"]
    elif feature.strand == "-":
        return feature.end + promoter_length <= contig_details["end"]


# Get the homopolymer count of the specific base within the promoter
def get_homopolymer_count(feature, contig_fasta, promoter_length, base, repeats):
    match_pattern = re.compile(rf"{base}{{{repeats},}}")
    if feature.strand == "+":
        promoter_sequence = contig_fasta[(feature.start - promoter_length - 1) : (feature.start - 1)]
    elif feature.strand == "-":
        promoter_sequence = contig_fasta[feature.end : (feature.end + promoter_length)].reverse.complement
    return len(match_pattern.findall(str(promoter_sequence).upper()))


# Parse arguments
def parse_arguments():
    def check_valid_gff3(gff_path):
        if not os.path.isfile(gff_path):
            raise argparse.ArgumentTypeError(f'{gff_path} does not exist.')
        
        with open(gff_path, "r") as gff_file:
            if gff_file.readline().strip() != "##gff-version 3":
                raise argparse.ArgumentTypeError(f'{gff_path} is not a valid GFF3 file.')
            
        return gff_path

    def check_int_more_than_one(val):
        try:
            val_int = int(val)
            if val_int <= 1:
                raise argparse.ArgumentTypeError(f'{val} is not a valid integer greater than 1.')
            return val_int
        except ValueError:
            raise argparse.ArgumentTypeError(f'{val} is not a valid integer.')

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        'gff',
        type=check_valid_gff3,
        help='path to a GFF3 file',
    )

    parser.add_argument(
        '-l', '--length',
        default=100,
        type=check_int_more_than_one,
        help='length of promoter region',
    )

    parser.add_argument(
        '-r', '--repeats',
        default=8,
        type=check_int_more_than_one,
        help='minimum repeats of the homopolymer',
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    main(args.gff, args.length, args.repeats)
