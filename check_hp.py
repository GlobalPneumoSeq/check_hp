#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
import re
import natsort
import gffutils
import pyfaidx
import tempfile


def main(gff_path, promoter_length, repeats):
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = str(Path(tmpdir) / "tmp.fasta")
        db_path = str(Path(tmpdir) / "tmp.gffutils_db")

        fasta = extract_fasta(gff_path, fasta_path)
        db = gffutils.create_db(gff_path, db_path)

        bases = ("A", "T", "C", "G")
        counter = {"promoter": 0, **{f"poly{base}": 0 for base in bases}}

        for contig in natsort.natsorted(db.seqids()):
            contig_details = get_contig_details(db, contig, fasta, featuretype="CDS")

            for i, feature in enumerate(contig_details["features"]):
                adjacent_features = get_adjacent_features(i, contig_details)
                
                if not has_promoter(feature, adjacent_features, promoter_length):
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


# Get upstream and downstream features if possible else returns None
def get_adjacent_features(cur_feature_index, contig_details):
    features = contig_details["features"]

    features_count = len(features)

    if cur_feature_index > 0:
        upstream_feature = features[cur_feature_index - 1]
        upstream_end = upstream_feature.end
    else:
        upstream_feature = None
        upstream_end = contig_details["start"] - 1

    if cur_feature_index + 1 < features_count: 
        downstream_feature = features[cur_feature_index + 1]
        downstream_start = downstream_feature.start
    else:
        downstream_feature = None
        downstream_start = contig_details["end"] + 1

    return {"upstream_feature": upstream_feature, "upstream_end": upstream_end, "downstream_feature": downstream_feature, "downstream_start": downstream_start}


# Get gap size between current feature and the relevant adjacent feature 
def get_gap(feature, adjacent_features):
    if feature.strand == "+":
        gap = feature.start - adjacent_features["upstream_end"] - 1
    elif feature.strand == "-":
        gap = adjacent_features["downstream_start"] - feature.end - 1
    else:
        raise ValueError("Malformatted GFF file with strand value neither + or -")
    
    return gap


# Check if gap size between current feature and the relevant adjacent feature is larger than promoter_length
def has_promoter(feature, adjacent_features, promoter_length):
    gap = get_gap(feature, adjacent_features)
    return gap >= promoter_length


# Get the homopolymer count of the specific base within the promoter
def get_homopolymer_count(feature, contig_fasta, promoter_length, base, repeats):
    match_pattern = re.compile(rf"{base}{{{repeats},}}")
    if feature.strand == "+":
        promoter_sequence = contig_fasta[(feature.start - promoter_length - 1) : (feature.start - 1)]
    elif feature.strand == "-":
        promoter_sequence = contig_fasta[feature.end : (feature.end + promoter_length)].reverse.complement
    else:
        raise ValueError("Malformatted GFF file with strand value neither + or -")
    
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
