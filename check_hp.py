#!/usr/bin/env python3

import argparse
import os
import re
import tempfile
import pandas as pd
import gffutils
import pyfaidx
from intervaltree import IntervalTree
import natsort


def main(gff_path, promoter_length, repeats):
    # Nucleotide bases to search for
    bases = ("A", "T", "C", "G")

    # To store promoters summary and homopolyers details
    promoters_summary_data = []
    homopolymers_details_data = []

    # Work within a temporary directory to handle extracted FASTA file and associated faidx index
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as workdir:
        all_contigs_details = get_all_contigs_details(gff_path, workdir, featuretype="CDS")

        for contig in all_contigs_details:
            contig_details = all_contigs_details[contig]

            features_intervaltree = get_features_intervaltree(contig_details["features"])

            for feature in contig_details["features"]:
                validate_feature_strand(feature, contig)

                if promoter_feature_overlapped(feature, features_intervaltree, promoter_length):
                    continue

                if not promoter_within_contig(feature, contig_details, promoter_length):
                    continue

                # Shared row data template for both promoters summary and homopolymers details
                row_data_template = {
                    "contig": contig,
                    "strand": feature.strand,
                    "cds_id": feature.id,
                    "cds_start_pos": feature.start,
                    "cds_end_pos": feature.end,
                    "gene": feature.attributes.get("gene", [""])[0],
                    "protein": feature.attributes.get("Name", [""])[0]
                }

                # To collect the homopolymer counts for each base in the promoter region
                promoters_summary_row_data = {}

                for base in bases:
                    homopolymers = get_homopolymers(feature, contig_details["fasta"], promoter_length, base, repeats)

                    promoters_summary_row_data[f"poly{base}"] = len(homopolymers)

                    for homopolymer in homopolymers:
                        homopolymers_details_data.append({**row_data_template, "repeated_base": base, **homopolymer})

                promoters_summary_data.append({**row_data_template, **promoters_summary_row_data})

    save_results(promoters_summary_data, homopolymers_details_data)


# Save details of all contigs from the GFF3 file to a nested dictionary
def get_all_contigs_details(gff_path, workdir, featuretype):
    contigs_details = {}

    db = gffutils.create_db(data=gff_path, dbfn=":memory:")

    with tempfile.NamedTemporaryFile(mode="w+t", dir=workdir, delete=False) as tmp_fasta_file:
        fasta = extract_fasta(gff_path, tmp_fasta_file.name)
        tmp_fasta_file.close()

        for contig in db.seqids():
            contigs_details[contig] = get_contig_details(db, contig, fasta, featuretype)

    return contigs_details


# Extract FASTA from the GFF3 file 
def extract_fasta(gff_path, fasta_path):
    with open(gff_path, "r") as gff_file, open(fasta_path, "w") as fasta_file:
        reached_fasta = False
       
        for line in gff_file:
            if line.strip() == "##FASTA":
                if reached_fasta:
                    raise ValueError("Multiple ##FASTA sections found in the provided GFF3 file")
                reached_fasta = True
                continue
            
            if reached_fasta:
                fasta_file.write(line)

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


# Create an interval tree for the features
def get_features_intervaltree(features):
    intervals = []

    for feature in features:
        # interval end is feature.end + 1 as IntervalTree is non-inclusive of the upper limit
        intervals.append((feature.start, feature.end + 1))
    
    return IntervalTree.from_tuples(intervals)


# Validate strand value of a feature
def validate_feature_strand(feature, contig):
    if feature.strand in ("+", "-"):
        return
    
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


# Get details of homopolymers of a specific base within the promoter of a feature
def get_homopolymers(feature, contig_fasta, promoter_length, base, repeats):
    match_pattern = re.compile(rf"{base}{{{repeats},}}")

    if feature.strand == "+":
        promoter_start_index = feature.start - promoter_length - 1
        promoter_end_index = feature.start - 1

        promoter_sequence = contig_fasta[promoter_start_index:promoter_end_index]

        matches = [
            {
                "repeats": len(match.group()),
                "repeat_start_pos": match.start() + promoter_start_index + 1,
                "repeat_end_pos": (repeat_end_pos := match.end() + promoter_start_index),
                "distance_from_cds": feature.start - repeat_end_pos
            }
            for match in match_pattern.finditer(str(promoter_sequence).upper())
        ]

    elif feature.strand == "-":
        promoter_start_index = feature.end
        promoter_end_index = feature.end + promoter_length

        promoter_sequence = contig_fasta[promoter_start_index:promoter_end_index].reverse.complement

        matches = [
            {
                "repeats": len(match.group()),
                "repeat_start_pos": (repeat_start_pos := promoter_start_index + promoter_length - match.end() + 1),
                "repeat_end_pos": promoter_start_index + promoter_length - match.start(),
                "distance_from_cds": repeat_start_pos - feature.end
            }
            for match in match_pattern.finditer(str(promoter_sequence).upper())
        ]

    return matches


# Save results to CSV files
# Rows sorted by contig then cds_id using natural sorting
def save_results(promoters_summary_data, homopolymers_details_data):
    results = [
        {   
            "data": promoters_summary_data,
            "filename": "promoters_summary.csv",  
            "columns":  ["contig", "strand", "cds_id", "cds_start_pos", "cds_end_pos", "gene", "protein", "polyA", "polyT", "polyC", "polyG" ]
        },
        {
            "data": homopolymers_details_data,
            "filename": "homopolymers_details.csv",
            "columns": ["contig", "strand", "cds_id", "cds_start_pos", "cds_end_pos", "gene", "protein", "repeated_base", "repeats", "repeat_start_pos", "repeat_end_pos", "distance_from_cds"]
        }
    ]

    for result in results:
        df = pd.DataFrame(result["data"], columns=result["columns"]).sort_values(by=["contig", "cds_id"], key=natsort.natsort_keygen())
        df.to_csv(result["filename"], index=False)


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
