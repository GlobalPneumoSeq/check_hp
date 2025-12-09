#!/usr/bin/env python3

import argparse
import gffutils
import os
import natsort


def main(gff_filename, promoter_region_length):
    try:
        fasta_filename = f"{'.'.join(gff_filename.split(".")[:-1])}.fasta"
        db_filename = f"{gff_filename}_db"

        extract_fasta(gff_filename, fasta_filename)
        db = gffutils.create_db(gff_filename, db_filename, force=True)

        promoter_region_count = 0

        for contig in natsort.natsorted(db.seqids()):
            contig_start, contig_end, cds_features = get_contig_details(db, contig)

            for i, feature in enumerate(cds_features):
                upstream_feature, downstream_feature = get_adjacent_features(i, cds_features)

                upstream_end = upstream_feature.end if upstream_feature else contig_start
                downstream_start = downstream_feature.start if downstream_feature else contig_end

                if feature.strand == "+":
                    if feature.start - upstream_end > promoter_region_length:
                        promoter_region_count += 1
                        # TODO: POLYCOUNT
                elif feature.strand == "-":
                    if downstream_start - feature.end > promoter_region_length:
                        promoter_region_count += 1
                        # TODO: POLYCOUNT
                else:
                    raise ValueError("Malformatted GFF file with strand value neither + or -")

        print(promoter_region_count)

    finally:
        for file in (db_filename, fasta_filename):
            if os.path.isfile(file):
                os.remove(file)


# Extract FASTA from GFF3
def extract_fasta(gff_filename, fasta_filename):
    with open(gff_filename, "r") as gff_file, open(fasta_filename, "w") as fasta_file:
        reached_fasta = False
        for line in gff_file:
            if reached_fasta:
                fasta_file.write(line)
            if line.strip() == "##FASTA":
                reached_fasta = True


# Extract details regarding the contig
def get_contig_details(db, contig):
    contig_start = db[contig].start
    contig_end = db[contig].end
    cds_features = [feature for feature in db.region(seqid=contig, featuretype="CDS")]
    
    return contig_start, contig_end, cds_features


# Get upstream and downstream features if possible else returns None
def get_adjacent_features(i, features):
    features_count = len(features)
    upstream_feature = features[i - 1] if i > 0 else None
    downstream_feature = features[i + 1] if i + 1 < features_count else None

    return upstream_feature, downstream_feature


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

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    main(args.gff, args.length)
