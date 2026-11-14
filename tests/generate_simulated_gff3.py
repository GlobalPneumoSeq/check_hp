# Generate deterministic Bakta-style GFF3 fixtures

from pathlib import Path
import pandas as pd
import random
import string
import textwrap


def main(fixtures_count=10):
    fixtures_path = get_fixtures_path()

    all_homopolymers = []

    for fixture_index in range(1, fixtures_count + 1):
        fixture_homopolymers = make_fixture(fixture_index, fixtures_path)
        all_homopolymers.extend(fixture_homopolymers)

    save_expected(fixtures_path, all_homopolymers)


# Generate and save a GFF3 fixture file, return its expected homopolymer records
def make_fixture(fixture_index, fixtures_path, contig_count_range=(2, 4), cds_count_range=(4, 8), contig_length_range=(5000, 7000), cds_length_range=(100, 400), promoter_length=100, homopolymer_length_range=(8, 15)):
    fixture_filename = f"SIM{fixture_index:06d}.gff3"
    bases = "ATCG"

    # Ensure deterministic output by seeding with the fixture index
    random_generator = random.Random(fixture_index)

    content_features = ["##gff-version 3", "# organism Simulated bacterium", "# Annotated with Bakta"]
    content_sequences = ['##FASTA']
    homopolymers = []
    cds_overall_counter = 0
    locus_tag_prefix = "".join(random_generator.choices(string.ascii_uppercase, k=6))

    contig_count = random_generator.randint(contig_count_range[0], contig_count_range[1])
    for contig_index in range(contig_count):
        contig, contig_content = make_contig(contig_index, bases, contig_length_range, random_generator)
        content_features.extend(contig_content)

        cds_count = random_generator.randint(cds_count_range[0], cds_count_range[1])
        for cds_index in range(cds_count):
            cds_overall_counter += 1
            cds, cds_content = make_cds(cds_index, cds_overall_counter, locus_tag_prefix, contig, cds_count, cds_length_range, promoter_length, random_generator)
            content_features.extend(cds_content)

            homopolymer_count = random_generator.randint(0, 3)
            for homopolymer_index in range(homopolymer_count):
                homopolymer = make_homopolymer(homopolymer_index, contig, cds, homopolymer_count, bases, homopolymer_length_range, promoter_length, random_generator)
                homopolymers.append({"file": fixture_filename, **homopolymer})

        # Save sequence only after changes made by make_homopolymer()
        content_sequences.extend([f">{contig["name"]}", *textwrap.wrap("".join(contig["sequence"]), width=60)])

    save_fixture(fixtures_path / fixture_filename, content_features, content_sequences)

    return homopolymers


# Generate contig, return contig details and content to add to fixture file
def make_contig(contig_index, bases, contig_length_range, random_generator):
    name = f"contig_{contig_index + 1}"
    base_count = len(bases)

    length = random_generator.randint(contig_length_range[0], contig_length_range[1])
    sequence = list(bases * (length // base_count) + bases[:(length % base_count)])

    contig_content = [
        f"##sequence-region {name} 1 {length}",
        "\t".join([
            name, "Bakta", "region", "1", str(length), ".", "+", ".", f"ID={name};Name={name}"
        ])
    ]

    contig_details = {
        "name": name,
        "length": length,
        "sequence": sequence
    }

    return contig_details, contig_content


# Generate CDS, return CDS details and content to add to fixture file
def make_cds(index, overall_counter, locus_tag_prefix, contig, count, length_range, promoter_length, random_generator):
    cds_id = f"{locus_tag_prefix}_{overall_counter:05d}"
    strand = random_generator.choice(["+", "-"])
    gene = f"sim_gene_{overall_counter}"
    protein = f"Simulated protein {overall_counter}"
    length = random_generator.randint(length_range[0], length_range[1])

    # Assign each CDS to a different fixed size segment in the contig
    segment_size = contig["length"] // count
    # +1 as GFF3 use 1-based index
    possible_start_range_start = segment_size * index + 1  
    possible_start_range_end = segment_size * (index + 1) - length + 1

    if strand == "+":
        start = random_generator.randint(possible_start_range_start + promoter_length, possible_start_range_end)
    elif strand == "-":
        start = random_generator.randint(possible_start_range_start, possible_start_range_end - promoter_length)
    end = start + length - 1

    cds_content = [
        "\t".join([
            contig["name"], "Pyrodigal", "CDS", str(start), str(end), ".", strand, "0",
            f"ID={cds_id};Name={protein};gene={gene};locus_tag={cds_id}"
        ])
    ]

    cds_details = {
        "strand": strand,
        "id": cds_id,
        "start": start,
        "end": end,
        "gene": gene,
        "protein": protein
    }

    return cds_details, cds_content

# Generate homopolymer and modify relevant sequence, return homopolymer details
def make_homopolymer(index, contig, cds, count, bases, length_range, promoter_length, random_generator):
    repeated_base = random_generator.choice(bases)

    length = random_generator.randint(length_range[0], length_range[1])

    # Assign each homopolymer to a different fixed size segment in the promoter
    segment_size = promoter_length // count
    possible_distance_range_start = segment_size * index
    possible_distance_range_end = segment_size * (index + 1) - length - 2 # -2 to ensure neighboring homopolymers are not using the same position for their flanking bases

    distance = random_generator.randint(possible_distance_range_start, possible_distance_range_end)

    if cds["strand"] == "+":
        ref_seq_base = repeated_base
        homopolymer_end = cds["start"] - distance - 1
        homopolymer_start = homopolymer_end - length + 1
    elif cds["strand"] == "-":
        ref_seq_base = repeated_base.translate(str.maketrans("ACGTU", "TGCAA"))
        homopolymer_start = cds["end"] + distance + 1
        homopolymer_end = homopolymer_start + length - 1

    # -1 as GFF3 use 1-based index 
    start_index, end_index = homopolymer_start - 1, homopolymer_end - 1

    contig["sequence"][start_index:end_index + 1] = [ref_seq_base] * length

    flanking_base = random_generator.choice([base for base in bases if base != ref_seq_base])
    # Safeguard from the case of contig["sequence"][-1]
    if start_index > 0:
        contig["sequence"][start_index - 1] = flanking_base
    # Safeguard from attempting to access non-existent index
    if end_index < len(contig["sequence"]) - 1:
        contig["sequence"][end_index + 1] = flanking_base

    homopolymer_details = {
        "contig": contig["name"],
        "strand": cds["strand"],
        "cds_id": cds["id"],
        "cds_start_pos": cds["start"],
        "cds_end_pos": cds["end"],
        "gene": cds["gene"],
        "protein": cds["protein"],
        "repeated_base": repeated_base,
        "repeats": length,
        "repeat_start_pos": homopolymer_start,
        "repeat_end_pos": homopolymer_end,
        "distance_from_cds": distance,
    }

    return homopolymer_details


def get_fixtures_path():
    fixtures_path = Path(__file__).parent / "fixtures"
    fixtures_path.mkdir(parents=True, exist_ok=True)
    return fixtures_path


def save_fixture(fixture_path, content_features, content_sequences):
    (fixture_path).write_text("\n".join(content_features + content_sequences))


def save_expected(fixtures_path, all_homopolymers, file_name="expected_homopolymers.csv"):
    pd.DataFrame(all_homopolymers).to_csv(fixtures_path / file_name, index=False)


if __name__ == "__main__":
    main()
