import pyfaidx
import pytest

from check_hp import extract_fasta


# Tests

def test_extract_fasta_typical(tmp_path):
    gff_path = tmp_path / "test.gff"
    fasta_path = tmp_path / "test.fasta"

    gff_path.write_text(
        "##gff-version 3\n"
        "##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo\n"
        "# organism Streptococcus pneumoniae\n"
        "# Annotated with Bakta\n"
        "# Software: v1.11.4\n"
        "# Database: v6.0, full\n"
        "# DOI: 10.1099/mgen.0.000685\n"
        "# URL: github.com/oschwengers/bakta\n"
        "#sequence-region contig_1 1 1000\n"
        "contig_1	Bakta	region	1	5	.	+	.	ID=contig_1;Name=contig_1\n"
        "contig_2	Bakta	region	1	6	.	+	.	ID=contig_2;Name=contig_2\n"
        "##FASTA\n"
        ">contig_1\n"
        "ATCGA\n"
        ">contig_2\n"
        "TTCGCA\n"
    )

    fasta = extract_fasta(gff_path, fasta_path)

    assert isinstance(fasta, pyfaidx.Fasta)
    assert set(fasta.keys()) == {"contig_1", "contig_2"}
    assert str(fasta["contig_1"]) == "ATCGA"
    assert str(fasta["contig_2"]) == "TTCGCA"


def test_extract_fasta_no_fasta(tmp_path):
    gff_path = tmp_path / "test.gff"
    fasta_path = tmp_path / "test.fasta"

    gff_path.write_text(
        "##gff-version 3\n"
        "##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo\n"
        "# organism Streptococcus pneumoniae\n"
        "# Annotated with Bakta\n"
        "# Software: v1.11.4\n"
        "# Database: v6.0, full\n"
        "# DOI: 10.1099/mgen.0.000685\n"
        "# URL: github.com/oschwengers/bakta\n"
        "#sequence-region contig_1 1 1000\n"
        "contig_1	Bakta	region	1	5	.	+	.	ID=contig_1;Name=contig_1\n"
        "contig_2	Bakta	region	1	6	.	+	.	ID=contig_2;Name=contig_2\n"
    )

    with pytest.raises(ValueError, match=r"No ##FASTA section"):
        extract_fasta(gff_path, fasta_path)


def test_extract_fasta_empty_fasta(tmp_path):
    gff_path = tmp_path / "test.gff"
    fasta_path = tmp_path / "test.fasta"

    gff_path.write_text(
        "##gff-version 3\n"
        "##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/so.obo\n"
        "# organism Streptococcus pneumoniae\n"
        "# Annotated with Bakta\n"
        "# Software: v1.11.4\n"
        "# Database: v6.0, full\n"
        "# DOI: 10.1099/mgen.0.000685\n"
        "# URL: github.com/oschwengers/bakta\n"
        "#sequence-region contig_1 1 1000\n"
        "contig_1	Bakta	region	1	5	.	+	.	ID=contig_1;Name=contig_1\n"
        "contig_2	Bakta	region	1	6	.	+	.	ID=contig_2;Name=contig_2\n"
        "##FASTA\n"
    )

    with pytest.raises(pyfaidx.FastaIndexingError, match=r"The FASTA file .+ does not contain a valid sequence. Check that sequence definition lines start with '>'."):
        extract_fasta(gff_path, fasta_path)
