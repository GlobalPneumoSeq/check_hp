from check_hp import get_all_contigs_details


# Tests

def test_get_all_contigs_details_loads_gff_features_and_fasta(tmp_path):
    gff_path = tmp_path / "annotation.gff3"
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
        "contig_1	Bakta	CDS	2	3	.	+	.	ID=contig_1_feature_1_id;Name=contig_1_feature_1_name\n"
        "contig_2	Bakta	region	1	4	.	+	.	ID=contig_2;Name=contig_2\n"
        "##FASTA\n"
        ">contig_1\n"
        "ATCGA\n"
        ">contig_2\n"
        "GCAT\n"
    )

    contigs = get_all_contigs_details(gff_path, tmp_path, featuretype="CDS")

    assert list(contigs) == ["contig_1", "contig_2"]
    assert contigs["contig_1"]["start"] == 1
    assert contigs["contig_1"]["end"] == 5
    assert contigs["contig_2"]["start"] == 1
    assert contigs["contig_2"]["end"] == 4
    assert str(contigs["contig_1"]["fasta"]) == "ATCGA"
    assert str(contigs["contig_2"]["fasta"]) == "GCAT"
    assert [feature.id for feature in contigs["contig_1"]["features"]] == ["contig_1_feature_1_id"]
    assert contigs["contig_2"]["features"] == []
