# Transcript Effects

### VCF INFO header lines

```
##SnpEffVersion="4.2 (build 2015-12-05)"
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele |
Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID |
Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length
| AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
```

```
##VEP="v96" time="2019-05-14 17:42:50" cache="/opt/vep/.vep/homo_sapiens/96_GRCh38"
ensembl=96.7a35428 ensembl-funcgen=96.9c3a0cd ensembl-variation=96.70d2777 ensembl-io=96.6e65b30
1000genomes="phase3" COSMIC="87" ClinVar="201901" ESP="V2-SSA137" HGMD-PUBLIC="20184"
assembly="GRCh38.p12" dbSNP="151" gencode="GENCODE 30" genebuild="2014-07" gnomAD="r2.1"
polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"
##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP.
Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|
HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|
DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
```

```
##VEP="v101" time="2020-10-28 14:12:58" cache="/opt/vep/.vep/homo_sapiens/101_GRCh38"
ensembl-funcgen=101.b918a49 ensembl-io=101.943b6c2 ensembl=101.856c8e8 ensembl-variation=101.50e7372
1000genomes="phase3" COSMIC="90" ClinVar="202003" ESP="V2-SSA137" HGMD-PUBLIC="20194"
assembly="GRCh38.p13" dbSNP="153" gencode="GENCODE 35" genebuild="2014-07" gnomAD="r2.1"
polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"
##INFO=<ID=ANN,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP.
Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|
HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|
DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID">
```

### Summary

| Field | bdg-formats                       | snpEff 4.2               | Ensembl VEP v96,v101  | Notes                                                                                                              |
| ----: | :-------------------------------- | :----------------------- | :-------------------- | :----                                                                                                              |
| 0     | `alternateAllele`                 | Allele                   | Allele                |                                                                                                                    |
| 1     | `effects`                         | Annotation               | Consequence           |                                                                                                                    |
| 2     | ``                                | Annotation_Impact        | IMPACT                | Putative_impact: A simple estimation of putative impact / deleteriousness : {HIGH, MODERATE, LOW, MODIFIER}        |
| 3     | `geneName`                        | Gene_Name                | SYMBOL                |                                                                                                                    |
| 4     | `geneId`                          | Gene_ID                  | Gene                  |                                                                                                                    |
| 5     | `featureType`                     | Feature_Type             | Feature_type          |                                                                                                                    |
| 6     | `featureId`                       | Feature_ID               | Feature               |                                                                                                                    |
| 7     | `biotype`                         | Transcript_BioType       | BIOTYPE               |                                                                                                                    |
| 8     | `rank`,`total`                    | Rank                     | EXON                  | Rank / total : Exon or Intron rank / total number of exons or introns.                                             |
| 9     | `transcriptHgvs`                  | HGVS.c                   | INTRON                | Rank / total : Exon or Intron rank / total number of exons or introns.                                             |
| 10    | `proteinHgvs`                     | HGVS.p                   | HGVSc                 | Map to `transcriptHgvs`                                                                                            |
| 11    | `cdnaPosition`,`cdnaLength`       | cDNA.pos / cDNA.length   | cDNA_position         | cDNA_position ​/ (cDNA_len optional)_ : Position in cDNA and trancript's cDNA length (one based).                   |
| 12    | `cdsPosition`,`cdsLength`         | CDS.pos / CDS.length     | CDS_position          | CDS_position ​/ (CDS_len optional)​: Position and number of coding bases (one based includes START and STOP codons). |
| 13    | `proteinPosition`,`proteinLength` | AA.pos / AA.length       | Protein_position      | Protein_position ​/ (Protein_len optional)​: Position and number of AA (one based, including START, but not STOP).   |
| 14    | `distance`                        | Distance                 | Amino_acids           | Amino_acids : Reference and variant amino acids                                                                    |
| 15    | `messages`                        | ERRORS / WARNINGS / INFO | Codons                | Codons : Reference and variant codon sequence                                                                      |
| 16    | ``                                |                          | Existing_variation    | Existing_variation : Identifier(s) of co-located known variants                                                    |
| 17    | ``                                |                          | DISTANCE              | Map to `distance`                                                                                                  |
| 18    | ``                                |                          | STRAND                | STRAND : Strand of the feature (1/-1)                                                                              |
| 19    | ``                                |                          | FLAGS                 | FLAGS : Transcript quality flags                                                                                   |
| 20    | ``                                |                          | SYMBOL_SOURCE         |                                                                                                                    |
| 21    | ``                                |                          | HGNC_ID               |                                                                                                                    |
