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

| bdg-formats                                      | Index | snpEff 4.2               | Index | Ensembl VEP v96,v101  | Notes                                                                                                                 |
| :--------------------------------                | ----: | :----------------------- | ----: | :-------------------- | :-------------------------------------------------------------------------------------------------------------------- |
| `alternateAllele`                                | 0     | Allele                   | 0     | Allele                |                                                                                                                       |
| `effects`                                        | 1     | Annotation               | 1     | Consequence           |                                                                                                                       |
| `impact`                                         | 2     | Annotation_Impact        | 2     | IMPACT                | New field `impact`                                                                                                    |
| `geneName`                                       | 3     | Gene_Name                | 3     | SYMBOL                |                                                                                                                       |
| `geneId`                                         | 4     | Gene_ID                  | 4     | Gene                  |                                                                                                                       |
| `featureType`                                    | 5     | Feature_Type             | 5     | Feature_type          |                                                                                                                       |
| `featureId`                                      | 6     | Feature_ID               | 6     | Feature               |                                                                                                                       |
| `biotype`                                        | 7     | Transcript_BioType       | 7     | BIOTYPE               |                                                                                                                       |
| `rank`, `total`                                  | 8     | Rank                     | 8     | EXON                  | Rank / total : Exon or Intron rank / total number of exons or introns; map to `rank`, `total`                         |
|                                                  |       |                          | 9     | INTRON                | Map to `rank`, `total`                                                                                                |
| `transcriptHgvs`                                 | 9     | HGVS.c                   | 10    | HGVSc                 |                                                                                                                       |
| `proteinHgvs`                                    | 10    | HGVS.p                   | 11    | HGVSp                 |                                                                                                                       |
| `cdnaPosition`, `cdnaLength`                     | 11    | cDNA.pos / cDNA.length   | 12    | cDNA_position         | cDNA_position _/ (cDNA_len optional)_ : Position in cDNA and trancript's cDNA length (one based).                     |
| `codingSequencePosition`, `codingSequenceLength` | 12    | CDS.pos / CDS.length     | 13    | CDS_position          | CDS_position _/ (CDS_len optional)_: Position and number of coding bases (one based includes START and STOP codons).  |
| `proteinPosition`, `proteinLength`               | 13    | AA.pos / AA.length       | 14    | Protein_position      | Protein_position _/ (Protein_len optional)_: Position and number of AA (one based, including START, but not STOP).    |
| `distance`                                       | 14    | Distance                 | 18    | Distance              |                                                                                                                       |
| `messages`                                       | 15    | ERRORS / WARNINGS / INFO |       |                       |                                                                                                                       |
| `referenceProteinSequence`                       |       |                          | 15    | Amino_acids           | Amino_acids : Reference and variant amino acids; new fields `referenceProteinSequence`, `alternateProteinSequence`    |
| `alternateProteinSequence`                       |       |                          | 15    | Amino_acids           | Amino_acids : Reference and variant amino acids; new fields `referenceProteinSequence`, `alternateProteinSequence`    |
| `referenceCodingSequence`                        |       |                          | 16    | Codons                | Codons : Reference and variant codon sequence; new fields `referenceCodingSequence`, `alternateCodingSequence`        |
| `alternateCodingSequence`                        |       |                          | 16    | Codons                | Codons : Reference and variant codon sequence; new fields `referenceCodingSequence`, `alternateCodingSequence`        |
|                                                  |       |                          | 17    | Existing_variation    | Existing_variation : Identifier(s) of co-located known variants                                                       |
| `strand`                                         |       |                          | 19    | STRAND                | STRAND : Strand of the feature (1/-1); new field `strand`                                                             |
|                                                  |       |                          | 20    | FLAGS                 | FLAGS: Transcript quality flags                                                                                       |
|                                                  |       |                          | 21    | SYMBOL_SOURCE         |                                                                                                                       |
|                                                  |       |                          | 22    | HGNC_ID               |                                                                                                                       |
