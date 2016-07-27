Variants
===

### Variant

| bdg-formats               | VCF 4.3 specification                                                         | GA4GH schema                                               |
|---------------------------|-------------------------------------------------------------------------------|------------------------------------------------------------|
| `contigName`              | column 1 "CHROM"                                                              | `Variant.reference_name`                                   |
| `start`                   | column 2 "POS", converted to 0-based coordinate system, closed-open intervals | `Variant.start`                                            |
| `end`                     | `start + referenceAllele.length()`                                            | `Variant.end`                                              |
| `names`                   | column 3 "ID", shared across all alleles in the same VCF record               | `Variant.names`                                            |
| `referenceAllele`         | column 4 "REF"                                                                | `Variant.reference_bases`                                  |
| `alternateAllele`         | column 5 "ALT", split for multi-allelic sites                                 | `Variant.alternate_bases`, split for multi-allelic sites   |
| `filtersApplied`          | column 7 "FILTER", true if any value other than the missing value             | `Variant.filters_applied`                                  |
| `filtersPassed`           | column 7 "FILTER", true if `"PASS"`                                           | `Variant.filters_passed`                                   |
| `filtersFailed`           | column 7 "FILTER", shared across all alleles in the same VCF record           | `Variant.filters_failed`                                   |
| `annotation`              |                                                                               |                                                            |


### VariantAnnotation

| bdg-formats                        | VCF 4.3 specification                                                                                      | GA4GH schema |
|------------------------------------|------------------------------------------------------------------------------------------------------------|--------------|
| `ancestralAllele`                  | INFO reserved key AA, Number=1, shared across all alternate alleles in the same VCF record                 | |
| `alleleCount`                      | INFO reserved key AC, Number=A, split for multi-allelic sites into a single integer value                  | |
| `readDepth`                        | INFO reserved key AD, Number=R, split for multi-allelic sites into a single integer value over two fields  | |
| `forwardReadDepth`                 | INFO reserved key ADF, Number=R, split for multi-allelic sites into a single integer value over two fields | |
| `reverseReadDepth`                 | INFO reserved key ADR, Number=R, split for multi-allelic sites into a single integer value over two fields | |
| `referenceReadDepth`               | INFO reserved key AD, Number=R, split for multi-allelic sites into a single integer value over two fields  | |
| `referenceForwardReadDepth`        | INFO reserved key ADF, Number=R, split for multi-allelic sites into a single integer value over two fields | |
| `referenceReverseReadDepth`        | INFO reserved key ADR, Number=R, split for multi-allelic sites into a single integer value over two fields | |
| `alleleFrequency`                  | INFO reserved key AF, Number=A, split for multi-allelic sites into a single float value                    | |
| `cigar`                            | INFO reserved key CIGAR, Number=A, split for multi-allelic sites into a single string value                | |
| `dbSnp`                            | INFO reserved key DB, Number=0, shared across all alternate alleles in the same VCF record                 | |
| `hapMap2`                          | INFO reserved key H2, Number=0, shared across all alternate alleles in the same VCF record                 | |
| `hapMap3`                          | INFO reserved key H3, Number=0, shared across all alternate alleles in the same VCF record                 | |
| `validated`                        | INFO reserved key VALIDATED, Number=0, shared across all alternate alleles in the same VCF record          | |
| `thousandGenomes`                  | INFO reserved key 1000G, Number=0, shared across all alternate alleles in the same VCF record              | |
| `somatic`                          | INFO reserved key SOMATIC, Number=0, shared across all alternate alleles in the same VCF record            | |
| `transcriptEffects`                | INFO key ANN, split for multi-allelic sites                                                                | |
| `attributes`                       | ^0 | |

Note 0:

* INFO key values with Number=., Number=0, Number=1, and Number=[n] are shared across all alternate alleles in the same VCF record
* INFO key values with Number=A are split for multi-allelic sites into a single value
* INFO key values with Number=R are split into an array of two values, [reference allele, alternate allele], separated by commas, e.g. `"0,1"`


### TranscriptEffect

| bdg-formats       | VCF ANN 1.0 specification                                | GA4GH schema                                                                             |
|-------------------|----------------------------------------------------------|------------------------------------------------------------------------------------------|
| `effects`         | Annotation (a.k.a. effect or consequence) field          | `TranscriptEffect.effects`                                                               |
| `geneName`        | Gene name field                                          | |
| `geneId`          | Gene identifier field                                    | |
| `featureType`     | Feature type field                                       | |
| `featureId`       | Feature identifier field                                 | `TranscriptEffect.feature_id`                                                            |
| `biotype`         | Transcript biotype field                                 | |
| `rank`            | Exon or intron rank, from Rank/total field               | |
| `total`           | Total number of exons or introns, from Rank/total field  | |
| `genomicHgvs`     |                                                          | `TranscriptEffect.hgvs_annotation` &rarr; `HGVSAnnotation.genomic`                       |
| `transcriptHgvs`  | HGVS.c field                                             | `TranscriptEffect.hgvs_annotation` &rarr; `HGVSAnnotation.transcript`                    |
| `proteinHgvs`     | HGVS.p field                                             | `TranscriptEffect.hgvs_annotation` &rarr; `HGVSAnnotation.protein`                       |
| `cdnaPosition`    | cDNA position from cDNA_position/cDNA_len field          | `TranscriptEffect.cdna_location` &rarr; `AlleleLocation.start`                           |
| `cdnaLength`      | cDNA length from cDNA_position/cDNA_len field            | `TranscriptEffect.cdna_location` &rarr; `(AlleleLocation.end - AlleleLocation.start)`    |
| `cdsPosition`     | CDS position from CDS_position/CDS_position field        | `TranscriptEffect.cds_location` &rarr; `AlleleLocation.start`                            |
| `cdsLength`       | CDS length from CDS_position/CDS_position field          | `TranscriptEffect.cds_location` &rarr; `(AlleleLocation.end - AlleleLocation.start)`     |
| `proteinPosition` | Protein position from Protein_position/Protein_len field | `TranscriptEffect.protein_location` &rarr; `AlleleLocation.start`                        |
| `proteinLength`   | Protein length from Protein_position/Protein_len field   | `TranscriptEffect.protein_location` &rarr; `(AlleleLocation.end - AlleleLocation.start)` |
| `distance`        | Distance field                                           | |
| `messages`        | Errors, warnings, or information messages field          | |


### Specification documents

 * [VCF 4.2 pdf](http://samtools.github.io/hts-specs/VCFv4.2.pdf)
 * [VCF 4.3 pdf](http://samtools.github.io/hts-specs/VCFv4.3.pdf)
 * [VCF 4.3 latest draft (tex)](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.tex)
 * [VCF ANN 1.0 pdf](http://snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf)
 * [Table with VCF INFO reserved keys (pull request)](https://github.com/samtools/hts-specs/pull/123)
 * [Issue regarding new AD/ADF/ADR INFO reserved keys in VCF 4.3](https://github.com/samtools/hts-specs/issues/78)
 * [TCGA VCF Specification](https://wiki.nci.nih.gov/display/TCGA/TCGA+Variant+Call+Format+%28VCF%29+Specification)
 * [NCBI human variation sets in VCF format](http://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)
 * [NCBI variation database glossary](http://www.ncbi.nlm.nih.gov/variation/docs/glossary/)
 * [GA4GH allele annotations schema](https://github.com/ga4gh/schemas/blob/master/src/main/proto/ga4gh/allele_annotations.proto)
 * [GA4GH variants schema](https://github.com/ga4gh/schemas/blob/master/src/main/proto/ga4gh/variants.proto)


## Draft notes

### Variant

| bdg-formats               | VCF 4.2 specification              | VCF 4.3 specification              | GA4GH variants specification   |
|---------------------------|------------------------------------|------------------------------------|--------------------------------|
| `contigName`              | column 1 "CHROM"                   | column 1 "CHROM"                   | `Variant.reference_name`       |
| `start`                   | column 2 "POS" ^1                  | column 2 "POS" ^1                  | `Variant.start`                |
| `end`                     | `start + referenceAllele.length()` | `start + referenceAllele.length()` | `Variant.end`                  |
| `names`                   | column 3 "ID" ^2                   | column 3 "ID" ^2                   | `Variant.names`                |
| `referenceAllele`         | column 4 "REF"                     | column 4 "REF"                     | `Variant.reference_bases`      |
| `alternateAllele`         | column 5 "ALT" ^3                  | column 5 "ALT" ^3                  | `Variant.alternate_bases` ^3   |
| `somatic`                 | INFO reserved key `SOMATIC`        | INFO reserved key `SOMATIC`        | |


#### Review questions

- [X] Is 1..1 `Variant` to `StructuralVariant` the right association? Yes
- [X] Include a mapping for `ID` on `Variant`
- [X] We should not include a mapping for `QUAL` on `Variant` or on `VariantAnnotation`
- [X] Should we include a mapping for `FILTER`, on `Variant` or on `VariantAnnotation`?  No keep it on`Genotype` and merge.
- [X] Are we implementing `variantErrorProbability` correctly?  Not useful, we will remove it
- [X] `Variant` is the best place for INFO key `SOMATIC` mapping?  Keep for now.
- [ ] Review how `QUAL` field is used in gVCF


Note 1:

VCF format uses 1-based positions, whereas bdg-formats uses 0-based, closed-open intervals.

Thus when reading from VCF format, `start` and `end` will be converted to 0-based, closed-open intervals.


Note 2:

Because bdg-formats splits one row in VCF into multiple `Variant`s, one per ALT allele, the value of
the `ID` column is copied to each `Variant`.


Note 3:

bdg-formats splits one row in VCF into multiple `Variant`s, one per ALT allele.



VCF 4.3 specification columns not present as fields of `Variant`:

```
QUAL - quality: Phred-scaled quality score for the assertion made in
ALT. i.e. −10log10 prob(call in ALT is wrong). If ALT is ‘.’ (no
variant) then this is −10log10 prob(variant), and if ALT is not ‘.’
this is −10log10 prob(no variant). If unknown, the missing value
should be specified.
```

Freebayes QUAL definition:
```
Of primary interest to most users is the QUAL field, which estimates
the probability that there is a polymorphism at the loci described by
the record.

In freebayes, this value can be understood as 1 - P(locus is homozygous given the data).
```


### StructuralVariant


Note `StructuralVariant` and `StructuralVariantType` have been removed.


| bdg-formats   | VCF 4.2 specification                     | VCF 4.3 specification                     |
|---------------|-------------------------------------------|-------------------------------------------|
| `type`        | INFO reserved key "SVTYPE" ^4             | INFO reserved key "SVTYPE" ^4             |
| `precise`     | negation of INFO reserved key "IMPRECISE" | negation of INFO reserved key "IMPRECISE" |
| `startWindow` | column 2 "POS" ^1                         | column 2 "POS" ^1                         |
| `endWindow`   | ^5                                        | ^5                                        |


#### Review questions

- [X] Is `assembly` field useful?  No, remove it.
- [X] Should we add `type`s for `BND` and `CNV`?  Yes.
- [X] Are there other structural variant INFO keys or FORMAT keys that should be fields?  No, shove into `StructuralVariant.attributes`.
- [X] `StructuralVariant` should have its own `attributes` key-value field.  [mlh] I'm less sure of this now, seems like much messy for little gain.
- [ ] Update `startWindow` as structural variant INFO reserved key "CIPOS" and `endWindow` as structural variant INFO reserved key "CIEND"?  Note both are `Number=2` fields.
- [ ] There isn't much here any more, is it worth keeping around?  (in fact, I'm not seeing any usage in ADAM code)

Note 4:

Mapping between INFO reserved key "SVTYPE" values and bdg-formats `type` field:

| value      | bdg-formats `type`                           |
|------------|----------------------------------------------|
| BND        | `StructuralVariantType.BREAKEND`             |
| CNV        | `StructuralVariantType.COPY_NUMBER_VARIABLE` |
| DEL        | `StructuralVariantType.DELETION`             |
| DEL:ME     | `StructuralVariantType.MOBILE_DELETION`      |
| DUP        | `StructuralVariantType.DUPLICATION`          |
| DUP:TANDEM | `StructuralVariantType.TANDEM_DUPLICATION`   |
| INS        | `StructuralVariantType.INSERTION`            |
| INS:ME     | `StructuralVariantType.MOBILE_INSERTION`     |
| INV        | `StructuralVariantType.INVERSION`            |


Note 5:

If `precise`, use `start + referenceAllele.length()`

Otherwise, use INFO reserved key "END", converted to 0-based, closed-open intervals


VCF 4.3 specification INFO reserved keys used for structural variants but not present as fields of `StructuralVariant`:
```
##INFO=<ID=NOVEL,Number=0,Type=Flag,Description="Indicates a novel structural variation">

##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">

##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
##INFO=<ID=METRANS,Number=4,Type=String,Description="Mobile element transduction info of the form CHR,START,END,POLARITY">
##INFO=<ID=DGVID,Number=1,Type=String,Description="ID of this element in Database of Genomic Variation">
##INFO=<ID=DBVARID,Number=1,Type=String,Description="ID of this element in DBVAR">
##INFO=<ID=DBRIPID,Number=1,Type=String,Description="ID of this element in DBRIP">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around the inserted material between breakends">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth of segment containing breakend">
##INFO=<ID=DPADJ,Number=.,Type=Integer,Description="Read Depth of adjacency">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number of segment containing breakend">
##INFO=<ID=CNADJ,Number=.,Type=Integer,Description="Copy number of adjacency">
##INFO=<ID=CICN,Number=2,Type=Integer,Description="Confidence interval around copy number for the segment">
##INFO=<ID=CICNADJ,Number=.,Type=Integer,Description="Confidence interval around copy number for the adjacency">
```

Note `DP` as defined here collides with its other definition, "combined depth across samples", which is mapped to
the `combinedDepth` field of `VariantAnnotation`.


VCF 4.3 specification FORMAT reserved keys used for structural variants but not present as fields of `StructuralVariant`:
```
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
##FORMAT=<ID=CNL,Number=G,Type=Float,Description="Copy number genotype likelihood for imprecise events">
##FORMAT=<ID=CNP,Number=G,Type=Float,Description="Copy number posterior probabilities">
##FORMAT=<ID=NQ,Number=1,Type=Integer,Description="Phred style probability score that the variant is novel">
##FORMAT=<ID=HAP,Number=1,Type=Integer,Description="Unique haplotype identifier">
##FORMAT=<ID=AHAP,Number=1,Type=Integer,Description="Unique identifier of ancestral haplotype">
```


### VariantAnnotation

| bdg-formats               | VCF 4.2 specification         | VCF 4.3 specification         | htsjdk VCFConstants | htsjdk VCFStandardHeaderLine |
|---------------------------|-------------------------------|-------------------------------|---------------------|------------------------------|
| `variant`                 |                               |                               | | |
| `ancestralAllele`         | INFO reserved key `AA`        | INFO reserved key `AA`        | `ANCESTRAL_ALLELE_KEY` | |
| `alleleCount`             | INFO reserved key `AC`        | INFO reserved key `AC`        | `ALLELE_COUNT_KEY` | X |
| `readDepth`               |                               | INFO reserved key `AD`        | | |
| `forwardReadDepth`        |                               | INFO reserved key `ADF`       | | |
| `reverseReadDepth`        |                               | INFO reserved key `ADR`       | | |
| `alleleFrequency`         | INFO reserved key `AF`        | INFO reserved key `AF`        | `ALLELE_FREQUENCY_KEY` | X |
| `cigar`                   | INFO reserved key `CIGAR`     | INFO reserved key `CIGAR`     | `CIGAR_KEY` | |
| `dbSnp`                   | INFO reserved key `DB`        | INFO reserved key `DB`        | `DBSNP_KEY` | X |
| `hapMap2`                 | INFO reserved key `H2`        | INFO reserved key `H2`        | `HAPMAP2_KEY`| |
| `hapMap3`                 | INFO reserved key `H3`        | INFO reserved key `H3`        | `HAPMAP3_KEY`| |
| `validated`               | INFO reserved key `VALIDATED` | INFO reserved key `VALIDATED` | `VALIDATED_KEY` | |
| `thousandGenomes`         | INFO reserved key `1000G`     | INFO reserved key `1000G`     | `THOUSAND_GENOMES_KEY` | |
| `transcriptEffects`       | INFO key `ANN`                | INFO key `ANN`                | | |
| `attributes`              |                               |                               | | |


Commented out attributes:

| bdg-formats                | VCF 4.2 specification         | VCF 4.3 specification         | htsjdk VCFConstants | htsjdk VCFStandardHeaderLine |
|----------------------------|-------------------------------|-------------------------------|---------------------|------------------------------|
| `allelesInCalledGenotypes` |                               | INFO reserved key `AN`        | | |
| `rmsBaseQuality`           | INFO reserved key `BQ`        | INFO reserved key `BQ`        | `RMS_BASE_QUALITY_KEY` | |
| `combinedDepth`            | INFO reserved key `DP`        | INFO reserved key `DP`        | `DEPTH_KEY` | X |
| `rmsMappingQuality`        | INFO reserved key `MQ`        | INFO reserved key `MQ`        | `RMS_MAPPING_QUALITY_KEY` | X |
| `mappingQualityZeroReads`  | INFO reserved key `MQ0`       | INFO reserved key `MQ0`       | `RMS_MAPPING_QUALITY_ZERO_KEY` | X |
| `samplesWithData`          | INFO reserved key `NS`        | INFO reserved key `NS`        | `SAMPLE_NUMBER_KEY` | |
| `strandBias`               | INFO reserved key `SB`        | INFO reserved key `SB`        | `STRAND_BIAS_KEY` | X |


#### Review questions

- [ ] `DP` vs. `AD/ADF/ADR` for read depth
- [ ] Are all of these INFO reserved key values useful as fields?
- [ ] Should we attempt to parse non-`ANN` spec INFO keys into transcript effects (e.g. `EFF`)?
- [ ] Should we define our own standard header lines in ADAM code where they are missing from htsjdk or push upstream?
- [ ] Are there other non-reserved INFO key values that should be fields?  (note there were some previously that have since been removed)
- [ ] Where INFO keys are `Number=R`, `Number=A`, or `Number=G`, should those be stored on `VariantAnnotation` or `GenotypeAnnotation`?
- [ ] Do any INFO key values need to be recalculated/split based on our split multi-allelic model?


Other commonly used VCF INFO keys:

TCGA:
```
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=GENE,Number=.,Type=String,Description="HUGO gene symbol or Unknown">
##INFO=<ID=RE,Number=0,Type=Flag,Description="Position known to have RNA-edits to occur">
##INFO=<ID=RGN,Number=.,Type=String,Description="Region where nucleotide variant occurs in relation to a gene">
##INFO=<ID=SID,Number=.,Type=String,Description="Unique identifier from gene annotation source or unknown">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=VLS,Number=1,Type=Integer,Description="Final validation status relative to non-adjacent Normal,0=none,1=germline,2=somatic,3=LOH,4=post transcriptional modification,5=unknown">
##INFO=<ID=VLSC,Number=1,Type=Integer,Description="Final somatic score between 0 and 255 when multiple lines of evidence are available"
##INFO=<ID=VT,Number=1,Type=String,Description="Variant type, can be SNP, INS , DEL, DNP, TNP, ONP or Consolidated">
```

NCBI:
```
ASS - VCF INFO. Indicates the variant is located in an acceptor splice site. The functional code (FxnCode) of such a variation = 73.

ASP - VCF INFO. Indicates the variant "is assembly specific". This flag is set if the variant maps to only one assembly.

CAF - VCF INFO. Comma delimited list of allele frequencies based on 1000Genomes.  The first frequency refers to the reference base,
  and alternate alleles follow in the order as in the ALT column.  Where 1000Genomes alternate allele is not in the dbSNPs alternate
  label set, the allele is added to the ALT column.  The minor allele is the second largest value in the list, and was previously
  reported as the GMAF in VCF, RefSNP and EntrezSNP pages and VariationReporter. This GMAF will be called 1000G MAF starting from b138.

CDA - VCF INFO. Indicates the variation is interrogated in a clinical diagnostic assay.

CFL - VCF INFO. Indicates the variant has an assembly conflict.  This flag is set for weight 1 and 2 variants that map to different
  chromosomes on different assemblies.

CLNACC - VCF INFO.  A string that is the accession and version number assigned by ClinVar to the genotype/phenotype relationship.

CLNALLE - VCF INFO. An integer that defines the alleles in the  REF or ALT columns.  0 is REF, 1 is the first ALT allele, etc. A
  value of -1 indicates that no allele was found to match a corresponding HGVS allele name.

CLNCUI - VCF INFO. A string that is the disease concept ID used in GTR and ClinVar for a phenotype associated with an allele.

CLNDBN - VCF INFO. A string that is the disease name used by the database specified by CLNSRC.

CLNHGVS - VCF INFO. A string that describes the variant names from HGVS. The order of these variants corresponds to the order of
  the info in the other clinical (CLN) INFO tags.

CLNORIGIN - VCF INFO. A string that describes the origin of the variant allele. One or more of the values may be added: 0 - unknown;
  1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental;
  256 - not-tested; 512 - tested-inconclusive; 1073741824 - other

CLNSIG - VCF INFO.  A string that describes the variant's clinical significance, where  0 - unknown, 1 - untested, 2 - non-pathogenic,
  3 - probable-non-pathogenic, 4 - probable-pathogenic, 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other.

CLNSRC - VCF INFO.  A string that describes the Variant's Clinical Sources.

CLNSRCID - VCF INFO. The identifier used for the allele from the source defined in CLNSRC.

dbSNPBuildID - VCF INFO. An integer that indicates the first dbSNP Build in which the rs was reported

DSS - VCF INFO. The variant is located in an donor splice site. The functional code (FxnCode) of such a variation = 75.

G5 - VCF INFO. Indicates the variant has a >5% minor allele frequency in one or more  populations.

G5A - VCF INFO. Indicates the variant has a >5% minor allele frequency in each and all populations.

GCF - VCF INFO. Indicates the variant has a Genotype Conflict  Same (rs, ind), different genotype.  N/N is not included.

GENEINFO - VCF INFO. Report of the gene symbol(s) and NCBI GeneID(s) at the location of the variation.  The gene symbol and
  ID are delimited by a colon (:) and each pair is delimited by a vertical bar (|).  Example: SYMBOl1:GeneID1|SYMBOl2:GeneID2| .

GMAF - VCF INFO. A floating decimal that indicates the Global Minor Allele Frequency [0, 0.5]. The global population is
  1000GenomesProject phase 1 genotype data from 629 individuals, released in the 11-23-2010 dataset.

GNO - VCF INFO. Indicates there are genotypes available for the variant and that the variant has an individual genotype (in SubInd table).

HD - VCF INFO. Indicates a marker for the variation is in a high density genotyping kit (50K density or greater);  The variant may
  have phenotype associations present in dbGaP.

INT - VCF INFO. Indicates the variant is located in an intron based on NCBI's annotation. The functional code (FxnCode) of
  such a variation = 6.

KGPhase1 - VCF INFO. Indicates the variation was in 1000 Genome phase 1 (including June Interim phase 1).

KGPilot123 - VCF INFO. Indicates the variation was in the 1000 Genome discovery -- all 2010 pilots (1,2,3).

KGPROD - VCF INFO. Indicates the variation was submitted as part of the 1000 Genomes Project.

KGValidated - VCF INFO. Indicates the variation was validated by 1000 Genomes.

LSD - VCF INFO. Indicates the variation was submitted from a locus-specific database (LSDB).

MTP - VCF INFO. Indicates the variation has microattribution/third-party annotation(TPA:GWAS,PAGE). x

MUT - VCF INFO. Indicates the allele has a low frequency and is cited in a journal or other reputable sources.

NOC - VCF INFO. Indicates the allele in the genome is not present in variant allele list. The reference sequence allele
  at the mapped position is not present in the variant allele list, adjusted for orientation.

NOV - VCF INFO. Indicates the rs cluster has non-overlapping allele sets. This is true when an rs set has more than 2 alleles
  from different submissions and these sets share no alleles in common.

NSF - VCF INFO. The consequence of the variation is a non-synonymous frameshift -- a coding region variation where one allele
  in the set changes all downstream amino acids. The functional class (FxnClass) of such a variation = 44.

NSM - VCF INFO. The consequence of the variation is a non-synonymous (missense) change -- it is a coding region variation where
  one allele in the set changes the amino acid, but translation continues. The functional class (FxnClass) of such a variation = 42.

NSN - VCF INFO. The consequence of the variation is a non-synonymous stop codon (nonsense) -- it is a coding region variation where
  one allele in the set changes to a STOP codon (TER or *). The functional class (FxnClass) of such a variation = 41.

OM - VCF INFO. Indicates the variation has a record in OMIM or OMIA.

OTH - VCF INFO. Indicates there is another variant with exactly the same set of mapped positions on NCBI reference assembly.

OTHERKG - VCF INFO. Indicates the variation was not included in the 1000 Genome submission.

PH3 - VCF INFO. Indicates the variation was HAP_MAP Phase 3 genotyped: filtered, non-redundant.

PM - VCF INFO. Indicates that the variant is from a clinical channel or is cited in PubMed (PM).

PMC - VCF INFO. Indicates that links exist from the variant's rs record to a PubMed Central article.

POPFREQ - VCF INFO. A string that gives the frequencies and count of the ALT alleles by population ID.  The form is
  p(na/ns):f(c1/c2)[|f(c1/c2)]...  where: p is the pop id; na is the number of alleles for the population; ns is the number of
  samples for the population; f is the frequency; c1 is the allele count; and c2 is sample count for that allele
  (c1 - homozygous count).  The populations ID, names, and handles shown above in dbSNP_POP_IDS, dbSNP_LOC_POP_IDS,
  and dbSNP_POP_HANDLES, respectively in corresponding order.

R3 - VCF INFO. Indicates the variant is located in a 3' gene region. The functional code (FxnCode) of such a variation = 13.

R5 - VCF INFO. Indicates the variant is located in a 5' gene region. The functional code (FxnCode) of such a variation = 15.

REF - VCF INFO. Indicates the variant "has reference". That is, it is a coding region variation where one allele in the set
  is identical to the reference sequence. The functional code (FxnCode) of such a variation = 8.

RSPOS - VCF INFO. Chromosome position as reported in dbSNP.

RV - VCF INFO. Indicates that the variation's "RS orientation is reversed".

S3D - VCF INFO. Indicates that the variant has 3D structure: SNP3D table.  Note: "S3D" is "SNP 3 dimensional".  We have
  changed the usage of "SNP" to the more inclusive term "variant"; the tag "S3D" remains in the vcf files.

SAO - VCF INFO. An integer that indicates variant allele origin.  The accepted values for this tag are: 0 - unspecified,
  1 - Germline, 2 - Somatic, 3 - Both.  Note: "SAO" is "SNP Allele Origin".  We have changed "SNP" to the more inclusive term
  "variant"; the tag "SAO" remains in the vcf files.

SLO - VCF INFO. Indicates that the variant's rs has a submitter provided "LinkOut" to the submitter's web site.

SSR - VCF INFO.  Tag that can be found in all records whose value is an integer that indicates the variant suspect reason code.
  The accepted values for this tag are: 0 - unspecified, 1 - Paralog, 2 - byEST, 3 - Para_EST, 4 - oldAlign, 5 - other.  Note:
  "SSR" is "SNP Suspect Reason".  We have changed "SNP" to the more inclusive term "variant"; the tag "SSR" remains in the vcf files.

SYN - VCF INFO. indicates the variant "has synonymous". That is, it is a coding region variation where one allele in the set
  does not change the encoded amino acid. The functional code (FxnCode) of such a variation = 3.

TPA - VCF INFO. The variant has provisional Third Party Annotation (TPA). This set is currently restricted to rs from PharmGKB,
  which provides phenotype data.

U3 - VCF INFO. The variant is located in a 3' untranslated region (UTR). The functional code (FxnCode) of such a variation = 53.

U5 - VCF INFO. The variant is located in a 5' untranslated region (UTR). The functional code (FxnCode) of such a variation = 55.

VLD - VCF INFO. The variant is validated. This flag is set if a variant has 2+ minor allele count based on frequency or genotype data.

VP - VCF INFO. A string that describes a "Variation Property" of the variant. 

VC - VCF INFO. A string that describes the "Variation Class" of the variant.

WGT - VCF INFO. A integer that indicates variant map weight. The accepted values for this tag are: 00 - unmapped, 1 - map weight 1;
  2 - map weight 2;  3 - map weight 3 or more.

WTD - VCF INFO. Indicates if one submitter of information about variation at this location was "withdrawn by submitter". The flag is
  set if one ss member of the rs cluster for a variant is withdrawn by the submitter.  If all ss members of the rs cluster are withdrawn,
  then the variant rs would be deleted to "SNPHistory" and not reported.
```


### TranscriptEffect

| bdg-formats       | VCF ANN 1.0 specification                                | GA4GH allele annotations specification                                                   |
|-------------------|----------------------------------------------------------|------------------------------------------------------------------------------------------|
| `effects`         | Annotation (a.k.a. effect or consequence) field          | `TranscriptEffect.effects`                                                               |
| `geneName`        | Gene name field                                          | |
| `geneId`          | Gene identifier field                                    | |
| `featureType`     | Feature type field                                       | |
| `featureId`       | Feature identifier field                                 | `TranscriptEffect.feature_id`                                                            |
| `biotype`         | Transcript biotype field                                 | |
| `rank`            | Exon or intron rank, from Rank/total field               | |
| `total`           | Total number of exons or introns, from Rank/total field  | |
| `genomicHgvs`     |                                                          | `TranscriptEffect.hgvs_annotation` &rarr; `HGVSAnnotation.genomic`                       |
| `transcriptHgvs`  | HGVS.c field                                             | `TranscriptEffect.hgvs_annotation` &rarr; `HGVSAnnotation.transcript`                    |
| `proteinHgvs`     | HGVS.p field                                             | `TranscriptEffect.hgvs_annotation` &rarr; `HGVSAnnotation.protein`                       |
| `cdnaPosition`    | cDNA position from cDNA_position/cDNA_len field          | `TranscriptEffect.cdna_location` &rarr; `AlleleLocation.start`                           |
| `cdnaLength`      | cDNA length from cDNA_position/cDNA_len field            | `TranscriptEffect.cdna_location` &rarr; `(AlleleLocation.end - AlleleLocation.start)`    |
| `cdsPosition`     | CDS position from CDS_position/CDS_position field        | `TranscriptEffect.cds_location` &rarr; `AlleleLocation.start`                            |
| `cdsLength`       | CDS length from CDS_position/CDS_position field          | `TranscriptEffect.cds_location` &rarr; `(AlleleLocation.end - AlleleLocation.start)`     |
| `proteinPosition` | Protein position from Protein_position/Protein_len field | `TranscriptEffect.protein_location` &rarr; `AlleleLocation.start`                        |
| `proteinLength`   | Protein length from Protein_position/Protein_len field   | `TranscriptEffect.protein_location` &rarr; `(AlleleLocation.end - AlleleLocation.start)` |
| `distance`        | Distance field                                           | |
| `messages`        | Errors, warnings, or information messages field          | |


#### Review questions

- [ ] Add analysis result per GA4GH and Ensembl VEP REST APIs?


GA4GH allele annotations specification fields not present as fields on `TranscriptEffect`:

```
message TranscriptEffect {
  ...

// Output from prediction packages such as SIFT.
  repeated AnalysisResult analysis_results = 9;
}

// An AnalysisResult record holds the output of a prediction package such as
// SIFT on a specific allele.
message AnalysisResult {
  // The ID of the analysis record for this result
  string analysis_id = 1;

  // The text-based result for this analysis
  string result = 2;

  // The numeric score for this analysis
  int32 score = 3;
}
```

### Genotype

| bdg-formats               | VCF 4.2 specification                       | VCF 4.3 specification                       | GA4GH variants specification        |
|---------------------------|---------------------------------------------|---------------------------------------------|-------------------------------------|
| `variant`                 | | | |
| `contigName`              | column 1 "CHROM"                            | column 1 "CHROM"                            | `Variant.reference_name`            |
| `start`                   | column 2 "POS" ^1                           | column 2 "POS" ^1                           | `Variant.start`                     |
| `end`                     | `start + variant.referenceAllele.length()`  | `start + variant.referenceAllele.length()`  | `Variant.end`                       |
| `genotypeAnnotation`      | | | |
| `sampleId`                | Header line sample identifier               | Header line sample identifier               | `CallSet.bio_sample_id` ^6          |
| `alleles`                 | Genotype field reserved key "GT" ^7         | Genotype field reserved key "GT" ^7         | `Call.genotype` ^7                  |
| `expectedAlleleDosage`    | | | |
| `referenceReadDepth`      | INFO reserved key "AD" ^8                   | INFO reserved key "AD" ^8                   | |
| `alternateReadDepth`      | INFO reserved key "AD" ^8                   | INFO reserved key "AD" ^8                   | |
| `readDepth`               | Genotype field reserved key "DP"            | Genotype field reserved key "DP"            | |
| `minReadDepth`            | | | |
| `genotypeQuality`         | Genotype field reserved key "GQ"            | Genotype field reserved key "GQ"            |
| `genotypeLikelihoods`     | Genotype field reserved key "PL" ^9         | Genotype field reserved key "PL" ^9         | `Call.genotype_likelihood`          |
| `nonReferenceLikelihoods` | | | |
| `strandBiasComponents`    | | | |
| `splitFromMultiAlleleic`  | `true` if this genotype was split from a multi-allelic variant | | |
| `phased`                  | `phaseSetId != null`                        | `phaseSetId != null`                        | |
| `phaseSetId`              | Genotype field reserved key "PS"            | Genotype field reserved key "PS"            | `Call.phaseset`                     |
| `phaseQuality`            | Genotype field reserved key "PQ"            | Genotype field reserved key "PQ"            | |


#### Review questions

- [X] Review sample mapping, implementation in ADAM code
- [ ] `DP` vs. `AD/ADF/ADR` for read depth
- [X] Mapping for `expectedAlleleDosage`, `strandBiasComponents`?  calculated during genotyping, add mapping to custom INFO keys
- [ ] Review `phased`, `phaseSetId` mappings, implementation in ADAM code  Find a good VCF with phased calls for unit tests.
- [ ] Are all of these FORMAT/genotype reserved key values useful as fields on `Genotype`?
- [ ] Should we define our own standard header lines where they are missing from htsjdk or push upstream?
- [ ] Are there other non-reserved FORMAT/genotype key values that should be fields?  (note there were some previously that have since been removed)
- [ ] Where INFO keys are `Number=R`, `Number=A`, or `Number=G`, should those be stored on `VariantAnnotation` or `GenotypeAnnotation`?
- [ ] Do any (other) FORMAT/genotype key values need to be recalculated/split based on our split multi-allelic model?


Note 6:

The join between `Call` and `CallSet` appears to be either on `Call.call_set_name` to `CallSet.name`
or `Call.call_set_id` to `CallSet.id`.

The join between bdg-formats `Genotype` and `Sample` is on `Genotype.sampleId` to `Sample.sampleId`.


Note 7:

Genotype field reserved key "GT" maps directly to GA4GH `Call.genotype`.  bdg-formats splits a multi-allelic
genotype into two genotype records.  For each index in the "GT" field array:

If the index is zero, use `GenotypeAllele.REF`

If the index matches the alternate allele of the variant for this genotype, use `GenotypeAllele.ALT`

If the index does not match the alternate allele of the variant for this genotype, use `GenotypeAllele.OTHER_ALT`

If the index is the no-call symbol `.`, use `GenotypeAllele.NO_CALL`

Ordering of the `GenotypeAllele` values in the `alleles` array is important if `phased` is true.

What about if the GT field shows phase (e.g. `0|1`) but `phaseSetId` is unset?  Per VCF spec,  "All phased genotypes
that do not contain a PS subfield are assumed to belong to the same phased set."


Note 8:

INFO reserved key field "AD" is currently split for multi-allelic genotypes

```scala
if (g.hasAD) {
  val ad = g.getAD
  gb.AD(Array(ad(0), ad(idx))
}
```
and then saved to `Genotype.referenceReadDepth` and `Genotype.alternateReadDepth`.

```scala
if (g.hasAD) {
  val ad = g.getAD
  genotype.setReferenceReadDepth(ad(0)).setAlternateReadDepth(ad(1))
}
```


Note 9:

bdg-formats uses genotype field reserved key "PL" whereas GA4GH uses genotype field reserved key "GL".

The "PL" field values are currently split for multi-allelic genotypes:

```scala
// Recompute PLs as needed to reflect stripped alleles.
// TODO: Collapse other alternate alleles into a single set of probabilities.
if (g.hasPL) {
  val oldPLs = g.getPL
  val maxIdx = oldPLs.length
  val newPLs = GenotypeLikelihoods.getPLIndecesOfAlleles(0, idx).filter(_ < maxIdx).map(oldPLs(_))
  // Normalize new likelihoods in log-space
  gb.PL(newPLs.map(_ - newPLs.min))
}
```

and then saved to `Genotype.genotypeLikelihoods` and `Genotype.nonReferenceLikelihoods` after
conversion from PHRED scores to log probabilities.

```scala
if (g.hasPL) b.setGenotypeLikelihoods(g.getPL.toList.map(p => jFloat(PhredUtils.phredToLogProbability(p))))
```

We should switch to using GL instead of PL.


### GenotypeAnnotation

| bdg-formats               | VCF 4.2 specification | VCF 4.3 specification |
|---------------------------|-----------------------|-----------------------|
| `filtersFailed`           | Genotype field reserved key "FT" ^10  | Genotype field reserved key "FT" ^10  |
| `filtersPassed`           | Genotype field reserved key "FT" ^10  | Genotype field reserved key "FT" ^10  |
| `downsampled`             | | |
| `baseQualityRankSum`      | | |
| `fisherStrandBiasPValue`  | | |
| `rmsMappingQuality`       | Genotype field reserved key "MQ"     | Genotype field reserved key "MQ"     |
| `mappingQualityZeroReads` | | |
| `mappingQualityRankSum`   | | |
| `readPositionRankSum`     | | |
| `genotypePriors`          | | |
| `genotypePosteriors`      | Genotype field reserved key "GP"     | Genotype field reserved key "GP"     |
| `vqslod`                  | | |
| `culprit`                 | | |
| `attributes`              | | |


#### Review questions

- [ ] Review filter mapping, implementation in ADAM code
- [ ] Mapping for `downsampled`, `baseQualityRankSum`, `fisherStrandBiasPValue`, ... ?
- [ ] Are all of these FORMAT/genotype reserved key values useful as fields on `GenotypeAnnotation`?
- [ ] Should we define our own standard header lines where they are missing from htsjdk or push upstream?
- [ ] Are there other non-reserved FORMAT/genotype key values that should be fields?  (note there were some previously that have since been removed)
- [ ] Where INFO keys are `Number=R`, `Number=A`, or `Number=G`, should those be stored on `VariantAnnotation` or `GenotypeAnnotation`?
- [ ] Do any (other) FORMAT/genotype key values need to be recalculated/split based on our split multi-allelic model?


Note 10:

If FT = ".", then both filtersFailed and filtersPassed are empty;

If FT = "PASS", then filtersPassed should contain all FILTER IDs from header
(but there isn't any distinction between FILTER filters and FT filters in header);

Otherwise filtersPassed is empty and filtersFailed contains failed FT values


VCF 4.3 specification genotype field reserved keys not present as fields of `Genotype` or `GenotypeAnnotation`:
```
GL : genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods
for all possible genotypes given the set of alleles defined in the REF and ALT fields.

GLE : genotype likelihoods of heterogeneous ploidy, used in presence of uncertain copy number.

GQ : conditional genotype quality, encoded as a phred quality −10log10 p(genotype call is wrong,
conditioned on the site’s being variant)

HQ : haplotype qualities, two comma separated phred qualities

EC : comma separated list of expected alternate allele counts for each alternate allele in the
same order as listed in the ALT field
```

Note `GLE` has been removed from the VCF 4.3 spec as of commit https://github.com/samtools/hts-specs/commit/1c5e2787daf6d7bf97d9e4d2d77c082716351750.


Other commonly used genotype field keys:

```
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depth of reads supporting alleles 0/1/2/3...">
##FORMAT=<ID=BQ,Number=.,Type=Integer,Description="Average base quality for reads supporting alleles">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases">
##FORMAT=<ID=DPN,Number=.,Type=Integer,Description="Strand specific depth of filtered reads supporting all reported alleles: fwd0,rev0,fwd1,rev1,fwd2,rev2,etc">
##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown">
##FORMAT=<ID=SSC,Number=1,Type=Integer,Description="Somatic score between 0 and 255">
##FORMAT=<ID=TE,Number=.,Type=String,Description="Translational effect of the variant in a codon">
##FORMAT=<ID=VAQ,Number=1,Type=Integer,Description="Variant allele quality">
```

Note `AD` as defined here collides with its other definition, "read depths for each allele; total (AD)",
which is mapped to the `readDepth` field of `VariantAnnotation`.


Other VCF headers for reference:

ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf.gz
```
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth summed across all datasets, excluding MQ0 reads">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Net Genotype quality across all datasets, defined as difference between most likely and next most likely genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Net Genotype across all datasets">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set in which this variant falls">

##INFO=<ID=DPSum,Number=1,Type=Integer,Description="Total read depth summed across all datasets, excluding MQ0 reads">
##INFO=<ID=filter,Number=1,Type=String,Description="Reason for filtering this genotype as uncertain">
##INFO=<ID=platforms,Number=1,Type=Integer,Description="Number of different platforms for which at least one callset called this genotype, whether filtered or not">
##INFO=<ID=platformnames,Number=.,Type=String,Description="Names of platforms for which at least one callset called this genotype, whether filtered or not">
##INFO=<ID=platformbias,Number=.,Type=String,Description="Names of platforms that have reads containing a variant at this location, but the high-confidence call is homozygous reference, indicating that there is a potential bias.">
##INFO=<ID=datasets,Number=1,Type=Integer,Description="Number of different datasets for which at least one callset called this genotype, whether filtered or not">
##INFO=<ID=datasetnames,Number=.,Type=String,Description="Names of datasets for which at least one callset called this genotype, whether filtered or not">
##INFO=<ID=datasetsmissingcall,Number=.,Type=Integer,Description="Names of datasets that are missing a call or have an incorrect call at this location, and the high-confidence call is a variant">
##INFO=<ID=callsets,Number=1,Type=Integer,Description="Number of different callsets that called this genotype, whether filtered or not">
##INFO=<ID=callsetnames,Number=.,Type=String,Description="Names of callsets that called this genotype, whether filtered or not">
##INFO=<ID=varType,Number=1,Type=String,Description="Type of variant">
##INFO=<ID=filt,Number=1,Type=String,Description="List of callsets that had this call filtered.">
##INFO=<ID=lowcov,Number=1,Type=String,Description="List of callsets that had this call in a region with low coverage of high MQ reads.">
##INFO=<ID=arbitrated,Number=1,Type=String,Description="TRUE if callsets had discordant calls so that arbitration was needed.">
```

#### VariantDB_Challenge

https://s3-us-west-2.amazonaws.com/mayo-bic-tools/variant_miner/vcfs/1KG.chr22.anno.infocol.vcf.gz
```
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Log-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GD,Number=1,Type=Float,Description="Genotype dosage.  Expected count of non-ref alleles [0,2]">
##FORMAT=<ID=OG,Number=1,Type=String,Description="Original Genotype input to Beagle">

##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=CB,Number=.,Type=String,Description="List of centres that called, UM (University of Michigan), BI (Broad Institute), BC (Boston College), NCBI">
##INFO=<ID=EUR_R2,Number=1,Type=Float,Description="R2 From Beagle based on European Samples">
##INFO=<ID=AFR_R2,Number=1,Type=Float,Description="R2 From Beagle based on AFRICAN Samples">
##INFO=<ID=ASN_R2,Number=1,Type=Float,Description="R2 From Beagle based on Asian Samples">
##INFO=<ID=SAVANT_IMPACT,Number=1,Type=String,Description="SAVANT_IMPACT">
##INFO=<ID=SAVANT_EFFECT,Number=1,Type=String,Description="SAVANT_EFFECT">
##INFO=<ID=ENST,Number=.,Type=String,Description="Ensembl transcript identifier">
##INFO=<ID=GENE,Number=.,Type=String,Description="Ensembl gene identifier">
##INFO=<ID=INFO,Number=.,Type=String,Description="strand/sequencelength/Exons/coding region size">
##INFO=<ID=LOC,Number=.,Type=String,Description="Exon number">
##INFO=<ID=HGVS,Number=.,Type=String,Description="HGVS Nomenclature">
##INFO=<ID=CLASS,Number=.,Type=String,Description="5PU: Variant in 5. untranslated region, 3PU: Variant in 3. untranslated region, INTR: Intronic variant that does not alter splice site bases, SS:  Intronic variant that alters a splice site base at -12 to +12 but not an ESS or , SS5 base10, ESS: Variant that alters essential splice site base (+1,+2,-1,-2), SS5: Variant that alters the +5 splice site base, but not an ESS base, SC: Synonymous change caused by a base substitution (i.e. does not alter amino , acid), NSC:  Nonsynonymous  change (missense)  caused  by  a  base  substitution (i.e. , alters amino acid), IF:  In-frame  insertion  and/or  deletion  (variant  alters  the  length  of  coding , sequence but not the frame), IM: Variant that alters the start codon, SG: Variant resulting in stop-gain (nonsense) mutation, SL: Variant resulting in stop-loss mutation, FC: Frameshifting insertion and/or deletion (variant alters the length and frame , of coding sequence), SCI: Synonymous change caused by a complex indel  that does not alter protein , length or amino acid, NCI:  Nonsynonymous  change  caused  by  a  complex  indel  that  does  not  alter , protein length but does change amino acid, SCSS:  Synonymous  coding  variant that  alters the  first  or last  three  bases  of  an , exon and may affect splicing (except the  first and last  three bases of  the whole , coding sequence)., NSCSS: Nonsynonymous coding variant that alters the first or last three bases of , an exon and may affect splicing (except the first and last three bases of the whole , coding sequence)., SCISS: Synonymous change by a complex indel that alters  the  first or last  three , bases of an exon. It does not change protein length or amino acid, but may affect , splicing (except the first and last three bases of the whole coding sequence)., NCISS: Nonsynonymous  change  by  a  complex  indel  that  alters  the  first  or  last , three  bases  of  an  exon. It  does  not  change  protein  length  but  does  change  an , amino acid and may affect  splicing (except  the  first and last  three  bases  of  the , whole coding sequence)., IFSS:  In-frame coding variant  that alters the  first or last  three bases of an exon , and may affect splicing (except the first and last three bases of the whole coding , sequence).">
##INFO=<ID=ALTANN,Number=.,Type=String,Description="None: variant has the same HGVS annotation regardless of its , left or right-alignment, HGVSNotClass: indel has an alternative HGVS but the same CLASS. HGVSAndClass: Multiple HGVS with different functional consequences.">
##INFO=<ID=ALTHGVS,Number=.,Type=String,Description="Alternate HGVS Nomenclature">
##INFO=<ID=ALTCLASS,Number=.,Type=String,Description="Alternate CLASS variable">
##INFO=<ID=DBSNP,Number=.,Type=String,Description="rsID from dbSNP">
##INFO=<ID=TYPE,Number=.,Type=String,Description="SNP or Indel">

// not sure why the duplication here...

##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Log-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GD,Number=1,Type=Float,Description="Genotype dosage.  Expected count of non-ref alleles [0,2]">
##FORMAT=<ID=OG,Number=1,Type=String,Description="Original Genotype input to Beagle">

##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=CB,Number=.,Type=String,Description="List of centres that called, UM (University of Michigan), BI (Broad Institute), BC (Boston College), NCBI">
##INFO=<ID=EUR_R2,Number=1,Type=Float,Description="R2 From Beagle based on European Samples">
##INFO=<ID=AFR_R2,Number=1,Type=Float,Description="R2 From Beagle based on AFRICAN Samples">
##INFO=<ID=ASN_R2,Number=1,Type=Float,Description="R2 From Beagle based on Asian Samples">
##INFO=<ID=hgmd_2014.2.CLASS,Number=1,Type=String,Description="CLASS. This field categorizes mutations and polymorphisms. There are seven possible values, DM, DM?, DP, DFP, FP, FTV and R. Details here -> http://bsiweb.mayo.edu/catalog/hgmd">
##INFO=<ID=hgmd_2014.2.PHEN,Number=1,Type=String,Description="The name for the disease or condition associated with the mutation. If the curator had some reservations about the adequacy of the evidence supporting the relationship to the disease, this name is followed by a question mark. Disease names are provided as reported in the literature.">
##INFO=<ID=hgmd_2014.2.OMIM_ID,Number=1,Type=String,Description="Identifier for the OMIM database, http://www.ncbi.nlm.nih.gov/omim.">
##INFO=<ID=hgmd_2014.2.PubMed,Number=1,Type=String,Description="Numeric PUBMED id for the reference in PubMed.">
##INFO=<ID=dbNSFP.FathmmPred,Number=1,Type=String,Description="If a FATHMM_score is <=-1.5 the corresponding NS is predicted as \"D(AMAGING)\"; otherwise it is predicted as \"T(OLERATED)\".">
##INFO=<ID=dbNSFP.FathmmRankscore,Number=1,Type=String,Description="0 is low impact, 1 is high">
##INFO=<ID=dbNSFP.FathmmScore,Number=1,Type=String,Description="FATHMM default score">
##INFO=<ID=dbNSFP.Gerp++Nr,Number=1,Type=String,Description="GERP++ neutral rate">
##INFO=<ID=dbNSFP.Gerp++Rs,Number=1,Type=String,Description="GERP++ RS score, the larger the score, the more conserved the site">
##INFO=<ID=dbNSFP.Gerp++RsRankscore,Number=1,Type=String,Description="0 is low impact, 1 is high">
##INFO=<ID=dbNSFP.InterproDomain,Number=1,Type=String,Description="domain or conserved site on which the variant locates. Domain annotations come from Interpro database. The number in the brackets following a specific domain is the count of times Interpro assigns the variant position to that domain, typically coming from different predicting databases">
##INFO=<ID=dbNSFP.LrtPred,Number=1,Type=String,Description="LRT prediction, D(eleterious), N(eutral) or U(nknown), which is not solely determined by score">
##INFO=<ID=dbNSFP.LrtScore,Number=1,Type=String,Description="The original LRT two-sided p-value (LRTori).">
##INFO=<ID=dbNSFP.MutationassessorPred,Number=1,Type=String,Description="utationAssessor's functional impact of a variant : predicted functional, i.e. high (\"H\") or medium (\"M\"), or predicted non-functional, i.e. low (\"L\") or neutral (\"N\"). The MAori score cutoffs between \"H\"and \"M\", \"M\"and \"L\", and \"L\"and \"N\", are 3.5, 1.9 and 0.8, respectively. The rankscore cutoffs between \"H\"and \"M\", \"M\"and \"L\", and \"L\"and \"N\", are 0.9416, 0.61387 and 0.26162, respectively.">
##INFO=<ID=dbNSFP.MutationassessorRankscore,Number=1,Type=String,Description="0 is low impact, 1 is high">
##INFO=<ID=dbNSFP.MutationassessorScore,Number=1,Type=String,Description="MutationAssessor functional impact combined score (MAori). The score ranges from -5.545 to 5.975 in dbNSFP. Please refer to Reva et al. (2011)Nucl. Acids Res. 39(17):e118 for details.">
##INFO=<ID=dbNSFP.MutationtasterConvertedRankscore,Number=1,Type=String,Description="The MTori scores were first converted: if the prediction is \"A\"or \"D\"MTnew=MTori; if the prediction is \"N\"or \"P\", MTnew=1-MTori. Then MTnew scores were ranked among all MTnew scores in dbNSFP. The rankscore is the ratio of the rank of the score over the total number of MTnew scores in dbNSFP. The scores range from 0.0931 to 0.80722.">
##INFO=<ID=dbNSFP.MutationtasterPred,Number=1,Type=String,Description="MutationTaster prediction, \"A\"(\"disease_causing_automatic\"),\"D\"(\"disease_causing\"), \"N\"(\"polymorphism\") or \"P\"(\"polymorphism_automatic\"). The score cutoff between \"D\"and \"N\"is 0.5 for MTori and 0.328 for the rankscore.">
##INFO=<ID=dbNSFP.MutationtasterScore,Number=1,Type=String,Description="MutationTaster p-value (MTori), ranges from 0 to 1.">
##INFO=<ID=dbNSFP.Phastcons46wayPrimate,Number=1,Type=String,Description="phastCons conservation score based on the multiple alignments of 10 primate genomes (including human). The larger the score, the more conserved the site">
##INFO=<ID=dbNSFP.Phastcons46wayPrimateRankscore,Number=1,Type=String,Description="0 is low impact, 1 is high">
##INFO=<ID=dbNSFP.Phylop46wayPrimate,Number=1,Type=String,Description="phastCons conservation score based on the multiple alignments of 10 primate genomes">
##INFO=<ID=dbNSFP.Phylop46wayPrimateRankscore,Number=1,Type=String,Description="0 is low impact, 1 is high">
##INFO=<ID=dbNSFP.Polyphen2HdivPred,Number=1,Type=String,Description="Polyphen2 prediction based on HumDiv, \"D\"(\"probably damaging\"), \"P\"(\"possibly damaging\") and \"B\"(\"benign\"). Multiple entries separated by \";\". Because the availability of multiple values, use expression such as 'D' in Polyphen2_HDIV_pred instead of 'D' = Polyphen2_HDIV_pred to filter variants that are probably damaging.">
##INFO=<ID=dbNSFP.Polyphen2HdivRankscore,Number=1,Type=String,Description="The maximum (most damaging) value of Polyphen2 score based on HumDiv, i.e. hdiv_prob. Use Polyphen2_HDIV_score to get a list of all scores.">
##INFO=<ID=dbNSFP.Polyphen2HdivScore,Number=1,Type=String,Description="Polyphen2 score based on HumDiv, i.e. hdiv_prob. The score ranges from 0 to 1, and the corresponding prediction is \"probably damaging\"if it is in [0.957,1]; \"possibly damaging\"if it is in [0.453,0.956]; \"benign\"if it is in [0,0.452]. Score cutoff for binary classification is 0.5, i.e. the prediction is \"neutral\"if the score is smaller than 0.5 and \"deleterious\"if the score is larger than 0.5.">
##INFO=<ID=dbNSFP.Polyphen2HvarPred,Number=1,Type=String,Description="Polyphen2 prediction based on HumVar, \"D\"(\"porobably damaging\"), \"P\"(\"possibly damaging\") and \"B\"(\"benign\"). Multiple entries separated by \";\". Because the availability of multiple values, use expression such as 'D' in Polyphen2_HVAR_pred instead of 'D' = Polyphen2_HVAR_pred to filter variants that are probably damaging.">
##INFO=<ID=dbNSFP.Polyphen2HvarRankscore,Number=1,Type=String,Description="The maximum (most damaging) value of all Polyphen2 score based on HumVar, i.e. hvar_prob. Use Polyphen2_HVAR_score_all to get a list of all scores.">
##INFO=<ID=dbNSFP.Polyphen2HvarScore,Number=1,Type=String,Description="Polyphen2 score based on HumVar, i.e. hvar_prob. The score ranges from 0 to 1, and the corresponding prediction is \"probably damaging\"if it is in [0.909,1]; \"possibly damaging\"if it is in [0.447,0.908]; \"benign\"if it is in [0,0.446]. Score cutoff for binary classification is 0.5, i.e. the prediction is \"neutral\"if the score is smaller than 0.5 and \"deleterious\"if the score is larger than 0.5.">
##INFO=<ID=dbNSFP.Pos,Number=1,Type=String,Description="physical position on the chromosome as to hg19">
##INFO=<ID=dbNSFP.RadialsvmPred,Number=1,Type=String,Description="Prediction of our SVM based ensemble prediction score,\"T(olerated)\"or \"D(amaging)\". The score cutoff between \"D\"and \"T\"is 0. The rankscore cutoff between \"D\"and \"T\"is 0.83357.">
##INFO=<ID=dbNSFP.SiftScore,Number=1,Type=String,Description="SIFT score, If a score is smaller than 0.05 the corresponding NS is predicted as \"D(amaging)\"otherwise it is predicted as \"T(olerated)\".">
##INFO=<ID=dbNSFP.Siphy29wayLogodds,Number=1,Type=String,Description="The estimated stationary distribution of A, C, G and T at the site, using SiPhy algorithm based on 29 mammals genomes.">
##INFO=<ID=dbNSFP.Siphy29wayLogoddsRankscore,Number=1,Type=String,Description="0 is low impact, 1 is high">
##INFO=<ID=dbNSFP.Vest3Rankscore,Number=1,Type=String,Description="0 is low impact, 1 is high">
##INFO=<ID=dbNSFP.Vest3Score,Number=1,Type=String,Description="Score ranges from 0 to 1. The larger the score the more likely the mutation may cause functional change">
##INFO=<ID=ExAC.Info.AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=ExAC.Info.AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=ExAC.Info.AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=ExAC.Info.AC_Het,Number=1,Type=Integer,Description="Adjusted Heterozygous Count (AC_Het)">
##INFO=<ID=ExAC.Info.AC_Hom,Number=1,Type=Integer,Description="Adjusted Homozygous Count">
##INFO=<ID=Clinvar.RCVaccession,Number=1,Type=String,Description="list of RCV accessions that report this variant">
##INFO=<ID=Clinvar.ReviewStatus,Number=1,Type=String,Description="highest review status for reporting this measure. For the key to the terms, and their relationship to the star graphics ClinVar displays on its web pages, see http://www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/#interpretation">
##INFO=<ID=Clinvar.ClinicalSignificance,Number=1,Type=String,Description="character, comma-separated list of values of clinical significance reported for this variation">
##INFO=<ID=Clinvar.OtherIDs,Number=1,Type=String,Description="list of other identifiers or sources of information about this variant">
##INFO=<ID=Clinvar.Guidelines,Number=1,Type=String,Description="ACMG only right now, for the reporting of incidental variation in a Gene">
##SnpEffVersion="3.5h (build 2014-04-01), by Pablo Cingolani"
##SnpEffCmd="SnpEff  -noStats -hgvs -lof GRCh37.74 1KG.chr22.vcf000.anno"
##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )'">
##INFO=<ID=Effect_Impact,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Functional_Class,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Codon_Change,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Amino_Acid_Change,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Amino_Acid_length,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Gene_Name,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Transcript_BioType,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Gene_Coding,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Transcript_ID,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Exon_Rank,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=Genotype_Number,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=ERRORS,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=WARNINGS,Number=1,Type=String,Description="Annotation from SNPEFF">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
```

https://s3-us-west-2.amazonaws.com/mayo-bic-tools/variant_miner/gvcfs/NA12878.chr22.g.vcf.gz
```
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">

##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
```
