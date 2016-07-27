Features
===

| bdg-formats     | GFF2/GTF                      | GFF3                          | BED                            | NarrowPeak            | IntervalList               |
|-----------------|-------------------------------|-------------------------------|--------------------------------|-----------------------|----------------------------|
| `featureId`     |                               | `ID` attribute key            |                                |                       | |
| `name`          |                               | `Name` attribute key          | optional column 4 "name" ^2    | column 4 "name" ^2    | column 5 "intervalName" ^2 |
| `source`        | column 2 "source"             | column 2 "source"             |                                |                       | |
| `featureType`   | column 3 "type"               | column 3 "type"               |                                |                       | |
| `contigName`    | column 1 "seqid"              | column 1 "seqid"              | column 1 "chrom"               | column 1 "chrom"      | column 1 "sequenceName"    |
| `start`         | column 4 "start" ^1           | column 4 "start" ^1           | column 2 "chromStart"          | column 2 "chromStart" | column 2 "start"           |
| `end`           | column 5 "end" ^1             | column 5 "end" ^1             | column 3 "chromEnd"            | column 3 "chromEnd"   | column 3 "end"             |
| `strand`        | column 7 "strand" ^3          | column 7 "strand" ^3          | optional column 6 "strand" ^3  | column 6 "strand"     | column 4 "strand" ^3       |
| `phase`         |                               | column 8 "phase"              |                                |                       | |
| `frame`         | column 8 "frame"              |                               |                                |                       | |
| `score`         | column 6 "score"              | column 6 "score"              | optional column 5 "score"      | column 5 "score"      | |
| `geneId`        | `gene_id` attribute key       | `gene_id` attribute key       |                                |                       | |
| `transcriptId`  | `transcript_id` attribute key | `transcript_id` attribute key |                                |                       | |
| `exonId`        | `exon_id` attribute key       | `exon_id` attribute key       |                                |                       | |
| `aliases`       |                               | `Alias` attribute key         |                                |                       | |
| `parentIds`     |                               | `Parent` attribute key        |                                |                       | |
| `target`        |                               | `Target` attribute key        |                                |                       | |
| `gap`           |                               | `Gap` attribute key           |                                |                       | |
| `derivesFrom`   |                               | `Derives_from` attribute key  |                                |                       | |
| `notes`         |                               | `Note` attribute key          |                                |                       | |
| `dbxrefs`       |                               | `Dbxref` attribute key        |                                |                       | |
| `ontologyTerms` |                               | `Ontology_term` attribute key |                                |                       | |
| `circular`      |                               | `Is_circular` attribute key   |                                |                       | |
| `attributes`    | column 9 "attributes"         | column 9 "attributes"         | ^4                             | ^5                    | |


Note 1:

GFF2/GTF, GFF3, and IntervalList formats use 1-based, fully closed intervals, whereas bdg-formats uses 0-based, closed-open intervals.

Thus when reading from these formats, `start` and `end` will be converted to 0-based, closed-open intervals.


Note 2:

When writing BED, NarrowPeak, or IntervalList formats, the "name" or "intervalName" column is written according to the following:
  * If `name` is set, use `name`
  * If `featureId` is set, use `featureId`
  * If `featureType` is "exon" or "SO:0000147" and `exonId` is set, use `exonId`
  * If `featureType` is "transcript" or "SO:0000673" and `transcriptId` is set, use `transcriptId`
  * If `featureType` is "gene" or "SO:0000704" and `geneId` is set, use `geneId`
  * If `featureType` is set but is none of the above, use `featureType`
  * Otherwise use "sequence_feature"


Note 3:

Mapping between strand column values and bdg-formats `strand` field:

| character     | bdg-formats `strand` |
|---------------|----------------------|
| `+`           | `Strand.FORWARD` |
| `-`           | `Strand.REVERSE` |
| `.`           | `Strand.INDEPENDENT` |
| `?`           | `Strand.UNKNOWN` |
| anything else | `null` |


Note 4:

BED format (BED12, specifically) optional columns 7 through 12 are stored in
`attributes` with attribute keys `thickStart`, `thickEnd`, `itemRgb`, `blockCount`,
`blockSizes`, and `blockStarts`, respectively.


Note 5:

NarrowPeak format optional columns 7 through 10 are stored in `attributes` with
attribute keys `signalValue`, `pValue`, `qValue`, and `peak`, respectively.
