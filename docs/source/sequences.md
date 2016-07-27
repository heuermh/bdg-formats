Sequences
===

### Sequence

| bdg-formats   | biojava 1.x                                                                  | biojava 1.x biojavax                   | biojava 4.x/5.x |
|---------------|------------------------------------------------------------------------------|----------------------------------------|-----------------|
| `name`        | `Sequence.getName()`                                                         | `RichSequence.getName()`               |  tbd            |
| `description` | `Sequence.getAnnotation().getProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE)` | `RichSequence.getDescription()`        |                 |
| `alphabet`    | `Sequence.getAlphabet()`                                                     | `RichSequence.getAlphabet()`           |                 |
| `sequence`    | `Sequence.seqString()`                                                       | `RichSequence.seqString()`             |                 |
| `length`      | `Long.valueOf(Sequence.length())`                                            | `Long.valueOf(RichSequence.length())`  |                 |
| `attributes`  |                                                                              |                                        |                 |


### Slice

| bdg-formats   | biojava 1.x                                                                  | biojava 1.x biojavax                   | biojava 4.x/5.x |
|---------------|------------------------------------------------------------------------------|----------------------------------------|-----------------|
| `name`        | `Sequence.getName()`                                                         | `RichSequence.getName()`               |  tbd            |
| `description` | `Sequence.getAnnotation().getProperty(FastaFormat.PROPERTY_DESCRIPTIONLINE)` | `RichSequence.getDescription()`        |                 |
| `alphabet`    | `Sequence.getAlphabet()`                                                     | `RichSequence.getAlphabet()`           |                 |
| `sequence`    | `Sequence.seqString()`                                                       | `RichSequence.seqString()`             |                 |
| `start`       | view start position, in 0-based coordinate system, closed-open intervals     |                                        |                 |
| `end`         | view end position, in 0-based coordinate system, closed-open intervals       |                                        |                 |
| `strand`      | defaults to `Strand.Independent`                                             | defaults to `Strand.Independent`       |                 |
| `attributes`  |                                                                              |                                        |                 |


### Read

| bdg-formats           | biojava 1.x                                          | biojava 4.x/5.x                                      |
|-----------------------|------------------------------------------------------|------------------------------------------------------|
| `name`                | `Fastq.getDescription()`, first token split by space | `Fastq.getDescription()`, first token split by space |
| `description`         | `Fastq.getDescription()`                             | `Fastq.getDescription()`                             |
| `alphabet`            | `Alphabet.DNA`                                       | `Alphabet.DNA`                                       |
| `sequence`            | `Fastq.getSequence()`                                | `Fastq.getSequence()`                                |
| `length`              | `Long.valueOf(Fastq.getSequence().length())`         | `Long.valueOf(Fastq.getSequence().length())`         |
| `qualityScores`       | `Fastq.getQuality()`                                 | `Fastq.getQuality()`                                 |
| `qualityScoreVariant` | `Fastq.getVariant()`                                 | `Fastq.getVariant()`                                 |
| `attributes`          |                                                      |                                                      |

