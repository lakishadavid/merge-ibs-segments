# merge-ibs-segments

Joins broken IBD segments in the output of Brian L. Browning's
[Refined IBD](https://faculty.washington.edu/browning/refined-ibd.html) while
keeping haplotypes apart: two segments of a pair of samples are joined only if
they are on the same haplotype of each sample. Each output segment also gets
its number of markers and how many haplotypes in the cohort carry it, and a
report lists every join and every segment end that needs a closer look.

It began as an exact Python translation of Browning's `merge-ibd-segments`
(version 17Jan20.102), kept as the tag `baseline-faithful` and checked against
the original program (see "Verification"). The current version replaces that
program's merging rule, adds three output columns, a report and a summary,
checks that every segment end is a VCF marker, every sample is in the VCF
file and the VCF records are in order, and writes the usage message with
argument errors to standard error. VCF and map reading, genetic-map
calculations, number formatting and the remaining messages are unchanged
from the baseline.

## What it does

The segments of each pair of samples are grouped by the haplotype of sample 1
and the haplotype of sample 2 (each 1 or 2, independently). Within each group
they are sorted and walked in order. Two consecutive segments of a group on
the same chromosome are joined if

- they share an end marker;
- no marker strictly between them is discordant, whatever the length of the
  gap; or
- at most `[discord]` markers strictly between them are discordant
  (`[discord]` is 0 or 1; any number of other markers may lie between them),
  and the gap is at most `[gap]` cM.

A marker is discordant when sample 1's allele on sample 1's haplotype differs
from sample 2's allele on sample 2's haplotype; the samples' other haplotypes
are not considered. The gap is measured with the genetic map from the first
segment's end marker to the second segment's start marker. Joins form chains,
so more than two segments can become one.

Segments of one group that overlap, including nested and repeated segments
and a one-marker segment at the other segment's end marker, stop the program
with an error. In Refined IBD each segment is a run of
identical alleles on the two haplotypes, extended until the next marker
differs, so two segments of one haplotype pair are not expected to overlap in
its output (reading its source, they can only after a stretch of about 25 cM
without markers); an overlap means the input is not what this program
assumes. Any other pair of segments (two or more discordant markers between
them, a discordant marker and a gap longer than `[gap]`, or one discordant
marker when `[discord]` is 0) is left as it is and not reported. A segment of
another haplotype pair lying between them does not prevent a join.

Each chain is written as one segment, with the group's haplotype numbers, the
largest LOD score in the chain (not recalculated for the joined segment), and
its length recomputed from the map. Beyond either end of the map, positions
are extrapolated along the line through the end point and the nearest map
point at least 5 cM inside it (the other end of the map if the map spans less
than 5 cM).

Every segment's start and end must be the position of exactly one VCF record,
both samples must be in the VCF file, each chromosome's VCF records must
form one block in position order (the markers of a segment and of a gap are
taken by record), and allele indices must be at most 255 (the program holds
alleles as bytes; biallelic markers use 0 and 1); otherwise the program stops
with an error naming the segment, sample or record.

### How the rule differs from Browning's merge-ibd-segments

| | Browning's program | This program |
|---|---|---|
| Segments that may be joined | any segments of the pair | only segments with the same haplotype of each sample |
| Overlapping or nested segments | joined without a test | the program stops with an error |
| Segments sharing an end marker | joined without a test | joined and flagged |
| Gap with no discordant marker | joined if at most `[gap]` cM | joined whatever its length, and flagged |
| Gap with discordant markers | joined if at most `[gap]` cM and at most `[discord]` discordant markers, counting the two end markers | joined if at most `[gap]` cM and at most `[discord]` (0 or 1) discordant markers strictly between the segments |
| Discordant marker | no allele of sample 1 equals either allele of sample 2 (the genotypes share no allele) | sample 1's allele on its haplotype differs from sample 2's allele on its haplotype |
| Haplotype 0 in the input | accepted | rejected |
| Haplotypes of a joined segment | written as 0 | kept |
| Output | 9 columns | 12 columns, a report and a summary |

The gap length calculation, the largest score and the length from the map are
Browning's.

## Use

No installation and no packages beyond the Python standard library are
needed; the code has been run with Python 3.12 and 3.13. Refined IBD writes
its `.ibd` output gzip-compressed, and standard input must be uncompressed.
From the repository root:

```
gunzip -c out.ibd.gz | python -m merge_ibs_segments phased.vcf.gz plink.map 0.6 1 report.tsv > merged.ibd
```

| Argument | Meaning |
|---|---|
| standard input | uncompressed Refined IBD `.ibd` segments: 8 or 9 white-space separated fields (sample 1, haplotype 1, sample 2, haplotype 2, chromosome, start, end, LOD score, optional length in cM) |
| `[vcf]` | the phased VCF file given to Refined IBD (plain text, or gzip/BGZF when the file name ends in `.gz`) |
| `[map]` | PLINK genetic map with cM positions (plain text, or gzip/BGZF when the file name ends in `.gz`) |
| `[gap]` | largest gap, in cM, to join across a discordant marker; also the length above which a joined gap with no discordant marker is reported as `joined_no_discordant_long_gap` |
| `[discord]` | largest number of discordant markers in a joined gap: 0 or 1 |
| `[report]` | optional: file for the report. Checked before any input is read: it must not be empty, a directory, the `[vcf]` or `[map]` file, or an existing file that is not writable, and its directory must exist. A symbolic link is followed. |

With no arguments the program prints its usage message to standard output.
Any argument error ends it with exit status 1, the usage message and the
error on standard error, and nothing on standard output.

The Refined IBD web page describes removing gaps that are shorter than 0.6 cM
and have at most one discordant homozygote, which corresponds approximately to
`0.6 1` in Browning's program. Those values were given for his genotype-based
discordance test; they have not been calibrated for this program's
haplotype-based test.

### Output

One line per chain, 12 tab-separated fields: Browning's nine (sample 1,
haplotype 1, sample 2, haplotype 2, chromosome, start, end, LOD score with at
most 2 decimals, length in cM with at most 3 decimals), then

10. **markers**: the number of VCF records from the segment's start to its
    end, both included;
11. **carriers**: the number of haplotypes in the VCF file that carry the
    segment's alleles, counting the pair's own two haplotypes. A haplotype
    carries them if it has the pair's allele at every marker of the segment
    where the pair's two haplotypes agree; a marker where they differ (the
    discordant marker of a joined gap) is skipped;
12. **frequency**: carriers divided by all haplotypes in the VCF file (twice
    the number of samples), rounded half to even to 6 significant digits and
    written without an exponent.

The carriers and the frequency depend on who is in the VCF file: they change
as samples are added or removed. The program holds one byte per haplotype per
VCF record in memory to count carriers.

### Report

If `[report]` is given, the program writes a tab-separated file with a header
line and one line per event, sorted by chromosome and position:

| Column | Meaning |
|---|---|
| `event` | see below |
| `flagged` | `yes` for every event except `joined_one_discordant` |
| `sample1` … `chromosome` | the segment's samples, haplotypes and chromosome |
| `start1`, `end1` | the (first) segment |
| `start2`, `end2` | the second segment of a join (`.` otherwise) |
| `snps_between` | markers strictly between the two segments of a join |
| `discordant_between` | discordant markers among them |
| `gap_cm` | the gap in cM, from the first segment's end marker to the second segment's start marker |

| Event | Meaning |
|---|---|
| `joined_one_discordant` | joined across a gap of at most `[gap]` cM with one discordant marker |
| `joined_no_discordant` | joined across a gap of at most `[gap]` cM with no discordant marker |
| `joined_no_discordant_long_gap` | joined across a gap longer than `[gap]` cM with no discordant marker |
| `joined_shared_end` | joined at a shared end marker (the gap fields are 0) |
| `start_matches_previous_snp` | the marker before the segment's start carries the same allele on both haplotypes |
| `end_matches_next_snp` | the marker after the segment's end carries the same allele on both haplotypes |

In Refined IBD output, a gap with no discordant marker, a shared end marker,
and a segment end whose neighbouring marker is not discordant are not
expected, except where a segment was cut at the edge of one of Refined IBD's
analysis windows; they are flagged so that they can be examined. The report
is written only when the program succeeds: it goes first to a temporary file
beside the named one, which replaces it after the segments have been written.
A failed run leaves no new report and leaves an earlier report at that path
unchanged, with two exceptions. If the report cannot be moved into place at
the very end (for example because another program holds the file open), the
segments have already been written and the program stops with an error, so
check the exit status. A run that is killed can leave its temporary file
(`.merge-ibs-report-*`) beside the report. If standard output fails (for
example a pipe whose reader has gone), the program stops with
`ERROR: cannot write the merged segments` and exit status 1.

A one-line summary of the counts goes to standard error.

## Verification

### Tests of the rule (`tests/rule`)

`tests/rule/test_rule.py` checks, on synthetic data:

- a join across one discordant marker, with `[discord]` 1 and not with 0;
  a gap of exactly `[gap]` cM is joined and is not joined when `[gap]` is one
  double smaller;
- two discordant markers, or one with a gap longer than `[gap]`: not joined
  and not reported;
- gaps with no discordant marker, short and 12.5 cM long, and a shared end
  marker: joined and flagged, with either `[discord]`;
- overlapping, nested and repeated segments of one haplotype pair, including
  one-marker segments at the other segment's end marker: the program stops
  with the expected message and writes neither output nor report; when the
  segments cannot be written, no new report or temporary file remains and an
  earlier report is unchanged;
- segments of other haplotype pairs neither join nor block a join; a change
  of haplotype of either sample, or both, prevents a join;
- only the respective haplotypes count as discordant, including a marker
  where the genotypes share an allele but the two haplotypes differ;
- chains of five segments, shuffled input, the output order, separate
  chromosomes;
- discordant markers are counted strictly between the segments, not at
  their end markers;
- the marker count (including two records at one position inside a
  segment), the carrier count (skipping a discordant marker inside a joined
  segment) and the frequency format;
- segment ends that are not a VCF marker or have two VCF records, sample 1
  or sample 2 not in the VCF file, VCF records out of position order or a
  chromosome split into two blocks, an allele index above 255, haplotype 0,
  `[discord]` values other than 0 and 1, and the number of arguments;
- `[report]` paths that are empty, a directory, in a missing directory,
  not writable, or the `[vcf]` or `[map]` file: stopped before anything is
  written;
- a sample name outside the Basic Multilingual Plane in the report;
- lengths extrapolated before and after the map and number rounding;
- the command line: output, report file and summary, and an output pipe
  closed by its reader (an error message and status 1);
- 30 seeded random data sets: haplotypes built as mosaics of three founder
  haplotypes with rare changes, every maximal identical run of 0.8 cM or more
  of each haplotype pair taken as a segment, some of them cut the ways Refined
  IBD can (an end cut short, a middle part missing, a split at a shared
  marker). The output, the report and the summary are compared with a
  separate implementation of the rule written for the test, and the run is
  repeated with the input in another order. Every kind of join and flag
  occurs.

All 35 tests pass with Python 3.12 and 3.13.

### Exact-conversion test of the baseline (`tests/parity`)

This test applies to the tag `baseline-faithful`, the exact translation of
Browning's program; the current version differs from the jar by its merging
rule, its checks, its output columns and its usage text.
`tests/parity/run_parity.py` runs the original jar and the baseline on the
same inputs and compares them. Successful runs must give the same output lines
(order aside). Failing runs must give the same exit status and the same
standard output and error after the normalisations listed in `run_parity.py`:
Java stack-trace frame lines, exceptions re-thrown from Java's parallel VCF
parsing, Java object identity hashes, the usage text's command name, segment
lines Java had already written, and blank lines. The inputs are 322 synthetic
cases generated by `tests/parity/cases.py`:

- every branch of the merging logic, threshold boundaries, shared endpoints,
  nested and duplicate segments, and changes of input order (29 cases);
- genetic-map interpolation and extrapolation, including repeated cM values
  (10);
- number parsing and formatting: 52,306 lengths and 48,813 scores printed by
  both programs, 55 further score strings and 21 cM values, and 124 unusual
  number strings given as scores (103) or haplotypes (21) (126);
- input formats: gzip, BGZF, line endings, blank lines, FORMAT fields,
  multi-allelic markers, Unicode sample names (19);
- command-line arguments (32) and single input errors (86);
- seeded random data sets (20).

With Java 25.0.4.1 on Linux (WSL2) and Python 3.12.14, 320 cases give
identical results and 2 show the known differences listed below. The full
record, with checksums of the jar and of every source file, is in
[`tests/parity/parity_record.md`](tests/parity/parity_record.md).

### Other differences from Browning's program

By design, in the baseline and the current version:

- Output lines are sorted: by chromosome in order of first appearance in the
  map, then start, end, score, samples (in VCF order; in the baseline,
  samples absent from the VCF follow, in order of first appearance in the
  input) and haplotypes. The Java program writes sample pairs in an
  arbitrary order.
- The whole output is written at the end. If the program stops with an error,
  no segments are written; the Java program may already have written some.
- Output lines end with `\n`. The Java program uses the platform's line
  separator, so the two agree where that is `\n` (Linux, where the
  exact-conversion test ran, and macOS), not on Windows.
- An error that Java re-throws from its parallel VCF parsing is reported once,
  as the original exception.
- Java object identity hashes (for example `@76ed5528`) are not printed.
- Java stack-trace frame lines (`at ...`) are not reproduced, and the usage
  message names this program's command. In the current version the usage
  message goes to standard error with any argument error (exit status 1);
  the Java program writes it to standard output, with exit status 0 for a
  wrong number of arguments.

Known and accepted, found by the exact-conversion test of the baseline:

- Some numbers of 2**53 (about 9.0 × 10^15) or more are printed with other
  digits. Both forms convert back to the same number. In the test, 169 of the
  1,794 such values printed differently, the smallest about 1.84 × 10^16. For
  example, the Java program prints 2**60 as `1152921504606846980` and this
  program as `1152921504606847000`. Real LOD scores and cM lengths are far
  smaller.
- A gap test for a sample that is not in the VCF file: the Java program stops
  at a bounds check or reads one element past its genotype array, then stops
  with an index exception; in the test both programs stopped with
  `IndexOutOfBoundsException`, naming index 9 (Java) and 8 (the baseline).
  The current version stops earlier, for any segment naming such a sample
  (as sample 1 or sample 2), with `ERROR: sample not in VCF file`.

Not compared with the jar:

- The messages for a damaged BGZF file.
- VCF files with more than one error: the Java program reads and parses a VCF
  file on several threads, so which error it reports can vary. This program
  reports the first in file order.
- Java versions other than 25.0.4.1; operating systems other than Linux;
  Python versions other than 3.12.

## Running the tests

The rule tests need only Python. From the repository root:

```
python -m unittest discover -s tests/rule
```

The exact-conversion test applies to the baseline, so check out its tag first
(`git checkout baseline-faithful`). It needs Java and the original program,
which is not distributed here. Download them from Browning's site:

- `https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar`
  (SHA-256 `f5a8e8d094e99fa8226e8489e2b55e4a2bc89b325c49de19501395b15769dc80`)
- source, for reference: `https://faculty.washington.edu/browning/refined-ibd/refined-ibd.17Jan20.102.zip`
  (SHA-256 `0b0abf48528ec53d3fec7f21a9742c6e5960857c5a225b0e49346179bc45e6ea`)

Download the jar into the repository root, then run from there:

```
python tests/parity/run_parity.py --jar merge-ibd-segments.17Jan20.102.jar --java-opt=-Xmx768m --work /tmp/mis_parity --report tests/parity/parity_record.md
```

The `--work` directory is emptied first. The script exits with status 0 when
no case shows an unexpected difference.

## Licence

Every file in this repository except LICENSE, including documentation and the
test code, is released under the GNU General Public License, version 3 or (at
your option) any later version (SPDX: GPL-3.0-or-later). LICENSE is the text
of that licence. The software comes with no warranty (LICENSE, sections
15-16).

This is a modified Python translation of merge-ibd-segments from Refined-IBD
17Jan20.102, Copyright (C) 2014-2016 Brian L. Browning, GPL-3.0-or-later. Each
translated module keeps his original notice and names the Java files it was
translated from. The BGZF reading in `merge_ibs_segments/input_it.py` is
translated from code Copyright (c) 2009 The Broad Institute under the MIT
licence, whose notice is kept in that file.

Copyright (C) 2026 LaKisha T. David.

## Not affiliated

This is an independent project. It is not written, maintained or endorsed by
Brian L. Browning or the University of Washington. Please report problems in
this repository's GitHub issues.

## Citation

If you use this software in published work, please cite both:

- this repository, with the version (tag or commit) you used, for example:
  David LT. merge-ibs-segments, commit `<hash>` [software].
  https://github.com/lakishadavid/merge-ibs-segments
- the Refined IBD paper, for Refined IBD and the original merge-ibd-segments
  program: Browning BL, Browning SR (2013). Improving the accuracy and
  efficiency of identity-by-descent detection in population data. *Genetics*
  194(2):459-471. doi:10.1534/genetics.113.150029

This is a request, not a condition of the licence.
