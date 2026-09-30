# merge-ibs-segments

Joins broken IBD segments in the output of Brian L. Browning's
[Refined IBD](https://faculty.washington.edu/browning/refined-ibd.html) while
keeping haplotypes apart: two segments of a pair of samples are joined only if
they are on the same haplotype of each sample.

It began as an exact Python translation of Browning's `merge-ibd-segments`
(version 17Jan20.102), kept as the tag `baseline-faithful` and checked against
the original program (see "Verification"). The current version replaces that
program's merging rule: haplotype 0 is rejected and the usage message
describes the new rule. Everything else (VCF and map reading, genetic-map
calculations, number formatting, other messages) is unchanged from the
baseline, except the index named when a sample absent from the VCF file takes
part in a gap test (see "Known and accepted").

## What it does

The segments of each pair of samples are grouped by the haplotype of sample 1
and the haplotype of sample 2. Within each group they are sorted and walked in
order. A segment joins the current chain if it is on the same chromosome and
either

- starts before the chain's running end (an overlapping or nested segment of
  the same haplotype pair joins without any test), or
- starts at or after the running end, and the gap passes both tests:
  - its length, measured with the genetic map from the marker at the chain's
    running end to the marker at the segment's start, is at most `[gap]` cM
    (a segment that starts at the running end has a gap of length zero and
    is always tested, even a one-marker segment inside the chain, which a
    minimum segment length rules out in practice);
  - it contains at most `[discord]` discordant markers, counting both of
    those markers. A marker is discordant when sample 1's allele on sample 1's
    haplotype differs from sample 2's allele on sample 2's haplotype. The
    samples' other haplotypes are not considered.

A segment of another haplotype pair lying in the gap does not prevent a join,
and segments of different haplotype pairs are never joined, even where they
overlap. Both gap positions must be markers in the VCF file. Otherwise the
program stops with "position missing from VCF file", which, as in Browning's
program, gives the chromosome as a 0-based index in map order and, for the
next segment, its end position.

Each chain is written as one segment, with the group's haplotype numbers, the
largest score in the chain, and its length recomputed from the map. Beyond
either end of the map, positions are extrapolated along the line through the
end point and the nearest map point at least 5 cM inside it (the other end of
the map if the map spans less than 5 cM).

### How the rule differs from Browning's merge-ibd-segments

| | Browning's program | This program |
|---|---|---|
| Segments that may be joined | any segments of the pair | only segments with the same haplotype of each sample |
| A segment that starts at the running end | joined without a test | tested as a gap of length zero (so that position must be a VCF marker) |
| Discordant marker | no allele of sample 1 equals either allele of sample 2 (the genotypes share no allele) | sample 1's allele on its haplotype differs from sample 2's allele on its haplotype |
| Haplotype 0 in the input | accepted | rejected |
| Haplotypes of a joined segment | written as 0 | kept |

The gap length test and its limits, the running-end walk, the largest score
and the length from the map are Browning's.

## Use

No installation and no packages beyond the Python standard library are
needed; the code has been run with Python 3.12 and 3.13. Refined IBD writes
its `.ibd` output gzip-compressed, and standard input must be uncompressed.
From the repository root:

```
gunzip -c out.ibd.gz | python -m merge_ibs_segments phased.vcf.gz plink.map 0.6 1 > merged.ibd
```

| Argument | Meaning |
|---|---|
| standard input | uncompressed Refined IBD `.ibd` segments: 8 or 9 white-space separated fields (sample 1, haplotype 1, sample 2, haplotype 2, chromosome, start, end, LOD score, optional length in cM) |
| `[vcf]` | the phased VCF file given to Refined IBD (plain text, or gzip/BGZF when the file name ends in `.gz`) |
| `[map]` | PLINK genetic map with cM positions (plain text, or gzip/BGZF when the file name ends in `.gz`) |
| `[gap]` | largest gap to join, in cM |
| `[discord]` | largest number of discordant markers in a gap |

The Refined IBD web page describes removing gaps that are shorter than 0.6 cM
and have at most one discordant homozygote, which corresponds approximately to
`0.6 1` in Browning's program. Those values were given for his genotype-based
discordance test; they have not been calibrated for this program's
haplotype-based test.

Output: one line per chain, 9 tab-separated fields in the input's order, with
the score written with at most 2 decimals and the length with at most 3.

## Verification

### Tests of the rule (`tests/rule`)

`tests/rule/test_rule.py` checks the rule on synthetic data:

- a pair keeps its haplotype numbers across a qualifying gap and is joined
  (for A1–B1 and A1–B2), and the joined segment keeps those numbers;
- a change of haplotype in sample 1, in sample 2, or in both, prevents a join;
- a segment of another haplotype pair spanning the gap does not prevent a
  join, and overlapping segments of different haplotype pairs are not joined;
- a gap of exactly `[gap]` cM is joined and is not joined when `[gap]` is one
  double smaller; exactly `[discord]` discordant markers are allowed and one
  more is not;
- only the respective haplotypes count as discordant, for sample 1 and for
  sample 2;
- a shared endpoint is tested as a gap of length zero, and must be a VCF
  marker; a segment starting at the running end is always tested;
- overlapping and nested segments of the same haplotype pair are joined, and
  the gap after a nested segment is measured from the chain's running end;
- chains of several segments, shuffled input, the output order, separate
  chromosomes, haplotype 0 (rejected), lengths extrapolated before and after
  the map, and number rounding;
- 30 seeded random data sets, compared with a separate implementation of the
  rule written for the test, and rerun in a second input order.

All 23 tests pass with Python 3.12 and 3.13.

### Exact-conversion test of the baseline (`tests/parity`)

This test applies to the tag `baseline-faithful`, the exact translation of
Browning's program; the current version differs from the jar by its merging
rule and usage text. `tests/parity/run_parity.py` runs the original jar and
the baseline on the same inputs and compares them. Successful runs must give
the same output lines (order aside). Failing runs must give the same exit
status and the same
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

These apply to the baseline and to the current version.

By design:

- Output lines are sorted: by chromosome in order of first appearance in the
  map, then start, end, score, samples (in VCF order; samples absent from the
  VCF after them, in order of first appearance in the input) and haplotypes.
  The Java program writes sample pairs in an arbitrary order.
- The whole output is written at the end. If the program stops with an error,
  no segments are written; the Java program may already have written some.
- Output lines end with `\n`. The Java program uses the platform's line
  separator, so the two agree where that is `\n` (Linux, where the
  exact-conversion test ran, and macOS), not on Windows.
- An error that Java re-throws from its parallel VCF parsing is reported once,
  as the original exception.
- Java object identity hashes (for example `@76ed5528`) are not printed.
- Java stack-trace frame lines (`at ...`) are not reproduced, and the usage
  message names this program's command.

Known and accepted:

- Some numbers of 2**53 (about 9.0 × 10^15) or more are printed with other
  digits. Both forms convert back to the same number. In the test, 169 of the
  1,794 such values printed differently, the smallest about 1.84 × 10^16. For
  example, the Java program prints 2**60 as `1152921504606846980` and this
  program as `1152921504606847000`. Real LOD scores and cM lengths are far
  smaller.
- If a sample that is not in the VCF file takes part in a gap test: depending
  on how the Java program stored the genotypes, it either stops at a bounds
  check or reads one element past its genotype array, and then stops with an
  index exception. In the test both programs stopped with
  `IndexOutOfBoundsException`, naming index 9 (Java) and 8 (the baseline).
  The current version names the index of the haplotype it reads
  (2 × sample + haplotype − 1).

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
