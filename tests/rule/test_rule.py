# Copyright (C) 2026 LaKisha T. David
#
# This file is part of merge-ibs-segments, a modified Python translation of
# merge-ibd-segments from Refined-IBD version 17Jan20.102 (Copyright (C)
# 2014-2016 Brian L. Browning).
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Tests of the merging rule, the added output columns and the report.

Run from the repository root: ``python -m unittest discover -s tests/rule``.

The hand-built cases each state their expected output, worked out from the
rule.  Unless a case says otherwise, their markers are 1 kb apart and
0.125 cM apart (an exact binary fraction, so gap lengths are exact), and
every genotype is 0|0 unless a case sets another allele.  The seeded random
cases build haplotypes as mosaics of a few founder haplotypes, take the
maximal identical runs of each haplotype pair as segments, cut some of them
the way Refined IBD can (an end cut short, a middle part missing, a split at
a shared marker), and compare the program with ``oracle``, a separate
implementation of the rule written for this test.  The oracle uses the
program's genetic map, its score parser and its formatting of scores and
lengths, which the exact-conversion test of the baseline checked against
Browning's jar; it formats the frequency on its own.

All data are synthetic and built here.
"""

from __future__ import annotations

import io
import math
import os
import random
import subprocess
import sys
import tempfile
import unittest
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from merge_ibs_segments import merge  # noqa: E402
from merge_ibs_segments.genetic_map import PlinkGenMap  # noqa: E402
from merge_ibs_segments.ids import Indexer  # noqa: E402
from merge_ibs_segments.java_compat import (  # noqa: E402
    JavaException, JavaTermination, java_decimal_format, java_encode,
    java_parse_float)

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
SAMPLES = ["A", "B", "C", "D"]                       # 8 haplotypes
POSITIONS = [1000 * i for i in range(1, 201)]       # 1000, 2000, ..., 200000
STEP_CM = 0.125                                     # cM between markers
MAP = [("1", repr(STEP_CM * i), p) for i, p in enumerate(POSITIONS)]
REPORT_HEADER = "\t".join(merge.REPORT_COLUMNS)


class Genotypes:
    """Phased genotypes, 0|0 unless set."""

    def __init__(self):
        self.calls: Dict[int, Dict[str, List[int]]] = {}

    def set(self, pos: int, sample: str, hap: int, allele: int = 1):
        self.calls.setdefault(pos, {}).setdefault(sample, [0, 0])[hap - 1] = allele
        return self

    def gt(self, pos: int, sample: str) -> str:
        return "%d|%d" % tuple(self.calls.get(pos, {}).get(sample, [0, 0]))


def write_inputs(directory: str, geno: Genotypes, chrom: str = "1",
                 positions: Sequence[int] = POSITIONS,
                 samples: Sequence[str] = SAMPLES,
                 map_entries: Sequence[Tuple[str, str, int]] = MAP,
                 extra_records: Sequence[str] = (),
                 after: Optional[Dict[int, str]] = None) -> Tuple[str, str]:
    """Writes a phased VCF and a PLINK map.  ``after`` puts an extra record
    straight after the record at a position; ``extra_records`` go at the
    end."""
    lines = ["##fileformat=VCFv4.2", HEADER + "".join("\t" + s for s in samples)]
    for pos in positions:
        lines.append("\t".join([chrom, str(pos), ".", "A", "G", ".", "PASS",
                                ".", "GT"] + [geno.gt(pos, s) for s in samples]))
        if after and pos in after:
            lines.append(after[pos])
    lines.extend(extra_records)
    vcf = os.path.join(directory, "in.vcf")
    with open(vcf, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    mp = os.path.join(directory, "in.map")
    with open(mp, "w", encoding="utf-8", newline="\n") as f:
        for c, cm, bp in map_entries:
            f.write("%s\t.\t%s\t%d\n" % (c, cm, bp))
    return vcf, mp


def run(rows: Sequence[str], geno: Optional[Genotypes] = None,
        gap: str = "0.6", discord: str = "1", **kw) -> merge.Result:
    """Runs the program on segment lines with a report file, writes its
    outputs, and checks that standard output and the file hold what the run
    returned."""
    with tempfile.TemporaryDirectory() as d:
        vcf, mp = write_inputs(d, geno or Genotypes(), **kw)
        report = os.path.join(d, "report.tsv")
        stdin = io.BytesIO(("\n".join(rows) + "\n").encode("utf-8"))
        result = merge.run([vcf, mp, gap, discord, report], stdin)
        stdout = io.BytesIO()
        merge.write_outputs(result, stdout)
        assert stdout.getvalue() == java_encode(result.output)
        with open(report, "rb") as f:
            assert f.read() == java_encode(result.report)
        assert sorted(os.listdir(d)) == ["in.map", "in.vcf", "report.tsv"]
        return result


def lines(rows, geno=None, **kw) -> List[str]:
    return run(rows, geno, **kw).output.splitlines()


def report_lines(rows, geno=None, **kw) -> List[str]:
    return run(rows, geno, **kw).report.splitlines()


def seg(s1, h1, s2, h2, start, end, score="3", chrom="1") -> str:
    return "\t".join(str(x) for x in (s1, h1, s2, h2, chrom, start, end, score))


def out(s1, h1, s2, h2, start, end, score, length, snps, carriers, freq,
        chrom="1") -> str:
    return "\t".join(str(x) for x in (s1, h1, s2, h2, chrom, start, end, score,
                                      length, snps, carriers, freq))


def rep(event, s1, h1, s2, h2, start1, end1, start2=".", end2=".", snps=".",
        disc=".", gap=".", chrom="1") -> str:
    flagged = "yes" if merge.FLAGGED[event] else "no"
    return "\t".join(str(x) for x in (event, flagged, s1, h1, s2, h2, chrom,
                                      start1, end1, start2, end2, snps, disc,
                                      gap))


def termination(rows, geno=None, gap="0.6", discord="1",
                **kw) -> JavaTermination:
    with tempfile.TemporaryDirectory() as d:
        vcf, mp = write_inputs(d, geno or Genotypes(), **kw)
        stdin = io.BytesIO(("\n".join(rows) + "\n").encode("utf-8"))
        try:
            merge.run([vcf, mp, gap, discord], stdin)
        except JavaTermination as t:
            return t
    raise AssertionError("the program did not stop")


ONE = merge.JOINED_ONE_DISCORDANT
NONE_SHORT = merge.JOINED_NO_DISCORDANT
NONE_LONG = merge.JOINED_NO_DISCORDANT_LONG_GAP
SHARED = merge.JOINED_SHARED_END
START = merge.START_MATCHES_PREVIOUS
END = merge.END_MATCHES_NEXT


class Joining(unittest.TestCase):

    def test_one_discordant_marker_joins(self):
        # A's haplotype 1 differs from B's haplotype 1 at 31000 (and at
        # 61000, just past the second segment).  C's haplotype 2 differs at
        # 20000, so it does not carry the segment; D's haplotype 1 differs
        # only at the skipped marker 31000, so it does.
        geno = (Genotypes().set(31000, "A", 1).set(61000, "A", 1)
                .set(20000, "C", 2).set(31000, "D", 1))
        rows = [seg("A", 1, "B", 1, 1000, 30000, "3"),
                seg("A", 1, "B", 1, 32000, 60000, "5")]
        got = run(rows, geno, discord="1")
        self.assertEqual(got.output.splitlines(),
                         [out("A", 1, "B", 1, 1000, 60000, "5", "7.375", 60, 7,
                              "0.875")])
        self.assertEqual(got.report.splitlines(),
                         [REPORT_HEADER,
                          rep(ONE, "A", 1, "B", 1, 1000, 30000, 32000, 60000,
                              1, 1, "0.25")])
        # With [discord] 0 the gap is not joined, and not reported.
        got = run(rows, geno, discord="0")
        self.assertEqual(got.output.splitlines(),
                         [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 7,
                              "0.875"),
                          out("A", 1, "B", 1, 32000, 60000, "5", "3.5", 29, 8,
                              "1")])
        self.assertEqual(got.report.splitlines(), [REPORT_HEADER])

    def test_gap_length_limit_for_one_discordant_marker(self):
        # One discordant marker (31000) among three in a 0.5 cM gap; the
        # identical markers next to the segments flag their ends.
        geno = Genotypes().set(31000, "A", 1).set(61000, "A", 1)
        rows = [seg("A", 1, "B", 1, 1000, 29000), seg("A", 1, "B", 1, 33000, 60000)]
        flags = [rep(END, "A", 1, "B", 1, 1000, 29000),
                 rep(START, "A", 1, "B", 1, 33000, 60000)]
        got = run(rows, geno, gap="0.5")
        self.assertEqual(got.output.splitlines(),
                         [out("A", 1, "B", 1, 1000, 60000, "3", "7.375", 60, 8, "1")])
        self.assertEqual(got.report.splitlines(),
                         [REPORT_HEADER, flags[0],
                          rep(ONE, "A", 1, "B", 1, 1000, 29000, 33000, 60000,
                              3, 1, "0.5"),
                          flags[1]])
        # The largest double below 0.5.
        got = run(rows, geno, gap=repr(math.nextafter(0.5, 0.0)))
        self.assertEqual(got.output.splitlines(),
                         [out("A", 1, "B", 1, 1000, 29000, "3", "3.5", 29, 8, "1"),
                          out("A", 1, "B", 1, 33000, 60000, "3", "3.375", 28, 8,
                              "1")])
        self.assertEqual(got.report.splitlines(), [REPORT_HEADER] + flags)

    def test_two_discordant_markers_are_ignored(self):
        geno = (Genotypes().set(31000, "A", 1).set(32000, "A", 1)
                .set(61000, "A", 1))
        rows = [seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 33000, 60000)]
        got = run(rows, geno)
        self.assertEqual(got.output.splitlines(),
                         [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 8, "1"),
                          out("A", 1, "B", 1, 33000, 60000, "3", "3.375", 28, 8,
                              "1")])
        self.assertEqual(got.report.splitlines(), [REPORT_HEADER])

    def test_discordance_with_a_long_gap_is_ignored(self):
        # One discordant marker in a 0.875 cM gap: not joined, and the join
        # is not reported (only the second segment's start is flagged).
        geno = Genotypes().set(31000, "A", 1).set(61000, "A", 1)
        rows = [seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 37000, 60000)]
        got = run(rows, geno)
        self.assertEqual(got.output.splitlines(),
                         [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 8, "1"),
                          out("A", 1, "B", 1, 37000, 60000, "3", "2.875", 24, 8,
                              "1")])
        self.assertEqual(got.report.splitlines(),
                         [REPORT_HEADER, rep(START, "A", 1, "B", 1, 37000, 60000)])

    def test_no_discordant_marker_in_a_short_gap_joins_and_flags(self):
        geno = Genotypes().set(61000, "A", 1)
        rows = [seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 33000, 60000)]
        for discord in ("0", "1"):
            with self.subTest(discord=discord):
                got = run(rows, geno, discord=discord)
                self.assertEqual(got.output.splitlines(),
                                 [out("A", 1, "B", 1, 1000, 60000, "3", "7.375",
                                      60, 8, "1")])
                self.assertEqual(got.report.splitlines(),
                                 [REPORT_HEADER,
                                  rep(END, "A", 1, "B", 1, 1000, 30000),
                                  rep(NONE_SHORT, "A", 1, "B", 1, 1000, 30000,
                                      33000, 60000, 2, 0, "0.375"),
                                  rep(START, "A", 1, "B", 1, 33000, 60000)])

    def test_no_discordant_marker_in_a_long_gap_joins_and_flags(self):
        # Refined IBD's case C: the middle of a run is missing, 12.5 cM here.
        geno = Genotypes().set(161000, "A", 1)
        rows = [seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 130000, 160000)]
        for discord in ("0", "1"):
            with self.subTest(discord=discord):
                got = run(rows, geno, discord=discord)
                self.assertEqual(got.output.splitlines(),
                                 [out("A", 1, "B", 1, 1000, 160000, "3", "19.875",
                                      160, 8, "1")])
                self.assertEqual(got.report.splitlines(),
                                 [REPORT_HEADER,
                                  rep(END, "A", 1, "B", 1, 1000, 30000),
                                  rep(NONE_LONG, "A", 1, "B", 1, 1000, 30000,
                                      130000, 160000, 99, 0, "12.5"),
                                  rep(START, "A", 1, "B", 1, 130000, 160000)])

    def test_shared_end_marker_joins_and_flags(self):
        # The marker beyond each end lies inside the other segment, so both
        # ends are flagged too.
        geno = Genotypes().set(61000, "A", 1)
        rows = [seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 30000, 60000)]
        got = run(rows, geno)
        self.assertEqual(got.output.splitlines(),
                         [out("A", 1, "B", 1, 1000, 60000, "3", "7.375", 60, 8, "1")])
        self.assertEqual(got.report.splitlines(),
                         [REPORT_HEADER,
                          rep(END, "A", 1, "B", 1, 1000, 30000),
                          rep(SHARED, "A", 1, "B", 1, 1000, 30000, 30000, 60000,
                              0, 0, 0),
                          rep(START, "A", 1, "B", 1, 30000, 60000)])

    def test_overlap_stops_the_program(self):
        cases = {
            "overlap": [seg("A", 1, "B", 1, 1000, 40000),
                        seg("A", 1, "B", 1, 30000, 60000)],
            "nested": [seg("A", 1, "B", 1, 1000, 60000),
                       seg("A", 1, "B", 1, 20000, 30000)],
            "duplicate": [seg("A", 1, "B", 1, 1000, 30000),
                          seg("A", 1, "B", 1, 1000, 30000)],
            # One-marker segments lie inside the other segment even where
            # they touch it only at its end marker.
            "one-marker duplicate": [seg("A", 1, "B", 1, 30000, 30000),
                                     seg("A", 1, "B", 1, 30000, 30000)],
            "one-marker at the end": [seg("A", 1, "B", 1, 1000, 30000),
                                      seg("A", 1, "B", 1, 30000, 30000)],
            "one-marker at the start": [seg("A", 1, "B", 1, 30000, 60000),
                                        seg("A", 1, "B", 1, 30000, 30000)],
        }
        names = {"overlap": ("1000-40000", "30000-60000"),
                 "nested": ("1000-60000", "20000-30000"),
                 "duplicate": ("1000-30000", "1000-30000"),
                 "one-marker duplicate": ("30000-30000", "30000-30000"),
                 "one-marker at the end": ("1000-30000", "30000-30000"),
                 "one-marker at the start": ("30000-30000", "30000-60000")}
        for case, rows in cases.items():
            with self.subTest(case=case):
                t = termination(rows)
                first, second = names[case]
                self.assertEqual(t.status, 1)
                self.assertEqual(t.stdout, "")
                self.assertEqual(
                    t.stderr,
                    "ERROR: segments of the same samples and haplotypes "
                    "overlap: A haplotype 1, B haplotype 1, chromosome 1, "
                    + first + " and A haplotype 1, B haplotype 1, chromosome "
                    "1, " + second + ".  Refined IBD is not expected to "
                    "report overlapping segments of one haplotype pair, so "
                    "this input is not what the program assumes.\n")

    def test_other_haplotype_groups_neither_join_nor_block(self):
        # Overlapping segments of different haplotype groups are fine.
        got = lines([seg("A", 1, "B", 1, 1000, 40000),
                     seg("A", 1, "B", 2, 30000, 60000)])
        self.assertEqual(got, [out("A", 1, "B", 1, 1000, 40000, "3", "4.875", 40, 8, "1"),
                               out("A", 1, "B", 2, 30000, 60000, "3", "3.75", 31, 8,
                                   "1")])
        # A segment of another group lying across the gap does not prevent
        # the join.
        got = lines([seg("A", 1, "B", 1, 1000, 30000), seg("A", 2, "B", 2, 20000, 45000),
                     seg("A", 1, "B", 1, 33000, 60000)])
        self.assertEqual(got, [out("A", 1, "B", 1, 1000, 60000, "3", "7.375", 60, 8, "1"),
                               out("A", 2, "B", 2, 20000, 45000, "3", "3.125", 26,
                                   8, "1")])
        # A change of haplotype of either sample, or both, prevents a join.
        for h1, h2 in ((2, 1), (1, 2), (2, 2)):
            with self.subTest(second=(h1, h2)):
                got = lines([seg("A", 1, "B", 1, 1000, 30000),
                             seg("A", h1, "B", h2, 30000, 60000)])
                self.assertEqual(len(got), 2)

    def test_only_the_respective_haplotypes_count(self):
        rows = [seg("A", 1, "B", 2, 1000, 30000), seg("A", 1, "B", 2, 33000, 60000)]
        joined = [out("A", 1, "B", 2, 1000, 60000, "3", "7.375", 60, 7, "0.875")]
        apart = [out("A", 1, "B", 2, 1000, 30000, "3", "3.625", 30, 8, "1"),
                 out("A", 1, "B", 2, 33000, 60000, "3", "3.375", 28, 8, "1")]
        # B's haplotype 2 differs from A's haplotype 1: discordant.
        self.assertEqual(lines(rows, Genotypes().set(31000, "B", 2), discord="0"),
                         apart)
        # B's haplotype 1 or A's haplotype 2 differs instead: not counted.
        # (That haplotype then does not carry the joined segment.)
        for sample, hap in (("B", 1), ("A", 2)):
            with self.subTest(other=(sample, hap)):
                self.assertEqual(
                    lines(rows, Genotypes().set(31000, sample, hap), discord="0"),
                    joined)
        # A 1|0 and B 0|1 share an allele, but A's haplotype 1 and B's
        # haplotype 1 differ: discordant (Browning's test would not count it).
        geno = Genotypes().set(31000, "A", 1).set(31000, "B", 2)
        self.assertEqual(
            lines([seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 33000, 60000)],
                  geno, discord="0"),
            [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 8, "1"),
             out("A", 1, "B", 1, 33000, 60000, "3", "3.375", 28, 8, "1")])

    def test_chain_of_five_and_shuffled_input(self):
        # Discordant markers at 21000, 41000, 61000 and 101000 split one run
        # into five segments, joined across one marker each; 131000 ends it.
        geno = Genotypes()
        for pos in (21000, 41000, 61000, 101000, 131000):
            geno.set(pos, "A", 2)
        rows = [seg("A", 2, "C", 2, 1000, 20000, "1"),
                seg("A", 2, "C", 2, 22000, 40000, "4"),
                seg("A", 2, "C", 2, 42000, 60000, "2"),
                seg("A", 2, "C", 2, 62000, 100000, "3"),
                seg("A", 2, "C", 2, 102000, 130000, "6")]
        # The four discordant markers inside are skipped when counting
        # carriers, so every haplotype carries the segment.
        expected = [out("A", 2, "C", 2, 1000, 130000, "6", "16.125", 130, 8,
                        "1")]
        self.assertEqual(lines(rows, geno), expected)
        for seed in range(5):
            shuffled = rows[:]
            random.Random(seed).shuffle(shuffled)
            self.assertEqual(lines(shuffled, geno), expected)

    def test_output_order(self):
        rows = [seg("B", 1, "C", 1, 50000, 90000), seg("A", 2, "B", 2, 50000, 90000),
                seg("A", 1, "C", 2, 1000, 30000), seg("A", 1, "B", 2, 1000, 30000)]
        self.assertEqual(lines(rows),
                         [out("A", 1, "B", 2, 1000, 30000, "3", "3.625", 30, 8, "1"),
                          out("A", 1, "C", 2, 1000, 30000, "3", "3.625", 30, 8, "1"),
                          out("A", 2, "B", 2, 50000, 90000, "3", "5", 41, 8, "1"),
                          out("B", 1, "C", 1, 50000, 90000, "3", "5", 41, 8, "1")])

    def test_different_chromosomes_not_joined(self):
        records = ["2\t%d\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0" % p
                   for p in POSITIONS]
        mp = MAP + [("2", repr(STEP_CM * i), p) for i, p in enumerate(POSITIONS)]
        got = lines([seg("A", 1, "B", 1, 1000, 30000),
                     seg("A", 1, "B", 1, 32000, 60000, chrom="2")],
                    extra_records=records, map_entries=mp)
        self.assertEqual(got, [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 8, "1"),
                               out("A", 1, "B", 1, 32000, 60000, "3", "3.5", 29, 8,
                                   "1", chrom="2")])


class Columns(unittest.TestCase):

    def test_carriers_and_frequency(self):
        # Only A's and B's haplotype 1 carry the segment's alleles.
        geno = (Genotypes().set(10000, "C", 1).set(10000, "C", 2)
                .set(12000, "A", 2).set(12000, "B", 2))
        row = [seg("A", 1, "B", 1, 1000, 30000)]
        self.assertEqual(lines(row, geno, samples=["A", "B", "C"]),
                         [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 2,
                              "0.333333")])
        # D's two haplotypes carry them too.
        self.assertEqual(lines(row, geno),
                         [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 4,
                              "0.5")])

    def test_frequency_format(self):
        for carriers, total, text in ((2, 96, "0.0208333"), (1, 3, "0.333333"),
                                      (8, 8, "1"), (1, 2, "0.5"), (1, 8, "0.125"),
                                      (2, 30000, "0.0000666667")):
            with self.subTest(carriers=carriers, total=total):
                self.assertEqual(merge.format_frequency(carriers, total), text)

    def test_discordance_counted_strictly_between_the_segments(self):
        # The input segments' own end markers are discordant here (input
        # Refined IBD would not write); only 31000, strictly between them,
        # counts, so the gap has one discordant marker and is joined.
        # (Browning's count includes both end markers, but counts a marker
        # only where the genotypes share no allele.)
        geno = (Genotypes().set(30000, "A", 1).set(31000, "A", 1)
                .set(32000, "A", 1).set(61000, "A", 1))
        rows = [seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 32000, 60000)]
        self.assertEqual(report_lines(rows, geno),
                         [REPORT_HEADER,
                          rep(ONE, "A", 1, "B", 1, 1000, 30000, 32000, 60000,
                              1, 1, "0.25")])

    def test_marker_count_includes_records_sharing_a_position(self):
        extra = "1\t30000\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0"
        self.assertEqual(lines([seg("A", 1, "B", 1, 1000, 40000)],
                               after={30000: extra}),
                         [out("A", 1, "B", 1, 1000, 40000, "3", "4.875", 41, 8,
                              "1")])

    def test_extrapolated_length_and_rounding(self):
        # Markers every 100 bp; 0 cM at 500 and 0.6 cM at 1100, 3.0 at 3100.
        # genPos(100) = (100 - 500) / 2600 * 3.0 = -0.4615...; genPos(1200)
        # = 0.72; length 1.1815... -> 1.182.  Score 2.675 is the float
        # 2.6749999523..., printed 2.67.
        old = dict(positions=list(range(100, 3100, 100)),
                   map_entries=[("1", "0.0", 500), ("1", "0.6", 1100),
                                ("1", "3.0", 3100)])
        self.assertEqual(lines([seg("A", 1, "B", 1, 100, 1200, "2.675")], **old),
                         [out("A", 1, "B", 1, 100, 1200, "2.67", "1.182", 12, 8,
                              "1")])

    def test_extrapolated_length_past_the_end_of_the_map(self):
        # The map ends at 2000 and spans less than 5 cM, so beyond it the
        # line runs from its first point to its last: genPos(2500) = 2.4.
        short = dict(positions=list(range(100, 3100, 100)),
                     map_entries=[("1", "0.0", 500), ("1", "0.6", 1100),
                                  ("1", "1.8", 2000)])
        self.assertEqual(lines([seg("A", 1, "B", 1, 1100, 2500)], **short),
                         [out("A", 1, "B", 1, 1100, 2500, "3", "1.8", 15, 8, "1")])

    def test_number_rounding(self):
        # Scores are float32 values rounded half to even: 0.045 is the float
        # 0.0450000017..., printed 0.05; 0.125 is exact and prints 0.12.
        # The length 0.0625 cM is exact and prints 0.062.
        tie = dict(positions=list(range(100, 3100, 100)),
                   map_entries=[("1", "0.0", 500), ("1", "0.0625", 600),
                                ("1", "3.0", 3100)])
        self.assertEqual(lines([seg("A", 1, "B", 1, 500, 600, "0.045"),
                                seg("A", 2, "B", 2, 500, 600, "0.125")], **tie),
                         [out("A", 1, "B", 1, 500, 600, "0.05", "0.062", 2, 8, "1"),
                          out("A", 2, "B", 2, 500, 600, "0.12", "0.062", 2, 8, "1")])


class Errors(unittest.TestCase):

    def test_discord_must_be_0_or_1(self):
        for arg, message in (("2", "ERROR: [discord] must be 0 or 1: 2\n"),
                             ("-1", "ERROR: [discord] < 0: -1\n"),
                             ("x", "ERROR: [discord] is not a nonnegative "
                                   "integer: x\n")):
            with self.subTest(arg=arg):
                t = termination([seg("A", 1, "B", 1, 1000, 30000)], discord=arg)
                # The usage text goes to standard error with the message, so
                # it never lands in [out].
                self.assertEqual(t.status, 1)
                self.assertEqual(t.stdout, "")
                self.assertEqual(t.stderr, merge.usage() + "\n" + message)

    def test_number_of_arguments(self):
        with self.assertRaises(JavaTermination) as ctx:
            merge.run([], io.BytesIO(b""))
        self.assertEqual((ctx.exception.status, ctx.exception.stdout,
                          ctx.exception.stderr), (0, merge.usage() + "\n", ""))
        for args in (["a", "b", "0.6"], ["a", "b", "0.6", "1", "r", "x"]):
            with self.subTest(n=len(args)):
                with self.assertRaises(JavaTermination) as ctx:
                    merge.run(args, io.BytesIO(b""))
                self.assertEqual(ctx.exception.status, 1)
                self.assertEqual(ctx.exception.stdout, "")
                self.assertEqual(ctx.exception.stderr,
                                 merge.usage() + "\nERROR: expected 4 or 5 "
                                 "arguments, found " + str(len(args)) + "\n")

    def test_vcf_records_in_position_order_in_one_block_per_chromosome(self):
        rows = [seg("A", 1, "B", 1, 1000, 30000)]
        swapped = POSITIONS[:]
        swapped[9], swapped[10] = swapped[10], swapped[9]
        t = termination(rows, positions=swapped)
        self.assertEqual(t.stderr, "ERROR: VCF records are not in position "
                                   "order: chromosome 1 position 10000 follows "
                                   "position 11000\n")
        # A chromosome 2 record between two chromosome 1 records.
        mp = MAP + [("2", "0.0", 5000), ("2", "1.0", 9000)]
        t = termination(rows, map_entries=mp,
                        after={100000: "2\t5000\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0"})
        self.assertEqual(t.stderr, "ERROR: the records of chromosome 1 are not "
                                   "in one block in the VCF file\n")

    def test_segment_end_must_be_one_vcf_marker(self):
        t = termination([seg("A", 1, "B", 1, 1000, 30500)])
        self.assertEqual(t.stderr,
                         "ERROR: segment end is not a VCF marker: position "
                         "30500 of A haplotype 1, B haplotype 1, chromosome 1, "
                         "1000-30500\n")
        extra = "1\t30000\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0\t0|0"
        t = termination([seg("A", 1, "B", 1, 1000, 30000)], after={30000: extra})
        self.assertEqual(t.stderr,
                         "ERROR: more than one VCF record at a segment end: "
                         "position 30000 of A haplotype 1, B haplotype 1, "
                         "chromosome 1, 1000-30000\n")

    def test_sample_must_be_in_the_vcf_file(self):
        for row in (seg("E", 1, "A", 1, 1000, 30000), seg("A", 1, "E", 1, 1000, 30000)):
            with self.subTest(row=row):
                t = termination([row])
                self.assertEqual(t.stderr, "ERROR: sample not in VCF file: E\n")

    def test_haplotype_zero_rejected(self):
        for row in (seg("A", 0, "B", 1, 1000, 30000), seg("A", 1, "B", 0, 1000, 30000)):
            with self.assertRaises(JavaException) as ctx:
                run([row])
            self.assertEqual(ctx.exception.to_string(),
                             "java.lang.IllegalArgumentException: invalid hap: 0")

    def test_no_output_and_no_report_after_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            vcf, mp = write_inputs(d, Genotypes())
            report = os.path.join(d, "report.tsv")
            rows = seg("A", 1, "B", 1, 1000, 40000) + "\n" + seg("A", 1, "B", 1, 30000, 60000)
            with self.assertRaises(JavaTermination):
                merge.run([vcf, mp, "0.6", "1", report],
                          io.BytesIO(rows.encode()))
            self.assertFalse(os.path.exists(report))

    def test_failed_output_leaves_no_new_report(self):
        # The segments cannot be written: the report is not written either,
        # no temporary file remains, and an earlier report is left as it was.

        class Broken(io.BytesIO):
            def write(self, data):
                raise OSError(28, "No space left on device")

        with tempfile.TemporaryDirectory() as d:
            vcf, mp = write_inputs(d, Genotypes())
            report = os.path.join(d, "report.tsv")
            for earlier in (None, b"earlier report\n"):
                with self.subTest(earlier=earlier):
                    if earlier is not None:
                        with open(report, "wb") as f:
                            f.write(earlier)
                    result = merge.run([vcf, mp, "0.6", "1", report], io.BytesIO(
                        (seg("A", 1, "B", 1, 1000, 30000) + "\n").encode()))
                    with self.assertRaises(JavaTermination) as ctx:
                        merge.write_outputs(result, Broken())
                    self.assertEqual(ctx.exception.stderr,
                                     "ERROR: cannot write the merged segments: "
                                     "[Errno 28] No space left on device\n")
                    expected = ["in.map", "in.vcf"] + ([] if earlier is None
                                                       else ["report.tsv"])
                    self.assertEqual(sorted(os.listdir(d)), expected)
                    if earlier is not None:
                        with open(report, "rb") as f:
                            self.assertEqual(f.read(), earlier)

    def test_sample_names_beyond_the_basic_plane_in_the_report(self):
        # A sample name outside the Basic Multilingual Plane is held as a
        # UTF-16 surrogate pair, as Java holds it, and written as UTF-8.
        name = "A\U0001F600"
        geno = Genotypes().set(31000, name, 1).set(61000, name, 1)
        got = run([seg(name, 1, "B", 1, 1000, 30000), seg(name, 1, "B", 1, 33000, 60000)],
                  geno, samples=[name, "B", "C", "D"])
        written = java_encode(got.report).decode("utf-8").splitlines()
        self.assertEqual(written[1:],
                         [rep(ONE, name, 1, "B", 1, 1000, 30000, 33000, 60000,
                              2, 1, "0.375"),
                          rep(START, name, 1, "B", 1, 33000, 60000)])


    def test_report_path_checked_before_anything_is_written(self):
        with tempfile.TemporaryDirectory() as d:
            vcf, mp = write_inputs(d, Genotypes())
            os.mkdir(os.path.join(d, "folder"))
            locked = os.path.join(d, "locked.tsv")
            with open(locked, "w") as f:
                f.write("earlier\n")
            os.chmod(locked, 0o444)
            cases = (("", "ERROR: [report] is empty"),
                     (os.path.join(d, "folder"), "ERROR: [report] is a directory: "),
                     (os.path.join(d, "missing", "r.tsv"),
                      "ERROR: the directory of [report] does not exist: "),
                     (locked, "ERROR: [report] is not writable: "),
                     (vcf, "ERROR: [report] is the same file as [vcf]: "),
                     (mp, "ERROR: [report] is the same file as [map]: "))
            try:
                for path, message in cases:
                    with self.subTest(message=message):
                        with self.assertRaises(JavaTermination) as ctx:
                            merge.run([vcf, mp, "0.6", "1", path], io.BytesIO(
                                (seg("A", 1, "B", 1, 1000, 30000) + "\n").encode()))
                        self.assertEqual(ctx.exception.status, 1)
                        self.assertEqual(ctx.exception.stdout, "")
                        self.assertTrue(ctx.exception.stderr.startswith(
                            merge.usage() + "\n" + message))
            finally:
                os.chmod(locked, 0o666)
            with open(vcf, encoding="utf-8") as f:
                self.assertTrue(f.read().startswith("##fileformat=VCFv4.2"))

    def test_allele_index_above_255(self):
        # REF A and 256 different ALT alleles, so that allele index 256 is
        # valid VCF.
        alleles = ",".join(["G"] + ["A" * (i + 2) for i in range(255)])
        record = "1\t150500\t.\tA\t%s\t.\tPASS\t.\tGT\t256|0\t0|0\t0|0\t0|0" % alleles
        t = termination([seg("A", 1, "B", 1, 1000, 30000)], after={150000: record})
        self.assertEqual(t.stderr, "ERROR: allele index above 255 at chromosome 1 "
                                   "position 150500\n")


class CommandLine(unittest.TestCase):

    def test_closed_output_pipe(self):
        # The reader of standard output has gone before the segments are
        # written: an error message and status 1, not a traceback.
        with tempfile.TemporaryDirectory() as d:
            vcf, mp = write_inputs(d, Genotypes())
            proc = subprocess.Popen(
                [sys.executable, "-m", "merge_ibs_segments", vcf, mp, "0.6", "1"],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, cwd=REPO)
            proc.stdout.close()
            proc.stdin.write((seg("A", 1, "B", 1, 1000, 30000) + "\n").encode())
            proc.stdin.close()
            err = proc.stderr.read()
            proc.stderr.close()
            proc.wait()
        self.assertEqual(proc.returncode, 1)
        self.assertTrue(err.decode().startswith(
            "ERROR: cannot write the merged segments: "), err.decode())
        self.assertNotIn("Traceback", err.decode())

    def test_output_report_and_summary(self):
        geno = Genotypes().set(61000, "A", 1)
        rows = [seg("A", 1, "B", 1, 1000, 30000), seg("A", 1, "B", 1, 33000, 60000)]
        with tempfile.TemporaryDirectory() as d:
            vcf, mp = write_inputs(d, geno)
            report = os.path.join(d, "report.tsv")
            done = subprocess.run(
                [sys.executable, "-m", "merge_ibs_segments", vcf, mp, "0.6", "1",
                 report],
                input=("\n".join(rows) + "\n").encode(), capture_output=True,
                cwd=REPO)
            with open(report, encoding="utf-8", newline="") as f:
                report_text = f.read()
        self.assertEqual(done.returncode, 0)
        self.assertEqual(done.stdout.decode(),
                         out("A", 1, "B", 1, 1000, 60000, "3", "7.375", 60, 8, "1")
                         + "\n")
        self.assertEqual(done.stderr.decode(),
                         "merge-ibs-segments: 2 segments read, 1 written; "
                         "joined: 0 across one discordant marker, 1 across a "
                         "gap with no discordant marker (0 of them longer than "
                         "[gap]), 0 at a shared end marker; segment ends whose "
                         "neighbouring marker matches: 2\n")
        self.assertEqual(report_text.splitlines()[0], REPORT_HEADER)
        self.assertEqual(len(report_text.splitlines()), 4)

    def test_report_is_optional(self):
        with tempfile.TemporaryDirectory() as d:
            vcf, mp = write_inputs(d, Genotypes())
            result = merge.run([vcf, mp, "0.6", "1"], io.BytesIO(
                (seg("A", 1, "B", 1, 1000, 30000) + "\n").encode()))
            stdout = io.BytesIO()
            merge.write_outputs(result, stdout)
            self.assertEqual(sorted(os.listdir(d)), ["in.map", "in.vcf"])
        self.assertEqual(stdout.getvalue().decode().splitlines(),
                         [out("A", 1, "B", 1, 1000, 30000, "3", "3.625", 30, 8, "1")])


# ---------------------------------------------------------------------------
# Seeded random cases against a separate implementation of the rule
# ---------------------------------------------------------------------------

def oracle(rows: List[Tuple], haps: Dict[str, List[List[int]]],
           positions: Dict[str, List[int]], gen_pos, gap: float, discord: int,
           n_samples: int, sample_col: Dict[str, int]):
    """The rule, implemented independently of merge.py: for each group
    (sample 1, haplotype 1, sample 2, haplotype 2, chromosome), keep a list
    of chains ordered by start and join the first neighbouring pair that
    qualifies until no pair does.  Returns the output lines and the report
    lines (without the header)."""
    groups = defaultdict(list)
    for s1, h1, s2, h2, chrom, start, end, score in rows:
        groups[(s1, h1, s2, h2, chrom)].append((start, end, java_parse_float(score)))
    out_lines, report = [], []
    n_haps = 2 * n_samples
    for (s1, h1, s2, h2, chrom), pieces in groups.items():
        a = haps[chrom][2 * sample_col[s1] + h1 - 1]
        b = haps[chrom][2 * sample_col[s2] + h2 - 1]
        index = {p: m for m, p in enumerate(positions[chrom])}
        pos = positions[chrom]
        for start, end, _ in pieces:
            i, j = index[start], index[end]
            if i > 0 and a[i - 1] == b[i - 1]:
                report.append((START, s1, h1, s2, h2, chrom, start, end,
                               ".", ".", ".", ".", "."))
            if j + 1 < len(pos) and a[j + 1] == b[j + 1]:
                report.append((END, s1, h1, s2, h2, chrom, start, end,
                               ".", ".", ".", ".", "."))
        chains = [[p] for p in sorted(pieces, key=lambda p: (p[0], p[1]))]
        joined = True
        while joined:
            joined = False
            for k in range(len(chains) - 1):
                left, right = chains[k][-1], chains[k + 1][0]
                if right[0] == left[1]:
                    event, fields = SHARED, ("0", "0", "0")
                else:
                    i, j = index[left[1]], index[right[0]]
                    bad = sum(1 for m in range(i + 1, j) if a[m] != b[m])
                    length = gen_pos(chrom, right[0]) - gen_pos(chrom, left[1])
                    fields = (str(j - i - 1), str(bad),
                              java_decimal_format(length, 3))
                    if bad == 0:
                        event = NONE_SHORT if length <= gap else NONE_LONG
                    elif bad <= discord and length <= gap:
                        event = ONE
                    else:
                        event = None
                if event is not None:
                    report.append((event, s1, h1, s2, h2, chrom, left[0], left[1],
                                   right[0], right[1]) + fields)
                    chains[k:k + 2] = [chains[k] + chains[k + 1]]
                    joined = True
                    break
        for chain in chains:
            start, end = chain[0][0], chain[-1][1]
            score = max(p[2] for p in chain)
            i, j = index[start], index[end]
            agree = [m for m in range(i, j + 1) if a[m] == b[m]]
            carriers = sum(1 for h in haps[chrom] if all(h[m] == a[m] for m in agree))
            length = gen_pos(chrom, end) - gen_pos(chrom, start)
            out_lines.append("\t".join((
                s1, str(h1), s2, str(h2), chrom, str(start), str(end),
                java_decimal_format(score, 2), java_decimal_format(length, 3),
                str(j - i + 1), str(carriers), "%.6g" % (carriers / n_haps))))
    report_lines = ["\t".join(str(x) for x in
                              (r[0], "yes" if merge.FLAGGED[r[0]] else "no") + r[1:])
                    for r in report]
    return out_lines, report_lines


class RandomAgainstOracle(unittest.TestCase):

    def test_seeded_random_cases(self):
        seen = Counter()
        for seed in range(1, 31):
            with self.subTest(seed=seed):
                seen.update(self._check(seed))
        # Every kind of join and flag occurred somewhere.
        for event, _ in merge.EVENTS:
            self.assertGreater(seen[event], 0, event)

    def _check(self, seed: int) -> Counter:
        rng = random.Random(seed)
        samples = ["r%d" % i for i in range(rng.randint(4, 6))]
        col = {s: k for k, s in enumerate(samples)}
        chroms = ["1", "2"]
        n_markers = 220
        positions = {c: sorted(rng.sample(range(1000, 2000000, 97), n_markers))
                     for c in chroms}
        haps: Dict[str, List[List[int]]] = {}
        for c in chroms:
            founders = [[rng.randint(0, 1) for _ in range(n_markers)]
                        for _ in range(3)]
            hs = []
            for _ in range(2 * len(samples)):
                f = rng.randrange(3)
                h = []
                for m in range(n_markers):
                    if rng.random() < 0.02:
                        f = rng.randrange(3)
                    allele = founders[f][m]
                    if rng.random() < 0.01:
                        allele = 1 - allele
                    h.append(allele)
                hs.append(h)
            haps[c] = hs
        map_entries = []
        for c in chroms:
            cm = 0.0
            for p in sorted(rng.sample(positions[c], 30)) + [positions[c][-1] + 500]:
                cm += rng.uniform(0.2, 1.5)
                map_entries.append((c, repr(round(cm, 4)), p))

        with tempfile.TemporaryDirectory() as d:
            mp = os.path.join(d, "in.map")
            with open(mp, "w", encoding="utf-8", newline="\n") as f:
                for c, cm, bp in map_entries:
                    f.write("%s\t.\t%s\t%d\n" % (c, cm, bp))
            chrom_ids = Indexer()
            gmap = PlinkGenMap.from_plink_map_file(mp, chrom_ids)

            def gen_pos(chrom, bp):
                return gmap.gen_pos(chrom_ids.get_index(chrom), bp)

            # Maximal identical runs of every haplotype pair of two samples,
            # at least 0.8 cM long, cut in some of the ways Refined IBD can.
            rows = []
            for c in chroms:
                pos = positions[c]
                hs = haps[c]
                for x in range(len(hs)):
                    for y in range(len(hs)):
                        if x // 2 >= y // 2:
                            continue
                        m = 0
                        while m < n_markers:
                            if hs[x][m] != hs[y][m]:
                                m += 1
                                continue
                            a = m
                            while m < n_markers and hs[x][m] == hs[y][m]:
                                m += 1
                            b = m - 1
                            if gen_pos(c, pos[b]) - gen_pos(c, pos[a]) < 0.8:
                                continue
                            pieces = [(a, b)]
                            r = rng.random()
                            if r < 0.15 and b - a >= 4:
                                pieces = [(a + rng.randint(0, 2), b - rng.randint(0, 2))]
                            elif r < 0.35 and b - a >= 4:
                                c1 = rng.randint(a, (a + b) // 2)
                                c2 = rng.randint(c1 + 1, b)
                                pieces = [(a, c1), (c2, b)]
                            elif r < 0.45 and b - a >= 2:
                                # Both pieces extend beyond the shared
                                # marker (a one-marker piece there would be
                                # nested, an error).
                                k = rng.randint(a + 1, b - 1)
                                pieces = [(a, k), (k, b)]
                            for i, j in pieces:
                                rows.append((samples[x // 2], x % 2 + 1,
                                             samples[y // 2], y % 2 + 1, c,
                                             pos[i], pos[j],
                                             "%.2f" % rng.uniform(0, 30)))
            rng.shuffle(rows)
            gap = rng.choice([0.05, 0.3, 0.6, 1.0])
            discord = rng.choice([0, 1])

            vcf = os.path.join(d, "in.vcf")
            with open(vcf, "w", encoding="utf-8", newline="\n") as f:
                f.write("##fileformat=VCFv4.2\n" + HEADER
                        + "".join("\t" + s for s in samples) + "\n")
                for c in chroms:
                    for m, p in enumerate(positions[c]):
                        f.write("\t".join([c, str(p), ".", "A", "G", ".", "PASS",
                                           ".", "GT"]
                                          + ["%d|%d" % (haps[c][2 * k][m],
                                                        haps[c][2 * k + 1][m])
                                             for k in range(len(samples))]) + "\n")

            def program(order):
                stdin = io.BytesIO("".join(
                    "\t".join(str(x) for x in r) + "\n" for r in order).encode())
                return merge.run([vcf, mp, repr(gap), str(discord)], stdin)

            got = program(rows)
            reordered = rows[:]
            random.Random(seed + 1000).shuffle(reordered)
            again = program(reordered)
        self.assertEqual(again, got)
        expected_out, expected_report = oracle(rows, haps, positions, gen_pos, gap,
                                               discord, len(samples), col)
        self.assertEqual(sorted(got.output.splitlines()), sorted(expected_out))
        got_report = got.report.splitlines()
        self.assertEqual(got_report[0], REPORT_HEADER)
        self.assertEqual(sorted(got_report[1:]), sorted(expected_report))
        c = Counter(line.split("\t")[0] for line in expected_report)
        self.assertEqual(
            got.summary,
            "merge-ibs-segments: %d segments read, %d written; joined: %d "
            "across one discordant marker, %d across a gap with no discordant "
            "marker (%d of them longer than [gap]), %d at a shared end marker; "
            "segment ends whose neighbouring marker matches: %d\n"
            % (len(rows), len(expected_out), c[ONE], c[NONE_SHORT] + c[NONE_LONG],
               c[NONE_LONG], c[SHARED], c[START] + c[END]))
        return c


if __name__ == "__main__":
    unittest.main()
