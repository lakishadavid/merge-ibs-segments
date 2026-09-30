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
"""Tests of the haplotype-specific merging rule.

Run from the repository root: ``python -m unittest discover -s tests/rule``.

The hand-built cases each state their expected output, worked out from the
rule.  The seeded random cases compare the program with ``oracle_merge``, a
separate implementation of the rule written for this test: it joins
neighbouring chains of a haplotype group until nothing more can be joined,
instead of walking the segments once.  Both use the program's genetic map and
number formatting, which the exact-conversion test of the baseline checked
against Browning's jar.

All data are synthetic and built here.
"""

from __future__ import annotations

import io
import os
import random
import sys
import tempfile
import unittest
from collections import defaultdict
from typing import Dict, List, Sequence, Tuple

REPO = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from merge_ibs_segments import merge  # noqa: E402
from merge_ibs_segments.genetic_map import PlinkGenMap  # noqa: E402
from merge_ibs_segments.ids import Indexer  # noqa: E402
from merge_ibs_segments.java_compat import (  # noqa: E402
    JavaException, java_decimal_format, java_parse_float)

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
SAMPLES = ["A", "B", "C"]
POSITIONS = list(range(100, 3100, 100))        # markers at 100, 200, ..., 3000
# 0 cM at 500 and exactly 0.6 cM at 1100, so the gap from 500 to 1100 is 0.6.
MAP = [("1", "0.0", 500), ("1", "0.6", 1100), ("1", "3.0", 3100)]


def write_inputs(directory: str, genotypes: Dict[int, Dict[str, str]],
                 chrom: str = "1", positions: Sequence[int] = POSITIONS,
                 samples: Sequence[str] = SAMPLES,
                 map_entries: Sequence[Tuple[str, str, int]] = MAP,
                 extra_records: Sequence[str] = ()) -> Tuple[str, str]:
    """Writes a phased VCF (every genotype 0|0 unless ``genotypes`` gives
    another for a position and sample) and a PLINK map."""
    lines = ["##fileformat=VCFv4.2", HEADER + "".join("\t" + s for s in samples)]
    for pos in positions:
        gts = [genotypes.get(pos, {}).get(s, "0|0") for s in samples]
        lines.append("\t".join([chrom, str(pos), ".", "A", "G", ".", "PASS",
                                ".", "GT"] + gts))
    lines.extend(extra_records)
    vcf = os.path.join(directory, "in.vcf")
    with open(vcf, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    mp = os.path.join(directory, "in.map")
    with open(mp, "w", encoding="utf-8", newline="\n") as f:
        for c, cm, bp in map_entries:
            f.write("%s\t.\t%s\t%d\n" % (c, cm, bp))
    return vcf, mp


def run(rows: Sequence[str], genotypes: Dict[int, Dict[str, str]] = None,
        gap: str = "0.6", discord: str = "1", **kw) -> List[str]:
    """Runs the program on segment lines; returns the output lines."""
    with tempfile.TemporaryDirectory() as d:
        vcf, mp = write_inputs(d, genotypes or {}, **kw)
        stdin = io.BytesIO(("\n".join(rows) + "\n").encode("utf-8"))
        return merge.run([vcf, mp, gap, discord], stdin).splitlines()


def seg(s1, h1, s2, h2, start, end, score="3", chrom="1") -> str:
    return "\t".join(str(x) for x in (s1, h1, s2, h2, chrom, start, end, score))


def out(s1, h1, s2, h2, start, end, score, length, chrom="1") -> str:
    return "\t".join(str(x) for x in (s1, h1, s2, h2, chrom, start, end,
                                      score, length))


# Lengths from the map: genPos is linear from 0 cM at 500 to 0.6 at 1100 and
# from 0.6 at 1100 to 3.0 at 3100 (0.0012 cM per bp); before 500 it is
# extrapolated along the line from 500 to 3100 (the first map position at
# least 5 cM inside is past the end, so the other end is used).


class HaplotypeRule(unittest.TestCase):

    def test_same_haplotypes_across_gap_merge_and_keep_numbers(self):
        for h2 in (1, 2):
            with self.subTest(b_haplotype=h2):
                got = run([seg("A", 1, "B", h2, 500, 800, "3"),
                           seg("A", 1, "B", h2, 1100, 1500, "5")])
                # From the running end 800 to 1100 is a 0.3 cM gap with no
                # discordant marker.
                self.assertEqual(got, [out("A", 1, "B", h2, 500, 1500, "5",
                                           "1.08")])

    def test_a_changes_not_merged(self):
        got = run([seg("A", 1, "B", 2, 500, 800), seg("A", 2, "B", 2, 1100, 1500)])
        self.assertEqual(got, [out("A", 1, "B", 2, 500, 800, "3", "0.3"),
                               out("A", 2, "B", 2, 1100, 1500, "3", "0.48")])

    def test_b_changes_not_merged(self):
        got = run([seg("A", 1, "B", 2, 500, 800), seg("A", 1, "B", 1, 1100, 1500)])
        self.assertEqual(got, [out("A", 1, "B", 2, 500, 800, "3", "0.3"),
                               out("A", 1, "B", 1, 1100, 1500, "3", "0.48")])

    def test_both_change_not_merged(self):
        got = run([seg("A", 1, "B", 1, 500, 800), seg("A", 2, "B", 2, 1100, 1500)])
        self.assertEqual(got, [out("A", 1, "B", 1, 500, 800, "3", "0.3"),
                               out("A", 2, "B", 2, 1100, 1500, "3", "0.48")])

    def test_other_haplotype_pair_across_gap_does_not_block(self):
        # A segment of another haplotype pair (A2-B2, 600-1300) spans the gap
        # 800-1100 between the two A1-B1 pieces; it does not prevent the join.
        got = run([seg("A", 1, "B", 1, 500, 800), seg("A", 2, "B", 2, 600, 1300),
                   seg("A", 1, "B", 1, 1100, 1500)])
        self.assertEqual(got, [out("A", 1, "B", 1, 500, 1500, "3", "1.08"),
                               out("A", 2, "B", 2, 600, 1300, "3", "0.74")])

    def test_segment_starting_at_running_end_is_always_gap_tested(self):
        # A one-marker segment at the running end lies inside the chain and
        # also starts at its end; it is tested as a zero-length gap.  (With a
        # 3 cM minimum segment length such a segment cannot occur, and with
        # [discord] of 1 or more the test always passes.)
        rows = [seg("A", 1, "B", 1, 500, 800), seg("A", 1, "B", 1, 800, 800)]
        self.assertEqual(run(rows, {800: {"A": "1|0"}}, discord="0"),
                         [out("A", 1, "B", 1, 500, 800, "3", "0.3"),
                          out("A", 1, "B", 1, 800, 800, "3", "0")])
        self.assertEqual(run(rows, {800: {"A": "1|0"}}, discord="1"),
                         [out("A", 1, "B", 1, 500, 800, "3", "0.3")])

    def test_overlapping_different_pairs_not_merged(self):
        got = run([seg("A", 1, "B", 1, 500, 1100), seg("A", 1, "B", 2, 800, 1500)])
        self.assertEqual(got, [out("A", 1, "B", 1, 500, 1100, "3", "0.6"),
                               out("A", 1, "B", 2, 800, 1500, "3", "0.78")])

    def test_gap_exactly_at_limit_merges_and_just_over_does_not(self):
        rows = [seg("A", 1, "B", 1, 500, 500), seg("A", 1, "B", 1, 1100, 1500)]
        self.assertEqual(run(rows, gap="0.6"),
                         [out("A", 1, "B", 1, 500, 1500, "3", "1.08")])
        # 0.5999999999999999 is the largest double below 0.6.
        self.assertEqual(run(rows, gap="0.5999999999999999"),
                         [out("A", 1, "B", 1, 500, 500, "3", "0"),
                          out("A", 1, "B", 1, 1100, 1500, "3", "0.48")])

    def test_discordance_threshold(self):
        # A's haplotype 1 carries 1 and B's haplotype 1 carries 0 at the
        # chosen gap markers (both boundary markers count).
        rows = [seg("A", 1, "B", 1, 500, 800), seg("A", 1, "B", 1, 1100, 1500)]
        one = {800: {"A": "1|0"}}
        two = {800: {"A": "1|0"}, 1100: {"A": "1|0"}}
        self.assertEqual(run(rows, one, discord="1"),
                         [out("A", 1, "B", 1, 500, 1500, "3", "1.08")])
        self.assertEqual(run(rows, two, discord="1"),
                         [out("A", 1, "B", 1, 500, 800, "3", "0.3"),
                          out("A", 1, "B", 1, 1100, 1500, "3", "0.48")])
        self.assertEqual(run(rows, one, discord="0"),
                         [out("A", 1, "B", 1, 500, 800, "3", "0.3"),
                          out("A", 1, "B", 1, 1100, 1500, "3", "0.48")])

    def test_only_the_respective_haplotypes_count(self):
        rows = [seg("A", 1, "B", 1, 500, 800), seg("A", 1, "B", 1, 1100, 1500)]
        # At 900 and 1000 A is 0|1 and B is 0|0: the second haplotypes
        # differ, but A's haplotype 1 equals B's haplotype 1.
        other = {900: {"A": "0|1"}, 1000: {"A": "0|1", "B": "0|0"}}
        self.assertEqual(run(rows, other, discord="0"),
                         [out("A", 1, "B", 1, 500, 1500, "3", "1.08")])
        # Opposite alleles on the respective haplotypes while the genotypes
        # share an allele (A 1|0, B 0|1): discordant here, although
        # Browning's genotype test would not count it.
        shared = {900: {"A": "1|0", "B": "0|1"}}
        self.assertEqual(run(rows, shared, discord="0"),
                         [out("A", 1, "B", 1, 500, 800, "3", "0.3"),
                          out("A", 1, "B", 1, 1100, 1500, "3", "0.48")])

    def test_shared_endpoint_is_tested_as_zero_gap(self):
        rows = [seg("A", 2, "B", 1, 500, 800), seg("A", 2, "B", 1, 800, 1100)]
        self.assertEqual(run(rows, discord="0"),
                         [out("A", 2, "B", 1, 500, 1100, "3", "0.6")])
        # The shared marker is discordant on the respective haplotypes.
        bad = {800: {"A": "0|1"}}
        self.assertEqual(run(rows, bad, discord="0"),
                         [out("A", 2, "B", 1, 500, 800, "3", "0.3"),
                          out("A", 2, "B", 1, 800, 1100, "3", "0.3")])

    def test_overlapping_and_nested_same_pair_merge_without_test(self):
        # Discordant markers inside the overlap do not matter.
        bad = {900: {"A": "1|1"}, 1000: {"A": "1|1"}}
        rows = [seg("A", 1, "B", 2, 500, 1000, "2"), seg("A", 1, "B", 2, 800, 1500, "4"),
                seg("A", 1, "B", 2, 600, 700, "9")]
        self.assertEqual(run(rows, bad, discord="0"),
                         [out("A", 1, "B", 2, 500, 1500, "9", "1.08")])

    def test_chain_of_four_and_shuffled_input(self):
        # Four pieces joined by gaps of 0.2, 0.2 and 0.24 cM, then a fifth
        # behind a gap of 0.72 cM (2000 -> 2600), which is not joined.
        rows = [seg("A", 2, "C", 2, 500, 600, "1"), seg("A", 2, "C", 2, 800, 900, "4"),
                seg("A", 2, "C", 2, 1100, 1300, "2"), seg("A", 2, "C", 2, 1500, 2000, "3"),
                seg("A", 2, "C", 2, 2600, 2900, "6")]
        expected = [out("A", 2, "C", 2, 500, 2000, "4", "1.68"),
                    out("A", 2, "C", 2, 2600, 2900, "6", "0.36")]
        self.assertEqual(run(rows), expected)
        for seed in range(5):
            shuffled = rows[:]
            random.Random(seed).shuffle(shuffled)
            self.assertEqual(run(shuffled), expected)

    def test_output_order(self):
        # Output is sorted by chromosome, start, end, score, samples and
        # haplotypes, whatever the input order and group.
        # [500,700] is nested in [500,800] of the same group and joins it.
        rows = [seg("B", 1, "C", 1, 1100, 1500), seg("A", 2, "B", 2, 1100, 1500),
                seg("A", 1, "C", 2, 500, 800), seg("A", 1, "B", 2, 500, 800),
                seg("A", 1, "B", 2, 500, 700)]
        self.assertEqual(run(rows),
                         [out("A", 1, "B", 2, 500, 800, "3", "0.3"),
                          out("A", 1, "C", 2, 500, 800, "3", "0.3"),
                          out("A", 2, "B", 2, 1100, 1500, "3", "0.48"),
                          out("B", 1, "C", 1, 1100, 1500, "3", "0.48")])

    def test_gap_measured_from_running_end_after_nested_segment(self):
        # [600,700] lies inside [500,1000]; the gap to [1100,1500] starts at
        # the running end 1000, so the discordant marker at 800 is not in it.
        # Measured from 700, it would be, and discord 0 would split them.
        rows = [seg("A", 1, "B", 1, 500, 1000), seg("A", 1, "B", 1, 600, 700),
                seg("A", 1, "B", 1, 1100, 1500)]
        self.assertEqual(run(rows, {800: {"A": "1|0"}}, discord="0"),
                         [out("A", 1, "B", 1, 500, 1500, "3", "1.08")])

    def test_sample_two_haplotype_two_is_the_one_compared(self):
        rows = [seg("A", 1, "B", 2, 500, 800), seg("A", 1, "B", 2, 1100, 1500)]
        # B's second allele differs from A's first allele: discordant.
        self.assertEqual(run(rows, {900: {"B": "0|1"}}, discord="0"),
                         [out("A", 1, "B", 2, 500, 800, "3", "0.3"),
                          out("A", 1, "B", 2, 1100, 1500, "3", "0.48")])
        # B's first allele differs instead: not counted.
        self.assertEqual(run(rows, {900: {"B": "1|0"}}, discord="0"),
                         [out("A", 1, "B", 2, 500, 1500, "3", "1.08")])

    def test_shared_endpoint_must_be_a_vcf_marker(self):
        # 850 is not a VCF marker; the zero-length gap at it cannot be tested.
        rows = [seg("A", 1, "B", 1, 500, 850), seg("A", 1, "B", 1, 850, 1100)]
        with self.assertRaises(merge.JavaTermination) as ctx:
            run(rows)
        self.assertEqual(ctx.exception.stderr,
                         "ERROR: position missing from VCF file: 0:850\n")

    def test_haplotype_zero_rejected(self):
        for row in (seg("A", 0, "B", 1, 500, 800), seg("A", 1, "B", 0, 500, 800)):
            with self.assertRaises(JavaException) as ctx:
                run([row])
            self.assertEqual(ctx.exception.to_string(),
                             "java.lang.IllegalArgumentException: invalid hap: 0")

    def test_extrapolated_length_and_rounding(self):
        # Before the map: genPos(100) = (100 - 500) / 2600 * 3.0 = -0.4615...
        # genPos(1200) = 0.6 + 100 * 0.0012 = 0.72; length 1.1815... -> 1.182.
        # Score 2.675 is the float 2.6749999523..., printed 2.67.
        self.assertEqual(run([seg("A", 1, "B", 1, 100, 1200, "2.675")]),
                         [out("A", 1, "B", 1, 100, 1200, "2.67", "1.182")])

    def test_extrapolated_length_past_the_end_of_the_map(self):
        # The map ends at 2000 and spans less than 5 cM, so beyond it the
        # line runs from its first point (500, 0 cM) to its last (2000,
        # 1.8 cM): genPos(2500) = 2000 / 1500 * 1.8 = 2.4, and the length
        # from 1100 (0.6 cM) is 1.8.
        short_map = [("1", "0.0", 500), ("1", "0.6", 1100), ("1", "1.8", 2000)]
        self.assertEqual(run([seg("A", 1, "B", 1, 1100, 2500)],
                             map_entries=short_map),
                         [out("A", 1, "B", 1, 1100, 2500, "3", "1.8")])

    def test_number_rounding(self):
        # Scores are float32 values rounded half to even: 0.045 is the float
        # 0.0450000017..., printed 0.05 (as a double it would print 0.04);
        # 0.125 is exact and prints 0.12 (half up would give 0.13).  The
        # length 0.0625 cM is exact and prints 0.062.
        tie_map = [("1", "0.0", 500), ("1", "0.0625", 600), ("1", "3.0", 3100)]
        self.assertEqual(run([seg("A", 1, "B", 1, 500, 600, "0.045"),
                              seg("A", 2, "B", 2, 500, 600, "0.125")],
                             map_entries=tie_map),
                         [out("A", 1, "B", 1, 500, 600, "0.05", "0.062"),
                          out("A", 2, "B", 2, 500, 600, "0.12", "0.062")])

    def test_different_chromosomes_not_merged(self):
        records = ["2\t%d\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0" % p
                   for p in POSITIONS]
        mp = MAP + [("2", "0.0", 500), ("2", "0.6", 1100), ("2", "3.0", 3100)]
        got = run([seg("A", 1, "B", 1, 500, 800), seg("A", 1, "B", 1, 1100, 1500, chrom="2")],
                  extra_records=records, map_entries=mp)
        self.assertEqual(got, [out("A", 1, "B", 1, 500, 800, "3", "0.3"),
                               out("A", 1, "B", 1, 1100, 1500, "3", "0.48", chrom="2")])


# ---------------------------------------------------------------------------
# Seeded random cases against a separate implementation of the rule
# ---------------------------------------------------------------------------

def oracle_merge(rows: List[Tuple], alleles: Dict[Tuple[str, int], Dict[str, Tuple[int, int]]],
                 gen_pos, gap: float, discord: int) -> List[str]:
    """The rule, implemented independently of merge.py: for each group
    (sample 1, haplotype 1, sample 2, haplotype 2, chromosome), keep a list
    of chains ordered by start and join the first neighbouring pair that
    qualifies, until no pair does."""
    groups = defaultdict(list)
    for s1, h1, s2, h2, chrom, start, end, score in rows:
        groups[(s1, h1, s2, h2, chrom)].append(
            [start, end, java_parse_float(score)])
    lines = []
    for (s1, h1, s2, h2, chrom), chains in groups.items():
        chains.sort(key=lambda c: (c[0], c[1]))
        joined = True
        while joined:
            joined = False
            for i in range(len(chains) - 1):
                left, right = chains[i], chains[i + 1]
                if right[0] < left[1]:
                    ok = True
                else:
                    length = gen_pos(chrom, right[0]) - gen_pos(chrom, left[1])
                    bad = sum(1 for p in alleles[chrom]
                              if left[1] <= p <= right[0]
                              and alleles[chrom][p][s1][h1 - 1]
                              != alleles[chrom][p][s2][h2 - 1])
                    ok = length <= gap and bad <= discord
                if ok:
                    chains[i:i + 2] = [[left[0], max(left[1], right[1]),
                                        max(left[2], right[2])]]
                    joined = True
                    break
        for start, end, score in chains:
            length = gen_pos(chrom, end) - gen_pos(chrom, start)
            lines.append("\t".join((s1, str(h1), s2, str(h2), chrom, str(start),
                                    str(end), java_decimal_format(score, 2),
                                    java_decimal_format(length, 3))))
    return lines


class RandomAgainstOracle(unittest.TestCase):

    def test_seeded_random_cases(self):
        for seed in range(1, 31):
            with self.subTest(seed=seed):
                self._check(seed)

    def _check(self, seed: int) -> None:
        rng = random.Random(seed)
        samples = ["r%d" % i for i in range(rng.randint(3, 7))]
        chroms = ["1", "2"]
        positions = {c: sorted(rng.sample(range(1000, 400000, 37), 150))
                     for c in chroms}
        alleles = {c: {p: {s: (rng.randint(0, 1), rng.randint(0, 1))
                           for s in samples} for p in positions[c]}
                   for c in chroms}
        # Long stretches where every sample carries the same alleles, so
        # that some gaps qualify.
        for c in chroms:
            for p in positions[c]:
                if (p // 20000) % 2 == 0 and rng.random() < 0.9:
                    same = (rng.randint(0, 1), rng.randint(0, 1))
                    for s in samples:
                        alleles[c][p][s] = same
        map_entries = []
        for c in chroms:
            cm = 0.0
            for p in sorted(rng.sample(positions[c], 20)) + [positions[c][-1] + 500]:
                cm += 0.0 if rng.random() < 0.1 else rng.uniform(0.01, 0.4)
                map_entries.append((c, repr(round(cm, 4)), p))
        rows = []
        for _ in range(rng.randint(20, 60)):
            i = rng.randrange(len(samples))
            j = rng.randrange(i, len(samples))
            h1, h2 = rng.randint(1, 2), rng.randint(1, 2)
            c = rng.choice(chroms)
            k = rng.randrange(len(positions[c]))
            for _ in range(rng.randint(1, 5)):
                e = min(k + rng.randint(0, 12), len(positions[c]) - 1)
                rows.append((samples[i], h1, samples[j], h2, c,
                             positions[c][k], positions[c][e],
                             "%.2f" % rng.uniform(0, 30)))
                k = max(0, min(e + rng.choice([-3, 0, 0, 1, 2, 4, 8]),
                               len(positions[c]) - 1))
        rng.shuffle(rows)
        gap = rng.choice([0.05, 0.2, 0.6, 1.0])
        discord = rng.choice([0, 1, 2])

        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, "in.vcf")
            with open(vcf, "w", encoding="utf-8", newline="\n") as f:
                f.write("##fileformat=VCFv4.2\n" + HEADER
                        + "".join("\t" + s for s in samples) + "\n")
                for c in chroms:
                    for p in positions[c]:
                        f.write("\t".join([c, str(p), ".", "A", "G", ".", "PASS",
                                           ".", "GT"]
                                          + ["%d|%d" % alleles[c][p][s]
                                             for s in samples]) + "\n")
            mp = os.path.join(d, "in.map")
            with open(mp, "w", encoding="utf-8", newline="\n") as f:
                for c, cm, bp in map_entries:
                    f.write("%s\t.\t%s\t%d\n" % (c, cm, bp))
            stdin = io.BytesIO("".join(
                "\t".join(str(x) for x in r) + "\n" for r in rows).encode())
            got = merge.run([vcf, mp, repr(gap), str(discord)], stdin).splitlines()
            # The output order does not depend on the input order.
            reordered = rows[:]
            random.Random(seed + 1000).shuffle(reordered)
            stdin = io.BytesIO("".join(
                "\t".join(str(x) for x in r) + "\n" for r in reordered).encode())
            again = merge.run([vcf, mp, repr(gap), str(discord)],
                              stdin).splitlines()
            self.assertEqual(again, got)
            chrom_ids = Indexer()
            gmap = PlinkGenMap.from_plink_map_file(mp, chrom_ids)

        def gen_pos(chrom, bp):
            return gmap.gen_pos(chrom_ids.get_index(chrom), bp)

        expected = oracle_merge(rows, {c: {p: alleles[c][p] for p in positions[c]}
                                       for c in chroms},
                                gen_pos, gap, discord)
        self.assertEqual(sorted(got), sorted(expected))
        self.assertLess(len(got), len(rows), "seed %d merged nothing" % seed)


if __name__ == "__main__":
    unittest.main()
