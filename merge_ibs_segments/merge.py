# merge-ibs-segments: a modified Python translation of merge-ibd-segments
# from Refined-IBD version 17Jan20.102.
#
# Copyright (C) 2014-2016 Brian L. Browning
# Copyright (C) 2026 LaKisha T. David
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
#
# Modified by LaKisha T. David, September 29-30, 2026: translated from Java to
# Python from src/ibdutil/MergeIbdSegments.java, src/blbutil/FileUtil.java,
# src/blbutil/Const.java, src/vcf/LowMafRefGT.java and
# src/vcf/LowMafRefDiallelicGT.java (including where it ends the program
# through Utilities.exit, which is translated in java_compat.py) in
# refined-ibd.17Jan20.102.zip (SHA-256
# 0b0abf48528ec53d3fec7f21a9742c6e5960857c5a225b0e49346179bc45e6ea). On
# September 30, 2026, Browning's merging rule was replaced by one that joins
# only segments on the same haplotype of each sample; on October 1, 2026, the
# joining conditions were changed, overlapping segments became an error, and
# three output columns and a report were added (see README).  Later changes
# are recorded in the Git history.
#
# Brian L. Browning's original notice, reproduced unchanged from those files:
#
# Copyright (C) 2014-2016 Brian L. Browning
#
# This file is part of Refined-IBD
#
# Refined-IBD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Beagle is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""The program: joins IBD segments of the same pair of samples that are on
the same haplotype of each sample.  Translated from ibdutil.MergeIbdSegments
(see the notice above for the other files it draws on), with its merging
rule replaced and three output columns and a report added.

Usage: ``cat [in] | python -m merge_ibs_segments [vcf] [map] [gap] [discord]
[report] > [out]``, where ``[report]`` is optional.

The IBD segments of each pair of samples are grouped by the haplotype of
sample 1 and the haplotype of sample 2 (1 or 2; haplotype 0 is rejected).
Within each group the segments are sorted as Browning's
StringIdSegment.compareTo sorts them and walked in order.  Two consecutive
segments of a group on the same chromosome are joined if

- they share an end marker;
- no marker strictly between them is discordant, whatever the length of the
  gap; or
- at most [discord] markers (0 or 1) strictly between them are discordant,
  and the gap, measured with the genetic map from the first segment's end
  marker to the second segment's start marker, is at most [gap] cM.

A marker is discordant when sample 1's allele on sample 1's haplotype
differs from sample 2's allele on sample 2's haplotype; the samples' other
haplotypes are not considered.  Segments of one group that overlap, or of
which one lies inside the other (including a one-marker segment at the
other's end marker, and repeated segments), stop the program with an error:
Refined IBD is not expected to report overlapping segments of one haplotype
pair, so an overlap means the input is not what this program assumes.
Every other pair of segments is left as it is and not reported.  A segment
of another haplotype group lying between them does not prevent a join.

Each chain of joined segments is written as one segment, with the group's
haplotypes, the largest score in the chain and its length from the genetic
map, followed by three columns: the number of VCF markers from its start to
its end, the number of haplotypes in the VCF file carrying the segment's
alleles (at every marker of the segment where the two haplotypes agree), and
that number divided by all haplotypes in the VCF file.

The report, if requested, lists every join and these flags: a join at a
shared end marker, a join across a gap with no discordant marker, and a
segment whose neighbouring marker beyond its start or end is not discordant
(Refined IBD ends a segment next to a discordant marker, except at the edge
of one of its analysis windows).  A summary goes to standard error.

Differences from Browning's merge-ibd-segments, of which the tag
baseline-faithful is an exact translation:

- Merging rule.  Browning's program joins segments of a pair whatever their
  haplotypes; joins overlapping and touching segments without a test; joins
  across a gap only if it is at most [gap] cM and has at most [discord]
  discordant markers, counting the two end markers; counts a marker as
  discordant only when no allele of sample 1 equals either allele of
  sample 2 (the genotypes share no allele); accepts haplotype 0 in its input;
  and writes haplotype 0 for every joined segment.
- Every segment's start and end must be the position of exactly one VCF
  record, both samples must be in the VCF file, each chromosome's VCF
  records must form one block in position order, and allele indices must be
  at most 255 (alleles are held as bytes).  Browning's program checks
  positions only in gap tests, does not check the record order, and reads
  outside its genotype arrays for a sample that is not in the VCF file.
- [report], if given, is checked before any input is read: not empty, not a
  directory, in an existing directory, writable if it exists, and not the
  [vcf] or [map] file.
- Output lines have three more columns, and there is a report and a summary.
- The usage message describes this program.  It goes to standard output
  only when there are no arguments; with an argument error it goes to
  standard error with the error message, exit status 1.  Browning's program
  writes it to standard output, with exit status 0 for a wrong number of
  arguments.
- Output lines are sorted by chromosome (in order of first appearance in the
  genetic map), start, end, score, sample 1, sample 2, haplotype 1 and
  haplotype 2.  The Java program writes pairs in java.util.HashMap order.
- When the program stops with an error, no segments are written and no
  report is written.  The Java program may already have written some
  segments before stopping.
- Java stack trace frame lines are not written.
- Output lines end with ``\\n``.  The Java program ends them with the
  platform's line separator.

Known difference found by the baseline's exact-conversion test
(tests/parity), outside the range of real scores and lengths: some numbers
of 2**53 (about 9.0e15) or more print with other digits than the jar's
(java_compat.java_decimal_format).
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
from decimal import ROUND_HALF_EVEN, Decimal, localcontext
from typing import BinaryIO, Dict, List, NamedTuple, Optional, Set, Tuple

from .genetic_map import PlinkGenMap
from .ids import Indexer
from .input_it import from_gzip_file, from_std_in
from .java_compat import (DOUBLE_MAX, JavaException, JavaTermination,
                          illegal_argument, java_decimal_format, java_encode,
                          java_compare_key, java_parse_double,
                          java_parse_int, java_trim, to_java_chars,
                          utilities_exit)
from .segments import Segment, parse_segment
from .vcf import VcfHeader, read_vcf

PROGRAM = "python -m merge_ibs_segments"

# Report events.  The second element says whether the event is flagged.
JOINED_ONE_DISCORDANT = "joined_one_discordant"
JOINED_NO_DISCORDANT = "joined_no_discordant"
JOINED_NO_DISCORDANT_LONG_GAP = "joined_no_discordant_long_gap"
JOINED_SHARED_END = "joined_shared_end"
START_MATCHES_PREVIOUS = "start_matches_previous_snp"
END_MATCHES_NEXT = "end_matches_next_snp"
EVENTS = ((JOINED_ONE_DISCORDANT, False), (JOINED_NO_DISCORDANT, True),
          (JOINED_NO_DISCORDANT_LONG_GAP, True), (JOINED_SHARED_END, True),
          (START_MATCHES_PREVIOUS, True), (END_MATCHES_NEXT, True))
FLAGGED = dict(EVENTS)

REPORT_COLUMNS = ("event", "flagged", "sample1", "hap1", "sample2", "hap2",
                  "chromosome", "start1", "end1", "start2", "end2",
                  "snps_between", "discordant_between", "gap_cm")


def usage() -> str:
    """The usage message, laid out as ``MergeIbdSegments.usage()``."""
    nl = "\n"
    return (nl
            + "usage: cat [in] | " + PROGRAM
            + " [vcf] [map] [gap] [discord] [report] > [out]" + nl
            + nl
            + "where" + nl
            + "  [in]      = IBD output file from Refined IBD analysis "
              "(uncompressed)" + nl
            + "  [vcf]     = Phased input VCF file from Refined IBD analysis"
            + nl
            + "  [map]     = PLINK format genetic map file with centimorgan "
              "(cM) distances" + nl
            + "  [gap]     = max length (cM) of a gap joined across a "
              "discordant marker;" + nl
            + "              a gap with no discordant marker is joined at any "
              "length" + nl
            + "  [discord] = max number of discordant markers in a gap (0 or "
              "1)" + nl
            + "  [report]  = optional file listing joins and flagged segments"
            + nl
            + "  [out]     = IBD output file with IBD segments after merging"
            + nl
            + nl
            + "IBD segments of a pair of samples are merged only if they are "
              "on the same" + nl
            + "haplotype of each sample and the same chromosome.  Two such "
              "segments are" + nl
            + "merged if they share an end marker, if no marker between them "
              "is discordant" + nl
            + "(the two haplotypes carry different alleles), or if at most "
              "[discord] markers" + nl
            + "between them are discordant and the gap is at most [gap] cM.  "
              "Overlapping" + nl
            + "or nested segments of the same samples and haplotypes stop "
              "the program with an" + nl
            + "error.  Merged segments keep their haplotype indices and have "
              "IBD score equal" + nl
            + "to the maximal score of the merged segments." + nl
            + "Each output line adds the number of markers in the segment, "
              "the number of" + nl
            + "haplotypes in [vcf] carrying the segment's alleles, and that "
              "number divided" + nl
            + "by all haplotypes in [vcf]." + nl)


def exit_with_error(message: str) -> JavaTermination:
    """``MergeIbdSegments.exitWithError``: the usage message, then the
    error message, exit status 1.  Browning's program writes the usage
    message to standard output; this program writes both to standard error,
    so that the usage text never lands in [out]."""
    return JavaTermination(stderr=usage() + "\n" + message + "\n", status=1)


def parse_and_check_max_gap(arg: str) -> float:
    """``MergeIbdSegments.parseAndCheckMaxGap``."""
    try:
        max_gap = java_parse_double(arg)
    except JavaException as error:
        if error.java_class == "java.lang.NumberFormatException":
            raise exit_with_error("ERROR: [gap] is not a numerical value: "
                                  + arg)
        raise
    if not math.isfinite(max_gap) or max_gap <= 0.0:
        raise exit_with_error("ERROR: [gap]=" + arg)
    return max_gap


def parse_and_check_max_discord(arg: str) -> int:
    """``MergeIbdSegments.parseAndCheckMaxDiscord``, limited to 0 and 1: a
    gap is joined across at most one discordant marker."""
    try:
        max_discord = java_parse_int(arg)
    except JavaException as error:
        if error.java_class == "java.lang.NumberFormatException":
            raise exit_with_error(
                "ERROR: [discord] is not a nonnegative integer: " + arg)
        raise
    if max_discord < 0:
        raise exit_with_error("ERROR: [discord] < 0: " + arg)
    if max_discord > 1:
        raise exit_with_error("ERROR: [discord] must be 0 or 1: " + arg)
    return max_discord


def _position_key(chrom_index: int, pos: int) -> int:
    return (chrom_index << 32) | (pos + (1 << 31))


def format_frequency(carriers: int, n_haplotypes: int) -> str:
    """``carriers / n_haplotypes`` rounded half to even to 6 significant
    digits, in plain decimal notation (``0.0208333``, ``0.5``, ``1``)."""
    with localcontext() as ctx:
        ctx.prec = 6
        ctx.rounding = ROUND_HALF_EVEN
        value = Decimal(carriers) / Decimal(n_haplotypes)
    return format(value, "f")


class Data:
    """``MergeIbdSegments.Data``: the genetic map, the VCF samples and the
    VCF records, with an index from (chromosome, position) to record and,
    built when first needed, each haplotype's alleles as bytes."""

    def __init__(self, vcf_path: str, map_path: str, chrom_ids: Indexer,
                 sample_ids: Indexer):
        self.chrom_ids = chrom_ids
        self.sample_ids = sample_ids
        self.gen_map = PlinkGenMap.from_plink_map_file(map_path, chrom_ids)
        lines = from_gzip_file(vcf_path)
        self.header, self.recs = read_vcf(vcf_path, chrom_ids, sample_ids,
                                          lines)
        # Segment ranges, the markers between segments and the neighbouring
        # markers are taken by record index, so each chromosome's records
        # must form one block in position order.
        index_map: Dict[int, int] = {}
        repeated: Set[int] = set()
        seen_chroms: Set[int] = set()
        prev = None
        for marker_index, rec in enumerate(self.recs):
            if prev is None or rec.chrom_index != prev.chrom_index:
                if rec.chrom_index in seen_chroms:
                    raise utilities_exit(
                        "ERROR: the records of chromosome "
                        + chrom_ids.id(rec.chrom_index)
                        + " are not in one block in the VCF file")
                seen_chroms.add(rec.chrom_index)
            elif rec.pos < prev.pos:
                raise utilities_exit(
                    "ERROR: VCF records are not in position order: "
                    "chromosome " + chrom_ids.id(rec.chrom_index)
                    + " position " + str(rec.pos) + " follows position "
                    + str(prev.pos))
            prev = rec
            key = _position_key(rec.chrom_index, rec.pos)
            if key in index_map:
                repeated.add(key)
            index_map[key] = marker_index
        self._index_map = index_map
        self._repeated = repeated
        self._haps: Optional[List[bytearray]] = None

    def describe(self, seg: Segment) -> str:
        """A segment as an error message names it."""
        return (self.sample_ids.id(seg.sample_index1) + " haplotype "
                + str(seg.hap1) + ", " + self.sample_ids.id(seg.sample_index2)
                + " haplotype " + str(seg.hap2) + ", chromosome "
                + self.chrom_ids.id(seg.chrom_index) + ", " + str(seg.start)
                + "-" + str(seg.end))

    def end_marker_index(self, seg: Segment, pos: int) -> int:
        """The index of the VCF record at a segment's start or end.  Stops
        the program if there is no record at that position or more than
        one."""
        key = _position_key(seg.chrom_index, pos)
        index = self._index_map.get(key)
        if index is None:
            raise utilities_exit("ERROR: segment end is not a VCF marker: "
                                 + "position " + str(pos) + " of "
                                 + self.describe(seg))
        if key in self._repeated:
            raise utilities_exit("ERROR: more than one VCF record at a "
                                 "segment end: position " + str(pos)
                                 + " of " + self.describe(seg))
        return index

    def haplotype(self, sample: int, hap: int) -> int:
        """The index of a sample's haplotype (1 or 2) in ``haplotypes()``.
        ``read_segments`` has checked that the sample is in the VCF file."""
        return 2 * self.header.sample_index(sample) + hap - 1

    def haplotypes(self) -> List[bytearray]:
        """Each haplotype's allele indices, one byte per VCF record:
        haplotype ``2s`` holds the first and ``2s + 1`` the second allele of
        the phased genotype of VCF sample ``s``."""
        if self._haps is None:
            n_samples = self.header.n_samples
            haps = [bytearray(len(self.recs)) for _ in range(2 * n_samples)]
            for m, rec in enumerate(self.recs):
                allele1 = rec.allele1
                allele2 = rec.allele2
                try:
                    for s in range(n_samples):
                        haps[2 * s][m] = allele1[s]
                        haps[2 * s + 1][m] = allele2[s]
                except ValueError:
                    raise utilities_exit(
                        "ERROR: allele index above 255 at chromosome "
                        + self.chrom_ids.id(rec.chrom_index) + " position "
                        + str(rec.pos))
            self._haps = haps
        return self._haps

    def gen_distance(self, start_index: int, end_index: int) -> float:
        """``Data.genDistance``."""
        if start_index > end_index:
            raise illegal_argument("startIndex=" + str(start_index)
                                   + " endIndex=" + str(end_index))
        start_marker = self.recs[start_index]
        end_marker = self.recs[end_index]
        if start_marker.chrom_index != end_marker.chrom_index:
            raise illegal_argument("inconsistent chromosomes: startIndex="
                                   + str(start_index) + " endIndex="
                                   + str(end_index))
        end_pos = self.gen_map.gen_pos(end_marker.chrom_index, end_marker.pos)
        start_pos = self.gen_map.gen_pos(start_marker.chrom_index,
                                         start_marker.pos)
        return end_pos - start_pos

    def alleles_equal(self, hap1: int, hap2: int, marker_index: int) -> bool:
        """Whether two haplotypes carry the same allele at a marker."""
        haps = self.haplotypes()
        return haps[hap1][marker_index] == haps[hap2][marker_index]

    def discordant_between(self, hap1: int, hap2: int, start_index: int,
                           end_index: int, limit: int) -> int:
        """The markers strictly between ``start_index`` and ``end_index`` at
        which the two haplotypes carry different alleles, counted up to
        ``limit``.  Replaces Browning's ``Data.ibdDiscordCnt``, which counts
        both end markers and counts a marker only when no allele of sample 1
        equals either allele of sample 2."""
        haps = self.haplotypes()
        a = haps[hap1]
        b = haps[hap2]
        if a[start_index + 1:end_index] == b[start_index + 1:end_index]:
            return 0
        count = 0
        for m in range(start_index + 1, end_index):
            if a[m] != b[m]:
                count += 1
                if count == limit:
                    break
        return count

    def carriers(self, hap1: int, hap2: int, start_index: int,
                 end_index: int) -> int:
        """The haplotypes in the VCF file carrying the segment's alleles:
        those equal to ``hap1`` at every marker from ``start_index`` to
        ``end_index`` where ``hap1`` and ``hap2`` agree.  Markers where they
        differ (a discordant marker inside a joined segment) are skipped.
        The two haplotypes themselves are counted."""
        haps = self.haplotypes()
        ref = haps[hap1]
        other = haps[hap2]
        stop = end_index + 1
        if ref[start_index:stop] == other[start_index:stop]:
            pieces = [(start_index, stop)]
        else:
            pieces = []
            piece_start = start_index
            for m in range(start_index, stop):
                if ref[m] != other[m]:
                    if piece_start < m:
                        pieces.append((piece_start, m))
                    piece_start = m + 1
            if piece_start < stop:
                pieces.append((piece_start, stop))
        wanted = [(a, b, ref[a:b]) for a, b in pieces]
        return sum(1 for hap in haps
                   if all(hap[a:b] == alleles for a, b, alleles in wanted))


class Placed(NamedTuple):
    """A segment with the VCF record indices of its start and end."""
    segment: Segment
    start_index: int
    end_index: int


class Result(NamedTuple):
    """What a successful run produces, and where the report goes (None if
    no report file was named)."""
    output: str
    report: str
    summary: str
    report_path: Optional[str]


def read_segments(lines, header: VcfHeader, sample_ids: Indexer,
                  chrom_ids: Indexer) -> Dict[Tuple[int, int, int, int],
                                              List[Segment]]:
    """``MergeIbdSegments.readSegments``, with segments grouped by sample
    pair and by the haplotype of each sample, in order of first appearance.
    Stops the program if a sample is not in the VCF file.  This is checked
    first: Browning's sample-order check would otherwise report a missing
    sample 2 as a sample-order error."""
    pairs: Dict[Tuple[int, int, int, int], List[Segment]] = {}
    for raw in lines:
        line = java_trim(raw)
        if len(line) > 0:
            segment = parse_segment(line, sample_ids, chrom_ids)
            s1 = segment.sample_index1
            s2 = segment.sample_index2
            for sample in (s1, s2):
                if header.sample_index(sample) == -1:
                    raise utilities_exit("ERROR: sample not in VCF file: "
                                         + sample_ids.id(sample))
            if header.sample_index(s1) > header.sample_index(s2):
                raise illegal_argument("Inconsistency in sample order:" + line)
            key = (s1, s2, segment.hap1, segment.hap2)
            group = pairs.get(key)
            if group is None:
                group = []
                pairs[key] = group
            group.append(segment)
    return pairs


def join_event(data: Data, last: Placed, nxt: Placed, hap1: int, hap2: int,
               max_gap: float,
               max_discord: int) -> Tuple[Optional[str], Tuple[str, ...]]:
    """Whether ``nxt`` joins a chain whose last segment is ``last`` (all of
    one haplotype group, sorted).  Returns the report event, or None if it
    does not join, and the report's gap fields (markers strictly between,
    discordant markers between, gap in cM).  Stops the program if the two
    segments overlap or one lies inside the other, including a one-marker
    segment at the other's end marker; only two segments that each extend
    beyond a shared end marker share an end."""
    a = last.segment
    b = nxt.segment
    if a.chrom_index != b.chrom_index:
        return None, ()
    if b.start < a.end or (b.start == a.end
                           and (a.start == a.end or b.start == b.end)):
        raise utilities_exit(
            "ERROR: segments of the same samples and haplotypes overlap: "
            + data.describe(a) + " and " + data.describe(b) + ".  Refined "
            "IBD is not expected to report overlapping segments of one "
            "haplotype pair, so this input is not what the program assumes.")
    if b.start == a.end:
        return JOINED_SHARED_END, ("0", "0", "0")
    discord = data.discordant_between(hap1, hap2, last.end_index,
                                      nxt.start_index, limit=2)
    gap_length = data.gen_distance(last.end_index, nxt.start_index)
    fields = (str(nxt.start_index - last.end_index - 1), str(discord),
              java_decimal_format(gap_length, 3))
    if discord == 0:
        if gap_length <= max_gap:
            return JOINED_NO_DISCORDANT, fields
        return JOINED_NO_DISCORDANT_LONG_GAP, fields
    if discord <= max_discord and gap_length <= max_gap:
        return JOINED_ONE_DISCORDANT, fields
    return None, ()


def end_flags(data: Data, placed: Placed, hap1: int,
              hap2: int) -> List[str]:
    """The flags for a segment whose neighbouring marker beyond its start or
    end, on the same chromosome, carries the same allele on both
    haplotypes."""
    flags = []
    seg = placed.segment
    before = placed.start_index - 1
    if (before >= 0 and data.recs[before].chrom_index == seg.chrom_index
            and data.alleles_equal(hap1, hap2, before)):
        flags.append(START_MATCHES_PREVIOUS)
    after = placed.end_index + 1
    if (after < len(data.recs)
            and data.recs[after].chrom_index == seg.chrom_index
            and data.alleles_equal(hap1, hap2, after)):
        flags.append(END_MATCHES_NEXT)
    return flags


def report_row(data: Data, event: str, first: Segment,
               second: Optional[Segment],
               gap_fields: Tuple[str, ...]) -> Tuple[Tuple, str]:
    """A report line, with its sort key."""
    if second is None:
        second_fields = (".", ".")
        gap_fields = (".", ".", ".")
    else:
        second_fields = (str(second.start), str(second.end))
    fields = ((event, "yes" if FLAGGED[event] else "no",
               data.sample_ids.id(first.sample_index1), str(first.hap1),
               data.sample_ids.id(first.sample_index2), str(first.hap2),
               data.chrom_ids.id(first.chrom_index), str(first.start),
               str(first.end)) + second_fields + gap_fields)
    key = (first.chrom_index, first.start, first.end,
           -1 if second is None else second.start, first.sample_index1,
           first.sample_index2, first.hap1, first.hap2, event)
    return key, "\t".join(fields) + "\n"


def write_chain(data: Data, chain: List[Placed], hap1: int, hap2: int,
                n_haplotypes: int) -> Tuple[Tuple, str]:
    """``MergeIbdSegments.printAndClearList``: the output line for a chain
    (with its sort key).  The chain keeps its group's haplotypes; Browning's
    program writes 0 for a joined segment.  Three columns follow Browning's
    nine: the number of VCF markers in the segment, the haplotypes carrying
    its alleles, and their frequency."""
    first = chain[0].segment
    score = -DOUBLE_MAX
    for placed in chain:
        if placed.segment.score > score:
            score = placed.segment.score
    start = first.start
    end = chain[-1].segment.end
    start_index = chain[0].start_index
    end_index = chain[-1].end_index
    gen_length = (data.gen_map.gen_pos(first.chrom_index, end)
                  - data.gen_map.gen_pos(first.chrom_index, start))
    carriers = data.carriers(hap1, hap2, start_index, end_index)
    line = "\t".join((data.sample_ids.id(first.sample_index1),
                      str(first.hap1),
                      data.sample_ids.id(first.sample_index2),
                      str(first.hap2), data.chrom_ids.id(first.chrom_index),
                      str(start), str(end), java_decimal_format(score, 2),
                      java_decimal_format(gen_length, 3),
                      str(end_index - start_index + 1), str(carriers),
                      format_frequency(carriers, n_haplotypes))) + "\n"
    key = (first.chrom_index, start, end, java_compare_key(score),
           first.sample_index1, first.sample_index2, first.hap1, first.hap2,
           line)
    return key, line


def summary_text(n_read: int, n_written: int, counts: Dict[str, int]) -> str:
    """The summary written to standard error."""
    return ("merge-ibs-segments: " + str(n_read) + " segments read, "
            + str(n_written) + " written; joined: "
            + str(counts[JOINED_ONE_DISCORDANT])
            + " across one discordant marker, "
            + str(counts[JOINED_NO_DISCORDANT]
                  + counts[JOINED_NO_DISCORDANT_LONG_GAP])
            + " across a gap with no discordant marker ("
            + str(counts[JOINED_NO_DISCORDANT_LONG_GAP])
            + " of them longer than [gap]), "
            + str(counts[JOINED_SHARED_END]) + " at a shared end marker; "
            + "segment ends whose neighbouring marker matches: "
            + str(counts[START_MATCHES_PREVIOUS] + counts[END_MATCHES_NEXT])
            + "\n")


def run(args: List[str], stdin: BinaryIO) -> Result:
    """Runs the program and returns the output, the report and the summary;
    ``write_outputs`` writes them.  Raises JavaTermination or JavaException
    where the program stops.  With no arguments the usage message goes to
    standard output (exit status 0); with another wrong number of arguments
    it is an error.  (Browning's program writes the usage message to
    standard output with status 0 for any wrong number.)"""
    if not args:
        raise JavaTermination(stdout=usage() + "\n", status=0)
    if len(args) not in (4, 5):
        raise exit_with_error("ERROR: expected 4 or 5 arguments, found "
                              + str(len(args)))
    vcf_path, map_path = args[0], args[1]
    max_gap = parse_and_check_max_gap(args[2])
    max_discord = parse_and_check_max_discord(args[3])
    report_path = (check_report_path(args[4], vcf_path, map_path)
                   if len(args) == 5 else None)
    chrom_ids = Indexer()
    sample_ids = Indexer()
    data = Data(vcf_path, map_path, chrom_ids, sample_ids)
    pairs = read_segments(from_std_in(stdin), data.header, sample_ids,
                          chrom_ids)
    n_haplotypes = 2 * data.header.n_samples
    counts = {event: 0 for event, _ in EVENTS}
    out: List[Tuple[Tuple, str]] = []
    report: List[Tuple[Tuple, str]] = []
    n_read = 0
    for (s1, s2, h1, h2), group in pairs.items():
        hap1 = data.haplotype(s1, h1)
        hap2 = data.haplotype(s2, h2)
        group.sort(key=Segment.sort_key)
        n_read += len(group)
        chain: List[Placed] = []
        for seg in group:
            placed = Placed(seg, data.end_marker_index(seg, seg.start),
                            data.end_marker_index(seg, seg.end))
            for flag in end_flags(data, placed, hap1, hap2):
                counts[flag] += 1
                report.append(report_row(data, flag, seg, None, ()))
            if chain:
                last = chain[-1]
                event, gap_fields = join_event(data, last, placed, hap1,
                                               hap2, max_gap, max_discord)
                if event is None:
                    out.append(write_chain(data, chain, hap1, hap2,
                                           n_haplotypes))
                    chain = []
                else:
                    counts[event] += 1
                    report.append(report_row(data, event, last.segment, seg,
                                             gap_fields))
            chain.append(placed)
        out.append(write_chain(data, chain, hap1, hap2, n_haplotypes))
    out.sort(key=lambda item: item[0])
    report.sort(key=lambda item: item[0])
    report_text = ("\t".join(REPORT_COLUMNS) + "\n"
                   + "".join(line for _, line in report))
    return Result("".join(line for _, line in out), report_text,
                  summary_text(n_read, len(out), counts), report_path)


def python_text(s: str) -> str:
    """A string from ``to_java_chars`` (characters beyond the Basic
    Multilingual Plane held as UTF-16 surrogate pairs) as an ordinary Python
    string, for use as a file name."""
    try:
        return s.encode("utf-16-le", "surrogatepass").decode("utf-16-le")
    except UnicodeDecodeError:
        return s


def check_report_path(arg: str, vcf_path: str, map_path: str) -> str:
    """The file the report goes to: ``[report]`` with symbolic links
    followed.  Checked before any input is read, so that a report path that
    cannot be used stops the program before anything is written: it must
    name a file, not a directory, in an existing directory; an existing file
    there must be writable; and it must not be ``[vcf]`` or ``[map]``."""
    if not arg:
        raise exit_with_error("ERROR: [report] is empty")
    path = os.path.realpath(python_text(arg))
    if os.path.isdir(path):
        raise exit_with_error("ERROR: [report] is a directory: " + arg)
    if not os.path.isdir(os.path.dirname(path)):
        raise exit_with_error("ERROR: the directory of [report] does not "
                              "exist: " + arg)
    if os.path.exists(path) and not os.access(path, os.W_OK):
        raise exit_with_error("ERROR: [report] is not writable: " + arg)
    for name, other in (("[vcf]", vcf_path), ("[map]", map_path)):
        if (os.path.normcase(path)
                == os.path.normcase(os.path.realpath(python_text(other)))):
            raise exit_with_error("ERROR: [report] is the same file as "
                                  + name + ": " + arg)
    return path


def write_outputs(result: Result, stdout: BinaryIO) -> None:
    """Writes the segments to ``stdout`` and, if a report file was named,
    the report.  The report is written first to a temporary file in the
    same directory, which replaces the named file only after the segments
    have been written and is removed if anything fails; so a failed run
    leaves no new or partial report, and leaves an earlier report at that
    path unchanged.  The one exception: if the report cannot be moved into
    place at the end (for example because another program holds the file
    open), the segments have already been written and the program stops
    with an error.  A run that is killed can leave the temporary file
    (``.merge-ibs-report-*``) behind."""
    temp_path = None
    try:
        if result.report_path is not None:
            directory = os.path.dirname(result.report_path)
            try:
                fd, temp_path = tempfile.mkstemp(prefix=".merge-ibs-report-",
                                                 dir=directory)
                with os.fdopen(fd, "wb") as f:
                    f.write(java_encode(result.report))
                # mkstemp creates the file readable by its owner only; give
                # it the permissions a newly created file would get.
                mask = os.umask(0)
                os.umask(mask)
                os.chmod(temp_path, 0o666 & ~mask)
            except OSError as error:
                raise utilities_exit("ERROR: cannot write [report] file "
                                     + result.report_path + ": "
                                     + str(error))
        # FileUtil.stdOutPrintWriter: segments go to standard output through
        # a writer in the default charset (UTF-8 in Java 18 and later).
        try:
            stdout.write(java_encode(result.output))
            stdout.flush()
        except OSError as error:
            termination = utilities_exit(
                "ERROR: cannot write the merged segments: " + str(error))
            termination.stdout_failed = True
            raise termination
        if temp_path is not None:
            try:
                os.replace(temp_path, result.report_path)
            except OSError as error:
                raise utilities_exit("ERROR: cannot write [report] file "
                                     + result.report_path + ": "
                                     + str(error))
            temp_path = None
    finally:
        if temp_path is not None:
            try:
                os.remove(temp_path)
            except OSError:
                pass


def main(argv: Optional[List[str]] = None) -> int:
    """Command-line entry point; returns the exit status."""
    args = [to_java_chars(a) for a in (sys.argv[1:] if argv is None
                                       else argv)]
    stdout = sys.stdout.buffer
    stderr = sys.stderr.buffer
    try:
        result = run(args, sys.stdin.buffer)
        write_outputs(result, stdout)
    except JavaTermination as termination:
        if getattr(termination, "stdout_failed", False):
            # Standard output is closed or failing (for example a pipe whose
            # reader has gone).  Point it at the null device, Python's
            # documented remedy, so that the interpreter's final flush of the
            # unwritten bytes cannot fail again; then report the error.
            try:
                devnull = os.open(os.devnull, os.O_WRONLY)
                os.dup2(devnull, sys.stdout.fileno())
            except (OSError, ValueError):
                pass
        elif termination.stdout:
            stdout.write(java_encode(termination.stdout))
            stdout.flush()
        stderr.write(java_encode(termination.stderr))
        stderr.flush()
        return termination.status
    except JavaException as error:
        stderr.write(java_encode('Exception in thread "main" '
                                 + error.to_string() + "\n"))
        stderr.flush()
        return 1
    stderr.write(java_encode(result.summary))
    stderr.flush()
    return 0
