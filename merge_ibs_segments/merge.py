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
# only segments on the same haplotype of each sample (see README). Later
# changes are recorded in the Git history.
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
rule replaced.

Usage: ``cat [in] | python -m merge_ibs_segments [vcf] [map] [gap] [discord]
> [out]``.

The IBD segments of each pair of samples are grouped by the haplotype of
sample 1 and the haplotype of sample 2 (1 or 2; haplotype 0 is rejected).
Within each group the segments are sorted as Browning's
StringIdSegment.compareTo sorts them and walked in order.  A segment joins
the current chain if it is on the same chromosome and either

- starts before the chain's running end (an overlapping or nested segment of
  the same haplotype pair), or
- starts at or after the running end, and the gap, measured with the genetic
  map from the marker at the running end to the marker at the segment's start
  (zero when they are the same position), is at most [gap] cM and contains at
  most [discord] discordant markers, both end markers included.  A marker is
  discordant when sample 1's allele on sample 1's haplotype differs from
  sample 2's allele on sample 2's haplotype; the samples' other haplotypes are
  not considered.  A segment that starts at the running end is always tested
  this way, even a one-marker segment lying inside the chain (which a minimum
  segment length of more than zero rules out in practice).

A segment of another haplotype group lying in the gap does not prevent a
join.  Each chain is written as one segment, with the group's haplotypes, the
largest score in the chain, and its length recomputed from the genetic map.

Differences from Browning's merge-ibd-segments, of which the tag
baseline-faithful is an exact translation:

- Merging rule.  Browning's program joins segments of a pair whatever their
  haplotypes; joins a segment that starts at the running end without a gap
  test; counts a marker as discordant only when no allele of sample 1 equals
  either allele of sample 2 (the genotypes share no allele); accepts
  haplotype 0 in its input; and writes haplotype 0 for every joined segment.
- The usage message describes this program's rule.
- Output lines are sorted by chromosome (in order of first appearance in the
  genetic map), start, end, score, sample 1, sample 2, haplotype 1 and
  haplotype 2.  The Java program writes pairs in java.util.HashMap order.
- When the program stops with an error, no segments are written.  The Java
  program may already have written some segments before stopping.
- Java stack trace frame lines are not written.
- Output lines end with ``\\n``.  The Java program ends them with the
  platform's line separator.

Known differences found by the baseline's exact-conversion test
(tests/parity), outside the range of real scores and lengths:

- Some numbers of 2**53 (about 9.0e15) or more print with other digits than
  the jar's (java_compat.java_decimal_format).
- A gap test for a sample that is in the IBD input but not in the VCF file:
  depending on how the jar stored the genotypes, it stops at a bounds check or
  reads one element past its genotype array.  This program stops with
  IndexOutOfBoundsException, naming the index of the haplotype it reads
  (2 x sample + haplotype - 1; the baseline named 2 x sample).
"""

from __future__ import annotations

import math
import sys
from typing import BinaryIO, Dict, List, Optional, Tuple

from .genetic_map import PlinkGenMap
from .ids import Indexer
from .input_it import from_gzip_file, from_std_in
from .java_compat import (DOUBLE_MAX, INT_MAX, INT_MIN, JavaException,
                          JavaTermination, illegal_argument,
                          java_decimal_format, java_encode,
                          java_compare_key, java_parse_double,
                          java_parse_int, java_trim, to_java_chars,
                          utilities_exit)
from .segments import Segment, parse_segment
from .vcf import VcfHeader, VcfRecord, read_vcf

PROGRAM = "python -m merge_ibs_segments"


def usage() -> str:
    """The usage message, laid out as ``MergeIbdSegments.usage()``."""
    nl = "\n"
    return (nl
            + "usage: cat [in] | " + PROGRAM
            + " [vcf] [map] [gap] [discord] > [out]" + nl
            + nl
            + "where" + nl
            + "  [in]      = IBD output file from Refined IBD analysis "
              "(uncompressed)" + nl
            + "  [vcf]     = Phased input VCF file from Refined IBD analysis"
            + nl
            + "  [map]     = PLINK format genetic map file with centimorgan "
              "(cM) distances" + nl
            + "  [gap]     = max length of gap between IBD segments (cM)" + nl
            + "  [discord] = max number of markers in the gap at which the "
              "two IBD haplotypes" + nl
            + "              carry different alleles" + nl
            + "  [out]     = IBD output file with IBD segments after merging"
            + nl
            + nl
            + "IBD segments for a pair of samples are merged only if they are "
              "on the same" + nl
            + "haplotype of each sample.  Such segments are merged if the "
              "later one starts" + nl
            + "before the earlier one ends, or if they are separated by a gap "
              "(or share an end" + nl
            + "position) having length <= [gap] cM and having <= [discord] "
              "markers at which" + nl
            + "the two haplotypes carry different alleles.  Merged segments "
              "keep their" + nl
            + "haplotype indices and have IBD score equal to the maximal "
              "score of the merged" + nl
            + "segments." + nl)


def exit_with_error(message: str) -> JavaTermination:
    """``MergeIbdSegments.exitWithError``: the usage message on standard
    output, then ``Utilities.exit(message)``."""
    termination = utilities_exit(message)
    termination.stdout = usage() + "\n"
    return termination


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
    """``MergeIbdSegments.parseAndCheckMaxDiscord``."""
    try:
        max_discord = java_parse_int(arg)
    except JavaException as error:
        if error.java_class == "java.lang.NumberFormatException":
            raise exit_with_error(
                "ERROR: [discord] is not a nonnegative integer: " + arg)
        raise
    if max_discord < 0:
        raise exit_with_error("ERROR: [discord] < 0: " + arg)
    return max_discord


def _position_key(chrom_index: int, pos: int) -> int:
    return (chrom_index << 32) | (pos + (1 << 31))


class Data:
    """``MergeIbdSegments.Data``: the genetic map, the VCF samples and the
    VCF records, with an index from (chromosome, position) to record."""

    def __init__(self, vcf_path: str, map_path: str, chrom_ids: Indexer,
                 sample_ids: Indexer):
        self.gen_map = PlinkGenMap.from_plink_map_file(map_path, chrom_ids)
        lines = from_gzip_file(vcf_path)
        self.header, self.recs = read_vcf(vcf_path, chrom_ids, sample_ids,
                                          lines)
        index_map: Dict[int, int] = {}
        for marker_index, rec in enumerate(self.recs):
            index_map[_position_key(rec.chrom_index, rec.pos)] = marker_index
        self._index_map = index_map

    def marker_index(self, chrom: int, pos: int) -> int:
        """The index of the last record at this chromosome and position, or
        -1."""
        return self._index_map.get(_position_key(chrom, pos), -1)

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

    def haplotype_discord_cnt(self, sample1: int, hap1: int, sample2: int,
                              hap2: int, start_index: int,
                              end_index: int) -> int:
        """The markers from ``start_index`` to ``end_index`` inclusive at
        which sample 1's allele on haplotype ``hap1`` differs from sample 2's
        allele on haplotype ``hap2`` (haplotypes 1 and 2 are the first and
        second alleles of the phased genotype).  Replaces Browning's
        ``Data.ibdDiscordCnt``, which counts a marker only when no allele of
        sample 1 equals either allele of sample 2."""
        if start_index > end_index:
            raise illegal_argument("startIndex=" + str(start_index)
                                   + " endIndex=" + str(end_index))
        n_samples = self.header.n_samples
        for sample, hap in ((sample1, hap1), (sample2, hap2)):
            if sample >= n_samples:
                # A sample that is not in the VCF file has no alleles.
                raise JavaException("java.lang.IndexOutOfBoundsException",
                                    str(2 * sample + hap - 1))
        discord_cnt = 0
        for m in range(start_index, end_index + 1):
            rec = self.recs[m]
            a = (rec.allele1 if hap1 == 1 else rec.allele2)[sample1]
            b = (rec.allele1 if hap2 == 1 else rec.allele2)[sample2]
            if a != b:
                discord_cnt += 1
        return discord_cnt


def read_segments(lines, header: VcfHeader, sample_ids: Indexer,
                  chrom_ids: Indexer) -> Dict[Tuple[int, int, int, int],
                                              List[Segment]]:
    """``MergeIbdSegments.readSegments``, with segments grouped by sample
    pair and by the haplotype of each sample, in order of first
    appearance."""
    pairs: Dict[Tuple[int, int, int, int], List[Segment]] = {}
    for raw in lines:
        line = java_trim(raw)
        if len(line) > 0:
            segment = parse_segment(line, sample_ids, chrom_ids)
            s1 = segment.sample_index1
            s2 = segment.sample_index2
            if header.sample_index(s1) > header.sample_index(s2):
                raise illegal_argument("Inconsistency in sample order:" + line)
            key = (s1, s2, segment.hap1, segment.hap2)
            group = pairs.get(key)
            if group is None:
                group = []
                pairs[key] = group
            group.append(segment)
    return pairs


def merge_with_list(to_merge: List[Segment], nxt: Segment, max_end: int,
                    data: Data, max_gap: float, max_discord: int) -> bool:
    """``MergeIbdSegments.mergeWithList`` with the haplotype rule.  All
    segments in ``to_merge`` and ``nxt`` belong to one haplotype group.  A
    segment that starts before the running end joins without a test; one
    that starts at the running end is tested as a zero-length gap."""
    if not to_merge:
        return True
    last = to_merge[-1]
    if last.chrom_index != nxt.chrom_index:
        return False
    if nxt.start < max_end:
        return True
    gap_start = data.marker_index(last.chrom_index, max_end)
    gap_end = data.marker_index(nxt.chrom_index, nxt.start)
    if gap_start == -1:
        raise exit_with_error("ERROR: position missing from VCF file: "
                              + str(last.chrom_index) + ":" + str(max_end))
    if gap_end == -1:
        # The Java program reports the segment's end position here.
        raise exit_with_error("ERROR: position missing from VCF file: "
                              + str(nxt.chrom_index) + ":" + str(nxt.end))
    gap_length = data.gen_distance(gap_start, gap_end)
    if gap_length > max_gap:
        return False
    discord_cnt = data.haplotype_discord_cnt(
        last.sample_index1, last.hap1, last.sample_index2, last.hap2,
        gap_start, gap_end)
    return discord_cnt <= max_discord


def print_and_clear_list(data: Data, to_merge: List[Segment],
                         sample_ids: Indexer,
                         chrom_ids: Indexer) -> Tuple[Tuple, str]:
    """``MergeIbdSegments.printAndClearList``: the output line for the chain
    (with its sort key), and the chain emptied.  The chain keeps its group's
    haplotypes; Browning's program writes 0 for a joined segment."""
    first = to_merge[0]
    score = -DOUBLE_MAX
    start = INT_MAX
    end = INT_MIN
    for seg in to_merge:
        if seg.start < start:
            start = seg.start
        if seg.end > end:
            end = seg.end
        if seg.score > score:
            score = seg.score
    start_pos = data.gen_map.gen_pos(first.chrom_index, start)
    end_pos = data.gen_map.gen_pos(first.chrom_index, end)
    gen_length = end_pos - start_pos
    hap1 = first.hap1
    hap2 = first.hap2
    line = "\t".join((sample_ids.id(first.sample_index1), str(hap1),
                      sample_ids.id(first.sample_index2), str(hap2),
                      chrom_ids.id(first.chrom_index), str(start), str(end),
                      java_decimal_format(score, 2),
                      java_decimal_format(gen_length, 3))) + "\n"
    key = (first.chrom_index, start, end, java_compare_key(score),
           first.sample_index1, first.sample_index2, hap1, hap2, line)
    to_merge.clear()
    return key, line


def run(args: List[str], stdin: BinaryIO) -> str:
    """Runs the program and returns its standard output.  Raises
    JavaTermination or JavaException where the Java program would stop."""
    if len(args) != 4:
        raise JavaTermination(stdout=usage() + "\n", status=0)
    vcf_path, map_path = args[0], args[1]
    max_gap = parse_and_check_max_gap(args[2])
    max_discord = parse_and_check_max_discord(args[3])
    chrom_ids = Indexer()
    sample_ids = Indexer()
    data = Data(vcf_path, map_path, chrom_ids, sample_ids)
    pairs = read_segments(from_std_in(stdin), data.header, sample_ids,
                          chrom_ids)
    out: List[Tuple[Tuple, str]] = []
    to_merge: List[Segment] = []
    for group in pairs.values():
        group.sort(key=Segment.sort_key)
        max_end = INT_MIN
        for nxt in group:
            if not merge_with_list(to_merge, nxt, max_end, data, max_gap,
                                   max_discord):
                out.append(print_and_clear_list(data, to_merge, sample_ids,
                                                chrom_ids))
                max_end = INT_MIN
            if nxt.end > max_end:
                max_end = nxt.end
            to_merge.append(nxt)
        out.append(print_and_clear_list(data, to_merge, sample_ids,
                                        chrom_ids))
    out.sort(key=lambda item: item[0])
    return "".join(line for _, line in out)


def main(argv: Optional[List[str]] = None) -> int:
    """Command-line entry point; returns the exit status."""
    args = [to_java_chars(a) for a in (sys.argv[1:] if argv is None
                                       else argv)]
    stdout = sys.stdout.buffer
    stderr = sys.stderr.buffer
    try:
        text = run(args, sys.stdin.buffer)
    except JavaTermination as termination:
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
    # FileUtil.stdOutPrintWriter: segments go to standard output through a
    # writer in the default charset (UTF-8 in Java 18 and later).
    stdout.write(java_encode(text))
    stdout.flush()
    return 0
