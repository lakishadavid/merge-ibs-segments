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
# Python from src/ibdutil/StringIdSegment.java in refined-ibd.17Jan20.102.zip
# (SHA-256 0b0abf48528ec53d3fec7f21a9742c6e5960857c5a225b0e49346179bc45e6ea).
# On September 30, 2026, haplotype 0 ceased to be accepted, for the
# haplotype-specific merging rule (see README). Later changes are recorded in
# the Git history.
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
"""IBD segments, converted from ibdutil.StringIdSegment.

A Refined IBD segment line has 8 or 9 white-space separated fields: sample 1,
haplotype of sample 1 (1 or 2), sample 2, haplotype of sample 2,
chromosome, start position, end position (both inclusive), LOD score and,
optionally, length in cM.  Browning's program also accepts haplotype 0
("unknown", which it writes for joined segments); this program does not,
because its merging rule needs each segment's haplotypes.
"""

from __future__ import annotations

import math
from typing import Tuple

from .ids import Indexer
from .java_compat import (get_fields_ws, illegal_argument,
                          java_compare_key, java_parse_float,
                          java_parse_int)


class Segment:
    """``ibdutil.StringIdSegment``: sample and chromosome identifiers are
    held as indices; the score and length are float32 values."""

    __slots__ = ("sample_index1", "hap1", "sample_index2", "hap2",
                 "chrom_index", "start", "end", "score", "length")

    def __init__(self, sample_index1: int, hap1: int, sample_index2: int,
                 hap2: int, chrom_index: int, start: int, end: int,
                 score: float, length: float):
        self.sample_index1 = sample_index1
        self.hap1 = hap1
        self.sample_index2 = sample_index2
        self.hap2 = hap2
        self.chrom_index = chrom_index
        self.start = start
        self.end = end
        self.score = score
        self.length = length

    def sort_key(self) -> Tuple:
        """``StringIdSegment.compareTo``: chromosome index, start, end,
        score (``Float.compare``), sample index 1, sample index 2,
        haplotype 1, haplotype 2."""
        return (self.chrom_index, self.start, self.end,
                java_compare_key(self.score), self.sample_index1,
                self.sample_index2, self.hap1, self.hap2)


def _check_haplotype(hap: int) -> None:
    """Accepts haplotypes 1 and 2 (Browning's check also accepts 0)."""
    if hap != 1 and hap != 2:
        raise illegal_argument("invalid hap: " + str(hap))


def parse_segment(line: str, sample_ids: Indexer,
                  chrom_ids: Indexer) -> Segment:
    """``StringIdSegment.parse(String)``."""
    fields = get_fields_ws(line)
    if len(fields) < 8 or len(fields) > 9:
        raise illegal_argument(
            "Line describing IBD segment does not contain 8 or 9 "
            "white-space delimited fields: " + line)
    hap1 = java_parse_int(fields[1])
    hap2 = java_parse_int(fields[3])
    start = java_parse_int(fields[5])
    end = java_parse_int(fields[6])
    score = java_parse_float(fields[7])
    length = math.nan if len(fields) == 8 else java_parse_float(fields[8])
    if start > end:
        raise illegal_argument("start=" + str(start) + " end=" + str(end))
    _check_haplotype(hap1)
    _check_haplotype(hap2)
    sample_index1 = sample_ids.get_index(fields[0])
    sample_index2 = sample_ids.get_index(fields[2])
    chrom_index = chrom_ids.get_index(fields[4])
    return Segment(sample_index1, hap1, sample_index2, hap2, chrom_index,
                   start, end, score, length)
