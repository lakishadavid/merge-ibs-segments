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
"""Synthetic test cases for the exact-conversion test.

Every case is built from code: no participant data is used or needed.  A case
is a VCF file, a PLINK map file, IBD segment lines for standard input, and the
[gap] and [discord] arguments.  The same seed always gives the same files.

The cases cover each branch of MergeIbdSegments and the classes it uses,
threshold boundaries, shared endpoints, duplicate and nested records, input
order, map interpolation and extrapolation, number parsing and rounding, and
seeded random data.  Error cases each contain a single error, because the
Java program reads and parses a VCF file on several threads and may report a
different one of several errors from run to run.
"""

from __future__ import annotations

import gzip
import random
import struct
import zlib
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"


@dataclass
class Case:
    """One run of both programs.  ``files`` maps a file name to its bytes;
    ``args`` refer to files by name; ``stdin`` is the IBD input."""

    name: str
    args: List[str]
    stdin: bytes
    files: Dict[str, bytes] = field(default_factory=dict)
    note: str = ""


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def vcf_text(samples: Sequence[str], records: Sequence[Tuple],
             meta: Sequence[str] = ("##fileformat=VCFv4.2",)) -> str:
    """records: (chrom, pos, genotypes) or (chrom, pos, genotypes, ref,
    alt, id, info, fmt)."""
    lines = list(meta)
    lines.append(HEADER + "".join("\t" + s for s in samples))
    for rec in records:
        chrom, pos, gts = rec[0], rec[1], rec[2]
        ref = rec[3] if len(rec) > 3 else "A"
        alt = rec[4] if len(rec) > 4 else "G"
        rid = rec[5] if len(rec) > 5 else "."
        info = rec[6] if len(rec) > 6 else "."
        fmt = rec[7] if len(rec) > 7 else "GT"
        lines.append("\t".join([chrom, str(pos), rid, ref, alt, ".", "PASS",
                                info, fmt] + list(gts)))
    return "\n".join(lines) + "\n"


def map_text(entries: Sequence[Tuple[str, float, int]]) -> str:
    return "".join(chrom + "\t.\t" + (cm if isinstance(cm, str) else repr(cm))
                   + "\t" + str(bp) + "\n" for chrom, cm, bp in entries)


def ibd_text(rows: Sequence[Sequence]) -> str:
    return "".join("\t".join(str(x) for x in row) + "\n" for row in rows)


def bgzip(data: bytes) -> bytes:
    """BGZF: blocks of at most 65280 input bytes, then the empty EOF
    block."""
    out = bytearray()
    for start in range(0, max(len(data), 1), 65280):
        block = data[start:start + 65280]
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        deflated = comp.compress(block) + comp.flush()
        bsize = 18 + len(deflated) + 8 - 1
        out += bytes([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0])
        out += struct.pack("<H", bsize)
        out += deflated
        out += struct.pack("<II", zlib.crc32(block) & 0xFFFFFFFF, len(block))
    # The BGZF end-of-file marker block defined in the SAM/BAM format
    # specification (hts-specs, SAMv1 section 4.1.2).
    out += bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    return bytes(out)


def simple_case(name: str, samples, records, map_entries, ibd_rows,
                gap: str = "0.6", discord: str = "1", note: str = "",
                vcf_name: str = "in.vcf", map_name: str = "in.map",
                vcf_bytes: Optional[bytes] = None,
                map_bytes: Optional[bytes] = None,
                stdin: Optional[bytes] = None) -> Case:
    files = {
        vcf_name: vcf_bytes if vcf_bytes is not None
        else vcf_text(samples, records).encode(),
        map_name: map_bytes if map_bytes is not None
        else map_text(map_entries).encode(),
    }
    return Case(name, [vcf_name, map_name, gap, discord],
                stdin if stdin is not None else ibd_text(ibd_rows).encode(),
                files, note)


# ---------------------------------------------------------------------------
# A small standard data set: 4 samples, chromosome 1, markers every 100 bp
# ---------------------------------------------------------------------------

S4 = ["A", "B", "C", "D"]


def _std_records(chrom: str = "1", n: int = 30, start: int = 100,
                 step: int = 100, discordant: Sequence[int] = ()):
    """A and B share an allele at every marker except the positions in
    ``discordant`` (A 0|0, B 1|1 there).  C and D vary."""
    recs = []
    for k in range(n):
        pos = start + k * step
        if pos in discordant:
            gts = ["0|0", "1|1", "0|1", "1|0"]
        else:
            gts = ["0|1", "1|0", "1|1", "0|0"] if k % 2 else \
                ["1|0", "0|1", "0|0", "1|1"]
        recs.append((chrom, pos, gts))
    return recs


def _std_map(chrom: str = "1"):
    # 0.01 cM per 100 bp from 100 to 3100: genPos(p) = (p - 100) / 10000.
    return [(chrom, 0.0, 100), (chrom, 0.3, 3100)]


def branch_cases() -> List[Case]:
    cases = []
    recs = _std_records()
    mp = _std_map()

    def add(name, rows, gap="0.6", discord="1", records=None, map_entries=None,
            note=""):
        cases.append(simple_case(name, S4, records or recs,
                                 map_entries or mp, rows, gap, discord, note))

    add("single_segment", [["A", 1, "B", 2, "1", 100, 500, 7.5, 0.04]])
    add("overlap_merge", [["A", 1, "B", 2, "1", 100, 500, 3.2, 0.04],
                          ["A", 2, "B", 1, "1", 400, 900, 5.1, 0.05]])
    add("shared_endpoint", [["A", 1, "B", 2, "1", 100, 500, 3.2, 0.04],
                            ["A", 1, "B", 2, "1", 500, 900, 5.1, 0.04]],
        note="start == running end: Java's overlap branch")
    add("nested", [["A", 1, "B", 1, "1", 100, 2000, 3.0, 0.19],
                   ["A", 2, "B", 2, "1", 500, 900, 9.0, 0.04],
                   ["A", 1, "B", 1, "1", 2500, 2900, 4.0, 0.04]], gap="0.1")
    add("duplicates", [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04]] * 3
        + [["A", 1, "B", 2, "1", 2500, 2900, 3.0, 0.04]] * 2, gap="0.05")
    add("gap_merge", [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
                      ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04]])
    add("gap_too_long", [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
                         ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04]],
        gap="0.02")
    # Map positions at the gap's end markers, so the gap length is
    # 1.1 - 0.5 in double arithmetic (0.6000000000000001); [gap] is set to
    # exactly that value, and to the double just below it.
    exact_map = [("1", 0.0, 100), ("1", 0.5, 500), ("1", 1.1, 1100),
                 ("1", 3.0, 3100)]
    add("gap_exactly_threshold",
        [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.5],
         ["A", 1, "B", 2, "1", 1100, 1500, 4.0, 0.4]],
        gap=repr(1.1 - 0.5), map_entries=exact_map,
        note="gap length equals [gap] exactly")
    add("gap_just_over_threshold",
        [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.5],
         ["A", 1, "B", 2, "1", 1100, 1500, 4.0, 0.4]],
        gap="0.59999999999999997", map_entries=exact_map)
    for n_disc, discord in ((0, "0"), (1, "1"), (2, "1"), (1, "0"), (2, "2")):
        disc = [600, 700][:n_disc]
        add("discord_%d_max_%s" % (n_disc, discord),
            [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
             ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04]],
            discord=discord, records=_std_records(discordant=disc))
    add("discord_at_boundary_markers",
        [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
         ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04]],
        discord="1", records=_std_records(discordant=[500, 800]),
        note="discordant markers at the gap's end markers are counted")
    add("chain_of_four", [["A", 1, "B", 2, "1", 100, 400, 3.0, 0.03],
                          ["A", 2, "B", 2, "1", 600, 900, 6.0, 0.03],
                          ["A", 1, "B", 1, "1", 1100, 1400, 2.0, 0.03],
                          ["A", 2, "B", 1, "1", 1600, 1900, 1.0, 0.03]])
    add("running_end_not_last", [["A", 1, "B", 2, "1", 100, 2000, 3.0, 0.19],
                                 ["A", 1, "B", 2, "1", 300, 600, 3.0, 0.03],
                                 ["A", 1, "B", 2, "1", 2300, 2600, 3.0, 0.03]],
        note="gap measured from the running end (2000), not the last segment")
    add("same_start_different_end", [["A", 1, "B", 2, "1", 100, 900, 3.0, 0.08],
                                     ["A", 2, "B", 2, "1", 100, 500, 3.0, 0.04],
                                     ["A", 1, "B", 1, "1", 100, 900, 2.5, 0.08]])
    add("many_pairs", [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
                       ["A", 1, "C", 2, "1", 100, 500, 3.0, 0.04],
                       ["A", 2, "C", 2, "1", 700, 900, 3.0, 0.04],
                       ["B", 1, "D", 2, "1", 100, 500, 3.0, 0.04],
                       ["C", 1, "D", 1, "1", 1500, 2500, 3.0, 0.04],
                       ["A", 1, "A", 2, "1", 100, 500, 3.0, 0.04],
                       ["A", 1, "A", 2, "1", 700, 900, 3.0, 0.04]])
    add("eight_fields", [["A", 1, "B", 2, "1", 100, 500, 3.0],
                         ["A", 1, "B", 2, "1", 600, 900, 3.0]])
    add("hap_zero_input", [["A", 0, "B", 0, "1", 100, 500, 3.0, 0.04],
                           ["A", 0, "B", 2, "1", 1500, 1900, 3.0, 0.04]])
    add("positions_not_in_vcf_single",
        [["A", 1, "B", 2, "1", 150, 450, 3.0, 0.03],
         ["A", 1, "C", 2, "1", 2950, 9999, 3.0, 0.03]],
        note="no gap test, so positions need not be VCF positions")
    add("positions_not_in_vcf_overlap",
        [["A", 1, "B", 2, "1", 150, 450, 3.0, 0.03],
         ["A", 1, "B", 2, "1", 250, 650, 3.0, 0.03]])
    add("sample_not_in_vcf_first",
        [["Z", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
         ["Z", 1, "B", 2, "1", 400, 900, 3.0, 0.04]],
        note="sample 1 absent from the VCF: allowed while no gap is tested")
    add("both_samples_not_in_vcf",
        [["Y", 1, "Z", 2, "1", 100, 500, 3.0, 0.04]])
    shuffled = [["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04],
                ["C", 1, "D", 2, "1", 100, 300, 1.0, 0.02],
                ["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
                ["A", 1, "B", 2, "1", 1300, 1400, 3.5, 0.01],
                ["C", 1, "D", 2, "1", 250, 700, 2.0, 0.04]]
    add("input_order_1", shuffled)
    add("input_order_2", list(reversed(shuffled)))
    add("input_order_3", shuffled[2:] + shuffled[:2])
    # Two chromosomes; map lists chromosome 2 first.
    recs2 = _std_records("1") + _std_records("2")
    map2 = _std_map("2") + _std_map("1")
    add("two_chromosomes", [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
                            ["A", 1, "B", 2, "2", 800, 1200, 4.0, 0.04],
                            ["A", 1, "B", 2, "2", 100, 500, 3.0, 0.04],
                            ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04]],
        records=recs2, map_entries=map2)
    return cases


def map_cases() -> List[Case]:
    """Interpolation and extrapolation: segments that begin before and end
    after the map, with map spans shorter and longer than 5 cM, repeated cM
    values and negative base positions."""
    cases = []
    recs = _std_records(n=30, start=100, step=100)
    maps = {
        "short_span": [("1", 1.0, 500), ("1", 2.0, 1000), ("1", 3.5, 2000)],
        "long_span": [("1", 0.0, 500), ("1", 2.0, 800), ("1", 6.0, 1200),
                      ("1", 9.0, 1500), ("1", 13.0, 2000), ("1", 20.0, 2400)],
        "exact_5cm": [("1", 0.0, 500), ("1", 5.0, 1000), ("1", 10.0, 1500)],
        "repeated_cm": [("1", 1.0, 500), ("1", 1.0, 700), ("1", 1.0, 900),
                        ("1", 6.0, 1300), ("1", 6.0, 1400), ("1", 11.0, 2000),
                        ("1", 11.0, 2100)],
        "negative_bp": [("1", 0.0, -5000), ("1", 4.0, 700), ("1", 8.0, 2000)],
        "negative_cm": [("1", -3.0, 500), ("1", -0.0, 800), ("1", 0.0, 900),
                        ("1", 2.5, 2000)],
        "fine_cm": [("1", 0.000123, 400), ("1", 0.000623, 800),
                    ("1", 0.001123, 1200), ("1", 0.0015, 1600),
                    ("1", 0.0025, 2000), ("1", 0.0035, 2400)],
        # Several map positions share the cM value that the 5 cM end rule
        # searches for, so the search's choice among equal values shows in
        # the extrapolated lengths (start side, then end side).
        "equal_cm_start": [("1", 0.0, 500), ("1", 5.0, 600), ("1", 5.0, 700),
                           ("1", 5.0, 800), ("1", 5.0, 900), ("1", 6.0, 1000),
                           ("1", 7.0, 1100), ("1", 8.0, 1200)],
        "equal_cm_end": [("1", 0.0, 500), ("1", 1.0, 600), ("1", 2.0, 700),
                         ("1", 3.0, 800), ("1", 3.0, 900), ("1", 3.0, 1000),
                         ("1", 3.0, 1100), ("1", 8.0, 1200)],
        "equal_cm_many": [("1", 0.0, 400)]
        + [("1", 5.0, 500 + 10 * k) for k in range(37)]
        + [("1", 10.0, 1300)] + [("1", 10.0, 1400 + 10 * k)
                                 for k in range(20)]
        + [("1", 15.0, 2800)],
    }
    rows = [["A", 1, "B", 2, "1", 100, 300, 3.0, 0.1],
            ["A", 1, "C", 2, "1", 100, 2900, 3.0, 0.1],
            ["A", 1, "D", 2, "1", 600, 1800, 3.0, 0.1],
            ["B", 1, "C", 2, "1", 2100, 3000, 3.0, 0.1],
            ["B", 1, "D", 2, "1", 500, 1000, 3.0, 0.1],
            ["C", 1, "D", 2, "1", 700, 700, 3.0, 0.1],
            ["A", 2, "A", 2, "1", -100000, 2000000000, 3.0, 0.1],
            ["B", 2, "B", 2, "1", -2147483648, 2147483647, 3.0, 0.1]]
    for name, mp in maps.items():
        cases.append(simple_case("map_" + name, S4, recs, mp, rows,
                                 gap="0.05"))
    return cases


SCORE_STRINGS = [
    "0.125", "0.375", "0.625", "0.875", "2.5", "3.5", "0.005", "0.015",
    "0.025", "1.005", "2.675", "99.995", "99.994999", "0.0049999", "-0.001",
    "-0.005", "-0.0", "0", "+3.5", ".5", "5.", "1e-50", "1.4e-45",
    "3.4028235e38", "3.4028236e38", "1e39", "NaN", "Infinity", "-Infinity",
    "0x1.8p1", "0x.8p1f", "1.5d", "2.25F", "16777217", "16777219",
    "1.00000005960464477539062501", "1.0000000596046448",
    "12345678.9", "1e7", "123456789012", "5e-5", "-5", "7.777", "1234.565",
    "0.994999", "0.995", "0.996", "9.995", "1.115", "-2.675",
    "3.4028234e38", "0.000001", "65504.004", "1E+2", "00012.500",
]


def number_cases() -> List[Case]:
    """Score parsing and DecimalFormat("0.##"); length DecimalFormat("0.###").

    Each score is on its own sample pair so that no segments merge.
    """
    n = 20
    samples = ["S%02d" % i for i in range(n)]
    pairs = [(i, j) for i in range(n) for j in range(i, n)]
    recs = []
    for k in range(40):
        recs.append(("1", 1000 + 1000 * k, ["0|1"] * n))
    rows = []
    for s, (i, j) in zip(SCORE_STRINGS, pairs):
        rows.append([samples[i], 1, samples[j], 2, "1", 1000, 2000, s, "1.5"])
    # Lengths: map positions with cM values chosen to make awkward
    # differences (DecimalFormat 0.### and its underflow edge).
    cms = ["0", "0.0005", "0.0015", "0.0025", "0.0035", "0.00049", "0.00051",
           "0.1234", "0.1235", "0.1245", "1.0005", "2.0015", "2.00149999",
           "12.3456", "12.3465", "99.9995", "100.0005", "123.4565", "250.0",
           "250.00050000000001", "300.1234567"]
    mp = []
    pos = 1000
    for cm in sorted(cms, key=float):
        mp.append(("1", cm, pos))
        pos += 1000
    mp_positions = [p for _, _, p in mp]
    cm_values = sorted(cms, key=float)
    k = len(rows)
    for a in range(len(mp_positions)):
        for b in range(a + 1, min(a + 4, len(mp_positions))):
            if k >= len(pairs):
                break
            i, j = pairs[k]
            rows.append([samples[i], 2, samples[j], 1, "1", mp_positions[a],
                         mp_positions[b], "3.0", "1"])
            k += 1
    cases = [simple_case("numbers", samples, recs, mp, rows, gap="0.001",
                         discord="0",
                         note="scores " + str(len(SCORE_STRINGS))
                         + ", cM values " + str(len(cm_values)))]
    return cases


def format_cases() -> List[Case]:
    """Input formats: FORMAT subfields, multi-digit and multi-allelic
    genotypes, compressed files, line endings, blank lines, white space."""
    cases = []
    rows = [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
            ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04],
            ["C", 1, "D", 2, "1", 100, 1200, 2.0, 0.11]]
    mp = _std_map()
    recs = _std_records()
    # FORMAT subfields.
    recs_fmt = [(c, p, [g + ":0.5:1,0" for g in gts], "A", "G", ".", ".",
                 "GT:DS:AD") for c, p, gts in recs]
    cases.append(simple_case("format_subfields", S4, recs_fmt, mp, rows))
    # Multi-allelic with multi-digit alleles.
    recs_multi = []
    for c, p, gts in recs:
        if p in (600, 700):
            recs_multi.append((c, p, ["10|2", "3|11", "0|10", "11|1"], "A",
                               "C,G,T,AC,AG,AT,CA,CC,CG,CT,GA"))
        else:
            recs_multi.append((c, p, gts))
    cases.append(simple_case("multi_allelic", S4, recs_multi, mp, rows,
                             discord="0"))
    # Extra sample fields beyond the header are ignored by Java.
    recs_extra = [(c, p, gts + ["1|1"]) for c, p, gts in recs]
    cases.append(simple_case("extra_sample_fields", S4, recs_extra, mp, rows))
    # IDs, INFO END, symbolic and star alleles, lower-case bases.
    recs_misc = []
    for c, p, gts in recs:
        if p == 300:
            recs_misc.append((c, p, gts, "a", "<DEL>", "rs1;.;rs2",
                              "SVTYPE=DEL;END=350;X=1"))
        elif p == 400:
            recs_misc.append((c, p, gts, "ACGTN", "*", ".;rs9"))
        elif p == 900:
            recs_misc.append((c, p, ["0|0", "0|0", "0|0", "0|0"], "T", "."))
        else:
            recs_misc.append((c, p, gts))
    cases.append(simple_case("vcf_field_variants", S4, recs_misc, mp, rows))
    # Compressed inputs.
    vcf = vcf_text(S4, recs).encode()
    mpb = map_text(mp).encode()
    cases.append(simple_case("vcf_gzip", S4, recs, mp, rows,
                             vcf_name="in.vcf.gz",
                             vcf_bytes=gzip.compress(vcf, mtime=0)))
    cases.append(simple_case("vcf_bgzip", S4, recs, mp, rows,
                             vcf_name="in.vcf.gz", vcf_bytes=bgzip(vcf)))
    cases.append(simple_case("map_gzip", S4, recs, mp, rows,
                             map_name="in.map.gz",
                             map_bytes=gzip.compress(mpb, mtime=0)))
    two_members = gzip.compress(vcf[:200], mtime=0) + gzip.compress(vcf[200:],
                                                                     mtime=0)
    cases.append(simple_case("vcf_gzip_two_members", S4, recs, mp, rows,
                             vcf_name="in.vcf.gz", vcf_bytes=two_members))
    # Line endings and blank lines.
    crlf = vcf.replace(b"\n", b"\r\n")
    cases.append(simple_case("crlf", S4, recs, mp, rows, vcf_bytes=crlf,
                             map_bytes=mpb.replace(b"\n", b"\r\n"),
                             stdin=ibd_text(rows).encode().replace(b"\n",
                                                                   b"\r\n")))
    cr = vcf.replace(b"\n", b"\r")
    cases.append(simple_case("cr_only", S4, recs, mp, rows, vcf_bytes=cr,
                             map_bytes=mpb.replace(b"\n", b"\r"),
                             stdin=ibd_text(rows).encode().replace(b"\n",
                                                                   b"\r")))
    lines = vcf.split(b"\n")
    blank_inside = b"\n".join(lines[:5] + [b"", b"  \t ", b""] + lines[5:])
    cases.append(simple_case("vcf_blank_lines_inside", S4, recs, mp, rows,
                             vcf_bytes=blank_inside))
    cases.append(simple_case("vcf_no_final_newline", S4, recs, mp, rows,
                             vcf_bytes=vcf.rstrip(b"\n")))
    header_ws = vcf.replace(b"##fileformat=VCFv4.2\n",
                            b"  ##fileformat=VCFv4.2  \n##x=y\n")
    cases.append(simple_case("vcf_meta_whitespace", S4, recs, mp, rows,
                             vcf_bytes=header_ws))
    ibd_ws = ("\n  A 1\tB  2 1 100 500 3.0 0.04  \n\n\t\n"
              "A\t1\tB\t2\t1\t800\t1200\t4.0\t0.04\x0b\n"
              "C\x011\x02D\x1f2 1 100 1200 2.0 0.11").encode()
    cases.append(simple_case("ibd_whitespace", S4, recs, mp, rows,
                             stdin=ibd_ws))
    map_ws = b"\n  1 . 0.0 100\n\n1\tsnp1\t0.3\t3100  \n\n"
    cases.append(simple_case("map_whitespace", S4, recs, mp, rows,
                             map_bytes=map_ws))
    cases.append(Case("empty_ibd_input", ["in.vcf", "in.map", "0.6", "1"],
                      b"", {"in.vcf": vcf, "in.map": mpb}))
    cases.append(Case("blank_ibd_input", ["in.vcf", "in.map", "0.6", "1"],
                      b"\n \n\t\n", {"in.vcf": vcf, "in.map": mpb}))
    # Duplicate VCF positions: the last record at a position is used.
    recs_dup = recs[:6] + [("1", 600, ["0|0", "1|1", "0|0", "1|1"])] + recs[6:]
    cases.append(simple_case("vcf_duplicate_position", S4, recs_dup, mp,
                             rows, discord="0"))
    # Unicode sample identifiers (UTF-8).
    uni = ["Ä", "éè", "\U0001d7d9", "D"]
    rows_uni = [[uni[0], 1, uni[1], 2, "1", 100, 500, 3.0, 0.04],
                [uni[0], 1, uni[1], 2, "1", 800, 1200, 4.0, 0.04],
                [uni[2], 1, uni[3], 2, "1", 100, 1200, 2.0, 0.11]]
    cases.append(simple_case("unicode_ids", uni, recs, mp, rows_uni))
    return cases


def argument_cases() -> List[Case]:
    vcf = vcf_text(S4, _std_records()).encode()
    mpb = map_text(_std_map()).encode()
    ibd = ibd_text([["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
                    ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04]]).encode()
    files = {"in.vcf": vcf, "in.map": mpb}
    cases = [Case("args_none", [], b"", files),
             Case("args_three", ["in.vcf", "in.map", "0.6"], b"", files),
             Case("args_five", ["in.vcf", "in.map", "0.6", "1", "x"], b"",
                  files)]
    gaps = ["abc", "inf", "Infinity", "-1", "0", "-0", "NaN", "1e400",
            "0x1p-2", " 0.6 ", "0.6f", "", "1.2.3", "1e", ".", "0x1p",
            "4.9e-324", "0.6e0d", "+.6", "٠.6"]
    for k, g in enumerate(gaps):
        cases.append(Case("gap_arg_%02d" % k, ["in.vcf", "in.map", g, "1"],
                          ibd, files, note=repr(g)))
    discords = ["-1", "1.5", "", "+1", "2147483648", "2147483647", " 1",
                "١", "0"]
    for k, d in enumerate(discords):
        cases.append(Case("discord_arg_%02d" % k, ["in.vcf", "in.map", "0.6",
                                                   d], ibd, files,
                          note=repr(d)))
    return cases


def error_cases() -> List[Case]:
    """Each case has exactly one error."""
    cases = []
    recs = _std_records()
    mp = _std_map()
    rows = [["A", 1, "B", 2, "1", 100, 500, 3.0, 0.04],
            ["A", 1, "B", 2, "1", 800, 1200, 4.0, 0.04]]
    vcf = vcf_text(S4, recs)
    mpt = map_text(mp)
    ibd = ibd_text(rows)

    def add(name, vcf_s=None, map_s=None, ibd_s=None, args=None,
            vcf_name="in.vcf", map_name="in.map", raw_vcf=None, raw_map=None):
        files = {}
        files[vcf_name] = raw_vcf if raw_vcf is not None else \
            (vcf if vcf_s is None else vcf_s).encode()
        files[map_name] = raw_map if raw_map is not None else \
            (mpt if map_s is None else map_s).encode()
        cases.append(Case(name, args or [vcf_name, map_name, "0.6", "1"],
                          (ibd if ibd_s is None else ibd_s).encode(), files))

    # Files that cannot be read.
    add("missing_vcf_file", args=["nosuch.vcf", "in.map", "0.6", "1"])
    add("missing_map_file", args=["in.vcf", "nosuch.map", "0.6", "1"])
    add("vcf_not_gzip", vcf_name="in.vcf.gz")
    add("vcf_gz_empty", vcf_name="in.vcf.gz", raw_vcf=b"")
    add("vcf_gz_truncated", vcf_name="in.vcf.gz",
        raw_vcf=gzip.compress(vcf.encode(), mtime=0)[:120])
    # Genetic map errors.
    add("map_three_fields", map_s="1 . 0.0\n1 . 0.3 3100\n")
    add("map_five_fields", map_s="1 . 0.0 100 x\n1 . 0.3 3100\n")
    add("map_bad_bp", map_s="1 . 0.0 100.5\n1 . 0.3 3100\n")
    add("map_bad_cm", map_s="1 . zero 100\n1 . 0.3 3100\n")
    add("map_nan_cm", map_s="1 . NaN 100\n1 . 0.3 3100\n")
    add("map_inf_cm", map_s="1 . -Infinity 100\n1 . 0.3 3100\n")
    add("map_duplicate_bp", map_s="1 . 0.0 100\n1 . 0.1 100\n1 . 0.3 3100\n")
    add("map_descending_bp", map_s="1 . 0.0 200\n1 . 0.1 100\n1 . 0.3 3100\n")
    add("map_descending_cm", map_s="1 . 0.2 100\n1 . 0.1 200\n1 . 0.3 3100\n")
    add("map_one_position", map_s="1 . 0.0 100\n")
    add("map_equal_cm", map_s="1 . 0.3 100\n1 . 0.3 3100\n")
    add("map_second_chrom_bad", map_s=mpt + "2 . 0.0 100\n2 . 0.0 200\n")
    add("map_missing_chrom_at_output", map_s="2 . 0.0 100\n2 . 0.3 3100\n")
    # VCF header errors.
    add("vcf_no_header", vcf_s="##fileformat=VCFv4.2\n")
    add("vcf_header_not_hash", vcf_s="CHROM\tPOS\n")
    add("vcf_header_bad_prefix", vcf_s=vcf.replace("#CHROM\tPOS", "#CHROM POS"))
    add("vcf_meta_no_equals", vcf_s="##fileformat\n" + vcf)
    add("vcf_meta_equals_last", vcf_s="##fileformat=\n" + vcf)
    add("vcf_empty_sample_id", vcf_s=vcf.replace("\tA\tB\tC\tD", "\tA\t\tC\tD"))
    add("vcf_duplicate_sample", vcf_s=vcf.replace("\tA\tB\tC\tD", "\tA\tB\tC\tA"))
    add("vcf_short_header_with_record",
        vcf_s=vcf.replace("\tFORMAT\tA\tB\tC\tD", ""))
    header_only = vcf.split("\n")[0] + "\n" + vcf.split("\n")[1] + "\n"
    add("vcf_no_records", vcf_s=header_only)
    add("vcf_header_then_blank_line", vcf_s=header_only + "\n\n")
    add("vcf_header_then_spaces", vcf_s=header_only + "   \n")
    add("vcf_trailing_blank_line", vcf_s=vcf + "\n")
    # VCF record errors.
    good = vcf.split("\n")[2]

    def with_record(line: str) -> str:
        parts = vcf.split("\n")
        parts.insert(5, line)
        return "\n".join(parts)

    add("rec_no_tab", vcf_s=with_record("1 900 . A G"))
    add("rec_seven_fields", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS"))
    add("rec_missing_chrom_short", vcf_s=with_record(good.replace("1\t100", ".\t950", 1)))
    add("rec_missing_chrom_long",
        vcf_s=with_record(".\t950\t" + "x" * 90 + "\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1"))
    add("rec_chrom_colon_short", vcf_s=with_record(good.replace("1\t100", "1:2\t950", 1)))
    add("rec_chrom_colon_long",
        vcf_s=with_record("1:2\t950\t" + "x" * 90 + "\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1"))
    add("rec_bad_pos", vcf_s=with_record(good.replace("\t100\t", "\t95x\t", 1)))
    add("rec_bad_pos_long",
        vcf_s=with_record("1\t95x\t" + "x" * 90 + "\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t0|1"))
    add("rec_empty_pos", vcf_s=with_record(good.replace("\t100\t", "\t\t", 1)))
    add("rec_pos_overflow", vcf_s=with_record(good.replace("\t100\t", "\t99999999999\t", 1)))
    add("rec_empty_id", vcf_s=with_record(good.replace("\t100\t.\t", "\t950\t\t", 1)))
    add("rec_id_space", vcf_s=with_record(good.replace("\t100\t.\t", "\t950\trs 1\t", 1)))
    add("rec_bad_ref", vcf_s=with_record(good.replace("\t100\t.\tA\t", "\t950\t.\tX\t", 1)))
    add("rec_empty_ref", vcf_s=with_record(good.replace("\t100\t.\tA\t", "\t950\t.\t\t", 1)))
    add("rec_empty_alt", vcf_s=with_record(good.replace("\tA\tG\t", "\tA\t\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_bad_alt", vcf_s=with_record(good.replace("\tA\tG\t", "\tA\tG,.\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_bad_symbolic_alt", vcf_s=with_record(good.replace("\tA\tG\t", "\tA\t<D,EL>\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_duplicate_allele", vcf_s=with_record(good.replace("\tA\tG\t", "\tA\tG,G\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_end_last_info_key", vcf_s=with_record(good.replace("\tPASS\t.\t", "\tPASS\tX=1;END=990\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_end_first_info_key", vcf_s=with_record(good.replace("\tPASS\t.\t", "\tPASS\tEND=10;X=1\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_end_before_pos_2", vcf_s=with_record(good.replace("\tPASS\t.\t", "\tPASS\tX=1;END=10;Y=2\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_end_empty", vcf_s=with_record(good.replace("\tPASS\t.\t", "\tPASS\tX=1;END=;Y=2\t", 1).replace("\t100\t", "\t950\t", 1)))
    add("rec_eight_tabs", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT"))
    add("rec_too_few_samples", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0"))
    add("rec_no_separator", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0\t0|1\t0|1"))
    add("rec_last_no_separator", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t1"))
    add("rec_colon_before_separator", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0:1|1\t0|1\t0|1"))
    add("rec_empty_genotype_last", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|1\t0|1\t"))
    add("rec_unphased", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0/1\t0|1\t0|1"))
    add("rec_missing_allele", vcf_s=with_record("1\t950\trs5;.\tA\tG\t.\tPASS\t.\tGT\t0|1\t.|1\t0|1\t0|1"))
    add("rec_missing_allele_ids_moved", vcf_s=with_record("1\t950\t.;rs5;;rs6\tA\tG,T\t.\tPASS\t.\tGT\t0|1\t1|.\t0|1\t0|1"))
    add("rec_invalid_allele", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t2|1\t0|1\t0|1"))
    add("rec_invalid_allele_letter", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\tA|1\t0|1\t0|1"))
    add("rec_invalid_allele_multi", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1a|1\t0|1\t0|1"))
    add("rec_negative_allele", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t-1|1\t0|1\t0|1"))
    add("rec_empty_first_allele", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t|1\t0|1\t0|1"))
    add("rec_empty_second_allele", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|\t0|1\t0|1"))
    add("rec_unphased_missing", vcf_s=with_record("1\t950\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t./.\t0|1\t0|1"))
    # IBD input errors.
    add("ibd_seven_fields", ibd_s="A 1 B 2 1 100 500\n")
    add("ibd_ten_fields", ibd_s="A 1 B 2 1 100 500 3.0 0.04 x\n")
    add("ibd_bad_hap", ibd_s="A x B 2 1 100 500 3.0 0.04\n")
    add("ibd_hap_three", ibd_s="A 3 B 2 1 100 500 3.0 0.04\n")
    add("ibd_hap_negative", ibd_s="A 1 B -1 1 100 500 3.0 0.04\n")
    add("ibd_start_after_end", ibd_s="A 1 B 2 1 600 500 3.0 0.04\n")
    add("ibd_bad_start", ibd_s="A 1 B 2 1 1e2 500 3.0 0.04\n")
    add("ibd_bad_score", ibd_s="A 1 B 2 1 100 500 high 0.04\n")
    add("ibd_score_two_points", ibd_s="A 1 B 2 1 100 500 1.2.3 0.04\n")
    add("ibd_bad_length", ibd_s="A 1 B 2 1 100 500 3.0 0.0.4\n")
    add("ibd_sample_order", ibd_s="B 1 A 2 1 100 500 3.0 0.04\n")
    add("ibd_second_sample_absent", ibd_s="A 1 Z 2 1 100 500 3.0 0.04\n")
    add("ibd_chrom_not_in_map", ibd_s=ibd + "A 1 C 2 7 100 500 3.0 0.04\n")
    add("ibd_gap_end_missing", ibd_s="A 1 B 2 1 100 500 3.0 0.04\nA 1 B 2 1 750 1200 3.0 0.04\n")
    add("ibd_gap_start_missing", ibd_s="A 1 B 2 1 100 550 3.0 0.04\nA 1 B 2 1 800 1200 3.0 0.04\n")
    add("ibd_gap_chrom_not_in_vcf", map_s=mpt + "5 . 0.0 100\n5 . 1.0 3100\n",
        ibd_s="A 1 B 2 5 100 500 3.0 0.04\nA 1 B 2 5 800 1200 3.0 0.04\n")
    add("ibd_absent_sample_gap_test",
        ibd_s="Z 1 B 2 1 100 500 3.0 0.04\nZ 1 B 2 1 800 1200 3.0 0.04\n")
    unsorted = vcf_text(S4, [recs[0], recs[1], recs[2], recs[3], recs[4],
                             recs[9], recs[8], recs[7], recs[6], recs[5]]
                        + recs[10:])
    add("vcf_unsorted_gap", vcf_s=unsorted,
        ibd_s="A 1 B 2 1 100 900 3.0 0.04\nA 1 B 2 1 1000 1500 3.0 0.04\n")
    return cases


def random_case(seed: int, n_samples: int = 30, n_markers: int = 600,
                chroms: Sequence[str] = ("1", "2", "X"),
                n_clusters: int = 250) -> Case:
    """Seeded random data: phased genotypes with common and rare alleles,
    some multi-allelic markers and FORMAT subfields; a genetic map at a
    subset of positions with repeated cM values; clusters of IBD segments
    that overlap, nest, touch or leave gaps."""
    rng = random.Random(seed)
    samples = ["s%d_%d" % (seed, i) for i in range(n_samples)]
    records = []
    positions: Dict[str, List[int]] = {}
    use_fmt = rng.random() < 0.5
    for chrom in chroms:
        pos = rng.randint(1, 5000)
        plist = []
        for _ in range(n_markers):
            pos += rng.randint(1, 3000)
            plist.append(pos)
            n_alleles = 2 if rng.random() < 0.9 else rng.randint(3, 4)
            if rng.random() < 0.15:
                freqs = [0.999] + [0.001 / (n_alleles - 1)] * (n_alleles - 1)
            else:
                raw = [rng.random() for _ in range(n_alleles)]
                total = sum(raw)
                freqs = [r / total for r in raw]
            gts = []
            for _ in range(n_samples):
                a = rng.choices(range(n_alleles), freqs)[0]
                b = rng.choices(range(n_alleles), freqs)[0]
                gt = "%d|%d" % (a, b)
                gts.append(gt + ":%.2f" % rng.random() if use_fmt else gt)
            alts = ",".join("CGT"[:n_alleles - 1])
            records.append((chrom, pos, gts, "A", alts, ".", ".",
                            "GT:DS" if use_fmt else "GT"))
        positions[chrom] = plist
    map_entries = []
    for chrom in chroms:
        plist = positions[chrom]
        cm = rng.uniform(-1.0, 3.0)
        chosen = sorted(set([plist[0] - rng.randint(0, 5000)]
                            + rng.sample(plist, len(plist) // 7)
                            + [plist[-1] + rng.randint(0, 5000)]))
        digits = rng.choice([3, 4, 6])
        for p in chosen:
            cm += 0.0 if rng.random() < 0.1 else rng.uniform(0.0, 0.8)
            map_entries.append((chrom, round(cm, digits), p))
    rows = []
    for _ in range(n_clusters):
        i = rng.randrange(n_samples)
        j = rng.randrange(n_samples) if rng.random() < 0.97 else i
        i, j = min(i, j), max(i, j)
        chrom = rng.choice(chroms)
        plist = positions[chrom]
        k = rng.randrange(len(plist))
        for _ in range(rng.randint(1, 6)):
            length = rng.randint(0, 40)
            s = min(k, len(plist) - 1)
            e = min(s + length, len(plist) - 1)
            h1 = rng.choice([1, 2, 1, 2, 0])
            h2 = rng.choice([1, 2, 1, 2, 0])
            score = rng.choice(["%.2f", "%.3f", "%.1f", "%.4f"]) % \
                rng.uniform(-1, 60)
            row = [samples[i], h1, samples[j], h2, chrom, plist[s], plist[e],
                   score]
            if rng.random() < 0.8:
                row.append("%.3f" % rng.uniform(0, 5))
            rows.append(row)
            k = e + rng.choice([-10, -3, -1, 0, 0, 1, 2, 3, 5, 8, 15])
            k = max(0, k)
    rng.shuffle(rows)
    gap = rng.choice(["0.6", "1", "0.25", "2.5", "0.05"])
    discord = rng.choice(["0", "1", "1", "2", "5"])
    return simple_case("random_seed_%d" % seed, samples, records,
                       map_entries, rows, gap, discord,
                       note="gap %s discord %s, %d segment lines"
                       % (gap, discord, len(rows)))


def _float32(x: float) -> float:
    return struct.unpack("<f", struct.pack("<f", x))[0]


def _length_values(rng: random.Random) -> List[float]:
    """Doubles for DecimalFormat("0.###"): a dense decimal grid, exact
    binary ties at the fourth decimal, values near the rounding edges below
    0.001, large magnitudes and seeded random values."""
    vals = set()
    for k in range(0, 20001):
        vals.add(k / 10000)
    for k in range(0, 4001):
        vals.add(k / 16)
        vals.add(k / 2048)
    for k in range(1, 2000):
        vals.add(k / 1000000)
        vals.add(k / 100000)
    for s in ("0.0005", "0.00049999999999999", "0.00050000000000001",
              "0.00045", "0.00055", "0.0004", "0.00005", "0.0001",
              "0.0009999", "0.00095", "0.0015", "0.0025", "0.00150000000001",
              "5e-5", "5e-6", "0.00051", "0.000501", "0.0006"):
        vals.add(float(s))
    for e in range(0, 23):
        vals.add(10.0 ** e)
        vals.add(10.0 ** e + 0.0005)
        vals.add(1.5 * 10.0 ** e)
    for x in (123456789.0625, 2.0 ** 53, 2.0 ** 53 + 2, 2.0 ** 60, 1e22, 1e23,
              12345678901234.5625, 99999999999.9995, 4503599627370495.5):
        vals.add(x)
    # [2**43, 2**53): values whose shortest decimal form has few fraction
    # digits although their exact binary value has more.
    for x in (8796093022208.03, 8796093022208.0625, 70368744177664.1,
              1125899906842624.2, 2.0 ** 46 + 2.0 ** -6, 2.0 ** 50 + 0.25,
              2.0 ** 52 + 0.5, 4503599627370495.5, 9007199254740991.0):
        vals.add(x)
    for _ in range(3000):
        e = rng.uniform(43, 53)
        v = 2.0 ** e
        vals.add(float(repr(round(v, rng.choice([1, 2, 3])))))
        vals.add(v)
    for _ in range(6000):
        vals.add(round(rng.uniform(0, 3), rng.choice([4, 5, 6, 8, 17])))
        vals.add(rng.uniform(0, 400))
        vals.add(10 ** rng.uniform(-6, 6))
    return sorted(v for v in vals if v >= 0.0)


def _score_strings(rng: random.Random) -> List[str]:
    """Strings for Float.parseFloat and DecimalFormat("0.##")."""
    out = []
    for k in range(-2000, 20001):
        out.append("%.3f" % (k / 1000))
    for k in range(-400, 4001):
        out.append(repr(k / 8))
        out.append(repr(k / 64))
    for _ in range(6000):
        bits = rng.getrandbits(32)
        f = struct.unpack("<f", struct.pack("<I", bits))[0]
        if f == f and abs(f) != float("inf"):
            out.append(repr(f))
        out.append(repr(_float32(rng.uniform(-100, 1000))))
        out.append("%.4f" % rng.uniform(-10, 100))
    out += ["1e-50", "1.4e-45", "7e-46", "1.401298464324817e-45",
            "3.4028235e38", "3.40282356779733661637539395458142568448e38",
            "3.4028235677973366e38", "0.0000000000000000000000001",
            "0x1.fffffep127", "0x1.ffffffp127", "0x1p-149", "0x1p-150",
            "0x1.8p-150", "+0x1P+3", "-0X.8P-1f", "1.D", ".5F", "5.e-1d",
            "0001.5000", "-0000", "+.5", "1E+2", "1e-0", "99.995", "99.985",
            "0.005", "0.015", "0.025", "0.035", "-0.005", "-0.015",
            "1.00000005960464477539062501", "1.0000000596046448",
            "16777217", "16777219", "33554435", "0.1", "0.2", "0.3"]
    return out


def vector_cases() -> List[Case]:
    """Record the jar's number formatting and parsing.

    The first case prints DecimalFormat("0.##") of many float scores and
    DecimalFormat("0.###") of many doubles: each output line is one segment
    of its own sample pair, running from a map position at 0 cM to a map
    position whose cM value is the double to be formatted.  The other cases
    each give one unusual number string to Float.parseFloat or
    Integer.parseInt, to record whether the jar accepts it and how it fails.
    """
    rng = random.Random(20260929)
    lengths = _length_values(rng)
    scores = _score_strings(rng)
    n_rows = max(len(lengths), len(scores))
    n = 1
    while n * (n + 1) // 2 < n_rows:
        n += 1
    samples = ["v%03d" % i for i in range(n)]
    pairs = [(i, j) for i in range(n) for j in range(i, n)]
    recs = [("1", 1, ["0|1"] * n), ("1", 2, ["1|0"] * n)]
    mp = [("1", "0.0", 1)]
    for k, v in enumerate(lengths):
        mp.append(("1", repr(v), 100 + k))
    rows = []
    for k in range(n_rows):
        i, j = pairs[k]
        end = 100 + (k % len(lengths))
        rows.append([samples[i], 1, samples[j], 2, "1", 1, end,
                     scores[k % len(scores)]])
    cases = [simple_case("vectors_decimal_format", samples, recs, mp, rows,
                         note="%d lengths, %d scores" % (len(lengths),
                                                         len(scores)))]
    recs4 = _std_records()
    mp4 = _std_map()
    float_strings = [
        "1.2.3", "..", "1..", ".", "-", "+", "e5", "1e", "1e+", "1e-", "1ee5",
        "1e5.5", "1.5f5", "1.5ff", "1.5fd", "0x", "0x1", "0x1p", "0x.p1",
        "0xp1", "0x1.8", "0x1.8p1x", "NaN1", "nan", "NAN", "Infinit",
        "infinity", "+-1", "--1", "1-", "1+", "1,5", "1_000", "١.٥",
        "٣", "１", "¹", "1.0e1000000000000", "1e2147483648",
        "0x1p2147483648", "0x1p-2147483649", "1d", "1D", "1F", "-.5e-3f",
        "+.e1", ".e1", "0.", "00", "-0x0p0", "0X1P1", "1.e5", "Infinityf",
        "NaNd", "1.2.x", "12.3e4.5", "1x.2.3", "0x1.8p1.5", "+NaN", "-NaN",
        "-Infinity", "+Infinity", "0x1.8P1D", "1e+05", "1E-5f", "0.5e",
        "1.5E+", "0x1p+", "0x1p-", "0x1.p1", "0x.1p1", "0xg", "1x", "x1",
        ".5.", "5..", "0.0.0", "-.", "+.", "-e1", "1.5 ",
        "0x1..p1", "0x1.2.3p1", "0x..p1", "-0x1..p1", "-0X.1.2p3", "0X.1.2",
        "0x1.8.", "0x.", "0x..", "1" * 700 + "x" + "2" * 800,
        "1." + "5" * 1500 + ".", "1" * 5000, "0." + "0" * 3000 + "15",
        "0x" + "f" * 2000 + "p-8000", "1e" + "9" * 5000]
    # Invalid strings around 1000 characters, with position-revealing
    # content, to record how the jar shortens long strings in its message.
    for length in (999, 1000, 1001, 1002, 1003, 1004, 2001):
        counter = "".join("%05d" % i for i in range(length))
        float_strings.append(counter[:length - 1] + "x")
    for k, s in enumerate(float_strings):
        cases.append(simple_case(
            "float_string_%02d" % k, S4, recs4, mp4,
            [["A", 1, "B", 2, "1", 100, 500, s, 0.04]], note=repr(s)))
    int_strings = ["+1", "01", "0001", "١", "１", "१",
                   "\U0001d7d9", "1_0", "x", "2147483648", "-2147483649",
                   "+", "-", "1.0", "+-1", "4294967297", "-0", "+2",
                   "0" * 5000 + "1", "1" * 5000, "-" + "0" * 20 + "2"]
    for k, s in enumerate(int_strings):
        cases.append(simple_case(
            "int_string_%02d" % k, S4, recs4, mp4,
            [["A", s, "B", 2, "1", 100, 500, 3.0, 0.04]], note=repr(s)))
    return cases


def all_cases(n_random: int = 20) -> List[Case]:
    cases = (branch_cases() + map_cases() + number_cases() + format_cases()
             + argument_cases() + error_cases() + vector_cases())
    for seed in range(1, n_random + 1):
        big = seed % 5 == 0
        cases.append(random_case(seed, n_samples=60 if big else 30,
                                 n_markers=2000 if big else 600,
                                 n_clusters=1500 if big else 250))
    names = [c.name for c in cases]
    assert len(names) == len(set(names)), "duplicate case names"
    return cases
