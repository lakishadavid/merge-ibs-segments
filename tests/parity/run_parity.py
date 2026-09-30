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
#
# This script runs Brian L. Browning's merge-ibd-segments jar, which is not
# distributed with merge-ibs-segments.  The jar
# (merge-ibd-segments.17Jan20.102.jar, SHA-256
# f5a8e8d094e99fa8226e8489e2b55e4a2bc89b325c49de19501395b15769dc80) and its
# source (refined-ibd.17Jan20.102.zip) are published at
# https://faculty.washington.edu/browning/refined-ibd.html.
"""Exact-conversion test: runs Browning's merge-ibd-segments jar and this
Python conversion on the same synthetic inputs and compares them.

Usage (inside a container that has Java and the reference jar):

    python tests/parity/run_parity.py --jar /path/merge-ibd-segments.17Jan20.102.jar \
        --work /tmp/mis_parity --report tests/parity/parity_record.md

Comparison rules:

- Both runs succeed (exit status 0): the complete output lines are compared
  as a multiset, with duplicate counts.  The Python output is sorted; order
  is the only accepted difference.  For the usage message (wrong number of
  arguments) the Java command is replaced by the Python command before
  comparing.
- Either run fails: the exit statuses must be equal, and so must the
  substantive diagnostic, which is standard output and standard error after
  these normalisations: Java stack-frame lines (``\\tat ...``,
  ``\\t... N more``) are removed; an exception that Java re-threw from a
  parallel stream is reduced to the original exception; Java object
  identity hashes (``@1b6d3586``) are removed; the usage message is
  replaced by a marker; segment lines the Java program had already written
  before an error are removed (the Python conversion writes no segments
  when it stops with an error).

Two differences are known and accepted (September 30, 2026); a case showing
one of them is reported as KNOWN, and only if nothing else differs:

- Formatted numbers of 2**53 and above: the jar prints digits the
  conversion does not reproduce.  A successful case is KNOWN when every
  differing output line differs only in the score or length column, and
  every differing value is at least 2**53 in magnitude.
- A gap test for a sample absent from the VCF file: the jar reads outside
  its genotype arrays.  The case ibd_absent_sample_gap_test is KNOWN when
  both programs exit with status 1 and report
  java.lang.IndexOutOfBoundsException.
"""

from __future__ import annotations

import argparse
import collections
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from typing import Dict, List, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

from cases import Case, all_cases  # noqa: E402

PY_COMMAND = "python -m merge_ibs_segments"
JAVA_COMMAND = re.compile(r"java -jar merge-ibd-segments\.[^ ]*\.jar")
USAGE_START = "usage: cat [in] | "
SOURCE_ZIP_SHA256 = ("0B0ABF48528EC53D3FEC7F21A9742C6E"
                     "5960857C5A225B0E49346179BC45E6EA")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd: List[str], cwd: str, stdin: bytes) -> Tuple[int, str, str, float]:
    start = time.time()
    p = subprocess.run(cmd, cwd=cwd, input=stdin, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, timeout=600)
    return (p.returncode, p.stdout.decode("utf-8", "replace"),
            p.stderr.decode("utf-8", "replace"), time.time() - start)


def _collapse_usage(text: str) -> str:
    """Replaces the usage message (it begins with an empty line and ends
    with the line about the maximal score) with a marker."""
    text = JAVA_COMMAND.sub(PY_COMMAND, text)
    pattern = re.compile(r"\nusage: cat \[in\] \| .*?maximal score of the "
                         r"merged segments\.\n\n?", re.S)
    return pattern.sub("<usage>\n", text)


_EXC_HEADER = re.compile(r'^(Exception in thread "main" |Caused by: )(.*)$')
_FRAME = re.compile(r"^\t(at |\.\.\. \d+ more)")


def _normalise_exceptions(text: str) -> str:
    """Removes stack frames and reduces a re-thrown exception to its
    innermost cause."""
    lines = text.split("\n")
    blocks: List[List[str]] = []   # exception headers with message lines
    out: List[str] = []
    current = None
    for line in lines:
        m = _EXC_HEADER.match(line)
        if m:
            current = [m.group(1), m.group(2)]
            blocks.append(current)
            out.append(("EXC", len(blocks) - 1))
            continue
        if _FRAME.match(line):
            current = None
            continue
        if current is not None:
            current[1] += "\n" + line
            continue
        out.append(line)
    result: List[str] = []
    main_done = False
    for item in out:
        if isinstance(item, tuple):
            index = item[1]
            header, body = blocks[index]
            if header.startswith("Caused by"):
                continue
            if main_done:
                continue
            main_done = True
            # The innermost cause, if the exception was re-thrown.
            causes = [b for b in blocks[index + 1:]
                      if b[0].startswith("Caused by")]
            if causes:
                body = causes[-1][1]
            else:
                # IllegalArgumentException(Throwable): "X: X: message".
                m = re.match(r"^([\w.$]+): \1: ", body)
                while m:
                    body = body[len(m.group(1)) + 2:]
                    m = re.match(r"^([\w.$]+): \1: ", body)
            result.append('Exception in thread "main" ' + body)
        else:
            result.append(item)
    return "\n".join(result)


_SEGMENT_LINE = re.compile(r"^[^\t\n]+(\t[^\t\n]*){8}$", re.M)


def normalise_failure(stdout: str, stderr: str) -> str:
    stdout = _SEGMENT_LINE.sub("", stdout)
    stdout = re.sub(r"\n{2,}", "\n", stdout) if stdout.strip() == "" else stdout
    text = "STDOUT:\n" + stdout + "\nSTDERR:\n" + stderr
    text = _collapse_usage(text)
    text = _normalise_exceptions(text)
    text = re.sub(r"@[0-9a-f]{4,}", "", text)
    text = re.sub(r"\n+", "\n", text)
    return text.strip()


LARGE = 2.0 ** 53
KNOWN_OUT_OF_RANGE_CASE = "ibd_absent_sample_gap_test"


def _only_large_numbers_differ(only_j: collections.Counter,
                               only_p: collections.Counter) -> bool:
    """True if the differing lines pair up by their first seven fields and
    differ only in the score and length columns, with every differing value
    at least 2**53 in magnitude in both outputs."""
    def by_key(counter):
        out = {}
        for line in counter.elements():
            fields = line.split("\t")
            if len(fields) != 9:
                return None
            key = tuple(fields[:7])
            if key in out:
                return None
            out[key] = fields
        return out
    j = by_key(only_j)
    p = by_key(only_p)
    if j is None or p is None or set(j) != set(p):
        return False
    for key, jf in j.items():
        pf = p[key]
        for col in (7, 8):
            if jf[col] != pf[col]:
                try:
                    if abs(float(jf[col])) < LARGE or abs(float(pf[col])) < LARGE:
                        return False
                except ValueError:
                    return False
    return True


def compare(java, py, name: str = "") -> Tuple[str, str]:
    j_code, j_out, j_err, _ = java
    p_code, p_out, p_err, _ = py
    if j_code != p_code:
        return "MISMATCH", "exit status java=%d python=%d" % (j_code, p_code)
    if j_code == 0:
        j_lines = collections.Counter(_collapse_usage(j_out).splitlines())
        p_lines = collections.Counter(_collapse_usage(p_out).splitlines())
        if j_lines != p_lines:
            only_j = j_lines - p_lines
            only_p = p_lines - j_lines
            detail = ("output differs: %d line(s) only in java, %d only in "
                      "python" % (sum(only_j.values()), sum(only_p.values())))
            if (not (j_err.strip() or p_err.strip())
                    and _only_large_numbers_differ(only_j, only_p)):
                return "KNOWN", detail + ", all values of 2**53 or more"
            return "MISMATCH", detail
        if j_err.strip() or p_err.strip():
            if normalise_failure("", j_err) != normalise_failure("", p_err):
                return "MISMATCH", "standard error differs"
        return "MATCH", "%d output lines" % sum(j_lines.values())
    jn = normalise_failure(j_out, j_err)
    pn = normalise_failure(p_out, p_err)
    if jn != pn:
        oob = 'Exception in thread "main" java.lang.IndexOutOfBoundsException'
        if (name == KNOWN_OUT_OF_RANGE_CASE and j_code == 1
                and jn.count(oob) == 1 and pn.count(oob) == 1):
            return "KNOWN", ("both report IndexOutOfBoundsException; the jar "
                             "reads outside its genotype arrays")
        return "MISMATCH", "diagnostic differs"
    return "MATCH", "exit %d, diagnostic equal" % j_code


def write_case(case: Case, directory: str) -> None:
    os.makedirs(directory, exist_ok=True)
    for name, data in case.files.items():
        with open(os.path.join(directory, name), "wb") as f:
            f.write(data)
    with open(os.path.join(directory, "stdin.ibd"), "wb") as f:
        f.write(case.stdin)


def environment(jar: str, java: str) -> Dict[str, str]:
    jv = subprocess.run([java, "-version"], stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT).stdout.decode()
    props = subprocess.run(
        [java, "-XshowSettings:properties", "-version"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.decode()
    keep = [line.strip() for line in props.splitlines()
            if re.match(r"\s*(file\.encoding|native\.encoding|stdout\.encoding"
                        r"|stderr\.encoding|sun\.jnu\.encoding|user\.language"
                        r"|user\.country|java\.version|java\.vendor|"
                        r"java\.runtime\.version|os\.name|os\.version|"
                        r"os\.arch) =", line)]
    return {
        "date": datetime.datetime.now(datetime.timezone.utc).isoformat(
            timespec="seconds"),
        "jar": os.path.basename(jar),
        "jar_sha256": sha256(jar),
        "source_zip_sha256": SOURCE_ZIP_SHA256.lower(),
        "java_version": jv.strip(),
        "java_properties": "\n".join(keep),
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "locale_env": " ".join("%s=%s" % (k, os.environ.get(k, ""))
                               for k in ("LANG", "LC_ALL", "LC_CTYPE")),
        "python_sources_sha256": "\n".join(
            "%s  %s" % (sha256(os.path.join(REPO, "merge_ibs_segments", n)), n)
            for n in sorted(os.listdir(os.path.join(REPO, "merge_ibs_segments")))
            if n.endswith(".py")),
        "test_sources_sha256": "\n".join(
            "%s  %s" % (sha256(os.path.join(HERE, n)), n)
            for n in sorted(os.listdir(HERE)) if n.endswith(".py")),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--jar", required=True)
    parser.add_argument("--java", default="java")
    parser.add_argument("--java-opt", action="append", default=[],
                        help="option for the Java virtual machine, for "
                             "example -Xmx768m (repeatable)")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--work", required=True,
                        help="directory for generated inputs and outputs "
                             "(emptied first)")
    parser.add_argument("--report", required=True, help="Markdown report")
    parser.add_argument("--random", type=int, default=20,
                        help="number of seeded random cases")
    parser.add_argument("--only", default="", help="regular expression of "
                        "case names to run")
    args = parser.parse_args()

    jar = os.path.abspath(args.jar)
    work = os.path.abspath(args.work)
    if os.path.isdir(work):
        shutil.rmtree(work)
    os.makedirs(work)
    env_info = environment(jar, args.java)
    env_info["java_options"] = " ".join(args.java_opt) or "(none)"
    py_env = dict(os.environ)
    py_env["PYTHONPATH"] = REPO + os.pathsep + py_env.get("PYTHONPATH", "")

    cases = all_cases(args.random)
    if args.only:
        cases = [c for c in cases if re.search(args.only, c.name)]
    results = []
    for case in cases:
        d = os.path.join(work, case.name)
        write_case(case, d)
        java = run([args.java] + args.java_opt + ["-jar", jar] + case.args,
                   d, case.stdin)
        p = subprocess.run([args.python, "-m", "merge_ibs_segments"]
                           + case.args, cwd=d, input=case.stdin,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           env=py_env, timeout=600)
        py = (p.returncode, p.stdout.decode("utf-8", "replace"),
              p.stderr.decode("utf-8", "replace"), 0.0)
        status, detail = compare(java, py, case.name)
        for label, r in (("java", java), ("python", py)):
            with open(os.path.join(d, label + ".stdout"), "w",
                      encoding="utf-8") as f:
                f.write(r[1])
            with open(os.path.join(d, label + ".stderr"), "w",
                      encoding="utf-8") as f:
                f.write(r[2])
        results.append({"case": case.name, "status": status,
                        "detail": detail, "java_exit": java[0],
                        "python_exit": py[0], "note": case.note})
        print("%-40s %-8s %s" % (case.name, status, detail), flush=True)

    n_match = sum(r["status"] == "MATCH" for r in results)
    n_known = sum(r["status"] == "KNOWN" for r in results)
    n_mismatch = len(results) - n_match - n_known
    with open(os.path.join(work, "results.json"), "w", encoding="utf-8") as f:
        json.dump({"environment": env_info, "results": results}, f, indent=1)
    lines = ["# Exact-conversion test record", "",
             "Synthetic inputs only (generated by `tests/parity/cases.py`).",
             "", "## Environment", ""]
    for key, value in env_info.items():
        if "\n" in value:
            lines += ["- %s:" % key, "", "```", value, "```", ""]
        else:
            lines.append("- %s: `%s`" % (key, value))
    lines += ["", "## Result", "",
              "%d of %d cases match; %d show only a known, accepted "
              "difference; %d mismatch." % (n_match, len(results), n_known,
                                            n_mismatch), "",
              "Known, accepted differences (see run_parity.py): numbers of "
              "2**53 and above are formatted with other digits by the jar; a "
              "gap test for a sample absent from the VCF file makes the jar "
              "read outside its genotype arrays.", "",
              "| Case | Result | Detail | Java exit | Python exit |",
              "|---|---|---|---|---|"]
    for r in results:
        lines.append("| %s | %s | %s | %d | %d |" % (
            r["case"], r["status"], r["detail"], r["java_exit"],
            r["python_exit"]))
    with open(args.report, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print("%d of %d cases match, %d known, %d mismatch"
          % (n_match, len(results), n_known, n_mismatch))
    return 0 if n_mismatch == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
