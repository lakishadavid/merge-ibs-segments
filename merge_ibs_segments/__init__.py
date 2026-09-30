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
"""merge-ibs-segments: a Python conversion of Brian Browning's
merge-ibd-segments (Refined-IBD version 17Jan20.102).

Modules: ``merge`` (the program), ``segments`` (IBD segments), ``vcf`` (VCF
reading), ``genetic_map`` (PLINK map), ``ids`` (identifier indexes),
``input_it`` (line input) and ``java_compat`` (Java behaviour the others
depend on).  The notice at the top of each module names the Java files it
was translated from.
"""
