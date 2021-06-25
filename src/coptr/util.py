"""
util.py
===============================
Miscellaneous utility functions.
"""

"""
This file is part of CoPTR.

CoPTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CoPTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CoPTR.  If not, see <https://www.gnu.org/licenses/>.
"""

import os.path


def get_fastq_name(fpath):
    """Remove the path and extension from a fastq file.

    Parameters
    ----------
        fpath : str
            path to .fastq, fastq.gz, .fq, or fq.gz file
    Returns
    -------
        fname : str
            the name of the file with its path and extension removeed
    """
    bn, ex = os.path.splitext(fpath)
    if ex == ".gz":
        bn, ex = os.path.splitext(bn)
    return bn
