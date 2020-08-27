"""
util.py
===============================
Miscellaneous utility functions.
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
    bn,ex = os.path.splitext(fpath)
    if ex == ".gz":
        bn,ex = os.path.splitext(bn)
    return bn