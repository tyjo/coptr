"""
print.py
======================
Print error, warning, and info messages.
"""

import sys
import time

def print_error(module, msg, exit=True):
    """Print an error message to the standard error.

    Parameters
    ----------
        module : str
            Which module generated the error message
        msg : str
            The error message itself
        exit :  bool
            If TRUE, the program terminates after printing the error
    """
    print("[ERROR] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr)
    if exit:
        exit(1)


def print_warning(module, msg):
    """Print a warning message to the standard error.
    
    Parameters
    ----------
        module : str
            Which module generated the message
        msg : str
            The message itself
    """

    print("[WARNING] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr)


def print_info(module, msg):
    """Print an info message to the standard error.
    
    Parameters
    ----------
        module : str
            Which module generated the message
        msg : str
            The message itself
    """

    print("[INFO] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr)