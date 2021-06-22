"""
print.py
======================
Print error, warning, and info messages.
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

import sys
import time


def print_error(module, msg, quit=True):
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
    print("[ERROR] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr, flush=True)
    if quit:
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

    print("[WARNING] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr, flush=True)


def print_info(module, msg):
    """Print an info message to the standard error.
    
    Parameters
    ----------
        module : str
            Which module generated the message
        msg : str
            The message itself
    """

    print("[INFO] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr, flush=True)