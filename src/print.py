import sys
import time

def print_error(module, msg, exit=True):
    """Print an error message and terminate.

    Parameters
    ----------
        module : str
            Which module generated the error message
        msg : str
            The error message itself
    """
    print("[ERROR] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr)
    if exit:
        exit(1)


def print_info(module, msg):
    print("[INFO] ({}) {}: {}".format(time.strftime("%b %d, %Y %X"), module, msg), file=sys.stderr)