import sys
from pkg_resources import parse_version

def check_requirements():
    if sys.version_info[0] < 3:
        raise Exception("Must be using Python 3")

    try:
        import matplotlib
        if parse_version(matplotlib.__version__) < parse_version('2.2.3'):
            raise Exception("Matplotlib requires version 2.2.3 or newer")
    except ModuleNotFoundError:
        raise Exception("Matplotlib not installed")

    try:
        import numpy
        if parse_version(numpy.__version__) < parse_version('1.15.1'):
            raise Exception("NumPy requires version 1.15.1 or newer")
    except ModuleNotFoundError:
        raise Exception("NumPy not installed")

    try:
        import scipy
        if parse_version(scipy.__version__) < parse_version('1.1.0'):
            raise Exception("SciPy requires version 1.1.0 or newer")
    except ModuleNotFoundError:
        raise Exception("SciPy not installed")
