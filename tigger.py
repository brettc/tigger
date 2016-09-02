"""
Tigger is an experiment in hybrid programming for phylogenetics. It implements
the TIGER algorithm.

Usage:
  tigger [--debug] <alignment_file>
  tigger --version


Options:
  -h --help     Show this screen.
  --version     Show version.
"""

__version__ = '1.0.2'

import logging
log = logging.getLogger("tigger.main")
logging.basicConfig(
    format="%(levelname)-8s | %(asctime)s | %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S"
)

from docopt import docopt
import sys
from tigger.alignment import Alignment, AlignmentError
from tigger._bounce import TigerDNA
from pathlib import Path
import numpy

def set_verbose():
    logging.getLogger("").setLevel(logging.DEBUG)
    # Enhance the format
    fmt = logging.Formatter(
        "%(levelname)-8s | %(asctime)s | %(name)-20s | %(message)s")
    logging.getLogger("").handlers[0].setFormatter(fmt)

def main(arguments):
    if arguments['--debug']:
        set_verbose()

    filepath = Path(arguments['<alignment_file>'])

    try:
        a = Alignment(filepath)
    except AlignmentError:
        return 1

    log.info("Species count %s", a.species_count)
    log.info("Alignment size %s", a.sequence_length)

    log.info("Beginning analysis --- ")
    t = TigerDNA()
    t.build_bitsets(a)
    rates = t.calc_rates()
    log.info("Finished analysis --- ")

    output_path = str(filepath.with_suffix('.tigger'))
    log.info("Saving file %s", output_path)
    numpy.savetxt(output_path, rates, fmt="%5f", delimiter='\n')

    # Everything ran okay
    return 0

if __name__ == "__main__":
    arguments = docopt(__doc__, version=__version__)
    sys.exit(main(arguments))



    

