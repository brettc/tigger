import logging
log = logging.getLogger("tigger.main")
logging.basicConfig(level=logging.INFO)

import sys
from tigger.alignment import Alignment
from tigger._bounce import TigerDNA
from pathlib import Path
import numpy

if __name__ == "__main__":
    a = Alignment()
    filepath = Path(sys.argv[1])
    log.info("Reading file %s", str(filepath))

    a.read(str(filepath))
    log.info("Species count %s", a.species_count)
    log.info("Alignment size %s", a.sequence_length)

    log.info("Beginning analysis --- ")
    t = TigerDNA()
    t.build_bitsets(a)
    rates = t.calc_rates()
    output_path = str(filepath.with_suffix('.tigger'))

    log.info("Finished analysis --- ")

    log.info("Saving file %s", output_path)
    numpy.savetxt(output_path, rates, fmt="%5f", delimiter='\n')



    

