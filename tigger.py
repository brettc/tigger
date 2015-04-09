import sys
from tigger.alignment import Alignment
from tigger._bounce import TigerDNA
from pathlib import Path
import numpy

if __name__ == "__main__":
    a = Alignment()
    filepath = Path(sys.argv[1])
    a.read(str(filepath))

    tigger = TigerDNA()
    tigger.build_bitsets(a)
    rates = tigger.calc_rates()
    output_path = str(filepath.with_suffix('.tigger'))
    numpy.savetxt(output_path, rates, fmt="%5f", delimiter='\n')


    

