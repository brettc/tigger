"""Loading, Saving, Parsing Alignment Files

    See the phyml details here:
    http://www.atgc-montpellier.fr/phyml/usersguide.php?type=command
"""
import logging
log = logging.getLogger("tigger.alignment")

import os
from pyparsing import (Word, OneOrMore, alphas, nums, Suppress, stringEnd,
                       ParseException, restOfLine, LineEnd, ZeroOrMore, Upcase)
import numpy
import array
import pathlib

class AlignmentError(RuntimeError):
    pass

class BaseParser(object):
    """Parses an alignment and returns species sequence tuples"""

    # I think this covers it...
    BASES = Upcase(Word(alphas + "?.-"))

    def __init__(self):
        self.sequence_length = None
        self.species_count = None
        self.sequences = []
        self.current_sequence = 0

        self.root_parser = self.init_parser() + stringEnd

    def parse(self, s):
        try:
            self.root_parser.parseString(s)

        except ParseException as p:
            log.error("Error in Alignment Parsing:" + str(p))
            log.error("A common cause of this error is having whitespace"
                      ", i.e. spaces or tabs, in the species names. Please check this and remove"
                      " all whitespace from species names, or replace them with e.g. underscores")

            raise AlignmentError

        # Check that all the sequences are equal length
        slen = None
        names = set()
        for nm, seq in self.sequences:
            if nm in names:
                log.error("Repeated species name '%s' is repeated "
                          "in alignment", nm)
                raise AlignmentError

            names.add(nm)
            if slen is None:
                # Use the first as the test case
                slen = len(seq)
            else:
                if len(seq) != slen:
                    log.error(
                        "Bad alignment file: Not all species have the same sequences length")
                    raise AlignmentError

        # Not all formats have a heading, but if we have one do some checking
        if self.sequence_length is None:
            self.sequence_length = len(self.sequences[0][1])
        else:
            if self.sequence_length != slen:
                log.error("Bad Alignment file: sequence length count in header does not match"
                          " sequence length in file, please check")
                raise AlignmentError

        if self.species_count is None:
            self.species_count = len(self.sequences)
        else:
            if len(self.sequences) != self.species_count:
                log.error("Bad Alignment file: species count in header does not match"
                          " number of sequences in file, please check")
                raise AlignmentError


class PhylipParser(BaseParser):

    def init_parser(self):

        INTEGER = Word(nums)
        INTEGER.setParseAction(lambda x: int(x[0]))

        header = INTEGER("species_count") + INTEGER("sequence_length") +\
            Suppress(restOfLine)
        header.setParseAction(self.set_header)

        sequence_name = Word(
            alphas + nums + "!#$%&\'*+-./;<=>?@[\\]^_`{|}~",
            max=100)

        # Take a copy and disallow line breaks in the bases
        bases = self.BASES.copy()
        bases.setWhitespaceChars(" \t")
        seq_start = sequence_name("species") + bases(
            "sequence") + Suppress(LineEnd())
        seq_start.setParseAction(self.set_seq_start)
        seq_start_block = OneOrMore(seq_start)
        seq_start_block.setParseAction(self.set_start_block)

        seq_continue = bases("sequence") + Suppress(LineEnd())
        seq_continue.setParseAction(self.set_seq_continue)

        seq_continue_block = Suppress(LineEnd()) + OneOrMore(seq_continue)
        seq_continue_block.setParseAction(self.set_continue_block)

        return header + seq_start_block + ZeroOrMore(seq_continue_block)

    def set_header(self, text, loc, tokens):
        self.sequence_length = tokens.sequence_length
        self.species_count = tokens.species_count

    def set_seq_start(self, text, loc, tokens):
        self.sequences.append([tokens.species, tokens.sequence])
        self.current_sequence += 1

    def set_start_block(self, tokens):
        # End of block
        # Reset the counter
        self.current_sequence = 0

    def set_seq_continue(self, text, loc, tokens):
        append_to = self.sequences[self.current_sequence]
        append_to[1] += tokens.sequence
        self.current_sequence += 1

    def set_continue_block(self, tokens):
        self.current_sequence = 0

class FastaParser(BaseParser):

    def init_parser(self):
        sequence_name = Word(
            alphas + nums + "!#$%&\'*+-./;?@[\\]^_`{|}~",
            max=100)

        sequence_name.setParseAction(self.set_name)
        name_block = Suppress(">") + sequence_name("name") + Suppress(LineEnd())

        # Take a copy and disallow line breaks in the bases
        bases = self.BASES.copy()
        seq = bases("sequence") + Suppress(LineEnd())
        seq.setParseAction(self.set_sequence)

        name_and_seq = name_block + seq
        return OneOrMore(name_and_seq)

    def set_name(self, tokens):
        self.current_name = tokens.name

    def set_sequence(self, tokens):
        self.sequences.append((self.current_name, tokens.sequence))


class Alignment(object):
    def __init__(self):
        self.species = []
        self.sequence_length = 0
        self.data = None

    @property
    def species_count(self):
        return len(self.species)

    def __str__(self):
        return "Alignment(%s species, %s codons)"\
               % (self.species_count, self.sequence_length)

    def parse(self, text, p):
        """Parse the sequence, then transfer data from the parser
        Note: parser returns tuples like ("dog", "GATC"), ("cat", "GATT")
        """
        # Allocate a numpy array using unsigned char
        p.parse(text)
        d = numpy.zeros((p.species_count, p.sequence_length), 'u1')
        self.sequence_length = p.sequence_length
        for i, (spec, codons) in enumerate(p.sequences):
            self.species.append(spec)
            # TODO: Is there a better way to do this directly into numpy?
            # Typecode "B" makes a unsigned int
            d[i] = array.array("B", codons)
        self.data = d

    def read(self, pth):
        path = pathlib.Path(pth)
        if not path.exists():
            log.error("Cannot find sequence file '%s'", pth)
            raise AlignmentError

        log.debug("Reading alignment file '%s'", pth)
        text = open(str(path), 'rU').read()
        suff = path.suffix.lower()
        if suff == '.phy':
            parser = PhylipParser()
        elif suff == '.fas':
            parser = FastaParser()
        else:
            log.error("Unknown file type: %s", str(path))
            raise AlignmentError

        self.parse(text, parser)


