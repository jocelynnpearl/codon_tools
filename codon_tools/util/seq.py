import Bio

from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import CodonUsage, seq3


# reverse translate AA seq to DNA seq
def back_translate(self):
    """Return the DNA sequence from an amino acid sequence by creating a new Seq object.

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_protein = Seq("MAIVMGR", IUPAC.protein)
    >>> my_protein
    Seq('MAIVMGR', IUPACProtein())
    >>> my_protein.back_translate()
    Seq('ATGGCCATTGTAATGGGCCGCTG', IUPACUnambiguousDNA())

    Trying to back-transcribe a DNA or RNA sequence raises an
    exception:

    >>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUG", IUPAC.unambiguous_rna)
    >>> messenger_rna.back_translate()
    Traceback (most recent call last):
    ...
    ValueError: Nucleic acids cannot be back translated!
    """
    base = Bio.Alphabet._get_base_alphabet(self.alphabet)
    if not isinstance(base, Bio.Alphabet.ProteinAlphabet):
        raise ValueError("Nucleic acids cannot be back translated!")

    # right now this just uses the most-prevalent codon for each AA
    # TODO: select codons with a weighted average using random.choice
    return Seq(
        "".join([CodonUsage.SynonymousCodons[seq3(AA).upper()][0] for AA in str(self)]),
        IUPAC.unambiguous_dna,
    )


Seq.back_translate = back_translate


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    from itertools import zip_longest

    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def codons(self):
    base = Bio.Alphabet._get_base_alphabet(self.alphabet)
    if isinstance(base, Bio.Alphabet.ProteinAlphabet):
        raise ValueError("Protein sequences cannot be split into codons!")
    for codon in grouper(self, 3):
        yield "".join(codon)


Seq.codons = codons
