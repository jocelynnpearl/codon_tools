from Bio.SeqUtils import CodonUsage

from . import Seq, logging
from ..data.codon_use_by_host import codon_tables, latin_name_to_taxID

logger = logging.getLogger(__name__)

# returns dictionary, counts the number of times each codon is used
def count_codons(dna_sequence):
    logger.info("===== COUNTING CODONS =====")

    codons_dict = CodonUsage.CodonsDict.copy()
    for codon in dna_sequence.codons():
        codons_dict[codon] += 1

    return codons_dict


# returns dictionary, calculates the % usage of each AA's codons
def calc_profile(codons_count):
    logger.info("===== CALCULATING PROFILE =====")
    codons_freq = CodonUsage.CodonsDict.copy()

    # loop through all amino acids
    for _, synonymous_codons in CodonUsage.SynonymousCodons.items():
        # add up number of times each AA appears
        tot_usage = sum([codons_count[codon] for codon in synonymous_codons])

        # calculate the distribution of synonymous codon use
        if tot_usage == 0:
            continue
        codons_freq.update(
            {codon: codons_count[codon] / tot_usage for codon in synonymous_codons}
        )

    return codons_freq


def _load_host_table(host):
    logger.info("===== PROCESSING HOST TABLE: {0} =====".format(host))
    table = CodonUsage.CodonsDict.copy()

    if host not in codon_tables and host in latin_name_to_taxID:
        host = latin_name_to_taxID[host]

    try:
        raw_table = codon_tables[host]
    except KeyError:
        raise ValueError(
            '"{}" is not a valid host id. '.format(host)
            + "Currently supported hosts (Latin and NCBI taxonomy IDs) are {}".format(
                host, ", ".join(latin_name_to_taxID.keys() + codon_tables.keys())
            )
        )
    else:
        for codon, usage in raw_table.items():
            # codons are stored with RNA alphabet
            table[str(Seq(codon).back_transcribe())] = usage["count"]

    return calc_profile(table)


# returns a calc_profile for a given host table
def process_host_table(host, threshold=0.10):

    table = _load_host_table(host)

    logger.info("HOST THRESHOLD SET TO: {0}".format(threshold))
    logger.detail("pre-threshold host table:")

    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        logger.detail("== {0} ==".format(AA))
        for codon in synonymous_codons:
            logger.detail("{0}: {1}".format(codon, table[codon]))
            if table[codon] < threshold:
                table[codon] = 0

    # recalculate profile after dropping rare codons
    table = calc_profile(table)

    if logger.isEnabledFor(logging.DETAIL):
        logger.detail("post-threshold host table:")
        for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
            logger.detail("== {0} ==".format(AA))
            for codon in synonymous_codons:
                logger.detail("{0}: {1}".format(codon, table[codon]))

    return table


def host_codon_usage(host):

    host_profile = _load_host_table(host)

    codon_use_by_aa = {}
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        codon_use_by_aa[AA] = [
            synonymous_codons,
            [host_profile[codon] for codon in synonymous_codons],
        ]

    return codon_use_by_aa
