#! /usr/bin/env python
import argparse
import random
import re
import sys

import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import CodonUsage, GC, seq3

from codon_tools.util import Seq, MutableSeq, codon_use, logging, log_levels
from codon_tools.data import GC_content, RestrictionEnz, RibosomeBindingSites


# options
parser = argparse.ArgumentParser(
    description="Optimize your AA or DNA sequence to harmonize with a host's condon usage.",
    epilog="2017-12-04, v0.47 (contact yhsia@uw.edu if stuff does not work)",
)
parser.add_argument("--input", type=str, required=True, help="input file with sequence")
parser.add_argument(
    "--host",
    type=str,
    default="413997",
    help='host table code: http://www.kazusa.or.jp/codon/, default is "Escherichia coli B"',
)
parser.add_argument(
    "--host_threshold",
    type=float,
    default="0.10",
    help="lowest codon fraction per AA in the host that is allowed",
)
parser.add_argument(
    "--verbose",
    type=int,
    default=0,
    choices=[0, 1, 2, 3],
    help="verbose output level (0=only result, 1=standard output, 2=extra output 3=debugging)",
)
parser.add_argument(
    "--local_homopolymer_threshold",
    type=int,
    default="4",
    help="number of consecutive NT repeats allowed",
)
parser.add_argument(
    "--cycles",
    type=int,
    default=1000,
    help="max number of cycles to run optimization, 0=unlimited",
)
parser.add_argument(
    "--max_relax",
    type=float,
    default="0.1",
    help="maximum % deviation from host profile",
)

args = parser.parse_args()


logging.basicConfig(level=log_levels[args.verbose])
logger = logging.getLogger(__name__)
random.seed()


# returns dictionary of comparison between two profiles
def compare_profiles(codons_count, host, relax):
    logger.info("===== COMPARING PROFILES =====")
    table = {}
    # loop AAs
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        logger.detail(AA)
        temp_table = {}

        # calculate total usage of codon in input
        tot_usage = sum([codons_count[codon] for codon in synonymous_codons])

        # calculate ideal usage of codon in host
        tot_ideal = 0
        for codon in synonymous_codons:
            ideal_usage_abs = int(round(host[codon] * tot_usage, 0))
            ideal_usage = int(round(host[codon] * relax * tot_usage, 0))
            logger.detail("{0}: {1}".format(codon, ideal_usage))
            tot_ideal += ideal_usage
            temp_table[codon] = {
                "input_count": codons_count[codon],
                "input_perc": codons_count[codon],
                "ideal_usage_abs": ideal_usage_abs,
                "ideal_usage": ideal_usage,
                "host_perc": float(host[codon]),
            }

        # if ideal number is too high, subtract from lowest host codon (that is not 0 already )
        while tot_ideal > tot_usage:
            lowest_table = temp_table
            lowest_codon = min(
                temp_table, key=lambda k: float(temp_table[k]["host_perc"])
            )
            while temp_table[lowest_codon]["ideal_usage"] == 0:
                lowest_table = removekey(lowest_table, lowest_codon)
                lowest_codon = min(
                    lowest_table, key=lambda k: float(lowest_table[k]["host_perc"])
                )

            temp_table[lowest_codon]["ideal_usage"] -= 1
            tot_ideal -= 1

        # if ideal number is too low, add to highest host codon
        while tot_ideal < tot_usage:
            highest_codon = max(
                temp_table, key=lambda k: float(temp_table[k]["host_perc"])
            )
            temp_table[highest_codon]["ideal_usage"] = (
                temp_table[highest_codon]["ideal_usage"] + 1
            )
            tot_ideal += 1

        # populate return table
        for codon in synonymous_codons:
            table[codon] = {
                "input_count": temp_table[codon]["input_count"],
                "ideal_usage_abs": temp_table[codon]["ideal_usage_abs"],
                "difference": temp_table[codon]["input_count"]
                - temp_table[codon]["ideal_usage"],
                "difference_abs": temp_table[codon]["input_count"]
                - temp_table[codon]["ideal_usage_abs"],
            }

    # calculate difference value
    resi_total = 0
    diff_total = 0
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        for codon in synonymous_codons:
            resi_total += table[codon]["ideal_usage_abs"]
            diff_total += abs(table[codon]["difference_abs"])

    logger.info(
        "CURRENT DIFFERENCE TOTAL: {0} of {1}".format(int(diff_total / 2), resi_total)
    )

    logger.info("CURRENT DIFFERENCE %: {0}".format(diff_total / resi_total / 2))
    diff = diff_total / resi_total / 2

    return table, diff


# returns mutated sequence after given a difference profile
def resample_codons(dna_sequence, codon_use_by_aa):
    """[summary]

    Args:
        dna_sequence ([type]): [description]
        codon_use_by_aa ([type]): [description]

    Returns:
        [type]: [description]
    """

    resampled_dna = "".join(
        [
            random.choices(*codon_use_by_aa[seq3(AA).upper()]).pop()
            for AA in dna_sequence.translate()
        ]
    )

    return Seq(resampled_dna, IUPAC.unambiguous_dna)


# check for local homopolymers
def remove_local_homopolymers(dna_sequence, n_codons=2):
    logger.info("===== REMOVE LOCAL HOMOPOLYMERS =====")
    mutable_seq = dna_sequence.tomutable()

    # look at each 6-mer
    keep_looping = True
    while keep_looping:
        for i in range(0, len(mutable_seq), 3):
            window = slice(
                i,
                i + (n_codons * 3)
                if i + (n_codons * 3) < len(mutable_seq)
                else len(mutable_seq),
            )

            seq = str(mutable_seq[window])
            nt_counts = {letter: seq.count(letter) for letter in set(seq)}
            letter = max(nt_counts, key=lambda letter: nt_counts[letter])

            if nt_counts[letter] <= args.local_homopolymer_threshold:
                keep_looping = False
                continue

            logger.detail("position: {0}: {1}".format(i, seq))
            logger.detail("{0}, count={1}".format(letter, nt_counts[letter]))

            for j in range(n_codons):
                codon_idx = slice(i + (j * 3), i + ((j + 1) * 3))
                mutable_seq[codon_idx] = mutate_codon(
                    mutable_seq[codon_idx], codon_use_table
                )
            keep_looping = True

    return mutable_seq.toseq()


# check for unwanted restriction sites
def remove_restriction_sites(dna_sequence, restrict_sites):
    logger.info("===== REMOVE RESTRICTION SITES =====")

    # check each unwanted restriction site
    analysis = Analysis(restrictionbatch=restrict_sites, sequence=dna_sequence)
    result = analysis.full()

    mutable_seq = dna_sequence.tomutable()
    for enz, cuts in result.items():
        for cut in cuts:
            logger.info(
                "The restriction enzyme {0} can cut the sequence before position {1}!".format(
                    str(enz), cuts
                )
            )
            # map sequence position to codon position
            # subtract 1 from `cut` to convert from sequence to string indexing
            codon_pos, offset = divmod((cut - 1) - (enz.size // 2), 3)

            # ensure the whole codon we mutate is being recognized by the restriction enzyme
            if offset:
                codon_pos += 1
            codon_idx = slice(codon_pos * 3, (codon_pos + 1) * 3)

            new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
            mutable_seq[codon_idx] = new_codon

    return mutable_seq.toseq()


# check for alternative start sites
def remove_start_sites(dna_sequence, ribosome_binding_sites, table_name="Standard"):

    codon_table = CodonTable.unambiguous_dna_by_name[table_name]
    logger.info(
        "===== REMOVE START SITES: {0} =====".format(
            ", ".join(codon_table.start_codons)
        )
    )

    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(dna_sequence))
    ]

    if not len(start_codon_positions):
        logger.info("No start codons found in sequence")
        return dna_sequence

    logger.info(
        "Found {0} start codons. Checking for upstream RBSs...".format(
            len(start_codon_positions)
        )
    )

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    _rbs_offset = 18
    rbs_positions = [
        pos - _rbs_offset for pos in start_codon_positions if pos >= _rbs_offset
    ]
    mutable_seq = dna_sequence.tomutable()

    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + _rbs_offset - 3
        rbs_query_seq = str(mutable_seq[rbs_start:rbs_stop])

        logger.detail(
            "checking sequence: {0}.{1}".format(
                rbs_query_seq, mutable_seq[rbs_stop : rbs_stop + 6]
            )
        )

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            logger.detail("checking start site: {0}, {1}".format(rbs, site))
            search = rbs_query_seq.find(site)

            count = 0  # counter to prevent infinite loop
            while search != -1 and count < 10:
                # mutate residues if site is found
                codon_pos = (search + rbs_start) // 3
                for i in range(2):
                    codon_idx = slice((codon_pos + i) * 3, (codon_pos + i + 1) * 3)
                    new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
                    mutable_seq[codon_idx] = new_codon

                # reset sequence and search again
                rbs_query_seq = str(mutable_seq[rbs_start : rbs_stop + 3])
                search = rbs_query_seq.find(site)
                count += 1

    return mutable_seq.toseq()


# mutate single codon
def mutate_codon(codon_in, codon_use_table):
    AA = seq3(CodonTable.standard_dna_table.forward_table[str(codon_in)]).upper()

    synonymous_codons, codon_use_freq = codon_use_table[AA]
    if len(synonymous_codons) == 1:
        return codon_in

    # pick new codon
    codon_out = codon_in
    while codon_in == codon_out:
        codon_out = random.choices(synonymous_codons, codon_use_freq).pop()

    logger.detail(
        "mutating [{0}] codon from {1} to {2}".format(AA, codon_in[1], codon_out)
    )

    return codon_out


# check GC content
def gc_scan(dna_sequence, window_size, low, high):
    logger.info(
        "===== GC CONTENT SCAN IN WINDOW: {0} bps, threshold: {1} < x < {2}=====".format(
            window_size, low, high
        )
    )

    # some windows may be expressed as function of the sequence length
    # for example, "x0.5"
    if isinstance(window_size, str) and window_size.startswith("x"):
        window_size = int(float(window_size[1:]) * len(dna_sequence))

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    overlap = codon_window // 2
    mutable_seq = dna_sequence.tomutable()

    # iterate by codon, but map back to sequence-based indices
    for i in range(0, len(mutable_seq) // 3, (codon_window - overlap) * 3):
        window = slice(
            i * 3,
            (i + codon_window) * 3
            if (i + codon_window) * 3 < len(mutable_seq)
            else len(mutable_seq),
        )
        logger.debug("Current segment: {0}".format(mutable_seq[window]))

        gc_percent = GC(mutable_seq[window]) / 100
        count = 0  # counter to prevent infinite loop
        # check gc_percent of current segment
        while (gc_percent < low or gc_percent > high) and count < codon_window * 2:
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice((i * 3) + position, ((i + 1) * 3) + position)

            init_codon = mutable_seq[codon_idx]
            new_codon = mutate_codon(init_codon, codon_use_table)

            if (GC(new_codon) < GC(init_codon) and gc_percent > high) or (
                GC(new_codon) > GC(init_codon) and gc_percent < low
            ):
                mutable_seq[codon_idx] = new_codon
                logger.debug("Mutating position: {0}".format(position))
                gc_percent = GC(mutable_seq[window]) / 100

            count += 1

    return mutable_seq.toseq()


# check for repeat segments
def remove_repeating_sequences(dna_sequence, window_size):
    logger.info(
        "===== REPEAT FRAGMENT SCAN FOR SIZE: {0} bps =====".format(window_size)
    )

    def _mutate_and_keep_looping(mutable_seq, window, offset):
        num_mutations = random.randint(1, 2)
        logger.debug("Mutating {0} codons".format(num_mutations))
        for _ in range(num_mutations):
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice(offset + position, (offset + 3) + position)
            new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
            mutable_seq[codon_idx] = new_codon

        return True

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    overlap = codon_window - 1
    mutable_seq = dna_sequence.tomutable()

    current_cycle = 0  # prevent infinite loops (caused by poly-TRP or poly-MET)
    keep_looping = True
    # `keep_looping` if any mutation is made,
    # i.e. continue until both checks pass without mutations
    while keep_looping and (current_cycle < (codon_window * 10)):

        keep_looping = False

        # iterate by codon, but map back to sequence-based indices
        for i in range(0, len(mutable_seq) // 3, (codon_window - overlap) * 3):
            window = slice(
                i * 3,
                (i + codon_window) * 3
                if (i + codon_window) * 3 < len(mutable_seq)
                else len(mutable_seq),
            )

            # make each mutable codon immutable so it can be hashed later
            codons = [
                mutable_seq[window][i : i + 3].toseq()
                for i in range(0, len(mutable_seq[window]), 3)
            ]

            # TODO: do we only want to check the first three? the first n? or all?
            # in this script the window is 9 bp, so it's only ever three codons

            # check if all codons in the window are identical
            if len(set(codons)) == 1:
                logger.detail("All codons in window are identical: {0}".format(codons))
                keep_looping = _mutate_and_keep_looping(mutable_seq, window, (i * 3))

            # check if the segment is found in the full sequence
            non_overlapping_matches = re.findall(
                str(mutable_seq[window]), str(mutable_seq)
            )
            if len(non_overlapping_matches) > 1 and len(mutable_seq[window]) > 3:
                logger.debug("Current window is repeated in the sequence")
                keep_looping = _mutate_and_keep_looping(mutable_seq, window, (i * 3))

        current_cycle += 1

    return mutable_seq.toseq()


##########################################################
#
# 	ACTUAL SCRIPT STARTS HERE
#
##########################################################

logger.info("===== SCRIPT START =====")

# read the input sequence file and parse lines
records = list(SeqIO.parse(args.input, "fasta", IUPAC.protein))

logger.info("Total number of sequences: {0}".format(len(records)))
[logger.detail(record) for record in records]

# generate host profile
codon_use_table = codon_use.host_codon_usage(args.host)

# process through all supplied sequences
for count, record in enumerate(records):
    logger.info("===== PROCESSING SEQUENCE {0} ===== {1}".format(count + 1, record.id))

    # check input seq style
    logger.info("INPUT IS AA SEQUENCE")
    dna = record.seq.back_translate()

    logger.detail("===== DUMPING SEQUENCE =====")
    logger.detail(str(dna))

    # TODO: compare input and host profiles
    # count codons and current profile
    count_table = codon_use.count_codons(dna)
    input_profile = codon_use.calc_profile(count_table)

    # run optimization
    difference, current_cycle, relax = 1.0, 1, 1.0

    # keep running while there are cycles AND difference between current and host is less than the % relax allowed
    while ((current_cycle <= args.cycles) or (args.cycles == 0)) and (
        difference >= (relax - 1)
    ):
        logger.info(
            "~~~~~~~~~~ Current cycle: {0}/{1} ~~~~~~~~~~".format(
                current_cycle, args.cycles
            )
        )

        # TODO: measure the deviation from the host profile
        # probably good to record the sequence that is closest
        # but at the end, after all the clean up is done.
        dna = resample_codons(dna, codon_use_table)

        # identify and remove undesirable features
        for _, gc_content in GC_content.items():
            # check various GC content requirements
            dna = gc_scan(dna, **gc_content)
        dna = remove_restriction_sites(dna, RestrictionEnz)
        dna = remove_start_sites(dna, RibosomeBindingSites, "Standard")
        dna = remove_repeating_sequences(dna, 9)
        dna = remove_local_homopolymers(dna)

        # count codons and current profile
        count_table = codon_use.count_codons(dna)
        input_profile = codon_use.calc_profile(count_table)

        # compare input and host profiles
        # this should probably use the CAI metric
        """
        # determine how much to relax harmonization
        relax = 1 + (args.max_relax * ((current_cycle - 1) / args.cycles))
        logger.info("Relax coeff: {0}".format(relax))

        mutation_table, difference = compare_profiles(
            input_profile, host_profile, relax
        )
        """

        # tick cycle
        current_cycle += 1

    # hit the max number of cycles?
    if current_cycle > args.cycles:
        logger.info("You hit the max number of cycles: {0}".format(args.cycles))

    # check GC content
    logger.info("===== GC CONTENT =====")
    gc_percent = round(GC(dna) / 100, 3)
    if gc_percent < 0.3 or gc_percent > 0.65:
        logger.warning("Overall GC content is {0}!".format(gc_percent))

    # display result name and sequence
    logger.output("===== SEQUENCE NAME =====")
    logger.output("{0}".format(record.id))

    # display final codon-use difference between host and current sequence (0.00 is ideal)
    logger.output("({0})".format(round(difference, 2)))
    logger.output("===== DUMPING SEQUENCE =====")
    logger.output(str(dna))
