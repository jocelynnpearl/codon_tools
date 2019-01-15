#! /usr/bin/env python
import argparse
import os.path
import random
import re
import sys

import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import CodonUsage, GC, seq3

from codon_tools.util import Seq, MutableSeq, logging, log_levels
from codon_tools.data import GC_content, RestrictionEnz, RibosomeBindingSites

##########################################################
#
# 	Usage:
# 		$python codon_tools.py --input INPUT_LIST.list
#
# 	Requires:
# 		Biopython, Python 3.4
#
# 	You can install Biopython by:
# 		$pip install biopython
#
# 	INPUT_LIST.list:
# 		>SEQ_1
# 		ACDEFGHIKLMNPQRSTVWY
# 		>SEQ_2
# 		ACDEFGHIKLMNPQRSTVWY
#
# 	What it does:
# 		1) reverse translates your input AA into an arbitrary DNA sequence OR translates your input DNA into AA.
# 		2) calculates the host's per-AA codon usage profile -- if the codon is used less than 10% (variable) of the time it is considered 0 instead.
# 		3) compares the current DNA sequence to the host's profile and determine which codons are overused/underused.
# 		4) stochastically mutates overused codons to underused codons.
# 		5) runs a series of other checks and optimizations:
# 			a) checks sequence for GC content and mutates codons to fall into a reasonable GC % (monte carlo)
# 			b) checks for unwanted restriction sites, stochastically mutate the codons if found.
# 			c) looks for all ATG/GTG/TTG sequences, then checks 18bp upstream of it for GA -rich regions. If found, mutate the codons.
# 			d) looks for 3-consecutive identical codons and 9-mer repeat chunks (in frame to speed stuff up) and mutate them away.
# 			e) checks for "local homopolymers", if there are areas with more than 4 (variable) consecutive identical bps, stochastically mutate the codons.
# 		6) repeat from step 3.
# 		7) cycles through until (variable) cycles are hit OR there the per-AA codon profile of current DNA and host profile matches, given a maximum deviation (and passes the checks).
#
# 	To do:
# 		1) remove RNA structure from sequence
#
##########################################################

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
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        # add up number of times each AA appears
        tot_usage = sum([codons_count[codon] for codon in synonymous_codons])

        # calculate the distribution of synonymous codon use
        if tot_usage == 0:
            continue
        codons_freq.update(
            {codon: codons_count[codon] / tot_usage for codon in synonymous_codons}
        )
    return codons_freq


# returns a calc_profile for a given host table
def process_host_table():
    logger.info("===== PROCESSING HOST TABLE: {0} =====".format(args.host))
    table = CodonUsage.CodonsDict.copy()

    dirname = os.path.dirname(__file__)
    filename = os.path.join(
        dirname, "template_files/codon_tables/{0}.txt".format(args.host)
    )
    with open(filename, "r") as inputfile:
        for line in inputfile:
            codon, _, count = line.split()

            if codon == "\0":
                continue
            # codons are stored with RNA alphabet
            table[str(Seq(codon).back_transcribe())] = int(count)

    calculated_table = calc_profile(table)

    # avoid iteration if the log level is disabled
    if logger.isEnabledFor(logging.DETAIL):
        logger.detail("pre-threshold host table:")
        for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
            logger.detail("== {0} ==".format(AA))
            for codon in synonymous_codons:
                logger.detail(
                    "{0}: {1}, {2}".format(codon, table[codon], calculated_table[codon])
                )

    logger.info("HOST THRESHOLD SET TO: {0}".format(args.host_threshold))
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        for codon in synonymous_codons:
            if calculated_table[codon] < args.host_threshold:
                table[codon] = 0

    # recalculate profile after threshold applied
    calculated_table = calc_profile(table)

    if logger.isEnabledFor(logging.DETAIL):
        logger.detail("post-threshold host table:")
        for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
            logger.detail("== {0} ==".format(AA))
            for codon in synonymous_codons:
                logger.detail(
                    "{0}: {1}, {2}".format(codon, table[codon], calculated_table[codon])
                )
    return calculated_table


def host_codon_usage(host=args.host):
    logger.info("===== PROCESSING HOST TABLE: {0} =====".format(host))
    host_profile = CodonUsage.CodonsDict.copy()

    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, "template_files/codon_tables/{0}.txt".format(host))
    with open(filename, "r") as inputfile:
        for line in inputfile:
            codon, frequency_per_1000, _ = line.split()

            if codon == "\0":
                continue
            # codons are stored with RNA alphabet
            host_profile[str(Seq(codon).back_transcribe())] = (
                float(frequency_per_1000) * 10
            )

    codon_use_by_aa = {}
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        codon_use_by_aa[AA] = [
            synonymous_codons,
            [host_profile[codon] for codon in synonymous_codons],
        ]
    return codon_use_by_aa


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
def remove_local_homopolymers(input):
    logger.info("===== REMOVE LOCAL HOMOPOLYMERS =====")
    check = 0
    while check != 1:
        # look at each 6-mer
        for x in range(0, (len(input) - 1)):
            mer = list(input[x][1] + input[x + 1][1])

            for base in ["T", "C", "A", "G"]:
                counter = 0
                found = 0
                max_found = 0
                while counter < len(mer):
                    if base == mer[counter]:
                        found += 1
                        if found > max_found:
                            max_found = found
                    else:
                        found = 0
                    counter += 1

                if max_found > args.local_homopolymer_threshold:
                    logger.detail("position: {0}: {1}".format(x * 3, mer))
                    logger.detail("{0}, count={1}".format(base, max_found))
                    apply_mutation(input, x)
                    apply_mutation(input, x + 1)
                    check = 0
                else:
                    check = 1
    return input


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
            # ensure the whole codon we mutate is being recognized
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

    dna_sequence = Seq("ATGGCCATTGTAATGGGCCGCTG", IUPAC.unambiguous_dna)
    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(dna_sequence))
    ]

    if len(start_codon_positions) == 0:
        logger.info("No start codons found in sequence")
        return dna_sequence
    else:
        logger.info(
            "Found {0} start codons. Checking for upstream RBSs...".format(
                len(start_codon_positions)
            )
        )

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    rbs_positions = [pos - 18 for pos in start_codon_positions if pos >= 18]
    mutable_seq = dna_sequence.tomutable()

    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + 18 - 3
        rbs_query_seq = mutable_seq[rbs_start:rbs_stop]

        if logger.isEnabledFor(logging.DETAIL):
            seq_frag_with_start = mutable_seq[rbs_stop : rbs_stop + 6]
            logger.detail(
                "checking sequence: {0}.{1}".format(rbs_query_seq, seq_frag_with_start)
            )

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            logger.detail("checking start site: {0}, {1}".format(rbs, site))
            search = rbs_query_seq.find(site)
            # test search
            count = 0
            while search != -1:
                position = (search + rbs_start + 1) // 3

                # mutate residues if site is found

                apply_mutation(mutable_seq, position)
                apply_mutation(mutable_seq, position + 1)

                # reset sequence and search again
                seq_frag = mutable_seq[rbs_start : rbs_stop + 3]
                search = seq_frag.find(site)

                # have a counter so it doesn't get infinite looped.
                count += 1
                if count >= 10:
                    break
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
    codon_window = int(window_size / 3)
    overlap = int(codon_window / 2)
    mutable_seq = dna_sequence.tomutable()

    # for seg_no, segment in enumerate(slidingWindow):
    for i in range(0, len(mutable_seq), codon_window - overlap):
        window = slice(
            i,
            i + codon_window
            if i + codon_window < len(mutable_seq)
            else len(mutable_seq),
        )
        logger.debug("Current segment: {0}".format(mutable_seq[window]))

        count = 0
        gc_percent = GC(mutable_seq[window]) / 100

        # check gc_percent of current segment
        while (gc_percent < low or gc_percent > high) and count < codon_window * 2:
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice(i + position, i + position + 3)

            init_codon = mutable_seq[codon_idx]
            new_codon = mutate_codon(init_codon, codon_use_table)

            if (GC(new_codon) < GC(init_codon) and gc_percent > high) or (
                GC(new_codon) > GC(init_codon) and gc_percent < low
            ):
                mutable_seq[codon_idx] = new_codon
                logger.debug("Mutating position: {0}".format(position))
                gc_percent = GC(mutable_seq[window]) / 100

            # prevent infinite loop
            count += 1
    return mutable_seq.toseq()


# check for repeat segments
def repeat_scan(input, frag_size):
    logger.info("===== REPEAT FRAGMENT SCAN FOR SIZE: {0} bps =====".format(frag_size))
    random.seed()

    # determine window and overlap size
    window = int((frag_size / 3))
    overlap = window - 1
    # dummy used so the fragment won't find itself
    dummy = "".zfill(window * 3)
    splices = [
        input[i : i + window] for i in range(0, len(input) - 2, window - overlap)
    ]

    loop_this = 1
    # this to make sure poly-TRP or poly-MET won't just infinite loop
    loop_break = 0

    while (loop_this != 0) and (loop_break < (window * 10)):
        # loop this if any mutation is made, until it passes through both checks without mutations
        loop_this = 0
        # loop through each segment
        for segment in splices:
            # check if it's 3 consecutive same codons
            if segment[0][1] == segment[1][1] and segment[1][1] == segment[2][1]:
                logger.detail("first 3 codons are identical: {0}".format(segment))
                num_mut = random.randint(1, 2)
                logger.debug("roll dice, mutate {0} codons".format(num_mut))
                for f in range(0, num_mut):
                    position = random.randint(0, num_mut - 1)
                    apply_mutation(segment, position)
                loop_this += 1
                # check if the segment is found in the full sequence
            seq = "".join(codon for innerlist in input for codon in innerlist[1])
            target = "".join(codon for innerlist in segment for codon in innerlist[1])
            seq = seq.replace(target, dummy, 1)
            if seq.find(target) != -1:
                logger.detail("Repeat fragment found with segment: {0}".format(segment))
                num_mut = random.randint(1, (frag_size / 3) - 1)
                logger.debug("roll dice, mutate {0} codons".format(num_mut))
                for f in range(0, num_mut):
                    position = random.randint(0, (frag_size / 3) - 1)
                    apply_mutation(segment, position)
                loop_this += 1
        loop_break += 1
    return input


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
codon_use_table = host_codon_usage(args.host)

# process through all supplied sequences
for count, record in enumerate(records):
    logger.info("===== PROCESSING SEQUENCE {0} ===== {1}".format(count + 1, record.id))

    # check input seq style
    logger.info("INPUT IS AA SEQUENCE")
    input_dna = record.seq.back_translate()

    logger.detail("===== DUMPING SEQUENCE =====")
    logger.detail("".join([seq[1] for seq in input]))

    # TODO: compare input and host profiles
    """
    # count codons and current profile
    count_table = count_codons(input_dna)
    input_profile = calc_profile(count_table)
    """

    # run optimization
    difference = 1
    cycles_current = 0
    relax = 1

    # keep running while there are cycles AND difference between current and host is less than the % relax allowed
    while ((cycles_current < args.cycles) or (args.cycles == 0)) and (
        difference >= (relax - 1)
    ):
        logger.info(
            "~~~~~~~~~~ Current cycle: {0}/{1} ~~~~~~~~~~".format(
                cycles_current + 1, args.cycles
            )
        )
        # determine how much to relax harmonization
        relax = 1 + (args.max_relax * ((cycles_current) / args.cycles))
        logger.info("Relax coeff: {0}".format(relax))

        # TODO: measure the deviation from the host profile
        input_dna = resample_codons(input_dna, codon_use_table)

        # check GC content in window
        for _, gc_content in GC_content.items():
            gc_scan(input_dna, **gc_content)

        # check for unwanted restriction sites
        remove_restriction_sites(input_dna, RestrictionEnz)

        # check for alternative start sites
        remove_start_sites(input_dna, RibosomeBindingSites, "Standard")

        # check for repeat fragments
        repeat_scan(input_dna, 9)
        # check for local homopolymers
        remove_local_homopolymers(input_dna)

        # count codons and current profile
        count_table = count_codons(input_dna)
        input_profile = calc_profile(count_table)

        # compare input and host profiles
        # this should probably use the CAI metric
        """
        mutation_table, difference = compare_profiles(
            input_profile, host_profile, relax
        )
        """

        # tick cycle
        cycles_current += 1

        # hit the max number of cycles?
    if cycles_current == args.cycles:
        logger.info("You hit the max number of cycles: {0}".format(args.cycles))

    # check GC content
    logger.info("===== GC CONTENT =====")
    gc_percent = round(GC(input_dna) / 100, 3)
    if gc_percent < 0.3 or gc_percent > 0.65:
        logger.warning("Overall GC content is {0}!".format(gc_percent))

    # display result name and sequence
    logger.output("===== SEQUENCE NAME =====")
    logger.output("{0}".format(record.id))

    # display final codon-use difference between host and current sequence (0.00 is ideal)
    logger.output("({0})".format(round(difference, 2)))
    logger.output("===== DUMPING SEQUENCE =====")
    logger.output("".join([seq[1] for seq in input]))
