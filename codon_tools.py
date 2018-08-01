#! /usr/bin/env python
import sys
import argparse
import Bio
import random
import re
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import CodonUsage
from Bio.SeqUtils import seq3

##########################################################
#
#	Usage:
#		$python codon_tools.py --input INPUT_LIST.list
#	
#	Requires:
#		Biopython, Python 3.4
#	
#	You can install Biopython by:
#		$pip install biopython 
#	
#	INPUT_LIST.list:
# 		>SEQ_1
#		ACDEFGHIKLMNPQRSTVWY
#		>SEQ_2
#		ACDEFGHIKLMNPQRSTVWY
#	
#	What it does:
#		1) reverse translates your input AA into an arbitrary DNA sequence OR translates your input DNA into AA.
#		2) calculates the host's per-AA codon usage profile -- if the codon is used less than 10% (variable) of the time it is considered 0 instead.
#		3) compares the current DNA sequence to the host's profile and determine which codons are overused/underused.
#		4) stochastically mutates overused codons to underused codons.
#		5) runs a series of other checks and optimizations:
#			a) checks sequence for GC content and mutates codons to fall into a reasonable GC % (monte carlo)
#			b) checks for unwanted restriction sites, stochastically mutate the codons if found.
#			c) looks for all ATG/GTG/TTG sequences, then checks 18bp upstream of it for GA -rich regions. If found, mutate the codons.
#			d) looks for 3-consecutive identical codons and 9-mer repeat chunks (in frame to speed stuff up) and mutate them away.
#			e) checks for "local homopolymers", if there are areas with more than 4 (variable) consecutive identical bps, stochastically mutate the codons.
#		6) repeat from step 3.
#		7) cycles through until (variable) cycles are hit OR there the per-AA codon profile of current DNA and host profile matches, given a maximum deviation (and passes the checks).
#
#	To do:
#		1) remove RNA structure from sequence
#
##########################################################

#options
parser = argparse.ArgumentParser(
	description='Optimize your AA or DNA sequence to harmonize with a host\'s condon usage.',
	epilog='2017-12-04, v0.47 (contact yhsia@uw.edu if stuff does not work)' )

parser.add_argument('--input', type=str, required=True, help='input file with sequence' )
parser.add_argument('--type', type=str, default='AA', choices=['AA', 'DNA'], help='is your input AA or DNA?' )
parser.add_argument('--host', type=str, default='413997', help='host table code: http://www.kazusa.or.jp/codon/, default is "Escherichia coli B"' )
parser.add_argument('--host_threshold', type=float, default='0.10', help='lowest codon fraction per AA in the host that is allowed' )
#parser.add_argument('--n_terminal_rares', type=bool, default='true', help='optimize the first 11 AAs 
parser.add_argument('--verbose', type=int, default=0, choices=[0, 1, 2, 3], help='verbose output level (0=only result, 1=standard output, 2=extra output 3=debugging)' )
parser.add_argument('--local_homopolymer_threshold', type=int, default='4', help='number of consecutive NT repeats allowed' )
parser.add_argument('--cycles', type=int, default=1000, help='max number of cycles to run optimization, 0=unlimited' )
parser.add_argument('--max_relax', type=float, default='0.1', help='maximum % deviation from host profile' )

args = parser.parse_args()

#reverse translate AA seq to DNA seq
def reverse_translate( input ):
	if args.verbose >= 1: print( "===== REVERSE TRANSLATING AA SEQUENCE =====" )
	aa_list = list( input )
	trans=[]
	x=0
	while x < len( input ):
		#print( aa_list[ x ] )
		codon = CodonUsage.SynonymousCodons[ str(seq3( aa_list[ x ] )).upper() ][0]
		trans.append( [ aa_list[ x ], codon, x+1 ] )
		x += 1
	if args.verbose >= 3: print( trans )
	return trans

#translate and format DNA seq
def translate_input( input ):
	if args.verbose >= 1: print( "===== TRANSLATING AA SEQUENCE =====" )
	input_triplet = store_triplets( input )
	#print( input_triplet )
	trans=[]
	x=0
	while x < ( len( input )/3 ):
		#print( str(Seq( input_triplet[x], IUPAC.unambiguous_dna ).translate) )
		if len( input_triplet[x] ) == 3:
			trans.append( [ str(Seq( input_triplet[x], IUPAC.unambiguous_dna ).translate()), input_triplet[x], x+1 ] )
		x += 1
	if args.verbose >= 3: print( trans )
	return trans

#returns list, splits an input dna sequence into a list of triplet codons
def store_triplets( input ):
	list=[]
	#make sure input DNA is divisible by 3
	assert len( input ) % 3 == 0, "Input DNA is not divisible by 3!"
		
	for start in range(0, len( input ), 3):
		list.append( input[start:start+3] )
	return list

#returns dictionary, counts the number of times each codon is used
def count_codons( input ):
	if args.verbose >= 1: print( "===== COUNTING CODONS =====" )
	codon_list=[]; x=0
	while x < len ( input ):
		codon_list.append( input[x][1] )
		x += 1
		
	table={}
	for base1 in ["T","C","A","G" ]:
		for base2 in ["T","C","A","G" ]:
			for base3 in ["T","C","A","G" ]:
					codon_count=codon_list.count( '{0}{1}{2}'.format( base1, base2, base3 ) )
					#print( codon_count )
					#codon_fraction=codon_count/len( input )
					#print( codon_fraction )
					#table.append( {'{0}{1}{2}'.format( base1, base2, base3 ): codon_count} )
					table['{0}{1}{2}'.format( base1, base2, base3 )] = [codon_count,0]
					#table.append( [codon_count, codon_fraction] )
					#print( table[-1] )
	return table

#returns dictionary, calculates the % usage of each AA's codons
def calc_profile( input ):
	if args.verbose >= 1: print( "===== CALCULATING PROFILE =====" )
	#loop through all amino acids
	for AA in CodonUsage.SynonymousCodons:
		#print( '== {0} =='.format( AA ) )
		#print( len( CodonUsage.SynonymousCodons[ AA ] ) )
		
		#add up number of times each AA(total codons) is used
		tot_usage=0
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			tot_usage += input[ syn_codon ][0]
		#print( tot_usage )
		
		#calculate how the fraction of times each codon is used per AA
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			usage=0
			if tot_usage == 0:
				usage = 0
			else:
				usage = input[ syn_codon ][0] / tot_usage
			#input.setdefault( syn_codon, [] )
			input[ syn_codon ][1] = usage
	
	#print( input )		
	return input
	
#returns a calc_profile for a given host table
def process_host_table():
	if args.verbose >= 1: print( "===== PROCESSING HOST TABLE: {0} =====".format( args.host ) )
	table={}
	dir = os.path.dirname(__file__)
	filename = os.path.join(dir, 'template_files/codon_tables/{0}.txt'.format( args.host ) )
	with open( filename, "r" ) as inputfile:
		for line in inputfile:
			tok = line.split()
			if tok[0] == '\0':
				continue
			codon_rna = Seq( tok[0] )
			codon_dna = codon_rna.back_transcribe()
			table[ '{0}{1}{2}'.format( codon_dna[0], codon_dna[1], codon_dna[2] ) ] = [ int( tok[2] ),0 ]
	#print( table )
	
	calculated_table = calc_profile( table )
	
	if args.verbose >= 2:
		print( "pre-threshold host table:" )
		for AA in CodonUsage.SynonymousCodons:
			print( '== {0} =='.format( AA ) )
			for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
				print( '{0}: {1}, {2}'.format( syn_codon, calculated_table[ syn_codon ][0], calculated_table[ syn_codon ][1] ) )
		#print(calculated_table)
	
	if args.verbose >= 1: print( "HOST THRESHOLD SET TO: {0}".format( args.host_threshold ) )
	for AA in CodonUsage.SynonymousCodons:
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			if calculated_table[ syn_codon ][1] < args.host_threshold:
				calculated_table[ syn_codon ][0] = 0
	
	#recalculate profile after threshold applied
	calculated_table = calc_profile( calculated_table )

	if args.verbose >= 2:
		print( "post-threshold host table:" )
		for AA in CodonUsage.SynonymousCodons:
			print( '== {0} =='.format( AA ) )
			for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
				print( '{0}: {1}, {2}'.format( syn_codon, calculated_table[ syn_codon ][0], calculated_table[ syn_codon ][1] ) )
		#print(calculated_table)
		
	return calculated_table	

#returns dictionary after a key gets removed (this is for reference issues).
def removekey(d, key):
    r = dict(d)
    del r[key]
    return r
    
#returns dictionary of comparison between two profiles
def compare_profiles( input, host, relax ):
	if args.verbose >= 1: print( "===== COMPARING PROFILES =====" )
	table={}
	#loop AAs
	for AA in CodonUsage.SynonymousCodons:
		if args.verbose >= 2: print( AA )
		temp_table={}
		#calculate total usage of codon in input
		tot_usage=0
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			tot_usage += input[ syn_codon ][0]
		
		#calculate ideal usage of codon in host
		tot_ideal=0
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			ideal_usage_abs = int( round( host[ syn_codon ][1] * tot_usage, 0 ) )
			ideal_usage = int( round( host[ syn_codon ][1] * relax * tot_usage, 0 ) )
			if args.verbose >= 2: print( "{0}: {1}".format(syn_codon,ideal_usage) )
			tot_ideal += ideal_usage
			temp_table[ syn_codon ] = { 'input_count': input[ syn_codon ][0], 'input_perc': input[ syn_codon ][1], 'ideal_usage_abs': ideal_usage_abs,'ideal_usage': ideal_usage, 'host_perc': float( host[ syn_codon ][1] ) }
		
		#if ideal number is too high, subtract from lowest host codon (that is not 0 already )
		while tot_ideal > tot_usage:
			lowest_table = temp_table
			#print( lowest_table )
			lowest_codon = min( temp_table, key=lambda k:float(temp_table[k]['host_perc']) )
			while temp_table[ lowest_codon ]['ideal_usage'] == 0 :
				#print( lowest_table )
				lowest_table = removekey( lowest_table, lowest_codon )
				#print( lowest_codon )
				lowest_codon = min( lowest_table, key=lambda k:float(lowest_table[k]['host_perc']) )
				#print( lowest_codon )
			temp_table[ lowest_codon ]['ideal_usage'] = temp_table[ lowest_codon ]['ideal_usage'] - 1
			#print( temp_table[ lowest_codon ]['ideal_usage'] )
			tot_ideal -= 1
		
		#if ideal number is too low, add to highest host codon
		while tot_ideal < tot_usage:
			highest_codon = max( temp_table, key=lambda k:float(temp_table[k]['host_perc']) )
			temp_table[ highest_codon ]['ideal_usage'] = temp_table[ highest_codon ]['ideal_usage'] + 1		
			tot_ideal += 1
							
		#populate return table
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			table[ syn_codon ] = { 'input_count': temp_table[ syn_codon ]['input_count'], 'ideal_usage_abs': temp_table[ syn_codon ]['ideal_usage_abs'], 'difference': temp_table[ syn_codon ]['input_count'] - temp_table[ syn_codon ]['ideal_usage'], 'difference_abs': temp_table[ syn_codon ]['input_count'] - temp_table[ syn_codon ]['ideal_usage_abs']}

	#calculate difference value
	resi_total=0
	diff_total=0
	for AA in CodonUsage.SynonymousCodons:
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			resi_total += table[ syn_codon ]['ideal_usage_abs']
			diff_total += abs( table[ syn_codon ]['difference_abs'] )
	if args.verbose >= 1: print( 'CURRENT DIFFERENCE TOTAL: {0} of {1}'.format( int( diff_total/2 ), resi_total ) )
	if args.verbose >= 1: print( 'CURRENT DIFFERENCE %: {0}'.format( diff_total/resi_total/2 ) )		
	diff = diff_total/resi_total/2
	
	return table, diff

#returns mutated sequence after given a difference profile
def optimize_sequence( input, mutation_profile ):
	if args.verbose >= 1: print( "===== OPTIMIZING SEQENCE =====" )
	random.seed()
	for AA in CodonUsage.SynonymousCodons:
		#figure out the index of codons in input sequence
		mutation_table={}
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			mutation_table[ syn_codon ] = { 'difference': mutation_profile[ syn_codon ][ 'difference' ], 'pos': None }
			#print( syn_codon )
			pos_list=[]
			for pos, item in enumerate( input ):
				if item[1] == syn_codon:
					pos_list.append( pos )
					#print( pos_list )
			mutation_table[ syn_codon ].update( {'pos': pos_list } )
		
		#print( AA )
		#print( mutation_table )
		#check if this AA even needs to be optimized
		tot_diff=0
		for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
			tot_diff += abs( mutation_table[ syn_codon ]['difference'] )
				
		while tot_diff >= 1 :
			#pick random codon1 and codon2
			
			##BLOCK HERE TO SKIP LOCKED CODONS			
			
			codon1=random.choice( list( mutation_table.keys() ) )
			while ( mutation_table[ codon1 ]['difference'] == 0 or
					mutation_table[ codon1 ]['difference'] < 0 ):
				codon1=random.choice( list( mutation_table.keys() ) )
				#print( codon1 )
			codon2=codon1
			#print( codon1 )
			
			while ( codon1 == codon2 or
					mutation_table[ codon2 ]['difference'] == 0 or
					mutation_table[ codon2 ]['difference'] > 0 ): #make sure codon2 is always receiving
				codon2=random.choice( list( mutation_table.keys() ) )
			#print( codon2 )

			#change codons positions
			#print( AA)
			#print( list( mutation_table[ codon1 ]['pos'] ) )
			codon_pos = random.choice( list( mutation_table[ codon1 ]['pos'] ) )
			#print( codon_pos )
			codon_pos_index = mutation_table[ codon1 ]['pos'].index( codon_pos )
			#delete from codon1
			#print( mutation_table[ codon1 ][ 'pos' ] )
			del mutation_table[ codon1 ][ 'pos' ][ codon_pos_index ]
			mutation_table[ codon1 ]['difference'] -= 1
			#print( mutation_table[ codon1 ][ 'pos' ] )
			#add to codon2
			mutation_table[ codon2 ]['pos'].append( codon_pos )
			mutation_table[ codon2 ]['difference'] += 1
			#print( input[ codon_pos ][1] )
			input[ codon_pos ][1] = codon2
			#print( input[ codon_pos ][1] )
			#print( mutation_table[ codon2 ][ 'pos' ] )

			#update difference
			tot_diff=0
			for syn_codon in CodonUsage.SynonymousCodons[ AA ]:
				tot_diff += abs( mutation_table[ syn_codon ]['difference'] )			
			
		#print( mutation_table ) #this should result in difference=0 for all
	return input

#check for local homopolymers
def remove_local_homopolymers( input ):
	if args.verbose >= 1: print ( "===== REMOVE LOCAL HOMOPOLYMERS =====" )
	check=0
	while ( check != 1 ):
		#look at each 6-mer
		for x in range(0, ( len( input ) - 1 ) ):
			mer = list( input[x][1] + input[x+1][1] )
	
			for base in ["T","C","A","G" ]:
				counter=0
				found=0
				max_found=0
				while counter < len( mer ):
					if base == mer[ counter ]:
						found += 1
						if found > max_found:
							max_found=found
					else:
						found=0
					counter += 1
				
				if max_found > args.local_homopolymer_threshold:
					if args.verbose >= 2:
						print( "position: {0}: {1}".format( x*3, mer ) )
						print( "{0}, count={1}".format( base, max_found ) )
					apply_mutation( input, x )
					apply_mutation( input, x+1 )
					check = 0
				else:
					check = 1
	return input

#parse unwanted sites file
def parse_unwanted_site_file( path_to_file, type ):	
	unwanted_sites = []
	with open( path_to_file, 'r' ) as sites_file:
		for count, line in enumerate( sites_file ):
			if count % 2 == 0:
				site_name = line
				site_name = site_name.replace('\n', '')
			else:
				seq = line
				seq = seq.replace('\n', '')
				unwanted_sites.append( [ site_name, seq ] )
	if args.verbose >= 1: print( 'Total number of unwanted {1} sites: {0}'.format( len(unwanted_sites), type ) )
	if args.verbose >= 2: print( unwanted_sites )
	return unwanted_sites

#check for unwanted restriction sites
def remove_restriction_sites( input, restrict_sites ):
	if args.verbose >= 1: print ( "===== REMOVE RESTRICTION SITES =====" )

	#check each unwanted restriction site
	for site_pair in restrict_sites:
		if args.verbose >= 1: print ( 'checking restriction site: {0}, {1}'.format( site_pair[0], site_pair[1] ) )

		#create sequence in a single string
		seq = ''.join( codon for innerlist in input for codon in innerlist[1] )
		#search for the restriction site in the sequence
		search = seq.find( site_pair[1] )
		#test search
		count=0
		while search != -1:
			position = int( search / 3 )
			#mutate residues if site is found
			apply_mutation( input, position )
			apply_mutation( input, position+1 )
			#reset sequence and search again
			seq = ''.join( codon for innerlist in input for codon in innerlist[1] )
			search = seq.find( site_pair[1] )
			
			#have a counter so it doesn't get infinite looped.
			count += 1
			if count >= 10: break
	return input
	
#check for alternative start sites
def remove_start_sites( input, start_sites, start_codon ):
	if args.verbose >= 1: print ( "===== REMOVE START SITES: {0} =====".format( start_codon ) )

	#create sequence in a single string
	seq = ''.join( codon for innerlist in input for codon in innerlist[1] )
	#find all start codon sites (xTG)
	find_start = [ m.start() for m in re.finditer( start_codon, seq ) ]
	if len( find_start ) == 0:
		if args.verbose >= 1: print( "No start codon found in sequence" )
		return input
	else:
		if args.verbose >= 1: print( "Start codons found: {0}".format( len( find_start ) ) )
	
	#check each start site for RBS (18 base pairs upstream of each ATG, ignore 3 bp closest to ATG )
	for start_pos in find_start:
		rbs_start = start_pos - 18
		seq_frag = seq[ rbs_start : start_pos - 3 ]
		
		if args.verbose >= 2:
			seq_frag_with_start = seq[ start_pos - 3 : start_pos + 3 ]
			print( 'checking sequence: {0}.{1}'.format( seq_frag, seq_frag_with_start ) )
		
		#check each unwanted RBS in each potential fragment
		for site_pair in start_sites:
			if args.verbose >= 2: print ( 'checking start site: {0}, {1}'.format( site_pair[0], site_pair[1] ) )
			search = seq_frag.find( site_pair[1] )
			#test search
			count=0
			while search != -1:
				position = int( ( search + rbs_start + 1 ) / 3 )
				#mutate residues if site is found
				apply_mutation( input, position )
				apply_mutation( input, position+1 )
				#reset sequence and search again
				seq = ''.join( codon for innerlist in input for codon in innerlist[1] )
				seq_frag = seq[ rbs_start : start_pos ]
				search = seq_frag.find( site_pair[1] )
				
				#have a counter so it doesn't get infinite looped.
				count += 1
				if count >= 10: break
	return input

#mutate single codon
def mutate_codon( codon_in ):
	random.seed()
	AA=str(seq3( str( Seq( codon_in[1], IUPAC.unambiguous_dna ).translate() ) )).upper() 
	num_codons=len( CodonUsage.SynonymousCodons[ AA ] )
		
	#pick new codon
	codon_out=random.choice( CodonUsage.SynonymousCodons[ AA ] )
	while ( codon_in[1] == codon_out ) and ( num_codons != 1 ):
		codon_out=random.choice( CodonUsage.SynonymousCodons[ AA ] )
	
	if args.verbose >= 2:
		print( 'mutating [{0}] codon at position {3} from {1} to {2}'.format( AA, codon_in[1], codon_out, codon_in[2] ) )
	return codon_out
	
#apply mutation to position x
def apply_mutation( input, position, old_codon=[] ):
	old_codon.clear()
	old_codon.append( input[ position ][1] )
	#print( old_codon )
	new_codon = mutate_codon( input[ position ] )
	input[ position ][1] = new_codon
	#print( input[ position ][1] )
	return input, position, old_codon
	
#sync mutation from template to mirror
#def sync_mutations( template, mirror ):
#	print( "sync_mutations" )
#	for template_codon in template:
#		abs_posi = template_codon[2]
#		if template_codon[1] == mirror[ abs_posi - 1 ][1]:
#			print( "template: {0}".format( template_codon ) )
#			print( "mirror: {0}".format( mirror[ abs_posi - 1 ] ) )
#		#print( template_codon )
#	return template, mirror

#dumps name of sequence
def dump_name( input ):
	if args.verbose >= 1: print( "===== SEQUENCE NAME =====" )
	print( "{0}".format( sequences[ count ][0] ), end=' ' )
	if args.verbose >= 1: print( '' )
	
#dumps string from triplet sequence
def dump_sequence( input ):
	if args.verbose >= 1: print( "===== DUMPING SEQUENCE =====" )
	for x in range(0, len( input ) ):
		print( input[x][1], end='')
	print()

#check GC content
def gc_scan( input, abs_window, low, high ):
	if args.verbose >= 1: print( "===== GC CONTENT SCAN IN WINDOW: {0} bps, threshold: {1} < x < {2}=====".format( abs_window, low, high ) )
	random.seed()

	#splice sequence into overlapping chunks
	window = int( abs_window / 3 )
	overlap = int( window / 2 )
	splices = [ input[ i: i + window ] for i in range( 0, len( input ),  window - overlap ) ] #this is somehow a reference to "input"?
	for segment in splices:
		gc_percent = gc_content( segment )
		if args.verbose >= 3: print( "Current segment: {0}".format( segment ) )
		count=0
		old_codon=[]
		#check gc_percent of current segment
		while gc_percent < low or gc_percent > high:
			position=random.randint( 0, len( segment ) - 1 )
			apply_mutation( segment, position, old_codon )
			gc_percent_new = gc_content( segment )
			if gc_percent_new >= gc_percent:
				segment[ position ][1] = old_codon[0]
				if args.verbose >= 3: print( "reverting position: {0}".format( segment[ position ] ) )
			else:
				gc_percent = gc_percent_new
			#have a counter so it doesn't get infinite looped. (= 2x segment length)
			count += 1
			if count >= ( len( segment ) * 2 ): break
		#sync_mutations( segment, input ) #don't need this since "splices" references "input" somehow
	return input

#count GC content
def gc_content( input ):
	seq = ''.join( codon for innerlist in input for codon in innerlist[1] )
	g_count = seq.count('G')
	c_count = seq.count('C')
	if args.verbose >= 1: print( ( g_count + c_count ) / len( seq ) )
	return( ( g_count + c_count ) / len( seq ) )

#check for repeat segments
def repeat_scan( input, frag_size ):
	if args.verbose >= 1: print( "===== REPEAT FRAGMENT SCAN FOR SIZE: {0} bps =====".format( frag_size ) )
	random.seed()
	
	#determine window and overlap size
	window = int( ( frag_size / 3 ) )
	overlap = ( window - 1 )
	#dummy used so the fragment won't find itself
	dummy = ''.zfill( window * 3 )
	splices = [ input[ i: i + window ] for i in range( 0, len( input ) - 2,  window - overlap ) ]
	
	loop_this = 1
	#this to make sure poly-TRP or poly-MET won't just infinite loop
	loop_break = 0
	
	while ( loop_this != 0 ) and ( loop_break < ( window * 10 ) ):
		#loop this if any mutation is made, until it passes through both checks without mutations
		loop_this = 0
		#loop through each segment
		for segment in splices:
			#check if it's 3 consecutive same codons
			if segment[0][1] == segment[1][1] and segment[1][1] == segment[2][1]:
				if args.verbose >= 2: print( "first 3 codons are identical: {0}".format( segment ) )
				num_mut = random.randint( 1, 2 )
				if args.verbose >= 3: print( "roll dice, mutate {0} codons".format( num_mut ) )
				for f in range( 0, num_mut ):
					position = random.randint( 0, num_mut - 1 )
					apply_mutation( segment, position )
				loop_this += 1
			#check if the segment is found in the full sequence
			seq = ''.join( codon for innerlist in input for codon in innerlist[1] )
			target = ''.join( codon for innerlist in segment for codon in innerlist[1] )
			seq = seq.replace( target, dummy, 1 )
			if seq.find( target ) != -1:
				if args.verbose >= 2: print( "Repeat fragment found with segment: {0}".format( segment ) )
				num_mut = random.randint( 1, ( frag_size / 3 ) - 1 )
				if args.verbose >= 3: print( "roll dice, mutate {0} codons".format( num_mut ) )
				for f in range( 0, num_mut ):
					position = random.randint( 0, ( frag_size / 3 ) - 1 )
					apply_mutation( segment, position )
				loop_this += 1
		loop_break += 1
	return input

##########################################################
#
#	ACTUAL SCRIPT STARTS HERE
#
##########################################################

if args.verbose >= 1: print( "===== SCRIPT START =====" )

#read the input sequence file and parse lines
sequences = []
with open( args.input, 'r' ) as input_file:
	for count, line in enumerate( input_file ):
		if count % 2 == 0:
			seq_name = line
			seq_name = seq_name.replace('\n', '')
		else:
			seq = line
			seq = seq.replace('\n', '')
			sequences.append( [ seq_name, seq ] )
if args.verbose >= 1: print( 'Total number of sequences: {0}'.format( len(sequences) ) )
if args.verbose >= 2: print( sequences )

#generate host profile
host_profile = process_host_table()
#parse data input files
dir = os.path.dirname(__file__)
filename = os.path.join(dir, 'template_files/restriction_sites.txt'.format( args.host ) )
if os.path.exists( filename ) == 1: restrict_sites = parse_unwanted_site_file( filename, 'restriction' )
else: restrict_sites = []
filename = os.path.join(dir, 'template_files/start_sites.txt'.format( args.host ) )
if os.path.exists( filename ) == 1: start_sites = parse_unwanted_site_file( filename, 'start' )
else: start_sites = []

#process through all supplied sequences
for count, seq in enumerate( sequences ):
	input_seq = sequences[ count ][ 1 ]

	if args.verbose >= 1: print( '===== PROCESSING SEQUENCE {0} ===== {1}'.format( count+1, sequences[ count ][0] ) )

	#check input seq style
	if args.type == "AA":
		if args.verbose >= 1: print( "INPUT IS AA SEQUENCE" )
		input_dna = reverse_translate( input_seq )
	elif args.type == "DNA":
		if args.verbose >= 1: print( "INPUT IS DNA SEQUENCE" )
		input_dna = translate_input( input_seq )
	
	if args.verbose >= 2: dump_sequence( input_dna )
	
	#count codons and current profile
	count_table = count_codons( input_dna )
	input_profile = calc_profile( count_table )
	
	#compare input and host profiles
	mutation_table, diff = compare_profiles( input_profile, host_profile, 1 )

	#run optimization
	difference=1
	cycles_current=0
	relax = 1	
	
	#keep running while there are cycles AND difference between current and host is less than the % relax allowed
	while ( ( cycles_current < args.cycles ) or ( args.cycles == 0 ) ) and ( difference >= ( relax - 1 ) ):
		if args.verbose >= 1: print( "~~~~~~~~~~ Current cycle: {0}/{1} ~~~~~~~~~~".format( cycles_current+1, args.cycles ) )
		#determine how much to relax harmonization
		relax = 1+( args.max_relax * ( ( cycles_current ) / args.cycles ) )		
		if args.verbose >=1: print( "Relax coeff: {0}".format( relax ) )
		
		#mutate residues to match host profile
		input_dna = optimize_sequence( input_dna, mutation_table )
		#check GC content in window
		gc_scan( input_dna, 20, 0.15, 0.90 )	#IDT values
		gc_scan( input_dna, 50, 0.15, 0.80 )	#twist values
		gc_scan( input_dna, 100, 0.28, 0.76 )	#IDT values, but close to twist values (100, 0.25, 0.75)
		gc_scan( input_dna, len( input_dna )*3, 0.3, 0.65 )	#twist values
		#check for unwanted restriction sites
		if len( restrict_sites ) != 0: remove_restriction_sites( input_dna, restrict_sites )
		#check for alternative start sites
		if len( start_sites ) != 0: remove_start_sites( input_dna, start_sites, 'ATG' )
		if len( start_sites ) != 0: remove_start_sites( input_dna, start_sites, 'GTG' )
		if len( start_sites ) != 0: remove_start_sites( input_dna, start_sites, 'TTG' )
		#check for repeat fragments
		repeat_scan( input_dna, 9 )
		#check for local homopolymers
		remove_local_homopolymers( input_dna )
		
		#count codons and current profile
		count_table = count_codons( input_dna )
		input_profile = calc_profile( count_table )
		
		#compare input and host profiles
		mutation_table, diff = compare_profiles( input_profile, host_profile, relax )
		difference=diff
			
		#tick cycle
		cycles_current += 1
	
	#hit the max number of cycles?
	if cycles_current == args.cycles:
		if args.verbose >= 1: print( "You hit the max number of cycles: {0}".format( args.cycles ) )

	#check GC content
	if args.verbose >= 1: print( "===== GC CONTENT =====" )
	gc_percent = round( gc_content( input_dna ), 3 )
	if gc_percent < 0.3 or gc_percent > 0.65: print( "WARNING: total GC content is {0}!".format( gc_percent ) )

	#dump result name and sequence
	dump_name( input_dna )
	#dumps a decimal of final difference (0.00 is ideal)
	if args.verbose <= 0: print( "({0})".format( round( difference, 2 ) ), end=" " )
	dump_sequence( input_dna )
#end
