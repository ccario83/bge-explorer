#!/usr/bin/python
# Clinton Cario 6/20/2014
# System libraries
import os, sys, datetime, re
# DB libraries
from peewee import *
from BGC_models import *
# For pretty printing 
from pprint import pprint

# Set the script's directory to its directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))
cwd = os.getcwd()
#SPYDER2 override: cwd = "/Users/ipqb/Documents/Fischbach Rotatation/scripts"
# Some defined paths
DATA_DIR       = os.path.abspath(os.path.join(cwd, os.pardir, "data"))
ANNOT_FILE     = os.path.join(DATA_DIR, "BGC_annotation_final.out")
ORGS_FILE      = os.path.join(DATA_DIR, "AllOrgsTable.out")
BGCS_DIR       = os.path.join(DATA_DIR, "HumanBGCs")
BGCS_SEQS_DIR  = os.path.join(DATA_DIR, "HumanBGCs_seqs")
ERROR_LOG      = "BGC_parse_errors.log"

# If the 'reset' flag was given, drop and recreate base tables 
if len(sys.argv)>1 and sys.argv[1] == "reset":
	print "Resetting Database Tables..."
	try:
		os.remove(DB_FILE)
	except OSError:
		pass
	db.connect()
	# Create main tables
	Organism.create_table()
	BGC.create_table()
	Contig.create_table()
	Gene.create_table()
	Domain.create_table()
	Sequence.create_table()
	Similarity.create_table()
	# Create Source Tables
	Source.create_table()
	OrganismSource.create_table()
	BGCSource.create_table()
	ContigSource.create_table()
	GeneSource.create_table()
	DomainSource.create_table()
	SequenceSource.create_table()
	# Create Success tables and entries
	Success.create_table()
	Success.create(step='First')
	Success.create(step='Second')
	Success.create(step='Third')
	Success.create(step='Fourth')
	Success.create(step='Fifth')

# Open a handle for the error log file
errfile = open(ERROR_LOG,'w')


print "Parsing Data..."

### ORGS_FILE
### HEADER AND EXAMPLE
### ==================
# taxon_oid								644736411
# Domain 								Archaea
# Status 								Finished
# Study Name 							Thermococcus gammatolerans EJ3
# Genome Name / Sample Name 			Thermococcus gammatolerans EJ3
# Sequencing Center 		 			Genoscope
# Phylum 								Euryarchaeota
# Class 								Thermococci
# Order 								Thermococcales
# Family 								Thermococcaceae
# Genus 								Thermococcus
# Species 								gammatolerans
# IMG Genome ID  						644736411
# Genome Size  							2045438
# Gene Count  							2206
print "\tParsing organism information..."
# Open the file and strip the header
if (Success.get(step='First').successful == False):
	try:
		Organism.drop_table()
	except:
		pass
	Organism.create_table()
	errfile.write("Organism Parsing errors\n=======================")
	infile = open(ORGS_FILE,'r')
	infile.next()
	line_no = 0
	with db.transaction():
		for line in infile:
			entry = line.strip("\n").split("\t")
			line_no = line_no + 1
			try:
				source = Source.create(location = "%s|%s" %(os.path.basename(ORGS_FILE), line_no))
				# These are unique, don't have to check if already exists
				organism = Organism.create(
											_domain				= entry[1],
											_phylum				= entry[6],
											_class				= entry[7],
											_order				= entry[8],
											_family				= entry[9],
											_genus				= entry[10],
											_species			= entry[11],
											taxon_id			= entry[0],
											status				= entry[2],
											sequencing_center	= entry[5],
											img_genome_id		= entry[12],
											genome_name			= entry[4],
											genome_name_key		= entry[4].replace(" ","_"),
											genome_size			= entry[13],
											gene_count			= entry[14],
										   )
				OrganismSource.create(source=source, organism=organism)

			except IndexError:
				#print "\t\tError parsing entry on line #%d, beginning '%s...'" % (line_no, line[0:20].replace("\n",""))
				errfile.write("\n[%d]: %s" % (line_no, line.replace("\n","")))
	infile.close()
	step = Success.get(step='First')
	step.successful = True
	step.save()
else:
	print "\t(Table marked as finished from previous run, skipping...)"


### ANNOT_FILE
### HEADER AND EXAMPLE
### ==================
#	[contig_id]_[cluster_idx]			637000276_7
#	kind								saccharide
#	organism_name						Bacillus_thuringiensis_sv_konkukian_97-27
#	start_locus_tag						BT9727_0525
#	end_locus_tag						BT9727_0525
print "\tParsing cluster annotation information..."
if (Success.get(step='Second').successful == False):
	try:
		BGC.drop_table()
	except:
		pass
	BGC.create_table()
	errfile.write("\n\nCluster Annotation Parsing errors\n=================================")
	# Open the file and strip the header
	infile = open(ANNOT_FILE,'r')
	line_no = 0
	with db.transaction():
		for line in infile:
			line_no = line_no + 1
			if (line_no % 1000)==1: print "\t\tOn line %d of ~14,000" % (line_no)
			entry = line.strip("\n").split("\t")
			try:
				source = Source.create(location = "%s|%s" %(os.path.basename(ANNOT_FILE), line_no))
				contig_id, cluster_idx  = entry[0].split("_")
				# Make a new contig with this id
				contig = None
				
				try: 
					contig = Contig.get(Contig.contig_id == contig_id)
				except Contig.DoesNotExist:
					contig = Contig(contig_id=contig_id, source=source)
					contig.save(force_insert=True)
					ContigSource.create(source=source, contig=contig)
				# Look up the organism
				organism = None
				
				try:
					organism = Organism.get(Organism.genome_name_key == entry[2])
				except Organism.DoesNotExist:
					#print "\t\tError parsing entry on line #%d, organism  '%s' doesn't exist" % (line_no, entry[2][0:20].replace("\n",""))
					errfile.write("\n[%d]: %s" % (line_no, line.replace("\n","")))
					continue
				
				# Create the BGC (these are unique, don't have to check if already exists)
				bgc = BGC.create(
									cluster_id 			= entry[0],
									contig				= contig,
									cluster_idx			= cluster_idx,
									organism			= organism,
									kind				= entry[1],
									start_locus_tag		= entry[3],
									end_locus_tag		= entry[4],
								)
				BGCSource.create(source=source, bgc=bgc)
			except ValueError:
				try:
					#print "\t\tError parsing entry on line #%d, beginning '%s...'" % (line_no, line[0:20].replace("\n",""))
					errfile.write("\n[%d]: %s" % (line_no, line.replace("\n","")))
				except:
					print sys.exc_info()[0]
	infile.close()
	step = Success.get(step='Second')
	step.successful = True
	step.save()
else:
	print "\t(Table marked as finished from previous run, skipping...)"



### BGCS_DIR .out files
### HEADER AND EXAMPLE
### ==================
#	contig ID (IMG JGI)		2502205815
#	img_gene ID (IMG JGI) 	2502290047
#	gene start (nt) 		106
#	gene end (nt) 			864
#	domain start (aa) 		61
#	domain end (aa) 		238
#	pfam ID (Pfam) 			PF01695 (Q for proteins with no Pfam annotations)
#	ClusterFinder prob 		0.520073176752
print "\tParsing Gene and Domain information..."
if (Success.get(step='Third').successful == False):
	try:
		Gene.drop_table()
		Domain.drop_table()
	except:
		pass
	Gene.create_table()
	Domain.create_table()
	errfile.write("\n\nGene/Domain (*.out) Parsing errors\n==================================")
	file_no = 0
	total_file_no = len(os.listdir(BGCS_DIR))/2
	for _file in os.listdir(BGCS_DIR):
		if _file.endswith("clusters.out"):
			file_no = file_no+1
			if (file_no % 1000)==1: print "\t\tOn file %s... [%d of %d]" % (_file, file_no, total_file_no)
			with open(os.path.join(BGCS_DIR,_file),'r') as infile:
				errfile.write("\n\n\tFile %s\n\t--------------------" % _file)
				line_no = 0
				with db.transaction():
					for line in infile:
						entry = line.strip("\n").split("\t")
						line_no = line_no + 1
						try:
							source = Source.create(location = "%s|%s" %(_file, line_no))
							# Get the contig or create (also saving the source)
							contig = None
							try: 
								contig = Contig.get(Contig.contig_id == entry[1])
							except Contig.DoesNotExist:
								contig = Contig(contig_id=entry[1])
								contig.save(force_insert=True)
								ContigSource.create(source=source, contig=contig)
							# Get the cluster
							bgc = None
							try: 
								bgc = BGC.get(BGC.cluster_id == entry[0])
							except BGC.DoesNotExist:
								# Create an entry, but doesn't have any annotation
								bgc = BGC.create(
												cluster_id 			= entry[0],
												contig				= contig,
												cluster_idx			= entry[0].split("_")[1],
												organism			= None,
												kind				= None,
												start_locus_tag		= None,
												end_locus_tag		= None,
											)
								BGCSource.create(source=source, bgc=bgc)

							# Get the gene or create it
							gene = None
							try:
								gene = Gene.get(
													Gene.img_gene_id	== entry[2],
													Gene.contig			== contig,
													Gene.start			== int(entry[3]),
													Gene.end			== int(entry[4])
												)
							except Gene.DoesNotExist:
								gene = Gene.create(
													img_gene_id			= entry[2],
													contig				= contig,
													bgc					= bgc,
													start				= int(entry[3]),
													end					= int(entry[4])
												  )
								GeneSource.create(source=source, gene=gene)
							
							# Get or create the domain
							try:
								Domain.get(
													Domain.gene					== gene,
													Domain.pfam_id				== entry[7],
													Domain.start				== int(entry[5]),
													Domain.end					== int(entry[6]),
													Domain.cluster_finder_prob 	== float(entry[8])
										  )
							except Domain.DoesNotExist:
								domain = Domain.create(
													gene						= gene,
													pfam_id						= entry[7],
													bgc							= bgc,
													start						= int(entry[5]),
													end							= int(entry[6]),
													cluster_finder_prob 		= float(entry[8])
											 )
								DomainSource.create(source=source, domain=domain)
						except:
							try:
								print "\t\t\tError parsing entry on line #%d of %s, beginning '%s'" % (line_no, _file, line[0:20].replace("\n",""))
								errfile.write("\n\t[%d]: %s" % (line_no, line.replace("\n","")))
							except:
								print sys.exc_info()[0]
		else:
			# A more comprehensive list of genes and domains can be found in the .out files
			continue
	step = Success.get(step='Third')
	step.successful = True
	step.save()
else:
	print "\t(Table marked as finished from previous run, skipping...)"


# select * from bgc, organism where bgc.kind = "saccharide" and organism._genus = "Bacteroides" now works
# > 836280 results


### BGCS_SEQS_DIR [X].fasta files
### HEADER AND EXAMPLE
### ==================
#	img_gene_id 		= 2500069248
#	locus_tag 			= Sugar phosphate permease
#	organism			= [Lactobacillus reuteri 100-23]
#...X					= cluster_id ([contig_id]_[contig_idx])
print "\tParsing Sequence information..."
if (Success.get(step='Fourth').successful == False):
	try:
		Sequence.drop_table()
	except:
		pass
	Sequence.create_table()
	errfile.write("\n\nSequence (*.fasta) Parsing errors\n=================================")
	file_no = 0
	total_file_no = len(os.listdir(BGCS_SEQS_DIR))
	for _file in os.listdir(BGCS_SEQS_DIR):
		if _file.endswith(".fasta"):
			file_no = file_no+1
			if (file_no % 1000)==1: print "\t\tOn file %s... [%d of %d]" % (_file, file_no, total_file_no)
			with open(os.path.join(BGCS_SEQS_DIR,_file),'r') as infile:
				errfile.write("\n\n\tFile %s\n\t--------------------" % _file)
				line_no = 0
				line = " "
				with db.transaction():
					# Define variables
					sequence 	= None
					gene 		= None
					locus_tag 	= None
					organism 	= None
					# Now cycle through the file
					while line:
						line_no = line_no + 1
						line = infile.readline()
						# See if this line is a fasta header (it is initially)
						fasta_regex = re.match('^>(\d+) (.*) \[(.*)\]$', line)
						# The line will either match this or be a sequence
						if fasta_regex:
							#print "\t\t\t----------------------------"
							#print "\t\t\tHeader match:\t%s" % fasta_regex.group()[0:40]
							# If a fasta line, try to get the gene and organism or skip
							source = Source.create(location = "%s|%s" %(_file, line_no))
							try:
								gene = Gene.get(Gene.img_gene_id == fasta_regex.group(1))
							except Gene.DoesNotExist:
								#print "\t\t\tError parsing gene on header line #%d of %s: '%s'" % (line_no, _file, fasta_regex.group(1))
								errfile.write("\n\t[%d]: %s" % (line_no, line.replace("\n","")))
								continue
							try:
								organism = Organism.get(Organism.genome_name == fasta_regex.group(3))
							except Organism.DoesNotExist:
								#print "\t\t\tError parsing organism on header line #%d of %s: '%s'" % (line_no, _file, fasta_regex.group(3))
								errfile.write("\n\t[%d]: %s" % (line_no, line.replace("\n","")))
								continue
							#print "\t\t\tGene:       \t%5d %s\n\t\t\tOrganism:   \t%5d %s" % (gene.gene_id, gene.img_gene_id, organism.organism_id, organism.genome_name)
							# Get the sequence
							sequence = ""
							line = infile.readline()
							while not re.match('^>(\d+) (.*) \[(.*)\]$', line) and line:
								sequence = sequence + line.replace("\n","")
								line = infile.readline()
							# Create a new sequence
							try:
								Sequence.get(
															gene 		= gene,
															locus_tag 	= fasta_regex.group(2),
															organism	= organism,
															sequence 	= sequence
											)
							except Sequence.DoesNotExist:
								sequence = Sequence.create(
															gene 		= gene,
															locus_tag 	= fasta_regex.group(2),
															organism	= organism,
															sequence 	= sequence
														  )
								SequenceSource.create(source=source, sequence=sequence)
							#print "\t\t\tSequence:  \t%s...%s" % (sequence[0:5], sequence[-6:-1])
							# Reset variables
							sequence 	= None
							gene 		= None
							locus_tag 	= None
							organism 	= None
	step = Success.get(step='Fourth')
	step.successful = True
	step.save()
else:
	print "\t(Table marked as finished from previous run, skipping...)"


# Add domain sequences to the DB
print "\tGenerating domain sequences..."
if (Success.get(step='Fifth').successful == False):
	errfile.write("\n\nDomain sequence computation errors\n==================================")
	entry_no = 0
	total_entry_no = Domain.select().count()
	with db.transaction():
		for domain in Domain.select():
			entry_no = entry_no+1
			if (entry_no % 1000)==1: print "\t\tOn entry: %d of %d" % (entry_no, total_entry_no)
			try:
				domain.sequence = domain.gene.protein_seq.get().sequence[(domain.start-1):(domain.end)]
				domain.save()
			except Sequence.DoesNotExist:
				pass
			if domain.gene.bgc:
				domain.bgc=domain.gene.bgc
				domain.save()
	step = Success.get(step='Fifth')
	step.successful = True
	step.save()
else:
	print "\t(Table marked as finished from previous run, skipping...)"


## Import CalculateDistanceSeq.py related data here?
errfile.close()

