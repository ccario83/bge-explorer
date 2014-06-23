from peewee import *
import os
import datetime

os.chdir(os.path.dirname(os.path.abspath(__file__)))
cwd = os.getcwd()
DATA_DIR	= os.path.abspath(os.path.join(cwd, "data"))
DB_FILE 	= os.path.join(DATA_DIR, "BGCs.db")
db = SqliteDatabase(DB_FILE)

class BGCModel(Model):
	class Meta:
		commit_select = True
		database = db


class Organism(BGCModel):
	organism_id			= PrimaryKeyField()
	_domain				= CharField()
	_phylum				= CharField()
	_class				= CharField()
	_order				= CharField()
	_family				= CharField()
	_genus				= CharField(index=True)
	_species			= CharField()
	taxon_id			= CharField()
	status				= CharField()
	sequencing_center	= CharField()
	img_genome_id		= CharField()
	genome_name			= CharField()
	genome_name_key		= CharField()
	genome_size			= IntegerField()
	gene_count			= IntegerField()
	#created_at			= DateTimeField(default=datetime.datetime.now)


class Contig(BGCModel):
	contig_id			= CharField(primary_key=True)
	#created_at			= DateTimeField(default=datetime.datetime.now)

# The ids are a bit confusing here...
# bgc_id is how the clusters are identified and referenced in this database
# cluster_id is how clusters are identified in the dataset which consists of ([contig_id]_[cluster_idx])
class BGC(BGCModel):
	bgc_id				= PrimaryKeyField()
	cluster_id 			= CharField()
	contig				= ForeignKeyField(Contig, related_name="contigs")
	cluster_idx			= IntegerField()
	organism			= ForeignKeyField(Organism, related_name="clusters", null=True)
	kind				= CharField(index=True, null=True)
	start_locus_tag		= CharField(null=True)
	end_locus_tag		= CharField(null=True)
	#created_at			= DateTimeField(default=datetime.datetime.now)


class Gene(BGCModel):
	gene_id				= PrimaryKeyField()
	img_gene_id 		= CharField()
	contig				= ForeignKeyField(Contig, related_name="genes", null=True)
	bgc					= ForeignKeyField(BGC, related_name="genes", null=True)
	start				= IntegerField()
	end					= IntegerField()
	#created_at			= DateTimeField(default=datetime.datetime.now)


class Domain(BGCModel):
	domain_id			= PrimaryKeyField()
	gene				= ForeignKeyField(Gene, related_name="domains")
	bgc					= ForeignKeyField(BGC, related_name="domains", null=True)
	pfam_id				= CharField(index=True)
	start				= IntegerField()
	end					= IntegerField()
	sequence 			= TextField(null=True)
	cluster_finder_prob = FloatField()
	#created_at			= DateTimeField(default=datetime.datetime.now)


class Sequence(BGCModel):
	sequence_id			= PrimaryKeyField()
	gene 				= ForeignKeyField(Gene, related_name="protein_seq")
	locus_tag 			= CharField()
	organism			= ForeignKeyField(Organism, related_name="organism")
	sequence 			= TextField()
	created_date 		= DateTimeField(default=datetime.datetime.now)


# Table for similarities
class Similarity(BGCModel):
	similarity_id		= PrimaryKeyField()
	comparing			= CharField()
	id1					= IntegerField()
	id2					= IntegerField()
	measure				= CharField(index=True)
	value				= CharField(null=True)
	created_date 		= DateTimeField(default=datetime.datetime.now)
	
	class Meta:
		indexes = (
					(('id1', 'id2', 'comparing'), False),
				  )

# Table for extra NCBI genbank annotations
class NCBIseq(BGCModel):
	ncbiseq_id 			= PrimaryKeyField()
	bgc 				= ForeignKeyField(BGC, related_name="extras", null=True)
	description 		= CharField(null=True)
	dbxrefs 			= CharField(null=True)
	sequence_id 		= CharField()
	sequence 			= TextField(null=True)
	# A json object
	annotations  		= TextField(null=True)
	# keys:
	#comment				= 
	#sequence_version    = 
	#source
	#taxonomy
	#keywords
	#references
	#accessions
	#data_file_division
	#date
	#organism
	#gi
	# 

class NCBIcds(BGCModel):
	nbcicds_id 			= PrimaryKeyField()
	ncbiseq 			= ForeignKeyField(NCBIseq, related_name="CDS", null=True)
	bgc 				= ForeignKeyField(BGC, related_name="ncbi_domains", null=True)
	gene 				= ForeignKeyField(Gene, related_name="fuzzy_ncbi_match", null=True)
	start				= IntegerField()
	end 				= IntegerField()
	strand 				= IntegerField()
	sequence 			= TextField(null=True)
	locus 				= CharField(null=True)
	description			= TextField(null=True)
	notes 				= TextField(null=True)

# Source tables indicate where data was originally found
class Source(BGCModel):
	source_id			= PrimaryKeyField()
	location			= CharField()

class OrganismSource(BGCModel):
	organism			= ForeignKeyField(Organism)
	source				= ForeignKeyField(Source)

class ContigSource(BGCModel):
	contig				= ForeignKeyField(Contig)
	source				= ForeignKeyField(Source)

class BGCSource(BGCModel):
	bgc					= ForeignKeyField(BGC)
	source				= ForeignKeyField(Source)

class GeneSource(BGCModel):
	gene				= ForeignKeyField(Gene)
	source				= ForeignKeyField(Source)

class DomainSource(BGCModel):
	domain				= ForeignKeyField(Domain)
	source				= ForeignKeyField(Source)

class SequenceSource(BGCModel):
	sequence			= ForeignKeyField(Sequence)
	source				= ForeignKeyField(Source)



# Used to keep track of population progress and to avoid repeating successful steps
class Success(BGCModel):
	success_id			= PrimaryKeyField()
	step 				= CharField()
	successful 			= BooleanField(default = False)





