#!/usr/bin/python
# Clinton Cario 6/20/2014
# OS, file/path manipulation libraries
import os, subprocess, tempfile, threading, shutil
# Maths
import math
import numpy as np
from munkres import Munkres
# DB libraries
from peewee import *
from BGC_models import *
# Webby stuff
import json
import webbrowser
import SimpleHTTPServer
import SocketServer

# Arrower (modified version from Peter)
from Arrower import *

# To augment from NCBI
from Bio import Entrez, SeqIO
Entrez.email = "nothanks@notgiven.org"

# Debug, print
from pprint import pprint

### Sample filter usage (see filters.txt)
## NRPS   = BGCs(filters={'Kind':"NRPS", 'Species':["%acnes%","%epi%"]})
## BACSAC = BGCs(filters={'Kind':"saccharide", 'Genus':"Bacteroides"})


# Global table definition (used to generate SQL query from filters)
TABLES        = ['Sequence', 'Domain', 'Gene', 'BGC', 'Organism']
LED           = ['domain','phylum','class','order','family','genus','species']
REPLACEMENTS  = {'domain_sequence': ('Domain','sequence'), 'gene_sequence':('Sequence','sequence')}

# Global Paths
OUTDIR        = os.path.abspath(os.path.join(cwd, "data/"))
FASTA_OUTDIR  = os.path.abspath(os.path.join(cwd, "data/BGCs"))
WEB_OUTDIR    = os.path.abspath(os.path.join(cwd, "data/web"))
ARROW_OUTDIR  = os.path.abspath(os.path.join(cwd, "data/arrows"))

# Query class to build SQL database queries from filters and to return results
class BGCQuery():
	table   = None
	columns = "*"
	filters = None
	query   = ""
	verbose = None

	def __init__(self, table='BGC', columns="*", filters=None, verbose=True):
		self.table = table
		if not self.table in TABLES:
			raise Exception('Invalid Table!')
		self.columns = columns
		self.filters = filters
		self.verbose = verbose

	# Returns a join clause syntax (to add to query)
	def join(self, parent, child): 
		return "\nINNER JOIN {parent} ON {parent}.{parent}_id = {child}.{parent}_id".format(parent=parent.lower(), child=child.lower())

	# Recursively joins tables up the BGC table hierarchy 
	def join_up(self, from_child):
		phrase = ""
		idx = TABLES.index(from_child)
		while idx+1 < len(TABLES):
			phrase += self.join(TABLES[idx+1], TABLES[idx])
			idx += 1
		return phrase

	# Turns humanized table column names to their actual database column forms
	def columnize(self, x):
		x = x.lower().replace(" ","_")
		if x in LED:
			x = '_' + x
		owner = None
		for table in TABLES:
			if hasattr(eval(table),x):
				owner = table
		if x in REPLACEMENTS.keys():
			owner = REPLACEMENTS[x][0]
			x     = REPLACEMENTS[x][1]
		if owner is None or x is None:
			raise Exception('Invalid column')
		return owner.lower() + "." + x

	# Returns a where clause syntax based on the filters keys and the type of their values
	def where(self):
		filters=self.filters
		phrase = ""
		lead   = "WHERE"
		for key, vals in filters.iteritems():
			column = self.columnize(key)
			if type(vals) == int:
				phrase += "\n{lead} {column} = {value}".format(lead=lead, column=column, value=vals)
			elif type(vals) == str:
				if vals[0:2] == "IS":
					phrase += "\n{lead} {column} {value}".format(lead=lead, column=column, value=vals)
				else:
					phrase += "\n{lead} {column} LIKE '{value}'".format(lead=lead, column=column, value=vals)
			elif type(vals) == list:
				if any([True for val in vals if '%' in str(val)]):
					like_phrase = "\" OR {column} LIKE \"".format(column=column).join([str(val) for val in vals])
					phrase += "\n{lead} ({column} LIKE \"{like}\")".format(lead=lead, column=column, like=like_phrase)
				else:
					phrase += "\n{lead} {column} IN ('{vals}')".format(lead=lead, column=column, vals="','".join([str(val) for val in vals]))
			lead = "AND"
		return phrase

	# Creates a select clause based on requested columns
	def select(self):
		table   = self.table
		cols    = self.columns
		phrase  = "SELECT {cols} FROM {table}".format(cols=cols, table=table.lower())
		return phrase

	# Builds the SQL query using above functions
	def build(self, verbose=None):
		verbose = verbose if verbose != None else self.verbose
		table   = self.table
		query   = self.select()
		query  += self.join_up(table)
		if self.filters: query += self.where()
		if self.verbose: print "Generated query:\n" + query
		self.query = query

	# Executes the SQL query generated above
	def execute(self, just_ids=False, as_dict=True, verbose=None):
		verbose = verbose if verbose != None else self.verbose
		if self.query == "": self.build(verbose=verbose)
		results = eval(self.table).raw(self.query).dicts().execute()
		if just_ids:
			results = [result[self.table.lower()+'_id'] for result in results]
		elif as_dict:
			results = dict([(result[self.table.lower()+'_id'],result) for result in results])
		else:
			results = [result for result in results]
		if verbose: print "Query returned %d results" % len(results)
		return results 

	# Convenience functions
	def get_bgcs(self, just_ids=False, as_dict=False, verbose=False):
		verbose = verbose if verbose != None else self.verbose
		''' Get all BGCs using the given filters '''
		self.table = "BGC"
		return self.execute(just_ids=just_ids, as_dict=as_dict, verbose=verbose)

	def get_genes(self, just_ids=False, as_dict=False, verbose=False):
		verbose = verbose if verbose != None else self.verbose
		''' Get all genes using the given filters '''
		self.table = "Gene"
		return self.execute(just_ids=just_ids, as_dict=as_dict, verbose=verbose)

	def get_domains(self, just_ids=False, as_dict=False, verbose=False):
		verbose = verbose if verbose != None else self.verbose
		''' Get all genes using the given filters '''
		self.table = "Domain"
		return self.execute(just_ids=just_ids, as_dict=as_dict, verbose=verbose)


# This is the man BGC class, which filters the database based on given parameters, 
# compares BGCS, visualizes results, and much more. 
# Descriptions of each function can be found in their definitions
class BGCs():
	# Whether to display output to the screen (True by default)
	verbose      = None
	# A handle to the actual DB filter (which is a BGCQuery object, as defined above)
	db           = None
	# A structure of the network after clustering ('nodes' and 'edges' are in a format liked by d3)
	network      = None  # {'nodes': [{ 'name' => bgc_id, 'attribs' => bgc }], 'edges': [{ 'source' => node idx, 'target' => node idx, 'value': edge_weight }], 'mapping': {bgc_id => node_idx} }
	# A structure of the clusters after clustering
	clusters     = None  # {'nodes': [{ 'cluster' => group_idx, 'attribs' => bgc, 'radius' => 10}], 'groups': [[{bgc_id => bgc}]], 'mapping': [{bgc_id => group_idx}] }

	def __init__(self, filters=None, verbose=True):
		self.verbose = verbose
		self.db      = BGCQuery(table="BGC", filters=filters, verbose=verbose)
		if verbose: print "If this is a new filter, please be sure to run .compare_bgcs(). This can take a while, but only has to be done once!"
		if verbose: print "You can then cluster with .cluster() and visualize with .visualize_network() or .visualize_clusters()"

	# Using the filters, performs the clustering using the given parameters
	# Results are stored in self.network, and self.clusters
	def cluster(self, measure="similarity", cutoff=0.35, use_mcl=True, I=3.5, prefilter=True, keep_singletons=False, verbose=True):
		''' Using the query to filter the database (stored in self.db.query), this will get all BGCs and 
		    cluster them using the given parameters. This only has to be run once per filter set.

		       measure:          Can be 'similarity' (combined), 'jaccard', or 'dds' measures 
		       cutoff:           If not 0, will remove any edges from the network with weight less than this value (prefilter = False, this gets set to 0)
		       use_mcl:          Whether to use mcl to perform clustering
		       I:                The mcl inflation value used if use_mcl is used
		       prefilter:        Whether to remove edges with weights less than 'cutoff' before mcl
		       keep_singletons:  Whether clusters of n=1 are kept for the visualization
		       verbose:          Display clustering information to the screen?
		'''
		# Check inputs
		if measure not in ['similarity', 'jaccard', 'dds']:
			print "Please use similarity, jaccard, or dds as measures"
			return
		if cutoff<0:
			print "Please use a cutoff greater than 0"
			return
		if type(use_mcl) != type(True):
			print "Please use True or False for the use_mcl parameter"
			return
		if I < 0:
			print "Please use a positive I value for mcl"
			return
		if type(prefilter) != type(True):
			print "Please use True or False for the prefilter parameter"
			return
		if type(keep_singletons) != type(True):
			print "Please use True or False for the keep_singletons parameter"
			return
		# Test verbose for True/False?? Eh,... no.

		verbose = verbose if verbose != None else self.verbose
		# Get all bgcs using the filters and define variables used to start building the json (for webpage viz) object
		bgcs           = self.db.get_bgcs(as_dict=True)
		bgc_ids        = self.db.get_bgcs(just_ids=True)
		self.network   = {'nodes': [], 'edges': [], 'mapping': {}}
		self.clusters  = {'nodes': [], 'groups': [], 'mapping': {}}
		if verbose: print "Found %d nodes" % len(bgc_ids)
		
		
		# Get and filter edges
		# --------------------------
		# Get all pairwise similarity scores for these ids (edges)
		if use_mcl and not prefilter:
			cutoff = 0.0
		if verbose: print "Getting edges with weight > %.2f" % cutoff
		edges = self.get_bgc_edges(measure, cutoff, verbose=False)
		if verbose: print "Found %d edges" % len(edges)
		if len(edges) == 0:
			print "Warning: No edges found. Please relax cutoffs and check filters"
			return 
		# Use mcl to trim edges if requested
		if use_mcl:
			if verbose: print "Pruning edges with MCL using I = %.2f" % I
			orig = len(edges)
			edges = self._mcl(bgcs, edges, I=I, verbose=False)
			if verbose: print "Removed %d edges" % (orig - len(edges))
			# At this point, edges looks like this: (bgc1,bgc2) => edge weight


		# Initialize the mappings from bgc id to what cluster/network number they are in (to None)
		for id_ in bgcs.keys():
			#self.network['mapping'][id_] = None
			self.clusters['mapping'][id_] = None

		# Iterate over bgc pairs and create edges between them if one exists
		# For all found edges, also build networks and clusters. 
		if verbose: print "Building network and clusters"
		cur_clust_no = -1
		cur_node_no   = 0
		bgc_pairs    = self._pairwise(bgc_ids, self_compare=False)
		for b1, b2 in bgc_pairs:
			# b1 always has the smaller ID in the edge object
			b1, b2 = sorted([b1, b2])
			edge = None
			# try to find the edge between these two bgcs or continue
			try:
				edge = edges[(b1,b2)]
			except:
				pass

			if (edge == None and keep_singletons) or edge != None:
				# Append the nodes and their attributes to the network graph, if needed 
				for id_ in [b1,b2]:
					try:
						self.network['mapping'][id_]
					except KeyError:
						self.network['nodes'].append({'name': id_, 'attribs': bgcs[id_]})
						self.network['mapping'][id_] = cur_node_no
						cur_node_no += 1

			if edge == None:
				continue

			# Add the edge to the network, using source/target indexes instead of IDs (as required by D3). 'value' is the edge weight
			self.network['edges'].append({'source': self.network['mapping'][b1], 'target': self.network['mapping'][b2], 'value': edge})

			# Determine cluster assignment. 'mapping' is a dictionary of node id => cluster number
			# If neither bgc belongs to a cluster (both none), then 
			if self.clusters['mapping'][b1] == None and self.clusters['mapping'][b2] == None:
				# this is a new cluster
				cur_clust_no += 1
				# and both ids should map to this cluster number
				self.clusters['mapping'][b1] = cur_clust_no
				self.clusters['mapping'][b2] = cur_clust_no
				# also add both bgcs {bgc_id => attribs} to 'groups' (a convenience object)
				self.clusters['groups'].append({b1: bgcs[b1], b2: bgcs[b2]})
				# One of the nodes (or both) are in a known cluster
			elif self.clusters['mapping'][b1] != None or self.clusters['mapping'][b1] != None:
				# so add the new node to it
				if self.clusters['mapping'][b1] == None:
					# The first bgc has no cluster
					# so add cluster idx to the mapping for b1, and then add b1 to b2's cluster
					b2_cluster = self.clusters['mapping'][b2]
					self.clusters['mapping'][b1] = b2_cluster
					self.clusters['groups'][b2_cluster][b1] = bgcs[b1]
				elif self.clusters['mapping'][b2] == None:
					# Same as above, but for b2
					b1_cluster = self.clusters['mapping'][b1]
					self.clusters['mapping'][b2] = b1_cluster
					self.clusters['groups'][b1_cluster][b2] = bgcs[b2]
				elif self.clusters['mapping'][b1] != None and self.clusters['mapping'][b1] != None:
					# Both bgcs belong to sub clusters, merge them
					b1_cluster = self.clusters['mapping'][b1]
					b2_cluster = self.clusters['mapping'][b2]
					# Change the mapping of every bgc in cluster2 to cluster1
					for b in self.clusters['groups'][b2_cluster]:
						self.clusters['mapping'][b] = b1_cluster
					# Update b1 cluster to include all of b2's bgcs
					self.clusters['groups'][b1_cluster].update(self.clusters['groups'][b2_cluster])
				else:
					raise Exception('Invalid cluster membership assignment and/or transfer. This should never happen.')
		
		# If requested, add singleton nodes to the mapping and build new cluster group objects
		if keep_singletons:
			for node, cluster in self.clusters['mapping'].iteritems():
					if cluster == None:
						cur_clust_no += 1
						self.clusters['mapping'][node] = cur_clust_no
						self.clusters['groups'].append({node: bgcs[node]})
		# Or remove network nodes with no edges
		elif verbose:
			print "Discarded %d singleton nodes from cluster visualization" % len([n for n,c in self.clusters['mapping'].iteritems() if c==None])


		# Build clusters node object for D3 cluster layout
		for node, cluster in self.clusters['mapping'].iteritems():
			if cluster == None: continue
			self.clusters['nodes'].append({'cluster': cluster, 'attribs': bgcs[node], 'radius': 1})


	def visualize_network(self, annotation="_species", verbose=None):
		''' Creates a json file that can be accessed in a web browser at 'http://localhost:8000/network_vis.html' after
			'python -m SimpleHTTPServer' is run by self._start_webserver

			annotation: The initial annotation column to use to color nodes when the page is launched
		'''
		if self.network['nodes']==[]:
			print "There are no nodes remaining after clustering to visualize. Please relax cutoffs and check filters"
			return
		verbose = verbose if verbose != None else self.verbose
		# Copy the network and add a 'group' key to each node corresponding to the attribute to group (color) by
		data = self.network.copy()
		data.pop('mapping',None)
		#sfor node in data['nodes']:
		#s	group = node['attribs'][annotation]
		#s	node['group'] = group
		#s	#node.pop('attribs',None)
		attribs = sorted(data['nodes'][0]['attribs'].keys())

		# Write the cluster to disk 
		with open(os.path.join(WEB_OUTDIR,'d3_network_data.json'), 'w') as outfile:
			json.dump({'attribs':attribs, 'inital_group':annotation, 'data':data}, outfile)
		
		# Start a webserver and serve from the data directory (to visualize cluster JSON file with d3)
		self._start_webserver()
		webbrowser.open("http://localhost:8000/data/web/network_vis.html", new=0)
		return


	def visualize_clusters(self, verbose=None, show_domains=True, color_key='pfams', group_nodes=True):
		''' Creates a json file that can be accessed in a web browser at 'http://localhost:8000/cluster_vis.html' after
			'python -m SimpleHTTPServer' is run by self._start_webserver

			show_domains:   Use light gray blocks to indicate domains in the genes shown in the arrower images
			color_key:      A BGC gene attribute (gene table column names) to color arrows by, or 'pfam' to use a in-order-list of PFAM ids as a color key
			group_nodes:    Show all nodes in a clustered group as one bubble, with radius proportional to the number of members of that cluster
		'''
		if self.clusters['nodes']==[]:
			print "There are no nodes remaining after clustering to visualize. Please relax cutoffs and check filters"
			return
		verbose = verbose if verbose != None else self.verbose
		# Make the arrows
		if verbose: print "Generating arrows (this can take a few seconds)"
		self.generate_arrows(show_domains=show_domains, color_key=color_key)
		# Copy the network and add a 'group' key to each node corresponding to the attribute to group (color) by
		if verbose: print "Preparing data structures"
		data = self.clusters.copy()
		data.pop('mapping',None)
		data.pop('groups',None)
		attribs = sorted(data['nodes'][0]['attribs'].keys())

		# Write the cluster to disk 
		with open(os.path.join(WEB_OUTDIR,'d3_cluster_data.json'), 'w') as outfile:
			json.dump({'attribs':attribs, 'data':data}, outfile)
		
		self._start_webserver()
		if group_nodes:
			webbrowser.open("http://localhost:8000/data/web/cluster_vis.html", new=0)
		else:
			webbrowser.open("http://localhost:8000/data/web/cluster_vis_old.html", new=0)
		return


	def augment(self, flush=False, verbose=None):
		verbose = verbose if verbose != None else self.verbose
		''' Augments the database with NCBI information (updates BGCs/domains) 
			This only has to be run once per filter set. Useful for getting more domain info like strand orientation
		'''
		bgcs = self.db.get_bgcs(as_dict=True)
		filters = self.db.filters.copy()

		# Create table if required
		if flush:
			try:
				NCBIseq.drop_table()
				NCBIcds.drop_table()
			except:
				pass
		try:
			NCBIseq.create_table()
			NCBIcds.create_table()
		except:
			pass

		# Some convenience functions
		getCDS = lambda x: [y for y in x.features if y.type == 'CDS'][0]
		getIDs = lambda x: (x.annotations['db_source'].split(" ")[1], x.annotations['accessions'][0])
		getPos = lambda x: [ int(filter(type(y).isdigit,y)) for y in x.qualifiers['coded_by'][0].strip(')').split(":")[1].split("..")]

		# Debug
		#a = [(a,b) for a,b in bgcs.iteritems()]
		#key,local_bgc = [z for z in a if z[0]==10545][0]
		for key,local_bgc in bgcs.iteritems():
			if verbose: print "Processing bgc with id %d" % key
			# Determine if this entry exists
			try:
				bgc_record = NCBIseq.get(NCBIseq.bgc == key)
				if verbose: print "Record exists, skipping (use flush=True to flush and re-enter these records)"
				continue
			except NCBIseq.DoesNotExist:
				pass

			# Get the protein entries from NCBI for the start/stop locus tags
			try:
				record = Entrez.read(Entrez.esearch(db="protein", retmax=1, term=local_bgc['start_locus_tag']))
				prot1  = SeqIO.read(Entrez.efetch(db="protein", id=record['IdList'][0], rettype="gb"), 'genbank')
				record = Entrez.read(Entrez.esearch(db="protein", retmax=1, term=local_bgc['end_locus_tag']))
				prot2  = SeqIO.read(Entrez.efetch(db="protein", id=record['IdList'][0], rettype="gb"), 'genbank')
			except IndexError:
				continue

			# Get the start position of the first locus tag and the end position of the second locus tag
			try:
				CDS1 = getCDS(prot1)
				IDs1 = getIDs(prot1)
				Pos1 = getPos(CDS1)
				CDS2 = getCDS(prot2)
				IDs2 = getIDs(prot2)
				Pos2 = getPos(CDS2)
			except IndexError:
				continue

			# Get the nucleotide genbank record for this region
			try:
				record = Entrez.read(Entrez.esearch(db="nucleotide", retmax=1, term=IDs2[0]))
				bgc    = SeqIO.read(Entrez.efetch(db="nucleotide", id=record['IdList'][0], rettype="gb", seq_start=Pos1[0], seq_stop=Pos2[1]), 'genbank')
			except: 
				continue


			# Get some information about this region and add to the database
			try:
				bgc.annotations.pop('references')
			except KeyError:
				pass
			bgc_record = NCBIseq.create(
										bgc            = local_bgc['bgc_id'],
										description    = bgc.description,
										dbxrefs        = bgc.dbxrefs,
										sequence_id    = bgc.id,
										sequence       = bgc.seq,
										annotations    = json.dumps(bgc.annotations),
									   )


			# Get local domains
			filters['bgc_id'] = key
			gq = BGCQuery(table="Gene", filters=filters, columns="gene_id AS id, gene.start AS start, gene.end AS end", verbose=False)
			genes = gq.execute(as_dict = False)
			offset = genes[0]['start']
			#print gq.query
			#pprint([(d['id'],d['start']-offset,d['end']-offset) for d in genes])
			# Define a lambda to find a fuzzy ncbi match based on start/end domain positions
			fuzzy_matches = lambda start,end: [m['id'] for m in genes if ((abs(m['start']-start-offset)+abs(m['end']-end-offset)) < 25)]
			total = count = 0
			# Iterating over each feature
			for feature in bgc.features:
				# If it is a coding sequence
				if feature.type == "CDS":
					total += 1
					# Determine if the entry exists
					start = feature.location.start.real
					end   = feature.location.end.real
					#print start, end
					try:
						first_fuzzy_match = fuzzy_matches(start,end)[0]
						count += 1
					except IndexError:
						first_fuzzy_match = None
					try:
						gene = NCBIcds.get((NCBIcds.bgc == local_bgc['bgc_id']) & (NCBIcds.start == start) & (NCBIcds.end == end))
					except NCBIcds.DoesNotExist:
						try:
							gene = NCBIcds.create(
														ncbiseq				= bgc_record,
														bgc					= local_bgc['bgc_id'],
														gene 				= first_fuzzy_match,
														start				= start,
														end					= end,
														strand				= feature.location.strand,
														sequence			= feature.qualifiers['translation'][0],
														locus 				= feature.qualifiers['locus_tag'][0],
														description			= feature.qualifiers['product'][0],
														notes				= json.dumps(feature.qualifiers),
													  )
						except KeyError:
							print "Missing a NCBI key for a gene, skipping"
							pass
			print "Associated %d of %d genes w/ NCBI CDS entries" % (count, total)
			#try:
			#	input("Press Enter to continue...")
			#except:
			#	pass


	def generate_arrows(self, use_augmented=True, verbose=None, show_domains=False, color_key='pfams'):
		''' Uses the arrower script to make bgc domain arrow images for each bgc in all clusters '''
		verbose = verbose if verbose != None else self.verbose
		found_some_supp = False
		# Remove the old output directory if exists
		try:
			shutil.rmtree(ARROW_OUTDIR)
		except OSError:
			pass
		# Remake it
		os.makedirs(ARROW_OUTDIR)
		# Debug
		# cluster_idx,bgcs = [(a,b) for a,b in enumerate(self.clusters['groups'])][1]
		for cluster_idx, bgcs in enumerate(self.clusters['groups']):
			# Get bgcs for this cluster and update to include domain information
			# bgc_id, attribs = [(a,b) for a,b in bgcs.iteritems() if b['bgc_id']==10545][0]
			for bgc_id, attribs in bgcs.iteritems():
				#if bgc_id==10545: print cluster_idx, bgc_id
				#if self.verbose: print "Making arrows for bgc_id: %d" %bgc_id
				# Copy DB filters
				filters = self.db.filters.copy()
				filters['bgc_id'] = bgc_id
				gq = BGCQuery(table="Gene", filters=filters, columns="gene.start AS start, gene.end AS end, *", verbose=False)
				genes = gq.execute(as_dict=False)
				# Add augmented NCBI data
				# gene = genes[0]
				if use_augmented:
					for gene in genes:
						try:
							supp = [s for s in self.get_augmented_ncbi([gene['gene_id']])][0]
						except:
							supp = {}
						gene.update(supp)
						if len(supp) != 0: 
							found_some_supp  = True
				# Add gene info to the bgc
				attribs['genes'] = genes

				# Add some domain specific filters
				#filters['domain_sequence'] = "IS NOT NULL"
				filters['pfam_id']         = "IS NOT 'Q'"
				# Peewee is responsible for the following stupidity...
				#  when returning results as a dict, keys are overwritten if table columns have same the name.
				#  Also, not sure why a standard SQL sum function doesn't work with peewee.
				dq = BGCQuery(table="Domain", \
							  filters=filters, \
							  columns=("gene.start AS gene_start, gene.end AS gene_end,\n  "
									    "domain.start AS domain_start, domain.end AS domain_end,\n  "
									    "gene.start+domain.start AS start, gene.start+domain.end AS end, *"), \
							  verbose=False)
				domains = dq.execute(as_dict=False)
				# Add domain information to the object
				attribs['domains'] = domains
			# Generate arrows for this cluster
			with open(os.path.join(ARROW_OUTDIR,"Cluster_%d_arrows.svg"%cluster_idx), 'w') as outfile:
				outfile.write(arrower(bgcs, show_domains=show_domains, color_key=color_key))
		# Warn the user if no augmentation data was found
		if verbose and use_augmented and not found_some_supp: 
			print "No augmented data (eg. DNA strand/description) found for any genes for one or more BGCs. Try running the augment() function"


	def get_augmented_ncbi(self, gene_ids):
		''' A simple helper function to see if NCBI augmented annotation information exists for a list of gene ids '''
		# Get ncbi cds matching these gene ids
		ncbi_cds = NCBIcds.select().where(NCBIcds.gene << gene_ids).dicts().execute()
		# For each one, yeild the results
		for ncbi_cd in ncbi_cds:
			yield(ncbi_cd)


	def _start_webserver(self):
		''' Starts a simple webserver to serve files used to visualize results (from the data directory; to visualize JSON cluster file with d3) '''
		try:
			Handler = SimpleHTTPServer.SimpleHTTPRequestHandler
			httpd = SocketServer.TCPServer(("", 8000), Handler)
			print "Started SimpleHTTPServer on port", 8000
			t1 = threading.Thread(target=httpd.serve_forever)
			t1.start()
		except SocketServer.socket.error:
			pass


	def get_bgc_edges(self, measure, cutoff, verbose=None):
		''' Return all pairwise edges between all bgcs (given by ids) for a measure with cutoff value '''
		verbose = verbose if verbose != None else self.verbose
		bgc_ids = self.db.get_bgcs(just_ids=True)
		edges = dict()
		query = '''
				SELECT id1, id2, value
				FROM similarity
				WHERE comparing="bgcs"
					AND measure=?
					AND value > ? 
					AND id1 IN(%s)
				'''
		query = query % ",".join([str(i) for i in bgc_ids])
		if verbose:
			print "\nRunning gene query:" + query + "\n[%s, %.2f]" % (measure, cutoff)
		results = Similarity.raw(query, measure, cutoff).execute()
		for result in results:
			edges[(result.id1, result.id2)] = result.value
		return edges
		# A cleaner alternative (keeping raw SQL in the Query class) is to use the peewee interface equivalent (not tested)
		#results = (Similarity
		#			.select(Similarity.id1, Similarity.id2, Similarity.value)
		#			.where( (Similarity.comparing=="bgcs") &
		#					(Similarity.measure == measure) & 
		#					(Similarity.value > cutoff) & 
		#					(Similarity.id1 << bgc_ids) )
		#			).execute()
		#for result in results:
		#	edges[(result.id1, result.id2)] = result.value
		#return edges


	def compare_bgcs(self, debug=False):
		''' Calls compare_two_bgcs for each pairwise combination of BGCS '''
		bgc_ids = self.db.get_bgcs(just_ids=True)
		num = len(bgc_ids)
		num_comarisons = math.factorial(num)/(math.factorial(num-2)*2)
		
		i = 0
		bgc_pairs = self._pairwise(bgc_ids, self_compare=False)
		#with db.transaction():
		for b1, b2 in bgc_pairs:
			i += 1
			if (i % 200 == 1): print "Currently comparing BGCs with ids %d and %d <%d of %d comparisons>" % (b1, b2, i, num_comarisons)
			self.compare_two_bgcs(b1, b2, debug)


	def compare_two_bgcs(self, bgc1_id, bgc2_id, debug=False, forced=False):
		''' Compares two BGC (bgc1_id, bgc2_id) using the Peter/Mohamed similarity measure.
			Results are returned and also stored in the DB for faster subsequent lookups.
			(also stores jaccard and dds measures in the DB)

			debug:   useful to show what the algorithm is doing (prints to screen for debug purposes)
			forced:  force the comparison (don't use value in the database)
		'''
		# Only needed for print statements
		if debug: binfo = self.db.get_bgcs(as_dict=True)

		# Sort the ids (id1 is always smaller value in database)
		# And get the domains for both BGCs
		bgc1_id, bgc2_id = sorted([bgc1_id, bgc2_id])

		# Return if the comparison exists in the database unless forced
		if not forced:
			try:
				Similarity.get((Similarity.id1==bgc1_id) & (Similarity.id2==bgc2_id) & (Similarity.comparing=="bgcs"))
				return
			except Similarity.DoesNotExist:
				pass

		# Copy the filters and add some domain specific ones 
		filters = self.db.filters.copy()
		filters['domain_sequence'] = "IS NOT NULL"
		filters['pfam_id']         = "IS NOT 'Q'"
		filters['bgc_id']          = bgc1_id
		bgc1_domains = BGCQuery(table="Domain", filters=filters, verbose=False).execute(as_dict = False)
		filters['bgc_id']          = bgc2_id
		bgc2_domains = BGCQuery(table="Domain", filters=filters, verbose=False).execute(as_dict = False)
		#bgc1_domains = self.get_bgc_domains(bgc1_id)
		#bgc2_domains = self.get_bgc_domains(bgc2_id)

		# Create a dictionary mapping pfam => list of domains (a more convient form)
		# And get pfam_ids for both BGCs
		get_pfam_ids  = lambda domains: [domain['pfam_id'] for domain in domains]
		bgc1_pfam_ids = get_pfam_ids(bgc1_domains)
		bgc2_pfam_ids = get_pfam_ids(bgc2_domains)
		#shared_pfams  = list(set(bgc1_pfam_ids) & set(bgc1_pfam_ids))

		# This function creates a list of domains for a pfam id and a cluster
		get_pfam_domain_ids = lambda pfam, domains: [domain['domain_id'] for domain in domains if domain['pfam_id']==pfam]
		# Set the overall dds and s variables
		dds, s = 0.0, 0.0
		# For each pfam found in either cluster
		for pfam in set(bgc1_pfam_ids + bgc2_pfam_ids):
			# Get the domains of each pfam occurrence for both clusters
			domains1 = get_pfam_domain_ids(pfam, bgc1_domains)
			domains2 = get_pfam_domain_ids(pfam, bgc2_domains)
			# Get their pairwise similarities (order doesn't matter, dms is a symmetric double hash)
			dms = self.get_bgc_similarity(domains1, domains2)

			# Whichever BGC has more domains for this pfam is seta (the list of domain_ids)
			# Flipped indicates if this reverse occurred (useful for printing and debugging)
			if len(domains2) > len(domains1):
				seta = domains2
				setb = domains1
			else:
				seta = domains1
				setb = domains2
			flipped = domains2==seta

			# General case of there being just one domain between both BGCs (really only have to test seta)
			# Just increase lin and domain duplication measure
			if len(seta) == 0 or len(setb) == 0:
				s += 1.
				dds += 1.
			# General case of 1 domain in each BGC for this pfam
			# Munkres not needed to find best pairs, just return the distance of the only possible pairing
			elif len(seta) == 1 and len(setb) == 1:
				# Get the similarity score for these two, or set distance to .9
				try:
					distance = 1 - (dms[seta[0]][setb[0]]/100.)
				except KeyError:
					distance = 1 - 0.1
				except TypeError:
					distance = 1 - 0.1
				s += 1.
				dds += distance
			# All other cases, use munkres to find the best possible pairing (most similar/least distant)
			else:
				if debug:
					statement = "\n" + "-"*100
					statement += "\n     %10s %15s %30s %30s"   % ("BGC ID","Cluster ID", "Start Locus Tag", "End Locus Tag")
					statement += "\nBGC1 %10d %15s %30s %30s"   % (bgc2_id, binfo[bgc2_id]['cluster_id'], binfo[bgc2_id]['start_locus_tag'], binfo[bgc2_id]['end_locus_tag'])
					statement += "\nBGC2 %10d %15s %30s %30s\n" % (bgc1_id, binfo[bgc1_id]['cluster_id'], binfo[bgc1_id]['start_locus_tag'], binfo[bgc1_id]['end_locus_tag'])
					statement += "-"*100
					statement += "\nPFAM          %s"           % (pfam)
					statement += "\nBGC1 domains |--%s--|"     % "--".join([str(i) for i in domains1])
					statement += "\nBGC2 domains |--%s--|\n"   % "--".join([str(i) for i in domains2])

				# Create an n by n pairwise matrix
				n = max(len(seta),len(setb))
				distances = np.zeros((n,n))
				doprint=False
				# Generate a distance matrix for all possible pairs
				for i in xrange(0,len(seta)):
					for j in xrange(0,len(setb)):
						distance = None
						if flipped and debug:
							statement += "\nindex%3d:\tBGC2 Domain: %d\nindex%3d:\tBGC1 Domain: %d (FLIPPED)" % (i, seta[i], j, setb[j])
						elif debug:
							statement += "\nindex%3d:\tBGC1 Domain: %d\nindex%3d:\tBGC2 Domain: %d" % (i, seta[i], j, setb[j])
						try:
							distance = 1 - (dms[seta[i]][setb[j]]/100.)
							if debug: statement += "\nDistance:\t%.2f\n" % (distance)
							doprint=True
						except KeyError:
							distance = 1 - 0.1
							if debug: statement += "\nDistance:\t%.2f (missing key! This shouldn't happen)\n" % (distance)
						except TypeError:
							distance = 1 - 0.1
							if debug: statement += "\nDistance:\t%.2f (no similarity match)\n" % (distance)
						distances[i,j] = distance
						# A more efficient loop does an internal loop of i to len(setb), in which case we can mirror over the diagonal
						#distances[j,i] = distance
				# Use munkres to find the best possible pairings from the distance matrix (least distant)
				m = Munkres()
				best_indices = m.compute(np.copy(distances))
				distance = sum([distances[bi] for bi in best_indices])
				# And compute the updated lin and domain duplication index
				s += max(len(seta),len(setb))
				dds += abs(len(seta)-len(setb)) - distance

				# print the debug statement if munkres was used (doprint) and debug is on
				if doprint and debug:
					print statement
					print "Distance Matrix:"
					print distances
					print "\nMunkres index pairings:"
					print best_indices
					for idx1, idx2 in best_indices:
						# Munkres pairs 'dummy' domains when the number of domains b/w BGCs is not =
						# These pairings have distance of 0 and don't increase dds or s and are ignored
						# But would be reported in this loop, so are passed if they are attempted to be looked-up
						try:
							print "Pairing Domains %d and %d" % (seta[idx1], setb[idx2])
						except: 
							pass

		# Calculate measures
		jaccard, similarity = (0, 0)
		jaccard = len(set(bgc1_pfam_ids) & set(bgc2_pfam_ids)) / \
				  float( len(set(bgc1_pfam_ids)) + len(set(bgc2_pfam_ids)) \
				  - len(set(bgc1_pfam_ids) & set(bgc2_pfam_ids)) )
		dds /= float(s)
		dds = math.exp(-dds)
		distance = 1 - (0.36 * jaccard) - (0.64 * dds)
		similarity = 1 - distance

		# Store measures in the database
		Similarity.create(
							id1			= bgc1_id,
							id2			= bgc2_id,
							measure		= "similarity",
							comparing	= "bgcs",
							value		= similarity
						 )
		Similarity.create(
							id1			= bgc1_id,
							id2			= bgc2_id,
							measure		= "dds",
							comparing	= "bgcs",
							value		= dds
						 )
		Similarity.create(
							id1			= bgc1_id,
							id2			= bgc2_id,
							measure		= "jaccard",
							comparing	= "bgcs",
							value		= jaccard
						 )
		return


	def precompute_similarities(self):
		''' Precomputes all BGC PFAM pairwise domain similarities by populating missing entries in the database
		    for faster access at some later point 
		'''
		bgc_ids = self.db.get_bgcs(just_ids=True)
		num = len(bgc_pairs)
		num_comarisons = math.factorial(num)/(math.factorial(num-2)*2)
		bgc_pairs = self._pairwise(bgc_ids, self_compare=False)
		i = 0
		for bgc1_id, bgc2_id in bgc_pairs:
			i += 1
			if (i % 200 == 1): print "Comparing BGCs with ids %d and %d <%d of %d>" % (bgc1_id, bgc2_id, i, num_comarisons)
			# Sort the ids (id1 is always smaller value in database)
			# And get the domains for both BGCs
			bgc1_id, bgc2_id = sorted([bgc1_id, bgc2_id])
			bgc1_domains = BGCQuery(table="Domain", filters=dict({'bgc_id':bgc1_id}, **self.db.filters)).execute()
			bgc2_domains = BGCQuery(table="Domain", filters=dict({'bgc_id':bgc2_id}, **self.db.filters)).execute()

			# Create a dictionary mapping pfam => list of domains (more convient form)
			# And get pfam_ids for both BGCs
			get_pfam_ids  = lambda domains: [domain['pfam_id'] for domain in domains]
			bgc1_pfam_ids = get_pfam_ids(bgc1_domains)
			bgc2_pfam_ids = get_pfam_ids(bgc2_domains)
			#shared_pfams  = list(set(bgc1_pfams) & set(bgc1_pfams))

			# This function creates a list of domains for a pfam id and a cluster
			get_pfam_domain_ids = lambda pfam, domains: [domain['domain_id'] for domain in domains if domain['pfam_id']==pfam]
			# Set the overall dds and s variables
			dds, s = 0.0, 0.0
			# For each pfam found in either cluster
			for pfam in set(bgc1_pfam_ids + bgc2_pfam_ids):
				domains1 = get_pfam_domain_ids(bgc1, pfam)
				domains2 = get_pfam_domain_ids(bgc2, pfam)
				# Get their pairwise similarities (order doesn't matter, dms is a symmetric double hash)
				self.get_bgc_similarity(domains1, domains2)
		return


	def get_bgc_similarity(self, bgc1_id, bgc2_id):
		''' Returns the similarities for a given list of domain_ids (also storing in the db if not present)'''
		similarities = dict()
		bgc1_domain_ids = BGCQuery(table="Domain", filters=dict({'bgc_id':bgc1_id}, **self.db.filters), verbose=False).execute(just_ids=True)
		bgc2_domain_ids = BGCQuery(table="Domain", filters=dict({'bgc_id':bgc2_id}, **self.db.filters), verbose=False).execute(just_ids=True)
		# Get all the entries we can from the database or use usearch to get the values
		for (d1, d2) in self._pairwise2(bgc1_domain_ids, bgc2_domain_ids):
			identity = None
			try:
				# Using peewee calls
				identity = Similarity.get((Similarity.id1==d1) & (Similarity.id2==d2) & (Similarity.comparing=="domains")).value
				identity = float(identity) if identity != None else None
			except Similarity.DoesNotExist:
				# Get all the missing ids pairwise identities via usearch
				usearch_sims = self._usearch([d1,d2])
				# populate the structure by
				# First seeing if usearch found an identity (or use None)
				try:
					identity = usearch_sims[d1][d2]
				except KeyError:
					identity = None
				# and in the database
				Similarity.create(
						id1			= d1,
						id2			= d2,
						measure		= "usearch_identity",
						comparing	= "domains",
						value		= identity
					   )
			# Set the identity in the final result structure (symmetric)
			try:
				similarities[d1][d2] = identity
			except KeyError:
				similarities[d1] = { d2: identity }
			try:
				similarities[d2][d1] = identity
			except KeyError:
				similarities[d2] = { d1: identity }
		return similarities


	def usearch_domains(self): self._usearch(BGCQuery(table="Domain", filters=dict({'bgc_id':self.db.get_bgcs(just_ids=True)}, **self.db.filters)).execute(just_ids=True))
	def _usearch(self, domain_ids, delete=True, verbose=None):
		'''Computes pairwise similarity calculation for given domains using usearch '''
		verbose = verbose if verbose != None else self.verbose
		domain_ids = self._force_int_list(domain_ids)
		# Generate some tempfiles to mediate usearch execution, and 
		# get fasta information using the make_fasta function
		results = None
		# use tempfile.NamedTemporaryFile(delete=False) to keep files
		with tempfile.NamedTemporaryFile(delete=delete) as intemp:
			with tempfile.NamedTemporaryFile(delete=delete) as outtemp:
				if delete==False and verbose:
					print intemp.name
					print outtemp.name
				self.make_fasta(domain_ids, location=intemp)
				intemp.flush()
				intemp.seek(0)
				command = "usearch " + \
						  "-acceptall " + \
						  "-allpairs_global \"%s\" " % intemp.name + \
						  "-userout \"%s\" " % outtemp.name + \
						  "-userfields query+target+id+evalue+alnlen+mism+opens+qlo+qhi+tlo+thi+bits"
				if verbose:
					print "Running usearch..."
				else:
					command += " > /dev/null 2>&1"
					#print "Usearch command: %s" % "".join(command)
				subprocess.call(command, shell=True)
				outtemp.flush()
				outtemp.seek(0)
				results = outtemp.read()
		identities = dict()
		for line in results.split('\n'):
			try:
				domain1, domain2, identity = line.split('\t')[:3]
				domain1  = int(domain1)
				domain2  = int(domain2)
				identity = float(identity)
			except ValueError:
				continue
			try:
				identities[domain1][domain2] = identity
			except KeyError:
				identities[domain1] = { domain2: identity }
			try:
				identities[domain2][domain1] = identity
			except KeyError:
				identities[domain2] = { domain1: identity }
		return identities


	def _mcl(self, bgcs, edges, I=2.0, delete=True, verbose=None):
		''' Uses MCL to automatically cluster bgcs 

			bgcs:    A dictionary of BGC objects (BGCQuery.get_bgcs(as_dict=True))
			edges:   A dictionary of BGC edges (BGCs.get_bgc_edges(measure, cutoff))
			delete:  delete input/output files used by mcl when done
		'''
		verbose = verbose if verbose != None else self.verbose
		# Generate some tempfiles to mediate mcl execution, and 
		results = None
		# use tempfile.NamedTemporaryFile(delete=False) to keep files
		with tempfile.NamedTemporaryFile(delete=delete) as intemp:
			with tempfile.NamedTemporaryFile(delete=delete) as outtemp:
				if delete==False and verbose:
					print intemp.name
					print outtemp.name
				self.make_cytoscape(bgcs, edges, location=intemp, annotation=False)
				intemp.flush()
				intemp.seek(0)
				command = "mcl " + \
						  "\"%s\" " % intemp.name + \
						  "--abc" + \
						  " -I %s" % str(I) + \
						  " -p 1.0E-15" + \
						  " -o \"%s\" " % outtemp.name
				if verbose:
					print "Running mcl..."
				else:
					command += " > /dev/null 2>&1"
				subprocess.call(command, shell=True)
				outtemp.flush()
				outtemp.seek(0)
				results = outtemp.read()
		out_edges = dict()
		for cluster in results.split('\n'):
			nodes = cluster.split()
			nodes = [int(node) for node in nodes]
			nodes.sort()
			for pair in self._pairwise(nodes, self_compare=False):
				if pair in edges:
					out_edges[pair] = edges[pair]
		return out_edges


	def make_fasta(self, ids, location=None, verbose=None):
		''' Generates a fasta file for the domains generated from a given ID list (ids)
			and corresponding query that generates sequences for each ID
			An optional location parameter specifies a file handle to 
			write the results ( or returns by default)
		'''
		verbose = verbose if verbose != None else self.verbose
		# Copy the current db filter and add ones to request unique, present sequences 
		filters = self.db.filters.copy()
		filters['domain_sequence'] = "IS NOT NULL"
		filters['pfam_id']         = "IS NOT 'Q'"
		
		# Requery the database for additional information to put in the fasta file
		if verbose:
			print "Writing to file" if type(location) in ['file', 'instance'] else ""
		fasta_info_headers = ["domain_id", "pfam_id", "img_gene_id", "cluster_id", "img_genome_id", \
							  "_genus", "_species", "kind","cluster_finder_prob"]
		results = None
		for id_ in ids:
			# Add a filter for this bgc's id and query for the sequences
			filters['bgc_id'] = id_
			seqs = BGCQuery(table="Domain", filters=filters, verbose=False).execute(as_dict=False)
			# For each pfam_id append results to the output file or generate a new one
			if seqs != None:
				for seq in seqs:
					# Prepare header line and sequence
					seq['domain_id'] = str(seq['domain_id']) + " "
					header = ">" + "|".join([str(seq[head]) for head in fasta_info_headers]) + "\n"
					seq    = seq['sequence']+"\n"
					# If no location is given, add the results list
					if location == None:
						results.append({'pfam_id': {'header': header, 'sequence': seq}})
					# Otherwise determine if the file handle is valid and writable, and write to it
					#elif type(location) == file:
					else:
						location.write(header)
						location.write(seq)
		return results


	def make_cytoscape(self, bgcs, edges, location=None, annotation=True, sep="\t", verbose=None):
		''' Generates a tab-delimited annotated output table for the given BGCs (given as list of ids)
			This file is compatible with both cytoscape and mcl
		'''
		verbose = verbose if verbose != None else self.verbose
		if verbose:
			print "Writing to file" if type(location) in ['file', 'instance'] else ""
		for pair, weight in edges.iteritems():
			bgc1 = bgcs[pair[0]]
			bgc2 = bgcs[pair[1]]
			# Get similar key/value pairs
			#same_annots = list(set([(key,val) for key,val in bgc1.items()]) & set([(key,val) for key,val in bgc2.items()]))
			line = "{}{sep}{}{sep}{:.2}".format(int(pair[0]), int(pair[1]), float(weight), sep=sep)
			if annotation:
				line += "{sep}{}{sep}{}{sep}{}{sep}{}".format(bgc1['cluster_id'], bgc2['cluster_id'], bgc1['_species'], bgc2['_species'], sep=sep)
			location.write(line+"\n")


	def write_clusters(self, outfile=None, sep="\t", header=True, verbose=None):
		''' Generates a tab-delimited (default, or use sep="delimiter") annotated output file for the clustered BGCs
		'''
		verbose = verbose if verbose != None else self.verbose
		if outfile == None:
			print "Please specifiy an output location with outfile=\"path/to/file\""
			return

		with open(outfile,'w') as out:
			if header:
				out.write("{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}".format("Cluster Group", "Domain",          "Phylum",        "Class",        "Order",              "Family",        "Genus",      "Species", s=sep))
				out.write("{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}".format("BGC Kind",      "BGC Start Locus", "BGC End Locus", "Clusterer ID", "BGC ID",             "Organism ID",   "Taxon ID",   "IMG Genome ID", s=sep))
				out.write("{}{s}{}{s}{}{s}{}{s}{}{s}{}\n".format(           "Genome Size",   "Gene Count",      "Genome Name",   "Genome Key",   "Sequencing Center",  "Status", s=sep))
			for cluster in self.clusters['nodes']:
				cluster_no = cluster['cluster']
				attribs    = cluster['attribs']
				out.write("{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}".format(cluster_no, attribs['_domain'], attribs['_phylum'], attribs['_class'], attribs['_order'], attribs['_family'], attribs['_genus'], attribs['_species'], s=sep))
				out.write("{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}{}{s}".format(attribs['kind'], attribs['start_locus_tag'], attribs['end_locus_tag'], attribs['cluster_id'], attribs['bgc_id'], attribs['organism'], attribs['taxon_id'], attribs['img_genome_id'], s=sep))
				out.write("{}{s}{}{s}{}{s}{}{s}{}{s}{}\n".format(           attribs['genome_size'], attribs['gene_count'], attribs['genome_name'], attribs['genome_name_key'], attribs['sequencing_center'], attribs['status'], s=sep))
		print "Done"

	# Some handy little functions
	def _pairwise(self, _list, self_compare=True):
		for i in xrange(len(_list)):
			for j in xrange(i,len(_list)):
				if i==j and self_compare==False:
					continue
				yield (_list[i],_list[j])

	def _pairwise2(self, _lista, _listb):
		for i in xrange(len(_lista)):
			for j in xrange(len(_listb)):
				yield (_lista[i],_listb[j])

	def _flatten(self, lol):
		return list(set([val for subl in lol for val in subl]))

	def _force_int_list(self, vals):
		if type(vals) in  [int, float, str]:
			vals = [vals]
		if type(vals) != list:
			print "Please specify an integer or integer list"
			return 
		try:
			vals = [int(val) for val in vals]
		except ValueError:
			print "Please specify an integer or integer list"
			return
		return vals

	# The following two weren't tested.... use with caution
	def _flush(self):
		''' Perform the directory flush and recreate structure '''
		response = None
		while not (response == 'n' or response == 'y'):
			response = raw_input("This will purge %s. Are you sure? (y/N) " % FASTA_OUTDIR).lower()
		if response == 'n':
			return
		if os.path.exists(FASTA_OUTDIR):
			shutil.rmtree(FASTA_OUTDIR)
		self._create_dir_tree()

	def _create_dir_tree(self):
		if not os.path.exists(FASTA_OUTDIR):
			os.makedirs(FASTA_OUTDIR)
			os.makedirs(os.path.join(FASTA_OUTDIR, "fastas"))
			os.makedirs(os.path.join(FASTA_OUTDIR, "usearch"))
			os.makedirs(ARROW_OUTDIR)
