from BGCExplorer import *

## NRPS = BGCs(filters={'Kind':"NRPS"})  # Put your filters here (see filters.txt)
## NRPS   = BGCs(filters={'Kind':"NRPS", 'Species':["%acnes%","%epi%"]})
BACSAC = BGCs(filters={'Kind':"saccharide", 'Genus':"Bacteroides"})
#NRPS.compare_bgcs()                   # This will compare all filtered BGCs and store values in the DB
#NRPS.augment()                        # Recommend but not required. Will get supplemental data on BGCs from NCBI
#NRPS.cluster(cutoff=.7)                        # Perform the clustering based on filtering criteria and similarity comparisons
#NRPS.visualize_clusters()             # Show the clustering results as a bubble plot in a web browser
#	-- or --
#NRPS.visualize_network()              # Show the clustering results as a network in a web browser
BACSAC.cluster(cutoff=.6)
BACSAC.visualize_network()
BACSAC.write_clusters()