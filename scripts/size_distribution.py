# Access to KEGG API
from bioservices.kegg import KEGG

import ora_msc

import matplotlib.pyplot as plt

# Define the path of metabolomics data
DATA_PATH = './data/'

# Stating the annotation files & modzscore files
pos_annot = DATA_PATH + 'annotation_pos.txt'
pos_mod = DATA_PATH + 'modzscore_pos_annotated.tsv'
neg_annot = DATA_PATH + 'annotation_neg.txt'
neg_mod = DATA_PATH + 'modzscore_neg_annotated.tsv'
# Initialise KEGG instance
kegg_instance = KEGG()
kegg_instance.organism = "eco"
# Initialise both backgrounds
test_compounds = ora_msc.get_all_compounds('eco')
zamboni_bg = ora_msc.loadTsv(DATA_PATH + 'annotation_all.txt')
# Remove metabolites detected in Zamboni but not in any E.coli pathway
zamboni_bg = zamboni_bg & test_compounds
# build {pathway: compounds} dictionary for E.coli
ecoli_pathways = kegg_instance.pathwayIds
pathway_2_compounds = dict()
for pathway in ecoli_pathways:
    parsed_output = kegg_instance.parse(kegg_instance.get(pathway)) # parsed_ouput has lots of information about the pathway
    try:
        compounds = set(parsed_output['COMPOUND'].keys())
        pathway_2_compounds[pathway] = compounds
    except KeyError: # Some pathways do not have defined compounds
        pass

# Translate KO number to gene name
sample_id_all = DATA_PATH + 'sample_id_modzscore.tsv'
all_knockouts = []# End product
fh_sample_id_all = open(sample_id_all, 'r')
for knockout in fh_sample_id_all:
    all_knockouts.append(knockout.rstrip())
fh_sample_id_all.close()


size_dist = []
for pathway in pathway_2_compounds:
    #if len(pathway_2_compounds[pathway]) == 1:
    #    print(pathway)
    size_dist.append(len(pathway_2_compounds[pathway]))

zamboni_size_dist = []
for pathway in pathway_2_compounds:
    compounds = pathway_2_compounds[pathway]
    cmpd_count = 0
    for compound in compounds:
        if compound in zamboni_bg:
            cmpd_count += 1
    zamboni_size_dist.append(cmpd_count)

plt.subplot(211)
plt.hist(zamboni_size_dist, bins=range(0, 145, 5))
plt.ylim(0, 40)
plt.xlabel('Pathway size')
plt.ylabel('Number of pathways')
plt.title('Pathway size distribution (Zamboni background)')

plt.subplot(212)
plt.hist(size_dist, bins=range(0, 145, 5))
plt.ylim(0, 40)
plt.xlabel('Pathway size')
plt.ylabel('Number of pathways')
plt.title('Pathway size distribution (all compounds)')

plt.tight_layout()
plt.show()