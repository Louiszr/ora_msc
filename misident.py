# Access to KEGG API
from bioservices.kegg import KEGG

import ora_msc

import random

import numpy as np

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


# Generate random knockout numbers
random_knockouts = np.random.randint(3717, size=50)

# Random metabolite mutation

# Change misidentification & number of simulations rate here
MISIDENT_RATE = 0.1
NUMBER_SIMULATION = 20
# 

# Change results file name here
filename = './mrate' + str(MISIDENT_RATE * 100) + '.tsv'

fh = open(filename, 'w')
for ko_number in random_knockouts:
    fp = 0
    fn = 0
    ora_results = []
    for i in range(0, NUMBER_SIMULATION + 1):
        ora_results.append([])
    (pvals, pathwayids, pathsizes) = ora_msc.oras_ko(ko_number, ecoli_pathways, zamboni_bg, pathway_2_compounds, 
                                  pos_annot, pos_mod, neg_annot, neg_mod, 2, 
                                  True, False, 0, [])
    for ind in range(0, len(pvals)):
        if pvals[ind] < 0.05:
            ora_results[0].append(pathwayids[ind])

    for k in range(0, NUMBER_SIMULATION): # Number of mutations per ko
        (pvals_mut, pathwayids_mut, pathsizes_mut) = ora_msc.oras_ko(ko_number, ecoli_pathways, zamboni_bg, pathway_2_compounds, 
                                              pos_annot, pos_mod, neg_annot, neg_mod, 2, 
                                              True, False, int(MISIDENT_RATE * len(zamboni_bg)), test_compounds)
        for ind in range(0, len(pvals_mut)):
            if pvals_mut[ind] < 0.05:
                ora_results[k+1].append(pathwayids_mut[ind])
    # write ora_result to a file
    for i in range(1, len(ora_results)):
        result = ora_results[i]
        fp += len(set(result) - set(ora_results[0]))
        fn += len(set(ora_results[0]) - set(result))
    fh.write('\t'.join([str(len(ora_results[0])), str(fp), str(fn), str(ko_number)]))
    fh.write('\n')
fh.close()