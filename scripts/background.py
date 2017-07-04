# Access to KEGG API
from bioservices.kegg import KEGG

import ora_msc

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


# Background Analysis
for ko_number in range(2406, 3717):
    nobg_pval, nobg_pathway_id, nobg_sizes = ora_msc.oras_ko(ko_number, ecoli_pathways, test_compounds, pathway_2_compounds, 
                        pos_annot, pos_mod, neg_annot, neg_mod, 2, True, True, 0, [])
    zamboni_pval, zamboni_pathway_id, zamboni_sizes = ora_msc.oras_ko(ko_number, ecoli_pathways, zamboni_bg, pathway_2_compounds, 
                           pos_annot, pos_mod, neg_annot, neg_mod, 2, True, True, 0, [])
    # Define where to save results here
    result_file = './Backgrounds/KO' + str(ko_number) + '.tsv'
    fh = open(result_file, 'w')
    for i in range(0, len(nobg_pathway_id)):
        fh.write('{}\t{}\t{}\t{}\t{}\n'.format(nobg_pathway_id[i][5:], nobg_pval[i], nobg_sizes[i], zamboni_pval[i], zamboni_sizes[i]))
    fh.close()