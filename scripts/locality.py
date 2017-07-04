# Directly from Notbook, modifications are needed before this script can work independently...

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


# Locality Analysis

# Find all enzymatic genes
eco_enzymes = []
for i in range(0, len(all_knockouts)):
    try:
        if 'Enzymes' in kegg_instance.parse(kegg_instance.get(kegg_instance.find('eco', all_knockouts[i])[:9])).get('BRITE', []):
            eco_enzymes.append(i)
    except AttributeError:
        pass

outfile = './locality/allKO.tsv'
fh = open(outfile, 'w')
for ko_number in eco_enzymes:
    ko_gene = all_knockouts[ko_number]
    gene_id = kegg_instance.find("eco", ko_gene)[:9]
    try:
        gene_pathways = kegg_instance.parse(kegg_instance.get(gene_id))['PATHWAY'].keys()
    except KeyError:
        gene_pathways = []
    except TypeError:
        gene_pathways = []
    if len(gene_pathways) > 0:
        fh.write(str(ko_number) + '\t' + ko_gene + '\t')
        fh.write(' '.join(list(gene_pathways)))
        fh.write('\n')
fh.close()

# Get KO number
# Load rankings for that KO
# Get significant/not related/average rank
kegg_results = './locality/allKO.tsv'
fh = open(kegg_results, 'r')
result_lines = fh.readlines()
fh.close()
out_fh = open('./locality/outresult_sig.tsv', 'w')
out_fh.write('KO_number\tGene\tSigpathways\tNot_KEGG\tSigKEGG\tRank\tPath_count\n')
for line in result_lines:
    fields = line.rstrip().split('\t')
    ko_number = fields[0]
    ko_name = fields[1]
    ko_kegg_pathways = fields[2].split()
    
    ko_ora_results = './allresult/KO' + ko_number + '.tsv'
    fh = open(ko_ora_results, 'r')
    ora_lines = fh.readlines()
    fh.close()
    ora_pvals = []
    ora_pathways = []
    ora_sigpathways = []
    for oraline in ora_lines:
        ora_pathway_result = oraline.rstrip().split('\t')
        ora_pvals.append(ora_pathway_result[1])
        ora_pathways.append(ora_pathway_result[0])
        if float(ora_pathway_result[1]) < 0.05:
            ora_sigpathways.append(ora_pathway_result[0])
    pathway_rank = dict(zip(ora_pathways, list(dup_argsort(ora_pvals))))
    # Sigpathway
    not_related_pathways = len(set(ora_sigpathways) - set(ko_kegg_pathways))
    missing_pathways = len(set(ko_kegg_pathways) - set(ora_sigpathways))
    # average rank
    ranksum = 0
    path_count = 0
    for path in ko_kegg_pathways:
        try:
            rank = pathway_rank[path]
            ranksum += rank
            path_count += 1
        except KeyError:
            pass
    if path_count != 0:
        rankavg = ranksum/path_count
    else:
        rankavg = 'N'
    output_str = '\t'.join([ko_number, ko_name, 
                  str(len(ora_sigpathways)), str(not_related_pathways), str(len(ora_sigpathways) - not_related_pathways), str(rankavg), str(path_count)])
    if len(ora_sigpathways) > 0:
        out_fh.write(output_str+'\n')
out_fh.close()