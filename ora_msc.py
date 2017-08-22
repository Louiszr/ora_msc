# Access to KEGG API
from bioservices.kegg import KEGG

import random
import numpy as np
import readline

# Hypergeometric tests
from scipy.stats import hypergeom

# FDR (Benjamini-Hochberg) in R
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector


def read_annotations(annotation_file):
    annotation_fh = open(annotation_file, 'r')
    annotations = annotation_fh.readlines()
    annotations = list(map(str.rstrip, annotations))
    annotation_fh.close()
    return annotations

def filter_zscores(ko_index, modzscore_file, annotations, zscore_threshold):
    '''
    Due to RAM restrictions, must generate input list on the fly
    '''
    metabolite_hits = []
    modz_fh = open(modzscore_file, 'r')
    modz_lines = modz_fh.readlines()
    modz_fh.close()
    for index in range(0, len(modz_lines)):
        line = modz_lines[index]
        zscores = line.rstrip().split()
        exp_zscore = float(zscores[ko_index]) # The zscore of that knockout
        if abs(exp_zscore) >= zscore_threshold:
            kegg_id = annotations[index]
            if kegg_id.startswith('C'): # Verifying KEGG ids
                metabolite_hits.append(kegg_id)
    metabolite_hits = set(metabolite_hits)
    return metabolite_hits

def build_metabo_input(ko_index, 
                       pos_annotation_file, pos_modzscore_file, 
                       neg_annotation_file, neg_modzscore_file, 
                       zscore_threshold):
    pos_annotations = read_annotations(pos_annotation_file)
    neg_annotations = read_annotations(neg_annotation_file)

    pos_hits = filter_zscores(ko_index, pos_modzscore_file, pos_annotations, zscore_threshold)
    neg_hits = filter_zscores(ko_index, neg_modzscore_file, neg_annotations, zscore_threshold)
    
    metabolite_hits = pos_hits | neg_hits
    return metabolite_hits

# Get ALL the compounds within a species's pathway db
def get_all_compounds(species):
    # Initiate KEGG instance
    kegg_inst = KEGG()
    kegg_inst.organism = species
    # Get all compounds
    all_compounds = set()
    species_pathways = kegg_inst.pathwayIds
    for pathway in species_pathways:
        parsed_output = kegg_inst.parse(kegg_inst.get(pathway)) # parsed_ouput has lots of information about the pathway
        try:
            compounds = set(parsed_output['COMPOUND'].keys())
            all_compounds = all_compounds | compounds
        except KeyError: # Some pathways do not have defined compounds
            pass
    return all_compounds

# ORA (only for E.coli)
def ora(in_metabolites, pathway_id, bg_metabolites, pathway_2_compounds, least_num_metabolites=1):
    '''
    Specifying a KEGG instance is way faster than creating one on the fly
    Need to specify organism for that instance though
    '''
    # Get compounds
    try:
        compounds = pathway_2_compounds[pathway_id]
    except KeyError:
        return 'No compounds (DB)'
    # Background filtering
    test_pathway_compounds = compounds & bg_metabolites
    if len(test_pathway_compounds) == 0:
        return 'No compounds (BG)'
    # Hypergeometric test
    test_in_metabolites = in_metabolites & test_pathway_compounds
    if len(test_in_metabolites) < least_num_metabolites:
        return 'Low metabolites'
    hyperg_test = hypergeom(len(bg_metabolites), len(test_pathway_compounds), len(in_metabolites & bg_metabolites))
    ora_raw_pval = 1 - hyperg_test.cdf(len(test_in_metabolites)) + hyperg_test.pmf(len(test_in_metabolites))
    return ora_raw_pval, pathway_id, len(test_pathway_compounds)

# Run ORAs for all pathways in a single knockout
def oras_ko(ko_number, testing_pathways, background_metabolites, pathway_2_compounds,  
            pos_annotation_file, pos_modzscore_file, 
            neg_annotation_file, neg_modzscore_file, 
            zscore_threshold, 
            multiple_testing_correction, minus_log_trans, max_mutation, mutation_pool):
    

	'''
	Only usable for the Zamboni dataset
	Could replace it with oras_allpaths in the future
	'''

    ko_metabolites = build_metabo_input(ko_number, 
                       pos_annotation_file, pos_modzscore_file, 
                       neg_annotation_file, neg_modzscore_file, 
                       zscore_threshold)
    # Random Mutation
    if max_mutation:
        tmp_metabolites = ko_metabolites # Copy ko_metabolites
        (background_metabolites, old_met, new_met) = misidentify(background_metabolites, 
                                                                 max_mutation, mutation_pool, True)
        met_translate = dict(zip(old_met, new_met)) # Make a dictionary for translation
        ko_metabolites = set() # Empty ko_metabolites
        for met in tmp_metabolites: # For each member in tmp, translate them and add them back to ko_metabolites
            ko_metabolites.add(met_translate.get(met, met))
    # Generating raw p values
    pval = []
    pathwayid = []
    pathwaysize = []
    for pathway in testing_pathways:
        ora_res = ora(ko_metabolites, pathway, background_metabolites, pathway_2_compounds)
        if len(ora_res) == 3: # if both ora_raw_pval and pathway_id are returned
            pval.append(ora_res[0])
            pathwayid.append(ora_res[1])
            pathwaysize.append(ora_res[2])
    # Multiple testing correction
    if multiple_testing_correction:
        pval = list(importr('stats').p_adjust(FloatVector(pval), method = 'BH'))
    # -log transformation
    if minus_log_trans:
        pval = list(map(np.log10, pval))
        pval = list(map(np.negative, pval))
    return pval, pathwayid, pathwaysize


def oras_allpaths(in_metabolites, testing_pathways, background_metabolites, pathway_2_compounds,  
            multiple_testing_correction, minus_log_trans, max_mutation, mutation_pool):
    
    '''
    Run ORAs for all pathways
    in_metabolites must be type set()
    '''

    # Random Mutation
    if max_mutation:
        tmp_metabolites = in_metabolites # Copy in_metabolites
        (background_metabolites, old_met, new_met) = misidentify(background_metabolites, 
                                                                 max_mutation, mutation_pool, True)
        met_translate = dict(zip(old_met, new_met)) # Make a dictionary for translation
        in_metabolites = set() # Empty in_metabolites
        for met in tmp_metabolites: # For each member in tmp, translate them and add them back to in_metabolites
            in_metabolites.add(met_translate.get(met, met))
    # Generating raw p values
    pval = []
    pathwayid = []
    pathwaysize = []
    for pathway in testing_pathways:
        ora_res = ora(in_metabolites, pathway, background_metabolites, pathway_2_compounds)
        if len(ora_res) == 3: # if both ora_raw_pval and pathway_id are returned
            pval.append(ora_res[0])
            pathwayid.append(ora_res[1])
            pathwaysize.append(ora_res[2])
    # Multiple testing correction
    if multiple_testing_correction:
        pval = list(importr('stats').p_adjust(FloatVector(pval), method = 'BH'))
    # -log transformation
    if minus_log_trans:
        pval = list(map(np.log10, pval))
        pval = list(map(np.negative, pval))
    return pval, pathwayid, pathwaysize

def misidentify(in_metabolites, times, pool_metabolites, return_changes = False):
    '''
    in_metabolites: set of in metabolites
    times: number of misidentification
    pool_metabolites: set pool of metabolites which the new identification can be
    '''
    in_metabolites = list(in_metabolites)
    pool_metabolites = list(pool_metabolites)
    # Adjusting times so times <= length of in metabolties
    if times > len(in_metabolites):
        times = len(in_metabolites)
    elif times > (len(pool_metabolites) - len(in_metabolites)):
        times = (len(pool_metabolites) - len(in_metabolites))
    if times == 0:
        return in_metabolites
    # Generating new identifications
    new_ident = []
    for i in range(0, times):
        random.shuffle(pool_metabolites)
        new_metabolite = pool_metabolites[0]
        # If new metab is already in the input list, discard and generate new ones
        while new_metabolite in in_metabolites:
            random.shuffle(pool_metabolites)
            new_metabolite = pool_metabolites[0]
        new_ident.append(new_metabolite)
    # Swapping existed metabolites for new metabolites
    old_ident = []
    random.shuffle(in_metabolites)
    for i in range(0, times):
        old_ident.append(in_metabolites.pop())
    in_metabolites += new_ident
    out_metabolites = set(in_metabolites)
    # Added an option to return identity changes
    if return_changes:
        return out_metabolites, old_ident, new_ident
    else:
        return out_metabolites

def loadTsv(filename):
    fh = open(filename, 'r')
    outset = set()
    for line in fh:
        if line.startswith('C'):
            outset = outset | set([line.rstrip()])
    fh.close()
    return outset

def dup_argsort(in_val):
    u, v = np.unique(in_val, return_inverse=True)
    out_ind = (np.cumsum(np.concatenate(([0], np.bincount(v)))))[v]
    return out_ind

# Generating null models
def make_null_model(original_model, metabolite_pool, overlap=False):
    '''
    original_model -- the dictionary with pathways as keys and metabolites as values
    metabolite_pool -- the collection of metabolites which the null model can include
    overlap -- True: overlap between pathways will not be randomised; False: overlap will be randomised
    '''
    if overlap:
        # If keeping the metabolite overlap, we only shuffle the labels of the metabolites
        
        # Getting the original metabolites
        original_metabolites = set()
        for value in original_model.values():
            original_metabolites = original_metabolites | value
        original_metabolites = list(original_metabolites)
        
        random.shuffle(metabolite_pool)
        metabolite_translate = dict(zip(original_metabolites, metabolite_pool[:len(original_metabolites)]))
    new_model = dict()
    for pathway in original_model:
        pathwaysize = len(original_model[pathway])
        if overlap:
            new_model[pathway] = set()
            for metabolite in original_model[pathway]:
                new_model[pathway].add(metabolite_translate[metabolite])
        else:
            for i in range(0, pathwaysize):
                random.shuffle(metabolite_pool)
                new_model[pathway] = set(metabolite_pool[:pathwaysize])
    return new_model

# Save null models and load null models
def save_null_model(model, filename):
    with open(filename, 'w') as fh:
        for pathway in model:
            metabolites = list(model[pathway])
            fh.write(pathway+'\t')
            fh.write(','.join(metabolites))
            fh.write('\n')
            
def load_null_model(filename):
    model = dict()
    with open(filename, 'r') as fh:
        for line in fh.readlines():
            fields = line.rstrip().split('\t')
            pathway = fields[0]
            metabolites = set(fields[1].split(','))
            model[pathway] = metabolites
    return model