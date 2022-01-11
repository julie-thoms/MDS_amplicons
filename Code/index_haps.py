#Module for assigning haplotypes to single cell amplicon data

#Authors:
#	Julie Thoms
#	Annatina Schnegg-Kaufmann
#	Fabio Zanini

#Version 1
#12th August 2021

#Version 2
#10th January 2022

#Import dependent modules
import os
import sys
import math
import numpy as np
import numpy.linalg
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import ternary
from scipy.stats import multinomial

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

#Define functions

def about():

	'''
	The functions in this module are used to assign haplotypes to single cells genotyped with PCR amplicons.

    The callable functions are;
        data_retrieval2(sourcefile, metadata, pt_id)
	    call_haps(data, pt_id, haps, reads,  cutoff)
		plot_hap_dist_sort_type(data, pt_id, haps, reads, cutoff, save = False)
        plot_index_heatmap_3(data, title, haps, reads, cutoff, save = False)
        calc_scVAF_mod_resc(data, pt_id, reads, draw_plot = False)

	
	'''
	print('This module contains functions to calculate allele frequencies and assign haplotypes. For more info type index_haps.about?')



def data_retrieval2(sourcefile, metadata, pt_id):

    '''
    This function reads data for a single patient from a master spreadsheet with amplicon data for all plates/patients.
    This works with clean data from the a_new_data_cleanup notebook.
    Input is the sourcefile with the readcounts (allele_counts_anon), and a metadata file which contains cell type for each plate (Amplicon_metadata_fixed_anon.xlsx).
    The function returns a dataframe containing just the data for the specified patient, ready to merge with index data and then plot.
    '''

    df = pd.read_csv(sourcefile, index_col = [0, 1, 2, 3], sep = '\t')
    df.columns = ['Reads']
    df['Plate'] = df.index.get_level_values(0)  #These lines send indexes to columns
    df['Well'] = df.index.get_level_values(1)
    df['Amplicon'] = df.index.get_level_values(2)
    df['Genotype'] = df.index.get_level_values(3)
    df[['Patient', 'one', 'two']] = df['Amplicon'].str.split('_', expand = True)
    df = df.drop(columns = ['one', 'two'])

    #Import information about plate cell type and patient
    key = pd.read_excel(metadata, sheet_name = 'PlateID')
    key = key.drop(['Cell Origin', 'Plate Nr', 'Plate Name','Nr of cells'], axis=1)
    key.rename(columns = {'Comments2':'Plate'}, inplace = True)
    key.rename(columns = {'Cell-group':'Celltype'}, inplace = True)

    #Make a dictionary to associate plates with patients and plate with cell type
    plate_pt_dict = dict(zip(key.Plate, key.Patient))
    plate_cell_dict = dict(zip(key.Plate, key.Celltype))

    #Now just look at data from selected patient, and apply filters to identify cells with enough reads/amplicon
    pt_allele_plate = df.loc[df['Patient'].isin([pt_id])] 
    pt_allele_plate = pt_allele_plate.drop(columns = 'Patient') #Drop the Patient ID column and other unwanted cols
    pt_allele_plate['Cell_type'] = pt_allele_plate['Plate'].replace(plate_cell_dict)
    pt_allele_plate['Plate_Well'] = pt_allele_plate['Plate'].astype(str) + '_' + pt_allele_plate['Well'].astype(str)

    return pt_allele_plate
    

def call_haps(data, metadata, pt_id, haps, reads,  cutoff):
    
    cond = f'{pt_id}_{haps}'
    print(cond)
    
    if cond == 'JP001_2':
        cols = ['JP001_RUNX1_c', 'JP001_RUNX1_g']
        allcols = ['JP001_RUNX1_c','JP001_RUNX1_g','JP001_SRSF2','JP001_TET2a','JP001_TET2b_c','JP001_TET2b_g']
    elif cond == 'JP001_3':
        cols = ['JP001_RUNX1_g', 'JP001_SRSF2', 'JP001_TET2a']
        allcols = ['JP001_RUNX1_c','JP001_RUNX1_g','JP001_SRSF2','JP001_TET2a','JP001_TET2b_c','JP001_TET2b_g']
    elif cond == 'JP001_4':
        cols = ['JP001_RUNX1_g', 'JP001_SRSF2', 'JP001_TET2a', 'JP001_TET2b_g']
        allcols = ['JP001_RUNX1_c','JP001_RUNX1_g','JP001_SRSF2','JP001_TET2a','JP001_TET2b_c','JP001_TET2b_g']
    elif cond == 'PD7153_3':
        cols = ['PD7153_SRSF2', 'PD7153_TET2a', 'PD7153_TET2b']
        allcols = ['PD7153_CUX1', 'PD7153_SRSF2', 'PD7153_TET2a', 'PD7153_TET2b', 'PD7153_TGFB3_c', 'PD7153_TGFB3_g']
    elif cond == 'PD7153_4':   
        cols = ['PD7153_SRSF2', 'PD7153_TET2a', 'PD7153_TET2b',  'PD7153_CUX1']
        allcols = ['PD7153_CUX1', 'PD7153_SRSF2', 'PD7153_TET2a', 'PD7153_TET2b', 'PD7153_TGFB3_c', 'PD7153_TGFB3_g']
    elif cond == 'PD7151_2': 
        cols = ['PD7151_TET2a', 'PD7151_TET2b']
        allcols = ['PD7151_TET2a', 'PD7151_TET2b']

    elif cond == 'PD7151_3': 
        cols = ['PD7151_TET2a', 'PD7151_TET2b', 'PD7151_SRSF2']
        allcols = ['PD7151_TET2a', 'PD7151_TET2b', 'PD7151_SRSF2']
    else:
        print('For JP001 enter 2/3/4 haplotypes, for PD7153 enter 3/4 haplotypes, for PD7151 enter 2 haplotypes')
    
    #Import information about plate cell type and patient
    key = pd.read_excel(metadata, sheet_name = 'PlateID')
    key = key.drop(['Cell Origin', 'Plate Nr', 'Plate Name','Nr of cells'], axis=1)
    key.rename(columns = {'Comments2':'Plate'}, inplace = True)
    key.rename(columns = {'Cell-group':'Celltype'}, inplace = True)

    
    #Make a dictionary to associate plates with patients and plate with cell type
    plate_pt_dict = dict(zip(key.Plate, key.Patient))
    plate_cell_dict = dict(zip(key.Plate, key.Celltype))
    
    #Group the data and apply filters
    df = data.copy()
    df = df.groupby(['Plate', 'Well', 'Amplicon']).sum().unstack()
    df.columns = df.columns.droplevel(0)
    
    df = df.loc[(df[cols] >= reads).all(axis=1)] #df1 contains just the rows with cells we want - use this to create a filter or key
    df['Plate'] = df.index.get_level_values(0)  #These lines send indexes to columns
    df['Well'] = df.index.get_level_values(1)
    df['Plate_Well'] = df['Plate'].astype(str) + '_' + df['Well'].astype(str)
    wells = df['Plate_Well'].drop_duplicates().to_list() 
    print(f'Cells with {reads} reads for {haps} genes = ', len(wells))

    df2 = data.copy()
    df2 = df2[df2['Plate_Well'].isin(wells)]
    df2 = df2[df2['Amplicon'].isin(cols)]
    
    #Calculate the allele frequency
    df2 = df2.iloc[:, 0:1].unstack(level = 3)
    df2['Total'] = df2.iloc[: , 0] + df2.iloc[: , 1]
    df2['Mut_freq'] = df2.iloc[:, 0]/df2['Total']
    
    #Assign Wt or MT to each allele
    df2 = df2.drop(columns = ['Reads', 'Total'])

    conditions = [(df2['Mut_freq'] <= cutoff), (df2['Mut_freq']) > cutoff ]
    values = ['w', 'm']
    df2['Genotype'] = np.select(conditions, values)
    df2 = df2.drop(columns = ['Mut_freq']).unstack(2)
    df2.columns = df2.columns.droplevel([0,1])
    
    if 'JP001_RUNX1_g' in df2.columns:
        df2.loc[:,'JP001_RUNX1_g'].replace({'w':'R','m':'r' }, inplace = True)
        
    if 'JP001_SRSF2' in df2.columns:  
        df2.loc[:,'JP001_SRSF2'].replace({'w':'S','m':'s' }, inplace = True)
        
    if 'JP001_TET2a' in df2.columns:     
        df2.loc[:,'JP001_TET2a'].replace({'w':'A','m':'a' }, inplace = True)
        
    if 'JP001_TET2b_g' in df2.columns:
        df2.loc[:,'JP001_TET2b_g'].replace({'w':'B','m':'b' }, inplace = True)
        
    if 'JP001_RUNX1_c' in df2.columns:   
        df2.loc[:,'JP001_RUNX1_c'].replace({'w':'C','m':'c' }, inplace = True)
        
    if 'PD7153_SRSF2' in df2.columns:   
        df2.loc[:,'PD7153_SRSF2'].replace({'w':'S','m':'s' }, inplace = True)
        
    if 'PD7153_TET2a' in df2.columns:   
        df2.loc[:,'PD7153_TET2a'].replace({'w':'A','m':'a' }, inplace = True)
        
    if 'PD7153_TET2b' in df2.columns:   
        df2.loc[:,'PD7153_TET2b'].replace({'w':'B','m':'b' }, inplace = True)
        
    if 'PD7153_TGFB3_g' in df2.columns:   
        df2.loc[:,'PD7153_TGFB3_g'].replace({'w':'T','m':'t' }, inplace = True)

    if 'PD7153_CUX1' in df2.columns:   
        df2.loc[:,'PD7153_CUX1'].replace({'w':'C','m':'c' }, inplace = True)    
        
    if 'PD7151_TET2a' in df2.columns:   
        df2.loc[:,'PD7151_TET2a'].replace({'w':'A','m':'a' }, inplace = True)
        
    if 'PD7151_TET2b' in df2.columns:   
        df2.loc[:,'PD7151_TET2b'].replace({'w':'B','m':'b' }, inplace = True)

    if 'PD7151_SRSF2' in df2.columns:   
        df2.loc[:,'PD7151_SRSF2'].replace({'w':'S','m':'s' }, inplace = True)

    
    df2['Haplotype'] = 'x'

    for idx, row in df2.iterrows():
        
        if cond == 'JP001_3':
            a = row['JP001_SRSF2'] + row['JP001_TET2a'] + row['JP001_RUNX1_g']
        elif cond == 'JP001_4':
            a = row['JP001_SRSF2'] + row['JP001_TET2a'] + row['JP001_TET2b_g'] + row['JP001_RUNX1_g']
        elif cond == 'JP001_2':
            a = row['JP001_RUNX1_c'] + row['JP001_RUNX1_g']   
        elif cond == 'PD7153_3':
            a = row['PD7153_TET2b'] + row['PD7153_TET2a'] + row['PD7153_SRSF2']
        elif cond == 'PD7153_4':
            #a = row['PD7153_TET2b'] + row['PD7153_TET2a'] + row['PD7153_SRSF2'] + row['PD7153_TGFB3_g']
            a = row['PD7153_TET2b'] + row['PD7153_TET2a'] + row['PD7153_SRSF2'] + row['PD7153_CUX1']
        elif cond == 'PD7151_2':
            a = row['PD7151_TET2b'] + row['PD7151_TET2a']
        elif cond == 'PD7151_3':
            a = row['PD7151_TET2b'] + row['PD7151_TET2a'] + row['PD7151_SRSF2']
        
        row['Haplotype'] = row['Haplotype'].replace('x', a)   

    df2['Sort_cell_type'] = df2.index.get_level_values(0)
    df2['Sort_cell_type'] = df2['Sort_cell_type'].replace(plate_cell_dict)
    df2['Plate'] = df2.index.get_level_values(0)
    df2['Well'] = df2.index.get_level_values(1)
    df2['Plate_Well'] = df2['Plate'].astype(str) + '_' + df2['Well'].astype(str)
    df2 = df2.drop(columns = cols)
    df2 = df2.drop(columns = ['Plate', 'Well'])
    
    return df2


def plot_hap_dist_sort_type(data, pt_id, haps, reads, cutoff, save = False): #plot based on cell type (need to merge output from above)
      
    #rename the input data and work out how many haplotypes it has
    df3 = data.copy()
    haps = len(df3.iloc[0,0]) #TODO: why isn't haps jsut coming from the input variable?
    cond = f'{pt_id}_{haps}'
    sortcells = sorted(df3['Sort_cell_type'].drop_duplicates().to_list(), reverse = True)
    cellnumb = len(sortcells)
    plotlen = int(math.ceil((cellnumb +1)/2))
    
    #Plot two haplotype data for 3 gene 100 amplicon set - second method to add colour for each haplotype
    fig, axes = plt.subplots(plotlen,2, figsize = (16,plotlen*2))
    fig.subplots_adjust(hspace = 1.2, wspace=.3)
    ax = axes.ravel()
    count = 0

    c = df3['Haplotype'].value_counts().rename_axis('hap').reset_index(name='counts')

    #haplotype dictionary
    hap_dict = {'CR': 'wt', 'Cr': 'r', 'cR': 'c', 'cr': 'cr', 'SAR': 'wt', 'SAr': 'r', 'SaR': 'a', 'Sar': 'ar', 'sAR': 's', 'sAr': 'sr', 'saR': 'sa', 'sar': 'sar', 'SABR': 'wt', 'SABr': 'r', 'SAbR': 'b', 'SAbr': 'br', 'SaBR': 'a', 'SaBr': 'ar', 'SabR': 'ab', 'Sabr': 'abr', 'sABR': 's', 'sABr': 'sr', 'sAbR': 'sb', 'sAbr': 'sbr', 'saBR': 'sa', 'saBr': 'sar', 'sabR': 'sab', 'sabr': 'sabr', 'BA': 'wt', 'Ba': 'a', 'bA': 'b', 'ba': 'ba', 'BAS': 'wt', 'BAs': 's', 'BaS': 'a', 'Bas': 'as', 'bAS': 'b', 'bAs': 'bs', 'baS': 'ba', 'bas': 'bas', 'BAST': 'wt', 'BASt': 't', 'BAsT': 's', 'BAst': 'st', 'BaST': 'a', 'BaSt': 'at', 'BasT': 'as', 'Bast': 'ast', 'bAST': 'b', 'bASt': 'bt', 'bAsT': 'bs', 'bAst': 'bst', 'baST': 'ba', 'baSt': 'bat', 'basT': 'bas', 'bast': 'bast', 'BASC': 'wt', 'BASc': 'c', 'BAsC': 's', 'BAsc': 'sc', 'BaSC': 'a', 'BaSc': 'ac', 'BasC': 'as', 'Basc': 'asc', 'bASC': 'b', 'bASc': 'bc', 'bAsC': 'bs', 'bAsc': 'bsc', 'baSC': 'ba', 'baSc': 'bac', 'basC': 'bas', 'basc': 'basc'}


    #set up correct variables for the number of input haplotypes
    
    # if cond == 'JP001_2':
    #     hap_poss = ['CR', 'Cr', 'cR', 'cr']
        
    # elif cond == 'JP001_3':
    #     hap_poss = ['SAR', 'SAr', 'SaR', 'Sar', 'sAR', 'sAr', 'saR', 'sar']

    # elif cond == 'JP001_4':
    #     hap_poss = ['SABR', 'SABr', 'SAbR', 'SAbr', 'SaBR', 'SaBr', 'SabR', 'Sabr', 'sABR', 'sABr', 'sAbR', 'sAbr', 'saBR', 'saBr', 'sabR', 'sabr']
      
    # elif cond == 'PD7151_2':
    #     hap_poss = ['BA', 'Ba', 'bA', 'ba']

    # elif cond == 'PD7153_3':
    #     hap_poss = ['BAS', 'BAs', 'BaS', 'Bas', 'bAS', 'bAs', 'baS', 'bas']

    # elif cond == 'PD7153_4':
    #     #hap_poss = ['BAST', 'BASt', 'BAsT', 'BAst', 'BaST', 'BaSt', 'BasT', 'Bast', 'bAST', 'bASt', 'bAsT', 'bAst', 'baST', 'baSt', 'basT', 'bast']
    #     hap_poss = ['BASC', 'BASc', 'BAsC', 'BAsc', 'BaSC', 'BaSc', 'BasC', 'Basc', 'bASC', 'bASc', 'bAsC', 'bAsc', 'baSC', 'baSc', 'basC', 'basc']

    if cond == 'JP001_2':
        hap_poss = ['CR', 'Cr', 'cR', 'cr']
        
    elif cond == 'JP001_3':
        hap_poss = ['SAR', 'sAR', 'SaR', 'SAr', 'saR', 'sAr', 'Sar', 'sar']

    elif cond == 'JP001_4':
        hap_poss = ['SABR', 'sABR', 'SaBR', 'SAbR', 'SABr', 'saBR', 'sAbR', 'sABr', 'SabR', 'SaBr','SAbr', 'sabR', 'saBr',  'sAbr','Sabr',  'sabr']
      
    elif cond == 'PD7151_2':
        hap_poss = ['BA', 'bA', 'Ba', 'ba']

    elif cond == 'PD7151_3':
        hap_poss = ['BAS', 'bAS', 'BaS', 'BAs', 'baS', 'bAs', 'Bas', 'bas']

    elif cond == 'PD7153_3':
        hap_poss = ['BAS', 'bAS', 'BaS', 'BAs', 'baS', 'bAs', 'Bas', 'bas']

    elif cond == 'PD7153_4':
        # hap_poss = ['BAST', 'bAST', 'BaST', 'BAsT',  'BASt', 'baST', 'bAsT', 'bASt',  'BasT', 'BaSt', 'BAst', 'basT', 'bAst', 'baSt', 'Bast', 'bast']
        hap_poss = ['BASC', 'bASC', 'BaSC', 'BAsC',  'BASc', 'baSC', 'bAsC', 'bASc',  'BasC', 'BaSc', 'BAsc', 'basC', 'bAsc', 'baSc', 'Basc', 'basc']


    num_col = len(hap_poss)
    cols = sns.color_palette("husl", num_col)     
    color = dict(zip(hap_poss, cols))
    hap_order = {}
    for i, j in enumerate(hap_poss):
        #hap_order[i] = j
        hap_order[j] = i
        
    #if any haplotype is not present, add it into the frame with freq 0 

    for h in hap_poss:
        if h not in str(c['hap']):  #for some reason this needs to be called as a string, wasn't needed outside function
            dfh = pd.DataFrame([[h, 0]], columns= ['hap', 'counts'])
            c = c.append(dfh)
            
    c['order'] = c['hap']
    c = c.replace({'order': hap_order})
    c = c.sort_values(by=['order'])
    d = c['counts'].sum()
    c['proportion'] = c['counts']/d

    #Plot all cells
    sns.barplot(x='hap', y='counts', data = c, palette = color, ax = ax[0], ci = None) #for scatter add  hue = 'hap'
    ax[0].set_title('All cells') 
    ax[0].set_ylabel('Number of cells', fontsize = 11)
    ax[0].set_xlabel('Haplotype', fontsize = 11)
    ax[0].tick_params(axis='x', labelrotation = 90)
    lbs = [hap_dict[l.get_text()] for l in ax[0].get_xticklabels()] #replace haplotype labels for figure
    ax[0].set_xticklabels(lbs)

    
    #Plot by cell type
    for cell in sortcells:
        count += 1
    
        if df3.loc[df3['Sort_cell_type'].isin([cell])].empty == False:

            a = df3.loc[df3['Sort_cell_type'].isin([cell])]['Haplotype'].value_counts().rename_axis('hap').reset_index(name='counts')

            #if any haplotype is not present, add it into the frame with freq 0 hap_3gene_poss has the possibilities
            for h in hap_poss:
                if h not in str(a['hap']):  #for some reason this needs to be called as a string, wasn't needed outside function
                    dfh = pd.DataFrame([[h, 0]], columns= ['hap', 'counts'])
                    a = a.append(dfh)    

            a['order'] = a['hap']
            a = a.replace({'order': hap_order})
            a = a.sort_values(by=['order'])
            b = a['counts'].sum()
            a['proportion'] = a['counts']/b

            sns.barplot(x='hap', y='counts', data = a, palette = color,  ax = ax[count], ci = None) #for scatter add  hue = 'hap'
            ax[count].set_title(str(cell)) 
            ax[count].set_ylabel('Number of cells', fontsize = 11)
            ax[count].set_xlabel('Haplotype', fontsize = 11)
            ax[count].tick_params(axis='x', labelrotation = 90)
            lbs = [hap_dict[l.get_text()] for l in ax[count].get_xticklabels()] #replace haplotype labels for figure
            ax[count].set_xticklabels(lbs)

        else:
            continue
    
    plt.rcParams['svg.fonttype'] = 'none'          
    if save == True:
        fig.savefig(f'../Results/{pt_id}_{haps}_{reads}_{cutoff}_haplotype_by_sortcell.svg',dpi=600)


def plot_index_heatmap_3(data, pt_id, haps, reads, cutoff, save = False):

    '''
    This version puts the scale bar on the right and orders haplotypes from least to most mutated. 
    It also adds a 0 column for any haplotypes not present in the data so the plots shows all possible haplotypes
    '''

    cond = f'{pt_id}_{haps}'
    print(cond)
    
    df = data.copy()
    a = df.groupby(['Haplotype', 'celltype']).size().unstack(fill_value = 0)
    a = a.drop(columns = ['unassigned']) #drop unassigned cells, no need to plot these
    alltypes = ['HSC','MPP','HSC_MPP','MDS_SC', 'CMP',  'GMP','GMP2', 'MEP',  'Neut', 'Mono','nBC']    
    col_order = {}            
    for i, typ in enumerate(alltypes):
         col_order[typ] = i

    #Haplotype plot order
    if cond == 'JP001_2':
        hap_poss = ['CR', 'Cr', 'cR', 'cr']
        
    elif cond == 'JP001_3':
        hap_poss = ['SAR', 'sAR', 'SaR', 'SAr', 'saR', 'sAr', 'Sar', 'sar']

    elif cond == 'JP001_4':
        hap_poss = ['SABR', 'sABR', 'SaBR', 'SAbR', 'SABr', 'saBR', 'sAbR', 'sABr', 'SabR', 'SaBr','SAbr', 'sabR', 'saBr',  'sAbr','Sabr',  'sabr']
      
    elif cond == 'PD7151_2':
        hap_poss = ['BA', 'bA', 'Ba', 'ba']

    elif cond == 'PD7151_3':
        hap_poss = ['BAS', 'bAS', 'BaS', 'BAs', 'baS', 'bAs', 'Bas', 'bas']

    elif cond == 'PD7153_3':
        hap_poss = ['BAS', 'bAS', 'BaS', 'BAs', 'baS', 'bAs', 'Bas', 'bas']

    elif cond == 'PD7153_4':
        hap_poss = ['BASC', 'bASC', 'BaSC', 'BAsC',  'BASc', 'baSC', 'bAsC', 'bASc',  'BasC', 'BaSc', 'BAsc', 'basC', 'bAsc', 'baSc', 'Basc', 'basc']


    hap_order = {}
    for i, j in enumerate(hap_poss):
        hap_order[j] = i
    
    #haplotype dictionary
    hap_dict = {'CR': 'wt', 'Cr': 'r', 'cR': 'c', 'cr': 'cr', 'SAR': 'wt', 'SAr': 'r', 'SaR': 'a', 'Sar': 'ar', 'sAR': 's', 'sAr': 'sr', 'saR': 'sa', 'sar': 'sar', 'SABR': 'wt', 'SABr': 'r', 'SAbR': 'b', 'SAbr': 'br', 'SaBR': 'a', 'SaBr': 'ar', 'SabR': 'ab', 'Sabr': 'abr', 'sABR': 's', 'sABr': 'sr', 'sAbR': 'sb', 'sAbr': 'sbr', 'saBR': 'sa', 'saBr': 'sar', 'sabR': 'sab', 'sabr': 'sabr', 'BA': 'wt', 'Ba': 'a', 'bA': 'b', 'ba': 'ba', 'BAS': 'wt', 'BAs': 's', 'BaS': 'a', 'Bas': 'as', 'bAS': 'b', 'bAs': 'bs', 'baS': 'ba', 'bas': 'bas', 'BAST': 'wt', 'BASt': 't', 'BAsT': 's', 'BAst': 'st', 'BaST': 'a', 'BaSt': 'at', 'BasT': 'as', 'Bast': 'ast', 'bAST': 'b', 'bASt': 'bt', 'bAsT': 'bs', 'bAst': 'bst', 'baST': 'ba', 'baSt': 'bat', 'basT': 'bas', 'bast': 'bast', 'BASC': 'wt', 'BASc': 'c', 'BAsC': 's', 'BAsc': 'sc', 'BaSC': 'a', 'BaSc': 'ac', 'BasC': 'as', 'Basc': 'asc', 'bASC': 'b', 'bASc': 'bc', 'bAsC': 'bs', 'bAsc': 'bsc', 'baSC': 'ba', 'baSc': 'bac', 'basC': 'bas', 'basc': 'basc'}

    a = a.T

    #Add in 0 counts for any haplotypes not present in the data
    for h in hap_poss:
        if h not in a.columns:
            print(hap_dict[h], 'not found in dataset, adding 0 for all cell types')
            a = a.assign(h = 0)
            a.rename(columns={'h':h}, inplace=True)

    a['ct'] = a.index.get_level_values(0)
    a = a.replace({'ct': col_order})
    a = a.sort_values(by=['ct'])
    a = a.drop(columns = ['ct'])
    a = a.T
    a['hap'] = a.index.get_level_values(0)
    a['order'] = a['hap'].replace(hap_order)
    a = a.sort_values(by=['order'])
    a = a.drop(columns = ['hap', 'order'])
    b = a.copy()
    a = a * 100 /a.sum(axis = 0)

    fig, (ax2,ax) = plt.subplots(2, 1, figsize = (8,4), gridspec_kw = dict(height_ratios = [0.8,5.2]))
    
    cbar_ax = fig.add_axes([1, 0.22, .03, 0.565])
    cbar_ax.tick_params(size=0)
    
    sns.heatmap(data = a, ax = ax, robust = True, vmax = 80, cmap = 'YlGnBu_r', cbar_ax=cbar_ax, cbar_kws=dict(ticks=[20, 40, 60, 80]))
    ax.tick_params(axis='y', labelrotation = 0)
    ax.tick_params(axis='x', labelrotation = 90)
    ax.set_xlabel('')
    ax.set_ylabel('')
    lbs = [hap_dict[l.get_text()] for l in ax.get_yticklabels()] #replace haplotype labels for figure
    ax.set_yticklabels(lbs)

    x = df.groupby(['Haplotype', 'celltype']).size().unstack(fill_value = 0)
    x = x.drop(columns = ['unassigned'])
    x = x.sum(axis = 0)
    x = x.to_frame()
    x['order'] = x.index.get_level_values(0)
    x = x.replace({'order': col_order})
    x = x.sort_values(by=['order'])
    x = x.drop(columns = ['order'])
    x['ct'] = x.index.get_level_values(0)
    x.columns = ['number', 'ct']

    
    sns.scatterplot(x = 'ct', y = 'number', data = x, ax = ax2, marker = 's', s = 120, zorder = 10, color = 'black')
    ax2_xlabels = x['number'].to_list()
    ax2.tick_params(axis='x', labelrotation = 0)
    ax2.set_xlabel('')
    ax2.set_ylabel('')
    ax2.axhline(10, ls = '--', lw = 0.5, c = 'gray')
    ax2.axhline(100, ls = '--', lw = 0.5, c = 'gray')
    ax2.axhline(1000, ls = '--', lw = 0.5, c = 'gray')
    ax2.spines['top'].set_visible(False)
    #ax2.spines['right'].set_visible(False)
    ax2.set_ylim(1,1200)
    ax2.set_yticks([1,10, 100, 1000])
    ax2.set_yticklabels(['1','10', '100','1000'])  
    ax2.set_xticklabels(ax2_xlabels)
    ax2.set_yscale('log') 
    
    fig.tight_layout(h_pad = 1) 
    
    plt.rcParams['svg.fonttype'] = 'none'  
    if save == True:
        fig.savefig(f'../Results/{pt_id}_{haps}_{reads}_{cutoff}_haplotype_by_indexcell_cbar_r.svg',bbox_inches='tight',dpi=600)
        fig.savefig(f'../Results/{pt_id}_{haps}_{reads}_{cutoff}_haplotype_by_indexcell_cbar_r.png',bbox_inches='tight',dpi=600)
    
    return b



def calc_scVAF_mod_resc(data, pt_id, reads, draw_plot = False):
    
    '''
    This function takes amplicon read counts for mt and wt and calculates the proportion of mutated alleles in each cell 
    that meets the specified read count for that amplicon.
    All CD34pos cells are treated as a single sample.
    The function returns a plot and data suitable for plotting (eg/ after merging all samples).
    
    '''
    cond = pt_id
    print(cond)


    if cond == 'JP001':  #cols_order is the columns actually being used, can be easily tweaked
        cols_order = ['JP001_SRSF2','JP001_TET2a','JP001_TET2b_g', 'JP001_RUNX1_g']
    elif cond == 'PD7153':
        cols_order = ['PD7153_TET2b','PD7153_SRSF2', 'PD7153_TET2a', 'PD7153_CUX1' ]
    elif cond == 'PD7151': 
        cols_order = ['PD7151_TET2b', 'PD7151_TET2a']

    else:
        print('Enter JP001,  PD7153, or PD7151 as pt_id')

    #Import information about plate cell type and patient
    key = pd.read_excel('../Data/Amp_data/Amplicon_metadata_fixed_anon.xlsx', sheet_name = 'PlateID')
    key = key.drop(['Cell Origin', 'Plate Nr', 'Plate Name','Nr of cells'], axis=1)
    key.rename(columns = {'Comments2':'Plate'}, inplace = True)
    key.rename(columns = {'Cell-group':'Celltype'}, inplace = True)

    #Make a dictionary to associate plates with patients and plate with cell type
    plate_pt_dict = dict(zip(key.Plate, key.Patient))
    plate_cell_dict = dict(zip(key.Plate, key.Celltype))

    #Group the data and apply filters
    df = data.copy()
    df = df.groupby(['Plate', 'Well', 'Amplicon']).sum().unstack()
    df.columns = df.columns.droplevel(0) #remove multi index column names
    df['Plate'] = df.index.get_level_values(0)  #These lines send indexes to columns
    df['Well'] = df.index.get_level_values(1)
    df['Plate_Well'] = df['Plate'].astype(str) + '_' + df['Well'].astype(str)
    df['Sort_cell_type'] = df['Plate'].replace(plate_cell_dict)
    rename = {'CD34+halfCD38-': 'CD34', 'CD34+/38-':'CD34', 'CD34+':'CD34', 'NEs':'Neut', 'Monocytes': 'Mono', 'nBCs': 'nBC'}
    df['Sort_cell_type'].replace(rename, inplace = True)#df now contains cell type as well

    result = {}
    sem_result = {}
    numcell = {}
    cols =[]

    for c in cols_order:

        actual_cols = df.columns

        if c not in actual_cols:
            continue
        else:
            cols.append(c)
            #Work out which wells fit the read criteria, and how many cells there are of each type
            df1 = df.loc[(df[c] >= reads)] #df1 contains just the rows with cells we want - use this to create a filter or key
            wells = df1['Plate_Well'].drop_duplicates().to_list()  #cells that meet the criteria
            cut = len(wells)
            print(f'Cells with at least {reads} reads for amplicon {c}  = ', len(wells))

            if cut == 0:
                print(c, 'has no reads, do not plot')
                continue
            else:

                dfcells = ['Mono', 'Neut', 'nBC', 'CD34']

                for dfc in dfcells:
                    cellno = df1.loc[df1['Sort_cell_type'].isin([dfc])].shape[0]
                    print(dfc, cellno)
                    sample_id = c + ',' + dfc
                    numcell[sample_id] = [cellno] #cellno dictionary now contains counts

                cell_counts = pd.DataFrame.from_dict(numcell, orient = 'Index', columns = ['cell_count'])
                #cell_counts gets created in each iteration from numcell, but numcell gets the counts added in each iteration
                #so the final iteration cell_counts has all the counts. This is what is used in the final output
                cell_counts['labels'] = cell_counts.index.get_level_values(0)
                cell_counts[['Amplicon', 'sort_celltype']] = cell_counts['labels'].str.split(',', expand = True)
                cell_counts = cell_counts.drop(columns = 'labels')

                df2 = data.copy()
                df2 = df2.loc[df2['Plate_Well'].isin(wells)]
                df2 = df2.loc[df2['Amplicon'].isin([c])]


                #Calculate the allele frequency
                df2 = df2.iloc[:, 0:1].unstack(level = 3)
                df2['Total'] = df2.iloc[: , 0] + df2.iloc[: , 1]
                df2['Mut_freq'] = df2.iloc[:, 0]/df2['Total']

                df2 = df2.drop(columns = ['Reads', 'Total'])
                df2 = df2.unstack(2)
                df2.columns = df2.columns.droplevel(0)
                df2['Sort_cell_type'] = df2.index.get_level_values(0)
                df2['Sort_cell_type'] = df2['Sort_cell_type'].replace(plate_cell_dict)

                df2['Sort_cell_type'] = df2['Sort_cell_type'].replace(rename)
                df2.sort_values(by=['Sort_cell_type'], inplace = True)

                #Get means for each amplicon/cell type
                x = df2.copy().groupby(by = 'Sort_cell_type').mean()
                x = x.unstack().to_frame()
                x['celltype'] = x.index.get_level_values(2)
                x['Amplicon'] = x.index.get_level_values(1)
                co = ['VAF', 'sort_celltype', 'Amplicon']
                x.columns = co   


                #Get sem for each amplicon/cell type
                x2 = df2.copy().groupby(by = 'Sort_cell_type').sem()        
                x2 = x2.unstack().to_frame()
                x2['celltype'] = x2.index.get_level_values(2)
                x2['Amplicon'] = x2.index.get_level_values(1)
                co2 = ['sem', 'sort_celltype', 'Amplicon']
                x2.columns = co2  

                result[c] = x.loc[x['Amplicon'].isin([c])]
                sem_result[c] = x2.loc[x2['Amplicon'].isin([c])]

    x_result = pd.concat(result.values(), axis = 0)
    x2_result = pd.concat(sem_result.values(), axis = 0)
    x_x2 = x_result.merge(x2_result, how = 'left', on = ['Amplicon', 'sort_celltype'] )

    final = x_x2.merge(cell_counts, how = 'left', on = ['Amplicon', 'sort_celltype'])

    order = {}
    for a, d in enumerate(cols_order):
        order[d] = a

    x_result['order'] =  x_result['Amplicon'].replace(order)
    x_result.sort_values(by=['order'], inplace = True)

    all_amps = ['PD7153_TET2a',
                'PD7151_TET2a',
                'JP001_TET2a',
                'PD7153_TET2b',
                'PD7151_TET2b',
                'JP001_TET2b_g',
                'PD7153_SRSF2',
                'JP001_SRSF2',
                'PD7153_CUX1',
                'JP001_RUNX1_g'
               ]
    colors = sns.color_palette('husl', n_colors = len(all_amps))
    allVAFcols = dict(zip(all_amps, colors))

    #This plots the mean only
    if draw_plot == True:
        fig, ax = plt.subplots(figsize = (2.5,4))
        sns.scatterplot(x = 'sort_celltype', y = 'VAF', data = x_x2, s = 80,  hue = 'Amplicon', palette = allVAFcols, alpha = 0.5, ax = ax)
        ax.legend(loc = 'upper left', bbox_to_anchor = [1,1], title = 'scVAFs')
        ax.set_ylim(0,0.6)
        ax.axhline(0.1, ls = '--', c = 'silver', zorder = 0)
        ax.axhline(0.2, ls = '--', c = 'silver', zorder = 0)
        ax.axhline(0.3, ls = '--', c = 'silver', zorder = 0)
        ax.axhline(0.4, ls = '--', c = 'silver', zorder = 0)
        ax.axhline(0.5, ls = '--', c = 'silver', zorder = 0)
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.tick_params(axis='x', labelrotation = 90)
        ax.margins(x=0.1)

    return final    


if __name__ == "__main__":

	access_module()