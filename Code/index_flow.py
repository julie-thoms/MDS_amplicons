#Module for importing and plotting fcs files, mostly for indexed files

#Authors:
#	Julie Thoms
#	Annatina Schnegg-Kaufmann
#	Fabio Zanini

#Version 1
#12th August 2021

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

import fcsparser


#Define functions

def about():

	'''
	The functions in this module are used to import and QC indexed fcs files, then assign cell types for normal or MDS stem cells, or neutrophils, naive B cells, Monocytes.

    The callable functions are;
        count_files(directory)
        get_channels(directory)
        plate_qc(directory, data_name, save = False)
        get_comp_data(directory, plate_key, channel_key, plot = False)
        get_comp_data_autolabel(directory, plate_key, plot = False)
        flowplot_byplate(compdata, plot_list, logs, gates,  data_name, plot = True, save = False)
        flowplot_bycelltype(assigndata, plot_list, logs, gates,  data_name, plot = True, save = False)
        MDS_BM_celltype_assign(source, gates, data_name, save = False)
        BM_celltype_assign(source, gates, data_name, save = False)
        PB_celltype_assign(source, gates, data_name, save = False)



	For cell type assignment to work, the following flow channels and labels should be used:
	For CD34pos cells;
		channel_key = {'YG582/15-A': 'CD34-PE', 
		   'YG670/30-A': 'Lin-PE-Cy5', 
		   'YG780/60-A': 'CD123-PE-Cy7', 
		   'V450/50-A': 'CD90-BV421', 
		   'V610/20-A': 'Zombie', 
		   'B530/30-A': 'CD45RA-FITC', 
		   'R660/20-A': 'IL1RAP-APC', 
		   'R780/60-A': 'CD38-APC-cy7'

	For mature cell types;
		channel_key = {
			'YG582/15-A': 'CD16-PE', 
			'YG670/30-A': 'CD14-Pe-Cy5', 
			'YG780/60-A': 'CD56-PE-Cy7', 
			'V450/50-A': 'CD66b-BV421', 
			'V610/20-A': 'Zombie', 
			'B530/30-A': 'CD45-FITC', 
			'B695/40-A': 'IgD-BB700', 
			'R660/20-A': 'Cd27-APC'
}


	'''
	print('This module contains functions to import and QC indexed fcs files. For more info type index_flow.about?')


def count_files(directory): 
	
	'''
	This function counts the number of files in an input directory.
	To call - count_files(directory)
	The function returns the number of files.
	'''
	
	count = 0
	for filename in os.listdir(directory):
		count += 1
		
	return count


def get_inx(meta : dict):  
	
	'''
	This function retrieves the wells which contain sort data from an already imported indexed fcs file.
	To call - get_inx(meta : dict)
	meta is an output from meta, data =  fcsparser.parse(filename, reformat_meta=True)
	The function returns the sort locations as a list, using the default software assignment x(0-23), y(0-15) for Aria
	'''
	
	i = 1
	key = f'INDEX SORTING LOCATIONS_{i}'
	sort_locs = [] #Sort locations from the metadata, numerical row then column starting at 0
	while key in meta:
		sort_locs.append(meta[key])
		i += 1
		key = f'INDEX SORTING LOCATIONS_{i}'
	sort_locs = ''.join(sort_locs).split(';')
	return sort_locs

def get_channels(directory):
    
    '''
    This functions returns a dictionary of the channels used for each plate
    
    '''
    out_dict = {}
    for count, filename in enumerate(os.listdir(directory)):
        fn = os.path.join(directory, filename)
        plateid = filename

        meta, data = fcsparser.parse(fn, reformat_meta=True)

        comp_fields = meta['SPILL'].split(',')
        n = int(comp_fields[0])
        channels = comp_fields[1: n+1]

        out_dict[plateid] = channels   
    
    return out_dict    


def plate_qc(directory, data_name, save = False): 
	
	'''
	This function gives a visual QC output of which wells in each plate contain a sort event.
	To call - plate_qc(directory, data_name)
	data_name is a string that will be used to label graphs and output files
	Sort locations are plotted using the default software assignment x(0-23), y(0-15) for Aria, with alphanumeric values applied to the graph only.
	When save = True, save location is ../Results/
	
	'''
	
	#count the plates
	plotlen = int(math.ceil(count_files(directory)/2))  #math.ceil rounds up to account for an odd number of plates

	fig, ax = plt.subplots(plotlen, 2, figsize = (16,(plotlen*5)))
	ax = ax.ravel()
	fig.subplots_adjust(hspace = 0.3, wspace=.3)

	xwell = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
	ywell = [0,1,2,3,4,5,6,7,8,9,0,11,12,13,14,15]
	allwells = []

	for a in ywell: #this creates a list of all possible locations
		for b in xwell:
			allwells.append(str(a) + ','+ str(b))

	#Read in files and plot the data for each

	for count, filename in enumerate(os.listdir(directory)):
		fn = os.path.join(directory, filename)
		plateid = filename

		meta, data = fcsparser.parse(fn, reformat_meta=True)
		
		sort_locs = get_inx(meta)

		for well in allwells:
			if well in sort_locs:  #if the well is in the list of wells with a sorted cell the colour will be darker
				alpha = 0.8
			else:
				alpha = 0.1        #plot empty wells in a lighter shade

			y, x = well.split(',')
			ax[count].scatter(x, y, alpha=alpha, color='gray')      
		ax[count].set_title(plateid)
		ax[count].invert_yaxis()  #flip the axis so the plate order looks natural
		ax[count].set_xticklabels(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'])
		ax[count].set_yticklabels(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'])

	fig.suptitle(f'{data_name}', fontsize=14)

	plt.rcParams['svg.fonttype'] = 'none'  
	if save == True:   
		fig.savefig(f'../Results/{data_name}_plot_sorted_by_plate.svg',dpi=600)     
	
	return


def get_comp_data(directory, plate_key, channel_key, plot = False):
	
	'''
	This function retrieves and applies compensation to indexed fcs files with an option to plot the compensation matrix.
	To call - get_comp_data(directory, plate_key, channel_key, plot = False)
	plate_key is a dictionary that matches fcs filename and user-defined plate name
	channel_key is a dictionary that matches laser/filter with the user-defined channel names
	The function ignores any user-applied channel labels in the metadata
	The function returns a dataframe containing the compensated values for each plate/well/channel
	
	'''
	

	data_dict = {}
	data_dict_comp = {}
	
	#count the plates
	plotlen = int(math.ceil(count_files(directory)/2))  #math.ceil rounds up to account for an odd number of plates
	
	#Set up comp plots
	
	fig, ax = plt.subplots(plotlen, 2, figsize = (12,(plotlen*4)))
	fig.subplots_adjust(hspace = 0.3, wspace=1.5)
	ax = ax.ravel()

	for count, filename in enumerate(os.listdir(directory)):
		fn = os.path.join(directory, filename)
		plateid = filename

		meta, data = fcsparser.parse(fn, reformat_meta=True)

		sort_locs = get_inx(meta)

		wells = []
		for loc in sort_locs:
			if loc == '':
				continue
			row_index, col = loc.split(',')
			col = str(int(col)+1)
			row = chr(65 + int(row_index)) 
			well = row+col
			wells.append(well)   #Wells is a list of well locations with data derived from the index file

		# Get antibodies and rename with well names
		channel_idx = [int(x[2:-1]) for x in meta if x.startswith('$P') and x.endswith('S')]

		data.index = pd.Index(wells, name='Sorted well') #Renames index with well name

		# Load compensation
		comp_fields = meta['SPILL'].split(',')
		n = int(comp_fields[0])
		channels = comp_fields[1: n+1]
		matrix = np.asarray(comp_fields[n+1:]).astype(np.float64).reshape((n, n)) 

		matrix = numpy.linalg.inv(matrix)

		spill_matrix = pd.DataFrame(
			matrix,
			index=channels,
			columns=channels,
		)  #spill_matrix is the comp matrix

		#Plots compensation matrix
		# Reorder the dyes by wavelength
		wls = [int(x.split('/')[0][-3:]) for x in channels]
		idx = np.argsort(wls)
		spill_by_wls = spill_matrix.iloc[idx].T.iloc[idx].T

		sns.heatmap(spill_by_wls, ax=ax[count])
		ax[count].set_title(plateid)

		#Apply comp and replace column names with antibodies
		data_comp = data.copy()
		compensation = spill_matrix

		for channel in channels:
			data_comp[channel] = compensation.loc[:,channel].values @ data[channels].values.T #@ for matrix multiplication

		#Store df and compdf for this iteration as a unique variable - 
		data_comp.rename(columns = channel_key, inplace = True)
		plate = plate_key.get(filename) #pull the plate name for this file and use to label output dictionaries
		data_comp['Well'] = data_comp.index.get_level_values(0)
		data_comp['Plate'] = plate  #add new column with plate name
		data_dict_comp[plate] = data_comp #output df into a dictionary

	fig.tight_layout()  
	
	if plot == False:
		plt.close() 

	alldata_comp = pd.concat(data_dict_comp.values(), axis = 0)
	alldata_comp['Plate_Well'] = alldata_comp['Plate'].astype(str) + '_' + alldata_comp['Well'].astype(str)
	
	return alldata_comp


def get_comp_data_autolabel(directory, plate_key, plot = False):
	
	'''
	This function retrieves and applies compensation to indexed fcs files with an option to plot the compensation matrix. 
	Channel names applied during fcs generation are used (but will cause issues if naming is inconsistent between plates) - if so use get_comp_data() instead.
	To call - get_comp_data(directory, plate_key, channel_key, plot = False)
	plate_key is a dictionary that matches fcs filename and user-defined plate name
	The function ignores any user-applied channel labels in the metadata
	The function returns a dataframe containing the compensated values for each plate/well/channel
	
	'''
	

	data_dict = {}
	data_dict_comp = {}
	
	#count the plates
	plotlen = int(math.ceil(count_files(directory)/2))  #math.ceil rounds up to account for an odd number of plates
	
	#Set up comp plots
	
	fig, ax = plt.subplots(plotlen, 2, figsize = (12,(plotlen*4)))
	fig.subplots_adjust(hspace = 0.3, wspace=1.5)
	ax = ax.ravel()

	for count, filename in enumerate(os.listdir(directory)):
		fn = os.path.join(directory, filename)
		plateid = filename

		meta, data = fcsparser.parse(fn, reformat_meta=True)

		sort_locs = get_inx(meta)

		wells = []
		for loc in sort_locs:
			if loc == '':
				continue
			row_index, col = loc.split(',')
			col = str(int(col)+1)
			row = chr(65 + int(row_index)) 
			well = row+col
			wells.append(well)   #Wells is a list of well locations with data derived from the index file

		# Get antibodies and rename with well names
		channel_idx = [int(x[2:-1]) for x in meta if x.startswith('$P') and x.endswith('S')]
		channeld = {meta['_channel_names_'][i-1]: meta[f'$P{i}S'] for i in channel_idx} #automatically retrieve labels

		data.index = pd.Index(wells, name='Sorted well') #Renames index with well name

		# Load compensation
		comp_fields = meta['SPILL'].split(',')
		n = int(comp_fields[0])
		channels = comp_fields[1: n+1]
		matrix = np.asarray(comp_fields[n+1:]).astype(np.float64).reshape((n, n)) 

		matrix = numpy.linalg.inv(matrix)

		spill_matrix = pd.DataFrame(
			matrix,
			index=channels,
			columns=channels,
		)  #spill_matrix is the comp matrix

		#Plots compensation matrix
		# Reorder the dyes by wavelength
		wls = [int(x.split('/')[0][-3:]) for x in channels]
		idx = np.argsort(wls)
		spill_by_wls = spill_matrix.iloc[idx].T.iloc[idx].T

		sns.heatmap(spill_by_wls, ax=ax[count])
		ax[count].set_title(plateid)

		#Apply comp and replace column names with antibodies
		data_comp = data.copy()
		compensation = spill_matrix

		for channel in channels:
			data_comp[channel] = compensation.loc[:,channel].values @ data[channels].values.T #@ for matrix multiplication

		#Store df and compdf for this iteration as a unique variable - 
		data_comp.rename(columns = channeld, inplace = True)
		plate = plate_key.get(filename) #pull the plate name for this file and use to label output dictionaries
		data_comp['Well'] = data_comp.index.get_level_values(0)
		data_comp['Plate'] = plate  #add new column with plate name
		data_dict_comp[plate] = data_comp #output df into a dictionary

	fig.tight_layout()  
	
	if plot == False:
		plt.close() 

	alldata_comp = pd.concat(data_dict_comp.values(), axis = 0)
	alldata_comp['Plate_Well'] = alldata_comp['Plate'].astype(str) + '_' + alldata_comp['Well'].astype(str)
	
	return alldata_comp


def flowplot_byplate(compdata, plot_list, logs, gates,  data_name, plot = True, save = False):

	'''
	This function plots retrieved and compensatd indexed fcs files, with colour applied according to the sorted plate.
	To call = flowplot_byplate(compdata, plot_list, logs, gates,  data_name, plot = True, save = False)
	compdata is the dataframe containing the compensated index data, eg/ the output from get_comp_data() or get_comp_data_autolabel()
	plot_list is a list of channel pairs to plot eg/ plot_list = [['FSC-A', 'SSC-A'],['FSC-A', 'FSC-W'],['SSC-A', 'SSC-H']]
	logs is a list of channels which should be plotted on a log10 scale (eg/ all antibody channels)
	gates = a dictionary of gate locations for all relevant channels, these will be added to plots  eg/ gates = {'Lin-PE-Cy5': 1500}
	data_name is a string which will be used to customise the name of the output file {data_name}_flowplot_by_plate.png
	When save = True, the plot will be saved in ../Results/

	'''
	
	sourcedata = compdata.copy()
	plotlen = int(math.ceil(len(plot_list)/2))

	plates = sourcedata['Plate'].drop_duplicates().to_list()
	cols = sns.color_palette('husl', n_colors = len(plates))
	palette = dict(zip(plates, cols))

	sourcedata['Colour'] = sourcedata['Plate'].map(palette)

	fig, axs = plt.subplots(plotlen,2, figsize = (12,(plotlen*4)))
	axs = axs.ravel()
	fig.subplots_adjust(hspace = 0.3, wspace=1.5)

	for ax,y in zip(axs, plot_list):
		x_label = y[0]
		y_label = y[1]
		ax.scatter(sourcedata[x_label] + 11, sourcedata[y_label] + 11, alpha = 0.3, c = sourcedata['Colour'], s = 5)
		ax.set_xlabel(x_label)
		ax.set_ylabel(y_label)
	   #ax.legend()
		if x_label in logs:
			ax.set_xscale('log')
			ax.set_xlim(left = 10)
		if y_label in logs:
			ax.set_yscale('log')
			ax.set_ylim(bottom = 10)
		if x_label in gates:
			ax.axvline(gates[x_label], ls = '--', c = 'k')
		if y_label in gates:
			ax.axhline(gates[y_label], ls = '--', c = 'k')   
		
		for hap in palette:
			point = ax.scatter([], [], color=palette[hap], s = 20, alpha = 0.3, label=hap)
			ax.add_artist(point)
			ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,title='Plate')
		
		ax.autoscale_view()  

	fig.suptitle(f'{data_name}', fontsize=16) 

	plt.rcParams['svg.fonttype'] = 'none'  
	if save == True:
		fig.savefig(f'../Results/{data_name}_flowplot_by_plate.svg',dpi=600) 

	if plot == False:
		plt.close()


def flowplot_bycelltype(assigndata, plot_list, logs, gates,  data_name, plot = True, save = False):

	'''
	This function plots retrieved and compensatd indexed fcs files, with colour applied according to the assigned cell type.
	To call = flowplot_bycelltype(assigndata, plot_list, logs, gates,  data_name, plot = True, save = False)
	assigndata is the dataframe containing the compensated index data, eg/ the output from MDS_BM_celltype_assign(), BM_celltype_assign(), PB_celltype_assign()
	plot_list is a list of channel pairs to plot eg/ plot_list = [['FSC-A', 'SSC-A'],['FSC-A', 'FSC-W'],['SSC-A', 'SSC-H']]
	logs is a list of channels which should be plotted on a log10 scale (eg/ all antibody channels)
	gates = a dictionary of gate locations for all relevant channels, these will be added to plots  eg/ gates = {'Lin-PE-Cy5': 1500}
	data_name is a string which will be used to customise the name of the output file {data_name}_flowplot_by_celltype.png
	When save = True, the plot will be saved in ../Results/
    '''

	


	sourcedata = assigndata.copy()
	plotlen = int(math.ceil(len(plot_list)/2))

	plates = sourcedata['celltype'].drop_duplicates().to_list()
	cols = sns.color_palette('husl', n_colors = len(plates))
	palette = dict(zip(plates, cols))

	sourcedata['Colour'] = sourcedata['celltype'].map(palette)

	fig, axs = plt.subplots(plotlen,2, figsize = (12,(plotlen*4)))
	axs = axs.ravel()
	fig.subplots_adjust(hspace = 0.3, wspace=1.5)

	for ax,y in zip(axs, plot_list):
		x_label = y[0]
		y_label = y[1]
		ax.scatter(sourcedata[x_label] + 11, sourcedata[y_label] + 11, alpha = 0.3, c = sourcedata['Colour'], s = 5)
		ax.set_xlabel(x_label)
		ax.set_ylabel(y_label)
		#ax.legend()
		if x_label in logs:
			ax.set_xscale('log')
			ax.set_xlim(left = 10)
		if y_label in logs:
			ax.set_yscale('log')
			ax.set_ylim(bottom = 10)
		if x_label in gates:
			ax.axvline(gates[x_label], ls = '--', c = 'k')
		if y_label in gates:
			ax.axhline(gates[y_label], ls = '--', c = 'k')   
		
		for hap in palette:
			point = ax.scatter([], [], color=palette[hap], s = 20, alpha = 0.3, label=hap)
			ax.add_artist(point)
			ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,title='Cell type')
		
		ax.autoscale_view()  

	fig.suptitle(f'{data_name}', fontsize=16)   

	plt.rcParams['svg.fonttype'] = 'none'  
	if save == True:
		fig.savefig(f'../Results/{data_name}_flowplot_by_celltype.svg',dpi=600) 

	if plot == False:
		plt.close()        



def MDS_BM_celltype_assign(source, gates, data_name, save = False):

	'''
	This function applies gates defined externally and assigns CD34pos cells to subsets (healthy stem cells, mds stem cells, CMP, GMP, MEP)
	To call - MDS_BM_celltype_assign(source, gates, data_name, save = False)
	source is a dataframe containing indexed flow data eg/ the output from get_comp_data() or get_comp_data_autolabel()
	gates = a dictionary of gate locations for all relevant channels, these will be added to plots  eg/ gates = {'Lin-PE-Cy5': 1500}
	data_name is a string which will be used to customise the name of the output file {data_name}_MDS_BMcelltype_assigned.png
	The cell assignment plot will be automatically be saved in ../Results/
	When save = True, the input dataframe with a new column 'celltype' containing the assigned will be saved. The function always returns this dataframe.
	The function assumes sort gates used for sorting actually dumped dead and lin pos cells

	'''

	data_in = source.copy()

	#Create boolean matrix based on gates in new columns

	data_in['CD34_pos'] = data_in['CD34-PE'] >= gates['CD34-PE']
	data_in['CD38_pos'] = data_in['CD38-APC-cy7'] >= gates['CD38-APC-cy7']
	data_in['CD38_neg'] = data_in['CD38-APC-cy7'] < gates['CD38-APC-cy7']
	data_in['IL1RAP_pos'] = data_in['IL1RAP-APC'] >= gates['IL1RAP-APC']
	data_in['IL1RAP_neg'] = data_in['IL1RAP-APC'] < gates['IL1RAP-APC']
	data_in['CD45RA_pos'] = data_in['CD45RA-FITC'] >= gates['CD45RA-FITC']
	data_in['CD45RA_neg'] = data_in['CD45RA-FITC'] < gates['CD45RA-FITC']
	data_in['CD123_pos'] = data_in['CD123-PE-Cy7'] >= gates['CD123-PE-Cy7']
	data_in['CD123_neg'] = data_in['CD123-PE-Cy7'] < gates['CD123-PE-Cy7']
	data_in['CD90_pos'] = data_in['CD90-BV421'] >= gates['CD90-BV421']
	data_in['CD90_neg'] = data_in['CD90-BV421'] < gates['CD90-BV421']

	#Define each cell type - doing this explicitly to make code easier to follow. Tweaked to include IL1RAP neg for normal 38- cells

	mds_sc1 = ['CD34_pos','CD38_neg','IL1RAP_pos'] #MDS sc are postivie for any of IL1RAP, CD45RA, CD123
	mds_sc2 = ['CD34_pos','CD38_neg','CD45RA_pos']
	mds_sc3 = ['CD34_pos','CD38_neg','CD123_pos']
	sc = ['CD34_pos','CD38_neg','IL1RAP_neg','CD45RA_neg','CD123_neg']
	CMP = ['CD34_pos','CD38_pos','CD45RA_neg','CD123_pos']
	gmp = ['CD34_pos','CD38_pos','CD45RA_pos']
	mep = ['CD34_pos','CD38_pos','CD45RA_neg','CD123_neg']

	col_names = ['mds_sc1', 'mds_sc2','mds_sc2','sc', 'cmp', 'gmp','mep']
	celltypes = ['MDS_SC', 'MDS_SC','MDS_SC','SC',  'CMP', 'GMP', 'MEP'] #will rename all mds_sc to a single label
	alltypes = [mds_sc1, mds_sc2, mds_sc3, sc, CMP, gmp, mep]

	#Store list of which columns are true for each cell type in a dictionary

	criteria = dict(zip(col_names, alltypes))

	#Use calculated criteria to add cell type true/false columns

	for x in col_names:
		data_in[x] = (data_in[criteria[x]] == True).all(axis = 1)

	data_in['unassigned'] = (data_in[col_names] == False).all(axis = 1)

	#Create a new column that contains the cell type

	col_names.append('unassigned')  #account for cells that were not assigned
	celltypes.append('unassigned')

	for a, b in zip(col_names, celltypes):

		data_in.loc[data_in[a] == True, 'celltype'] = b 

	#Define colour palette here, and make a new column for it

	ct = data_in['celltype'].drop_duplicates().to_list()
	col = sns.color_palette('husl', n_colors = len(ct))
	palette = dict(zip(ct, col))
	data_in['Colour'] = data_in['celltype'].map(palette)

	#Create a file containing assigned cell type for each well
	if save == True:
		data_in.to_csv(f'../Data/{data_name}_index_refined.tsv', sep = '\t')    #Save the big df to a file
		#wellID = data_in[['Plate_Well', 'celltype']] #don't need this and it make sthe fucntion fail for bulk cell

	#Make a plot of cell type distributions

	fig, ax = plt.subplots()

	c = data_in['celltype'].value_counts().rename_axis('cell').reset_index(name='counts')
	ax = sns.barplot(x='cell', y='counts', data = c, palette = palette)

	fig.suptitle(f'{data_name}', fontsize=16)
	plt.rcParams['svg.fonttype'] = 'none'  
	fig.savefig(f'../Results/{data_name}_MDS_BMcelltype_assigned.svg',dpi=600) 
	
	return data_in

def BM_celltype_assign(source, gates, data_name, save = False):

	'''
	This function applies gates defined externally and assigns CD34pos cells to subsets (HSC, MPP, LMPP, CMP, GMP, MEP)
	To call - BM_celltype_assign(source, gates, data_name, save = False)
	source is a dataframe containing indexed flow data eg/ the output from get_comp_data() or get_comp_data_autolabel()
	gates = a dictionary of gate locations for all relevant channels, these will be added to plots  eg/ gates = {'Lin-PE-Cy5': 1500}
	data_name is a string which will be used to customise the name of the output file {data_name}_BMcelltype_assigned.png
	The cell assignment plot will be automatically be saved in ../Results/
	When save = True, the input dataframe with a new column 'celltype' containing the assigned will be saved. The function always returns this dataframe.
	The function assumes sort gates used for sorting actually dumped dead and lin pos cells

	'''

	data_in = source.copy()

	#Create boolean matrix based on gates in new columns

	data_in['CD34_pos'] = data_in['CD34-PE'] >= gates['CD34-PE']
	data_in['CD38_pos'] = data_in['CD38-APC-cy7'] >= gates['CD38-APC-cy7']
	data_in['CD38_neg'] = data_in['CD38-APC-cy7'] < gates['CD38-APC-cy7']
	data_in['CD45RA_pos'] = data_in['CD45RA-FITC'] >= gates['CD45RA-FITC']
	data_in['CD45RA_neg'] = data_in['CD45RA-FITC'] < gates['CD45RA-FITC']
	data_in['CD123_pos'] = data_in['CD123-PE-Cy7'] >= gates['CD123-PE-Cy7']
	data_in['CD123_neg'] = data_in['CD123-PE-Cy7'] < gates['CD123-PE-Cy7']
	data_in['CD90_pos'] = data_in['CD90-BV421'] >= gates['CD90-BV421']
	data_in['CD90_neg'] = data_in['CD90-BV421'] < gates['CD90-BV421']

	#Define each cell type - doing this explicitly to make code easier to follow. 


	hsc = ['CD34_pos','CD38_neg','CD45RA_neg','CD123_neg','CD90_pos'] 
	mpp = ['CD34_pos','CD38_neg','CD45RA_neg','CD123_neg','CD90_neg']
	lmpp = ['CD34_pos','CD38_neg','CD45RA_pos','CD123_neg']
	CMP = ['CD34_pos','CD38_pos','CD45RA_neg','CD123_pos']
	gmp = ['CD34_pos','CD38_pos','CD45RA_pos']
	mep = ['CD34_pos','CD38_pos','CD45RA_neg','CD123_neg']

	col_names = ['hsc', 'mpp','lmpp', 'cmp', 'gmp','mep']
	celltypes = ['HSC', 'MPP','LMPP',  'CMP', 'GMP', 'MEP']
	alltypes = [hsc, mpp, lmpp,  CMP, gmp, mep]

	#Store list of which columns are true for each cell type in a dictionary

	criteria = dict(zip(col_names, alltypes))

	#Use calculated criteria to add cell type true/false columns

	for x in col_names:
		data_in[x] = (data_in[criteria[x]] == True).all(axis = 1)

	data_in['unassigned'] = (data_in[col_names] == False).all(axis = 1)

	#Create a new column that contains the cell type

	col_names.append('unassigned')  #account for cells that were not assigned
	celltypes.append('unassigned')

	for a, b in zip(col_names, celltypes):

		data_in.loc[data_in[a] == True, 'celltype'] = b 

	#Define colour palette here, and make a new column for it

	ct = data_in['celltype'].drop_duplicates().to_list()
	col = sns.color_palette('husl', n_colors = len(ct))
	palette = dict(zip(ct, col))
	data_in['Colour'] = data_in['celltype'].map(palette)

	#Create a file containing assigned cell type for each well
	if save == True:
		data_in.to_csv(f'../Data/{data_name}_index_refined.tsv', sep = '\t')    #Save the big df to a file

	#Make a plot of cell type distributions

	fig, ax = plt.subplots()

	c = data_in['celltype'].value_counts().rename_axis('cell').reset_index(name='counts')
	ax = sns.barplot(x='cell', y='counts', data = c, palette = palette)

	plt.rcParams['svg.fonttype'] = 'none'  
	fig.suptitle(f'{data_name}', fontsize=16)
	fig.savefig(f'../Results/{data_name}_BMcelltype_assigned.svg',dpi=600) 
	
	return data_in


def PB_celltype_assign(source, gates, data_name, save = False):

	'''
	This function applies gates defined externally and assigns mature cell types(neutrophils (NE), monocytes (Mono), naive B cells (nBC))
	To call - PB_celltype_assign(source, gates, data_name, save = False)
	source is a dataframe containing indexed flow data eg/ the output from get_comp_data() or get_comp_data_autolabel()
	gates = a dictionary of gate locations for all relevant channels, these will be added to plots  eg/ gates = {'Lin-PE-Cy5': 1500}
	data_name is a string which will be used to customise the name of the output file {data_name}_PBcelltype_assigned.png
	The cell assignment plot will be automatically be saved in ../Results/
	When save = True, the input dataframe with a new column 'celltype' containing the assigned will be saved. The function always returns this dataframe.
	The function assumes sort gates used for sorting actually dumped dead cells

	'''

	data_in = source.copy()
	
	#Create boolean matrix based on gates in new columns
	
	data_in['CD45_pos'] = data_in['CD45-FITC'] >= gates['CD45-FITC']
	data_in['IgD_pos'] = data_in['IgD-BB700'] >= gates['IgD-BB700']
	data_in['IgD_neg'] = data_in['IgD-BB700'] < gates['IgD-BB700']
	data_in['CD27_pos'] = data_in['Cd27-APC'] >= gates['Cd27-APC']
	data_in['CD27_neg'] = data_in['Cd27-APC'] < gates['Cd27-APC']
	data_in['CD66b_pos'] = data_in['CD66b-BV421'] >= gates['CD66b-BV421']
	data_in['CD66b_neg'] = data_in['CD66b-BV421'] < gates['CD66b-BV421']
	data_in['CD16_pos'] = data_in['CD16-PE'] >= gates['CD16-PE']
	data_in['CD16_neg'] = data_in['CD16-PE'] < gates['CD16-PE']
	data_in['CD14_pos'] = data_in['CD14-Pe-Cy5'] >= gates['CD14-Pe-Cy5']
	data_in['CD14_neg'] = data_in['CD14-Pe-Cy5'] < gates['CD14-Pe-Cy5']


	#Define each cell type - doing this explicitly to make code easier to follow. Tweaked to include IL1RAP neg for normal 38- cells

	nBC = ['CD45_pos','IgD_pos','CD27_neg']
	Neut = ['CD45_pos','IgD_neg','CD66b_pos','CD16_pos']
	Mono = ['CD45_pos','IgD_neg','CD66b_neg']

	col_names = ['nBC', 'Neut', 'Mono']
	celltypes = ['nBC', 'Neut', 'Mono']
	alltypes = [nBC, Neut, Mono]

	#Store list of which columns are true for each cell type in a dictionary

	criteria = dict(zip(col_names, alltypes))


	#Use calculated criteria to add cell type true/false columns

	for x in col_names:
		data_in[x] = (data_in[criteria[x]] == True).all(axis = 1)

	data_in['unassigned'] = (data_in[col_names] == False).all(axis = 1)

	#Create a new column that contains the cell type

	col_names.append('unassigned')  #account for cells that were not assigned
	celltypes.append('unassigned')

	for a, b in zip(col_names, celltypes):

		data_in.loc[data_in[a] == True, 'celltype'] = b 

	#Define colour palette here, and make a new column for it

	ct = data_in['celltype'].drop_duplicates().to_list()
	col = sns.color_palette('husl', n_colors = len(ct))
	palette = dict(zip(ct, col))
	data_in['Colour'] = data_in['celltype'].map(palette)

	#Create a file containing assigned cell type for each well
	if save == True:
		data_in.to_csv(f'../Data/{data_name}_index.tsv', sep = '\t')    #Save the big df to a file

	#Make a plot of cell type distributions

	fig, ax = plt.subplots()

	c = data_in['celltype'].value_counts().rename_axis('cell').reset_index(name='counts')
	ax = sns.barplot(x='cell', y='counts', data = c, palette = palette)

	plt.rcParams['svg.fonttype'] = 'none'  
	fig.suptitle(f'{data_name}', fontsize=16)
	fig.savefig(f'../Results/{data_name}_PBcelltype_assigned.svg',dpi=600) 
	
	return data_in




if __name__ == "__main__":

	access_module()