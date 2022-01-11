# MDS_amplicons
Amplicon sequencing analysis

Create a local folder /MDS_amplicons/Results/ to hold output files (this file path is in .gitignore)

Create another local folder /MDS_amplicons/Data/ to hold input fcs files and any code-generated datafiles (this file path is in .gitignore)

This folder requires some sub-folders;
/Amp_data   #read count data
/VAF_data  #bulk VAF data
/JP001_BM   #fcs files
/JP001_PB    #fcs files
/PD7151_BM   #fcs files
/PD7151_PB    #fcs files
/PD7153_BM   #fcs files
/PD7153_PB    #fcs files

#fcs files are available at https://flowrepository.org/id/FR-FCM-Z4PR


#Module: index_flow

The functions in this module are used to import and QC indexed fcs files, then assign cell types for normal or MDS stem cells, or neutrophils, naive B cells, Monocytes.

The callable functions are;
    about()
    count_files(directory)
    get_inx(meta : dict)
    get_channels(directory)
    plate_qc(directory, data_name, save = False)
    get_comp_data(directory, plate_key, channel_key, plot = False)
    get_comp_data_autolabel(directory, plate_key, plot = False)
    flowplot_byplate(compdata, plot_list, logs, gates,  data_name, plot = True, save = False)
    flowplot_bycelltype(assigndata, plot_list, logs, gates,  data_name, plot = True, save = False)
    flowplot_bycelltype_gating(assigndata, logs, gates, data_name, plot = True, save = False)
    MDS_BM_celltype_assign(source, gates, data_name, save = False)
    BM_celltype_assign(source, gates, data_name, save = False)
    PB_celltype_assign(source, gates, data_name, save = False)


#Module: index_haps

The functions in this module are used to assign haplotypes to single cells genotyped with PCR amplicons.

The callable functions are;
    about()
    data_retrieval2(sourcefile, metadata, pt_id)
    call_haps(data, pt_id, haps, reads,  cutoff)
    plot_hap_dist_sort_type(data, pt_id, haps, reads, cutoff, save = False)
    plot_index_heatmap_3(data, title, haps, reads, cutoff, save = False)
    calc_scVAF_mod_resc(data, pt_id, reads, draw_plot = False)