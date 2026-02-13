#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHOSR_MAXQUANT_TMT } from './modules/local/phosr-maxquant-tmt/main'

workflow PIPELINE {
    main:
    PHOSR_MAXQUANT_TMT (
        params.input_file ? Channel.fromPath(params.input_file).collect() : Channel.of([]),
        params.fasta_file ? Channel.fromPath(params.fasta_file).collect() : Channel.of([]),
        params.annotation_file ? Channel.fromPath(params.annotation_file).collect() : Channel.of([]),
        params.annotation_protein_file ? Channel.fromPath(params.annotation_protein_file).collect() : Channel.of([]),
        Channel.value(params.feature_id_col ?: ''),
        Channel.value(params.site_col ?: ''),
        Channel.value(params.protein_col ?: ''),
        Channel.value(params.probability_col ?: ''),
        Channel.value(params.min_probability ?: ''),
        params.protein_file ? Channel.fromPath(params.protein_file).collect() : Channel.of([]),
        Channel.value(params.protein_feature_id_col ?: ''),
        Channel.value(params.protein_id_col ?: ''),
        params.comparison_file ? Channel.fromPath(params.comparison_file).collect() : Channel.of([]),
        Channel.value(params.log2_transform ?: ''),
        Channel.value(params.exclude_conditions ?: ''),
        Channel.value(params.remove_norm_channel ?: ''),
        Channel.value(params.col_filter ?: ''),
        Channel.value(params.row_filter ?: ''),
        Channel.value(params.impute_order ?: ''),
        Channel.value(params.impute ?: ''),
        Channel.value(params.normalize_method ?: ''),
        Channel.value(params.aggregation_order ?: ''),
        Channel.value(params.summarization_method ?: ''),
        Channel.value(params.adjust_method ?: ''),
        Channel.value(params.alpha ?: ''),
        Channel.value(params.lfc_threshold ?: ''),
        Channel.value(params.organism ?: ''),
        Channel.value(params.perform_kinase_analysis ?: ''),
        Channel.value(params.kinase_top_substrates ?: ''),
        Channel.value(params.kinase_score_threshold ?: ''),
        Channel.value(params.kinase_num_motifs ?: ''),
        Channel.value(params.top_pathways_plot ?: ''),
        Channel.value(params.top_ptmset_plot ?: ''),
        Channel.value(params.top_diff_sites_heatmap ?: ''),
        Channel.value(params.top_phosphosite_heatmap ?: ''),
        Channel.value(params.network_plot_top_interactions ?: ''),
    )
}

workflow {
    PIPELINE ()
}
