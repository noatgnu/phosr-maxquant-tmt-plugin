process PHOSR_MAXQUANT_TMT {
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://cauldron/phosr-maxquant-tmt:1.0.0' :
        'cauldron/phosr-maxquant-tmt:1.0.0' }"

    input:
    path input_file
    path fasta_file
    path annotation_file
    path annotation_protein_file
    val feature_id_col
    val site_col
    val protein_col
    val probability_col
    val min_probability
    path protein_file
    val protein_feature_id_col
    val protein_id_col
    path comparison_file
    val log2_transform
    val exclude_conditions
    val remove_norm_channel
    val col_filter
    val row_filter
    val impute_order
    val impute
    val normalize_method
    val aggregation_order
    val summarization_method
    val adjust_method
    val alpha
    val lfc_threshold
    val organism
    val perform_kinase_analysis
    val kinase_top_substrates
    val kinase_score_threshold
    val kinase_num_motifs
    val top_pathways_plot
    val top_ptmset_plot
    val top_diff_sites_heatmap
    val top_phosphosite_heatmap
    val network_plot_top_interactions

    output:
    
    path "dpa_results.txt", emit: dpa_results, optional: true
    path "dpu_results_phosr_format.txt", emit: dpu_results_phosr_format, optional: true
    path "phosphosite_intensities.txt", emit: phosphosite_intensities, optional: true
    path "kinase_substrate_scores.txt", emit: kinase_substrate_scores, optional: true
    path "kinase_activities.txt", emit: kinase_activities, optional: true
    path "differential_kinase_activities.txt", emit: differential_kinase_activities, optional: true
    path "top_differential_sites.txt", emit: top_differential_sites_data, optional: true
    path "pathway_enrichment_all_comparisons.txt", emit: pathway_enrichment_all_comparisons, optional: true
    path "ptmset_enrichment_all_comparisons.txt", emit: ptmset_enrichment_all_comparisons, optional: true
    path "kinase_substrate_network_all_comparisons.txt", emit: kinase_substrate_network_all_comparisons, optional: true
    path "qc_plots.pdf", emit: qc_plots, optional: true
    path "volcano_plots.pdf", emit: volcano_plots, optional: true
    path "kinase_heatmap.pdf", emit: kinase_heatmap, optional: true
    path "phosphosite_heatmap.pdf", emit: phosphosite_heatmap, optional: true
    path "site_class_distribution.pdf", emit: site_class_distribution, optional: true
    path "top_differential_sites.pdf", emit: top_differential_sites_heatmap, optional: true
    path "kinase_family_analysis.pdf", emit: kinase_family_analysis, optional: true
    path "pathway_enrichment.pdf", emit: pathway_enrichment_plots, optional: true
    path "ptmset_enrichment.pdf", emit: ptmset_enrichment_plots, optional: true
    path "kinase_substrate_network.pdf", emit: kinase_substrate_network_plots, optional: true
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Build arguments dynamically to match CauldronGO PluginExecutor logic
    ARG_LIST=()

    
    # Mapping for row_filter
    VAL="$row_filter"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--row_filter" "\$VAL")
    fi
    
    # Mapping for aggregation_order
    VAL="$aggregation_order"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--aggregation_order" "\$VAL")
    fi
    
    # Mapping for perform_kinase_analysis
    VAL="$perform_kinase_analysis"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        if [ "\$VAL" = "true" ]; then
            ARG_LIST+=("--perform_kinase_analysis")
        fi
    fi
    
    # Mapping for kinase_score_threshold
    VAL="$kinase_score_threshold"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--kinase_score_threshold" "\$VAL")
    fi
    
    # Mapping for top_phosphosite_heatmap
    VAL="$top_phosphosite_heatmap"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--top_phosphosite_heatmap" "\$VAL")
    fi
    
    # Mapping for input_file
    VAL="$input_file"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--input_file" "\$VAL")
    fi
    
    # Mapping for protein_file
    VAL="$protein_file"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--protein_file" "\$VAL")
    fi
    
    # Mapping for top_diff_sites_heatmap
    VAL="$top_diff_sites_heatmap"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--top_diff_sites_heatmap" "\$VAL")
    fi
    
    # Mapping for organism
    VAL="$organism"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--organism" "\$VAL")
    fi
    
    # Mapping for protein_feature_id_col
    VAL="$protein_feature_id_col"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--protein_feature_id_col" "\$VAL")
    fi
    
    # Mapping for probability_col
    VAL="$probability_col"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--probability_col" "\$VAL")
    fi
    
    # Mapping for min_probability
    VAL="$min_probability"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--min_probability" "\$VAL")
    fi
    
    # Mapping for exclude_conditions
    VAL="$exclude_conditions"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--exclude_conditions" "\$VAL")
    fi
    
    # Mapping for remove_norm_channel
    VAL="$remove_norm_channel"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        if [ "\$VAL" = "true" ]; then
            ARG_LIST+=("--remove_norm_channel")
        fi
    fi
    
    # Mapping for normalize_method
    VAL="$normalize_method"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--normalize_method" "\$VAL")
    fi
    
    # Mapping for lfc_threshold
    VAL="$lfc_threshold"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--lfc_threshold" "\$VAL")
    fi
    
    # Mapping for kinase_num_motifs
    VAL="$kinase_num_motifs"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--kinase_num_motifs" "\$VAL")
    fi
    
    # Mapping for annotation_file
    VAL="$annotation_file"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--annotation_file" "\$VAL")
    fi
    
    # Mapping for site_col
    VAL="$site_col"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--site_col" "\$VAL")
    fi
    
    # Mapping for comparison_file
    VAL="$comparison_file"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--comparison_file" "\$VAL")
    fi
    
    # Mapping for log2_transform
    VAL="$log2_transform"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        if [ "\$VAL" = "true" ]; then
            ARG_LIST+=("--log2_transform")
        fi
    fi
    
    # Mapping for alpha
    VAL="$alpha"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--alpha" "\$VAL")
    fi
    
    # Mapping for fasta_file
    VAL="$fasta_file"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--fasta_file" "\$VAL")
    fi
    
    # Mapping for annotation_protein_file
    VAL="$annotation_protein_file"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--annotation_protein_file" "\$VAL")
    fi
    
    # Mapping for col_filter
    VAL="$col_filter"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--col_filter" "\$VAL")
    fi
    
    # Mapping for impute_order
    VAL="$impute_order"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--impute_order" "\$VAL")
    fi
    
    # Mapping for impute
    VAL="$impute"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--impute" "\$VAL")
    fi
    
    # Mapping for adjust_method
    VAL="$adjust_method"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--adjust_method" "\$VAL")
    fi
    
    # Mapping for kinase_top_substrates
    VAL="$kinase_top_substrates"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--kinase_top_substrates" "\$VAL")
    fi
    
    # Mapping for top_pathways_plot
    VAL="$top_pathways_plot"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--top_pathways_plot" "\$VAL")
    fi
    
    # Mapping for feature_id_col
    VAL="$feature_id_col"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--feature_id_col" "\$VAL")
    fi
    
    # Mapping for protein_col
    VAL="$protein_col"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--protein_col" "\$VAL")
    fi
    
    # Mapping for protein_id_col
    VAL="$protein_id_col"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--protein_id_col" "\$VAL")
    fi
    
    # Mapping for summarization_method
    VAL="$summarization_method"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--summarization_method" "\$VAL")
    fi
    
    # Mapping for top_ptmset_plot
    VAL="$top_ptmset_plot"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--top_ptmset_plot" "\$VAL")
    fi
    
    # Mapping for network_plot_top_interactions
    VAL="$network_plot_top_interactions"
    if [ -n "\$VAL" ] && [ "\$VAL" != "null" ] && [ "\$VAL" != "[]" ]; then
        ARG_LIST+=("--network_plot_top_interactions" "\$VAL")
    fi
    
    Rscript --vanilla --slave /app/phosr_maxquant_tmt_modular.R --args \
        "\${ARG_LIST[@]}" \
        --output_folder . \
        \${args:-}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PhosR MaxQuant TMT Analysis: 1.0.0
    END_VERSIONS
    """
}
