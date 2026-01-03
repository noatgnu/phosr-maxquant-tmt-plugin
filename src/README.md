# PhosR MaxQuant TMT - Modular Structure

This plugin has been refactored into a modular structure for better maintainability and code organization.

## Directory Structure

```
phosr-maxquant-tmt/
├── phosr_maxquant_tmt_modular.R  # Main entry point
├── plugin.yaml                    # Plugin configuration
├── r-packages.txt                 # R package dependencies
├── README.md                      # User documentation
└── src/                           # Modular components
    ├── README.md                  # This file
    ├── utils.R                    # Utility functions
    ├── ptm_processing.R           # PTM site processing
    ├── data_loading.R             # Data loading and preparation
    ├── normalization.R            # Filtering, normalization, imputation
    ├── statistical_analysis.R     # Differential analysis with limma
    ├── msqrob2_analysis.R         # MSqRob2 model fitting and hypothesis testing
    ├── visualization.R            # All plotting functions
    ├── kinase_analysis.R          # PhosR kinase activity analysis
    ├── pathway_analysis.R         # Pathway enrichment analysis
    └── signaling_analysis.R       # Signaling network analysis
```

## Module Descriptions

### utils.R
- `detect_delimiter()`: Auto-detect file delimiter (CSV/TSV)

### ptm_processing.R
- `process_ptm_with_fasta()`: FASTA-based PTM site mapping with probability thresholds

### data_loading.R
- `load_and_prepare_ptm_data()`: Load and process PTM evidence file
- `load_and_prepare_protein_data()`: Load and process protein data
- `create_se_from_data()`: Create SummarizedExperiment objects

### normalization.R
- `filter_and_normalize_se()`: Filter missing values, normalize, and impute PTM data
- `filter_and_normalize_protein()`: Same for protein data

### statistical_analysis.R
- `perform_differential_analysis()`: Limma-based differential analysis
- `save_intensity_matrix()`: Save processed intensity matrices

### msqrob2_analysis.R
- `perform_msqrob2_model_fitting()`: Fit MSqRob2 models to protein and PTM data
  - Handles run effects
  - Calculates DPU (Differential Protein Use)
  - Saves protein intensities
- `perform_msqrob2_hypothesis_testing()`: Perform hypothesis testing for all comparisons
  - Generates DPA (Differential Peptide Abundance) results
  - Generates DPU (Differential Protein Use) results
  - Generates protein-level differential results
- `generate_msqrob2_plots()`: Generate MSqRob2-style QC and volcano plots
  - Protein and PTM boxplots
  - PCA plots for protein and PTM data
  - Volcano plots for protein, DPA, and DPU results

### visualization.R
- `generate_qc_plots()`: QC plots (boxplots, PCA, correlation)
- `generate_volcano_plots()`: Volcano plots for differential results
- `generate_site_class_distribution()`: pS/pT/pY distribution plots
- `generate_top_differential_sites_heatmap()`: Heatmap of top differential sites
- `generate_kinase_family_analysis()`: Kinase family enrichment and activity
- `generate_phosphosite_heatmap()`: Phosphosite-kinase score heatmap
- `generate_kinase_activity_plots()`: Kinase activity heatmaps and bar charts

### kinase_analysis.R
- `perform_kinase_analysis()`: Complete PhosR kinase analysis workflow
  - Sequence window extraction
  - Kinase-substrate scoring
  - Kinase activity calculation
  - Differential kinase activity
  - Orchestrates all kinase-related visualizations

### pathway_analysis.R
- `perform_pathway_enrichment()`: Pathway enrichment analysis using PhosR
  - Reactome pathway enrichment
  - Gene Ontology enrichment
  - Visualization of enriched pathways

### signaling_analysis.R
- `perform_signaling_analysis()`: Signaling network analysis
  - PTM-SET enrichment analysis
  - Kinase-substrate network analysis
  - Signaling network visualization

## Working Directory

The backend executor (`script_executor.go`) automatically sets the working directory to the plugin location (`plugins/phosr-maxquant-tmt/`) before executing the R script. This ensures all relative imports (e.g., `source("src/utils.R")`) work correctly.

The script can safely use relative paths for:
- Sourcing module files from `src/`
- Any other file operations within the plugin directory

## Execution Flow

1. Backend sets working directory to plugin location
2. Load libraries
3. Source all module files from `src/`
4. Parse command-line arguments
5. Execute analysis workflow:
   - Load and prepare PTM data
   - Create SummarizedExperiment
   - Normalize and impute
   - Process protein data (optional, for DPU analysis)
   - Differential analysis (MSqRob2 or limma)
   - Generate QC plots and volcano plots
   - Kinase analysis (optional)
     - Kinase-substrate scoring
     - Kinase activity calculation
     - Differential kinase activity
   - Pathway enrichment analysis (optional)
     - Reactome pathways
     - Gene Ontology terms
   - Signaling network analysis (optional)
     - PTM-SET enrichment
     - Kinase-substrate networks

## Advantages of Modular Structure

- **Maintainability**: Each module has a clear responsibility
- **Reusability**: Functions can be easily reused in other contexts
- **Testing**: Individual modules can be tested independently
- **Readability**: Smaller, focused files are easier to understand
- **Collaboration**: Multiple developers can work on different modules
