#!/usr/bin/env bash

echo "Plotting..."
Rscript "/home/sofia/Projects/etmr/defaria/scripts/03_plot.R" /home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/paper_data_filt_comvar_int.h5ad /home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/merged_ref_paper_data_filt_comvar.h5ad /home/sofia/Projects/etmr/defaria/scRNA/paper_data/data/processed/all_data_int.h5ad
echo "Done."

