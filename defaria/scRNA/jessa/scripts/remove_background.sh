#! /usr/bin/env bash
# i halved learning rate because the learning curve looked strange and the documentation suggested that 
cellbender remove-background \
                 --input /home/sofia/Projects/etmr/defaria/scRNA/jessa/data/raw/raw_feature_bc_matrix.h5 \
                 --output /home/sofia/Projects/etmr/defaria/scRNA/jessa/data/processed/jessa_background.h5ad \
                 --learning-rate 5e-5 \ 
                 --cuda \

