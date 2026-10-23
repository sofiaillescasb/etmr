library(qs2)
library(Seurat)
library(clustree)


healthy <- qs_read(
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/human_dev_downsampled_bayesian.qs2", 
  nthreads = 6
)

clus_tree_unint <- clustree(healthy, prefix = "unintegrated_snn_res.")

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/unintegrated/sct/bayesian/clustree.png", 
  clus_tree_unint, 
  units = "in", 
  width = 14, 
  height = 14
  )

healthy_int <- qs_read(
  "/home/sofia/Projects/etmr/defaria/healthy/data/processed/human_dev_downsampled_bayesian_integrated.qs2", 
  nthreads = 6
)

clus_tree_int <- clustree(healthy_int, prefix = "integrated_snn_res.")

ggsave(
  "/home/sofia/Projects/etmr/defaria/healthy/plots/transfered_labels/integrated/sct/bayesian/clustree.png", 
  clus_tree_int, 
  units = "in", 
  width = 14, 
  height = 14
)