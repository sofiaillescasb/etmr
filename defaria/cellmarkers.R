cellmarkers <- list(
  
  # Tumor stem / radial glia-like
  RDG = c(
    "PAX6","SOX2","VIM","NES",
    "HES1","HES5","FABP7","SLC1A3",
    "PTN","HMGA2"
  ),
  
  # Neuroepithelial progenitor
  NEC = c(
    "MDK",
    "MIF",
    "PSMA4",
    "SMS",
    "HSPA5",
    "RANBP1",
    "TPM1",
    "IGFBP2",
    "PFN1",
    "OSTC",
    "HMGA1",
    "LAPTM4A",
    "MYL12A",
    "PDIA6",
    "PRELID1",
    "ARF4",
    "OCIAD2",
    "MAD2L1",
    "PDIA3",
    "KPNA2",
    "AKR1B1",
    "EIF4E",
    "TMED9",
    "MRPL42"
  ),
  
  # IPC
  IPC= c(
    "EOMES", "HES6", "NEUROD4"
  ),
  
  # Immature neurons
  Neuroblasts = c(
    "DCX","NCAM1",
    "NEUROD1","NEUROD2",
    "STMN2","TUBB3"
  ),
  
  # Excitatory neurons
  Glut = c(
    "SLC17A6", "SLC17A7", "SLC17A8"
  ),
  
  # Inhibitory neurons
  GABA = c(
    "GAD1","GAD2","DLX1","DLX2"
  ),
  
  # Mature neurons
  Neurons = c(
    "RBFOX3","SNAP25",
    "SYT1","ELAVL4","NEFM", "MAP2", "GRIN1"
  ),
  
  # Astrocytes
  Astrocytes = c(
    "ALDH1L1","AQP4",
    "GFAP","GLUL","SLC1A2"
  ),
  
  # OPC
  OPC = c(
    "PDGFRA","CSPG4",
    "OLIG1","OLIG2","SOX10"
  ),
  
  # Oligodendrocytes
  Oligodendrocytes = c(
    "MBP","MOG","PLP1",
    "CNP","MOBP"
  ),
  
  # Microglia
  Microglia = c(
    "AIF1","TMEM119",
    "P2RY12","CX3CR1",
    "TYROBP","HEXB"
  ),
  
  # Endothelial
  Endothelial = c(
    "CLDN5","FLT1",
    "KDR","EMCN","PECAM1"
  ),
  
  # Pericytes
  Pericytes = c(
    "PDGFRB","RGS5",
    "CSPG4","DES","MCAM"
  ),
  
  # Fibroblasts
  Fibroblasts = c(
    "COL1A1","COL1A2",
    "DCN","LUM"
  ),
  
  # T cells
  Tcell = c(
    "CD2","CD3D",
    "CD3E","PTPRC"
  ),
  
  # Cell cycle
  Cycling = c(
    "MKI67","TOP2A",
    "AURKA","UBE2C"
  )
)
