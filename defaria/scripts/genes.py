s_genes =  ["MCM5",     "PCNA",     "TYMS",     "FEN1",     "MCM2",     "MCM4",    
"RRM1",     "UNG",      "GINS2",    "MCM6",     "CDCA7",    "DTL",     
"PRIM1",    "UHRF1",    "MLF1IP",   "HELLS",    "RFC2",     "RPA2",    
"NASP",     "RAD51AP1", "GMNN",     "WDR76",    "SLBP",     "CCNE2",   
"UBR7",     "POLD3",    "MSH2",     "ATAD2",    "RAD51",    "RRM2",    
"CDC45",    "CDC6",     "EXO1",     "TIPIN",    "DSCC1",    "BLM",
"CASP8AP2", "USP1",     "CLSPN",    "POLA1",    "CHAF1B",   "BRIP1",
"E2F8"]




g2m_genes = ["HMGB2",   "CDK1",   "NUSAP1",  "UBE2C",   "BIRC5",   "TPX2",   
"TOP2A",  "NDC80",   "CKS2",    "NUF2",    "CKS1B",   "MKI67",  
"TMPO",    "CENPF",   "TACC3",   "FAM64A",  "SMC4",    "CCNB2",  
"CKAP2L",  "CKAP2",   "AURKB",   "BUB1",    "KIF11",   "ANP32E", 
"TUBB4B",  "GTSE1",   "KIF20B",  "HJURP",   "CDCA3",   "HN1",    
"CDC20",   "TTK",     "CDC25C",  "KIF2C",   "RANGAP1", "NCAPD2", 
"DLGAP5",  "CDCA2",   "CDCA8",   "ECT2",    "KIF23",   "HMMR",
"AURKA",   "PSRC1",   "ANLN",    "LBR" ,    "CKAP5",     "CENPE",
"CTCF",    "NEK2",    "G2E3",    "GAS2L3",  "CBX5",    "CENPA"]

etmr1_markers = {
  "Astrocytes": ["APOE", "SLC1A3", "GFAP", "SLC1A2"], # ASTROCITOS
  "Microglia": ["AIF1", "TYROBP"],  # MICROGLIA
  "IPC": ["ASCL1", "GADD45G", "BTG2", "CCND2", "EOMES", "NEUROD1"], # IPC
  "RDG": ["PAX6", "SOX2", "VIM", "NES", "CCND1", "PTN", "HMGA2", "HES1", "HES5", "FABP7", "SLC1A3"], # RDG
  "Neurons (general)": ["SNAP25", "DCX", "DNER", "RBFOX3"], # NEURONA
  "GABAergic": ["GAD1", "GAD2", "SLC6A1", "SLC32A1", "DLX1", "DLX2"], # GABAÉRGICAS
  "Cycling": ["TOP2A", "MKI67", "AURKA"], # Cycling
  "Glutamatergic": ["SLC17A6", "EMX1", "TBR1", "SLC17A7", "GRIN1", "CAMK2A"] # GLUTAMATÉRGICAS
}



prov_ann = {
    "Endothelial": ["CLDN5", "MCAM"],
    "Red blood cells": ["HEMGN"],
    "VLMC": ["COL1A1", "VTN"],
    "Pericytes": ["CSPG4", "PDGFRB"],
    "Cajal_Retzius": ["RELN", "LHX5"],
    "Astrocytes": ["APOE", "ALDH1L1", "SLC1A3", "GFAP", "SLC1A2"],
    "Oligodendrocytes": ["OLIG1", "OLIG2", "PDGFRA", "CSPG4", "NKX2-2", "SOX10"],
    "Microglia": ["AIF1", "TMEM119", "TYROBP"],
    "Neurons (general)": ["SNAP25", "DCX", "DNER"],
    "GABA": ["GAD1", "GAD2", "DLX2", "DLX1"],
    "Glutamtergic": ["NEUROD2", "NEUROD6", "SLC17A6", "GRIN1"],
    "Cycling": ["TOP2A", "MKI67", "AURKA"],
    "IPC": ["GADD45G", "BTG2", "CCND2", # general
         "ASCL1",                    # inhibitorias
         "NEUROG2"]  ,
    "Neuroepithelial-like RDG": ["PAX6", "SOX2"],
    "RDG": ["NES", "CCND1", "PTN", "HMGA2", "HES1", "HES5", "FABP7", "SLC1A3", "VIM"],
    "Neuroepithelial-like": ["PRTG", "WNT8B", "ID3"],
    "Neuroepithelial_sig" : [
"MDK",
"MIF",
# "AC090498.1",
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
]
}

