"""
Gene list module for GASLIT-AF Variant Analysis.
Defines and processes gene lists for analysis.
"""

import logging

# Configure logging
log = logging.getLogger("gaslit-af")

# Comprehensive GASLIT-AF gene list including new additions
GASLIT_AF_GENES_TEXT = """
# Immune & Inflammatory System
IDO2 AHR AHRR IL36RN CFH MBL2 NLRP3 IL1B IL6 IL17 IL13 IL4 HLA-DQB1 PTPN22 CTLA4 ASXL1 CBL DNMT3B ETV6 IDH1 IL6R

# Autonomic & Neurotransmitter System
COMT CHRM2 DRD2 GABRA1 CHRNA7 ADRB1 ADRB2 NOS3 GNB3 SLC6A2 NET EZH2 SLC6A4 HTR2A TAAR1 OPRM1 GCH1 TRPV2 MYT1L NRXN3

# Structural & Connective Tissue
TNXB ADAMTS10 SELENON NEB MYH7 MAPRE1 ADGRV1 PLXNA2 COL3A1 FBN1 FLNA COL5A1 FKBP14 PLOD1 CDON SULF2

# Metabolic System
APOE PCSK9 UGT1A1 HNF1A ABCC8 TFAM C19orf12 MT-ATP6 MT-ATP8 PDHA1 SDHB NAMPT NMRK1 PGC1A PRKAA2

# Endocannabinoid System
CNR1 CNR2 FAAH MGLL

# Calcium & Ion Channels
ITPR1 KCNJ5 RYR2 KCNA5 KCND3 KCNE1 KCNQ1 HCN4 CAMK2B

# Mast Cell Activation
TPSAB1 KIT HNMT TET2

# Kynurenine Pathway
IDO1 KMO KYNU TDO2 HAAO ARNT BECN1 ATG5

# Vascular & RAS System
ROCK1 ROCK2 ARG1 "Ang-(1-7)" ACE ACE2 "ANG I" "ANG II" TGFβ1 TGFβ2 TGFβ3 GDF-15 "Activin B" Follistatin "Hif-1α"

# Mitochondrial & Cellular Stress
DRP1 "PINK-1" SIRT1 IFNα IFNβ IFNγ IFNL1 PGE2 "α-NAGA" ATG13 NEFL S100B TWEAK S100PBP AKAP1 USP6NL

# Cardiac Development & Conduction
PITX2 SPEN KIAA1755 GATA4 GATA5 GATA6 TBX3 TBX5 NKX2-5 ZFHX3 GREM2 NPPA SCN5A SH3PXD2A MYL4 LMNA

# ME/CFS & Post-Viral Syndromes
S100PBP AKAP1 USP6NL CDON SULF2
"""

# Known SNPs of interest based on the gene variant map
KNOWN_SNPS = {
    # Cognition & Brain Function
    "CHRM2": ["rs8191992", "rs2350780"],
    "DRD2": ["rs6277"],
    "TFAM": ["rs1937"],
    "BCL2": ["rs956572"],
    "ST8SIA6": [],  # Add specific SNPs if available
    "CHRNA5": [],
    "NRG1": [],
    "MAPRE1": [],
    "GYPC": [],
    "CABP5": [],
    
    # Sleep Traits
    "ADA": ["rs73598374"],
    "VRK1": [],
    "RALYL": [],
    "FOXO6": [],
    
    # Cardiovascular Health
    "VKORC1": [],
    "CYP2C9": [],
    "CYP2C19": [],
    "GP6": [],
    "PTGS1": [],
    
    # Atrial Fibrillation Related
    "PITX2": [],  # Transcription factor implicated in cardiac development
    "IL6R": [],   # Immune responses related to AF
    "SPEN": [],   # Hormone-inducible transcriptional coregulator
    "KIAA1755": [], # Associated with heart rate variability
    "GATA4": [],  # Cardiac morphogenesis transcription factor
    "GATA5": [],  # Cardiac morphogenesis transcription factor
    "GATA6": [],  # Cardiac morphogenesis transcription factor
    "KCNA5": [],  # Potassium channel subunit
    "KCND3": [],  # Potassium channel subunit
    "KCNE1": [],  # Potassium channel subunit
    "KCNQ1": [],  # Potassium channel gene
    "HCN4": [],   # Pacemaker channel protein
    "PRKAA2": [], # AMP-activated protein kinase
    "CAMK2B": [], # Calcium/calmodulin-dependent protein kinase II
    "TBX3": [],   # Transcription factor for cardiac conduction
    "TBX5": [],   # Transcription factor for cardiac conduction
    "NKX2-5": [], # Homeobox-containing transcription factor
    "ZFHX3": [],  # Cardiac muscle function and development
    "GREM2": [],  # Cardiac laterality and atrial rhythm regulation
    "NPPA": [],   # Atrial natriuretic peptide
    "SCN5A": [],  # Cardiac sodium channel
    "SH3PXD2A": [], # Cellular signaling pathways
    "MYL4": [],   # Atrial-specific myosin light chain
    "LMNA": [],   # Nuclear lamina proteins
    
    # ME/CFS Related
    "S100PBP": [], # S100P-binding protein
    "AKAP1": [],  # A-kinase anchoring protein 1, mitochondrial function
    "USP6NL": [], # USP6 N-terminal like protein, GTPase regulation
    "CDON": [],   # Cell adhesion associated, oncogene regulated
    "SULF2": [],  # Sulfatase 2, extracellular sulfatase
    
    # Rare & Neurological Conditions
    "ADGRV1": ["rs575602255", "rs555466095"],
    "C19orf12": ["rs146170087"],
    "PRSS1": ["rs202003805", "rs1232891794"],
    "ATM": ["rs531617441"]
}

# Build a reverse lookup dictionary for SNP to gene mapping
SNP_TO_GENE = {}
for gene, snps in KNOWN_SNPS.items():
    for snp in snps:
        SNP_TO_GENE[snp] = gene

def parse_gene_list():
    """
    Parse the GASLIT-AF gene list text into a set of gene symbols.
    Handles quoted multi-word genes correctly.
    
    Returns:
        set: Set of gene symbols
    """
    genes = set()
    in_quotes = False
    current_gene = ""
    
    for char in GASLIT_AF_GENES_TEXT:
        if char == '"':
            in_quotes = not in_quotes
            if not in_quotes and current_gene:  # End of quoted gene
                genes.add(current_gene.strip())
                current_gene = ""
        elif in_quotes:
            current_gene += char
        elif char.isspace():
            if current_gene:
                genes.add(current_gene.strip())
                current_gene = ""
        else:
            current_gene += char
    
    # Add the last gene if there is one
    if current_gene:
        genes.add(current_gene.strip())
    
    log.info(f"Parsed {len(genes)} GASLIT-AF genes")
    return genes

# Initialize the gene set
GASLIT_AF_GENES = parse_gene_list()
