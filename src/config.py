"""Configuration. Edit the paths at the top; everything else is the grid."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ----------------------------------------------------------------- binaries
# Built per INSTALL.md. Verified: IQ-TREE 2.4.0, SimPhy from adamallo/SimPhy.
IQTREE_BIN = ROOT.parent / "iqtree-2.4.0-Linux-intel" / "bin" / "iqtree2"
SIMPHY_BIN = ROOT.parent / "SimPhy" / "bin" / "simphy"

# --------------------------------------------------------------- directories
DATA_DIR   = ROOT / "data"
WORK_DIR   = ROOT / "work"
RESULT_DIR = ROOT / "results"

# =============================================================== PART A
# The real S100 benchmark (Zhang, Rabiee, Sayyari & Mirarab 2018; also the
# evaluation set for wASTRAL, Zhang & Mirarab 2022). NOT bundled -- see
# INSTALL.md for the download. Set S100_DIR to the extracted root.
S100_DIR = DATA_DIR / "S100"

# Filenames inside each replicate, as distributed. run_s100.py --check reports
# exactly what it finds rather than assuming any of this is present.
S100_SPECIES_TREE = "s_tree.trees"
S100_TRUE_GENE_GLOB = "g_trees*.trees"
S100_EST_GENE_TEMPLATE = (
    "bestMLestimatedgenetree/"
    "estimatedgenetre_{length}.gtr.rerooted.final.contracted.non{suffix}")
S100_SUPPORT = {
    # suffix -> (support maximum, human label)
    ".abayes": (1.0, "aBayes"),
    "":        (100.0, "bootstrap"),
}
S100_LENGTHS = [200, 400, 800, 1600]
S100_REPLICATES = [f"{i:02d}" for i in range(1, 11)]
S100_MAX_GENES = 200          # per replicate/length; raise once it runs

# =============================================================== PART B
# Regenerating S100-protocol data with SimPhy when you need conditions the
# published dataset does not contain (e.g. model misspecification).
#
# These values approximate the published S100 protocol. VERIFY THEM AGAINST
# Zhang et al. (2018) before putting any of them in a write-up -- they are
# starting points, not quotations.
SIM_N_TAXA      = 51
SIM_BIRTH_RATE  = 1e-7          # per generation
SIM_TREE_HEIGHT = 2_500_000     # generations
SIM_POP_SIZE    = 400_000       # haploid effective population size
SIM_SUBST_RATE  = 1e7           # e:rate notation
SIM_N_LOCI      = 50
SIM_SEQ_LENGTHS = [200, 1600]

# (generating model, analysis model, label).
# "correct" is the consistency control: ML is consistent there, so error is
# stochastic. "misspecified" drops rate heterogeneity from the analysis model,
# which makes part of the error systematic.
SIM_MODEL_PAIRS = [
    ("GTR{1.0,2.0,1.0,1.0,2.0}+G4{0.15}", "GTR+G4", "correct"),
    ("GTR{1.0,2.0,1.0,1.0,2.0}+G4{0.15}", "JC",     "misspecified"),
]
SIM_REPLICATES = 1

# =============================================================== PART C
# A real biological gene-tree set, bundled. 3,679 UCE loci, 48 taxa, RAxML
# with rapid bootstrap, from the Avian Phylogenomics Project (Jarvis et al.
# 2014, Science 346:1320-1331), redistributed in the TREE-QMC tutorial
# (Han & Molloy 2023, Genome Res 33:1042-1052).
AVIAN_TREES = DATA_DIR / "avian" / "avian_uce_trees_3679.tre"
AVIAN_SUPPORT_MAX = 100.0

MASTER_SEED = 20260804
CONFIDENT_FRACTION = 0.95      # "confident" = support >= this * support_max
