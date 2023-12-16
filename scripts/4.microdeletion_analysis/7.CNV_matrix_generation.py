##### Script purpose: Generate CNV matrix and plot profiles
##### Author: Ege Ulgen
##### Date: Dec 2023

from pathlib import Path

from SigProfilerAssignment.Analyzer import cosmic_fit
from SigProfilerMatrixGenerator.scripts.CNVMatrixGenerator import generateCNVMatrix
from sigProfilerPlotting import plotCNV

# install desired reference genome - no need to run more than once
# from SigProfilerMatrixGenerator import install as genInstall
# genInstall.install(GENOME_BUILD)

GENOME_BUILD = "GRCh38"
INPUT_FILE = "data/TCGA_ASCAT.tsv"
PROJECT = "PMN_project"
OUT_DIR = Path("output")

SPA_OUT_DIR = str(OUT_DIR / "SigProfilerAssignment")
PLOT_OUT_DIR = str(OUT_DIR / "sigProfilerPlotting")

COUNT_MATRIX_PATH = str(OUT_DIR / f"{PROJECT}.CNV48.matrix.tsv")

# Generate CNV matrix
count_mat = generateCNVMatrix(file_type="ASCAT", input_file=INPUT_FILE, project=PROJECT, output_path=str(OUT_DIR))

# assign CN signatures
cosmic_fit(samples=COUNT_MATRIX_PATH, output=SPA_OUT_DIR, input_type="matrix", genome_build=GENOME_BUILD, cosmic_version=3.4, collapse_to_SBS96 = False)

# plotting CNV counts
plotCNV(COUNT_MATRIX_PATH, PLOT_OUT_DIR, f"{PROJECT}_by_sample", percentage=False, aggregate=False)
plotCNV(COUNT_MATRIX_PATH, PLOT_OUT_DIR, f"{PROJECT}_aggregate", percentage=False, aggregate=True)
