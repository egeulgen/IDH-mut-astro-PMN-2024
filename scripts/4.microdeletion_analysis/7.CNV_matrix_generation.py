##### Script purpose: Generate CNV matrices and plot profiles
##### Author: Ege Ulgen
##### Date: Dec 2023

from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna
import sigProfilerPlotting as sigPlt

input_file = "data/TCGA_ASCAT.tsv" 
output_path = "output/"
project = "PMN-latest"
scna.generateCNVMatrix("ASCAT", input_file, project, output_path)


matrix_path = f"output/{project}.CNV48.matrix.tsv"
output_path = "output/sigProfilerPlotting/"

sigPlt.plotCNV(matrix_path, output_path, f"{project}_by_sample", percentage=False, aggregate=False)
sigPlt.plotCNV(matrix_path, output_path, f"{project}_aggregate", percentage=False, aggregate=True)
