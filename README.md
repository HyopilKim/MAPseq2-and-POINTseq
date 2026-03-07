# MAPseq2-and-POINTseq
Codes and scripts related to MAPseq2 (https://www.biorxiv.org/content/10.1101/2025.06.23.661165v1) and POINTseq (https://www.biorxiv.org/content/10.1101/2025.06.25.661405v1).

1. System requirements

MATLAB (MathWorks, version R2022a)
R (version 4.4.1 or later)
No non-standard hardware is required.



2. Installation guide

MATLAB
No installation is required. Download functions from GitHub
https://github.com/HyopilKim/MAPseq2-and-POINTseq/MAPseq2/MATLAB/helperfunctions
and add the folder to the MATLAB path.

R
Install the required R package: MetaNeighbor (version 1.25.0)
Typical installation time on a normal desktop computer: < 1 minute.



3. Demo

The dataset included in the repository can be used to run the full analysis.

MATLAB analysis
Instructions for generating raw barcode matrices from FASTQ files, creating projection matrices, and clustering MAPseq1 and MAPseq2 projection matrices are provided in:
https://github.com/HyopilKim/MAPseq2-and-POINTseq/MAPseq2/MATLAB/MAPseq1andAMPseq2clustering.mlx

All raw barcode matrices and projection matrices are provided in the same directory, which also can be used as demo.

Expected output:
* Raw barcode matrices
* Projection matrices
* Clustering results for MAPseq1 and MAPseq2 projection matrices
Expected run time on a normal desktop computer: < 1 hour

R analysis
MetaNeighbor comparison between MAPseq1 and MAPseq2 can be performed using:
https://github.com/HyopilKim/MAPseq2-and-POINTseq/MAPseq2/R/metaneighbor.R

Expected output:
* MetaNeighbor comparison of MAPseq1 and MAPseq2
Expected run time on a normal desktop computer: < 2 hours



4. Instructions for use

All raw sequencing datasets are available in the NCBI Sequence Read Archive (SRA) under the following BioSample accession numbers:
MAPseq virus library (SAMN53393094)
Barrel cortex (SAMN49783975)
Motor cortex diluted replicates (SAMN49783974)
MAPseq2-1 (SAMN53373349)
MAPseq2-2 (SAMN53373327)

Processed MAPseq data are available in Supplementary File 2 and also in the GitHub repository:
https://github.com/HyopilKim/MAPseq2-and-POINTseq/MAPseq2

To start from FASTQ files, follow the instructions in:
https://github.com/HyopilKim/MAPseq2-and-POINTseq/MAPseq2/Bash/Bash_processing.txt
and
https://github.com/HyopilKim/MAPseq2-and-POINTseq/MAPseq2/MATLAB/MAPseq1andAMPseq2clustering.mlx

To start from the provided raw barcode matrices or projection matrices, follow the instructions in:
https://github.com/HyopilKim/MAPseq2-and-POINTseq/MAPseq2/MATLAB/MAPseq1andAMPseq2clustering.mlx
