#!/bin/bash -l
#SBATCH --job-name=bulkRNAseq
#SBATCH --partition=tier1q
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load gcc/12.1.0
module load miniconda3/24.4.0

conda activate /gpfs/data/bioinformatics/software/conda_envs/cri-bulk-rnaseq-report-v1.0

multiqc results -c /gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/app/multiqc_config.yml -f --no-data-dir

R -e "shiny::runApp('/gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/app', host = '0.0.0.0', port = 3838)"


