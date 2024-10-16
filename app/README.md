# RNAseq analysis routine workflow 

## Step 1. Transfer raw data to the project folder on Randi
Please check data completeness with `sha256sum`

## Step 2. Run nf-core/rnaseq pipeline on Randi

***Use `tmux` to avoid accidental terminal disconnection.***

<span style="color:blue">***A demo can be found on Randi at /gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/demo. You can copy the scripts from demo to your working directory and start from there:***</span>

```
cd /path/to/your/working/dir
cp /gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/demo/run.sh .
cp /gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/demo/nextflow.config .
```

`run.sh` contains the following nextflow command:

```
nextflow run nf-core/rnaseq -r 3.16.0 \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCm39M27 \
    --gencode \
    --extra_salmon_quant_args="--gcBias" \
    -resume
```

To prepare the `samplesheet.csv` file, please check [https://nf-co.re/rnaseq/3.16.0/docs/usage/](https://nf-co.re/rnaseq/3.16.0/docs/usage/) for the format requirements.

We use Gencode genome annotation references. The reference genome parameter set by `--genome` can be customized in the `nextflow.config` file. *Please check for the latest release each time you run the pipeline.*

The `--gcBias` option is recommended for estimating correction factors to address systematic biases commonly found in RNA-seq data (Love, Hogenesch, & Irizarry, 2016; Patro et al., 2017), unless you are confident that your data are free of such biases.

Once all parameters and configurations are set, run the pipeline for QC and read mapping:

```
sh run.sh
```

<span style="color:blue">***So far, the downstream DE analysis is only compatible with the default `STAR-Salmon` route in nf-core/rnaseq pipeline. We will make it more flexible in the future.***</span>

## Step 3. Generate the DE analysis report on Randi

Before we run the DE analysis program, please make sure the following files are available in your working directory:
1) A metadata table containing sample information in TXT format. 
2) A `start.slurm` script to submit for running jobs.
3) The nf-core/rnaseq output folder.

### Metadata

The metadata table must meet the following requirements: 
- The first column should match sample names used in nf-core/rnaseq.
- A column named "group" must be included to define the sample groups for comparison.
- A column named "batch" must be included to define the batch information if batch correction is needed.

Here is an example:
|sample|timepoint|treatment|group|
|------|---------|---------|-----|
|RW06|D3|Control|D3.C|
|RW05|D3|Control|D3.C|
|RW07|D3|Control|D3.C|
|RW13|D5|IR|D5.IR|
|RW14|D5|IR|D5.IR|
|RW15|D5|IR|D5.IR|
|RW16|D5|Control|D5.C|
|RW17|D5|Control|D5.C|
|RW18|D5|Control|D5.C|

### Slurm job

Copy the `start.slurm` script from demo to your working directory: 
```
cp /gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/demo/start.slurm .
```

Here is the content of the script:
```
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

multiqc results -c /gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/app/multiqc_config.yml

R -e "shiny::runApp('/gpfs/data/bioinformatics/Pipelines/Production/Internal/CRI-BulkRNAseq-Pipeline/app', host = '0.0.0.0', port = 3838)"
```

The process starts by activating the Conda environment with all required dependencies installed. Next, it runs MultiQC to regenerate a customized QC report which is better than the default one from nf-core/rnaseq. Finally, the Shiny app is launched, providing an interface for configuring parameters.

To submit the job:
```
sbatch start.slurm
```

*If the job fails due to an unavailable port, try changing the default port (3838) to an alternative one.*

---

To access the interface in your **local** browser, run the following command in your **local** terminal:
```
ssh -N -f -L 8080:<running_node>:3838 <username>@randi.cri.uchicago.edu
```

Replace `<username>` with your user name.

Replace `<running_node>` with node that your `start.slurm` is running on, i.e., `cri22cn096`.

Then open [http://localhost:8080/](http://localhost:8080/) in your local browser. Fill out the parameter form and submit to run the program following the instructions. The report will be generated in minutes depending on how large your dataset is. 

<span style="color:blue">***MultiQC usually takes at least 10 minutes to be done. Please wait until the multiQC report appears in your working directory and refresh [http://localhost:8080/](http://localhost:8080/). You should be able to see a webpage like below:***</span>

![app](figures/app.png)
