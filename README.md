# m00a-snakemake
Snakemake workflow for executing the entirety of GATK-SV Module00a on a local Slurm HPC facility as well as generating GVCF files with Haplotype Caller for input into Module00c. 

# Table of Contents
1. [Obtaining](#obtaining)
2. [Requirements](#requirements)
   * [How to make a Google Cloud account](#how-to-make-a-Google-Cloud-service-account)
   * [How to create a Snakemake job execution profile](#how-to-create-a-Snakemake-job-execution-profile)
3. [Usage](#usage)
   * [Sub-Module Descriptions](#sub-module-descriptions)
   * [Steps after running each module](#steps-after-running-each-module)
4. [DAG of Workflow Modules with Subcommand Memory Allocation](#dag-of-workflow-modules-with-subcommand-memory-allocation) 
   * [counts](#counts)
   * [variants](#variants)
   * [fixvariants](#fixvariants)
   * [haplotype](#haplotype)

## Obtaining
Clone this repository using Git:
``` git clone https://github.com/hevmarriott/m00a-snakemake```

## Requirements
- Snakemake installed via Conda in a standalone environment - [link to recommended installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
   * NOTE - this environment takes up 1.1 Gigabytes of disk storage
- A Google Cloud service account to access the public GATK-SV resource buckets and a dedicated Google Cloud Storage Bucket
   * NOTE - with Cromwell, underscores in the storage bucket name are not permitted
- Snakemake profile for submitting many parallel jobs via the SLURM cluster
   * NOTE - this should be created in the m00a-snakemake/ directory

### How to make a Google Cloud service account
If you do not already have a service account, navigate to [this page](https://cloud.google.com/iam/docs/creating-managing-service-accounts) for a detailed walkthrough on how to create an account using the gcloud command-line tool, which needs to be installed beforehand. After following all of the instructions, you can check if the account has been credentialed with:
```
gcloud auth list
``` 
If you have multiple accounts (i.e. your own personal gmail account and service account) and the wrong one is active, you can activate the service account by running:
```
gcloud config set account `NAME_OF_SERVICE_ACCOUNT@PROJECT_NAME.iam.gserviceaccount.com`
```

Next, ensure that [billing is enabled](https://cloud.google.com/billing/docs/how-to/modify-project) for your project that the service account is linked to, and that the Cloud Life Sciences, Compute Engine and Cloud Storage APIs are enabled.

### How to create a Snakemake job execution profile
First, cookiecutter needs to be installed via conda in the activated snakemake environment:
```
conda install -c conda-forge cookiecutter
```
Then, follow the instructions [here](https://github.com/Snakemake-Profiles/slurm) using either Example 1 or Example 2 depending on your needs. The result will be a SLURM profile directory (slurm.account_name) in the m00a-snakemake directory with three different files - slurm-submit.py, slurm-jobscript.sh, slurm-status.py.

## Usage
Before running the workflow, make sure that you are in the m00a-snakemake directory and have logged into your Google cloud account by running:
``` gcloud auth application-default login```

NOTE - you only have to do this for the first pass of the entire pipeline, as after that, the files stay in the local directory as gatk-sv-resources-public/ and gcp-public-data--broad-references/

Next, go to the config.yaml file and enter the locations for the input_data, out and MELT directories, sample names and MELT parameters. 

Then, invoke the pipeline using the following code:
``` 
snakemake *COMMAND* --profile *NAME_OF_SNAKEMAKE_PROFILE_DIR* --use-conda --use-singularity -j *NUMBER_OF_PARALLEL_JOBS_TO_SUBMIT*
```
Where ```--use singularity``` is only relevant for the *variants* command as Manta and Wham docker containers are used. 

### Sub-Module Descriptions
Below are descriptions for each command in the workflow, with the parallel job allocation and approximate storage usage listed when running the pipeline with 10 samples at a time. This will enable you to easily scale the workflow in accordance to any storage and job scheduler limits. 

| Command | Description | Job allocation with 10 Samples | Storage usage with 10 Samples | Clock Time with 10 Samples | 
| --- | --- | --- | --- | --- | 
| bam | Converts cram input to bam format necessary for functioning of the entire workflow | -j 10 | 550 GB | 1 hour | 
| counts | Performs split-read and paired-read evidence collection on the bam files as well as collecting the binned read counts and generating count interval lists | -j 20 | 12 GB | 1 hour
| variants | Runs Delly, Manta, MELT and Whamg variant callers on the bam files | -j 40 | 46.5 MB | 18-21 hours
| fixvariants | Converts Delly to VCF output and reformats Whamg and Melt headers | -j 30 | N/A | < 1 minute
| haplotype | Runs GATK Haplotype Caller on the bam files to obtain GVCFs which can be used to generate B allele frequency in Module00c | -j 10 | 120 GB | 24-30 hours

### After running each sub-module
You will find that after passing the first batch of samples through the pipeline, new subdirectories appear in the out_directory location in addition to those in the results directory. These include the aforementioned gatk-sv resource folders which are there for future running of the workflow without need to authenticate, but also a benchmarks folder. Within this, there is a folder for each step of the workflow, and subsequent files for each sample which documents the wall and CPU time of the job as well as various memory usage  and request metrics. More importantly, these are only created after the job has completed, which means that the creation of files is also one way to measure if the files for each sample are fully complete. 

When you have all the appropriate files, you then need to upload them to Google Cloud Storage in order to be able downstream GATK-SV steps. This can be gone via the gsutil command line tool. 

If you do not have a bucket, you can create one under your service account by performing:
``` 
gsutil mb gs://my-bucket 
```
Once the bucket is created, you can then move the files i.e. delly.vcf.gz and index to your bucket:
```
cd out_dir/results/
gsutil mv delly/* gs://my-bucket
```

## DAG of Workflow Modules with Subcommand Memory Allocation
The workflows below specify the input dependencies and output of each command, as well as the memory and thread allocation in the Snakefile for each subcommand. The advantage to using snakemake is that if a file did not generate output due to a failure with the cluster etc., then by invoking the command again, Snakemake will automatically decipher on which remaining files need to be generated, without duplicating outputs.

### counts
![alt text](https://github.com/hevmarriott/m00a-snakemake/blob/main/counts.PNG)

* CramToBam: 16GB, 4 threads
* CollectCounts: 12GB, 1 thread
* CountsToIntervals: 1GB, 1 thread
* PESRCollection: 4GB, 1 thread

### variants
![alttext](https://github.com/hevmarriott/m00a-snakemake/blob/main/variants.PNG)

* runDelly: 16GB, 1 thread
* runManta: 4GB, 8 threads
* runMELT: 32GB, 1 thread
* runWhamg: 32GB, 8 threads
### fixvariants
![alttext](https://github.com/hevmarriott/m00a-snakemake/blob/main/fixvariants.PNG)

* Dellybcf2vcf: 4GB, 1 thread
* MELTFixOutput: 1GB, 1 thread
* WhamgOutput: 4GB, 4 threads
* WhamgFixOutput: 1GB, 1 thread

### haplotype
![alttext](https://github.com/hevmarriott/m00a-snakemake/blob/main/haplotype.PNG)

* Module00cGVCF: 16GB, 16 threads
