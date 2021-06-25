# run SV tools on grid

Warning: These scripts are tailored to run on the ProjectMinE user interface. These scripts do NOT run on another system.


## Shell scripts
There are two type of scripts:
- Scripts run the tools on workernode: delly.sh,manta.sh,melt.sh,whamg.sh.
- Script that submit all jobs jobs (submit_all.sh)

Submitting the file should be done with
`submit_all.sh sample_file.txt`
where sample_file.txt is tab separate file with sample ID with next to it the path to a cram file.

## readonly file system (/cvmfs)

To access the files under /cvmfs directory, you must mount /cvmfs/softdrive on your system

## Environment

The environment files for conda are made by:

```
conda activate /cvmfs/softdrive.nl/kooyman/sv
conda env export  > env_full.yml
conda env export --from-history > env_from_history.yml
```

Restoring the env can be done with:
```
conda env create -f env_full.yml
```
