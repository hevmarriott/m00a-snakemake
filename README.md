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

# Run time

Statistics of core hours and memory used per samplein a pilot run.

|Tool   |       N|Core hours mean|Core hours max|Core hours stddev| memory max(MB)| memory mean(MB)|
|-------|--------|---------------|---|---|---|---|
|delly  |30.0    |14.0  |  26.5   | 3.5     |5727  |5457|
|manta  |29.0    |11.8  | 27.8   | 3.7     |468   |308|
|melt   |30.0    |9.4   | 14.5   | 2.1     |3066  |2487|
|whamg  |30.0    |3.8   | 4.4    | 0.3     |14404 |7728|
