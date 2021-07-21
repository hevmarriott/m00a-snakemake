#GRID output=/projectmine-nfs/Disk/User/kooyman/delly/
set +xeu
conda init bash
source ~/.bashrc
conda activate /cvmfs/softdrive.nl/kooyman/sv
set -xue
samtools index ${input1}

delly call -g /cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa -o ${sample}.delly.bcf -x /cvmfs/softdrive.nl/kooyman/sv/resource/delly_human.hg38.excl.tsv ${input1}
rm ${input1}.crai
