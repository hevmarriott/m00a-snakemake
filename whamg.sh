#GRID output=/projectmine-nfs/Disk/User/kooyman/whamg/
set +xeu
conda init bash
source ~/.bashrc
conda activate /cvmfs/softdrive.nl/kooyman/sv
set -xue
/usr/bin/time samtools view -xXA -1 --write-index ${input1} -o ${sample}.bam##idx##${sample}.bam.bai

whamg_c="chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY"
threads=1

singularity exec /cvmfs/softdrive.nl/kooyman/sv/wham.sif whamg -c ${whamg_c} -x ${threads} -a /cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa -f ${sample}.bam | bgzip -c > ${sample}.wham_bad_header_bad_tags.vcf.gz
tabix -p vcf ${sample}.wham_bad_header_bad_tags.vcf.gz
rm -f ${sample}.ba*
