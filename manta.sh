#GRID output=/projectmine-nfs/Disk/User/kooyman/manta/
set +xeu
conda init bash
source ~/.bashrc
conda activate /cvmfs/softdrive.nl/kooyman/sv

set -xue
samtools index ${input1}

configManta.py --runDir "manta" --reference "/cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa" --bam "${input1}"  --callRegions /cvmfs/softdrive.nl/kooyman/sv/resource/mantaregions.bed.gz
  

cd manta  
./runWorkflow.py --quiet -m local -j 1    
cd ..
cp manta/results/variants/* .
/cvmfs/softdrive.nl/kooyman/sv/share/manta-1.6.0-1/libexec/convertInversion.py /cvmfs/softdrive.nl/kooyman/sv/bin/samtools /cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa diploidSV.vcf.gz | bcftools reheader -s <(echo "${sample}") > diploidSV.vcf
bgzip -c diploidSV.vcf > ${sample}.manta.vcf.gz
tabix -p vcf ${sample}.manta.vcf.gz

rm -rf manta diploidSV.vcf*
rm ${input1}.crai
