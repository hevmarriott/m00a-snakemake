#GRID output=/projectmine-nfs/Disk/User/kooyman/melt/
set +xeu
conda init bash
source ~/.bashrc
conda activate /cvmfs/softdrive.nl/kooyman/sv

set -xue
samtools index ${input1}

melt_dir="/cvmfs/softdrive.nl/kooyman/sv/resource/MELTv2.2.2" 
mean_chrom_length=40000000
vcfsort="/cvmfs/softdrive.nl/kooyman/sv/resource/vcf-sort.pl"

insertsize=$(samtools stats ${input1}  chr1:200000-400000 |grep ^SN |grep "insert size average:"|grep -oP "[0-9.]+"|cut -f 1 -d ".")
readlength=$(samtools stats ${input1}  chr1:200000-400000 |grep ^SN |grep -P "average length"|grep -oP "[0-9.]+"|cut -f 1 -d ".")


mkdir -p result
ls ${melt_dir}/me_refs/Hg38/*zip | sed 's/\*//g' > result/transposon_reference_list         
java -Xmx12000m -jar ${melt_dir}/MELT.jar Single \
-bamfile ${input1} \
-h /cvmfs/softdrive.nl/projectmine_sw/resources/Build38/hs38DH/hs38DH.fa \
-r $readlength \
-e $insertsize \
-d $mean_chrom_length \
-t result/transposon_reference_list \
-n ${melt_dir}/add_bed_files/Hg38/Hg38.genes.bed \
-w result

cat result/SVA.final_comp.vcf | grep "^#" > "result/${sample}.header.txt"
cat result/SVA.final_comp.vcf | grep -v "^#" > "result/${sample}.sva.vcf"
cat result/LINE1.final_comp.vcf | grep -v "^#" > "result/${sample}.line1.vcf"
cat result/ALU.final_comp.vcf | grep -v "^#" > "result/${sample}.alu.vcf"
cat result/${sample}.header.txt result/${sample}.sva.vcf result/${sample}.line1.vcf result/${sample}.alu.vcf | perl $vcfsort -c | bgzip -c > ${sample}.melt_fix.vcf.gz
tabix -p vcf ${sample}.melt_fix.vcf.gz

rm -rf result
rm ${input1}.crai *.disc.bai *.disc *.fq

