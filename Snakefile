configfile: "config.yaml"


import os
import os.path
import glob
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider


# google_bucket_resources
GS = GSRemoteProvider()
GS_REFERENCE_PREFIX = "gcp-public-data--broad-references"
GS_RESOURCES_PREFIX = "gatk-sv-resources-public"

# common parameters
primary_contigs_list = GS_REFERENCE_PREFIX + "/hg38/v0/sv-resources/resources/v1/primary_contigs.list"
reference_fasta = GS_REFERENCE_PREFIX + "/hg38/v0/Homo_sapiens_assembly38.fasta"
reference_index = GS_REFERENCE_PREFIX + "/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
reference_dict = GS_REFERENCE_PREFIX + "/hg38/v0/Homo_sapiens_assembly38.dict"

# coverage inputs
preprocessed_intervals = (
    GS_RESOURCES_PREFIX
    + "/hg38/v0/sv-resources/resources/v1/preprocessed_intervals.interval_list"
)

# manta inputs
manta_region_bed = (
    GS_REFERENCE_PREFIX
    + "/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz"
)
manta_region_bed_index = (
    GS_REFERENCE_PREFIX
    + "/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz.tbi"
)

# PESR inputs 
sd_locs_vcf = (
    GS_REFERENCE_PREFIX
    + "/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
)

# MELT inputs
melt_standard_vcf_header = (
    GS_RESOURCES_PREFIX
    + "/hg38/v0/sv-resources/resources/v1/melt_standard_vcf_header.txt"
) 

#melt_metrics_intervals = config[]
#insert_size = config[]
#read_length = config[]
#coverage = config[]
#metrics_intervals = config[]
#pct_chimeras = config[]
#total_reads = config[]
#pf_reads_improper_pairs = config[]

# WHAM inputs
wham_include_list_bed_file = (
    GS_REFERENCE_PREFIX 
    + "/hg38/v0/sv-resources/resources/v1/wham_whitelist.bed"
) 

# module metrics parameters 
primary_contigs_fai = GS_REFERENCE_PREFIX + "/hg38/v0/sv-resources/resources/v1/contig.fai"
# NEED THE BASELINE VCFs - OPTIONAL FOR METRICS 

# checking_directory_structure
cram_dir = config["cram_dir"]
out_dir = config["out_dir"]
melt_dir = config["melt_dir"]

# checking sample names and cram files can be found and read
sample_name = config["sample_name"]
cram_files = expand(cram_dir + "{sample}.cram", sample=sample_name)


# creating output_directories
bam_dir = out_dir + "bam/"
counts_dir = out_dir + "counts/"
pesr_dir = out_dir + "pesr/"
whamg_dir = out_dir + "whamg/"
multiple_metrics_dir = out_dir + "multiple_metrics/"
final_melt_dir = out_dir + "melt"
module00cgvcf_dir = out_dir + "module00c_gvcf/"


# pseudo-rule to collect all target files
rule bam:
    input:
        expand(bam_dir + "{sample}.bam", sample=sample_name) + expand(
            bam_dir + "{sample}.bam.bai", sample=sample_name
        ),


rule counts:
    input:
        expand(counts_dir + "{sample}.counts.tsv.gz", sample=sample_name) + expand(
            counts_dir + "condensed_counts.{sample}.tsv.gz", sample=sample_name
        ) + expand(
            counts_dir + "{sample}.counts_intervals.interval_list", sample=sample_name
        ),


rule pesr:
    input:
        expand(pesr_dir + "{sample}.pe.txt.gz", sample=sample_name) + expand(
            pesr_dir + "{sample}.pe.txt.gz.tbi", sample=sample_name
        ) + expand(
            pesr_dir + "{sample}.sr.txt.gz", sample=sample_name
        ) + expand(
            pesr_dir + "{sample}.sr.txt.gz.tbi", sample=sample_name
        ) + expand(
            pesr_dir + "{sample}.sd.txt.gz", sample=sample_name
        ) + expand(
            pesr_dir + "{sample}.sd.txt.gz.tbi", sample=sample_name),


rule variants:
    input:
        expand(out_dir + "manta/{sample}.manta.vcf.gz", sample=sample_name) + expand(
            out_dir + "manta/{sample}.manta.vcf.gz.tbi", sample=sample_name
        ) + expand(
            whamg_dir + "{sample}.wham.vcf.gz", sample=sample_name
        ) + expand(
            whamg_dir + "{sample}.wham.vcf.gz.tbi", sample=sample_name
        ) + expand(
            final_melt_dir + "/{sample}/{sample}.melt_fix.vcf.gz", sample=sample_name
        ) + expand(
            final_melt_dir + "/{sample}/{sample}.melt_fix.vcf.gz.tbi",
            sample=sample_name,
        ),

rule fixvariants:
    input:
        expand(final_melt_dir + "/{sample}.melt.vcf.gz", sample=sample_name) + expand(
            final_melt_dir + "/{sample}.melt.vcf.gz.tbi", sample=sample_name
        )

rule haplotype:
    input:
        expand(module00cgvcf_dir + "{sample}.g.vcf.gz", sample=sample_name) + expand(
            module00cgvcf_dir + "{sample}.g.vcf.gz.tbi", sample=sample_name
        ),

# need to create a rule for scramble, do module metrics module and rule to include bams without cram conversion
            
rule CramToBam:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        cram_file=cram_dir + "{sample}.final-gatk.cram",
    output:
        bam_file=bam_dir + "{sample}.bam",
        bam_index=bam_dir + "{sample}.bam.bai",
    benchmark:
        "benchmarks/CramToBam/{sample}.tsv"
    conda:
        "envs/samtools.yaml"
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        """
        samtools view -b -h -@ {threads} -T {input[0]} -o {output.bam_file} {input.cram_file}
        samtools index -@ {threads} {output.bam_file}
        """

rule CollectCounts:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        GS.remote(reference_dict, keep_local=True),
        GS.remote(preprocessed_intervals, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        bam_index=rules.CramToBam.output.bam_index,
    output:
        counts_file=counts_dir + "{sample}.counts.tsv.gz",
    benchmark:
        "benchmarks/CollectCounts/{sample}.tsv"
    conda:
        "envs/gatk.yaml"
    params:
        sample="{sample}",
        temp=counts_dir + "{sample}.counts.tsv",
    threads: 1
    resources:
        mem_mb=12000,
    shell:
        """
        gatk --java-options "-Xmx10024m" CollectReadCounts --input {input.bam_file} --read-index {input.bam_index} --reference {input[0]} -L {input[3]} --format TSV --interval-merging-rule OVERLAPPING_ONLY -O {params.temp} --disable-read-filter MappingQualityReadFilter 
        sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:{params.sample}/g" {params.temp}
        bgzip {params.temp}
        """

rule CondenseReadCounts:
    input:
        counts = rules.CollectCounts.output.counts_file
    output: 
        condensed_counts = counts_dir + "condensed_counts.{sample}.tsv.gz"
    benchmark:
        "benchmarks/CondenseReadCounts/{sample}.tsv"
    params:
        sample="{sample}",
        temp_in_rd=DNAscan_results_dir + "counts/{sample}.in.rd.txt.gz",
        temp_out_rd=DNAscan_results_dir + "counts/{sample}.out.rd.txt.gz",
        temp_ref_dict=DNAscan_results_dir + "counts/{sample}.ref.dict"
    threads: 1
    resources:
        mem_mb = 3000
    conda:
        "envs/gatk.yaml"
    shell:
        """
        zcat {input.counts} | grep '^@' | grep -v '@RG' > {params.temp_ref_dict}
        zcat {input.counts} | grep -v '^@' | sed -e 1d | awk 'BEGIN{{FS=OFS="\t";print "#Chr\tStart\tEnd\tNA21133"}}{{print $1,$2-1,$3,$4}}' | bgzip > {params.temp_in_rd}
        tabix -0 -s1 -b2 -e3 {params.temp_in_rd}
        gatk --java-options -Xmx2g CondenseDepthEvidence -F {params.temp_in_rd} -O {params.temp_out_rd} --sequence-dictionary {params.temp_ref_dict} --max-interval-size 2000 --min-interval-size 101
        cat {params.temp_ref_dict} <(zcat {params.temp_out_rd} | awk 'BEGIN{{FS=OFS="\t";print "@RG\tID:GATKCopyNumber\tSM:{params.sample}\\nCONTIG\tSTART\tEND\tCOUNT"}}{{if(NR>1)print $1,$2+1,$3,$4}}') | bgzip > {output.condensed_counts}
        rm {params.temp_ref_dict}
        rm {params.temp_in_rd}*
        rm {params.temp_out_rd}*
        """

rule CountsToIntervals:
    input:
        counts=rules.CollectCounts.output.counts_file,
    output:
        counts_intervals=counts_dir + "{sample}.interval_list",
    benchmark:
        "benchmarks/CountsToIntervals/{sample}.tsv"
    shell:
        """
        zgrep "^@" {input.counts} > {output.counts_intervals}
        zgrep -v "^@" {input.counts} | sed -e 1d | awk -F "\t" -v OFS="\t" '{{print $1,$2,$3,"+","."}}' >> {output.counts_intervals}
        """

rule PESRCollection:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_dict, keep_local=True),
        GS.remote(sd_locs_vcf, keep_local=True),
        GS.remote(preprocessed_intervals, keep_local=True),
        GS.remote(primary_contigs_list, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
    output:
        PE_file=pesr_dir + "{sample}.pe.txt.gz",
        PE_file_index=pesr_dir + "{sample}.pe.txt.gz.tbi",
        SR_file=pesr_dir + "{sample}.sr.txt.gz",
        SR_file_index=pesr_dir + "{sample}.sr.txt.gz.tbi",
        SD_file=pesr_dir + "{sample}.sd.txt.gz",
        SD_index=pesr_dir + "{sample}.sd.txt.gz.tbi",
    benchmark:
        "benchmarks/PESRCollection/{sample}.tsv"
    conda:
        "envs/gatk.yaml"
    resources:
        mem_mb=4000,
    params:
        sample="{sample}",
    shell:
        """
        gatk --java-options "-Xmx3250m" CollectSVEvidence -I {input.bam_file} --sample-name {params.sample} -F {input[2]} -SR {output.SR_file} -PE {output.PE_file} -SD {output.SD_file} --site-depth-min-mapq 6 --site-depth-min-baseq 10 -R {input[0]} -L {input[4]}
        """

rule runMantastep1:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        GS.remote(manta_region_bed, keep_local=True),
        GS.remote(manta_region_bed_index, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        bam_index=rules.CramToBam.output.bam_index,
    output:
        diploid_vcf=out_dir + "manta/{sample}/results/variants/diploidSV.vcf.gz"
    benchmark:
        "benchmarks/runMantastep1/{sample}.tsv"
    singularity:
        "docker://hevmarriott/manta:v1.6"
    threads: 8
    resources:
        mem_mb=16000,
    params:
        memGb=16,
        sample="{sample}",
        manta_run_dir=out_dir + "manta/{sample}",
    shell:
        """
        /usr/local/bin/manta/bin/configManta.py --bam {input.bam_file} --referenceFasta {input[0]} --runDir {params.manta_run_dir} --callRegions {input[2]}
        {params.manta_run_dir}/runWorkflow.py --mode local --jobs {threads} --memGb {params.memGb}
        python2 /usr/local/bin/manta/libexec/convertInversion.py /usr/local/bin/samtools {input[0]} {params.manta_run_dir}/results/variants/diploidSV.vcf.gz | bcftools reheader -s <(echo "{params.sample}") > {params.manta_run_dir}/results/variants/diploidSV.vcf
        bgzip -c {params.manta_run_dir}/results/variants/diploidSV.vcf > {output.manta_vcf}
        tabix -p vcf {output.manta_vcf}
        """

rule runMantastep2:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        rules.runMantastep1.output.diploid_vcf,
    output:
        manta_vcf=out_dir + "manta/{sample}.manta.vcf.gz",
        manta_index=out_dir + "manta/{sample}.manta.vcf.gz.tbi",
    benchmark:
        "benchmarks/runMantastep2/{sample}.tsv"
    singularity:
        "docker://hevmarriott/manta:v1.6"
    threads: 8
    resources:
        mem_mb=16000,
    params:
        sample="{sample}",
        manta_run_dir=out_dir + "manta/{sample}",
    shell:
        """
        python2 /usr/local/bin/manta/libexec/convertInversion.py /usr/local/bin/samtools {input[0]} {input[1]} | bcftools reheader -s <(echo "{params.sample}") > {params.manta_run_dir}/results/variants/diploidSV.vcf
        bgzip -c {params.manta_run_dir}/results/variants/diploidSV.vcf > {output.manta_vcf}
        tabix -p vcf {output.manta_vcf}
        """

rule runWhamg:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        GS.remote(wham_include_list_bed_file, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        bam_index=rules.CramToBam.output.bam_index,
    output:
        whamg_vcf=whamg_dir + "{sample}.wham.vcf.gz",
        whamg_index=whamg_dir + "{sample}.wham.vcf.gz.tbi"
    benchmark:
        "benchmarks/runWhamg/{sample}.tsv"
    threads: 8
    resources:
        mem_mb=32000,
    singularity:
        "docker://gatksv/wham:8645aa"
    params:
        sample = "{sample}",
        whamg_c="chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY",
    shell:
        """
        mkdir {whamg_dir}/{params.sample}_temp
        cd {whamg_dir}/{params.sample}_temp
        awk 'BEGIN{FS=OFS="\t"}{printf("%07d\t%s\n",NR,$1":"$2"-"$3)}' {input[2]} |\
          while read -r line interval; do
            vcfFile="$line.wham.vcf.gz"
            whamg -c "{params.whamg_c} -x {threads} -a {input[0]} -f {input.bam_file} -r $interval | bgzip -c > $vcfFile
            bcftools index -t $vcfFile
          done

        ls -1 *.wham.vcf.gz > vcf.list
        bcftools concat -a -Ov -f vcf.list | sed -e 's/^#CHROM\t.*/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{params.sample}/' -e 's/;TAGS=[^;]*;/;TAGS={params.sample};/' | bgzip -c > {output.whamg_vcf}
        cd ..
        bcftools index -t {output.whamg_vcf}

        df -h
        ls -l
        """

#GOT TO HERE - THINK I MIGHT NEED COLLECT COVERAGE FILES
#coverage = value of MEAN_COVERAGE from CollectWgsMetrics
#read_length = MEAN_READ_LENGTH from CollectAlignmentSummaryMetrics
#insert_size = MEAN_INSERT_SIZE from CollectInsertSizeMetrics

rule runMELTInputMetrics:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
    output:
        alignment_file = multiple_metrics_dir + "{sample}.alignment_summary_metrics",
        insert_file = multiple_metrics_dir + "{sample}.insert_size_metrics",
        quality_file = multiple_metrics_dir + "{sample}.quality_distribution_metrics",
        wgs_metrics_file = multiple_metrics_dir + "{sample}_wgs_metrics.txt"
    benchmark:
        "benchmarks/runMELTInputMetrics/{sample}.tsv"
    resources:
        mem_mb=4000, 
    conda:
        "envs/gatk.yaml"
    params:
        sample = "{sample}",
        metrics_base = out_dir + "multiple_metrics"
        read_length=config["mean_read_length"],
    shell:
        """
        gatk --java-options -Xmx3250m CollectMultipleMetrics -I {input.bam_file} -O {params.metrics_base}/{params.sample} -R {input[0]} --ASSUME_SORTED true --PROGRAM null --PROGRAM CollectAlignmentSummaryMetrics \
        --PROGRAM CollectInsertSizeMetrics --PROGRAM CollectSequencingArtifactMetrics --PROGRAM CollectGcBiasMetrics --PROGRAM QualityScoreDistribution --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE
        gatk --java-options -Xmx3250m CollectWgsMetrics --INPUT {input.bam_file} --VALIDATION_STRINGENCY SILENT --REFERENCE_SEQUENCE {input[0]} --READ_LENGTH {params.mean_read_length} --INCLUDE_BQ_HISTOGRAM true --OUTPUT {params.metrics_base}_{params.sample}_wgs_metrics.txt --USE_FAST_ALGORITHM true
        """

rule runMELT:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        melt_bed_file=melt_dir + "/add_bed_files/Hg38/Hg38.genes.bed",
    output:
        melt_fix_vcf=final_melt_dir + "/{sample}/{sample}.melt_fix.vcf.gz",
        melt_fix_index=final_melt_dir + "/{sample}/{sample}.melt_fix.vcf.gz.tbi",
    benchmark:
        "benchmarks/runMELT/{sample}.tsv"
    resources:
        mem_mb=32000,
    conda:
        "envs/manta_melt.yaml"
    params:
        coverage=config["mean_coverage"],
        read_length=config["mean_read_length"],
        insert_size=config["mean_insert_size"],
        mean_chrom_length=40000000,
        sample="{sample}",
        melt_results_dir=final_melt_dir + "/{sample}",
        vcf_sort="resources/vcf-sort.pl",
    shell:
        """
        ls {melt_dir}/me_refs/Hg38/*zip | sed 's/\*//g' > {params.melt_results_dir}/transposon_reference_list         
        java -Xmx12000m -jar {melt_dir}/MELT.jar Single -bamfile {input.bam_file} -h {input[0]} -c {params.coverage} -r {params.read_length} -e {params.insert_size} -d {params.mean_chrom_length} -t {params.melt_results_dir}/transposon_reference_list -n {input.melt_bed_file} -w {params.melt_results_dir}
        cat {params.melt_results_dir}/SVA.final_comp.vcf | grep "^#" > "{params.melt_results_dir}/{params.sample}.header.txt"
        cat {params.melt_results_dir}/SVA.final_comp.vcf | grep -v "^#" > "{params.melt_results_dir}/{params.sample}.sva.vcf"
        cat {params.melt_results_dir}/LINE1.final_comp.vcf | grep -v "^#" > "{params.melt_results_dir}/{params.sample}.line1.vcf"
        cat {params.melt_results_dir}/ALU.final_comp.vcf | grep -v "^#" > "{params.melt_results_dir}/{params.sample}.alu.vcf"
        cat {params.melt_results_dir}/{params.sample}.header.txt {params.melt_results_dir}/{params.sample}.sva.vcf {params.melt_results_dir}/{params.sample}.line1.vcf {params.melt_results_dir}/{params.sample}.alu.vcf | perl {params.vcf_sort} -c | bgzip -c > {output.melt_fix_vcf}
        tabix -p vcf {output.melt_fix_vcf} 
        """


rule MELTFixOutput:
    input:
        melt_vcf_header="resources/melt_standard_vcf_header.txt",
        bam_file=rules.CramToBam.output.bam_file,
        melt_fix_vcf=rules.runMELT.output.melt_fix_vcf,
    output:
        melt_vcf=final_melt_dir + "/{sample}.melt.vcf.gz",
        melt_index=final_melt_dir + "/{sample}.melt.vcf.gz.tbi",
    params:
        sample="{sample}",
        melt_results_dir=final_melt_dir + "/{sample}",
        temp_header="/{sample}_temp_header.txt",
        temp_vcf="/{sample}.temp_melt_vcf.gz",
    benchmark:
        "benchmarks/MELTFixOutput/{sample}.tsv"
    conda:
        "envs/manta_melt.yaml"
    shell:
        """
        vcf_text=$(bgzip -cd {input.melt_fix_vcf})
        grep '^#CHR' <<<"$vcf_text" | sed 's|{{basename({input.bam_file}, ".bam")}}|{params.sample}|g' > {params.melt_results_dir}/{params.temp_header}
        grep -v '^#' <<<"$vcf_text" | sed 's/No Difference/No_Difference/g' >> {params.melt_results_dir}/{params.temp_header}
        FILEDATE=$(grep -F 'fileDate=' <<<"$vcf_text")
        cat "{input.melt_vcf_header}" {params.melt_results_dir}/{params.temp_header} | sed "2i$FILEDATE" | bgzip -c > {input.melt_fix_vcf}
        mv {input.melt_fix_vcf} {params.melt_results_dir}/{params.temp_vcf}
        bcftools reheader -s <( echo "{params.sample}") {params.melt_results_dir}/{params.temp_vcf} > {output.melt_vcf}  
        rm {params.melt_results_dir}/{params.temp_vcf}
        tabix -p vcf {output.melt_vcf}
        """

rule Module00cGVCF:
    input:
        GS.remote(reference_fasta, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        interval_file=rules.CountsToIntervals.output.counts_intervals,
    output:
        GVCF=module00cgvcf_dir + "{sample}.g.vcf.gz",
        GVCF_index=module00cgvcf_dir + "{sample}.g.vcf.gz.tbi",
    benchmark:
        "benchmarks/Module00cGVCF/{sample}.tsv"
    conda:
        "envs/gatk.yaml"
    resources:
        mem_mb=16000,
    threads: 16
    shell:
        """
        gatk --java-options "-Xmx16G" HaplotypeCaller -R {input[0]} -I {input.bam_file} -L {input.interval_file} -O {output.GVCF} -ERC GVCF --native-pair-hmm-threads {threads} --interval-padding 100 
        """

