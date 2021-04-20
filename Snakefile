configfile: "config.yaml"


import os
import os.path
import glob
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

# google_bucket_resources
GS = GSRemoteProvider()
GS_REFERENCE_PREFIX = "gcp-public-data--broad-references"
GS_RESOURCES_PREFIX = "gatk-sv-resources-public"

reference_fasta = GS_REFERENCE_PREFIX + "/hg38/v0/Homo_sapiens_assembly38.fasta"
reference_index = GS_REFERENCE_PREFIX + "/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
reference_dict = GS_REFERENCE_PREFIX + "/hg38/v0/Homo_sapiens_assembly38.dict"
preprocessed_intervals = (
    GS_RESOURCES_PREFIX
    + "/hg38/v0/sv-resources/resources/v1/preprocessed_intervals.interval_list"
)
delly_exclude_intervals = (
    GS_RESOURCES_PREFIX + "/hg38/v0/sv-resources/resources/v1/delly_human.hg38.excl.tsv"
)
manta_region_bed = (
    GS_REFERENCE_PREFIX
    + "/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz"
)
manta_region_bed_index = (
    GS_REFERENCE_PREFIX
    + "/hg38/v0/sv-resources/resources/v1/primary_contigs_plus_mito.bed.gz.tbi"
)
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
pesr_dir = out_dir + "pesr_files/"
final_delly_dir = out_dir + "delly/"
final_whamg_dir = out_dir + "whamg/"
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
            counts_dir + "{sample}.counts_intervals.interval_list", sample=sample_name
        ) + expand(pesr_dir + "{sample}.disc.txt.gz", sample=sample_name) + expand(
            pesr_dir + "{sample}.disc.txt.gz.tbi", sample=sample_name
        ) + expand(
            pesr_dir + "{sample}.split.txt.gz", sample=sample_name
        ) + expand(
            pesr_dir + "{sample}.split.txt.gz.tbi", sample=sample_name
        ),


rule variants:
    input:
        expand(out_dir + "manta/{sample}.manta.vcf.gz", sample=sample_name) + expand(
            out_dir + "manta/{sample}.manta.vcf.gz.tbi", sample=sample_name
        ) + expand(
            final_delly_dir + "{sample}/{sample}.delly.bcf", sample=sample_name
        ) + expand(
            final_whamg_dir + "{sample}/{sample}.wham_bad_header_bad_tags.vcf.gz",
            sample=sample_name,
        ) + expand(
            final_whamg_dir + "{sample}/{sample}.wham_bad_header_bad_tags.vcf.gz.tbi",
            sample=sample_name,
        ) + expand(
            final_melt_dir + "/{sample}/{sample}.melt_fix.vcf.gz", sample=sample_name
        ) + expand(
            final_melt_dir + "/{sample}/{sample}.melt_fix.vcf.gz.tbi",
            sample=sample_name,
        ),


rule fixvariants:
    input:
        expand(final_delly_dir + "{sample}.delly.vcf.gz", sample=sample_name) + expand(
            final_delly_dir + "{sample}.delly.vcf.gz.tbi", sample=sample_name
        ) + expand(
            final_whamg_dir + "{sample}/{sample}.tags_annotation_file.tsv.gz",
            sample=sample_name,
        ) + expand(
            final_whamg_dir + "{sample}/{sample}.tags_annotation_file.tsv.gz.tbi",
            sample=sample_name,
        ) + expand(
            final_whamg_dir + "{sample}/{sample}.wham_bad_header.vcf.gz",
            sample=sample_name,
        ) + expand(
            final_whamg_dir + "{sample}/{sample}.wham_bad_header.vcf.gz.tbi",
            sample=sample_name,
        ) + expand(
            final_whamg_dir + "{sample}.wham.vcf.gz", sample=sample_name
        ) + expand(
            final_whamg_dir + "{sample}.wham.vcf.gz.tbi", sample=sample_name
        ) + expand(
            final_melt_dir + "/{sample}.melt.vcf.gz", sample=sample_name
        ) + expand(
            final_melt_dir + "/{sample}.melt.vcf.gz.tbi", sample=sample_name
        ),


rule haplotype:
    input:
        expand(module00cgvcf_dir + "{sample}.g.vcf.gz", sample=sample_name) + expand(
            module00cgvcf_dir + "{sample}.g.vcf.gz.tbi", sample=sample_name
        ),
            
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
        gatk --java-options "-Xmx10024m" CollectReadCounts -I {input.bam_file} --read-index {input.bam_index} -R {input[0]} -L {input[3]} --format TSV --interval-merging-rule OVERLAPPING_ONLY -O {params.temp}
        sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:{params.sample}/g" {params.temp}
        bgzip {params.temp}
        """


rule CountsToIntervals:
    input:
        counts=rules.CollectCounts.output.counts_file,
    output:
        counts_intervals=counts_dir + "{sample}.counts_intervals.interval_list",
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
        bam_file=rules.CramToBam.output.bam_file,
    output:
        PE_file=pesr_dir + "{sample}.disc.txt.gz",
        PE_file_index=pesr_dir + "{sample}.disc.txt.gz.tbi",
        SR_file=pesr_dir + "{sample}.split.txt.gz",
        SR_file_index=pesr_dir + "{sample}.split.txt.gz.tbi",
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
        gatk --java-options "-Xmx3250m" PairedEndAndSplitReadEvidenceCollection -I {input.bam_file} --pe-file {output.PE_file} --sr-file {output.SR_file} --sample-name {params.sample} -R {input[0]}
        tabix -f -s1 -b 2 -e 2 {output.PE_file}
        tabix -f -s1 -b 2 -e 2 {output.SR_file}
        """


rule runDelly:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        GS.remote(delly_exclude_intervals, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        bam_index=rules.CramToBam.output.bam_index,
    output:
        delly_bcf=final_delly_dir + "{sample}/{sample}.delly.bcf",
    benchmark:
        "benchmarks/runDelly/{sample}.tsv"
    conda:
        "envs/delly.yaml"
    resources:
        mem_mb=16000,
    shell:
        """
        delly call -x {input[2]} -o {output.delly_bcf} -g {input[0]} {input.bam_file}
        """


rule Dellybcf2vcf:
    input:
        delly_bcf=rules.runDelly.output.delly_bcf,
    output:
        delly_vcf=final_delly_dir + "{sample}.delly.vcf.gz",
        delly_index=final_delly_dir + "{sample}.delly.vcf.gz.tbi",
    params:
        temp=final_delly_dir + "{sample}/{sample}.delly.vcf",
        vcf_sort="resources/vcf-sort.pl",
    benchmark:
        "benchmarks/Dellybcf2vcf/{sample}.tsv"
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=4000,
    shell:
        """
        bcftools view {input.delly_bcf} > {params.temp}
        cat {params.temp} | perl {params.vcf_sort} -c | bgzip -c > {output.delly_vcf}
        tabix {output.delly_vcf}
        """


rule runManta:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        GS.remote(manta_region_bed, keep_local=True),
        GS.remote(manta_region_bed_index, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        bam_index=rules.CramToBam.output.bam_index,
    output:
        manta_vcf=out_dir + "manta/{sample}.manta.vcf.gz",
        manta_index=out_dir + "manta/{sample}.manta.vcf.gz.tbi",
    benchmark:
        "benchmarks/runManta/{sample}.tsv"
    singularity:
        "docker://gatksv/manta:8645aa"
    threads: 8
    resources:
        mem_mb=4000,
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


rule runWhamg:
    input:
        GS.remote(reference_fasta, keep_local=True),
        GS.remote(reference_index, keep_local=True),
        bam_file=rules.CramToBam.output.bam_file,
        bam_index=rules.CramToBam.output.bam_index,
    output:
        whamg_bad_tags=(
            final_whamg_dir + "{sample}/{sample}.wham_bad_header_bad_tags.vcf.gz"
        ),
        whamg_bad_tags_index=(
            final_whamg_dir + "{sample}/{sample}.wham_bad_header_bad_tags.vcf.gz.tbi"
        ),
    benchmark:
        "benchmarks/runWhamg/{sample}.tsv"
    threads: 8
    resources:
        mem_mb=32000,
    singularity:
        "docker://gatksv/wham:8645aa"
    params:
        whamg_c="chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY",
    shell:
        """
        whamg -c "{params.whamg_c}" -x {threads} -a {input[0]} -f {input.bam_file} | bgzip -c > {output.whamg_bad_tags}
        tabix -p vcf {output.whamg_bad_tags}
        """


rule WhamgOutput:
    input:
        whamg_bad_tags=rules.runWhamg.output.whamg_bad_tags,
        whamg_bad_tags_index=rules.runWhamg.output.whamg_bad_tags_index,
    output:
        annotation_file=(
            final_whamg_dir + "{sample}/{sample}.tags_annotation_file.tsv.gz"
        ),
        annotation_file_index=(
            final_whamg_dir + "{sample}/{sample}.tags_annotation_file.tsv.gz.tbi"
        ),
        whamg_bad_header=final_whamg_dir + "{sample}/{sample}.wham_bad_header.vcf.gz",
        whamg_bad_header_index=(
            final_whamg_dir + "{sample}/{sample}.wham_bad_header.vcf.gz.tbi"
        ),
    params:
        sample="{sample}",
    benchmark:
        "benchmarks/WhamgOutput/{sample}.tsv"
    conda:
        "envs/samtools.yaml"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        """
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t{params.sample}\n' {input.whamg_bad_tags} | bgzip -c > {output.annotation_file}
        tabix -f -s1 -b2 -e2 {output.annotation_file}
        bcftools annotate -a {output.annotation_file} -c CHROM,POS,REF,ALT,INFO/TAGS {input.whamg_bad_tags} | bgzip -c > {output.whamg_bad_header}
        tabix -p vcf {output.whamg_bad_header}
        """


rule WhamgFixOutput:
    input:
        GS.remote(reference_index, keep_local=True),
        whamg_bad_header=rules.WhamgOutput.output.whamg_bad_header,
    output:
        whamg_vcf=final_whamg_dir + "{sample}.wham.vcf.gz",
        wham_index=final_whamg_dir + "{sample}.wham.vcf.gz.tbi",
    params:
        sample="{sample}",
        contigs_pattern="^chr1\t|^chr2\t|^chr3\t|^chr4\t|^chr5\t|^chr6\t|^chr7\t|^chr8\t|^chr9\t|^chr10\t|^chr11\t|^chr12\t|^chr13\t|^chr14\t|^chr15\t|^chr16\t|^chr17\t|^chr18\t|^chr19\t|^chr20\t|^chr21\t|^chr22\t|^chrX\t|^chrY\t",
    benchmark:
        "benchmarks/WhamgFixOutput/{sample}.tsv"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        OLD_HEADER=$(bcftools view -h {input.whamg_bad_header} | grep -v '##contig=')
        CONTIGS_HEADER=$(grep "{params.contigs_pattern}" -P {input[0]} | awk '{{print "##contig=<ID=" $1 ",length=" $2 ">"}}')
        bcftools reheader -h <( echo "$OLD_HEADER" | sed \$d ; echo "$CONTIGS_HEADER" ; echo "$OLD_HEADER" | tail -n 1) -s <( echo "{params.sample}") {input.whamg_bad_header} > {output.whamg_vcf}
        tabix {output.whamg_vcf}
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

