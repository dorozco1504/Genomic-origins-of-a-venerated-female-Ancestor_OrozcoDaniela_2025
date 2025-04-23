#######################################
# General Configuration
#######################################
configfile: "config.yaml"

# Load paths and parameters from config
PATHS = config['Paths']
RAW = PATHS["data_dir"]
OUT = PATHS["output_dir"] 
ROUND = config["rounds"]
THREADS = config["threads"]

# Detect sample IDs based on FASTQ filenames
IND = glob_wildcards(f"{RAW}/Yucatan/{{ind,[^/no]+}}_1.fastq.gz").ind

#######################################
# Helper Functions to Define Final Outputs
#######################################
def final_AdapterRemoval():
    return expand(f"{OUT}/AdapterRemoval/{{round}}/{{ind}}_1-2.fastq.collapsed.gz", round=ROUND, ind=IND)

def final_bams():
    return expand(f"{OUT}/aln/{{round}}/{{ind}}.F4_q25.rmdup.uniq.bam", round=ROUND, ind=IND)

def final_stats():
    return expand(f"{OUT}/aln/{{round}}/{{ind}}.stats", round=ROUND, ind=IND)

def final_plink():
    return expand(f"{OUT}/varcalling/{{round}}/{{ind}}.bim", round=ROUND, ind=IND)

#######################################
# Main Rule
#######################################
rule all:
    input:
        final_plink()

#######################################
# Adapter Removal using AdapterRemoval
#######################################
rule Adapter_Removal:
    input:
        r1 = lambda wildcards: f"{RAW}/{wildcards.round}/{wildcards.ind}_1.fastq.gz",
        r2 = lambda wildcards: f"{RAW}/{wildcards.round}/{wildcards.ind}_2.fastq.gz"
    output:
        collapsed = f"{OUT}/AdapterRemoval/{{round}}/{{ind}}_1-2.fastq.collapsed.gz",
        settings = f"{OUT}/AdapterRemoval/{{round}}/{{ind}}_1-2.fastq.settings"
    threads: THREADS["adapter_removal"]
    shell:
        """
        AdapterRemoval \
            --threads {threads} \
            --trimns \
            --trimqualities \
            --qualitybase 33 \
            --minlength 30 \
            --collapse \
            --gzip \
            --file1 {input.r1} \
            --file2 {input.r2} \
            --basename {wildcards.ind}_1-2.fastq \
            --outputcollapsed {output.collapsed} \
            --settings {output.settings} \
            --outputcollapsedtruncated /dev/null \
            --output1 /dev/null \
            --output2 /dev/null \
            --singleton /dev/null \
            --discarded /dev/null
        """

#######################################
# Alignment and Filtering
#######################################
rule bwa_aln:
    input: rules.Adapter_Removal.output.collapsed
    output: f"{OUT}/aln/{{round}}/{{ind}}.sai"
    params:
        reference="/data/users/dorozco/hackaton/data/human_genome_ref/human_g1k_v37_rCRS.fasta.gz"
    threads: 10
    shell:
        """
        /data/programs/bwa.kit/bwa aln -l 500 -t {threads} {params.reference} {input} > {output}
        """

rule bwa_samse:
    input: 
        sai = rules.bwa_aln.output,
        collapsed = rules.Adapter_Removal.output.collapsed
    output: f"{OUT}/aln/{{round}}/{{ind}}.sam.gz" 
    params:
        reference="/data/users/dorozco/hackaton/data/human_genome_ref/human_g1k_v37_rCRS.fasta.gz"
    shell:
        """
        /data/programs/bwa.kit/bwa samse -r "@RG\\tID:{wildcards.ind}\\tSM:{wildcards.ind}" {params.reference} {input.sai} {input.collapsed} | gzip -c > {output}
        """

rule samtools_bam:
    input: rules.bwa_samse.output
    output: f"{OUT}/aln/{{round}}/{{ind}}.all.bam"
    threads: 8
    shell:
        "samtools view -@ {threads} -bSh {input} > {output}"

rule filter_quality:
    input: rules.samtools_bam.output
    output: f"{OUT}/aln/{{round}}/{{ind}}.F4_q25.bam"
    shell:
        "samtools view -bh -F4 -q 25 {input} > {output}"

rule filter_unmapped:
    input: rules.samtools_bam.output
    output: f"{OUT}/aln/{{round}}/{{ind}}.unmapped.bam"
    shell:
        "samtools view -bh -f4 {input} > {output}"

rule filter_rmdup:
    input: rules.filter_quality.output
    output: 
        sort= f"{OUT}/aln/{{round}}/{{ind}}.F4_q25.sort.bam" ,
        rmdup= f"{OUT}/aln/{{round}}/{{ind}}.F4_q25.rmdup.bam"
    shell:
        """
        samtools sort -o {output.sort} {input}
        samtools rmdup -s {output.sort} {output.rmdup}
        """

rule filter_uniquemap:
    input: rules.filter_rmdup.output.rmdup
    output: f"{OUT}/aln/{{round}}/{{ind}}.F4_q25.rmdup.uniq.bam" 
    shell:
        """
        samtools index {input}
        samtools view -h {input} | grep -v 'XT:A:R' | grep -v 'XA:Z' | grep -v 'XT:A:M' | \
        awk '{{if($0~/X1:i:0/||$0~/^@/)print $0}}' | samtools view -bS - > {output}
        """

rule make_stats:
    input:
        all= rules.samtools_bam.output,
        qual_mapped= rules.filter_quality.output,
        rmdup= rules.filter_rmdup.output.rmdup,
        unique= rules.filter_uniquemap.output
    output: f"{OUT}/aln/{{round}}/{{ind}}.stats"
    shell:
        """
        echo -e ind_ID'\t'all_reads'\t'mapped_reads'\t'nodup_reads'\t'uniq_reads > {output}
        all=$(samtools view -c {input.all})
        qual_mapped=$(samtools view -c {input.qual_mapped})
        rmdup=$(samtools view -c {input.rmdup})
        uniq=$(samtools view -c {input.unique})
        echo -e {wildcards.ind}'\t'$all'\t'$qual_mapped'\t'$rmdup'\t'$uniq >> {output}
        """

#######################################
# Variant Calling and Format Conversion
#######################################
rule angsd_calling:
    input:
        bam= rules.filter_uniquemap.output,
        sites= PATHS['sites']
    output:
        haplo = f"{OUT}/varcalling/{{round}}/{{ind}}.haplo.gz",
        arg = f"{OUT}/varcalling/{{round}}/{{ind}}.arg"
    params:
        basename = f"{OUT}/varcalling/{{round}}/{{ind}}"
    threads: 10
    shell:
        "/data/programs/angsd/angsd -i {input.bam} -dohaplocall 1 -minMapQ 30 -minQ 20 -doCounts 1 -sites {input.sites} -doMajorMinor 3 -out {params.basename} -nThreads {threads}"       

rule haplotoplink:
    input: rules.angsd_calling.output.haplo
    output:
        tped = f"{OUT}/varcalling/{{round}}/{{ind}}.tped",
        tfam = f"{OUT}/varcalling/{{round}}/{{ind}}.tfam"
    params:
        basename = f"{OUT}/varcalling/{{round}}/{{ind}}"
    shell:
        "/data/programs/angsd/misc/haploToPlink {input} {params.basename}"

rule tped_bim:
    input:
        tped = rules.haplotoplink.output.tped,
        tfam = rules.haplotoplink.output.tfam
    output:
        bim = f"{OUT}/varcalling/{{round}}/{{ind}}.bim",
        fam = f"{OUT}/varcalling/{{round}}/{{ind}}.fam",
        bed = f"{OUT}/varcalling/{{round}}/{{ind}}.bed"
    params:
        basename = f"{OUT}/varcalling/{{round}}/{{ind}}"
    shell:
        """
        plink --tfile {params.basename} --make-bed --out {params.basename}
        echo ancient {wildcards.ind} 0 0 0 1 > {output.fam}
        """