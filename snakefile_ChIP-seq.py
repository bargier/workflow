workdir: '/homes/m2bsg_2018_2019/m2bsg_2019_user11/Run_231'

SAMPLE = ['S002593_Input-PF382-98098', 
	'S002594_H3K4me3-PF382-98099', 
	'S002595_Input-P12-98100', 
	'S002596_H3K4me3-P12-98101',]

CHIP = ['S002596_H3K4me3-P12-98101', 'S002594_H3K4me3-PF382-98099']

CHIP2INPUT = {'S002596_H3K4me3-P12-98101':'S002595_Input-P12-98100',
              'S002594_H3K4me3-PF382-98099':'S002593_Input-PF382-98098'}

def get_input(wildcards):
    return 'output/bowtie/' + CHIP2INPUT[wildcards.smp] + '.bam'

INDEX = '/home/m2bsg_2018_2019/genomes/human/ucsc/bowtie2_index/hg38'


rule final: 
	input: expand('output/fastqc/{smp}/done', smp = SAMPLE),      \
		expand('output/bowtie/{smp}.bam.bai', smp = SAMPLE),  \
		expand('output/bamCoverage/{smp}.bw', smp = SAMPLE),  \
                expand('output/macs2/{smp}/{smp}_peaks.bed', smp=CHIP)
	threads: 1


rule fastqc: #quality control
	input: r1 = '/home/m2bsg_2018_2019/data/run_231/fastq/{smp}_R1.fq.gz', 
		r2 = '/home/m2bsg_2018_2019/data/run_231/fastq/{smp}_R2.fq.gz'
	output: 'output/fastqc/{smp}/done'
	threads: 1
	shell: '''fastqc - o output/fastqc/{wildcards.smp} {input.r1} {input.r2}
		touch {output}
		'''

rule sickle: #trim
	input: r1 = '/home/m2bsg_2018_2019/data/run_231/fastq/{smp}_R1.fq.gz', 
		r2 = '/home/m2bsg_2018_2019/data/run_231/fastq/{smp}_R2.fq.gz'
	output: r1 = 'output/sickle/{smp}_R1.fq.gz',
		r2 = 'output/sickle/{smp}_R2.fq.gz',
		sing = 'output/sickle/{smp}_sing.fq.gz'
	threads: 1
	shell: 'sickle pe -f {input.r1} -r {input.r2} -o {output.r1} -p {output.r2} -s{output.sing} -t sanger -g'

rule bowtie: #map
        input: r1 = 'output/sickle/{smp}_R1.fq.gz', r2 = 'output/sickle/{smp}_R2.fq.gz',
        output: 'output/bowtie/{smp}.bam',
	params: ind = INDEX
	threads: 25
        shell:''' bowtie2 -p 20 -x {params.ind} -1 {input.r1} -2 {input.r2} | samtools view -bS -q 30 - | samtools sort -@5 -m 20G - output/bowtie/{wildcards.smp}'''

rule samtools_index: #index
	input:'output/bowtie/{smp}.bam'
	output:'output/bowtie/{smp}.bam.bai'
	threads: 1
	shell: 'samtools index {input}'

rule bigwig:
	input: bam = 'output/bowtie/{smp}.bam', bai = 'output/bowtie/{smp}.bam.bai'
	output: 'output/bamCoverage/{smp}.bw'
	threads: 1
	shell:'bamCoverage -bs 25 -b {input.bam} -o {output}'



rule macs2:
	input: bam = 'output/bowtie/{smp}.bam', bai = 'output/bowtie/{smp}.bam.bai', inp=get_input
	output: 'output/macs2/{smp}/{smp}_peaks.xls'
	threads: 1
	shell:'macs2 callpeak -t {input.bam} -c {input.inp} --broad -g hs --broad-cutoff 0.1 \
               --outdir output/macs2/{wildcards.smp} -n {wildcards.smp} '
rule macs_to_bed:
	input: 'output/macs2/{smp}/{smp}_peaks.xls'
	output: 'output/macs2/{smp}/{smp}_peaks.bed'
	threads: 1
	shell:""" cat {input}  | grep -v "^#" | grep -v "start"| awk 'BEGIN{{FS=OFS="\\t"}}{{print $1,$2,$3,$9,$8,"+"}}' > {output}"""