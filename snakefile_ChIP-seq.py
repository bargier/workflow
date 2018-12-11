workdir: '/gpfs/tgml/bargier/projet/run_235'

FASTQ_DIR = '/gpfs/tgml/reads/quality/qc-ngs-tgml/runs/Run_235_NS500-142_06.12.2018_SK_BSoS/data/fastq/'
file_handler = open("sample.txt")
CHIP2INPUT = dict()


for line in file_handler:
	line = line.rstrip("\r").rstrip("\n")
	token = line.split("\t")
	CHIP2INPUT[token[1]] = token[0]

SAMPLE = list(CHIP2INPUT.keys()) + list(CHIP2INPUT.values())

def get_input(wildcards):
    return 'output/bowtie/' + CHIP2INPUT[wildcards.smp] + '.bam'

def get_input_bai(wildcards):
    return 'output/bowtie/' + CHIP2INPUT[wildcards.smp] + '.bam.bai'

INDEX = '/gpfs/tgml/reads/quality/qc-ngs-tgml/Genomes/Homo_sapiens/GRCh38.primary_assembly.genome'


rule final: 
	input: expand('output/fastqc/{smp}/done', smp = SAMPLE),      \
		expand('output/bowtie/{smp}.bam.bai', smp = SAMPLE),  \
		expand('output/bamCoverage/{smp}.bw', smp = SAMPLE),  \
                expand('output/macs2/{smp}/{smp}_peaks.bed', smp=CHIP2INPUT.keys())
	threads: 1

rule fastqc: #quality control
	input: r1 = FASTQ_DIR + '{smp}_R1_001.fastq.gz', 
		r2 = FASTQ_DIR + '{smp}_R2_001.fastq.gz'
	output: 'output/fastqc/{smp}/done'
	threads: 1
	conda: "env/env.yaml"
	shell: '''fastqc -o output/fastqc/{wildcards.smp} {input.r1} {input.r2}
		touch {output}
		'''

rule sickle: #trim
	input: r1 = FASTQ_DIR + '{smp}_R1_001.fastq.gz', 
		r2 = FASTQ_DIR + '{smp}_R2_001.fastq.gz'
	output: r1 = 'output/sickle/{smp}_R1_001.fastq.gz',
		r2 = 'output/sickle/{smp}_R2_001.fastq.gz',
		sing = 'output/sickle/{smp}_sing.fastq.gz'
	threads: 1
	conda: "env/env.yaml"
	shell: 'sickle pe -f {input.r1} -r {input.r2} -o {output.r1} -p {output.r2} -s{output.sing} -t sanger -g'

rule bowtie: #map
        input: r1 = 'output/sickle/{smp}_R1_001.fastq.gz', r2 = 'output/sickle/{smp}_R2_001.fastq.gz',
        output: 'output/bowtie/{smp}.bam',
	params: ind = INDEX
	threads: 15
	conda: "env/env.yaml"
	shell:''' bowtie2 -p 10 -x {params.ind} -1 {input.r1} -2 {input.r2} | samtools view -bS -q 30 - | samtools sort -@5 -m 20G - output/bowtie/{wildcards.smp}'''

rule samtools_index: #index
	input:'output/bowtie/{smp}.bam'
	output:'output/bowtie/{smp}.bam.bai'
	threads: 1
	conda: "env/env.yaml"
	shell: 'samtools index {input}'

rule bigwig:
	input: bam = 'output/bowtie/{smp}.bam', bai = 'output/bowtie/{smp}.bam.bai'
	output: 'output/bamCoverage/{smp}.bw'
	threads: 1
	conda: "env/env.yaml"
	shell:'bamCoverage -bs 25 -b {input.bam} -o {output}'



rule macs2:
	input: bam = 'output/bowtie/{smp}.bam', bai = 'output/bowtie/{smp}.bam.bai', inp=get_input, inp_bai=get_input_bai
	output: 'output/macs2/{smp}/{smp}_peaks.xls'
	threads: 1
	conda: "env/env.yaml"
	shell:'macs2 callpeak -t {input.bam} -c {input.inp} --broad -g hs --broad-cutoff 0.1 \
               --outdir output/macs2/{wildcards.smp} -n {wildcards.smp} '
rule macs_to_bed:
	input: 'output/macs2/{smp}/{smp}_peaks.xls'
	output: 'output/macs2/{smp}/{smp}_peaks.bed'
	threads: 1
	conda: "env/env.yaml"
	shell:""" cat {input}  | grep -v "^#" | grep -v "start"| awk 'BEGIN{{FS=OFS="\\t"}}{{print $1,$2,$3,$9,$8,"+"}}' > {output}"""