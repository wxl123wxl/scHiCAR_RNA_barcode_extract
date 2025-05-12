shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")
import collections
configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON'])) ##  fastq for each lane
SAMPLES = sorted(FILES.keys())
BWA_INDEX = config['BWA_INDEX']
chromsizes = config['chromsizes']
cool_bin       = config['cool_bin']
TARGETS = []
genome = config['genome']
tmp='/work/xw171/tmp'
barcode = json.load(open('barcode.json')) # json file to include barcode file location.

genome_version = 'hs'
smooth_window  = 150
shiftsize      = -75
pval_thresh    = 0.01


def load_group(grouptxt): # group per sample key
    groupid = collections.defaultdict(list)
    with open(grouptxt,'r') as f:
        for line in f:
            line = line.strip()
            if not line =='':
                group, value = line.split("\t")
                groupid[group].append(value)
    return groupid

groupid = load_group('group.txt')



# frag_path = config['frag_path']

## constructe the target if the inputs are fastqs
# bam = expand("01_bam/{sample}.bam", sample = SAMPLES)
# raw_pairsam = expand("pairs-hg38/{sample}/{sample}.raw.pairsam.gz", sample = SAMPLES)

# peak_pairs = expand("peaks-{genome}/{sample}.final.pairs.gz", sample = SAMPLES, genome = genome)
# all_pairs = expand("filtered-{genome}/{sample}.valid.pairs.gz", sample = SAMPLES, genome = genome)

## update the fastq ## add the gzip later
#TARGETS.extend(expand("00_raw_fq_update/{sample}_cut_L001_R1_001.fastq.gz",sample = SAMPLES ))
#TARGETS.extend(expand("00_raw_fq_update/{sample}_adapter23_L001_R1_001.fastq",sample = SAMPLES ))
#TARGETS.extend(expand("00_raw_fq_update/{sample}_cutHicar1_L001_R1_001.fastq",sample = SAMPLES ))
# TARGETS.extend(expand("02_barcode_info/{sample}.barcode_final_map",sample = SAMPLES ))
# TARGETS.extend(expand("02_barcode_info/{sample}.barcode_final_map",sample = SAMPLES ))
TARGETS.extend(expand("03_corrected_fq/{sample}_all_L001_R1_001.fastq.gz",sample = SAMPLES ))
TARGETS.extend(expand("03_corrected_fq/{sample}_all_L001_R2_001.fastq.gz",sample = SAMPLES ))
TARGETS.extend(expand("03_corrected_fq/{sample}_ME_L001_R1_001.fastq.gz",sample = SAMPLES ))
TARGETS.extend(expand("03_corrected_fq/{sample}_ME_L001_R2_001.fastq.gz",sample = SAMPLES ))
# TARGETS.extend(expand("03_ME_fq/{sample}_cutME_L001_R1_001.fastq.gz",sample = SAMPLES ))
#TARGETS.extend(expand("05_ME_fq/{sample}_cutbind_cutME_L001_R2_001.fastq.gz",sample = SAMPLES ))
#TARGETS.extend(expand("05_ME_fq/{sample}_cutbind_cutRNAunknown_L001_R2_001.fastq.gz",sample = SAMPLES ))
# TARGETS.extend(expand("03_ME_fq/{sample}_cutMeKpnI1_cutME_cutBind_L001_R2_001.fastq.gz",sample = SAMPLES ))

## looking for the cellular deduplicated pairs. 
# TARGETS.extend(expand("bulk_pairs_summary-{genome}/{sample}.read_summary", sample = SAMPLES, genome = genome ))  ## cis/trans ratio for each fastq pairs.
# TARGETS.extend(expand("R2-{genome}/{sample}.tsv.gz",  sample = groupid.keys(), genome = genome))


# TARGETS.extend(expand("merged_pairs-{genome}/{sample}.merged.pairs.gz",  sample = groupid.keys(), genome = genome))
# TARGETS.extend(expand("merged_peaks-{genome}/{sample}.merged.peaks.pairs.gz",  sample = groupid.keys(), genome = genome))
# TARGETS.extend(expand("merged_peaks-{genome}/{sample}.R2.end.gz"
# # "merged_peaks-{genome}/{sample}.merged.peaks.pairs.gz",
# # "merged_peaks-{genome}/{sample}.R2_end.gz"
# # "merged_pairs-{genome}/{sample}.merged.pairs.gz",

# ## update after umap result ready
# # TARGETS.extend(expand("split_pairs-{genome}/{sample}/{sample}.splited.pairs.gz", genome = genome, sample = barcode.keys()))
# TARGETS.extend(expand("split_coolers-{genome}/{sample}.{cool_bin}.cool",  genome = genome, sample = barcode.keys(), cool_bin = cool_bin))

# TARGETS.extend(expand("macs2_peak_{sample}/{sample}_{genome}_peaks.narrowPeak", genome = genome, sample = barcode.keys()))
# TARGETS.extend(expand("longreads/{sample}-{cool_bin}-{genome}/{sample}.chr17.long.intra.bedpe", genome = genome, sample = barcode.keys(), cool_bin =  cool_bin))
# TARGETS.extend(expand("short_reads/{sample}-{genome}/{sample}.chr1.shrt.vip.bed",  genome = genome, sample = barcode.keys()))
localrules: targetfiles, read_info_summary
rule targetfiles:
    input: TARGETS




## run for each fastq pairs individually. 

rule raw_fq_adapter:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: 
        r1 = "00_raw_fq_update/{sample}_adapter_L001_R1_001.fastq",
        r2 = "00_raw_fq_update/{sample}_adapter_L001_R2_001.fastq",
        r3 = "00_raw_fq_update/{sample}_NOadapter_L001_R1_001.fastq",
        r4 = "00_raw_fq_update/{sample}_NOadapter_L001_R2_001.fastq"
    shell:
        """
        cutadapt -j 11 --action=retain  -g  'NNNNNNNNNNNNNNCCATTCCAGCAGCGTGTGCGAACTCAGACCNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG;min_overlap=82; max_error_rate=0.086;' --minimum-length 124:146  -o {output[0]} -p {output[1]} --untrimmed-output {output[2]} --untrimmed-paired-output {output[3]}  {input[0]} {input[1]}
        """

rule raw_fq_adapter123:
    input:
        r1 = "00_raw_fq_update/{sample}_adapter_L001_R1_001.fastq"
    output: 
        r1 = "00_raw_fq_update/{sample}_adapter123_umi_L001_R1_001.fastq"
    shell:
        """
        export TMPDIR=/datacommons/ydiaolab/diaolab_group/xiaolin/scratch/tmp/
        awk  '{{if(NR%2==0) {{print  substr($0,82,6)  substr($0,46,6) substr($0,1,6)  substr($0,7,8) substr($0,length($0)-10+1,10) substr($0,88)}} else  {{print $0 }} }}' {input}  > {output} 
        """

rule cut_R2_tso:
    input:
        r1 = "00_raw_fq_update/{sample}_adapter123_umi_L001_R1_001.fastq",
        r2 = "00_raw_fq_update/{sample}_adapter_L001_R2_001.fastq"
        #r3 = "00_raw_fq_update/{sample}_adapter1ME_adapter23_L001_R1_001.fastq",
        #r4 = "00_raw_fq_update/{sample}_adapter1ME_adapter23_L001_R2_001.fastq"
    output: 
        r1 = temp("00_raw_fq_update/{sample}_cut_tso2_L001_R1_001.fastq.gz"),
        r2 = temp("00_raw_fq_update/{sample}_cut_tso2_L001_R2_001.fastq.gz")
        #r1 = temp("01_raw_fq_update/{sample}_index_L001_R1_001.fastq"),
        #r2 = temp("01_raw_fq_update/{sample}_index_L001_R2_001.fastq")
    shell:
        "cutadapt  -e 0.15 --action=trim --minimum-length 27  -G 'XCAGTGGTATCAACGCAGAGTACATGGG;min_overlap=10' -o {output[0]} -p {output[1]} {input[0]} {input[1]}"


rule cut_R2_polyA:
    input:
        r1 = "00_raw_fq_update/{sample}_cut_tso2_L001_R1_001.fastq.gz",
        r2 = "00_raw_fq_update/{sample}_cut_tso2_L001_R2_001.fastq.gz"
        #r3 = "00_raw_fq_update/{sample}_adapter1ME_adapter23_L001_R1_001.fastq",
        #r4 = "00_raw_fq_update/{sample}_adapter1ME_adapter23_L001_R2_001.fastq"
    output: 
        r1 = "00_raw_fq_update/{sample}_cut_tso2_polyA_L001_R1_001.fastq",
        r2 = "00_raw_fq_update/{sample}_cut_tso2_polyA_L001_R2_001.fastq",
        r3 = "00_raw_fq_update/{sample}_cut_R2adapter_L001_R1_001.fastq",
        r4 = "00_raw_fq_update/{sample}_cut_R2adapter_L001_R2_001.fastq"
        #r1 = temp("01_raw_fq_update/{sample}_index_L001_R1_001.fastq"),
        #r2 = temp("01_raw_fq_update/{sample}_index_L001_R2_001.fastq")
    shell:
        """
        cutadapt -e 0.15 --action=trim --minimum-length 27  -A 'AAAAAAAAAAAAAAAA;min_overlap=14' -o {output[0]} -p {output[1]} {input[0]} {input[1]}
        cutadapt -e 0.15 --action=trim --minimum-length 27  -A 'GCTGTANNNNNNCGAATGCTCTGGCCT;min_overlap=20' -o {output[2]} -p {output[3]} {output[0]} {output[1]}
        """
rule raw_fq_ME_oligdT:
    input:
        r1 = "00_raw_fq_update/{sample}_cut_R2adapter_L001_R1_001.fastq",
        r2 = "00_raw_fq_update/{sample}_cut_R2adapter_L001_R2_001.fastq"
    output:
        #r1 = "00_raw_fq_update/{sample}_adapter1noME_L001_R1_001.fastq",
        #r2 = "00_raw_fq_update/{sample}_adapter1noME_L001_R2_001.fastq",
        "00_raw_fq_update/{sample}_adapter1ME_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1ME_L001_R2_001.fastq",
        "00_raw_fq_update/{sample}_adapter1RNA_oligodT_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1RNA_oligodT_L001_R2_001.fastq",
        "00_raw_fq_update/{sample}_adapter1RNA_randomN_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1RNA_randomN_L001_R2_001.fastq"
        #r1 = temp("01_raw_fq_update/{sample}_index_L001_R1_001.fastq"),
        #r2 = temp("01_raw_fq_update/{sample}_index_L001_R2_001.fastq")
    shell:
        """
        cutadapt -Z -j {threads} -e 0.2 --action=retain -g file:ME_index --minimum-length 34:27\
        -o 00_raw_fq_update/{wildcards.sample}_adapter1{{name}}_L001_R1_001.fastq \
        -p 00_raw_fq_update/{wildcards.sample}_adapter1{{name}}_L001_R2_001.fastq \
        {input[0]}  {input[1]}
        """




rule barcode_QC: ## extract total barcodes list
    input: 
        "00_raw_fq_update/{sample}_adapter1RNA_oligodT_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1RNA_randomN_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1ME_L001_R1_001.fastq"
    #input: "01_raw_fq_update/{sample}_index_L001_R1_001.fastq"
    output: 
        "02_barcode_info/{sample}_oligodT_raw_barcode_count.txt",
        "02_barcode_info/{sample}_randomN_raw_barcode_count.txt",
        "02_barcode_info/{sample}_ME_raw_barcode_count.txt"
    threads: 8 
    shell:
        """
        export TMPDIR=/work/xw171/tmp/
        # awk  '{{if(NR%4==1) print substr($0,2,30)}}' {input} | sort --parallel={threads} --temporary-directory={tmp} -S 20G | uniq -c | sort -nr   > {output} 
        awk  '{{if(NR%4==2) print substr($0,1,18) }}' {input[0]} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output[0]} 
        awk  '{{if(NR%4==2) print substr($0,1,18) }}' {input[1]} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output[1]}
        awk  '{{if(NR%4==2) print substr($0,1,18) }}' {input[2]} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output[2]}
        """

        # update the temp dir for sort
        # the merged_barcode sequence is at reads name
        # double check the reads formate and the part of barcode to included.  
        #  awk  '{{if(NR%4==1) print substr($0,2,18)}}' {input} | sort | uniq -c | sort -nr   > {output} 

rule find_right_barcodes: # mismatch correction for each barcode
    input: 
        "02_barcode_info/{sample}_oligodT_raw_barcode_count.txt",
    output: sum = "02_barcode_info/{sample}_oligodT.barcode_final_summary",
            map = "02_barcode_info/{sample}_oligodT.barcode_final_map",
            log = "02_barcode_info/{sample}_oligodT.barcode_log"
    script:
        "script/barcode_hash_v2.py"

rule find_right_barcodes_randomN: # mismatch correction for each barcode
    input: 
        "02_barcode_info/{sample}_randomN_raw_barcode_count.txt"
    output: sum = "02_barcode_info/{sample}_randomN.barcode_final_summary",
            map = "02_barcode_info/{sample}_randomN.barcode_final_map",
            log = "02_barcode_info/{sample}_randomN.barcode_log"
    script:
        "script/barcode_hash_v2.py"

rule find_right_barcodes_ME: # mismatch correction for each barcode
    input: 
        "02_barcode_info/{sample}_ME_raw_barcode_count.txt"
    output: sum = "02_barcode_info/{sample}_ME.barcode_final_summary",
            map = "02_barcode_info/{sample}_ME.barcode_final_map",
            log = "02_barcode_info/{sample}_ME.barcode_log"
    script:
        "script/barcode_hash_v2.py"


rule read1_barcode_correction:
    input :
        "00_raw_fq_update/{sample}_adapter1RNA_oligodT_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1RNA_oligodT_L001_R2_001.fastq",
        "02_barcode_info/{sample}_oligodT.barcode_final_map"
        #"01_raw_fq_update/{sample}_index_L001_R1_001.fastq",
        #"01_raw_fq_update/{sample}_index_L001_R2_001.fastq",
    output :
        "03_corrected_fq/{sample}_oligodT_L001_R1_001.fastq",
        "03_corrected_fq/{sample}_oligodT_L001_R2_001.fastq"
    log: "00_log/{sample}_L001_R1_corrected.log"
    script:
        "script/fq_barcode_correction_R1.py"

rule read1_barcode_correction_randomN:
    input :
        "00_raw_fq_update/{sample}_adapter1RNA_randomN_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1RNA_randomN_L001_R2_001.fastq",
        "02_barcode_info/{sample}_randomN.barcode_final_map"
        #"01_raw_fq_update/{sample}_index_L001_R1_001.fastq",
        #"01_raw_fq_update/{sample}_index_L001_R2_001.fastq",
    output :
        "03_corrected_fq/{sample}_randomN_L001_R1_001.fastq",
        "03_corrected_fq/{sample}_randomN_L001_R2_001.fastq"
    log: "00_log/{sample}_L001_R1_corrected.log"
    script:
        "script/fq_barcode_correction_R1.py"

rule read1_barcode_correction_ME:
    input :
        "00_raw_fq_update/{sample}_adapter1ME_L001_R1_001.fastq",
        "00_raw_fq_update/{sample}_adapter1ME_L001_R2_001.fastq",
        "02_barcode_info/{sample}_ME.barcode_final_map"
        #"01_raw_fq_update/{sample}_index_L001_R1_001.fastq",
        #"01_raw_fq_update/{sample}_index_L001_R2_001.fastq",
    output :
        "03_corrected_fq/{sample}_ME_L001_R1_001.fastq",
        "03_corrected_fq/{sample}_ME_L001_R2_001.fastq"
    log: "00_log/{sample}_L001_R1_corrected.log"
    script:
        "script/fq_barcode_correction_R1.py"

#rule read2_barcode_correction:
#    input :
#        "01_raw_fq_update/{sample}_index_L001_R2_001.fastq",
#        "02_barcode_info/{sample}.barcode_final_map"
#    output :
#        #"03_corrected_fq/{sample}_L001_R2_001.fastq"
#    log:"00_log/{sample}_L001_R2_corrected.log"   
#    script:
#        "script/fq_barcode_correction.py"
#
rule r1_zip:
    input  : "03_corrected_fq/{sample}_oligodT_L001_R1_001.fastq"
    output : "03_corrected_fq/{sample}_oligodT_L001_R1_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"

rule r2_zip:
    input  : "03_corrected_fq/{sample}_oligodT_L001_R2_001.fastq"
    output : "03_corrected_fq/{sample}_oligodT_L001_R2_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"

rule r1_zip_randomN:
    input  : 
        "03_corrected_fq/{sample}_randomN_L001_R1_001.fastq",
        "03_corrected_fq/{sample}_oligodT_L001_R1_001.fastq.gz"
    output : 
        "03_corrected_fq/{sample}_randomN_L001_R1_001.fastq.gz",
        "03_corrected_fq/{sample}_all_L001_R1_001.fastq.gz"
    threads: 11
    shell:
        """
        pigz -p {threads} {input[0]}
        cat {input[1]} {output[0]} > {output[1]}
        """


rule r2_zip_randomN:
    input  : 
        "03_corrected_fq/{sample}_randomN_L001_R2_001.fastq",
        "03_corrected_fq/{sample}_oligodT_L001_R2_001.fastq.gz"
    output : 
        "03_corrected_fq/{sample}_randomN_L001_R2_001.fastq.gz",
        "03_corrected_fq/{sample}_all_L001_R2_001.fastq.gz"
    threads: 11
    shell:
        """
        pigz -p {threads} {input[0]}
        cat {input[1]} {output[0]} > {output[1]}
        """

rule r1_zip_ME:
    input  : "03_corrected_fq/{sample}_ME_L001_R1_001.fastq"
    output : "03_corrected_fq/{sample}_ME_L001_R1_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"

rule r2_zip_ME:
    input  : "03_corrected_fq/{sample}_ME_L001_R2_001.fastq"
    output : "03_corrected_fq/{sample}_ME_L001_R2_001.fastq.gz"
    threads: 11
    shell:
        "pigz -p {threads} {input}"

rule index_read2:
    input:
        r1 = "03_corrected_fq/{sample}_L001_R1_001.fastq",
        r2 = "03_corrected_fq/{sample}_L001_R2_001.fastq"
    output: 
        r1 = "04_index_fq/{sample}_index2_L001_R1_001.fastq",
        r2 = "04_index_fq/{sample}_index2_L001_R2_001.fastq"
        #r1 = temp("01_raw_fq_update/{sample}_index_L001_R1_001.fastq"),
        #r2 = temp("01_raw_fq_update/{sample}_index_L001_R2_001.fastq")
    script:
        "script/index_update.py"

#rule cutBind:
#    input:
#        "03_ME_fq/{sample}_cutME_L001_R1_001.fastq.gz",
#        "03_ME_fq/{sample}_cutME_L001_R2_001.fastq.gz"
#    output:
#        "03_ME_fq/{sample}_cutME_cutBind_L001_R1_001.fastq.gz",
#        "03_ME_fq/{sample}_cutME_cutBind_L001_R2_001.fastq.gz"
#    threads: 16
#    shell:
#        """
#        cutadapt -Z -j {threads} -e 0.2 --quality-cutoff 10,10 -g file:bind_index \
#        -o 03_ME_fq/{wildcards.sample}_cutME_cut{{name}}_L001_R1_001.fastq.gz \
#        -p 03_ME_fq/{wildcards.sample}_cutME_cut{{name}}_L001_R2_001.fastq.gz \
#        {input[0]}  {input[1]}  
#        """

rule cutTSO:
    input:
        "03_ME_fq/{sample}_cutME_cutBind_L001_R1_001.fastq.gz",
        "03_ME_fq/{sample}_cutME_cutBind_L001_R2_001.fastq.gz"
    output:
        #"03_ME_fq/{sample}_cutMeKpnI2_cutME_cutBind_L001_R1_001.fastq.gz",
        "03_ME_fq/{sample}_cutTSO_L001_R1_001.fastq.gz",
        "03_ME_fq/{sample}_cutTSO_L001_R2_001.fastq.gz"
        #"03_ME_fq/{sample}_cutMeKpnI3_cutME_cutBind_L001_R2_001.fastq.gz"
    threads: 16
    shell:
        """
        cutadapt -Z -j {threads} -e 0.25 --quality-cutoff 10,10 -a "ACTCTGCGTTGATACCACTGCTT;min_overlap=13" -G "AAGCAGTGGTATCAACGCAGAGTNNNNNNNNNNNNN;min_overlap=13"\
        -o {output[1]} \
        -p {output[0]} \
        {input[1]}  {input[0]}  
        """

rule bwa_align:
    input:
        r1 = "03_ME_fq/{sample}_cutTSO_L001_R1_001.fastq.gz",
        r2 = "03_ME_fq/{sample}_cutTSO_L001_R2_001.fastq.gz"
        #r1 = "03_corrected_fq/{sample}_L001_R2_001.fastq.gz1",
        #r2 = "03_corrected_fq/{sample}_L001_R2_001.fastq.gz1"
    output: "04_bam/{sample}.bam"
    threads: 16
    message: "bwa {input}: {threads} threads"
    log:
         "00_log/{sample}.bwa"
    shell:
        """
        module load BWA ## dcc module
        module load samtools
        bwa mem  -SP -t {threads} {BWA_INDEX} {input} | samtools view -bS - > {output}  2> {log}
        """

rule prase_bam_to_pairs: ## no flip to keep the R1 R2 position for the peak calling, extract 28 bp cell barcode 
    input:  "04_bam/{sample}.bam"
    output: "pairs-{genome}/{sample}.select.pairsam.gz",  "pairs-{genome}/{sample}.raw.pairsam.stat"
    message: "prase bam {input} "
    threads: 10
    shell:
        """
        module load samtools ## set the min-mapq as 10
        pairtools parse -c {chromsizes}  \
        --assembly {genome} --min-mapq 10 \
        --max-molecule-size 2000 --max-inter-align-gap 20 \
        --walks-policy mask  --no-flip --drop-seq --drop-sam  \
        --output-stats {output[1]}  {input} | pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' | \
        awk  -F $"\\t" 'BEGIN {{OFS=FS}} ; {{ if($1 ~ /^#/) {{ print $0}} else {{ print substr($1,1,28),$2,$3,$4,$5,$6,$7,$8}}}}' | \
        pbgzip -n {threads} -c > {output[0]}
        """

rule flip_pairs_unique_map_sort_by_position: ## flip and select valid pairs
    input:  "pairs-{genome}/{sample}.select.pairsam.gz"
    output: "pairs-{genome}/{sample}.fliped.select.pairs.gz"
    message: "flip and sort {input} "
    threads: 8
    shell:
        """
         pairtools flip -c {chromsizes} {input} -o {output}
        """


rule flip_pairs_sort_by_cells: # sort pairs by cell and position
    input:  "pairs-{genome}/{sample}.fliped.select.pairs.gz"
    output: "pairs-{genome}/{sample}.cell.sorted.pairs"
    message: " sort by cells {input} "
    threads: 8
    shell:
        """
        # zcat {input} | head -n 2000 | grep "^#" > {output} || true
        export LC_COLLATE=C; export LANG=C;
        zcat {input} | grep -v "^#" | \
        sort -V -k1,1 -k2,2 -k4,4 -k3,3n -k5,5n --stable  --parallel={threads} --temporary-directory={tmp} -S 20G  > {output}
        """

rule gzip_pairs: # sort pairs by cell and position
    input:  "pairs-{genome}/{sample}.cell.sorted.pairs"
    output: "pairs-{genome}/{sample}.cell.sorted.pairs.gz"
    threads: 12
    shell: 
        "pbgzip -n {threads} {input}"

rule cell_pairs_dedup:
    input:  "pairs-{genome}/{sample}.cell.sorted.pairs",  "pairs-{genome}/{sample}.select.pairsam.gz"
    output: "filtered-{genome}/{sample}.dedup.pairs.gz" ,"filtered-{genome}/{sample}.dedup.pairs.stat", "filtered-{genome}/{sample}.dedup.pairs.head"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
        zcat {input[1]} | head -n 2000 | grep "^#" > {output[2]} || true
        pairtools dedup --max-mismatch 1 --method max -o {output[0]}  --output-stats  {output[1]} \
         < (cat {output[2]} {input[0]})
        """

rule read_info_summary:  # get basic hic quality for each fastq files
    input:  "pairs-{genome}/{sample}.raw.pairsam.stat", "filtered-{genome}/{sample}.dedup.pairs.stat"
    output: "bulk_pairs_summary-{genome}/{sample}.read_summary"
    threads: 1
    script:
        "script/read_summary.R"



rule peaks_file_select: ## select valid peaks pairs 
    input:  "pairs-{genome}/{sample}.raw.pairsam.gz"
    output: "peaks-{genome}/{sample}.selected.pairs.gz"
    threads: 8
    shell:
        """
        pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' {input}  -o {output}
        """

## merge R2 reads from all fastq pairs. use them for the arcR analysis.  


def get_R2_input(wildcards):  ## group functin to define the input sample
    return expand("peaks-{genome}/{sample}.selected.pairs.gz", genome = wildcards.genome,
                                                               sample = groupid[wildcards.sample] )



rule merge_peaks:  ## merge peaks pair, then split by cell id for distance selection
    input:  get_R2_input  ## require select pairs for each fastq pairs
    output: "merged_peaks-{genome}/{sample}.merged.peaks.pairs"
    threads: 8
    shell:
        """
        zcat {input[0]} | head -n 1000 | grep "^#" > {output} || true
        zcat {input} |  grep -v '^#'  >> {output} 
        """

rule merge_peaks_gzip:  ## merge peaks pair, then split by cell id for distance selection
    input:  "merged_peaks-{genome}/{sample}.merged.peaks.pairs"
    output: "merged_peaks-{genome}/{sample}.merged.peaks.pairs.gz"
    threads: 11
    shell: 
        "pbgzip -n {threads} {input}"

rule extract_R2_reads:  ## extract R2 ends reads for cell type classification
    input:  "merged_peaks-{genome}/{sample}.merged.peaks.pairs.gz"
    output: "R2-{genome}/{sample}.tsv.gz","R2-{genome}/{sample}.tsv.gz.tbi"
    threads: 8
    shell:
        """
        export TMPDIR=/work/xw171/tmp
        zcat {input} | grep -v "^#" | \
        awk  -F $"\\t" 'BEGIN {{OFS=FS}} ; {{ {{ if ($7 == "+") {{$5 = $5 + 4}} else if ($7 == "-") {{$5 = $5 - 5}}  print $4, $5, $5+1, "*", "*", $7}} }} '  | \
        sort    -k4,4 -k1,1 -k2,2n --stable  --parallel={threads} --temporary-directory={tmp} -S 20G | uniq | \
        sort -V -k1,1 -k2,2n --stable  --parallel={threads} --temporary-directory={tmp} -S 20G | \
        pbgzip -n {threads} -c  > {output[0]} 
        tabix -p bed {output[0]}
        """


# rule extract_R2_reads_single_lane:  ## extract R2 ends reads for cell type classification
#     input:  "peaks-{genome}/{sample}.selected.pairs.gz"
#     output: "peaks-{genome}/{sample}.tsv.gz", "peaks-{genome}/{sample}.tsv.gz.tbi"
#     threads: 8
#     shell:  ## first round sort remove duplicate in each cells. then sort by reads for archR input.
#         """
#         export TMPDIR=/datacommons/ydiaolab/diaolab_group/xiaolin/scratch/tmp
#         zcat {input} | grep -v "^#" | \
#         awk  -F $"\\t" 'BEGIN {{OFS=FS}} ;  {{ print $4,$5,$5+1,$1}} ' | \
#         sort  -k4,4 -k1,1 -k2,2n --stable  --parallel={threads} --temporary-directory={tmp} -S 20G | uniq | \
#         sort -V -k1,1 -k2,2n --stable  --parallel={threads} --temporary-directory={tmp} -S 20G | \
#         pbgzip -n {threads} -c  > {output[0]} 
#         tabix -p bed {output[0]}
#         """


def get_pairs_input(wildcards):  ## group functin to define the input sample
    return expand("filtered-{genome}/{sample}.dedup.pairs.gz",  
                  genome = wildcards.genome,
                  sample = groupid[wildcards.sample])

rule merge_pairs:
    input: get_pairs_input ## all the dedup pairs based on the groupid dictionary.
    output: "merged_pairs-{genome}/{sample}.merged.pairs.gz",
     temp("merged_pairs-{genome}/{sample}.head"),
     temp("merged_pairs-{genome}/{sample}.body")
    threads: 10
    shell:
        """
        zcat {input[0]} | head -n 1000 |grep "^#" > {output[1]} || true
        zcat {input} | grep -v "^#"  > {output[2]}
        cat  {output[1]}  {output[2]} |  pbgzip -n {threads} -c > {output[0]}
        """

# split the merged pairs after get cell type barcode from UMP clustering.

rule split_pairs_for_each_celltype:  ## the sample in the merged sample 
    input: expand("merged_pairs-{{genome}}/{sample}.merged.pairs.gz", sample = groupid.keys()), ## assuam only one merged pairs!
            lambda wildcards: barcode[wildcards.sample]
    output: "split_pairs-{genome}/{sample}/{sample}.splited.pairs.gz",
     temp("split_pairs-{genome}/{sample}/{sample}.splited.head"),
     temp("split_pairs-{genome}/{sample}/{sample}.splited.body")
    threads: 10
    shell:
        """
        zcat {input[0]} | head -n 1000 |grep "^#" > {output[1]} || true
        awk  -F $"\\t" 'BEGIN {{OFS=FS}} ;  FNR==NR {{f1[$1];next}} ($1 in f1)' \
        {input[1]} <(gzip -dc  {input[0]}) > {output[2]}
        cat  {output[1]}  {output[2]} |  pbgzip -n {threads} -c > {output[0]}
        """
## sort and load into cooler 

rule split_pairs_sort_deduplicate:
    input:  "split_pairs-{genome}/{sample}/{sample}.splited.pairs.gz"
    output: "split_pairs-{genome}/{sample}/{sample}.valid.pairs.gz" ,"split_pairs-{genome}/{sample}/{sample}.valid.pairs.stat"
    message: "sort {input} "
    threads: 8
    shell:
        """
        pairtools sort  --nproc 8  --memory 15G  --tmpdir {tmp} {input} |
        pairtools dedup --max-mismatch 1 --method max -o {output[0]}  \
        --output-stats  {output[1]}
        """

rule make_index:
    input:  "split_pairs-{genome}/{sample}/{sample}.valid.pairs.gz" 
    output: "split_pairs-{genome}/{sample}/{sample}.valid.pairs.gz.px2"
    threads: 2
    shell:
        """
         pairix -p pairs  {input}
        """

rule HiC_contact_matrices_Cooler:
    input:  "split_pairs-{genome}/{sample}/{sample}.valid.pairs.gz" ,"split_pairs-{genome}/{sample}/{sample}.valid.pairs.gz.px2"
    output: "split_coolers-{genome}/{sample}.{cool_bin}.cool"
    message: "cooler {input} "
    # params: res = {cool_bin} 
    threads: 10
    shell:
        """
        cooler cload pairix --assembly {wildcards.genome} --nproc {threads} \
        --max-split 2 {chromsizes}:{wildcards.cool_bin} {input[0]} {output}
        """


## cooler is ready for further maps call

rule dump_long_reads:
    input: "split_coolers-{genome}/{sample}.{cool_bin}.cool"
    output: "longreads/{sample}-{cool_bin}-{genome}/{sample}.chr17.long.intra.bedpe"
    threads: 1
    params: "longreads/{sample}-{cool_bin}-{genome}"
    shell:
        """
        mkdir  -p {params}
        cooler dump -t pixels -H --join {input} | \
        awk -v setname={wildcards.sample} -v outdir={params}   -F $"\\t" '{{if($1 == $4) {{print > outdir"/"setname"."$1".long.intra.bedpe"}} }}'
        """



## prepare the 1D peaks track. 

        # """
        # export TMPDIR=/hpc/home/xw171/scratch/ntmp
        # zcat {input} | grep -v "^#" | \
        # awk  -F $"\\t" 'BEGIN {{OFS=FS}} ;  {{ print $4,$5,$5+1,substr($1,1,28)}} ' | \
        # sort  -k4,4 -k1,1 -k2,2n --stable  --parallel={threads}  -S 20G | uniq | \
        # sort -V  -k1,1 -k2,2n  --stable  --parallel={threads}  -S 20G | bgzip -@ {threads} -c  > {output}
        # """


# run after umap result ready
rule split_peaks_for_each_celltype:
    input: expand("merged_peaks-{{genome}}/{sample}.merged.peaks.pairs.gz", sample = groupid.keys()), ## assuam only one merged pairs!
            lambda wildcards: barcode[wildcards.sample]
    output: "split_peaks-{genome}/{sample}/{sample}.splited.pairs.gz",
     temp("split_peaks-{genome}/{sample}/{sample}.splited.head"),
     temp("split_peaks-{genome}/{sample}/{sample}.splited.body")
    threads: 8
    shell:
        """
        zcat {input[0]} | head -n 1000 |grep "^#" > {output[1]} || true
        awk  -F $"\\t" 'BEGIN {{OFS=FS}} ;  FNR==NR {{f1[$1];next}} ($1 in f1)' \
        {input[1]} <(gzip -dc  {input[0]}) > {output[2]}
        cat  {output[1]}  {output[2]} |  pbgzip -n {threads} -c > {output[0]}
        """


rule R2_for_each_celltypes:
    input:  "split_peaks-{genome}/{sample}/{sample}.splited.pairs.gz"
    output: "split_peaks-{genome}/{sample}/{sample}.longRange_Trans.pairs.gz", "split_peaks-{genome}/{sample}/{sample}.short.pairs.gz"
    message: "flip to filted {input} "
    threads: 6
    shell:
        """
         pairtools select '(chrom1==chrom2) and (abs(pos1 - pos2) < 1e4)'  -o {output[1]}  --output-rest {output[0]}  {input}
        """


rule ATAC_reads_Tn5_shifting_duplicate_remove:  
    input:  "split_peaks-{genome}/{sample}/{sample}.longRange_Trans.pairs.gz"
    output: "split_peaks-{genome}/{sample}/{sample}.R2.ATAC.bed.gz"
    threads: 8
    shell:
        """
        export TMPDIR={tmp}
        zcat {input} | \
        awk ' BEGIN {{OFS="\\t"}} ;  /^[^#]/ {{ {{ if ($7 == "+") {{$5 = $5 + 4}} else if ($7 == "-") {{$5 = $5 - 5}}  print $4, $5, $5+1, "*", "*", $7}} }} ' |\
        sort -k1,1 -k2,2n  --stable  --parallel={threads}  -S 20G  | uniq  | pbgzip -n {threads} -c  > {output}
        """

## dump short reads for maps2 

rule dump_short_reads:
    input: "split_peaks-{genome}/{sample}/{sample}.R2.ATAC.bed.gz"
    output: "short_reads/{sample}-{genome}/{sample}.chr1.shrt.vip.bed"
    threads: 1
    params: "short_reads/{sample}-{genome}"
    shell:
        """
        mkdir -p {params}
        zcat {input[0]} | awk -v setname={wildcards.sample} -v outdir={params}   -F $"\\t"  'BEGIN {{OFS=FS}};{{print $1,$2,$3 > outdir"/"setname"."$1".shrt.vip.bed" }}' 
        """
        # cooler dump -t pixels -H --join {input} | \
        # awk -v setname={wildcards.sample} -v outdir={params} '{{if($1 == $4) {{print > outdir"/"setname"."$1".long.intra.bedpe"}} }}'


rule ATAC_macs2_peaks:
    input:  "split_peaks-{genome}/{sample}/{sample}.R2.ATAC.bed.gz"
    output: "macs2_peak_{sample}/{sample}_{genome}_peaks.narrowPeak"
    threads:1
    params: name = "{sample}_{genome}"
    shell:
        """
        macs2 callpeak -t {input} -f BED -n {params.name}  -g {genome_version} --qval {pval_thresh} \
        --shift {shiftsize} --extsize {smooth_window} --nomodel -B --SPMR --keep-dup all --call-summits \
        --outdir macs2_peak_{wildcards.sample}
        """
