## This pipeline was used for extract the cell barcode, trim adaptor sequence from RNA fastq file and write the three barcode 5 end of read1

snakemake --latency-wait 50 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &


  
#### 1. trim specific sequence in the 5 end of read1 
```bash

# bash code

cutadapt -j 11 --action=retain  -g  'NNNNNNNNNNNNNNCCATTCCAGCAGCGTGTGCGAACTCAGACCNNNNNNNATCCACGTGCTTGAGAGGCCAGAGCATTCG;min_overlap=82; max_error_rate=0.086;' --minimum-length 124:146  -o {output[0]} -p {output[1]} --untrimmed-output {output[2]} --untrimmed-paired-output {output[3]}  {input[0]} {input[1]}
```
#### 2. extract the three barcodes and write the three barcode at 5 end of read1
```

export TMPDIR=/datacommons/ydiaolab/diaolab_group/xiaolin/scratch/tmp/
awk  '{{if(NR%2==0) {{print  substr($0,82,6)  substr($0,46,6) substr($0,1,6)  substr($0,7,8) substr($0,length($0)-10+1,10) substr($0,88)}} else  {{print $0 }} }}' {input}  > {output}
```
#### 3. remove the TSO sequence in the 5 end of read2
```

cutadapt  -e 0.15 --action=trim --minimum-length 27  -G 'XCAGTGGTATCAACGCAGAGTACATGGG;min_overlap=10' -o {output[0]} -p {output[1]} {input[0]} {input[1]}
```
#### 4. remove the polyA seqence and adaptor sequence in the 3 end of read2
```

cutadapt -e 0.15 --action=trim --minimum-length 27  -A 'AAAAAAAAAAAAAAAA;min_overlap=14' -o {output[0]} -p {output[1]} {input[0]} {input[1]}
cutadapt -e 0.15 --action=trim --minimum-length 27  -A 'GCTGTANNNNNNCGAATGCTCTGGCCT;min_overlap=20' -o {output[2]} -p {output[3]} {output[0]} {output[1]}
```
#### 5. split the fastq file into two fastq. one was generated using oligo_dT primer and the other using random hexamer primer.
```

cutadapt -Z -j {threads} -e 0.2 --action=retain -g file:ME_index --minimum-length 34:27\
        -o 00_raw_fq_update/{wildcards.sample}_adapter1{{name}}_L001_R1_001.fastq \
        -p 00_raw_fq_update/{wildcards.sample}_adapter1{{name}}_L001_R2_001.fastq \
        {input[0]}  {input[1]}
```

#### 6. extract total barcodes list, count each barcode
```

export TMPDIR=/work/xw171/tmp/
# awk  '{{if(NR%4==1) print substr($0,2,30)}}' {input} | sort --parallel={threads} --temporary-directory={tmp} -S 20G | uniq -c | sort -nr   > {output}
awk  '{{if(NR%4==2) print substr($0,1,18) }}' {input[0]} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output[0]}
awk  '{{if(NR%4==2) print substr($0,1,18) }}' {input[1]} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output[1]}
awk  '{{if(NR%4==2) print substr($0,1,18) }}' {input[2]} | sort --parallel={threads} --temporary-directory={tmp}  | uniq -c | sort -nr   > {output[2]}
```
#### 7. compare the extracted barcodes with the whitelist
```

script/barcode_hash_v2.py
```
#### 8. correct the barcode which barcode only have one mismatch
```

script/fq_barcode_correction_R1.py
```
#### 9. compress fastq file and combine the two fastq file：one generated using oligo_dT primer and the other using random hexamer primer.
```

pigz -p {threads} {input}
cat {input[1]} {output[0]} > {output[1]}
```
