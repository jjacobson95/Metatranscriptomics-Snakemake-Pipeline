# Metatranscriptomics-Snakemake-Pipeline  
  
  

Metatranscriptomics pipeline to output data ready for differential expression analysis (DEA) or pathway enrichment analysis (PEA).  
  
Pipeline Now Functional.  
  
Steps:  
1)‎ Deinterleave with **BBmap**   
2) Quality Check and Trimming with **fastp**   
3) Remove rRNA using **SortMeRNA** aligner and rRNA reference files  
4) Concatenate all genome assembly and all annotation files (seperately)  
5) Create STAR database using **STAR --runMode genomeGenerate**  
6) Use STAR to align reads using **STAR --runMode alignReads**  
7) Sort files using **samtools**  
8) Generate files ready for DEA or PEA using **htseq-count** with **stranded=reverse**    
  
  <br/>
  
  
### Prepare files for R using bash:  
This step may be added to pipeline.    
To title and join all htseq files into a compact tsv prior to R analysis, consider using the following bash code:  
  
Add header to htseq files. This header is the file name (ie. sample name).   
Must be in directory with no other files.  
  
```
for f in *; do sed -i "1i File\t$f" $f; done
```
  
Join all files together use:  
  
```
function multijoin() { 
        out=$1 
        shift 1 
        cat $1 | awk '{print $1}' > $out 
        for f in $*; do join $out $f > tmp; mv tmp $out; done \
            }
multijoin ../OD_DESeq_input/DESeq2_counts_file.tsv *
```
  
