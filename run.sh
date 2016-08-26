#!/bin/bash
#ENCODE SPP, MACS2, IDR pipeline
#usage ./run.sh genome_location_prefix hs/gs/etc chromosize_location

threads=$(( $(grep -c ^processor /proc/cpuinfo) * 2 ));
genome=$1
export gsize=$2
export csize=$3

#alignments, and sam file generation
for f in *.fastq; do
 bwa aln -I -B 0 -t $threads $genome $f > $f.sai;
 bwa samse $genome $f.sai $f.fastq > $f.sam;
done;

#q30 filter
find . -name "*sam" | xargs -n 1 -P $threads -iFILES sh -c 'samtools view -F 1548 -h -S -q 30 -o FILES_q30.sam FILES;';

#generate tag densities for q30
find . -name "*_q30.sam" | xargs -n 1 -P $threads -iFILES sh -c 'samtools view -Sb FILES | bedtools bamtobed -i stdin | awk '\''BEGIN{FS="\t";OFS="\t"}{$4="N";print$0}'\'' | gzip -c > FILES.tagAlign.gz;';

#generate bams
find . -name "*.sam" | xargs -n 1 -P $threads -iFILES sh -c 'samtools view -Sbh FILES > FILES.bam;';

#generate pooled tags from three biological replicates
find . -name "*rep1*.tagAlign.gz" | xargs -n 1 -P $threads -iFILES sh -c 'r1=FILES; r2=$(echo $r1 | sed "s/rep1/rep2"); r3=$(echo $r1 | sed "s/rep1/rep3"); r0=$(echo $r1 | sed "s/rep1/rep0"); zcat $r1 $r2 $r3 > $r0;';

#rename inputs to pooled format
find . -name "*input*.tagAlign.gz" | xargs -n 1 -P $threads -iFILES sh -c 'f=$(echo FILES | sed "s/input\_rep1/input_rep0/"); mv FILES $f;';

#shuffle tags and generate psuedo replicates
find . -name "*rep*tagAlign.gz" | xargs -n 1 -P $threads -iFILES sh -c 'nline=$(( ($(zcat FILES | wc -l) + 1) / 2 )); zcat FILES | shuf | split -d -l $nline - ./FILES.0; gzip ./FILES.000; gzip ./FILES.001; mv ./FILES.000.gz ./FILES.pr1.gz; mv ./FILES.001.gz ./FILES.pr2.gz;';

#spp peak calls
find . -name "*rep*.tagAlign.gz" | xargs -n 1 -P $threads -iFILES sh -c 'inp=$(echo FILES | sed "s/rep./input\.rep0/"); Rscript run_spp.R -c=FILES -i=$inp -npeak=30000 -x=-500:85 -odir=./ -p=10 -savr -savp -rf -out=FILES.cc;';

#macs2 peak calls
find . -name "*rep*.tagAlign.gz" | xargs -n 1 -P $threads -iFILES sh -c 'rep=$(echo FILES | sed "s/\.tagAlign.*//g"); inp=$(echo FILES | sed "s/rep./input\.rep0/"); ss=$(cat FILES.cc | grep FILES | cut -f3 | head -n 1 | cut -d "," -f1); macs2 callpeak -t FILES -c $inp -g $gsize --to-large --nomodel --extsize $ss -p 0.1 -n "$rep"_VS_input;';

#IDR self consistency analysis

#TODO

exit 0
