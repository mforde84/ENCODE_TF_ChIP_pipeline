#!/bin/bash
#ENCODE SPP, MACS2, IDR pipeline
#usage ./run.sh genome_fasta_prefix hs/gs/etc chromosize_location

threads=$(( $(grep -c ^processor /proc/cpuinfo) * 2 ));
genome=$1
export gsize=$2
export csize=$3

#alignments, and sam file generation
for f in *.fastq; do
 bwa aln -I -B 0 -t $threads $genome $f > $f.sai;
 bwa samse $genome $f.sai $f > $f.sam;
done;

#q30 filter
find . -name "*sam" | xargs -n 1 -P $threads -iFILES sh -c 'samtools view -F 1548 -h -S -q 30 -o FILES_q30.sam FILES;';

#generate tag densities for q30
find . -name "*_q30.sam" | xargs -n 1 -P $threads -iFILES sh -c 'samtools view -Sb FILES | bedtools bamtobed -i stdin | awk '\''BEGIN{FS="\t";OFS="\t"}{$4="N";print$0}'\'' | gzip -c > FILES.tagAlign.gz;';

#generate bams
find . -name "*.sam" | xargs -n 1 -P $threads -iFILES sh -c 'samtools view -Sbh FILES > FILES.bam;';

#generate pooled tags from three biological replicates
find . -name "*rep1*.tagAlign.gz" | xargs -n 1 -P $threads -iFILES sh -c 'r1=FILES; r2=$(echo $r1 | sed "s/rep1/rep2/"); r3=$(echo $r1 | sed "s/rep1/rep3/"); r0=$(echo $r1 | sed "s/rep1/rep0/"); zcat $r1 $r2 $r3 | gzip -c > $r0;';

#shuffle tags and generate psuedo replicates
find . -name "*tagAlign.gz"  -not -name "*input*" | xargs -n 1 -P $threads -iFILES sh -c 'nline=$(( ($(zcat FILES | wc -l) + 1) / 2 )); zcat FILES | shuf | split -d -l $nline - ./FILES.0; gzip ./FILES.000; gzip ./FILES.001; mv ./FILES.000.gz ./FILES.pr1.gz; mv ./FILES.001.gz ./FILES.pr2.gz;';

#spp peak calls
for f in $(find . -name "*tagAlign*gz"  -not -name "*input*"); do
 inp=$(echo $f | sed "s/rep./input\.rep0/" | sed "s/\.pr[1-9]\.gz//g"); 
 Rscript run_spp.R -c=$f -i=$inp -npeak=30000 -x=-500:85 -odir=./ -p=10 -savr -savp -rf -out=$f.cc;
done;

#macs2 peak calls
find . -name "*tagAlign*gz"  -not -name "*input*" -not -name "*VS*" | xargs -n 1 -P $threads -iFILES sh -c 'rep=$(echo FILES | sed "s/\.tagAlign.*//g"); inp=$(echo FILES | sed "s/rep./input\.rep0/"); ss=$(cat FILES.cc | grep FILES | cut -f3 | head -n 1 | cut -d "," -f1); macs2 callpeak -t FILES -c $inp -g $gsize --to-large --nomodel --extsize $ss -p 0.1 -n "$rep"_VS_input;';

#IDR self consistency analysis

# sort MACS2 peaks
find . -name "*.narrowPeak" | xargs -n 1 -P $threads -iFILES sh -c 'sort -k 8nr,8nr FILES | head -n 30000 > FILES.sorted;';

# consistency analysis for self-psuedo pr1 pr2 MACS2
find . -name "*pr1*narrowPeak.sorted" | xargs -n 1 -P $threads -iFILES sh -c '
	pr2=$(echo FILES | sed "s/pr1/pr2/g");
	base=$(echo FILES | sed "s/pr1.*/pr/g";
	Rscript batch-consistency-analysis.r FILES $pr2 -1 "$base1_VS_$base2".macs 0 F p.value
	Rscript batch-consistency-plot.r 1 "$base1_VS_$base2".macs "$base1_VS_$base2".macs;
';

# consistency analysis for self-psuedo pr1 pr2 SPP
find . -name "*pr1*VS*input*regionPeak.gz" | xargs -n 1 -P $threads -iFILES sh -c '
	pr2=$(echo FILES | sed "s/pr1/pr2/g");
	base=$(echo FILES | sed "s/pr1.*/pr/g";
	Rscript batch-consistency-analysis.r FILES $pr2 -1 "$base1_VS_$base2" 0 F signal.value;
	Rscript batch-consistency-plot.r 1 "$base1_VS_$base2" "$base1_VS_$base2";
';
				
# consistency original reps SRR		
find . -name "*rep1*VS*input*regionPeak.gz" -not -name "*\.pr[1-2]\.*" | xargs -n 1 -P 16 -iFILES sh -c '
	rep2=$(echo FILES | sed "s/rep1/rep2/g"); 
	rep3=$(echo FILES | sed "s/rep1/rep3/g");
	base=$(echo FILES | sed "s/rep1.*/rep/g");
	Rscript batch-consistency-analysis.r FILES $rep2 -1 "$base1_VS_$base2" 0 F signal.value &
	Rscript batch-consistency-analysis.r FILES $rep3 -1 "$base1_VS_$base3" 0 F signal.value &
	Rscript batch-consistency-analysis.r $rep2 $rep3 -1 "$base2_VS_$base3" 0 F signal.value;
	Rscript batch-consistency-plot.r 1 "$base1_VS_$base2" "$base1_VS_$base2" & 
	Rscript batch-consistency-plot.r 1 "$base1_VS_$base3" "$base1_VS_$base3" &
	Rscript batch-consistency-plot.r 1 "$base2_VS_$base3" "$base2_VS_$base3" &
	Rscript batch-consistency-plot.r 3 "$base1_VS_$base2" "$base1_VS_$base3" "$base2_VS_$base3"; 
'	

# consistency original reps MACS2
find . -name "*rep1*VS*input*narrowPeak.sorted" -not -name "*\.pr[1-2]\.*" | xargs -n 1 -P 16 -iFILES sh -c '
	rep2=$(echo FILES | sed "s/rep1/rep2/g"); 
	rep3=$(echo FILES | sed "s/rep1/rep3/g");
	base=$(echo FILES | sed "s/rep1.*/rep/g");
	Rscript batch-consistency-analysis.r FILES $rep2 -1 "$base1_VS_$base2".macs 0 F signal.value &
	Rscript batch-consistency-analysis.r FILES $rep3 -1 "$base1_VS_$base3".macs 0 F signal.value &
	Rscript batch-consistency-analysis.r $rep2 $rep3 -1 "$base2_VS_$base3".macs 0 F signal.value;
	Rscript batch-consistency-plot.r 1 "$base1_VS_$base2".macs "$base1_VS_$base2".macs & 
	Rscript batch-consistency-plot.r 1 "$base1_VS_$base3".macs "$base1_VS_$base3".macs &
	Rscript batch-consistency-plot.r 1 "$base2_VS_$base3".macs "$base2_VS_$base3".macs &
	Rscript batch-consistency-plot.r 3 "$base1_VS_$base2".macs "$base1_VS_$base3".macs "$base2_VS_$base3".macs; 
'

# pooled optimal / conservative spp 
find . -name "*rep0*pr1*VS*rep0*pr2*overlap*" -not -name "*macs*" | xargs -n 1 -P $threads -iFILES sh -c '
 hold=FILES;
 optimal=$(wc -l FILES);
 conservative=0;
 export base=$(echo FILE | sed "s/.rep0.*//g");
 for f in $base*overlap*; do
  if [[ "$f" != "$hold" ]]; then
   $comp=$(wc -l $f);
   if [ $comp -gt $optimal ]; then
    $optimal=$comp;	
   fi;
   if [ $comp -gt $consevative ]; then
    $conservative=$comp;	
   fi;
  fi;
 done;
 peakfile=$(find . -name "$base.rep0*VS*input*regionPeak.gz" -not -name "*pr[1-2]*");
 cat "$peakfile" | sort -k 8nr,8nr | head -n $optimal | gzip -c > "$peakfile.spp.optimal.gz";  
 cat "$peakfile" | sort -k 8nr,8nr | head -n $conservative | gzip -c > "$peakfile.spp.conservative.gz";  
'
	
# pooled optimal / conservative macs2 
find . -name "*rep0*pr1*VS*rep0*pr2*macs*overlap*" | xargs -n 1 -P $threads -iFILES sh -c '
 hold=FILES;
 optimal=$(wc -l FILES);
 conservative=0;
 export base=$(echo FILE | sed "s/.rep0.*//g");
 for f in $base*overlap*; do
  if [[ "$f" != "$hold" ]]; then
   $comp=$(wc -l $f);
   if [ $comp -gt $optimal ]; then
    $optimal=$comp;	
   fi;
   if [ $comp -gt $consevative ]; then
    $conservative=$comp;	
   fi;
  fi;
 done;
 peakfile=$(find . -name "$base.rep0*VS*input*narrowPeak" -not -name "*pr[1-2]*");
 cat "$peakfile" | sort -k 8nr,8nr | head -n $optimal | gzip -c > "$peakfile.macs.optimal.gz";  
 cat "$peakfile" | sort -k 8nr,8nr | head -n $conservative | gzip -c > "$peakfile.macs.conservative.gz";  
'

exit 0
