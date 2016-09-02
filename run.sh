#!/bin/bash
#ENCODE SPP, MACS2, IDR pipeline
#usage ./run.sh genome_fasta raw_reads_directory output_dir macs2_genome_size(e.g., hs/gs/etc) chromo_size_path max_num_replicates

threads=$(( $(grep -c ^processor /proc/cpuinfo) * 2 ));
genome=$1;
export reads_dir=$2;
export out_dir=$3;
export gsize=$4;
export csize=$5;
export reps=$6;

rm -rf encode_temp_dir;
rm -rf "$out_dir";
mkdir encode_temp_dir;
mkdir "$out_dir";

#alignments and sam file generation - multithreaded aln, samse is memory heavy
for f in "$reads_dir"/*.fastq; do
 temp_loc=$(echo $f | sed "s/$reads_dir\///"); 
 bwa aln -I -B 0 -t $threads $genome $f | bwa samse $genome - $f > encode_temp_dir/"$temp_loc".sam;
done;

#q30 filter - single core
find encode_temp_dir/ -name "*sam" | xargs -n 1 -P $threads -iFILES bash -c '
 samtools view -F 1548 -h -S -q 30 -o FILES_q30.sam FILES;
';

#generate tag densities for q30 - single core
find encode_temp_dir/ -name "*_q30.sam" | xargs -n 1 -P $threads -iFILES bash -c '
 samtools view -Sb FILES | bedtools bamtobed -i stdin | awk '\''BEGIN{FS="\t";OFS="\t"}{$4="N";print$0}'\'' | gzip -c > FILES.tagAlign.gz;
';

#generate bams - single core
find encode_temp_dir/ -name "*_q30.sam" | xargs -n 1 -P $threads -iFILES bash -c '
 out="$out_dir"/$(echo FILES | sed -e "s/encode_temp_dir.//g" | sed "s/\.sam//g")
 samtools view -Sbh FILES > FILES.bam;
 samtools sort FILES.bam "$out".sort;
 samtools index "$out"sort.bam "$out"sort.bam.bai;
 perl bam2wiggle.perl -m -D "$out_dir" "$out".sort.bam;
 genomeCoverageBed -split -bg -ibam "$out".sort.bam -g "$csize" > FILES.bedgraph;
 sort -k1,1 -k2,2n FILES.bedGraph > FILES.sort.bedGraph;
 bedGraphToBigWig FILES.sort.bedGraph "$csize" FILES.bigWig;
 bigWigToWig FILES.bigWig "$out".wig
';

#generate pooled tags from biological replicates - single core
find encode_temp_dir/ -name "*rep1*.tagAlign.gz" | xargs -n 1 -P $threads -iFILES bash -c '
 pool=$(echo FILES | sed "s/rep1/rep0/");
 stringer="";
 for ((x=1; x<=$reps; x++)); do
  rep=$(echo FILES | sed "s/rep1/rep$x/g");
  stringer+="$rep"" ";
 done;
 zcat $stringer | gzip -c > $pool;
';

#shuffle tags and generate psuedo replicates - single core
find encode_temp_dir/ -name "*tagAlign.gz"  -not -name "*input*" | xargs -n 1 -P $threads -iFILES bash -c '
 nline=$(( ($(zcat FILES | wc -l) + 1) / 2 )); 
 zcat FILES | shuf | split -d -l $nline - ./FILES.0; 
 gzip ./FILES.000 &
 gzip ./FILES.001;
 mv ./FILES.000.gz ./FILES.pr1.gz & 
 mv ./FILES.001.gz ./FILES.pr2.gz &
';

#spp peak calls - multithreaded through R snow package support
for f in $(find . -name "*tagAlign*gz"  -not -name "*input*"); do
 inp="$(echo $f | sed "s/rep./input\.rep0/" | sed "s/\.pr[1-9]\.gz//g")"; 
 base="$out_dir"/"$(echo $f | sed "s/..encode\_temp\_dir\///g" | sed "s/.fastq.sam_q30.sam.tagAlign//" | sed "s/\.gz//g")";
 echo "$base".regionPeak;
 echo "$base".plot.pdf;
 Rscript run_spp.R -c=$f -i=$inp -npeak=30000 -x=-500:85 -odir=./ -p=10 -savr -savp -rf -out="$f".cc -savr="$base".regionPeak -savp="$base".plot.pdf
 cp "$out_dir"/*.regionPeak.gz encode_temp_dir/;
done;
find . -name "*regionPeak.gz" | xargs -n 1 -P $threads -iFILES bash -c 'gunzip FILES;';

# macs2 peak calls - macs2 single core
find encode_temp_dir/ -name "*tagAlign*gz"  -not -name "*input*" | xargs -n 1 -P $threads -iFILES bash -c '
 inp=$(echo FILES | sed "s/rep./input\.rep0/" | sed "s/.pr[1-2].gz//g"); 
 base=$(echo FILES | sed "s/^.*\///g" | sed "s/.fastq.sam_q30.sam.tagAlign.gz//g" | sed "s/\.gz//g");
 ss=$(cat FILES.cc | cut -f3 | head -n 1 | cut -d "," -f1); 
 macs2 callpeak -t FILES -c $inp -g $gsize --to-large --nomodel --extsize "$ss" -p 0.1 -n "$base" --outdir encode_temp_dir/ --tempdir encode_temp_dir/;
';

# sort MACS2 peaks - single core
find encode_temp_dir/ -name "*.narrowPeak" | xargs -n 1 -P $threads -iFILES bash -c '
 base="$out_dir"/"$(echo FILES | sed "s/^.*\///g")";
 sort -k 8nr,8nr FILES | head -n 30000 | gzip -c > FILES.gz;
 rename "s/\_peaks\.narrowPeak\.gz/\.narrowPeak\.gz/g" encode_temp_dir/*peaks.narrowPeak.gz;
 cp encode_temp_dir/*narrowPeak.gz "$out_dir";
';

# comparisons for IDR analysis, single core
rm -rf encode_temp_dir/comparisons.txt encode_temp_dir/stringer.txt;
find encode_temp_dir/ -name "*rep1*Peak*gz" -not -name "*\.pr[1-2]\.*" | xargs -n 1 -P 16 -iFILES bash -c '
 stringer="";	 
 sig="signal.value"; 
 if [[ FILES =~ .*narrow.* ]]; then 
  sig="p.value"; 
 fi;
 for ((x=0; x<=$reps; x++)); do
  pr1=$(echo FILES | sed "s/rep1\./rep$x\.pr1\./g");
  pr2=$(echo FILES | sed "s/rep1\./rep$x\.pr2\./g");
  base_pr2=$(echo $pr2 | sed "s/^.*\///g");
  echo -ne "$pr1\t$pr2\t$base_pr2\t$sig\n" >> encode_temp_dir/comparisons.txt;
  if [ $x -lt $reps ] && [ $x -gt 0 ]; then
   z=$(( $x + 1 ));
   for ((y=$z; y<=$reps; y++)); do
    f1=$(echo FILES | sed "s/rep1/rep$x/g");
    f2=$(echo FILES | sed "s/rep1/rep$y/g"); 
    base_f1=$(echo $f1 | sed "s/^.*\///g"); 
    base_f2=$(echo $f2 | sed "s/^.*\///g"); 
    stringer+="$base_f1""_VS_""$base_f2""~~~";
    echo -ne "$f1\t$f2\t$base_f2\t$sig\n" >> encode_temp_dir/comparisons.txt;
   done;
  fi; 
 done;
 echo -ne "$stringer\n" | sed "s/~~~$//" >> encode_temp_dir/stringer.txt;
'

#consistency analysis
cat encode_temp_dir/comparisons.txt | xargs -n 1 -P 16 -iLINES bash -c '
 pr1=$(echo LINES | cut -d " " -f 1); 
 pr2=$(echo LINES | cut -d " " -f 2); 
 base_pr2=$(echo LINES | cut -d " " -f 3); 
 sig=$(echo LINES | cut -d " " -f 4);
 Rscript batch-consistency-analysis.r "$pr1" "$pr2" -1 "$pr1"_VS_"$base_pr2" 0 F "$sig";
 Rscript batch-consistency-plot.r 1 "$pr1"_VS_"$base_pr2" "$pr1"_VS_"$base_pr2";
'

# pooled optimal / conservative macs2 - single core
find encode_temp_dir/ -name "*rep0*pr1*VS*rep0*pr2*overlap*" | xargs -n 1 -P $threads -iFILES bash -c '
 hold=FILES;
 type="region"; 
 if [[ "$hold" =~ .*narrow.* ]]; then 
  type="narrow"; 
 fi;
 optimal=$(wc -l FILES | cut -d " " -f 1);
 conservative=0;
 base=$(echo FILES | sed "s/.rep0.*//g" | sed "s/^.*\///g");
 for f in encode_temp_dir/"$base"*"$type"Peak*overlap*; do
  if [[ "$f" != "$hold" ]]; then 
   comp=$(wc -l $f | cut -d " " -f 1); 
   if [ "$comp" -gt "$optimal" ]; then
    optimal=$comp;	
   fi;
   if [ "$comp" -gt "$conservative" ]; then
    conservative=$comp;	
   fi;
  fi;
 done;
 peakfile=$(find encode_temp_dir/ -name "$base*rep0*$typePeak.gz" -not -name "*pr[1-2]*" -not -name "*VS*");
 zcat "$peakfile" | sort -k 8nr,8nr | head -n $optimal | gzip -c > "$peakfile".optimal.gz;  
 zcat "$peakfile" | sort -k 8nr,8nr | head -n $conservative | gzip -c > "$peakfile".conservative.gz;  
'

cd encode_temp_dir/;
cp *conservative* "$out_dir";
cp *optimal* "$out_dir";
cp *plot* "$out_dir";
cp *overlap* "$out_dir";

echo "Finished!";

exit 0;
