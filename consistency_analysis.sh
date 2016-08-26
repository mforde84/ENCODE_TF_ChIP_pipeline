#TODO


# self-psuedo - pr1 pr2

find . -name "*_r*pr1*regionPeak" | xargs -n 1 -P 16 -iFILES sh -c '
	f=FILES; 
	pr1=$(echo $f | sed '\''s/\.\///g'\'' | sed '\''s/\.regionPeak//g'\''); 
	pr2=$(echo $pr1 | sed '\''s/pr1/pr2/g'\''); 
	Rscript batch-consistency-analysis.r "$pr1.regionPeak" "$pr2.regionPeak" -1 $pr1"_VS_"$pr2 0 F signal.value; 
	Rscript batch-consistency-plot.r 1 $pr1"_VS_"$pr2 $pr1"_VS_"$pr2; numPeaks_Rep2_pr=$(awk '\''$11 <= 0.01 {print $0}'\'' $pr1"_VS_"$pr2-overlapped-peaks.txt | wc -l); 
	echo "MSG: $pr1 vs $pr1 $numPeaks_Rep2_pr" >> num_peaks.txt; 
	sort -k 8nr,8nr "$pr1.peaks.narrowPeak" | head -n 30000 > "$pr1.narrowPeak.sorted"; 
	sort -k 8nr,8nr "$pr2.peaks.narrowPeak" | head -n 30000 > "$pr2.narrowPeak.sorted"; 
	Rscript batch-consistency-analysis.r "$pr1.narrowPeak.sorted" "$pr2.narrowPeak.sorted" -1 $pr1"_VS_"$pr2.macs 0 F p.value; 
	Rscript batch-consistency-plot.r 1 $pr1"_VS_"$pr2.macs $pr1"_VS_"$pr2.macs; 
	numPeaks_Rep2_pr_macs2=$(awk '\''$11 <= 0.01 {print $0}'\'' $pr1"_VS_"$pr2".macs-overlapped-peaks.txt" | wc -l); 
'
			
			
# consistency rep1 and rep 2s
		
find . -name "*_r1.regionPeak" | xargs -n 1 -P 16 -iFILES sh -c '
	f=FILES;
	pr1=$(echo $f | sed '\''s/\.\///g'\'' | sed '\''s/.regionPeak//g'\''); 
	pr2=$(echo $pr1 | sed '\''s/_r1/_r2/g'\''); 
	Rscript batch-consistency-analysis.r "$pr1.regionPeak" "$pr2.regionPeak" -1 "$pr1""_VS_""$pr2" 0 F signal.value; 
	numPeaks_Rep2_Rep1=$(awk '\''$11 <= 0.02 {print $0}'\'' "$pr1""_VS_""$pr2-overlapped-peaks.txt" | wc -l); 
	if [ $numPeaks_Rep2_Rep1 -gt $max_numPeaks_Rep ]; then 
	 max_numPeaks_Rep=$numPeaks_Rep2_Rep1; 
	fi;
	sort -k 8nr,8nr "$pr1.peaks.narrowPeak" | head -n 30000 > "$pr1.peaks.narrowPeak.sorted"; 
	sort -k 8nr,8nr "$pr2.peaks.narrowPeak" | head -n 30000 > "$pr2peaks.narrowPeak.sorted"; 
	Rscript batch-consistency-analysis.r "$pr1.peaks.narrowPeak.sorted" "$pr2.peaks.narrowPeak.sorted" -1 "$pr1""_VS_""$pr2.macs" 0 F p.value; 
	numPeaks_Rep2_Rep1_macs2=$(awk '\''$11 <= 0.02 {print $0}'\'' "$pr1""_VS_""$pr2.macs-overlapped-peaks.txt" | wc -l); 
	if [ $numPeaks_Rep2_Rep1_macs2 -gt $max_numPeaks_Rep_macs2 ]; then 
	 max_numPeaks_Rep_macs2=$numPeaks_Rep2_Rep1_macs2; 
	fi;
'

# consistency: pool pooled-pseduoreplicates

find . -name "*pool.regionPeak.gz" | xargs -n 1 -P 16 -iFILES sh -c '
	f=FILES; pool=$(echo $f | sed '\''s/\.\///g'\'' | sed 's/\.regionPeak\.gz//'); 
	Rscript batch-consistency-analysis.r "$pool"".pr1.regionPeak.gz" "$pool"".pr2.regionPeak.gz" -1 "$pool"".pr1_VS_""$pool"".pr2" 0 F signal.value; 
	Rscript batch-consistency-plot.r 1 "$pool"".pr1_VS_""$pool"".pr2" "$pool"".pr1_VS_""$pool"".pr2"; 
	numPeaks_Rep0_pr=$(awk '\''$11 <= 0.01 {print $0}'\'' "$pool"".pr1_VS_""$pool"".pr2-overlapped-peaks.txt" | wc -l); 
	echo "MSG: numPeaks_Rep0_pr $numPeaks_Rep0_pr"; 
	if [ $numPeaks_Rep0_pr -gt $max_numPeaks_Rep ]; then 
	 optThresh=$numPeaks_Rep0_pr; 
	else 
	 optThresh=$max_numPeaks_Rep; 
	fi; 
	echo "MSG: optThresh $optThresh"; zcat "$pool"".regionPeak.gz" | sort -k7nr,7nr | head -n $optThresh | gzip -c > "spp.optimal.""$pool"".regionPeak.gz"; 
	echo "MSG: conThresh $max_numPeaks_Rep"; zcat "$pool"".regionPeak.gz" | sort -k7nr,7nr | head -n 3032 | gzip -c > "spp.conservative.""$pool"".regionPeak.gz"; 
	sort -k 8nr,8nr "$pool"".pr1.peaks.narrowPeak" | head -n 30000 > "$pool"".pr1.peaks.narrowPeak.sorted"; 
	sort -k 8nr,8nr "$pool"".pr2.peaks.narrowPeak" | head -n 30000 > "$pool"".pr2.peaks.narrowPeak.sorted"; 
	Rscript batch-consistency-analysis.r "$pool"".pr1.peaks.narrowPeak.sorted" "$pool"".pr2.peaks.narrowPeak.sorted" -1 "$pool"".pr1_VS_""$pool"".pr2.macs" 0 F p.value; 
	Rscript batch-consistency-plot.r 1 "$pool"".pr1_VS_""$pool"".pr2.macs" "$pool"".pr1_VS_""$pool"".pr2.macs"; 
	numPeaks_Rep0_pr_macs2=$(awk '\''$11 <= 0.01 {print $0}'\'' "$pool"".pr1_VS_""$pool"".pr2.macs-overlapped-peaks.txt" | wc -l); 
	if [ $numPeaks_Rep0_pr_macs2 -gt $max_numPeaks_Rep ]; then 
	 optThresh_macs2=$numPeaks_Rep0_pr_macs2; 
	else 
	 optThresh_macs2=$max_numPeaks_Rep_macs2; 
	fi; 
	cat "$pool"".peaks.narrowPeak" | sort -k 8nr,8nr | head -n $optThresh_macs2 | gzip -c > "macs.optimal.""$pool"".peaks.narrowPeak.gz"; 
	cat "$pool"".peaks.narrowPeak" | sort -k 8nr,8nr | head -n $max_numPeaks_Rep_macs2 | gzip -c > "macs.conservative.""$pool"".peaks.narrowPeak.gz";
'
