# ENCODE TF ChIP pipeline using SPP and MACS2 for peak calling, and IDR for consistency analysis

Includes two major scripts:

1) install_software.sh - if the VM is new will install the needed software for analysis minus the genome index

2) run.sh - main pipeline for alignments, peak calling, and IDR analysis. 

usage ./run.sh genome_fasta raw_reads_directory output_dir macs2_genome_size(e.g., hs/gs/etc) chromo_size_path max_num_replicates

e.g. nohup ./run.sh /mnt/genome/hg19.fa /mnt/reads /mnt/output hs /mnt/genome/chromo.size.txt 3 &
