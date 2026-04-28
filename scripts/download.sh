#!/bin/bash  
#SBATCH --mail-user=isacaste@iu.edu  
#SBATCH --nodes=2  
#SBATCH --mem=200G  
#SBATCH -p gpu  
#SBATCH --ntasks-per-node=2  
#SBATCH --gpus-per-node=2  
#SBATCH --time=1-23:59:00  
#SBATCH --mail-type=BEGIN,FAIL,END  
#SBATCH --job-name=fastq 
#SBATCH -o fastq.out  
#SBATCH -e fastq.err 
#SBATCH -A r00750  

module load sra-toolkit

SAMPLES=(
    SRR31972716
    SRR31972717
    SRR31972718
    SRR31972719
    SRR31972720
    SRR31972721
    SRR31972722
    SRR31972723
    SRR31972724
    SRR31972725
    SRR31972726
    SRR31972727
    SRR31972728
    SRR31972729
    SRR31972730
    SRR31972731
    SRR31972732
    SRR31972733
    SRR31972734
    SRR31972735
    SRR31972736
    SRR31972737
    SRR31972738
    SRR31972739
    SRR31972740
    SRR31972741
    SRR31972742
    SRR31972743
    SRR31972744
    SRR31972745
    SRR31972746
    SRR31972747
    SRR31972748
    SRR31972749
    SRR31972750
    SRR31972751
    SRR31972752
    SRR31972753
)

for sample in "${SAMPLES[@]}"; do
    echo "Downloading $sample..."
    fasterq-dump --split-files --gzip $sample
    echo "$sample done."
done

echo "All downloads complete."

