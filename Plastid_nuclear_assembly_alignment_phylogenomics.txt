# ---- 1. Plastid CDS data, extracted with PhyloSuite

# Assemble plastomes for samples US31 and US32

# Index the reference
bwa index /Data/cbarrett/004_2024_MH_EpiBase_seqcap/getorg_seeds/Stereosandra_javanica.fasta

# Output folder
mkdir -p mapped_bams

# Define input arrays
samples=("OS2" "US31" "US32")

for sample in "${samples[@]}"; do
  R1="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/fastp_polyg_merged/fastp-${sample}_Stereosandra_javanica_merged_R1.fastq.gz"
  R2="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/fastp_polyg_merged/fastp-${sample}_Stereosandra_javanica_merged_R2.fastq.gz"
  REF="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/getorg_seeds/Stereosandra_javanica.fasta"
  OUT="mapped_bams/${sample}.bam"

  bwa mem -t 8 "$REF" "$R1" "$R2" | samtools view -bS - | samtools sort -o "$OUT"
  samtools index "$OUT"
done

for bam in mapped_bams/*.bam; do
  sample=$(basename "$bam" .bam)
  echo "== $sample =="

  # 1. Total mapped reads
  total_reads=$(samtools view -c "$bam")
  mapped_reads=$(samtools view -c -F 4 "$bam")
  percent_mapped=$(echo "scale=2; 100 * $mapped_reads / $total_reads" | bc)
  echo "Total reads: $total_reads"
  echo "Mapped reads: $mapped_reads"
  echo "Percent mapped: $percent_mapped%"

  # 2. Mean coverage
  mean_cov=$(samtools depth "$bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
  echo "Mean coverage depth: $mean_cov"
  echo ""
done

mkdir -p consensus

for bam in mapped_bams/*.bam; do
  sample=$(basename "$bam" .bam)
  REF="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/getorg_seeds/Stereosandra_javanica.fasta"

  # Generate BCF
  bcftools mpileup -f "$REF" "$bam" | bcftools call -mv -Oz -o "consensus/${sample}.vcf.gz"
  bcftools index "consensus/${sample}.vcf.gz"

  # Generate consensus FASTA
  bcftools consensus -f "$REF" "consensus/${sample}.vcf.gz" > "consensus/${sample}_consensus.fasta"
done
# First align with mafft

# Directories
INDIR="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/2025_June_Nervilieae/CDS_NUC"
OUTDIR_MAFFT="$INDIR/mafft"
mkdir -p "$OUTDIR_MAFFT"

# MAFFT alignment using 80 threads
for file in "$INDIR"/*.fas; do
    base=$(basename "$file" .fas)
    mafft --thread 80 "$file" > "$OUTDIR_MAFFT/${base}_mafft.fas"
done


# Codon align with macse2

# Paths
MAFFT_DIR="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/2025_June_Nervilieae/CDS_NUC/mafft"
MACSE_OUT="/Data/cbarrett/004_2024_MH_EpiBase_seqcap/2025_June_Nervilieae/CDS_NUC/macse"
MACSE_JAR="/usr/local/bin/macse_v2.07.jar"   # <- update this path!

mkdir -p "$MACSE_OUT"

for aln in "$MAFFT_DIR"/*_mafft.fas; do
    base=$(basename "$aln" _mafft.fas)
    
    java -Xmx32g -jar "$MACSE_JAR" \
        -prog alignSequences \
        -seq "$aln" \
        -out_NT "$MACSE_OUT/${base}_macse_NT.fasta" \
        -out_AA "$MACSE_OUT/${base}_macse_AA.fasta"
done



# ---- 2. Nuclear (A353) data

# Download additional read files from SRA
fasterq-dump SRR26934839 SRR27827645 SRR6008325 --split-files -e 50

## run fastp

for f1 in *_R1.fastq.gz
	do
        f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
        fastp -i $f1 -I $f2 -w 16 --trim_poly_g --trim_poly_x -l 25 --cut_right -o "../fastp_polyg/fastp-$f1" -O "../fastp_polyg/fastp-$f2"
	done


# FastQC/MultiQC

conda activate /usr/local/src/conda_envs/fastqc

# Define input and output directories
dirs=("reads" "fastp_polyg")
for dir in "${dirs[@]}"; do
    # Create output directory
    outdir="fastqc_${dir}"
    mkdir -p "$outdir"

    # Run FastQC on all .fastq.gz files in the directory
    fastqc "${dir}"/*.fastq.gz -o "$outdir" --threads 8
done

# Run MultiQC on each FastQC output folder
for dir in "${dirs[@]}"; do
    multiqc "fastqc_${dir}" -o "fastqc_${dir}"
done



conda activate /usr/local/src/conda_envs/hybpiper

mkdir hybpiper_2025_06_03
sudo chmod 777 -R hybpiper_2025_06_03
cd hybpiper_2025_06_03

# Hybpiper loop
while read name; 
do hybpiper assemble -t_dna ../orchids_targetfile.fasta -r ../Nervilieae_keepers/$name*.fastq.gz --prefix $name --bwa; 
done < ../namelist.txt

hybpiper stats -t_dna ../orchids_targetfile.fasta gene ../namelist.txt

# Extract sequences
hybpiper retrieve_sequences -t_dna ../orchids_targetfile.fasta dna --sample_names ../namelist.txt
hybpiper paralog_retriever ../namelist.txt -t_dna ../orchids_targetfile.fasta dna


# Align loci with MAFFT

for i in *.FNA; do
mafft --adjustdirection --thread 32 ${i} > ${i%.*}_mafft.fasta;
done

# Loop to run mixturefinder for all alignments
conda activate /usr/local/src/conda_envs/binf

# Define directories
input_dir=/Data/cbarrett/004_2024_MH_EpiBase_seqcap/2025_June_Nervilieae/mafft  # Directory with alignment files
output_dir="mix_loop_out_nervilieae"     # Directory to save IQ-TREE outputs
mkdir -p "$output_dir"            # Create output directory if it doesn't exist

# Initialize a results file to store the best model and BIC scores
bic_file="$output_dir/best_model_bic_summary.txt"
echo -e "Alignment\tBest_Model\tBIC_Score" > "$bic_file"  # Header for BIC scores summary

# Function to run ModelFinder for each alignment and capture BIC score
run_modelfinder_bic() {
    alignment="$1"
    base=$(basename "$alignment" .fas)
    
	# run for 1-4 mixture classes
    iqtree2 -s "$alignment" -m MIX+MFP -cmax 4 -B 1000 --prefix "$output_dir/${base}_model_selection" -T 2 > /dev/null 2>&1

	best_model=$(grep "Best-fit model according to BIC:" "$output_dir/${base}_model_selection.iqtree" | awk '{print $6}')
    bic_score=$(grep "BIC score:" "$output_dir/${base}_model_selection.iqtree" | awk '{print $3}')
    echo -e "${base}\t${best_model}\t${bic_score}" >> "$bic_file"
}

# Export functions and variables for GNU parallel
export -f run_modelfinder_bic
export output_dir model_set bic_file

# Step #2: Run ModelFinder in parallel for each alignment
echo "Starting ModelFinder for each alignment..."
find "$input_dir" -name "*.fna" | parallel -j 45 run_modelfinder_bic

# Wait for all ModelFinder processes to complete
wait  # Ensures all parallel jobs are completed before moving to the next steps
echo "ModelFinder analysis complete for all alignments. Summary of best models and BIC scores saved to $bic_file."

### Sum the BIC scores and # parameters across alignments
### Pull values and create a table: alignment, NP, BIC, logL


#!/bin/bash

# Define output file for the results table
output_file="free_parameters_bic_logl_summary.tsv"
echo -e "Filename\tFree_Parameters\tLog_Likelihood\tBIC_Score" > "$output_file"

# Loop through each .iqtree file in the directory
for file in mix_loop_out_nervilieae/*.iqtree; do
  echo "Processing file: $file"  # Debugging output to track file processing

  # Initialize variables to store extracted data
  free_params="NA"
  log_likelihood="NA"
  bic_score="NA"

  # Parse the .iqtree file line by line
  while IFS= read -r line; do
    # Extract the number of free parameters
    if [[ "$line" =~ "Number of free parameters" ]]; then
      free_params=$(echo "$line" | awk '{print $NF}')
    # Extract the log-likelihood (capture the negative or positive number only)
    elif [[ "$line" =~ "Log-likelihood of the tree:" ]]; then
      log_likelihood=$(echo "$line" | awk '{for(i=1;i<=NF;i++) if ($i ~ /^-?[0-9]+\.[0-9]+$/) print $i}')
    # Extract the BIC score
    elif [[ "$line" =~ "Bayesian information criterion (BIC) score:" ]]; then
      bic_score=$(echo "$line" | awk '{print $NF}')
    fi
  done < "$file"

  # Check if any fields are still NA, indicating a parsing failure
  if [[ "$free_params" == "NA" || "$log_likelihood" == "NA" || "$bic_score" == "NA" ]]; then
    echo "Warning: Missing data for $file. Manual review may be needed."  # Debugging output
  fi

  # Append extracted data to the output file
  echo -e "$(basename "$file")\t$free_params\t$log_likelihood\t$bic_score" >> "$output_file"
done



# Concatenate gene trees and run wastral

cd /Data/cbarrett/004_2024_MH_EpiBase_seqcap/2025_June_Nervilieae/mix_loop_out_nervilieae

cat *.contree > mixfinder_allparalogs.treefile

wastral --root fastp-Tropidia_curculigoides -t 32 -i mixfinder_allparalogs.treefile -o mixfinder_wastral_nervilieae.tre

# Verbose output for plotting pie charts at nodes

wastral --root fastp-Tropidia_curculigoides -t 32 -u 2 -i mixfinder_allparalogs.treefile -o mixfinder_wastral_nervilieae_verbose.tre



# Missing data

# Full plastid data with Epipogium
total=$(seqkit seq -i Nervilieae_pt_concat.fasta | grep -v '^>' | tr -d '\n' | wc -c)
missing=$(seqkit seq -i Nervilieae_pt_concat.fasta | grep -v '^>' | tr -cd '\-\?N' | wc -c)

pct=$(echo "scale=2; 100 * $missing / $total" | bc)
echo "Missing data percentage: $pct%"

# Full plastid data withOUT Epipogium
total=$(seqkit seq -i Nervilieae_pt_concat_no_Epipogium.fasta | grep -v '^>' | tr -d '\n' | wc -c)
missing=$(seqkit seq -i Nervilieae_pt_concat_no_Epipogium.fasta | grep -v '^>' | tr -cd '\-\?N' | wc -c)

pct=$(echo "scale=2; 100 * $missing / $total" | bc)
echo "Missing data percentage: $pct%"

# Convert wastral coalescent tree from coalescent BL units to substitutions/site

astral4 -i <gene-tree-path> -C -c <species-tree-path> -o <output-path>

astral4 -i mixfinder_allparalogs.treefile -C -c mixfinder_wastral_nervilieae.tre -o Nervilieae_mega353_castles_bl.tre --genelength 1000 --root fastp-Tropidia_curculigoides

