#!/bin/bash

# Check if enough arguments are passed
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_dir>"
    echo "Please provide path to input directory"
    exit 1
fi

# Parse arguments
input_dir="${1%/}"  # First argument: Input directory containing sample subdirectories

# Ensure the input directory exists
if [ ! -d "$input_dir" ]; then
    echo "$input_dir is not a valid directory"
    exit 1
fi

# Loop through all directories inside the input directory
for dir in "$input_dir"/*/; do

    clean_dir="${dir%/}"
    echo "Processing directory: $dir"
    dir_name=$(basename "$dir")

    echo "Processing sample in directory: $dir"

    # Create the lambda_alignment directory if it doesn't exist
    mkdir -p $dir/lambda_alignment
    mkdir -p $dir/unaligned_output

    cat $clean_dir/*fastq > $clean_dir/merged.fastq
    echo "fastq merged into one"
    # Call the lambda_alignment_summary.sh script with correct arguments
    ./lambda_alignment_summary.sh \
        $dir/merged.fastq \
        /mnt/d/READ_DM1/DATA/REFERENCE/lambda.fasta \
        "$clean_dir/lambda_alignment" \
        "$dir_name"

    
    # Convert unaligned BAM to FASTQ and save it
     
    # Use an array to hold all unaligned BAM files
    unaligned_bam_files=$(ls $clean_dir/lambda_alignment/*_unaligned.bam | head -n 1)
    #echo "$unaligned_bam_files"
    echo "============================================"
    # Check if there are any unaligned BAM files found
    if [ ${#unaligned_bam_files[@]} -gt 0 ]; then
        samtools fastq "${unaligned_bam_files[@]}" > "$dir/unaligned_output/${dir_name}.fastq"
    else
        echo "No unaligned BAM files found in $dir/lambda_alignment."
    fi
    # Run the prepare script with specified parameters
    echo "Starting STRAT_prepare..."
    python3 ./strat_prepare.py \
        --config ./config.yaml \
        --input_path "$clean_dir/unaligned_output" \
        --output_path "$clean_dir/unaligned_output"
    
    #clean_dir="${dir%/}" 
    # Merge all .fastq.ontarget.tsv files in the current directory
    if ls "$clean_dir"/*.fastq.ontarget.tsv 1> /dev/null 2>&1; then
        cat "$clean_dir"/*.fastq.ontarget.tsv > "$clean_dir/merged.ontarget.tsv"
        rm "$clean_dir"/*.fastq.ontarget.tsv
        echo "Merged .ontarget.tsv file created: $clean_dir/merged.ontarget.tsv"
    else
        echo "No .fastq.ontarget.tsv files found in $clean_dir"
    fi

    # Call STRAT process
    echo "Starting STRAT process..."
    python3 strat_process.py \
        --motif_prim "CAG" \
        --motif_scnd "CAG" \
        --threshold 1 \
        --input_path "$clean_dir/merged.ontarget.tsv" \
        --output_path "$dir"

    # Run fastq2tsv
    python3 ./fastq2tsv.py \
        --fastq_path "$clean_dir" \
        --output_path "$clean_dir/fastq.tsv"

    # Call plot.py
    python3 ./plots.py \
        --fastq_tsv_path "$clean_dir" \
        --merged_ontarget_path "$clean_dir" \
        --raw_ontarget_merged_path "$clean_dir" \
        --images_path "$clean_dir" \
        --processed_path "$clean_dir"

    # Call TSV summary script (proba.py)
    echo "Running TSV summary..."
    python3 ./summarize.py \
        --lambda_alignment_summary "$clean_dir/lambda_alignment/lambda_alignment_summary.tsv" \
        --strat_process_output "$clean_dir/merged.ontarget.processed.tsv" \
        --output_dir_path "$dir"

done
