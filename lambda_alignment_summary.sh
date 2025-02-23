#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <fastq_file> <reference_genome> <output_path> <sample_name>"
    exit 1
fi

# Assign command-line arguments to variables
FASTQ_FILE=$1
REFERENCE_GENOME=$2
OUTPUT_FOLDER=$3
SAMPLE_NAME=$4

# Define output BAM file names with sample name and correct paths
OUTPUT_BAM="$OUTPUT_FOLDER/${SAMPLE_NAME}.bam"
ALIGNED_BAM="$OUTPUT_FOLDER/${SAMPLE_NAME}_aligned.bam"
UNALIGNED_BAM="$OUTPUT_FOLDER/${SAMPLE_NAME}_unaligned.bam"
SORTED_ALIGNED_BAM="$OUTPUT_FOLDER/${SAMPLE_NAME}_sorted_aligned.bam"
SORTED_UNALIGNED_BAM="$OUTPUT_FOLDER/${SAMPLE_NAME}_sorted_unaligned.bam"
echo "==================================="
echo "$OUTPUT_BAM"
echo "====================================="
# Align the FASTQ file to the reference genome using minimap2
minimap2 -ax map-ont "$REFERENCE_GENOME" "$FASTQ_FILE" | samtools view -Sb - > "$OUTPUT_BAM"

# Separate aligned and unaligned reads
samtools view -b -f 4 "$OUTPUT_BAM" > "$UNALIGNED_BAM"  # Unaligned reads
samtools view -b -F 4 "$OUTPUT_BAM" > "$ALIGNED_BAM"    # Aligned reads

# Sort the aligned BAM file
# SORTED_ALIGNED_BAM="${ALIGNED_BAM%.bam}_sorted.bam"
samtools sort "$ALIGNED_BAM" -o "$SORTED_ALIGNED_BAM"

# Sort the unaligned BAM file
#SORTED_UNALIGNED_BAM="${UNALIGNED_BAM%.bam}_sorted.bam"
samtools sort "$UNALIGNED_BAM" -o "$SORTED_UNALIGNED_BAM"

# Index the sorted aligned BAM file
samtools index "$SORTED_ALIGNED_BAM"
samtools index "$SORTED_UNALIGNED_BAM"

# Generate alignment summary
echo -e "Type\tCount" > "$OUTPUT_FOLDER/lambda_alignment_summary.tsv"
echo -e "Total\t$(samtools view -c $OUTPUT_BAM)" >> "$OUTPUT_FOLDER/lambda_alignment_summary.tsv"
echo -e "Aligned\t$(samtools view -c -F 4 $SORTED_ALIGNED_BAM)" >> "$OUTPUT_FOLDER/lambda_alignment_summary.tsv"
echo -e "Unaligned\t$(samtools view -c -f 4 $SORTED_UNALIGNED_BAM)" >> "$OUTPUT_FOLDER/lambda_alignment_summary.tsv"

# Report completion
echo "Alignment completed:"
echo "Aligned reads: $SORTED_ALIGNED_BAM"
echo "Unaligned reads: $SORTED_UNALIGNED_BAM"
