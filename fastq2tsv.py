# FASTQ to TSV file conversion
import argparse
import gzip
from os import listdir
from os.path import isfile, join

# parse fastq path 

def parse_arguments():
    parser = argparse.ArgumentParser(description='FASTQ to TSV file converter - Takes Raw .FASTQ files and makes one .TSV file ')
    
    # Define string parameters
    parser.add_argument('--fastq_path', type=str, required=True, help='Path to directory containing input FASTQ(.GZ) files (ex. "/data/fastqs/")')
    parser.add_argument('--output_path', type=str, required=True, help='Path to store output files (ex. "/data/outputs/file.tsv")')
    return parser


def read_fastq(fastq_path):  # Ulazi kroz komandnu liniju kao dir
    fastq_name = fastq_path.split('/')[-1]
    gzipped = fastq_name.endswith('.gz')
    openner = gzip.open if gzipped else open

    reads = []

    if not isfile(fastq_path):
        return reads

    with openner(fastq_path, 'rt') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i%4 == 0:
                if line.startswith('@'):
                    id = line.split(' ')[0]
                else:
                    print(f'Error in {fastq_path} line {i} - not an ID line')
                    raise
            elif i%4 == 1:
                seq = line
            elif i%4 == 2:
                if line.startswith('+'):
                    opt = line
                else:
                    print(f'Error in {fastq_path} line {i} - not a + line')
                    raise
            elif i%4 == 3:
                qual = line
                reads.append({
                    'id': id,
                    'seq': seq,
                    'qual': qual,
                })

    return reads


def write_reads(input_path, output_path):
    try:
        fastq_paths = sorted(join(input_path, f) for f in listdir(input_path) if 'merged.fastq' in f and isfile(join(input_path, f)))
    except FileNotFoundError:
        return {
            'fastqs': -1,
            'reads': -1,
        }

    with open(output_path, 'wt') as o:
    
        for fastq_path in fastq_paths:
            reads = read_fastq(fastq_path)
            for read in reads:
                row = f"{read['seq']}\t{read['id']}\n"
                o.write(row)



# Parse the command-line arguments
args = parse_arguments().parse_args()
fastq_path = args.fastq_path
output_path = args.output_path      # called as    dir/fastq.tsv     

#output_path_file = f'{output_path}/fastq.tsv'
print(write_reads(fastq_path, output_path))