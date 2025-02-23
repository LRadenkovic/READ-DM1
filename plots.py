from os.path import isfile, join
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import pandas as pd
import seaborn as sns

import common

pd.set_option('display.max_columns', 14)
pd.set_option('display.max_rows', 10)

COLORS = {
    '6': '#7777FF',  # light blue
    '5': '#5555FF',  # blue
    '4': '#FF5555',  # red
    '3': '#FF7777',  # light red
    '0': '#FFFFFF',  # white
    '_': '#333333',  # dark grey
}

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

def rev_comp(seq, comps):
    return ''.join(comps.get(n, n) for n in reversed(seq))

def plot_waterfall_processed(df, col_len, col_seq, stretch, grid, output_path):
    
    width = min(1500, df[col_len].max())
    
    cond = df['direction'] == 'fwd'
    fwd = df[cond][[col_seq, col_len]].sort_values([col_len, col_seq], ascending=[True, True]).reset_index()

    cond = df['direction'] == 'rev'
    rev = df[cond][[col_seq, col_len]].sort_values([col_len, col_seq], ascending=[False, False]).reset_index()

    fwd = list(fwd[col_seq])
    rev = list(rev[col_seq])
    inss = fwd + [width*'I'] + rev

    width = (width)*stretch+stretch
    height = len(inss)
    image = common.Image.new('RGB', (width, height), 'grey')
    draw = common.ImageDraw.Draw(image)
    bottom = 0
    for i, seq in enumerate(inss):
        y = i
        # for j, n in enumerate(reversed(seq)):
        for j, n in enumerate(seq):
            x = j + 1
            N = 'CAG'[j%3]
            if n == N:
                color = 'black'
            else:
                color = common.COLORS[n]
            draw.line([(stretch*x, y), (stretch*x, y+1)], width=stretch, fill=color)

    for i in range(width):
        y = stretch*i+stretch//2
        if i % 3 == 0:
            draw.line([(y, 0), (y, height)], width=grid, fill='#AAAAAA')
    
        if i % 30 == 0:
            draw.line([(y, 0), (y, height)], width=grid, fill='white')

        if i % 300 == 0:
            draw.line([(y, 0), (y, height)], width=grid, fill='black')

    image.save(output_path)


def plot_waterfalls_processed(df, output_path):
    col_seq = 'ins_ext_aln'
    col_len = 'len_' + col_seq

    cond = df[col_len] <= 150
    if len(df[cond]) > 0:
        df_sampled = df[cond].sample(n=1000, replace=True, random_state=42)
        output_path_50 = f'{output_path}.wtrf.50.png'
        plot_waterfall_processed(df_sampled, col_len, col_seq, 15, 2, output_path_50)

    cond = df[col_len] > 150
    if len(df[cond]) > 0:
        df_sampled = df[cond].sample(n=1000, replace=True, random_state=42)
        output_path_51 = f'{output_path}.wtrf.51.png'
        plot_waterfall_processed(df_sampled, col_len, col_seq, 7, 1, output_path_51)

def prepare_row_for_waterfall(row):
    plt_seq = row['seq']
    
    if not pd.isna(row['prefix_flank']):
        # print(row['prefix_flank'])
        plt_seq = plt_seq.replace(row['prefix_flank'], '0000000000')
        plt_seq = plt_seq.replace(row['suffix_flank'], '0000000000')
    
    plt_seq = plt_seq.replace('CTGCTG', '666666')
    plt_seq = plt_seq.replace('CAGCAG', '333333')
    plt_seq = plt_seq.replace('CTG', '555')
    plt_seq = plt_seq.replace('CAG', '444')

    plt_seq = plt_seq.replace('A', '_')
    plt_seq = plt_seq.replace('T', '_')
    plt_seq = plt_seq.replace('C', '_')
    plt_seq = plt_seq.replace('G', '_')

    return plt_seq


def prepare_df_for_waterfall(df, sample_size=4000, max_len_seq=3500):
    df['len_seq'] = df['seq'].str.len()
    cond = df['len_seq'] <= max_len_seq
    plt_df = df[cond].sample(sample_size, replace=True)
    plt_df['plt_seq'] = plt_df.apply(prepare_row_for_waterfall, axis=1)
    #plt_df['fwd'] = plt_df['plt_seq'].str.count('6')
    #plt_df['rev'] = plt_df['plt_seq'].str.count('3')
    
    plt_df['fwd'] = plt_df['plt_seq'].str.count('3')
    plt_df['rev'] = plt_df['plt_seq'].str.count('6')
    
    
    plt_df['dir'] = plt_df['fwd'] - plt_df['rev']
    cond = plt_df['dir'] > 0
    plt_df_fwd = plt_df[cond].sort_values(['len_seq', 'seq'])

    cond = plt_df['dir'] <= 0
    plt_df_rev = plt_df[cond].sort_values(['len_seq', 'seq'], ascending=[False, False])
    plt_df = pd.concat([plt_df_fwd, plt_df_rev])

    return plt_df


def plot_waterfall(df, col_len, col_seq, output_path):
    width = df[col_len].max()
    height = len(df)

    inss = list(df[col_seq])

    image = common.Image.new('RGB', (width, height), 'grey')
    draw = common.ImageDraw.Draw(image)
    bottom = 0
    for i, seq in enumerate(inss):
        y = i
        half = len(seq) / 2
        # left = 0
        left = (width - len(seq)) // 2
        # right = width - len(seq)
        for j, n in enumerate(seq):
            color = COLORS[n]
            draw.point([left+j, i], fill=color)
            # if j < half:
            #     draw.point([left+j, i], fill=color)
            # else:
            #     draw.point([right+j, i], fill=color)

    image.save(output_path)

def plot_histogram(df, x, hue, base, output_histogram):
    fig, ax = plt.subplots(figsize=(16, 10))

    # Create the histogram
    gfg = sns.histplot(df, x=x, discrete=True, hue=hue, multiple='stack')

    # Set major ticks at specified intervals (base)
    loc = plticker.MultipleLocator(base=base)
    ax.xaxis.set_major_locator(loc)  # Apply locator to the axis

    # Rotate the x-tick labels for better readability
    plt.xticks(rotation=90)

    # Save the figure
    fig.savefig(output_histogram, bbox_inches='tight')  # 'tight' ensures proper layout

    # Optionally, close the figure to release memory
    plt.close(fig)

def plot_histogram(df, x, hue, base, output_histogram):
    fig, ax = plt.subplots(figsize=(16, 10))
    gfg = sns.histplot(df, x=x, discrete=True, hue=hue, multiple='stack')
    # gfg.set_xlim(0, 1000)
    # gfg.set_yscale("log")
    loc = plticker.MultipleLocator(base=base)
    gfg.xaxis.set_major_locator(loc)
    gfg.set_xticks(gfg.get_xticks())  # Set the ticks first
    gfg.set_xticklabels(gfg.get_xticklabels(), rotation=90)
    fig.savefig(output_histogram)

def plot_histograms(df, x, output_path):
    cond = df[x] <= 50
    if len(df[cond]) > 0:
        plot_histogram(df[cond].sort_values('direction', ascending=False), x, 'direction', 5, f'{output_path}.hist.50.png')
    cond = df[x] > 50
    if len(df[cond]) > 0:
        plot_histogram(df[cond].sort_values('direction', ascending=False), x, 'direction', 5, f'{output_path}.hist.51.png')


def load_fastq_tsv(path):
    df = pd.read_csv(path, sep='\t', header=None)
    df.columns = ['seq', 'id']
    return df


def load_ontarget(path):
    return common.load_tsv(path, common.COLUMNS_PREPARED)


def load_kmers(path):
    cols = ['direction', 'id', 'prefix_flank', 'ins', 'suffix_flank', 'start_cnt', 'start_stdev', 'end_cnt', 'end_stdev']
    return common.load_tsv(path, cols)


def load_processed(path):
    return common.load_tsv(path)


def load_kmers_processed(path):
    return common.load_tsv(path)     

def parse_arguments():
    
    parser = argparse.ArgumentParser(description='FASTQ to TSV file converter - Takes Raw .FASTQ files and makes one .TSV file ')
    
    # Define string parameters
    #parser.add_argument('--motif', type=str, required=True, help='Motif of the repeated block in forward orientation')
    #parser.add_argument('--limit', type=int, required=True, help='Pathogenic repeat number limit')
    parser.add_argument('--fastq_tsv_path', type=str, required=True, help='Path to directory containing input FASTQ(.GZ) files (ex. "/data/fastqs/")')
    parser.add_argument('--merged_ontarget_path', type=str, required=True, help='Path to directory containing input FASTQ(.GZ) files (ex. "/data/fastqs/")')
    parser.add_argument('--raw_ontarget_merged_path', type=str, required=True, help='Path to directory containing input FASTQ(.GZ) files (ex. "/data/fastqs/")')
    parser.add_argument('--images_path', type=str, required=True, help='Path to directory containing input FASTQ(.GZ) files (ex. "/data/fastqs/")')
    parser.add_argument('--processed_path', type=str, required=True, help='Path to directory containing input FASTQ(.GZ) files (ex. "/data/fastqs/")')

    
    return parser


def load():

    # Parse the command-line arguments
    args = parse_arguments().parse_args()
    #global limit; limit = args.limit 
    #global motif; motif = args.motif
    fastq_tsv_path = args.fastq_tsv_path
    merged_ontarget_path = args.merged_ontarget_path              
    raw_ontarget_merged_path = args.raw_ontarget_merged_path
    processed_path = args.processed_path
    images_path = args.images_path

    path = f'{fastq_tsv_path}/fastq.tsv'   # ovo treba da se napravi i da se zove fastq.tsv
    df_fastq = load_fastq_tsv(path)
    reads = len(df_fastq)
    

    path = f'{merged_ontarget_path}/merged.ontarget.tsv' # merged ontarget tsv
    print(path)
    df_ontarget = load_ontarget(path)
    ontarget = len(df_ontarget)
    cond = df_ontarget['direction'] == 'fwd'
    ontarget_fwd = len(df_ontarget[cond])
    ontarget_rev = ontarget - ontarget_fwd
    
    # Merge FASTQ and on-target dataframes
    df = pd.merge(df_fastq, df_ontarget, how="outer", on=["id", "id"])
    # save merged
    path = f'{raw_ontarget_merged_path}/fastq.ontarget.tsv'
    df.to_csv(path, sep='\t', index=False)
    
    # Waterfall
    plt_df = prepare_df_for_waterfall(df)
    output_path = f'{images_path}/waterfall.png'
    plot_waterfall(plt_df, 'len_seq', 'plt_seq', output_path)

    limit = 150
    path = f'{processed_path}/merged.ontarget.processed.tsv'
    df_processed = load_processed(path)
    cond_fwd = df_processed['direction'] == 'fwd'
    cond_rev = df_processed['direction'] == 'rev'
    df_processed['len'] = df_processed['len_ins_ext_aln'] / 3
    cond = df_processed['len'] <= limit
    processed_50_fwd = sum(cond & cond_fwd)
    processed_50_rev = sum(cond & cond_rev)
    percentiles_50 = df_processed[cond]['len'].quantile([0.1, 0.5, 0.9])
    cond = df_processed['len'] > limit
    processed_51_fwd = sum(cond & cond_fwd)
    processed_51_rev = sum(cond & cond_rev)
    percentiles_51 = df_processed[cond]['len'].quantile([0.1, 0.5, 0.9])

    # Histograms

    output_path = f'{images_path}'
    plot_histograms(df_processed, 'len', output_path)

    #Waterfalls
    output_path = f'{images_path}'
    plot_waterfalls_processed(df_processed, output_path)
    print('completed!\n')
    print('*****************************************************************************')
    
    
load()