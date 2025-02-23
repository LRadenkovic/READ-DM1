import pandas as pd
import numpy as np
import argparse
import seaborn as sns 
import matplotlib.pyplot as plt 

def parse_arguments():
    parser = argparse.ArgumentParser(description='FASTQ to TSV file converter - Takes Raw .FASTQ files and makes one .TSV file ')
    
    parser.add_argument('--lambda_alignment_summary', type=str, help='Path to the summary.tsv file')
    parser.add_argument('--strat_process_output', type=str, help='Path to another TSV file')
    parser.add_argument('--output_dir_path', type=str, help='Path to output folder')
    
    return parser

# Define functions
def get_subset(df, condition):
    return df[condition]

# Main function to encapsulate script logic
def main(lambda_alignment_summpary, strat_process_output, output_dir_path):
    # Import data
    # Import lambda-alignment summary
    print('starting main, line 15, proba.pz')
    df_lambda = pd.read_csv(lambda_alignment_summpary, sep='\t', header=None)
    df_lambda.columns = ['Read_type', 'Count']
    df_lambda.set_index('Read_type', inplace=True)

    lambda_reads = int(df_lambda.loc['Aligned', 'Count'])
    non_lambda_reads = int(df_lambda.loc['Unaligned', 'Count'])

    if lambda_reads == 0 and non_lambda_reads ==0:
        lambda_sum = -1
    else:
        lambda_sum =   lambda_reads + non_lambda_reads  

    # Import df TSV file from STRAT process
    df = pd.read_csv(strat_process_output, sep='\t')
    #print(df.columns)

    directions = list(set(df['direction'])) 
    #print(directions)
    alleles = ['wt', 'exp']
    mix_combos = ['wt_fwd', 'wt_rev', 'exp_fwd', 'exp_rev']

    # wt and exp conditions
    wt_cond = df['len_ins_ext_aln'] < 150  #yaml parser
    exp_cond = df['len_ins_ext_aln'] >= 150 #yaml parser

    fwd_cond = df['direction']== 'fwd'
    rev_cond = df['direction'] == 'rev'

    # all reads counts
    directions_reads_counts = []
    total_reads_count = df.shape[0]
    directions_reads_counts.append(('total',total_reads_count)) 

    index = 0
    for cond in [fwd_cond, rev_cond]:
        subset_df = get_subset(df, cond)
        directions_reads_counts.append((directions[index], subset_df.shape[0]))
        index +=1
    directions_reads_counts    

    # all reads counts, wt plus exp
    alleles_reads_counts = []
    index = 0
    for cond in [wt_cond, exp_cond]:
        subset_df = get_subset(df, cond)
        alleles_reads_counts.append((alleles[index], subset_df.shape[0]))
        index +=1
    alleles_reads_counts 

    # wt and exp counts, fwd and rev
    allele_directions_reads_counts = []

    index = 0
    for cond_allele in [wt_cond, exp_cond]:
        for cond_direction in [fwd_cond, rev_cond]:
            subset_df = get_subset(df, cond_allele & cond_direction )
            allele_directions_reads_counts.append((mix_combos[index], subset_df.shape[0]))
            index+=1
    allele_directions_reads_counts

    # 10th, median, 90th percentile for wt and exp
    percentiles = []
    index = 0
    for cond in [wt_cond, exp_cond]:
        allele = alleles[index]
        for perc in [10, 50, 90]:
            subset_df = get_subset(df, cond)
            p_nth = float(np.percentile(subset_df['len_ins_ext_aln'], perc))
            percentiles.append((f'{allele}_{perc}', p_nth/3))
        index+=1

    # return master_dataframe
    # dodati total reads
    print(strat_process_output)
    sample_name = strat_process_output.split('/')[2]
    master_df = {'Sample':sample_name,
                'Total_reads': lambda_reads+non_lambda_reads,
                'Lambda_reads':  lambda_reads,
                'Non-lambda_reads': non_lambda_reads, 
                'Total on-target reads': directions_reads_counts[0][1], 
                'Fwd reads': directions_reads_counts[1][1], 
                'Rev reads': directions_reads_counts[2][1], 
                'Wt reads':  alleles_reads_counts[0][1], 
                'Exp reads': alleles_reads_counts[1][1], 
                'Wt_fwd': allele_directions_reads_counts[0][1],
                'Wt_rev': allele_directions_reads_counts[1][1],
                'Exp_fwd': allele_directions_reads_counts[2][1],
                'Exp_rev': allele_directions_reads_counts[3][1],
                'Wt_10th': percentiles[0][1],
                'Wt_50th': percentiles[1][1],
                'Wt_90th': percentiles[2][1],
                'Exp_10th': percentiles[3][1],
                'Exp_50th': percentiles[4][1],
                'Exp_90th': percentiles[5][1],
                }
    master_df = pd.DataFrame(master_df, index=[0])      
    save_path_csv = f"{strat_process_output.split('.')[0]}_output.csv"
    master_df.to_csv(save_path_csv, index=False)
    print('===============')
# return master_dataframe
# dodati total reads

    print(strat_process_output.split('/'))

    master_df_percents = {'Sample':sample_name,
                'Total_reads': lambda_reads+non_lambda_reads,
                'Lambda_reads':  f"{lambda_reads} ({round(lambda_reads/(lambda_sum),2)}%)",
                'Non-lambda_reads': f"{non_lambda_reads} ({round(non_lambda_reads/(lambda_sum),2)}%)", 
                'On-target reads': f"{directions_reads_counts[0][1]} (100%) ",
                'Wt reads':  f"{alleles_reads_counts[0][1]} ({round(alleles_reads_counts[0][1]/(directions_reads_counts[0][1]),2)}%)", 
                'Exp reads': f"{alleles_reads_counts[1][1]} ({round(alleles_reads_counts[1][1]/(directions_reads_counts[0][1]),2)}%)",  
                'Wt_fwd': f"{allele_directions_reads_counts[0][1]} ({round(allele_directions_reads_counts[0][1]/(directions_reads_counts[0][1]),2)}%)", 
                'Wt_rev': f"{allele_directions_reads_counts[1][1]} ({round(allele_directions_reads_counts[1][1]/(directions_reads_counts[0][1]),2)}%)", 
                'Exp_fwd': f"{allele_directions_reads_counts[2][1]} ({round(allele_directions_reads_counts[2][1]/(directions_reads_counts[0][1]),2)}%)", 
                'Exp_rev': f"{allele_directions_reads_counts[3][1]} ({round(allele_directions_reads_counts[3][1]/(directions_reads_counts[0][1]),2)}%)", 
                'Wt_10th': percentiles[0][1],
                'Wt_50th': percentiles[1][1],
                'Wt_90th': percentiles[2][1],
                'Exp_10th': percentiles[3][1],
                'Exp_50th': percentiles[4][1],
                'Exp_90th': percentiles[5][1],
                }
    master_df_percents = pd.DataFrame(master_df_percents, index=[0])      
    #master_df_percents.to_csv(f"{strat_process_output.split('.')[0]}_output_percents.csv", index=False)

    save_path_csv = f"{strat_process_output.split('.')[0]}_output_percents.csv"
    print(save_path_csv)
    #master_df.to_csv(f"{strat_process_output.split('.')[0]}_output.csv", index=False
    master_df.to_csv(save_path_csv, index=False)
    


    # Density plot for 'column1'
    #plt.figure(figsize=(13,8))
    index = 0
    for cond in [wt_cond, exp_cond]:
        allele = alleles[index]
        subset_df = get_subset(df, cond)
        plt.figure(figsize=(13,8))
        sns.kdeplot(subset_df['len_ins_ext_aln']/3, shade=True, linewidth=1.5)    
        plt.savefig(f'{output_dir_path}/{sample_name}_{allele}.png', dpi=300)    
        index+=1

   # def plot_waterfall_processed(df, col_len, col_seq, stretch, grid):
   #     width = min(1500, df[col_len].max())
   #     
   #     cond = df['direction'] == 'fwd'
   #     fwd = df[cond][[col_seq, col_len]].sort_values([col_len, col_seq], ascending=[True, True]).reset_index()#

   #     cond = df['direction'] == 'rev'
   #     rev = df[cond][[col_seq, col_len]].sort_values([col_len, col_seq], ascending=[False, False]).reset_index()#

   #     fwd = list(fwd[col_seq])
   #     rev = list(rev[col_seq])
   #     inss = fwd + [width*'I'] + rev#

   #     width = (width)*stretch+stretch
   #     height = len(inss)
   #     image = common.Image.new('RGB', (width, height), 'grey')
   #     draw = common.ImageDraw.Draw(image)
   #     bottom = 0
   #     for i, seq in enumerate(inss):
   #         y = i
   #         # for j, n in enumerate(reversed(seq)):
   #         for j, n in enumerate(seq):
   #             x = j + 1
   #             N = 'CAG'[j%3]
   #             if n == N:
   #                 color = 'black'
   #             else:
   #                 color = common.COLORS[n]
   #             draw.line([(stretch*x, y), (stretch*x, y+1)], width=stretch, fill=color)#

   #     for i in range(width):
   #         y = stretch*i+stretch//2
   #         if i % 3 == 0:
   #             draw.line([(y, 0), (y, height)], width=grid, fill='#AAAAAA')
   #     
   #         if i % 30 == 0:
   #             draw.line([(y, 0), (y, height)], width=grid, fill='white')#

   #         if i % 300 == 0:
   #             draw.line([(y, 0), (y, height)], width=grid, fill='black')#

   #     image.save(output_path)#
#

   # def plot_waterfalls_processed(df, output_path):
   #     col_seq = 'ins_ext_aln'
   #     col_len = 'len_' + col_seq#

   #     cond = df[col_len] <= 150
   #     if len(df[cond]) > 0:
   #         df_sampled = df[cond].sample(n=1000, replace=True, random_state=42)
   #         output_path_50 = f'{output_path}.wtrf.50.png'
   #         plot_waterfall_processed(df_sampled, col_len, col_seq, 15, 2, output_path_50)#

   #     cond = df[col_len] > 150
   #     if len(df[cond]) > 0:
   #         df_sampled = df[cond].sample(n=1000, replace=True, random_state=42)
   #         output_path_51 = f'{output_path}.wtrf.51.png'
   #         plot_waterfall_processed(df_sampled, col_len, col_seq, 7, 1, output_path_51)#
#
#
#


# Parse the command-line arguments
args = parse_arguments().parse_args()
lambda_alignment_summary = args.lambda_alignment_summary
strat_process_output = args.strat_process_output
output_dir_path = args.output_dir_path
    

# Call the main function with the provided input files
main(lambda_alignment_summary, strat_process_output, output_dir_path)

