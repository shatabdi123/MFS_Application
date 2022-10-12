import pandas as pd
import os, shutil, argparse

#Lots of options, most of which should never need to change
parser = argparse.ArgumentParser()
parser.add_argument('--distance_file',help="structural bed file containing information on the centomere, telomore and Knob information")
parser.add_argument('--gene_info_file',help="File containing information about the gene coordinates")
parser.add_argument('--output_dir',help='Directory to store the outputs containing information on distances of  the genes from the centomere, telomore and Knob')
args = parser.parse_args()

if not args.distance_file:
    print("You need to provide a distance file with the \"--distance_file\" option. Dying now.")

if not args.gene_info_file:
    print("You need to provide a gene coordinate file with the \"--gene_info_file\" option. Dying now.")

if os.path.exists(args.output_dir):
    shutil.rmtree(args.output_dir)

def repeat_split(x):
    x = x.split("=")[1]
    return x

def distance_calculator(distance_file = args.distance_file, gene_file = args.gene_info_file, output_directory= args.output_dir ):
    os.mkdir(output_directory)
    frame_header = ['chrName','start', 'end', 'length', 'type', 'assembly']
    dis_df = pd.read_csv(distance_file, sep='\t', names = frame_header)
    dis_df['type'] = dis_df['type'].apply(repeat_split)
    # Telomere datafrome
    TL_file_df = dis_df[dis_df['type'] == 'TR-1']
    TL_file_df['mid_telo'] = (TL_file_df['start'] + TL_file_df['end'])/2
    # Knob dataframe
    Knob_file_df = dis_df[dis_df['type'] == 'knob180']
    Knob_file_df['mid_knob'] = (Knob_file_df['start'] + Knob_file_df['end']) / 2
    # Centomere dataframe
    Cento_file_df = dis_df[dis_df['type'] == 'cenH3']
    Cento_file_df['mid_cento'] = (Cento_file_df['start'] + Cento_file_df['end']) / 2
    # Gene dataframe
    gene_file_df = pd.read_csv(gene_file,sep='\t')
    gene_file_df['mid_dis'] = (gene_file_df['start'] + gene_file_df['end'])/2
    #merge dataframes
    telo_merged_left = pd.merge(left=TL_file_df, right=gene_file_df, how='left', left_on='chrName', right_on='chrName')
    telo_merged_left["distance_to_telomere"] = abs(telo_merged_left['mid_telo'] - telo_merged_left["mid_dis"])
    telo_merged_left = (telo_merged_left.groupby('ID')['distance_to_telomere','chrName'].min().reset_index())

    knob_merged_left = pd.merge(left=Knob_file_df, right=gene_file_df, how='left', left_on='chrName', right_on='chrName')
    knob_merged_left["distance_to_knob"] = abs(knob_merged_left['mid_knob'] - knob_merged_left["mid_dis"])
    knob_merged_left = (knob_merged_left.groupby('ID')['distance_to_knob', 'chrName'].min().reset_index())

    cento_merged_left = pd.merge(left=Cento_file_df, right=gene_file_df, how='left', left_on='chrName', right_on='chrName')
    cento_merged_left["distance_to_centomere"] = abs(cento_merged_left['mid_cento'] - cento_merged_left["mid_dis"])
    cento_merged_left = (cento_merged_left.groupby('ID')['distance_to_centomere', 'chrName'].min().reset_index())

    final_merge = pd.merge(left=telo_merged_left, right=knob_merged_left, how='left', left_on='ID', right_on='ID')
    final_merge = pd.merge(left=final_merge, right=cento_merged_left, how='left', left_on='ID',
                           right_on='ID')
    final_merge = final_merge[['ID', 'chrName', 'distance_to_telomere', 'distance_to_knob', 'distance_to_centomere']]
    final_merge.to_csv(output_directory + '/gene_structural_distance' + '.txt', sep='\t', index=False, header=True)
distance_calculator()