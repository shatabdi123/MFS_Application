# importing csv module
import csv
import pandas as pd
import re
import os, shutil, argparse

#Lots of options, most of which should never need to change
parser = argparse.ArgumentParser()
parser.add_argument('--gff_file',help="GFF format file containing information on the genome of the organism being studied")
parser.add_argument('--meta_info_dir',help="Directory to storing the meta informations generated on the fly while processing the GFF files")
parser.add_argument('--output_dir',help='Directory to store the outputs from GFF files')
args = parser.parse_args()

if not args.gff_file:
    print("You need to provide a gff file with the \"--gff_file\" option. Dying now.")

if os.path.exists(args.meta_info_dir):
    shutil.rmtree(args.meta_info_dir)

if os.path.exists(args.output_dir):
    shutil.rmtree(args.output_dir)

def ID_conversion(x):
    x = x.split("_")[0]
    return x

def gff_file_calculator(gff_file = args.gff_file,split_directory = args.meta_info_dir,output_directory= args.output_dir):
# def gff_file_calculator(gff_file= "C:/Users/Shatabdi/OneDrive/NAM_Datasets/B97/B97.gff3", split_directory="C:/Users/Shatabdi/PycharmProjects/NAM_FeatureStore/GFF_Meta_Files/",
#                             output_directory="C:/Users/Shatabdi/PycharmProjects/NAM_FeatureStore/GFF_Files_Output"):
    os.mkdir(split_directory)
    os.mkdir(output_directory)
    frame_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'ID']
    rows = []
    # reading csv file
    with open(gff_file, 'r') as csvfile:
        # creating a csv reader object
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            if not row[0].startswith("#"):
                rows.append(row)
    frame = pd.DataFrame(rows, columns=frame_header, index=None)
    datas = {}
    for i, g in frame.groupby('type'):
        datas.update({'data_' + str(i): g.reset_index(drop=True)})
        datas['data_' + str(i)].to_csv(
            split_directory + '/' + str(
                i) + '.txt', sep='\t', index=False, header=False)
    for filename in os.listdir(split_directory):
        rows = []
        with open(split_directory + '/' + filename, "r") as csvfile:
            csvreader = csv.reader(csvfile, delimiter='\t')
            for row in csvreader:
                mod = re.split(";* *\w*=", row[8])[1:]
                row = row[0:8]
                row.extend(mod)
                rows.append(row)
        if filename.startswith("chromosome"):
            chrom_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'ID', 'Name']
            chrom_df = pd.DataFrame(rows, columns=chrom_header, index=None)
        elif filename.startswith("exon"):
            exon_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'Parent', 'Name',
                           'ensembl_end_phase', 'ensembl_phase', 'exon_id', 'rank']
            exon_df = pd.DataFrame(rows, columns=exon_header, index=None)
            # Exon num per transcript
            data = exon_df.groupby('Parent').apply(lambda x: sum(list(abs(pd.to_numeric(x.end) - pd.to_numeric(x.start) + 1)))).reset_index(name='exonLength')
            data['count'] = exon_df.groupby('Parent')['Parent'].count().reset_index(name="count")['count']
            exon_df = pd.DataFrame(columns=["parent","exonLength","exonNum"])
            exon_df['parent'] = data['Parent']
            exon_df['exonLength'] = data['exonLength']
            exon_df['exonNum'] = data['count']
            exon_df['ID'] = exon_df['parent'].apply(ID_conversion)
            exon_df.drop(['parent'], axis=1)
            exon_df["exonLength"] = exon_df.groupby("ID")["exonLength"].sum().reset_index(name="exonLength")["exonLength"]
            exon_df["exonNum"] = exon_df.groupby("ID")["exonNum"].sum().reset_index(name="exonNum")["exonNum"]
            exon_df['avgExonLength'] = round(exon_df["exonLength"]/exon_df["exonNum"],2)
            exon_df.to_csv(output_directory + '/exonNum' + '.txt' , sep='\t', index=False, header=True)
        elif filename.startswith("five"):
            five_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'Parent']
            five_df = pd.DataFrame(rows, columns=five_header, index=None)
            five_df = five_df.groupby('Parent').apply(lambda x: sum(list(abs(pd.to_numeric(x.end) - pd.to_numeric(x.start) + 1)))).reset_index(name='fiveUTRLen')
            five_df['ID'] = five_df['Parent'].apply(ID_conversion)
            five_df.drop(['Parent'], axis=1)
            five_df = five_df.groupby('ID')['fiveUTRLen'].sum().reset_index(name="fiveUTRLen")
            five_df.to_csv(output_directory + '/fiveUTRlen' + '.txt' , sep='\t', index=False, header=True)
        elif filename.startswith("three"):
            three_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'Parent']
            three_df = pd.DataFrame(rows, columns=three_header, index=None)
            three_df = three_df.groupby('Parent').apply(
                lambda x: sum(list(abs(pd.to_numeric(x.end) - pd.to_numeric(x.start) + 1)))).reset_index(
                name='threeUTRlen')
            three_df['ID'] = three_df['Parent'].apply(ID_conversion)
            three_df.drop(['Parent'], axis=1)
            three_df = three_df.groupby('ID')['threeUTRlen'].sum().reset_index(name="threeUTRlen")
            three_df.to_csv(output_directory + '/threeUTRlen' + '.txt', sep='\t', index=False, header=True)
        elif filename.startswith("gene"):
            gene_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'ID', 'biotype',
                           'logic_name']
            gene_df = pd.DataFrame(rows, columns=gene_header, index=None)
            # gene_df chr,id,start,end,strand
            gene_df = gene_df[['chrName', 'ID', 'start', 'end', 'strand']]
            # generate gene length output file
            gene_df['geneLength']  = pd.to_numeric(gene_df['end']) - pd.to_numeric(gene_df['start']) + 1
            gene_df.to_csv(output_directory + '/geneLength' + '.txt' , sep='\t', index=False, header=True)
            # calculate chromosomal distance
            gene_dis = gene_df[['chrName', 'ID', 'start', 'end']]
            temp_2 = gene_dis.groupby("chrName")['end','start'].apply(lambda x:  (abs((max(pd.to_numeric(x.end))-min(pd.to_numeric(x.start)))/2))).reset_index()
            merged_left = pd.merge(left=gene_dis, right=temp_2, how='left', left_on='chrName', right_on='chrName')
            merged_left['distance'] = merged_left[0]- pd.to_numeric(merged_left['start'])
            merged_left = merged_left[["ID","chrName","distance"]]
            merged_left.to_csv(output_directory + '/genedis' + '.txt' , sep='\t', index=False, header=True)

        elif filename.startswith("mRNA"):
            mRNA_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'ID', 'Parent',
                           'biotype', 'transcript_id', 'canonical_transcript']
            mRNA_df = pd.DataFrame(rows, columns=mRNA_header, index=None)
            mRNA_isoform = mRNA_df.groupby("Parent")["type"].count().reset_index(name="isoformCount")
            mRNA_isoform.to_csv(output_directory + '/isoform' + '.txt' , sep='\t', index=False, header=True)
            mRNA_df = mRNA_df.loc[mRNA_df["canonical_transcript"] == str(1), ["Parent","start","end"]]
            mRNA_df["mrnaLength"] = pd.to_numeric(mRNA_df['end']) - pd.to_numeric(mRNA_df['start']) + 1
            mRNA_df = mRNA_df[["Parent","mrnaLength"]]
            mRNA_df.to_csv(output_directory + '/canonicalMrnaLen' + '.txt' , sep='\t', index=False, header=True)
        elif filename.startswith("CDS"):
            CDS_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'ID', 'Parent',
                          'protein_id']
            CDS_df = pd.DataFrame(rows, columns=CDS_header, index=None)
            CDS_df = CDS_df.groupby('Parent').apply(lambda x: sum(list(abs(pd.to_numeric(x.end) - pd.to_numeric(x.start) + 1)))).reset_index(name='CDSLen')
            CDS_df['ID'] = CDS_df['Parent'].apply(ID_conversion)
            CDS_df.drop(['Parent'], axis=1)
            CDS_df = CDS_df.groupby('ID')['CDSLen'].sum().reset_index(name="CDSLen")
            CDS_df.to_csv(output_directory + '/CDSlen' + '.txt' , sep='\t', index=False, header=True)
        elif filename.startswith("scaffold"):
            scaffold_header = ['chrName', 'source', 'type', 'start', 'end', 'dot', 'strand', 'value', 'ID', 'Name']
            scaffold_df = pd.DataFrame(rows, columns=scaffold_header, index=None)
    return output_directory

gff_file_calculator()