import pandas as pd
import numpy as np
import json
from pandas import DataFrame
from bson.json_util import dumps
from src.common.database import Database
from src.models.charts import charts
from sklearn.utils import resample
from itertools import islice

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

def Convert(lst):
    res_dct = {lst[i]: 1 for i in range(0, len(lst), 1)}
    return res_dct

col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/colfile_new.txt", "r")
col_name = []
for line in col_file:
    col_name.append(line.strip())

class plot_seq():

    @staticmethod
    def structure_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('struc_dis', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())

        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)

        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

        for column in df.columns:
            if column == 'distance':
                df[column] = df[column].abs()
                # print(df[df['ID']=="Zm00001eb207600"])
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df)

        return plot_sns

    @staticmethod
    def downsampled_structure_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('struc_dis', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        print(df_downsampled.head())
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_plots'):
            plot_sns = charts.downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_scatter_plots'):
            plot_sns = charts.downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'downsampled_dendogram_plots'):
            plot_sns = charts.downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'downsampled_dendogram_with_cluster_plots'):
            plot_sns = charts.downsampled_cluster_plot(df_downsampled)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns


    @staticmethod
    def codon_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        # print(df[df['ID'] == 'Zm00001eb125200'])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'ENC_plot'):
            plot_sns = charts.ENC_plot(df)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df)
        return plot_sns

    @staticmethod
    def downsampled_codon_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())

        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())

        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'down_hierarchial_plots'):
            plot_sns = charts.downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'down_hierarchial_scatter_plots'):
            plot_sns = charts.downsampled_hierarchial_scatter_plot(df_downsampled, cluster)
        elif (analysis[0] == 'down_dendogram_plots'):
            plot_sns = charts.downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'ENC_plot'):
            plot_sns = charts.ENC_plot(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def ProteinStructure_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')
        if 'SignalP' in df.columns:
            df['SignalP'] = df['SignalP'].replace(1, 'SignalP')
            df['SignalP'] = df['SignalP'].replace(0, 'No SignalP')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        # print(df[df['ID'] == 'Zm00001eb125200'])
        if (analysis[0] == 'categorical_bar_chart'):
            plot_sns = charts.categorical_bar_chart(df)
        elif (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df)
        return plot_sns

    @staticmethod
    def downsampled_ProteinStructure_plot(select, analysis,cluster):
        cluster = cluster[0]
        print("cluster is: ", cluster)
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if 'SignalP' in df.columns:
            df['SignalP'] = df['SignalP'].replace(1, 'SignalP')
            df['SignalP'] = df['SignalP'].replace(0, 'No SignalP')
        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if (analysis[0] == 'categorical_bar_chart'):
            plot_sns = charts.categorical_bar_chart(df_downsampled)
        elif (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_plots'):
            plot_sns = charts.Count_downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_scatter_plots'):
            plot_sns = charts.Count_downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'Count_down_dendogram_plots'):
            plot_sns = charts.Count_downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def ProteinLocalization_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())

        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())

        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        # print(df[df['ID'] == 'Zm00001eb125200'])
        if (analysis[0] == 'categorical_bar_chart'):
            plot_sns = charts.categorical_bar_chart(df)
        elif (analysis[0] == 'Kmode_cluster_plots'):
            plot_sns = charts.Kmode_cluster_plot(df)
        return plot_sns

    @staticmethod
    def downsampled_ProteinLocalization_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if (analysis[0] == 'categorical_bar_chart'):
            plot_sns = charts.categorical_bar_chart(df_downsampled)
        elif (analysis[0] == 'Kmode_cluster_plots'):
            plot_sns = charts.Kmode_cluster_plot(df_downsampled)
        return plot_sns

    @staticmethod
    def GeneExp_plot(select, analysis):
        dict = {}
        print(select)

        print(select[0])
        if select[0] == 'Varotto_Lab':
            Vartobo_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Vartobo.txt", "r")
            Vartobo_name = []
            for line in Vartobo_file:
                Vartobo_name.append(line.strip())
            dict.update(Convert(Vartobo_name))
        elif select[0] == 'Scanlon_Lab':
            Schanlon_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Schanlon.txt",
                                 "r")
            Schanlon_name = []
            for line in Schanlon_file:
                Schanlon_name.append(line.strip())
            dict.update(Convert(Schanlon_name))
        elif select[0] == 'Pereira_Lab':
            Pereira_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Pereira.txt", "r")
            Pereira_name = []
            for line in Pereira_file:
                Pereira_name.append(line.strip())
            dict.update(Convert(Pereira_name))
        elif select[0] == 'Kaeppler_and_Walley':
            Kaepler_Walley_file = open(
                "C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Kaepler_Walley.txt", "r")
            Kaepler_Walley_name = []
            for line in Kaepler_Walley_file:
                Kaepler_Walley_name.append(line.strip())
            dict.update(Convert(Kaepler_Walley_name))
        elif select[0] == 'Springer_Lab':
            Springer_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Springer.txt",
                                 "r")
            Springer_name = []
            for line in Springer_file:
                Springer_name.append(line.strip())
            dict.update(Convert(Springer_name))
        elif select[0] == 'Hochholdinger_Lab':
            Holchholdiner_file = open(
                "C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Holchholdiner.txt", "r")
            Holchholdiner_name = []
            for line in Holchholdiner_file:
                Holchholdiner_name.append(line.strip())
            dict.update(Convert(Holchholdiner_name))
        elif select[0] == 'Fowler_Lab':
            fowler_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Fowler.txt", "r")
            fowler_name = []
            for line in fowler_file:
                fowler_name.append(line.strip())
            dict.update(Convert(fowler_name))
        else:
            dict.update(Convert(select))
        # id_dict = {}
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        new_dict = {}
        print(dict)
        if (analysis[0] == 'Heatmap_plots') or (analysis[0] == 'PCA_2D_samples_plots') \
                or (analysis[0] == 'PCA_3D_samples_plots') or (analysis[0] == 'PCA_2D_label_plots') \
                or (analysis[0] == 'Count_downsampled_Hierarchial_heatmap')\
                or (analysis[0] == 'Count_downsampled_dendogram_plot') :
            # id_dict["_id"] = 0
            cursor = Database.find_by_structure('gene_exp', dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)
            print(df.head())
        else:

            new_dict.update(take(5, dict.items()))
            new_dict['ID'] = 1
            new_dict['_id'] = 0
            print("new dict")
            print(new_dict)
            cursor = Database.find_by_structure('gene_exp', new_dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            if 'classical_label' in dict:
                del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            if 'core_label' in dict:
                del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            if 'Origin' in dict:
                del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        # print(df[df['ID'] == 'Zm00001eb125200'])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Heatmap_plots'):
            plot_sns = charts.Heatmap_plots(df)
        elif (analysis[0] == 'PCA_2D_samples_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df.columns)):
                df = df.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'PCA_3D_samples_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df.columns)):
                df = df.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.Gene_PCA_3D_samples_plots(df)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df.columns)):
                df = df.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.PCA_2D_label_plots(df)
        return plot_sns

    @staticmethod
    def downsampled_GeneExp_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = {}
        print(select)

        print(select[0])
        if select[0] == 'Varotto_Lab':
            Vartobo_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Vartobo.txt", "r")
            Vartobo_name = []
            for line in Vartobo_file:
                Vartobo_name.append(line.strip())
            dict.update(Convert(Vartobo_name))
        elif select[0] == 'Scanlon_Lab':
            Schanlon_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Schanlon.txt",
                                 "r")
            Schanlon_name = []
            for line in Schanlon_file:
                Schanlon_name.append(line.strip())
            dict.update(Convert(Schanlon_name))
        elif select[0] == 'Pereira_Lab':
            Pereira_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Pereira.txt", "r")
            Pereira_name = []
            for line in Pereira_file:
                Pereira_name.append(line.strip())
            dict.update(Convert(Pereira_name))
        elif select[0] == 'Kaeppler_and_Walley':
            Kaepler_Walley_file = open(
                "C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Kaepler_Walley.txt", "r")
            Kaepler_Walley_name = []
            for line in Kaepler_Walley_file:
                Kaepler_Walley_name.append(line.strip())
            dict.update(Convert(Kaepler_Walley_name))
        elif select[0] == 'Springer_Lab':
            Springer_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Springer.txt",
                                 "r")
            Springer_name = []
            for line in Springer_file:
                Springer_name.append(line.strip())
            dict.update(Convert(Springer_name))
        elif select[0] == 'Hochholdinger_Lab':
            Holchholdiner_file = open(
                "C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Holchholdiner.txt", "r")
            Holchholdiner_name = []
            for line in Holchholdiner_file:
                Holchholdiner_name.append(line.strip())
            dict.update(Convert(Holchholdiner_name))
        elif select[0] == 'Fowler_Lab':
            fowler_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Fowler.txt", "r")
            fowler_name = []
            for line in fowler_file:
                fowler_name.append(line.strip())
            dict.update(Convert(fowler_name))
        else:
            dict.update(Convert(select))
        # id_dict = {}
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        new_dict = {}
        print(dict)
        if (analysis[0] == 'Heatmap_plots') or (analysis[0] == 'PCA_2D_samples_plots') \
                or (analysis[0] == 'PCA_3D_samples_plots') or (analysis[0] == 'PCA_2D_label_plots') \
                or (analysis[0] == 'PCA_2D_biplot_plots') or (analysis[0] == 'Count_downsampled_Hierarchial_heatmap') \
                or (analysis[0] == 'PCA_3D_label_plots') or (analysis[0] == 'Count_downsampled_dendogram_plot'):
            # id_dict["_id"] = 0
            cursor = Database.find_by_structure('gene_exp', dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)
            print(df.head())
        else:

            new_dict.update(take(5, dict.items()))
            new_dict['ID'] = 1
            new_dict['_id'] = 0
            print("new dict")
            print(new_dict)
            cursor = Database.find_by_structure('gene_exp', new_dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            if 'classical_label' in dict:
                del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            if 'core_label' in dict:
                del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            if 'Origin' in dict:
                del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())

        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'Count_downsampled_Hierarchial_heatmap'):
            plot_sns = charts.Count_downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'Count_downsampled_hierarchial_scatter_plot'):
            plot_sns = charts.Count_downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'Count_downsampled_dendogram_plot'):
            plot_sns = charts.Count_downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Heatmap_plots'):
            plot_sns = charts.Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_samples_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_samples_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.Gene_PCA_3D_samples_plots(df_downsampled)
        # elif (analysis[0] == 'PCA_2D_label_plots'):
        #     if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df_downsampled.columns)):
        #         df_downsampled = df_downsampled.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
        #     plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_biplot_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.PCA_3D_biplot_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            if set(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_gene_abundance', 'gene_breadth', 'Max_FPKM'], axis=1)
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def ProteinExp_plot(select, analysis):
        dict = Convert(select)
        id_dict ={}
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        if (analysis[0] == 'Heatmap_plots') or (analysis[0] == 'PCA_2D_samples_plots') \
                or (analysis[0] == 'PCA_3D_samples_plots') or (analysis[0] == 'PCA_2D_label_plots'):
            id_dict["_id"] = 0
            cursor = Database.find_all('prot_exp',id_dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)
        else:

            cursor = Database.find_by_structure('prot_exp', dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Heatmap_plots'):
            plot_sns = charts.Heatmap_plots(df)
        elif (analysis[0] == 'PCA_2D_samples_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df.columns)):
                df = df.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.Protein_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'PCA_3D_samples_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df.columns)):
                df = df.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.Protein_PCA_3D_samples_plots(df)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df.columns)):
                df = df.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.PCA_2D_label_plots(df)
        return plot_sns

    @staticmethod
    def downsampled_ProteinExp_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        id_dict = {}
        if (analysis[0] == 'Heatmap_plots') or (analysis[0] == 'PCA_2D_samples_plots') \
                or (analysis[0] == 'PCA_3D_samples_plots')or (analysis[0] == 'PCA_2D_label_plots'):
            id_dict["_id"] = 0
            cursor = Database.find_all('prot_exp',id_dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)
        else:

            cursor = Database.find_by_structure('prot_exp', dict)
            list_cur = list(cursor)
            df = pd.DataFrame(list_cur)

        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'Count_downsampled_Hierarchial_heatmap'):
            plot_sns = charts.Count_downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'Count_downsampled_hierarchial_scatter_plot'):
            plot_sns = charts.Count_downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'Count_downsampled_dendogram_plot'):
            plot_sns = charts.Count_downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Heatmap_plots'):
            plot_sns = charts.Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_samples_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.Protein_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_samples_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.Protein_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_biplot_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.PCA_3D_biplot_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            if set(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF']).issubset(set(df_downsampled.columns)):
                df_downsampled = df_downsampled.drop(['avg_protein_abundance', 'protein_breadth', 'Max_dSNAF'], axis=1)
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def seq_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        label_dict = {}
        print(dict)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] =1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] =1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        if 'dinucleic' in str(select) and "count_dinucleic_box_plots" in str(analysis):
            del dict['dinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_2"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('di_nuc', dict)
        elif 'dinucleic' in str(select) and "frequency_dinucleic_box_plots" in str(analysis):
            del dict['dinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_2"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('norm_di_nuc', dict)
        elif 'trinucleic' in str(select) and "count_tri_nucleic_box_plots" in str(analysis):
            del dict['trinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_3"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('tri_nuc', dict)
        elif 'trinucleic' in str(select) and "frequency_tri_nucleic_box_plots" in str(analysis):
            del dict['trinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_3"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('norm_tri_nuc', dict)
        list_cur = list(cursor)
        print(list_cur)

        df = pd.DataFrame(list_cur)
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df.head())
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        for name in df.columns:
            if name.startswith("kmer_2") or name.startswith("kmer_3"):
                new_name = name.split("_")[2]
                df.rename(columns={name: new_name}, inplace=True)
        if "classical_label" in df.columns:
            df = df.melt(id_vars=["classical_label"],
                         var_name="variable",
                         value_name="value")
        elif "core_label" in df.columns:
            df = df.melt(id_vars=["core_label"],
                         var_name="variable",
                         value_name="value")
        elif "Origin" in df.columns:
            df = df.melt(id_vars=["Origin"],
                         var_name="variable",
                         value_name="value")
        else:
            df['None'] = 1
            df = df.melt(id_vars=["None"],
                         var_name="variable",
                         value_name="value")

        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'count_dinucleic_box_plots'):
            plot_sns = charts.dinucleic_box_plots(df)
        elif (analysis[0] == 'count_tri_nucleic_box_plots'):
            plot_sns = charts.trinucleic_box_plots(df)
        elif (analysis[0] == 'frequency_dinucleic_box_plots'):
            plot_sns = charts.frequency_dinucleic_box_plots(df)
        elif (analysis[0] == 'frequency_tri_nucleic_box_plots'):
            plot_sns = charts.frequency_trinucleic_box_plots(df)


        return plot_sns

    @staticmethod
    def downsampled_seq_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        label_dict = {}
        print(dict)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        if 'dinucleic' in str(select) and "count_dinucleic_box_plots" in str(analysis):
            del dict['dinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_2"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('di_nuc', dict)
        elif 'dinucleic' in str(select) and "frequency_dinucleic_box_plots" in str(analysis):
            del dict['dinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_2"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('norm_di_nuc', dict)
        elif 'trinucleic' in str(select) and "count_tri_nucleic_box_plots" in str(analysis):
            del dict['trinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_3"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('tri_nuc', dict)
        elif 'trinucleic' in str(select) and "frequency_tri_nucleic_box_plots" in str(analysis):
            del dict['trinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_3"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('norm_tri_nuc', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=label_df, right=df, how='left', left_on='ID', right_on='ID')
        print(df.head())
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        if "ID" in df_downsampled.columns:
            df_downsampled = df_downsampled.drop("ID", axis=1)
        for name in df_downsampled.columns:
            if name.startswith("kmer_2") or name.startswith("kmer_3"):
                new_name = name.split("_")[2]
                df_downsampled.rename(columns={name: new_name}, inplace=True)
        if "classical_label" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["classical_label"],
                                                 var_name="variable",
                                                 value_name="value")
        elif "core_label" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["core_label"],
                                                 var_name="variable",
                                                 value_name="value")
        elif "Origin" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["Origin"],
                                                 var_name="variable",
                                                 value_name="value")
        print(df_downsampled.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'count_dinucleic_box_plots'):
            plot_sns = charts.dinucleic_box_plots(df_downsampled)
        elif (analysis[0] == 'count_tri_nucleic_box_plots'):
            plot_sns = charts.trinucleic_box_plots(df_downsampled)
        elif (analysis[0] == 'frequency_dinucleic_box_plots'):
            plot_sns = charts.frequency_dinucleic_box_plots(df_downsampled)
        elif (analysis[0] == 'frequency_tri_nucleic_box_plots'):
            plot_sns = charts.frequency_trinucleic_box_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def same_len_downsampled_seq_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        label_dict = {}
        print(dict)
        cursor_len = Database.find_by_structure('struc_dis', {"ID":1,"CDSlength":1,"_id":0})
        list_cursor_len = list(cursor_len )
        len_df = DataFrame(list_cursor_len)
        len_df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        if 'dinucleic' in str(select) and "dinucleic_box_plots" in str(analysis):
            del dict['dinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_2"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('di_nuc', dict)
        elif 'trinucleic' in str(select) and "tri_nucleic_box_plots" in str(analysis):
            del dict['trinucleic']
            cols = []
            for name in col_name:
                if name.startswith("kmer_3"):
                    cols.append(name)
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('tri_nuc', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        df.drop_duplicates(subset="ID", keep=False, inplace=True)

        df = pd.merge(left=label_df, right=df, how='left', left_on='ID', right_on='ID')

        df.to_csv('C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/k_mer_len5' + '.txt', sep='\t',
                  index=False, header=True)
        df = pd.merge(left=df, right=len_df, how='left', left_on='ID', right_on='ID')
        print(df.head())
        df.to_csv('C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/k_mer_len6' + '.txt', sep='\t',
              index=False, header=True)
        df = df[(df['CDSlength'] >= 2000) & (df['CDSlength'] <= 3000)]
        df.to_csv('C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/k_mer_len7' + '.txt', sep='\t',
                  index=False, header=True)
        # df = df.loc[:, df.columns != 'genelength']

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            df_downsampled.to_csv('C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/k_mer_len3' + '.txt', sep='\t',
                      index=False, header=True)
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        if "ID" in df_downsampled.columns:
            df_downsampled = df_downsampled.drop("ID", axis=1)
        for name in df_downsampled.columns:
            if name.startswith("kmer_2") or name.startswith("kmer_3"):
                new_name = name.split("_")[2]
                df_downsampled.rename(columns={name: new_name}, inplace=True)
        if "classical_label" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["classical_label"],
                                                 var_name="variable",
                                                 value_name="value")
        elif "core_label" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["core_label"],
                                                 var_name="variable",
                                                 value_name="value")
        elif "Origin" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["Origin"],
                                                 var_name="variable",
                                                 value_name="value")
        print(df_downsampled.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'dinucleic_box_plots'):
            plot_sns = charts.dinucleic_box_plots(df_downsampled)
        elif (analysis[0] == 'tri_nucleic_box_plots'):
            plot_sns = charts.trinucleic_box_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def pep_seq_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        label_dict = {}
        print(dict)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        if 'DC' in str(select) and "count_dipeptide_box_plots" in str(analysis):
            del dict['DC']
            cols = []
            for name in col_name:
                if name.startswith("DC"):
                    cols.append(name)
            cols = cols[0:20]
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('di_pep', dict)
        elif 'TC' in str(select) and "count_tripeptide_box_plots" in str(analysis):
            del dict['TC']
            cols = []
            for name in col_name:
                if name.startswith("TC"):
                    cols.append(name)
            cols = cols[0:20]
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('tri_pep', dict)
        list_cur = list(cursor)
        print(list_cur)

        df = pd.DataFrame(list_cur)
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df.head())
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        for name in df.columns:
            if name.startswith("DC") or name.startswith("TC"):
                new_name = name.split("_")[1]
                df.rename(columns={name: new_name}, inplace=True)
        if "classical_label" in df.columns:
            df = df.melt(id_vars=["classical_label"],
                         var_name="variable",
                         value_name="value")
        elif "core_label" in df.columns:
            df = df.melt(id_vars=["core_label"],
                         var_name="variable",
                         value_name="value")
        elif "Origin" in df.columns:
            df = df.melt(id_vars=["Origin"],
                         var_name="variable",
                         value_name="value")
        else:
            df["None"] = 1
            df = df.melt(id_vars=["None"],
                         var_name="variable",
                         value_name="value")
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'count_dipeptide_box_plots'):
            plot_sns = charts.dipeptide_box_plots(df)
        elif (analysis[0] == 'count_tripeptide_box_plots'):
            plot_sns = charts.tripeptide_box_plots(df)

        return plot_sns

    @staticmethod
    def pep_downsampled_seq_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        label_dict = {}
        print(dict)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        if 'DC' in str(select) and "count_dipeptide_box_plots" in str(analysis):
            del dict['DC']
            cols = []
            for name in col_name:
                if name.startswith("DC"):
                    cols.append(name)
            cols = cols[0:20]
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('di_pep', dict)
        elif 'TC' in str(select) and "count_tripeptide_box_plots" in str(analysis):
            del dict['TC']
            cols = []
            for name in col_name:
                if name.startswith("TC"):
                    cols.append(name)
            cols = cols[0:20]
            dict1 = Convert(cols)
            dict.update(dict1)
            dict['_id'] = 0
            print(dict)
            cursor = Database.find_by_structure('tri_pep', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=label_df, right=df, how='left', left_on='ID', right_on='ID')
        print(df.head())
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        if "ID" in df_downsampled.columns:
            df_downsampled = df_downsampled.drop("ID", axis=1)
        for name in df_downsampled.columns:

            if name.startswith("DC") or name.startswith("TC"):
                new_name = name.split("_")[1]
                df_downsampled.rename(columns={name: new_name}, inplace=True)
        if "classical_label" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["classical_label"],
                                                 var_name="variable",
                                                 value_name="value")
        elif "core_label" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["core_label"],
                                                 var_name="variable",
                                                 value_name="value")
        elif "Origin" in df_downsampled.columns:
            df_downsampled = df_downsampled.melt(id_vars=["Origin"],
                                                 var_name="variable",
                                                 value_name="value")
        print(df_downsampled.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'count_dipeptide_box_plots'):
            plot_sns = charts.dipeptide_box_plots(df_downsampled)
        elif (analysis[0] == 'count_tripeptide_box_plots'):
            plot_sns = charts.tripeptide_box_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def ChromStates_plot(select, analysis):
        dict = Convert(select)
        id_dict = {}
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())

        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        # print(df[df['ID'] == 'Zm00001eb125200'])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df)
        elif (analysis[0] == 'Leaf_chromatin_state_plot'):
            plot_sns = "Leaf_emissions_11.png"
        elif (analysis[0] == 'ear_chromatin_state_plot'):
            plot_sns = 'ear_emissions_9.png'
        elif (analysis[0] == 'Leaf_fold_enrichment_plot'):
            plot_sns = 'Leaf_11_overlap.png'
        elif (analysis[0] == 'ear_fold_enrichment_plot'):
            plot_sns = 'ear_9_overlap.png'
        elif (analysis[0] == 'Leaf_TES_fold_enrichment_plot'):
            plot_sns = 'Leaf_11_RefSeqTES_neighborhood.png'
        elif (analysis[0] == 'Leaf_TSS_fold_enrichment_plot'):
            plot_sns = 'Leaf_11_RefSeqTSS_neighborhood.png'
        elif (analysis[0] == 'ear_TES_fold_enrichment_plot'):
            plot_sns = 'ear_9_RefSeqTES_neighborhood.png'
        elif (analysis[0] == 'ear_TSS_fold_enrichment_plot'):
            plot_sns = 'ear_9_RefSeqTSS_neighborhood.png'
        return plot_sns

    @staticmethod
    def downsampled_ChromStates_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_plots'):
            plot_sns = charts.Count_downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_scatter_plots'):
            plot_sns = charts.Count_downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'Count_down_dendogram_plots'):
            plot_sns = charts.Count_downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        elif (analysis[0] == 'Leaf_chromatin_state_plot'):
            plot_sns = "Leaf_emissions_11.png"
        elif (analysis[0] == 'ear_chromatin_state_plot'):
            plot_sns = 'ear_emissions_9.png'
        elif (analysis[0] == 'Leaf_fold_enrichment_plot'):
            plot_sns = 'Leaf_11_overlap.png'
        elif (analysis[0] == 'ear_fold_enrichment_plot'):
            plot_sns = 'ear_9_overlap.png'
        elif (analysis[0] == 'Leaf_TES_fold_enrichment_plot'):
            plot_sns = 'Leaf_11_RefSeqTES_neighborhood.png'
        elif (analysis[0] == 'Leaf_TSS_fold_enrichment_plot'):
            plot_sns = 'Leaf_11_RefSeqTSS_neighborhood.png'
        elif (analysis[0] == 'ear_TES_fold_enrichment_plot'):
            plot_sns = 'ear_9_RefSeqTES_neighborhood.png'
        elif (analysis[0] == 'ear_TSS_fold_enrichment_plot'):
            plot_sns = 'ear_9_RefSeqTSS_neighborhood.png'
        return plot_sns

    @staticmethod
    def ATACseq_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            if column == 'distance':
                df[column] = df[column].abs()
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df)

        return plot_sns

    @staticmethod
    def downsampled_ATACseq_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        print(df_downsampled.head())
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_plots'):
            plot_sns = charts.downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_scatter_plots'):
            plot_sns = charts.downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'downsampled_dendogram_plots'):
            plot_sns = charts.downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def miRNA_Count_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        for column in df.columns:
            if column == 'ID':
                pass
            else:
                df[column] = df[column].replace(0, 'No')
                df[column] = df[column].replace(1, 'Yes')
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')


        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        if (analysis[0] == 'categorical_bar_chart'):
            plot_sns = charts.categorical_bar_chart(df)
        return plot_sns

    @staticmethod
    def downsampled_miRNA_Count_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        for column in df.columns:
            if column == 'ID':
                pass
            else:
                df[column] = df[column].replace(0, 'No')
                df[column] = df[column].replace(1, 'Yes')

        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
        if (analysis[0] == 'categorical_bar_chart'):
            plot_sns = charts.categorical_bar_chart(df_downsampled)
        return plot_sns

    @staticmethod
    def count_plot(select, analysis):
        dict = Convert(select)
        id_dict = {}
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        # print(df[df['ID'] == 'Zm00001eb125200'])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df)

        return plot_sns

    @staticmethod
    def downsampled_count_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_plots'):
            plot_sns = charts.Count_downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_scatter_plots'):
            plot_sns = charts.Count_downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'Count_down_dendogram_plots'):
            plot_sns = charts.Count_downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def Varionomic_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('varionomic', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            if column == 'distance':
                df[column] = df[column].abs()
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df)

        return plot_sns

    @staticmethod
    def downsampled_Varionomic_plot(select, analysis, cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('varionomic', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        print(df_downsampled.head())
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_plots'):
            plot_sns = charts.varionomic_downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_scatter_plots'):
            plot_sns = charts.downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'downsampled_dendogram_plots'):
            plot_sns = charts.downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def others_plot(select, analysis):
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            if column == 'distance':
                df[column] = df[column].abs()
            df[column] = df[column].replace('none', np.nan)
            df[column] = df[column].replace('none', np.nan)
            df[column] = df[column].replace('NA16', np.nan)
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('#VALUE!', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df)

        return plot_sns

    @staticmethod
    def downsampled_others_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('many_col_data', dict)
        list_cur = list(cursor)
        # print(list_cur)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            del dict['Origin']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('none', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('NA16', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('#VALUE!', np.nan)

        print(df_downsampled.head())
        print(df_downsampled.dtypes)
        if (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_plots'):
            plot_sns = charts.downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'downsampled_hierarchial_scatter_plots'):
            plot_sns = charts.downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'downsampled_dendogram_plots'):
            plot_sns = charts.downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Codon_Heatmap_plots'):
            plot_sns = charts.Codon_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Codon_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Codon_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns

    @staticmethod
    def TF_count_plot(select, analysis):
        dict = {}
        if  select[0] == "ARF":
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/ARF_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif select[0] == "BZIP":
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/BZIP_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif select[0] == "EREB" :
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/EREB_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif select[0] == "LBD":
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/LBD_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif select[0] == "SBP":
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/SBP_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        else:
            dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            # del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            # del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            # del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        if "no_label" in str(select):
            df = df
        else:
            df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)
        for column in df.columns:
            df[column] = df[column].replace('None', np.nan)
            df[column] = df[column].replace('None', np.nan)
        print(df.head())
        # print(df[df['dis_to_AF'].notna()].head())
        # return render_template('base.html')
        # print(df[df['ID'] == 'Zm00001eb125200'])
        if ((select[0]=="ARF" or select[0]=="BZIP" or select[0]=="EREB" or select[0]=="LBD" or select[0]=="SBP" )
             and analysis[0] == 'histogram'):
            plot_sns = charts.TF_histogram(df)
        elif (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df)
        return plot_sns

    @staticmethod
    def downsampled_TF_count_plot(select, analysis,cluster):
        cluster = cluster[0]
        dict = {}
        print(str(select))
        # for names in select:
        #     print(names)
        if select[0] == "ARF" :

            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/ARF_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
            dict = Convert(col_name)
            print(dict)
        elif select[0] == "BZIP":
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/BZIP_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
            dict = Convert(col_name)
        elif select[0] == "EREB" :
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/EREB_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
            dict = Convert(col_name)
        elif select[0] == "LBD":
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/LBD_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
            dict = Convert(col_name)
        elif select[0] == "SBP" :
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/SBP_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
            dict = Convert(col_name)
        else:
            dict = Convert(select)
        dict['ID'] = 1
        dict['_id'] = 0
        label_dict = {}
        print(dict)
        cursor = Database.find_by_structure('count', dict)
        list_cur = list(cursor)
        df = pd.DataFrame(list_cur)
        if "classical_label" in str(select):
            # del dict['classical_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['classical_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            # print(list_cur)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "core_label" in str(select):
            # del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['core_label'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        elif "Origin" in str(select):
            # del dict['core_label']
            label_dict['_id'] = 0
            label_dict['ID'] = 1
            label_dict['Origin'] = 1
            cursor_label = Database.find_by_structure('label', label_dict)
            list_cur = list(cursor_label)
            label_df = DataFrame(list_cur)
            label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(label_df.head())
        df.drop_duplicates(subset="ID", keep=False, inplace=True)
        df = pd.merge(left=df, right=label_df, how='left', left_on='ID', right_on='ID')

        print(df)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace('None', 0)
            print(df['classical_label'].value_counts())
            print(df['classical_label'].value_counts().argmax())
            df_majority = df[df.classical_label == df['classical_label'].value_counts().argmax()]
            df_minority = df[df.classical_label != df['classical_label'].value_counts().argmax()]
            counts = df['classical_label'].value_counts().tolist()
            df_downsampled = pd.DataFrame()
            # Downsample majority class
            df_majority_downsampled = resample(df_majority,
                                               replace=False,  # sample without replacement
                                               n_samples=min(counts),  # to match minority class
                                               random_state=123)
            df_downsampled = pd.concat([df_majority_downsampled, df_minority])
            print(df_downsampled.classical_label.value_counts())
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace('Core Gene', 1)
            df['core_label'] = df['core_label'].replace('Near-Core Gene', 2)
            df['core_label'] = df['core_label'].replace('Dispensable Gene', 3)
            df['core_label'] = df['core_label'].replace('Private Gene', 4)

            print(df['core_label'].value_counts())
            print(df['core_label'].value_counts().argmin())
            label = 'core_label'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.core_label.value_counts())
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace('WGD', 1)
            df['Origin'] = df['Origin'].replace('Tandem', 2)
            df['Origin'] = df['Origin'].replace('Both', 3)

            print(df['Origin'].value_counts())
            print(df['Origin'].value_counts().argmin())
            label = 'Origin'
            g = df.groupby(label, group_keys=False)
            df_downsampled = pd.DataFrame(g.apply(lambda x: x.sample(g.size().min(),random_state=123))).reset_index(drop=True)
            print(df_downsampled.Origin.value_counts())
        for column in df_downsampled.columns:
            if column == 'distance':
                df_downsampled[column] = df_downsampled[column].abs()
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)
            df_downsampled[column] = df_downsampled[column].replace('None', np.nan)

        # print(df_downsampled.head())
        # print(df_downsampled[df_downsampled['ID'] == 'Zm00001eb125200'])
        # print(df_downsampled[df_downsampled['Nc']== 0])
        if ((select[0]=="ARF" or select[0]=="BZIP" or select[0]=="EREB" or select[0]=="LBD" or select[0]=="SBP" )
             and analysis[0] == 'histogram'):
            plot_sns = charts.TF_histogram(df_downsampled)
        elif (analysis[0] == 'histogram'):
            plot_sns = charts.histogram(df_downsampled)
        elif (analysis[0] == 'Count_and_distribution'):
            plot_sns = charts.plots(df_downsampled)
        elif (analysis[0] == 'Pair_plots'):
            plot_sns = charts.pair_plots(df_downsampled)
        elif (analysis[0] == 'Box_plots'):
            plot_sns = charts.box_plots(df_downsampled)
        elif (analysis[0] == 'Violin_plots'):
            plot_sns = charts.violin_plots(df_downsampled)
        elif (analysis[0] == 'Joint_plots'):
            plot_sns = charts.joint_plots(df_downsampled)
        elif (analysis[0] == 'Scatter_plots'):
            plot_sns = charts.scatter_plots(df_downsampled)
        elif (analysis[0] == 'Correlation_plots'):
            plot_sns = charts.correlation_plots(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_plots'):
            plot_sns = charts.Count_downsampled_Hierarchial_heatmap(df_downsampled)
        elif (analysis[0] == 'Count_down_hierarchial_scatter_plots'):
            plot_sns = charts.Count_downsampled_hierarchial_scatter_plot(df_downsampled,cluster)
        elif (analysis[0] == 'Count_down_dendogram_plots'):
            plot_sns = charts.Count_downsampled_dendogram_plot(df_downsampled)
        elif (analysis[0] == 'Count_Heatmap_plots'):
            plot_sns = charts.Count_Heatmap_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_2D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_2D_samples_plots(df_downsampled)
        elif (analysis[0] == 'Count_Gene_PCA_3D_samples_plots'):
            plot_sns = charts.Count_Gene_PCA_3D_samples_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_label_plots'):
            plot_sns = charts.PCA_2D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_3D_label_plots'):
            plot_sns = charts.PCA_3D_label_plots(df_downsampled)
        elif (analysis[0] == 'PCA_2D_biplot_plots'):
            plot_sns = charts.PCA_2D_label_biplot_plots(df_downsampled)
        return plot_sns