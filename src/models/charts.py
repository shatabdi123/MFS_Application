import io
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import ttest_ind
from scipy.stats import chi2
from scipy import stats
from src.models.dendo_cluster import get_clust_graph
import plotly.figure_factory as ff
import dash_bio as dashbio
import plotly.express as px
from kmodes.kmodes import KModes
from matplotlib import colors as mcolors
from plotly.subplots import make_subplots
from scipy.spatial.distance import pdist, squareform
import plotly.graph_objects as go
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import base64
# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.cluster import AgglomerativeClustering

K = []


def downsampled_cluster_plot(df):
    rid = list(df['ID'])
    if "ID" in df.columns:
        df = df[df.columns].set_index('ID')
        # df = df.drop("ID", axis=1)
    if "classical_label" in df.columns:
        df = df.drop(["classical_label"], axis=1)
    if "core_label" in df.columns:
        df = df.drop(["core_label"], axis=1)
    if "Origin" in df.columns:
        df = df.drop(["Origin"], axis=1)
    cols = df.columns
    df = df.fillna(0)
    print(df.head())
    std_scaler = StandardScaler()
    df_label = std_scaler.fit_transform(df)
    df_label = pd.DataFrame(df_label, columns=cols)
    cluster_sns = get_clust_graph(df_label, rid, 786, transpose=True, dataname="Structure")
    return cluster_sns

rid = []
df_label = pd.DataFrame()
def count_downsampled_cluster_plot(df):
    rid = list(df['ID'])
    if "ID" in df.columns:
        df = df[df.columns].set_index('ID')
        # df = df.drop("ID", axis=1)
    if "classical_label" in df.columns:
        df = df.drop(["classical_label"], axis=1)
    if "core_label" in df.columns:
        df = df.drop(["core_label"], axis=1)
    if "Origin" in df.columns:
        df = df.drop(["Origin"], axis=1)
    df = df.fillna(0)
    df_label = df
    cluster_sns = get_clust_graph(df_label, rid, 786, transpose=True, dataname="Structure")
    return cluster_sns


class charts():

    @staticmethod
    def histogram(df):
        print(df.shape)
        dataset_bin = pd.DataFrame()  # To contain our dataframe with our discretised continuous variables
        fig_h = 0
        fig_w = 0
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            dataset_bin['classical_label'] = df['classical_label']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 1250
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 1250
            else:
                fig_h = 3000
                fig_w = 1250
        elif "core_label" in df.columns:
            dataset_bin['core_label'] = df['core_label']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 2000
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 2000
            else:
                fig_h = 3000
                fig_w = 2000
        elif "Origin" in df.columns:
            dataset_bin['Origin'] = df['Origin']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 2000
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 2000
            else:
                fig_h = 3000
                fig_w = 2000
        else:
            if len(df.columns) == 2:
                fig_h = 1000
                fig_w = 1250
            elif len(df.columns) == 3:
                fig_h = 1800
                fig_w = 1250
            else:
                fig_h = 2800
                fig_w = 1250

        if "classical_label" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=2)
        elif "core_label" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=4)
        elif "Origin" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=3)
        else:
            fig = make_subplots(rows=len(df.columns), cols=1)
        for i, column in enumerate(df.columns, start=1):
            data_1_count = pd.DataFrame()
            data_0_count = pd.DataFrame()
            print(i)
            for j in range(1, 2):
                print(j)
                print([i], [j])
                if column.startswith("WGCNA") or column.startswith("classical_label") or column.startswith(
                        "core_label")  or column.startswith("Origin"):
                    pass
                    print(column)
                else:
                    if "classical_label" in dataset_bin.columns:
                        data_1 = df.loc[df['classical_label'] == 1.0][column].dropna()
                        data_1_count[column] = df.loc[df['classical_label'] == 1.0][column].dropna()
                        data_1_count = data_1_count.groupby(column)[column].count().reset_index(name="counts")

                        mu_1, std_1 = norm.fit(data_1)
                        data_0 = df.loc[df['classical_label'] == 0.0][column].dropna()
                        data_0_count[column] = df.loc[df['classical_label'] == 0.0][column].dropna()
                        data_0_count = data_0_count.groupby(column)[column].count().reset_index(name="counts")
                        print(data_1.head())
                        print(data_1_count.tail())
                        # hover_label = column + ' Count'
                        # data_1_mean = np.mean(data_1)
                        # data_0_mean = np.mean(data_0)
                        ttest, pval = ttest_ind(data_1, data_0)
                        print("p-value", pval)
                        hover_label = column
                        print(np.max(list(data_1_count["counts"])))
                        print(np.max(list(data_0_count["counts"])))
                        counts_1 = plt.hist(data_1, bins=50)
                        counts_0 = plt.hist(data_0, bins=50)
                        print(np.max(counts_1[0]))
                        if np.max(counts_1[0]) > np.max(counts_0[0]):
                            max = np.max(counts_1[0])
                        else:
                            max = np.max(counts_0[0])

                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         name='Classical genes',
                                         nbinsx=100,
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}',
                                         marker=dict(color='#1f77b4'),
                                         ),
                            row=i, col=j
                        )

                        label_1 = "<b>Classical genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        # fig.update_xaxes(range=[0, max(data_1[column])])
                        fig.update_xaxes(title_text=hover_label, title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),row=i, col=j)
                        print("the max is",max)
                        print("the log max",math.log(max,10))

                        fig.update_yaxes(type="log", range=[-1, math.log(max, 10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         title_font=dict(size=20, family='Arial Black', color='black'), row=i,col=j)

                        mu_0, std_0 = norm.fit(data_0)
                        label_0 = "<b>Other genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_0, std_0)
                        fig.add_trace(
                            go.Histogram(x=data_0,
                                         nbinsx=100,
                                         name='Other genes',
                                         marker=dict(color='#ff7f0e'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 1
                        )

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_0, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 1)
                        # fig.update_layout(hovermode="x unified")
                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 1)

                        fig.update_yaxes(type="log", range=[-1, math.log(max, 10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),row=i,col=j + 1)

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E  </b>" % (column,pval)
                            print("p_val len", len(p_val),type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 2.1
                            elif (int(len(p_val)) >= 130) and (len(p_val) < 150) :

                                x_len = 2.2
                            elif int(len(p_val)) >= 150:

                                x_len = 2.3
                            else:

                                x_len = 1.9
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=1700,legend=dict(title_font_family="Arial Black",
                              font=dict(size=18, family='Arial Black', color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=1700,legend=dict(title_font_family="Arial Black",
                              font=dict(size=18, family='Arial Black', color='black')))

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E  </b>" % (column,pval)

                            print("p_val len", len(p_val),type(len(p_val)))
                            if(len(p_val)> 110) and (len(p_val) < 130):

                                x_len = 2.1
                            elif int(len(p_val)) >= 130 and (len(p_val) < 150) :

                                x_len = 2.2
                            elif int(len(p_val)) >= 150:

                                x_len = 2.3
                            else:

                                x_len = 1.9

                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)

                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=1700,legend=dict(title_font_family="Arial Black",
                              font=dict(size=18, family='Arial Black', color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=1700,legend=dict(title_font_family="Arial Black",
                              font=dict(size=18, family='Arial Black', color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            print("p_val len", len(p_val),type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 130):

                                x_len = 1.9
                            elif (int(len(p_val)) >= 130) and (len(p_val) < 150) :

                                x_len = 1.9
                            elif int(len(p_val)) >= 150:

                                x_len = 1.9
                            else:

                                x_len = 1.9
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=1700,legend=dict(title_font_family="Arial Black",
                              font=dict(size=18, family='Arial Black', color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=1700,legend=dict(title_font_family="Arial Black",
                              font=dict(size=18, family='Arial Black', color='black')))

                    elif "core_label" in dataset_bin.columns:
                        data_1 = df.loc[df['core_label'] == 1.0][column].dropna()
                        data_2 = df.loc[df['core_label'] == 2.0][column].dropna()
                        data_3 = df.loc[df['core_label'] == 3.0][column].dropna()
                        data_4 = df.loc[df['core_label'] == 4.0][column].dropna()
                        # hover_label = column + ' Count'
                        hover_label = column
                        anovatest, pval = stats.f_oneway(data_1, data_2,data_3, data_4)
                        print("p-value", pval)
                        counts_1 = plt.hist(data_1, bins=50)
                        counts_2 = plt.hist(data_2, bins=50)
                        counts_3 = plt.hist(data_3, bins=50)
                        counts_4 = plt.hist(data_4, bins=50)
                        max_list = []
                        max_list.append(np.max(counts_1[0]))
                        max_list.append(np.max(counts_2[0]))
                        max_list.append(np.max(counts_3[0]))
                        max_list.append(np.max(counts_4[0]))
                        print(max_list)
                        maximum = np.max(max_list)
                        print(maximum)
                        mu_1, std_1 = norm.fit(data_1)
                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         nbinsx=100,
                                         name='Core genes',
                                         marker=dict(color='#1f77b4'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j
                        )
                        label_1 = "<b>Core genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black')
                                         ,row=i, col=j)

                        fig.update_yaxes(type="log", range=[-1,math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        mu_2, std_2 = norm.fit(data_2)
                        label_2 = "<b>Near-core genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_2, std_2)
                        fig.add_trace(
                            go.Histogram(x=data_2,
                                         nbinsx=100,
                                         name='Near-core genes',
                                         marker=dict(color='#ff7f0e'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 1
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_2, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 1)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 1)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 1)

                        mu_3, std_3 = norm.fit(data_3)
                        fig.add_trace(
                            go.Histogram(x=data_3,
                                         nbinsx=100,
                                         name='Dispensable genes',
                                         marker=dict(color='#2ca02c'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),

                            row=i, col=j + 2
                        )
                        label_3 = "<b>Dispensable genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_3, std_3)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_3, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 2)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 2)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 2)
                        mu_4, std_4 = norm.fit(data_4)
                        label_4 = "<b>Private genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_4, std_4)
                        fig.add_trace(
                            go.Histogram(x=data_4,
                                         nbinsx=100,
                                         name='Private Genes',
                                         marker=dict(color='#d62728'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 3
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_4, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 3)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 3)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 3)
                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in the different core genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            print("len of p_val",len(p_val))
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 3.3
                            elif int(len(p_val)) >= 130 and int(len(p_val)) < 150:

                                x_len = 3.6
                            elif int(len(p_val)) >= 150:
                                x_len = 3.7
                            else:
                                x_len = 3.00
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=3300, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=3300, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in the different core genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            print("len of p_val", len(p_val))
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 3.3
                            elif int(len(p_val)) >= 130 and int(len(p_val)) < 150:

                                x_len = 3.6
                            elif int(len(p_val)) >= 150:
                                x_len = 3.7
                            else:
                                x_len = 3.00
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=3300, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=3300, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                        else:
                            # print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            print("len of p_val",len(p_val))
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 3.3
                            elif int(len(p_val)) >= 130 and int(len(p_val)) < 150:

                                x_len = 3.6
                            elif int(len(p_val)) >= 150:
                                x_len = 3.7
                            else:
                                x_len = 3.00
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=3300, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=3300, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))


                    elif "Origin" in dataset_bin.columns:
                        data_1 = df.loc[df['Origin'] == 1.0][column].dropna()
                        data_2 = df.loc[df['Origin'] == 2.0][column].dropna()
                        data_3 = df.loc[df['Origin'] == 3.0][column].dropna()
                        # hover_label = column + ' Count'
                        hover_label = column
                        anovatest, pval = stats.f_oneway(data_1, data_2, data_3)
                        print("p-value", pval)
                        counts_1 = plt.hist(data_1, bins=50)
                        counts_2 = plt.hist(data_2, bins=50)
                        counts_3 = plt.hist(data_3, bins=50)
                        max_list = []
                        max_list.append(np.max(counts_1[0]))
                        max_list.append(np.max(counts_2[0]))
                        max_list.append(np.max(counts_3[0]))
                        print(max_list)
                        maximum = np.max(max_list)
                        print(maximum)
                        mu_1, std_1 = norm.fit(data_1)
                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         nbinsx=100,
                                         name='WGD genes',
                                         marker=dict(color='#1f77b4'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j
                        )
                        label_1 = "<b>WGD genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        fig.update_xaxes(title_text=hover_label ,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                        row=i, col=j)

                        fig.update_yaxes(type="log", range=[-1, math.log(maximum, 10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j)

                        mu_2, std_2 = norm.fit(data_2)
                        label_2 = "<b>Tandem genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_2, std_2)
                        fig.add_trace(
                            go.Histogram(x=data_2,
                                         nbinsx=100,
                                         name='Tandem genes',
                                         marker=dict(color='#ff7f0e'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 1
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_2, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 1)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 1)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum, 10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j + 1)
                        mu_3, std_3 = norm.fit(data_3)
                        fig.add_trace(
                            go.Histogram(x=data_3,
                                         nbinsx=100,
                                         name='Both genes',
                                         marker=dict(color='#2ca02c'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),

                            row=i, col=j + 2
                        )
                        label_3 = "<b>Both genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_3, std_3)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_3, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 2)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 2)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum, 10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j + 2)

                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            print("len of p_val", len(p_val))
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 2.6
                            elif int(len(p_val)) >= 130 and int(len(p_val)) < 150:

                                x_len = 2.8
                            elif int(len(p_val)) >= 150:
                                x_len = 2.9
                            else:
                                x_len = 2.4
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=2700, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=2700, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            print("len of p_val", len(p_val))
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 2.6
                            elif int(len(p_val)) >= 130 and int(len(p_val)) < 150:

                                x_len = 2.8
                            elif int(len(p_val)) >= 150:
                                x_len = 2.9
                            else:
                                x_len = 2.4
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)

                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=2700, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=2700, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None</b>"
                            print("len of p_val", len(p_val))
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 2.6
                            elif int(len(p_val)) >= 130 and int(len(p_val)) < 150:

                                x_len = 2.8
                            elif int(len(p_val)) >= 150:
                                x_len = 2.9
                            else:
                                x_len = 2.4
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                fig.update_layout(height=fig_h, width=2700, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))
                            else:
                                fig.update_layout(height=fig_h, width=2700, legend=dict(title_font_family="Arial Black",
                                                                                        font=dict(size=18,
                                                                                                  family='Arial Black',
                                                                                                  color='black')))

                    else:
                        data_1 = df[column].dropna()
                        data_1_count[column] = df[column].dropna()
                        data_1_count = data_1_count.groupby(column)[column].count().reset_index(name="counts")

                        mu_1, std_1 = norm.fit(data_1)
                        hover_label = column

                        print(data_1.head())
                        print(data_1_count.tail())


                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         nbinsx=100,
                                         name='Genes',
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}',
                                         marker=dict(color='#1f77b4'),
                                         ),
                            row=i, col=j
                        )

                        label_1 = "<b>All genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        # fig.update_xaxes(range=[0, max(data_1[column])])
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)

                        fig.update_yaxes(type="log", title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j)

                        fig.update_layout(showlegend=False)

                        if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                            fig.update_layout(height=fig_h, width=800)
                        else:
                            fig.update_layout(height=fig_h, width=800)
        return fig
        # fig.show()

    @staticmethod
    def TF_histogram(df):
        print(df.shape)
        dataset_bin = pd.DataFrame()  # To contain our dataframe with our discretised continuous variables
        fig_h = 0
        fig_w = 0
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "classical_label" in df.columns:
            dataset_bin['classical_label'] = df['classical_label']
            # if len(df.columns) == 2:
            fig_h = 1200
            fig_w = 1500

        elif "core_label" in df.columns:
            dataset_bin['core_label'] = df['core_label']
            # if len(df.columns) == 2:
            fig_h = 1200
            fig_w = 3000

        elif "Origin" in df.columns:
            dataset_bin['Origin'] = df['Origin']
            fig_h = 1200
            fig_w = 2800
        else:
            fig_h = 1200
            fig_w = 800
        for col in df.columns:
            if "ARF" in col:
                TF_col = "ARF"
            elif "BZIP" in col:
                TF_col = "BZIP"
            elif "EREB" in col:
                TF_col = "EREB"
            elif "LBD" in col:
                TF_col = "LBD"
            elif "SBP" in col:
                TF_col = "SBP"

        df[TF_col] = df.sum(axis=1)

        if "classical_label" in df.columns:
            df = df[[TF_col, "classical_label"]]
            fig = make_subplots(rows=len(df.columns), cols=2)
            print(len(df.columns))
        elif "core_label" in df.columns:
            df = df[[TF_col, "core_label"]]
            fig = make_subplots(rows=len(df.columns), cols=4)
        elif "Origin" in df.columns:
            df = df[[TF_col, "Origin"]]
            fig = make_subplots(rows=len(df.columns), cols=3)
        else:
            df = df[[TF_col]]
            fig = make_subplots(rows=len(df.columns), cols=1)


        for i, column in enumerate(df.columns, start=1):
            data_1_count = pd.DataFrame()
            data_0_count = pd.DataFrame()
            for j in range(1, 2):
                print(j)
                print([i], [j])
                if column.startswith("WGCNA") or column.startswith("classical_label") or column.startswith(
                        "core_label")  or column.startswith("Origin") :
                    pass
                    print(column)
                else:
                    if "classical_label" in dataset_bin.columns:
                        data_1 = df.loc[df['classical_label'] == 1.0][column].dropna()
                        data_1_count[column] = df.loc[df['classical_label'] == 1.0][column].dropna()
                        data_1_count = data_1_count.groupby(column)[column].count().reset_index(name="counts")
                        mu_1, std_1 = norm.fit(data_1)
                        data_0 = df.loc[df['classical_label'] == 0.0][column].dropna()
                        data_0_count[column] = df.loc[df['classical_label'] == 0.0][column].dropna()
                        data_0_count = data_0_count.groupby(column)[column].count().reset_index(name="counts")
                        # hover_label = column + ' Count'
                        ttest, pval = ttest_ind(data_1, data_0)
                        print("p-value", pval)
                        hover_label = column
                        print(np.max(list(data_1_count["counts"])))
                        print(np.max(list(data_0_count["counts"])))
                        counts_1 = plt.hist(data_1, bins=50)
                        counts_0 = plt.hist(data_0, bins=50)
                        print(np.max(counts_1[0]))
                        if np.max(counts_1[0]) > np.max(counts_0[0]):
                            max = np.max(counts_1[0])
                        else:
                            max = np.max(counts_0[0])
                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         name='Classical genes',
                                         nbinsx=100,
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}',
                                         marker=dict(color='#1f77b4'),
                                         ),
                            row=i, col=j
                        )

                        label_1 = "<b>Classical genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        # fig.update_xaxes(range=[0, max(data_1[column])])
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)
                        print("the max is", max)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        mu_0, std_0 = norm.fit(data_0)
                        label_0 = "<b>Other genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_0, std_0)
                        fig.add_trace(
                            go.Histogram(x=data_0,
                                         nbinsx=100,
                                         name='Other genes',
                                         marker=dict(color='#ff7f0e'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 1
                        )

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_0, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 1)
                        # fig.update_layout(hovermode="x unified")
                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 1)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 1)

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            fig.add_annotation(xref="x domain", yref="y domain", x=2.13, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            fig.add_annotation(xref="x domain", yref="y domain", x=2.13, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            fig.add_annotation(xref="x domain", yref="y domain", x=2.13, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))

                    elif "core_label" in dataset_bin.columns:
                        data_1 = df.loc[df['core_label'] == 1.0][column].dropna()
                        data_2 = df.loc[df['core_label'] == 2.0][column].dropna()
                        data_3 = df.loc[df['core_label'] == 3.0][column].dropna()
                        data_4 = df.loc[df['core_label'] == 4.0][column].dropna()
                        # hover_label = column + ' Count'
                        hover_label = column
                        anovatest, pval = stats.f_oneway(data_1, data_2, data_3, data_4)
                        print("p-value", pval)
                        counts_1 = plt.hist(data_1, bins=50)
                        counts_2 = plt.hist(data_2, bins=50)
                        counts_3 = plt.hist(data_3, bins=50)
                        counts_4 = plt.hist(data_4, bins=50)
                        max_list = []
                        max_list.append(np.max(counts_1[0]))
                        max_list.append(np.max(counts_2[0]))
                        max_list.append(np.max(counts_3[0]))
                        max_list.append(np.max(counts_4[0]))
                        print(max_list)
                        maximum = np.max(max_list)
                        print(maximum)

                        mu_1, std_1 = norm.fit(data_1)
                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         nbinsx=100,
                                         name='Core genes',
                                         marker=dict(color='#1f77b4'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j
                        )
                        label_1 = "<b>Core genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)
                        mu_2, std_2 = norm.fit(data_2)
                        label_2 = "<b>Near-core genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_2, std_2)
                        fig.add_trace(
                            go.Histogram(x=data_2,
                                         nbinsx=100,
                                         name='Near-core genes',
                                         marker=dict(color='#ff7f0e'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 1
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_2, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 1)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 1)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j + 1)

                        mu_3, std_3 = norm.fit(data_3)
                        fig.add_trace(
                            go.Histogram(x=data_3,
                                         nbinsx=100,
                                         name='Dispensable genes',
                                         marker=dict(color='#2ca02c'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),

                            row=i, col=j + 2
                        )
                        label_3 = "<b>Dispensable genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_3, std_3)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_3, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 2)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 2)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j + 2)

                        mu_4, std_4 = norm.fit(data_4)
                        label_4 = "<b>Private genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_4, std_4)
                        fig.add_trace(
                            go.Histogram(x=data_4,
                                         nbinsx=100,
                                         name='Private Genes',
                                         marker=dict(color='#d62728'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 3
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_4, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 3)
                        fig.update_xaxes(title_text=hover_label, row=i, col=j + 3)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j + 3)

                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))
                        if pval < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in the different core genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            fig.add_annotation(xref="x domain", yref="y domain", x=3.1, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in the different core genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            fig.add_annotation(xref="x domain", yref="y domain", x=3.1, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            fig.add_annotation(xref="x domain", yref="y domain", x=3.1, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))

                    elif "Origin" in dataset_bin.columns:
                        data_1 = df.loc[df['Origin'] == 1.0][column].dropna()
                        data_2 = df.loc[df['Origin'] == 2.0][column].dropna()
                        data_3 = df.loc[df['Origin'] == 3.0][column].dropna()
                        hover_label = column
                        anovatest, pval = stats.f_oneway(data_1, data_2, data_3)
                        print("p-value", pval)
                        counts_1 = plt.hist(data_1, bins=50)
                        counts_2 = plt.hist(data_2, bins=50)
                        counts_3 = plt.hist(data_3, bins=50)
                        max_list = []
                        max_list.append(np.max(counts_1[0]))
                        max_list.append(np.max(counts_2[0]))
                        max_list.append(np.max(counts_3[0]))
                        print(max_list)
                        maximum = np.max(max_list)
                        print(maximum)

                        mu_1, std_1 = norm.fit(data_1)
                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         nbinsx=100,
                                         name='WGD genes',
                                         marker=dict(color='#1f77b4'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j
                        )
                        label_1 = "<b>WGD genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=19
                            ), row=i, col=j)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        mu_2, std_2 = norm.fit(data_2)
                        label_2 = "<b>Tandem genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_2, std_2)
                        fig.add_trace(
                            go.Histogram(x=data_2,
                                         nbinsx=100,
                                         name='Tandem genes',
                                         marker=dict(color='#ff7f0e'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),
                            row=i, col=j + 1
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_2, font=dict(
                                family="Arial Black",
                                size=19
                            ), row=i, col=j + 1)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 1)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j + 1)
                        mu_3, std_3 = norm.fit(data_3)
                        fig.add_trace(
                            go.Histogram(x=data_3,
                                         nbinsx=100,
                                         name='Both genes',
                                         marker=dict(color='#2ca02c'),
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                         ),

                            row=i, col=j + 2
                        )
                        label_3 = "<b>Both genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_3, std_3)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_3, font=dict(
                                family="Arial Black",
                                size=19
                            ), row=i, col=j + 2)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 2)
                        fig.update_yaxes(type="log", range=[-1, math.log(maximum,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i,col=j + 2)


                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))
                        if pval < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            fig.add_annotation(xref="x domain", yref="y domain", x=2.4, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Rockwell",
                                    size=18,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E  </b>" % (column, pval)
                            fig.add_annotation(xref="x domain", yref="y domain", x=2.4, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Rockwell",
                                    size=18,
                                    color="Crimson"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            fig.add_annotation(xref="x domain", yref="y domain", x=2.4, y=1.17, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Rockwell",
                                    size=18,
                                    color="purple"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))

                    else:
                        data_1 = df[column].dropna()
                        data_1_count[column] = df[column].dropna()
                        data_1_count = data_1_count.groupby(column)[column].count().reset_index(name="counts")
                        mu_1, std_1 = norm.fit(data_1)


                        hover_label = column

                        fig.add_trace(
                            go.Histogram(x=data_1,
                                         name='Genes',
                                         nbinsx=100,
                                         hovertemplate=column + ': %{x} <br>Number of genes: %{y}',
                                         marker=dict(color='#1f77b4'),
                                         ),
                            row=i, col=j
                        )

                        label_1 = "<b>All genes fit results: mu = %.2f,  std = %.2f </b>" % (mu_1, std_1)

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.1, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=19
                            ), row=i, col=j)
                        # fig.update_xaxes(range=[0, max(data_1[column])])
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_layout(showlegend=False)
                        fig.update_yaxes(type="log", title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)
                        fig.update_layout(height=fig_h, width=fig_w)

        return fig


    @staticmethod
    def categorical_bar_chart(df):
        dataset_bin = pd.DataFrame()  # To contain our dataframe with our discretised continuous variables
        fig_h = 0
        fig_w = 0
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "classical_label" in df.columns:
            dataset_bin['classical_label'] = df['classical_label']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 1500
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 1500
            else:
                fig_h = 3400
                fig_w = 1500
        elif "core_label" in df.columns:
            dataset_bin['core_label'] = df['core_label']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 3000
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 3000
            else:
                fig_h = 3400
                fig_w = 3000
        elif "Origin" in df.columns:
            dataset_bin['Origin'] = df['Origin']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 2700
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 2700
            else:
                fig_h = 3400
                fig_w = 2700
        else:
            if len(df.columns) == 2:
                fig_h = 1000
                fig_w = 1250
            elif len(df.columns) == 3:
                fig_h = 1800
                fig_w = 1250
            else:
                fig_h = 2800
                fig_w = 1250
        if "classical_label" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=2)
        elif "core_label" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=4)
        elif "Origin" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=3)
        else:
            fig = make_subplots(rows=len(df.columns), cols=1)
        for i, column in enumerate(df.columns, start=1):
            data_0 = pd.DataFrame()
            data_1 = pd.DataFrame()
            data_2 = pd.DataFrame()
            data_3 = pd.DataFrame()
            data_4 = pd.DataFrame()
            print(i)
            for j in range(1,2):
                print(j)
                print([i],[j])
                if column.startswith("classical_label") or column.startswith(
                        "core_label") or column.startswith("Origin"):
                    pass
                    print(column)
                else:
                    if "classical_label" in dataset_bin.columns:
                        data_1[column] = df.loc[df['classical_label'] == 1.0][column].dropna()
                        data_1 = data_1.groupby(column)[column].count().reset_index(name="counts")
                        data_0[column] = df.loc[df['classical_label'] == 0.0][column].dropna()
                        data_0 = data_0.groupby(column)[column].count().reset_index(name="counts")
                        print(data_0)
                        print(data_1)
                        print("count : ", np.max([np.max(data_0['counts']),np.max(data_1['counts'])]))
                        # print("count : ", np.max(data_0_count['count']))
                        max = np.max([np.max(data_0['counts']),np.max(data_1['counts'])])
                        print("the max:",max)
                        contingency_table = pd.crosstab(df[column], df['classical_label'])
                        print('contingency_table :-\n', contingency_table)
                        # Observed Values
                        Observed_Values = contingency_table.values
                        print("Observed Values :-\n", Observed_Values)
                        b = stats.chi2_contingency(contingency_table)
                        Expected_Values = b[3]
                        print("Expected Values :-\n", Expected_Values)

                        no_of_rows = contingency_table.shape[0]
                        print("shape :", contingency_table.shape)
                        no_of_columns = contingency_table.shape[1]
                        ddof = (no_of_rows - 1) * (no_of_columns - 1)
                        print("Degree of Freedom:-", ddof)
                        alpha = 0.05
                        chi_square = sum([(o - e) ** 2. / e for o, e in zip(Observed_Values, Expected_Values)])
                        chi_square_statistic = chi_square[0] + chi_square[1]
                        print("chi-square statistic:-", chi_square_statistic)
                        critical_value = chi2.ppf(q=1 - alpha, df=ddof)
                        print('critical_value:', critical_value)
                        # p-value
                        p_value = 1 - chi2.cdf(x=chi_square_statistic, df=ddof)

                        hover_label = column
                        fig.add_trace(
                            go.Bar(x=data_1[column], y=data_1['counts'],
                                   name='Classical genes',
                                   marker=dict(color='#1f77b4'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),

                            row=i, col=j
                        )
                        label_1 = "<b>Classical genes</b>"

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_1, font=dict(
            family="Arial Black",
            size=20
            ), row=i, col=j)
                        # fig.update_xaxes(range=[0, max(data_1[column])])
                        fig.update_xaxes(title_text= hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        label_0 = "<b>Other genes</b>"
                        fig.add_trace(
                            go.Bar(x=data_0[column], y=data_0['counts'],
                                   name='Other genes',
                                   marker=dict(color='#ff7f0e'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),
                            row=i, col=j+1
                        )

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_0, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j+1)
                        # fig.update_layout(hovermode="x unified")
                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j+1)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 1)

                        if p_value < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E  </b>" % (column, p_value)
                            if (len(p_val) > 110) and (len(p_val) < 130):

                                x_len = 2.3
                            elif int(len(p_val)) >= 130:

                                x_len = 2.5
                            else:

                                x_len = 2.13
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                            fig.update_xaxes(tickangle=45)

                        elif p_value > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E  </b>" % (column, p_value)
                            if (len(p_val) > 110) and (len(p_val) < 130):

                                x_len = 2.3
                            elif int(len(p_val)) >= 130:

                                x_len = 2.5
                            else:

                                x_len = 2.13
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                            fig.update_xaxes(tickangle=45)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            if (len(p_val) > 110) and (len(p_val) < 130):

                                x_len = 2.3
                            elif int(len(p_val)) >= 130:

                                x_len = 2.5
                            else:

                                x_len = 2.13
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                            fig.update_xaxes(tickangle=45)

                    elif "core_label" in dataset_bin.columns:
                        data_1[column] = df.loc[df['core_label'] == 1.0][column].dropna()
                        data_1 = data_1.groupby(column)[column].count().reset_index(name="counts")
                        data_2[column] = df.loc[df['core_label'] == 2.0][column].dropna()
                        data_2 = data_2.groupby(column)[column].count().reset_index(name="counts")
                        data_3[column] = df.loc[df['core_label'] == 3.0][column].dropna()
                        data_3 = data_3.groupby(column)[column].count().reset_index(name="counts")
                        data_4[column] = df.loc[df['core_label'] == 4.0][column].dropna()
                        data_4 = data_4.groupby(column)[column].count().reset_index(name="counts")
                        hover_label = column + ' Count'
                        max = np.max([np.max(data_1['counts']), np.max(data_2['counts']),np.max(data_3['counts']),
                                      np.max(data_4['counts'])])
                        contingency_table = pd.crosstab(df[column], df['core_label'])
                        print('contingency_table :-\n', contingency_table)
                        # Observed Values
                        Observed_Values = contingency_table.values
                        print("Observed Values :-\n", Observed_Values)
                        b = stats.chi2_contingency(contingency_table)
                        Expected_Values = b[3]
                        print("Expected Values :-\n", Expected_Values)

                        no_of_rows = contingency_table.shape[0]
                        print(contingency_table.shape)
                        no_of_columns = contingency_table.shape[1]
                        print(no_of_rows,no_of_columns)
                        ddof = (no_of_rows - 1) * (no_of_columns - 1)
                        print("Degree of Freedom:-", ddof)
                        alpha = 0.05
                        chi_square = sum([(o - e) ** 2. / e for o, e in zip(Observed_Values, Expected_Values)])
                        chi_square_statistic = chi_square[0] + chi_square[1]
                        print("chi-square statistic:-", chi_square_statistic)
                        critical_value = chi2.ppf(q=1 - alpha, df=ddof)
                        print('critical_value:', critical_value)
                        # p-value
                        p_value = 1 - chi2.cdf(x=chi_square_statistic, df=ddof)
                        fig.add_trace(
                            go.Bar(x=data_1[column], y=data_1['counts'],
                                   name='Core genes',
                                   marker=dict(color='#1f77b4'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),

                            row=i, col=j
                        )
                        label_1 = "<b>Core genes </b>"

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        label_2 = "<b>Near-core genes</b>"
                        fig.add_trace(
                            go.Bar(x=data_2[column], y=data_2['counts'],
                                   name='Near-core genes',
                                   marker=dict(color='#ff7f0e'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),
                            row=i, col=j + 1
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_2, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 1)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 1)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j+1)
                        fig.add_trace(
                            go.Bar(x=data_3[column], y=data_3['counts'],
                                   name='Dispensable genes',
                                   marker=dict(color='#2ca02c'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),

                            row=i, col=j+2
                        )
                        label_3 = "<b>Dispensable genes </b>"

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_3, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j+2)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j+2)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 2)

                        label_4 = "<b>Private genes  </b>"
                        fig.add_trace(
                            go.Bar(x=data_4[column], y=data_4['counts'],
                                   name='Private Genes',
                                   marker=dict(color='#d62728'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),
                            row=i, col=j + 3
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_4, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 3)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 3)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 3)

                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))
                        # fig.update_layout(height=fig_h, width=fig_w)
                        if p_value < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in the different core genes" \
                                    " as p-value = %.2E  </b>" % (column, p_value)
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 3.3
                            elif int(len(p_val)) >= 130:

                                x_len = 3.6
                            else:
                                x_len = 3.00
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            fig.update_xaxes(tickangle=45)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))

                        elif p_value > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in the different core genes" \
                                    " as p-value = %.2E  </b>" % (column, p_value)
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 3.3
                            elif int(len(p_val)) >= 130:

                                x_len = 3.6
                            else:
                                x_len = 3.00
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                            fig.update_xaxes(tickangle=45)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None</b>"
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 3.3
                            elif int(len(p_val)) >= 130:

                                x_len = 3.6
                            else:
                                x_len = 3.00
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                            fig.update_xaxes(tickangle=45)

                    elif "Origin" in dataset_bin.columns:
                        data_1[column] = df.loc[df['Origin'] == 1.0][column].dropna()
                        data_1 = data_1.groupby(column)[column].count().reset_index(name="counts")
                        data_2[column] = df.loc[df['Origin'] == 2.0][column].dropna()
                        data_2 = data_2.groupby(column)[column].count().reset_index(name="counts")
                        data_3[column] = df.loc[df['Origin'] == 3.0][column].dropna()
                        data_3 = data_3.groupby(column)[column].count().reset_index(name="counts")
                        hover_label = column + ' Count'
                        max = np.max([np.max(data_1['counts']), np.max(data_2['counts']),np.max(data_3['counts'])])
                        contingency_table = pd.crosstab(df[column], df['Origin'])
                        print('contingency_table :-\n', contingency_table)
                        # Observed Values
                        Observed_Values = contingency_table.values
                        print("Observed Values :-\n", Observed_Values)
                        b = stats.chi2_contingency(contingency_table)
                        Expected_Values = b[3]
                        print("Expected Values :-\n", Expected_Values)

                        no_of_rows = contingency_table.shape[0]
                        print(contingency_table.shape)
                        no_of_columns = contingency_table.shape[1]
                        print(no_of_rows,no_of_columns)
                        ddof = (no_of_rows - 1) * (no_of_columns - 1)
                        print("Degree of Freedom:-", ddof)
                        alpha = 0.05
                        chi_square = sum([(o - e) ** 2. / e for o, e in zip(Observed_Values, Expected_Values)])
                        chi_square_statistic = chi_square[0] + chi_square[1]
                        print("chi-square statistic:-", chi_square_statistic)
                        critical_value = chi2.ppf(q=1 - alpha, df=ddof)
                        print('critical_value:', critical_value)
                        # p-value
                        p_value = 1 - chi2.cdf(x=chi_square_statistic, df=ddof)
                        fig.add_trace(
                            go.Bar(x=data_1[column], y=data_1['counts'],
                                   name='WGD genes',
                                   marker=dict(color='#1f77b4'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),

                            row=i, col=j
                        )
                        label_1 = "<b>WGD genes </b>"

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        label_2 = "<b>Tandem genes</b>"
                        fig.add_trace(
                            go.Bar(x=data_2[column], y=data_2['counts'],
                                   name='Tandem genes',
                                   marker=dict(color='#ff7f0e'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),
                            row=i, col=j + 1
                        )
                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_2, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j + 1)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j + 1)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j+1)

                        fig.add_trace(
                            go.Bar(x=data_3[column], y=data_3['counts'],
                                   name='Both genes',
                                   marker=dict(color='#2ca02c'),
                                   hovertemplate= column + ': %{x} <br>Number of genes: %{y}'
                                   ),

                            row=i, col=j+2
                        )
                        label_3 = "<b>Both genes </b>"

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_3, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j+2)
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j+2)
                        fig.update_yaxes(type="log", range=[-1, math.log(max,10)], title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j + 2)

                        fig.update_layout(
                            hoverlabel=dict(
                                bgcolor="white",
                                font_size=16,
                                font_family="Rockwell"
                            ))
                        # fig.update_layout(height=fig_h, width=fig_w)
                        if p_value < 0.05:
                            print("we reject null hypothesis")
                            p_val = "<b>There is significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E  </b>" % (column, p_value)
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 2.7
                            elif int(len(p_val)) >= 130:

                                x_len = 2.8
                            else:
                                x_len = 2.5
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                            fig.update_xaxes(tickangle=45)

                        elif p_value > 0.05:
                            print("we accept null hypothesis")
                            p_val = "<b>There is no significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E  </b>" % (column, p_value)
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 2.7
                            elif int(len(p_val)) >= 130:

                                x_len = 2.9
                            else:
                                x_len = 2.5
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                            fig.update_xaxes(tickangle=45)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            if (len(p_val) > 110) and (len(p_val) < 130):
                                x_len = 2.7
                            elif int(len(p_val)) >= 130:

                                x_len = 2.8
                            else:
                                x_len = 2.5
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.11, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                            fig.update_layout(height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                            font=dict(size=18, family='Arial Black', color='black')))
                            fig.update_xaxes(tickangle=45)
                    else:
                        data_1[column] = df[column].dropna()
                        data_1 = data_1.groupby(column)[column].count().reset_index(name="counts")

                        print(data_1)


                        hover_label = column
                        fig.add_trace(
                            go.Bar(x=data_1[column], y=data_1['counts'],
                                   name='Genes',
                                   marker=dict(color='#1f77b4'),
                                   hovertemplate=column + ': %{x} <br>Number of genes: %{y}'
                                   ),

                            row=i, col=j
                        )
                        label_1 = "<b>All genes</b>"

                        fig.add_annotation(xref="x domain", yref="y domain", x=0.5, y=1.06, showarrow=False,
                                           text=label_1, font=dict(
                                family="Arial Black",
                                size=20
                            ), row=i, col=j)
                        # fig.update_xaxes(range=[0, max(data_1[column])])
                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(type="log", title_text="Count of genes",
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        fig.update_layout(showlegend=False)
                        fig.update_layout(height=fig_h, width=fig_w)
                        fig.update_xaxes(tickangle=45)


        return fig

    @staticmethod
    def plots(df):
        dataset_bin = pd.DataFrame()  # To contain our dataframe with our discretised continuous variables
        dataset_con = pd.DataFrame()
        if "ID" in df.columns:
            df = df.drop("ID",axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            dataset_bin['classical_label'] = df['classical_label']
            dataset_con['classical_label'] = df['classical_label']
            fig, axs = plt.subplots(nrows=len(df.columns) - 1, ncols=1, figsize=(12, 25))
        elif "core_label" in df.columns:
            dataset_bin['core_label'] = df['core_label']
            dataset_con['core_label'] = df['core_label']
            fig, axs = plt.subplots(nrows=len(df.columns) - 1, ncols=1, figsize=(12, 25))
        elif "Origin" in df.columns:
            dataset_bin['Origin'] = df['Origin']
            dataset_con['Origin'] = df['Origin']
            fig, axs = plt.subplots(nrows=len(df.columns) - 1, ncols=1, figsize=(12, 25))
        else:
            fig, axs = plt.subplots(nrows=len(df.columns), ncols=1, figsize=(9, 20))

        img = io.BytesIO()

        print(len(df.columns) - 1)
        for i, column in enumerate(df.columns, start=1):
            for j in range(1):
                if column.startswith("WGCNA") or column.startswith("classical_label") or column.startswith("core_label")\
                        or column.startswith("Origin"):
                    pass
                else:
                    dataset_con[column] = df[column]
                    plt.style.use('seaborn-whitegrid')
                    if "classical_label" in dataset_con.columns:
                        if(len(df.columns)==2):
                            axis_val = axs
                        else:
                            axis_val = axs[i-1]
                        data = pd.DataFrame()
                        data2 = pd.DataFrame()
                        data["1"] = dataset_con.loc[dataset_con['classical_label'] == 1.0][column].dropna()
                        data2["0"] = dataset_con.loc[dataset_con['classical_label'] == 0.0][column].dropna()
                        ttest, pval = ttest_ind(data["1"], data2["0"])
                        print("p-value", pval)
                        sns.distplot(dataset_con.loc[dataset_con['classical_label'] == 1.0][column].dropna(),hist = False, kde_kws={'linewidth': 3,"label": "classical"},
                                     ax=axis_val);
                        sns.distplot(dataset_con.loc[dataset_con['classical_label'] == 0.0][column].dropna(),hist = False, kde_kws={'linewidth': 3,"label": "other"},
                                     ax=axis_val);
                        if pval < 0.05:
                            print("we reject null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "There is significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E " % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "There is significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E  " % (column, pval)

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "There is no significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "There is no significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E" % (column, pval)
                        else:
                            p_val = "Atleast one of the groups contain Null values, therefore p-value is None"
                        axis_val.set_title(p_val + '\n Probability Density plots',fontsize=15,fontweight='bold')
                    elif "core_label" in dataset_con.columns:
                        if (len(df.columns) == 2):
                            axis_val = axs
                        else:
                            axis_val = axs[i - 1]
                        data1 = pd.DataFrame()
                        data2 = pd.DataFrame()
                        data3 = pd.DataFrame()
                        data4 = pd.DataFrame()
                        data1["1"] = dataset_con.loc[dataset_con['core_label'] == 1.0][column].dropna()
                        data2["2"] = dataset_con.loc[dataset_con['core_label'] == 2.0][column].dropna()
                        data3["3"] = dataset_con.loc[dataset_con['core_label'] == 3.0][column].dropna()
                        data4["4"] = dataset_con.loc[dataset_con['core_label'] == 4.0][column].dropna()
                        anovatest, pval = stats.f_oneway(data1, data2, data3, data4)
                        print("p-value", pval)
                        sns.distplot(dataset_con.loc[dataset_con['core_label'] == 1.0][column].dropna(), hist = False,kde_kws={'linewidth': 3,"label": "Core Gene"},
                                     ax=axis_val);
                        sns.distplot(dataset_con.loc[dataset_con['core_label'] == 2.0][column].dropna(),hist = False,kde_kws={'linewidth': 3,"label": "Near-Core Gene"},
                                     ax=axis_val);
                        sns.distplot(dataset_con.loc[dataset_con['core_label'] == 3.0][column].dropna(),hist = False, kde_kws={'linewidth': 3,"label": "Dispensable Gene"},
                                     ax=axis_val);
                        sns.distplot(dataset_con.loc[dataset_con['core_label'] == 4.0][column].dropna(),hist = False,kde_kws={'linewidth': 3,"label": "Private Gene"},
                                     ax=axis_val);

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "There is significant difference of %s in the different core genes" \
                                        " as p-value = %.2E  " % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "There is significant difference of %s in the different core genes" \
                                        " as p-value = %.2E  " % (column, pval)

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "There is no significant difference of %s in the different core genes" \
                                        " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "There is no significant difference of %s in the different core genes" \
                                        " as p-value = %.2E" % (column, pval)
                        else:
                            p_val = "Atleast one of the groups contain Null values, therefore p-value is None"
                        axis_val.set_title(p_val + '\n Probability Density plots',fontsize=15,fontweight='bold')
                    elif "Origin" in dataset_con.columns:
                        if (len(df.columns) == 2):
                            axis_val = axs
                        else:
                            axis_val = axs[i - 1]
                        data1 = pd.DataFrame()
                        data2 = pd.DataFrame()
                        data3 = pd.DataFrame()
                        data1["1"] = dataset_con.loc[dataset_con['Origin'] == 1.0][column].dropna()
                        data2["2"] = dataset_con.loc[dataset_con['Origin'] == 2.0][column].dropna()
                        data3["3"] = dataset_con.loc[dataset_con['Origin'] == 3.0][column].dropna()
                        anovatest, pval = stats.f_oneway(data1, data2, data3)
                        print("p-value", pval)
                        sns.distplot(dataset_con.loc[dataset_con['Origin'] == 1.0][column].dropna(), hist = False,kde_kws={'linewidth': 3,"label": "WGD Gene"},
                                     ax=axis_val);
                        sns.distplot(dataset_con.loc[dataset_con['Origin'] == 2.0][column].dropna(),hist = False,kde_kws={'linewidth': 3,"label": "Tandem Gene"},
                                     ax=axis_val);
                        sns.distplot(dataset_con.loc[dataset_con['Origin'] == 3.0][column].dropna(),hist = False, kde_kws={'linewidth': 3,"label": "Both Gene"},
                                     ax=axis_val);


                        if pval < 0.05:
                            print("we reject null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "There is significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E  " % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "There is significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E  " % (column, pval)

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "There is no significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "There is no significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E" % (column, pval)
                        else:
                            p_val = "Atleast one of the groups contain Null values, therefore p-value is None"

                        axis_val.set_title(p_val + '\n Probability Density plots',fontsize=15,fontweight='bold')

                    else:
                        if(len(df.columns)==2):
                            axis_val = axs
                        else:
                            print(i-1)
                            axis_val = axs[i-1]

                        sns.distplot(dataset_con[column].dropna(),hist = False,kde_kws={'linewidth': 3},
                                     ax=axis_val);


                        axis_val.set_title('Probability Density plots',fontsize=15,fontweight='bold')

                    axis_val.legend(loc=0, ncol=1, bbox_to_anchor=(0, 0, 1, 1),
                                    prop={'size': 15,'weight':'bold'}, fancybox=True, shadow=False, title='Gene Type',
                                    title_fontsize=15).get_title().set_weight("bold")
                    axis_val.set_ylabel('Density',fontsize=15,fontweight='bold')
                    axis_val.xaxis.set_tick_params(color='b',width=4)
                    axis_val.yaxis.set_tick_params(color='b',width=4)
                    axis_val.xaxis.get_label().set_fontsize(15)
                    axis_val.xaxis.get_label().set_fontweight('bold')
                    axis_val.set_xticklabels(axis_val.get_xticks(),fontdict={'fontweight':'bold','fontsize':13})
                    axis_val.set_yticklabels(axis_val.get_yticks(),fontdict={'fontweight':'bold','fontsize':13})
                    axis_val.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    axis_val.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))
                    plt.tight_layout()

        fig.savefig('C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/plot.png')
        fig.savefig(img, format='png')
        img.seek(0)
        plot_url = base64.b64encode(img.getvalue()).decode()
        return plot_url

    @staticmethod
    def pair_plots(df):
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            var = [x for i,x in enumerate(df.columns) if x!='classical_label']
            fig = px.scatter_matrix(df,
                                    dimensions=var,
                                    color="classical_label", symbol="classical_label",
                                    opacity=0.5,
                                    title="<b>Scatter matrix of Classical/Other genes</b>")

        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1,'Core Genes')
            df['core_label'] = df['core_label'].replace(2,'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3,'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4,'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')
            var = [x for i, x in enumerate(df.columns) if x != 'core_label']
            fig = px.scatter_matrix(df,
                                    dimensions=var,
                                    color="core_label",symbol="core_label",
                                    opacity=0.5,
                                    title="<b>Scatter matrix of Core/Non-core/Dispensible/Private Genes</b>")
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1,'WGD Genes')
            df['Origin'] = df['Origin'].replace(2,'Tandem Genes')
            df['Origin'] = df['Origin'].replace(3,'Both Genes')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')
            var = [x for i, x in enumerate(df.columns) if x != 'Origin']
            fig = px.scatter_matrix(df,
                                    dimensions=var,
                                    color="Origin",symbol="Origin",
                                    opacity=0.5,
                                    title="<b>Scatter matrix of WGD/Tandem/Both Genes</b>")
        else:
            var = [x for i, x in enumerate(df.columns)]
            fig = px.scatter_matrix(df,
                                    dimensions=var,
                                    opacity=0.5,
                                    line={'color': 'black'},
                                    title="<b>Scatter matrix of All Genes</b>")
        fig.update_layout(
            dragmode='select',
            width=2300,
            height=2300,
            hovermode='closest',legend=dict(title_font_family="Arial Black",title_font_size=22,title_font_color="black",
                font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    @staticmethod
    def box_plots(df):
        print(df.shape)
        dataset_bin = pd.DataFrame()  # To contain our dataframe with our discretised continuous variables
        fig_h = 0
        fig_w = 0
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            dataset_bin['classical_label'] = df['classical_label']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 1550
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 1550
            else:
                fig_h = 3000
                fig_w = 1550
        elif "core_label" in df.columns:
            dataset_bin['core_label'] = df['core_label']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 4000
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 4000
            else:
                fig_h = 3200
                fig_w = 4000
        elif "Origin" in df.columns:
            dataset_bin['Origin'] = df['Origin']
            if len(df.columns) == 2:
                fig_h = 1200
                fig_w = 2700
            elif len(df.columns) == 3:
                fig_h = 2000
                fig_w = 2700
            else:
                fig_h = 3000
                fig_w = 2700
        else:
            if len(df.columns) == 2:
                fig_h = 1000
                fig_w = 800
            elif len(df.columns) == 3:
                fig_h = 1800
                fig_w = 800
            else:
                fig_h = 2800
                fig_w = 800
        if "classical_label" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=2)
        elif "core_label" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=4)
        elif "Origin" in df.columns:
            fig = make_subplots(rows=len(df.columns), cols=3)
        else:
            fig = make_subplots(rows=len(df.columns), cols=1)
        for i, column in enumerate(df.columns, start=1):
            print(i)
            for j in range(1, 2):
                print(j)
                print([i], [j])
                if column.startswith("WGCNA") or column.startswith("classical_label") or column.startswith(
                        "core_label") or column.startswith("Origin"):
                    pass
                    print(column)
                else:
                    if "classical_label" in dataset_bin.columns:
                        data_1 = df.loc[df['classical_label'] == 1.0][column].dropna()
                        data_0 = df.loc[df['classical_label'] == 0.0][column].dropna()
                        ttest, pval = ttest_ind(data_1, data_0)
                        print("p-value", pval)
                        hover_label = 'Gene type'
                        fig.add_trace(
                            go.Box(y=data_1,
                                         name='Classical genes',
                                            marker=dict(color='#1f77b4'),
                                         # hovertemplate=column + ': %{x} <br>: %{y}'
                                         ),
                            row=i, col=j
                        )

                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(title_text=column,showline=True, linewidth=2, linecolor='black',
                        title_font=dict(size=20, family='Arial Black', color='black'),
                        tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)

                        fig.add_trace(
                            go.Box(y=data_0,
                                   name='Other genes',
                                   marker=dict(color='#ff7f0e'),
                                   # hovertemplate=column + ': %{x} <br>: %{y}'
                                   ),
                            row=i, col=j
                        )

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "<b>There is significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E  </b>" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "<b>There is significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E  </b>" % (column, pval)
                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 2.04
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 2.08
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 2.18
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.24

                            elif (int(len(p_val)) > 128):

                                x_len = 2.4
                            else:

                                x_len = 2.0
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            # fig.update_layout(height=fig_h, width=fig_w)

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "There is no significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "There is no significant difference of %s in both Classical and Other genes" \
                                        " as p-value = %.2E" % (column, pval)

                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 2.04
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 2.08
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 2.18
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.24
                            elif int(len(p_val)) > 128:

                                x_len = 2.4
                            else:

                                x_len = 2.0
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)

                        else:
                            print("we reject null hypothesis")

                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None</b>"

                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 2.04
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 2.08
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 2.18
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.24
                            elif (int(len(p_val)) > 128):

                                x_len = 2.4
                            else:

                                x_len = 2.0
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                        fig.update_layout( title='Box Plots',
                                           height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                                        font=dict(size=18, family='Arial Black', color='black')))

                    elif "core_label" in dataset_bin.columns:
                        data_1 = df.loc[df['core_label'] == 1.0][column].dropna()
                        data_2 = df.loc[df['core_label'] == 2.0][column].dropna()
                        data_3 = df.loc[df['core_label'] == 3.0][column].dropna()
                        data_4 = df.loc[df['core_label'] == 4.0][column].dropna()
                        hover_label = 'Gene type'
                        anovatest, pval = stats.f_oneway(data_1, data_2, data_3, data_4)
                        print("p-value", pval)
                        fig.add_trace(
                            go.Box(y=data_1,
                                   name='Core genes',
                                   marker=dict(color='#1f77b4'),
                                   ),
                            row=i, col=j
                        )

                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(title_text=column,
                                         showline=True, linewidth=2, linecolor='black',
                                         title_font=dict(size=20, family='Arial Black', color='black'),
                                         tickfont=dict(size=18, family='Arial Black', color='black'),
                                         row=i, col=j)

                        fig.add_trace(
                            go.Box(y=data_2,
                                   name='Non-Core genes',
                                   marker=dict(color='#ff7f0e'),
                                   ),
                            row=i, col=j
                        )

                        fig.add_trace(
                            go.Box(y=data_3,
                                   name='Dispensable Genes',
                                   marker=dict(color='#2ca02c'),
                                   ),
                            row=i, col=j
                        )

                        fig.add_trace(
                            go.Box(y=data_4,
                                   name='Private genes',
                                   marker=dict(color='#d62728'),
                                   ),
                            row=i, col=j
                        )

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "<b>There is significant difference of %s in the different core genes" \
                                        " as p-value = %.2E  </b>" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "<b>There is significant difference of %s in the different core genes" \
                                        " as p-value = %.2E  </b>" % (column, pval)
                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 1.4
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 1.8
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 2.0
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.14
                            elif (int(len(p_val)) > 128):

                                x_len = 2.24
                            else:

                                x_len = 1.3
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            # fig.update_layout(height=fig_h, width=fig_w)

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "<b>There is no significant difference of %s in the different core genes" \
                                        " as p-value = %.2E  </b>" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "<b>There is no significant difference of %s in the different core genes" \
                                        " as p-value = %.2E  </b>" % (column, pval)

                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 1.4
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 1.8
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 2.0
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.14
                            elif int(len(p_val)) > 128:

                                x_len = 2.24
                            else:

                                x_len = 1.3
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                        else:
                            print("we reject null hypothesis")
                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None </b>"
                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 1.4
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 1.8
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 2.0
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.14
                            elif (int(len(p_val)) > 128):

                                x_len = 2.24
                            else:

                                x_len = 1.3
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)
                        fig.update_layout(title='<b>Box Plots</b>',
                                          height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                                            font=dict(size=18, family='Arial Black', color='black')))
                    elif "Origin" in dataset_bin.columns:
                        data_1 = df.loc[df['Origin'] == 1.0][column].dropna()
                        data_2 = df.loc[df['Origin'] == 2.0][column].dropna()
                        data_3 = df.loc[df['Origin'] == 3.0][column].dropna()
                        hover_label = 'Gene type'
                        anovatest, pval = stats.f_oneway(data_1, data_2, data_3)
                        print("p-value", pval)
                        fig.add_trace(
                            go.Box(y=data_1,
                                   name='WGD genes',
                                   marker=dict(color='#1f77b4'),
                                   ),
                            row=i, col=j
                        )

                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(title_text=column,showline=True, linewidth=2, linecolor='black',
                        title_font=dict(size=20, family='Arial Black', color='black'),
                        tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)

                        fig.add_trace(
                            go.Box(y=data_2,
                                   name='Tandem genes',
                                   marker=dict(color='#ff7f0e'),
                                   ),
                            row=i, col=j
                        )

                        fig.add_trace(
                            go.Box(y=data_3,
                                   name='Both Genes',
                                   marker=dict(color='#2ca02c'),
                                   ),
                            row=i, col=j
                        )

                        if pval < 0.05:
                            print("we reject null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "<b>There is significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E  </b>" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "<b>There is significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E  </b>" % (column, pval)
                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 1.6
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 1.7
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 1.9
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.0
                            elif (int(len(p_val)) > 128):

                                x_len = 2.1
                            else:

                                x_len = 1.5
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="LightSeaGreen"
                                ), row=i, col=j)
                            # fig.update_layout(height=fig_h, width=fig_w)

                        elif pval > 0.05:
                            print("we accept null hypothesis")
                            if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                                p_val = "<b>There is no significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E  </b>" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                            else:
                                p_val = "<b>There is no significant difference of %s in the different Origin genes" \
                                        " as p-value = %.2E  </b>" % (column, pval)

                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 1.6
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 1.7
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 1.9
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.0
                            elif int(len(p_val)) > 128:

                                x_len = 2.1
                            else:

                                x_len = 1.5
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="Crimson"
                                ), row=i, col=j)
                        else:

                            p_val = "<b>Atleast one of the groups contain Null values, therefore p-value is None</b>"
                            print("p_val len", len(p_val), type(len(p_val)))
                            if (len(p_val) > 110) and (len(p_val) < 115):

                                x_len = 1.6
                            elif (len(p_val) >= 115) and (len(p_val) < 120):

                                x_len = 1.7
                            elif (len(p_val) >= 120) and (len(p_val) < 125):

                                x_len = 1.9
                            elif (len(p_val) >= 125) and (len(p_val) <= 128):

                                x_len = 2.0
                            elif (int(len(p_val)) > 128):

                                x_len = 2.1
                            else:

                                x_len = 1.5
                            fig.add_annotation(xref="x domain", yref="y domain", x=x_len, y=1.09, showarrow=False,
                                               text=p_val, font=dict(
                                    family="Arial Black",
                                    size=19,
                                    color="purple"
                                ), row=i, col=j)

                        fig.update_layout(title='<b>Box Plots</b>',
                                          height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                                            font=dict(size=18, family='Arial Black', color='black')))
                    else:
                        data_1 = df[column].dropna()
                        hover_label = 'Gene type'
                        fig.add_trace(
                            go.Box(y=data_1,
                                            marker=dict(color='#1f77b4'),
                                            name='All Genes',
                                         # hovertemplate=column + ': %{x} <br>: %{y}'
                                         ),
                            row=i, col=j
                        )

                        fig.update_xaxes(title_text=hover_label,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_yaxes(title_text=column,title_font=dict(size=20, family='Arial Black', color='black'),
                                         showline=True, linewidth=2, linecolor='black',
                                         tickfont=dict(size=18, family='Arial Black', color='black'), row=i, col=j)
                        fig.update_layout( title='Box Plots',
                                           height=fig_h, width=fig_w,legend=dict(title_font_family="Arial Black",
                                            font=dict(size=18, family='Arial Black', color='black')))
                        fig.update_layout(showlegend=False)
                fig.update_layout(title_font=dict(size=20, family='Arial Black', color='black'))
        return fig

    @staticmethod
    def violin_plots(df):
        fig_w = 0
        fig_h = 0
        dataset_con = pd.DataFrame()
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            dataset_con['classical_label'] = df['classical_label']
            a = len(df.columns) - 1  # number of rows
            b = 1  # number of columns
            c = 1  # initialize plot counter
            if len(df.columns) == 2:
                fig_h = 12
                fig_w = 10
            elif len(df.columns) == 3:
                fig_h = 12
                fig_w = 19
            else:
                fig_h = 12
                fig_w = 21
        elif "core_label" in df.columns:
            dataset_con['core_label'] = df['core_label']
            a = len(df.columns) - 1  # number of rows
            b = 1  # number of columns
            c = 1  # initialize plot counter
            if len(df.columns) == 2:
                fig_h = 12
                fig_w = 10
            elif len(df.columns) == 3:
                fig_h = 12
                fig_w = 19
            else:
                fig_h = 12
                fig_w = 21
        elif "Origin" in df.columns:
            dataset_con['Origin'] = df['Origin']
            a = len(df.columns) - 1  # number of rows
            b = 1  # number of columns
            c = 1  # initialize plot counter
            if len(df.columns) == 2:
                fig_h = 12
                fig_w = 10
            elif len(df.columns) == 3:
                fig_h = 12
                fig_w = 19
            else:
                fig_h = 12
                fig_w = 20

        else:
            a = len(df.columns) # number of rows
            b = 1  # number of columns
            c = 1  # initialize plot counter
            if len(df.columns) == 2:
                fig_h = 7
                fig_w = 8
            elif len(df.columns) == 3:
                fig_h = 7
                fig_w = 17
            else:
                fig_h = 7
                fig_w = 19
        img = io.BytesIO()
        #  plot Numerical Data

        plt.style.use('seaborn-whitegrid')
        # fig = plt.figure(figsize=(12, 17))
        print(fig_h, fig_w)
        fig = plt.figure(figsize=(fig_h, fig_w))
        for i, column in enumerate(df.columns):
            if column.startswith("classical_label") or column.startswith("core_label") or column.startswith("Origin") :
                pass
            else:
                dataset_con[column] = df[column]
                plt.style.use('seaborn-whitegrid')
                if "classical_label" in dataset_con.columns:
                    data = pd.DataFrame()
                    data2 = pd.DataFrame()
                    data["1"] = dataset_con.loc[dataset_con['classical_label'] == 1.0][column].dropna()
                    data2["0"] = dataset_con.loc[dataset_con['classical_label'] == 0.0][column].dropna()
                    ttest, pval = ttest_ind(data["1"], data2["0"])
                    print("p-value", pval)
                    ax = fig.add_subplot(a, b, c)
                    sns.violinplot(x="classical_label", y=column, data=dataset_con);
                    plt.xlabel("Gene type")
                    plt.ylabel(column)
                    ax.xaxis.get_label().set_fontsize(15)
                    ax.yaxis.get_label().set_fontsize(15)
                    c = c + 1
                    if pval < 0.05:
                        print("we reject null hypothesis")
                        if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                            p_val = "There is significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                        else:
                            p_val = "There is significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E" % (column, pval)

                    elif pval > 0.05:
                        print("we accept null hypothesis")
                        if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                            p_val = "There is no significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                        else:
                            p_val = "There is no significant difference of %s in both Classical and Other genes" \
                                    " as p-value = %.2E" % (column, pval)
                    else:
                        print("we reject null hypothesis")

                        p_val = "Atleast one of the groups contain Null values, therefore p-value is None"
                    ax.set_title(p_val, fontsize=15,fontweight='bold')
                    ax.xaxis.set_tick_params(color='b', width=4)
                    ax.yaxis.set_tick_params(color='b', width=4)
                    ax.xaxis.get_label().set_fontweight('bold')
                    ax.yaxis.get_label().set_fontweight('bold')
                    ax.set_xticklabels(ax.get_xticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.set_yticklabels(ax.get_yticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(ticks=[0, 1], labels=['Classical', 'others'])
                elif "core_label" in dataset_con.columns:
                    data1 = pd.DataFrame()
                    data2 = pd.DataFrame()
                    data3 = pd.DataFrame()
                    data4 = pd.DataFrame()
                    data1["1"] = dataset_con.loc[dataset_con['core_label'] == 1.0][column].dropna()
                    data2["2"] = dataset_con.loc[dataset_con['core_label'] == 2.0][column].dropna()
                    data3["3"] = dataset_con.loc[dataset_con['core_label'] == 3.0][column].dropna()
                    data4["4"] = dataset_con.loc[dataset_con['core_label'] == 4.0][column].dropna()
                    anovatest, pval = stats.f_oneway(data1, data2, data3, data4)
                    print("p-value", pval)
                    ax = fig.add_subplot(a, b, c)
                    sns.violinplot(x="core_label", y=column, data=dataset_con);
                    plt.xlabel("Gene type")
                    plt.ylabel(column)
                    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
                    # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                    ax.xaxis.get_label().set_fontsize(15)
                    ax.yaxis.get_label().set_fontsize(15)
                    c = c + 1
                    if pval < 0.05:
                        print("we reject null hypothesis")
                        if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                            p_val = "There is significant difference of %s in the different core genes" \
                                    " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                        else:
                            p_val = "There is significant difference of %s in the different core genes" \
                                    " as p-value = %.2E" % (column, pval)

                    elif pval > 0.05:
                        print("we accept null hypothesis")
                        if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                            p_val = "There is no significant difference of %s in the different core genes" \
                                    " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                        else:
                            p_val = "There is no significant difference of %s in the different core genes" \
                                    " as p-value = %.2E" % (column, pval)
                    else:
                        print("we reject null hypothesis")

                        p_val = "Atleast one of the groups contain Null values, therefore p-value is None"
                    ax.set_title(p_val, fontsize=15, fontweight='bold')
                    ax.xaxis.set_tick_params(color='b', width=4)
                    ax.yaxis.set_tick_params(color='b', width=4)
                    ax.xaxis.get_label().set_fontweight('bold')
                    ax.yaxis.get_label().set_fontweight('bold')
                    ax.set_xticklabels(ax.get_xticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.set_yticklabels(ax.get_yticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(ticks=[0, 1, 2, 3],
                              labels=['Core Gene', 'Near-Core Gene', 'Dispensable Gene', 'Private Gene'])
                elif "Origin" in dataset_con.columns:
                    data1 = pd.DataFrame()
                    data2 = pd.DataFrame()
                    data3 = pd.DataFrame()
                    data1["1"] = dataset_con.loc[dataset_con['Origin'] == 1.0][column].dropna()
                    data2["2"] = dataset_con.loc[dataset_con['Origin'] == 2.0][column].dropna()
                    data3["3"] = dataset_con.loc[dataset_con['Origin'] == 3.0][column].dropna()
                    anovatest, pval = stats.f_oneway(data1, data2, data3)
                    print("p-value", pval)
                    ax = fig.add_subplot(a, b, c)
                    sns.violinplot(x="Origin", y=column, data=dataset_con);
                    plt.xlabel("Gene type")
                    plt.ylabel(column)
                    ax.xaxis.get_label().set_fontsize(15)
                    ax.yaxis.get_label().set_fontsize(15)
                    c = c + 1
                    if pval < 0.05:
                        print("we reject null hypothesis")
                        if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                            p_val = "There is significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                        else:
                            p_val = "There is significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E" % (column, pval)

                    elif pval > 0.05:
                        print("we accept null hypothesis")
                        if column == 'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':
                            p_val = "There is no significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E" % ('var_eff_5_UTR_premature_start_codon_gain', pval)
                        else:
                            p_val = "There is no significant difference of %s in the different Origin genes" \
                                    " as p-value = %.2E" % (column, pval)
                    else:
                        print("we reject null hypothesis")

                        p_val = "Atleast one of the groups contain Null values, therefore p-value is None"
                    ax.set_title(p_val, fontsize=15,fontweight='bold')
                    ax.xaxis.set_tick_params(color='b', width=4)
                    ax.yaxis.set_tick_params(color='b', width=4)
                    ax.xaxis.get_label().set_fontweight('bold')
                    ax.yaxis.get_label().set_fontweight('bold')
                    ax.set_xticklabels(ax.get_xticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.set_yticklabels(ax.get_yticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(ticks=[0, 1, 2], labels=['WGD Gene','Tandem Gene','Both Gene'])
                else:
                    ax = fig.add_subplot(a, b, c)
                    sns.violinplot(y=column, data=dataset_con);
                    plt.xlabel("All genes")
                    plt.ylabel(column)
                    # plt.xticks(ticks=[0], labels=['All genes'])
                    # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

                    ax.xaxis.get_label().set_fontsize(15)
                    ax.yaxis.get_label().set_fontsize(15)
                    ax.xaxis.set_tick_params(color='b', width=4)
                    ax.yaxis.set_tick_params(color='b', width=4)
                    ax.xaxis.get_label().set_fontweight('bold')
                    ax.yaxis.get_label().set_fontweight('bold')
                    ax.set_xticklabels(ax.get_xticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.set_yticklabels(ax.get_yticks(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    c = c + 1

        fig.suptitle('Violin plots', fontsize=17,fontweight='bold')
        fig.subplots_adjust(top=1.1)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.savefig('C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/plot.png')
        fig.savefig(img, format='png')
        img.seek(0)
        plot_url = base64.b64encode(img.getvalue()).decode()
        return plot_url

    @staticmethod
    def joint_plots(df):
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        plt.style.use('seaborn-whitegrid')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            var = [x for i, x in enumerate(df.columns) if x != 'classical_label']
            fig = px.scatter(df, x=var[0], y=var[1], marginal_x="histogram", marginal_y="rug",
                             color="classical_label", title="<b>Joint Plots</b><br>Click on the legend items!"
                             )

        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')
            var = [x for i, x in enumerate(df.columns) if x != 'core_label']
            fig = px.scatter(df, x=var[0], y=var[1], marginal_x="histogram", marginal_y="rug",
                             color="core_label", title="<b>Joint Plots</b><br>Click on the legend items!"
                             )

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')
            var = [x for i, x in enumerate(df.columns) if x != 'Origin']
            fig = px.scatter(df, x=var[0], y=var[1], marginal_x="histogram", marginal_y="rug",
                             color="Origin", title="<b>Joint Plots</b>"
                             )
        else:
            var = [x for i, x in enumerate(df.columns)]
            fig = px.scatter(df, x=var[0], y=var[1], marginal_x="histogram", marginal_y="rug",
                        title="<b>Joint Plots</b><br>Click on the legend items!"
                             )
        # fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
        # fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    @staticmethod
    def scatter_plots(df):
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        df = df.dropna()
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant': 'var_eff_5_UTR_SCG'},
                      inplace=True)
        if "core_label" in df.columns:
            # df = df.drop("core_label", axis=1)
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')
            var = [x for i, x in enumerate(df.columns) if x != 'core_label']
            fig = px.scatter(x=df[var[0]], y=df[var[1]], trendline="ols", color=df["core_label"],
                             title="<b>Scatter plots </b>"
                             )
        elif "Origin" in df.columns:
            # df = df.drop("Origin", axis=1)
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')
            var = [x for i, x in enumerate(df.columns) if x != 'Origin']
            fig = px.scatter(x=df[var[0]], y=df[var[1]], trendline="ols", color=df["Origin"],
                             title="<b>Scatter plots </b>"
                             )
        elif "classical_label" in df.columns:
            # df = df.drop("classical_label", axis=1)
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            var = [x for i, x in enumerate(df.columns) if x != 'classical_label']
            fig = px.scatter(x=df[var[0]], y=df[var[1]], trendline="ols", color=df["classical_label"],
                             title="<b>Scatter plots </b>"
                             )
        else:
            var = [x for i, x in enumerate(df.columns)]
            fig = px.scatter(x=df[var[0]], y=df[var[1]], trendline="ols",
                             title="<b>Scatter plots </b>"
                             )
        print(df.columns)

        print(len(df))

        fig.update_xaxes(title_text=var[0])
        fig.update_yaxes(title_text=var[1])
        # fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
        # fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )

        return fig


    @staticmethod
    def correlation_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        img = io.BytesIO()
        plt.style.use('seaborn-whitegrid')
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.labelweight"] = "bold"
        if "classical_label" in df.columns:
            # correlation - classical
            df_classical = df[(df["classical_label"] == 1)]
            df_classicalCorr = df_classical.drop(["classical_label"], axis=1).corr()
            # correlation - other
            df_other = df[(df["classical_label"] == 0)]
            df_otherCorr = df_other.drop(["classical_label"], axis=1).corr()
            fig = plt.figure(figsize=(15, 8))
            plt.subplot(121)  # subplot 1 - classical
            plt.title('Classical.corr(),  subplot: 121', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_classicalCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

            plt.subplot(122)  # subplot 2 - Other
            plt.title('Other.corr(),  subplot: 122', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_otherCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

        elif "core_label" in df.columns:
        # correlation - core
            df_core = df[(df["core_label"] == 1)]
            df_coreCorr = df_core.drop(["core_label"], axis=1).corr()
            # correlation - near-core
            df_near = df[(df["core_label"] == 2)]
            df_nearCorr = df_near.drop(["core_label"], axis=1).corr()
            # correlation - dispensible
            df_dispensible = df[(df["core_label"] == 3)]
            df_dispensibleCorr = df_dispensible.drop(["core_label"], axis=1).corr()
            # correlation - private
            df_private = df[(df["core_label"] == 4)]
            df_privateCorr = df_private.drop(["core_label"], axis=1).corr()
            fig = plt.figure(figsize=(16, 12))
            plt.subplot(221)  # subplot 1 - classical
            plt.title('Core.corr(),  subplot: 121', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_coreCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

            plt.subplot(222)  # subplot 2 - Other
            plt.title('Near-Core.corr(),  subplot: 122', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_nearCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

            plt.subplot(223)  # subplot 1 - classical
            plt.title('Dispensible.corr(),  subplot: 221', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_dispensibleCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

            plt.subplot(224)  # subplot 2 - Other
            plt.title('Private.corr(),  subplot: 222', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_privateCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

        elif "Origin" in df.columns:
        # correlation - core
            df_core = df[(df["Origin"] == 1)]
            df_coreCorr = df_core.drop(["Origin"], axis=1).corr()
            # correlation - near-core
            df_near = df[(df["Origin"] == 2)]
            df_nearCorr = df_near.drop(["Origin"], axis=1).corr()
            # correlation - dispensible
            df_dispensible = df[(df["Origin"] == 3)]
            df_dispensibleCorr = df_dispensible.drop(["Origin"], axis=1).corr()

            fig = plt.figure(figsize=(15, 12))
            plt.subplot(221)  # subplot 1 - classical
            plt.title('WGD.corr(),  subplot: 121', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_coreCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

            plt.subplot(222)  # subplot 2 - Other
            plt.title('Tandem.corr(),  subplot: 122', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_nearCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

            plt.subplot(223)  # subplot 1 - classical
            plt.title('Both.corr(),  subplot: 221', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_dispensibleCorr, annot=True, fmt='.2f', square=True, cmap='Blues_r')
        else:
            # correlation - classical

            df_Corr = df.corr()
            fig = plt.figure(figsize=(15, 8))
            plt.subplot(121)  # subplot 1 - classical
            plt.title('All Genes.corr(),  subplot: 121', fontname="Arial Black", size=15,fontweight="bold")
            sns.heatmap(df_Corr, annot=True, fmt='.2f', square=True, cmap='Blues_r')

        plt.suptitle('Correlation plots', fontname="Arial Black", size=20,fontweight="bold")
        plt.tight_layout()
        plt.savefig('C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/plot.png')
        plt.savefig(img, format='png')
        img.seek(0)
        plot_url = base64.b64encode(img.getvalue()).decode()
        return plot_url

    @staticmethod
    def downsampled_Hierarchial_heatmap(df):
        df=df[df.columns].set_index('ID')
        # label_ID = list(df['ID'])
        df = df.fillna(0)
        # if "ID" in df.columns:
        #     df = df.drop("ID", axis=1)
        if "classical_label" in df.columns:
            df_label= df.drop(["classical_label"], axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        elif "core_label" in df.columns:
            df_label = df.drop(["core_label"], axis=1)
        elif "Origin" in df.columns:
            df_label = df.drop(["Origin"], axis=1)
        cols = df_label.columns
        # df_label = normalize(df_label)
        # df_label = pd.DataFrame(df_label, columns=cols)
        std_scaler = StandardScaler()
        df_label = std_scaler.fit_transform(df_label)
        df_label = pd.DataFrame(df_label, columns=cols)
        print(df_label.head())
        columns = list(cols)
        print(columns)
        rows = list(df.index)
        print(rows)
        print(df_label.values.shape)
        color_palette = [
            'rgb(128, 0, 96)',
            'rgb(230, 115, 0)',
            'rgb(255, 191, 0)'
        ]
        component = dashbio.Clustergram(
            data=df_label.values,
            color_threshold={'row': 150,
                             'col': 700},
            color_map=color_palette,
            color_list={
                'row': [color_palette[0], color_palette[1], color_palette[2]],
                'col': [color_palette[1], color_palette[2], color_palette[0]],
                'bg': 'rgb(255,255,255)'
            },
            line_width=2,
            column_labels=list(df_label.columns.values),
            row_labels=list(df.index),
            hidden_labels=['row'],
            width=1300,
            height=800
        )
        component.update_layout(title="<b>Hierarchial heatmap of the downsampled genes</b>", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        # component.show()
        return component

    @staticmethod
    def downsampled_dendogram_plot(df):
        downsampled_cluster_plot(df)
        rid = list(df['ID'])
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        cols = df.columns
        df = df.fillna(0)
        print(cols)
        std_scaler = StandardScaler()
        df_label = std_scaler.fit_transform(df)
        df_label = pd.DataFrame(df_label, columns=cols)
        fig = ff.create_dendrogram(df_label, labels=rid,
                                   linkagefun=lambda x: linkage(df_label, 'ward', metric='euclidean'))
        fig.update_layout(width=1300, height=730, title="<b>Dendogram of the downsampled genes</b>",font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        return fig

    @staticmethod
    def downsampled_hierarchial_scatter_plot(df,cluster):
        K=int(cluster)
        rid = list(df['ID'])
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        # df = df.iloc[: , :5]
        cols = df.columns

        print(cols)
        df = df.fillna(0)
        # df_label = normalize(df)
        std_scaler = StandardScaler()
        df_label = std_scaler.fit_transform(df)
        df_label = pd.DataFrame(df_label, columns=cols)
        # print(df_label[df_label.classical_label !=0].head())
        X = df_label.iloc[:, :].values
        fig = ff.create_dendrogram(df_label, labels=rid,
                                   linkagefun=lambda x: linkage(df_label, 'ward', metric='euclidean'))
        # for i in fig.data:
        #     K.append(i.to_plotly_json()['marker']['color'])
        # A_set = set(K)
        # print(A_set)
        # print(len(A_set))
        print(K)
        cluster = AgglomerativeClustering(n_clusters=K, affinity='euclidean', linkage='ward')
        cluster.fit_predict(X)
        labels = cluster.labels_
        labels = labels.tolist()
        labels = [str(x) for x in labels]
        print(labels)
        df_label = df_label[cols[:5]]
        print(df_label.head())
        print(cols[:5])
        fig = px.scatter_matrix(df_label,
                                dimensions=cols[:5],
                                color=labels,
                                opacity = 0.5)
        fig.update_layout(width=2300, height=2300, title="<b>Hierarchial scatter plot of the downsampled genes</b>\n(shows a maximum of 5 variable at a time for better vizualization)" , font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        # fig.show()
        return fig



    @staticmethod
    def Count_downsampled_Hierarchial_heatmap(df):
        df = df[df.columns].set_index('ID')
        # label_ID = list(df['ID'])
        df = df.fillna(0)
        # if "ID" in df.columns:
        #     df = df.drop("ID", axis=1)
        if "classical_label" in df.columns:
            df_label = df.drop(["classical_label"], axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        elif "core_label" in df.columns:
            df_label = df.drop(["core_label"], axis=1)
        elif "Origin" in df.columns:
            df_label = df.drop(["Origin"], axis=1)
        cols = df_label.columns
        # df_label = normalize(df_label)
        # df_label = normalize(df_label)
        # df_label = pd.DataFrame(df_label, columns=cols)
        print(df_label.head())
        columns = list(cols)
        print(columns)
        rows = list(df.index)
        print(rows)
        df_label = df_label.T
        print(df_label.values.shape)
        print(df_label.values.max())

        color_palette = [
            'rgb(128, 0, 96)',
            'rgb(230, 115, 0)',
            'rgb(255, 191, 0)'
        ]
        component = dashbio.Clustergramm(
            data=df_label.values,
            column_labels=list(df_label.columns.values),
            row_labels=columns,
            cluster='all',
            display_ratio=[0.3, 0.1],
            hidden_labels=['col'],
            color_threshold=dict(row=150, col=700),
            color_map=color_palette,
            color_list={
                'row': [color_palette[0], color_palette[1], color_palette[2]],
                'col': [color_palette[1], color_palette[2], color_palette[0]],
                'bg': 'rgb(255,255,255)'
            },
            annotation_font=dict(
                color='white',
                size=10
            ),

            tick_font=dict(
                size=11,
                color='black'
            ),
            optimal_leaf_order=True,
            center_values=False,
            log_transform=False,
            width=1400,
            height=800
        )
        # data = [component]
        # print(component)
        # component.data.update(marker=dict(cmin=0, cmax=50))
        component.update_layout(width=1500,title="<b>Hierarchial heatmap of the downsampled genes</b>", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        return component

    @staticmethod
    def varionomic_downsampled_Hierarchial_heatmap(df):
        df = df[df.columns].set_index('ID')
        # label_ID = list(df['ID'])
        df = df.fillna(0)
        # if "ID" in df.columns:
        #     df = df.drop("ID", axis=1)
        if "classical_label" in df.columns:
            df_label = df.drop(["classical_label"], axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant': 'var_eff_5_UTR_SCG'},
                      inplace=True)
        elif "core_label" in df.columns:
            df_label = df.drop(["core_label"], axis=1)
        elif "Origin" in df.columns:
            df_label = df.drop(["Origin"], axis=1)
        cols = df_label.columns
        # df_label = normalize(df_label)
        # df_label = pd.DataFrame(df_label, columns=cols)
        std_scaler = StandardScaler()
        df_label = std_scaler.fit_transform(df_label)
        df_label = pd.DataFrame(df_label, columns=cols)
        print(df_label.head())
        columns = list(cols)
        print(columns)
        rows = list(df.index)
        print(rows)
        print(df_label.values.shape)
        color_palette = [
            'rgb(128, 0, 96)',
            'rgb(230, 115, 0)',
            'rgb(255, 191, 0)'
        ]
        component = dashbio.Clustergram_var(
            data=df_label.values,
            color_threshold={'row': 150,
                             'col': 700},
            column_labels=list(df_label.columns.values),
            color_map=color_palette,
            color_list={
                'row': [color_palette[0], color_palette[1], color_palette[2]],
                'col': [color_palette[1], color_palette[2], color_palette[0]],
                'bg': 'rgb(255,255,255)'
            },
            line_width=2,
            row_labels=list(df.index),
            hidden_labels=['row'],
            width=1500,
            height=800
        )
        component.update_layout(title="<b>Hierarchial heatmap of the downsampled genes</b>", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        # component.show()
        return component

    @staticmethod
    def Count_downsampled_dendogram_plot(df):
        count_downsampled_cluster_plot(df)
        rid = list(df['ID'])
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df_label = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        cols = df.columns
        df = df.fillna(0)
        print(cols)
        df_label = df
        fig = ff.create_dendrogram(df_label, labels=rid,
                                   linkagefun=lambda x: linkage(df_label, 'ward', metric='euclidean'))
        fig.update_layout(width=1300, height=730, title="<b>Dendogram of the downsampled genes</b>", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        return fig

    @staticmethod
    def Count_downsampled_hierarchial_scatter_plot(df,cluster):
        K = int(cluster)
        rid = list(df['ID'])
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        cols = df.columns
        print(cols)
        df = df.fillna(0)
        df_label = df
        # print(df_label[df_label.classical_label !=0].head())
        X = df_label.iloc[:, :].values
        fig = ff.create_dendrogram(df_label, labels=rid,
                                   linkagefun=lambda x: linkage(df_label, 'ward', metric='euclidean'))
        # for i in fig.data:
        #     K.append(i.to_plotly_json()['marker']['color'])

        cluster = AgglomerativeClustering(n_clusters=K, affinity='euclidean', linkage='ward')
        cluster.fit_predict(X)
        labels = cluster.labels_
        labels = labels.tolist()
        labels=[str(x) for x in labels]
        print(labels)
        df_label = df_label[cols[:5]]
        print(df_label.head())
        print(cols[:5])
        fig = px.scatter_matrix(df_label,
                                dimensions=cols[:5],
                                color=labels,
                                opacity = 0.5)

        fig.update_layout(height=2300,width=2300,title="<b>Hierarchial scatter plot of the downsampled genes</b> \n(shows a maximum of 5 variable at a time for better vizualization)", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        return fig



    @staticmethod
    def Kmode_cluster_plot(df):
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1,'Core Genes')
            df['core_label'] = df['core_label'].replace(2,'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3,'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4,'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')
        # df = df.fillna(0)
        img = io.BytesIO()
        plt.style.use('seaborn-whitegrid')
        fig = plt.figure(figsize=(12, 6))
        cost = []
        K = range(1, 5)
        cols = df.columns
        for num_clusters in list(K):
            kmode = KModes(n_clusters=num_clusters, init="random", n_init=5, verbose=1)
            kmode.fit_predict(df)
            cost.append(kmode.cost_)
        # plt.plot(K, cost, 'bx-')
        # plt.xlabel('No. of clusters')
        # plt.ylabel('Cost')
        # plt.title('Elbow Method For Optimal k')
        kmode = KModes(n_clusters=3, init="random", n_init=5, verbose=1)
        clusters = kmode.fit_predict(df)
        df.insert(0, "Cluster", clusters, True)
        print(df.head())
        cluster_0 = df[df['Cluster'] == 0]
        cluster_1 = df[df['Cluster'] == 1]
        cluster_2 = df[df['Cluster'] == 2]
        print(len(df.columns) - 1)
        print(df.columns)
        print(cols)
        fig, axs = plt.subplots(nrows=len(df.columns) - 1, ncols=1, figsize=(12, 27))
        for i, column in enumerate(cols, start=1):
            print(i)
            for j in range(1):
                print(j)
                print([i - 1], [j])
                sns.countplot(x=df[column],order=df[column].value_counts().index,hue=df['Cluster'],ax=axs[i-1])
                title = "Count Plot of the Categorical Variables in each cluster"
                axs[i-1].set_title(title, fontsize=17, fontweight='bold')
                axs[i-1].set_ylabel('Count', fontsize=15)
                axs[i-1].xaxis.get_label().set_fontsize(15)
                axs[i-1].set_xticklabels(axs[i-1].get_xticklabels(), rotation=45, ha='right',
                                            rotation_mode='anchor')
                axs[i-1].legend(loc='upper right', ncol=1, bbox_to_anchor=(0, 0, 1, 1),
                                prop={'size': 15, 'weight': 'bold'}, fancybox=True, shadow=False, title='Gene Type',
                                title_fontsize=15).get_title().set_weight("bold")
                axs[i-1].xaxis.set_tick_params(color='b', width=4)
                axs[i-1].yaxis.set_tick_params(color='b', width=4)
                axs[i-1].xaxis.get_label().set_fontsize(15)
                axs[i-1].xaxis.get_label().set_fontweight('bold')
                axs[i - 1].yaxis.get_label().set_fontsize(15)
                axs[i - 1].yaxis.get_label().set_fontweight('bold')
                axs[i-1].set_xticklabels(axs[i-1].get_xticklabels(), fontdict={'fontweight': 'bold', 'fontsize': 13})
                axs[i-1].set_yticklabels(axs[i-1].get_yticklabels(), fontdict={'fontweight': 'bold', 'fontsize': 13})
        plt.tight_layout()
        plt.savefig('C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/plot.png')
        plt.savefig(img, format='png')
        img.seek(0)
        plot_url = base64.b64encode(img.getvalue()).decode()
        return plot_url

    @staticmethod
    def Heatmap_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        # if "ID" in df.columns:
        #     df = df.drop("ID", axis=1)
        print(df.head())
        if "ID" in df.columns:

            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        # new_data_df = df.dropna()
        print("Number of NaNs: {}".format(df.shape[0] - new_data_df.shape[0]))
        print(new_data_df)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        gene_expression = new_data_df.values.tolist()
        fig = go.Figure(data=go.Heatmap(
            z=gene_expression[0:100],
            x=cid[0:100],
            y=rid[0:100],
            hoverongaps=False,
            zmin = 0,
            zmax=50
        ))
        fig.update_layout(height=1300,width=1300,title="<b>Heatmap plot of the first 100 downsampled genes or entire gene set</b>", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        return fig

    @staticmethod
    def Codon_Heatmap_plots(df):
        # if "ID" in df.columns:
        #     df = df.drop("ID", axis=1)
        print(df.head())
        if "ID" in df.columns:

            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'},
                      inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        # new_data_df = df.dropna()
        print("Number of NaNs: {}".format(df.shape[0] - new_data_df.shape[0]))
        print(new_data_df)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        std_scaler = StandardScaler()
        new_data_df = std_scaler.fit_transform(new_data_df)
        new_data_df = pd.DataFrame(new_data_df, columns=cid)
        gene_expression = new_data_df.values.tolist()
        fig = go.Figure(data=go.Heatmap(
            z=gene_expression[0:100],
            x=cid[0:100],
            y=rid[0:100],
            hoverongaps=False
        ))
        fig.update_layout(height=1300,width=1300,title="<b>Heatmap plot of the first 100 downsampled genes or entire gene set</b>", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        return fig

    @staticmethod
    def Codon_Gene_PCA_2D_samples_plots(df):
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        temp_xx = new_data_df.T
        print(temp_xx.head())
        xx = temp_xx.values
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=2)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()

        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
        fig = px.scatter_matrix(yy, labels=labels, dimensions=range(2), color=temp_xx.index)
        fig.update_layout(
            height=600, width=1200,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    @staticmethod
    def Codon_Gene_PCA_3D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        temp_xx = new_data_df.T
        print(temp_xx.head())
        xx = temp_xx.values
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=3)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        fig = px.scatter_3d(yy, x=0, y=1, z=2, color= temp_xx.index,
                            title=f'Total Explained Variance: {var}',
                            labels={'0': 'PC1', '1': 'PC2', '2': 'PC3'})
        fig.update_layout(
            height=800, width=1400,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=15
        )
        return fig

    @staticmethod
    def Count_Heatmap_plots(df):
        # if "ID" in df.columns:
        #     df = df.drop("ID", axis=1)
        if "SignalP" in df.columns:
            df['SignalP'] = df['SignalP'].replace('SignalP',1)
            df['SignalP'] = df['SignalP'].replace('No SignalP', 0)
        print(df.head())
        if "ID" in df.columns:

            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        # new_data_df = df.dropna()
        print("Number of NaNs: {}".format(df.shape[0] - new_data_df.shape[0]))
        print(new_data_df)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        # std_scaler = StandardScaler()
        # new_data_df = std_scaler.fit_transform(new_data_df)
        # new_data_df = pd.DataFrame(new_data_df, columns=cid)
        gene_expression = new_data_df.values.tolist()
        fig = go.Figure(data=go.Heatmap(
            z=gene_expression[0:100],
            x=cid[0:100],
            y=rid[0:100],
            hoverongaps=False
        ))
        fig.update_layout(height=1300,width=1300,title="<b>Heatmap plot of the first 100 downsampled genes or entire gene set</b>", font_size=18,font_family="Arial Black",
            font_color="black", title_font=dict(size=20, family='Arial Black', color='black'))
        return fig

    @staticmethod
    def Count_Gene_PCA_2D_samples_plots(df):
        if "SignalP" in df.columns:
            df['SignalP'] = df['SignalP'].replace('SignalP',1)
            df['SignalP'] = df['SignalP'].replace('No SignalP', 0)
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        temp_xx = new_data_df.T
        print(temp_xx.head())
        xx = temp_xx.values
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=2)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()

        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
        fig = px.scatter_matrix(yy, labels=labels, dimensions=range(2), color=temp_xx.index)
        fig.update_layout(
            height=600, width=1200,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    @staticmethod
    def Count_Gene_PCA_3D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        if "SignalP" in df.columns:
            df['SignalP'] = df['SignalP'].replace('SignalP',1)
            df['SignalP'] = df['SignalP'].replace('No SignalP', 0)
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        temp_xx = new_data_df.T
        print(temp_xx.head())
        xx = temp_xx.values
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=3)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        fig = px.scatter_3d(yy, x=0, y=1, z=2, color=temp_xx.index,
                            title=f'Total Explained Variance: {var}',
                            labels={'0': 'PC1', '1': 'PC2', '2': 'PC3'})
        fig.update_layout(
            height=800, width=1400,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=15
        )
        return fig

    @staticmethod
    def Gene_PCA_2D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        pca_label= pd.read_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/PCA_Label.txt",sep='\t')
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        gene_expression = new_data_df.values.tolist()
        temp_xx = new_data_df.T
        # temp_xx = temp_xx.rename(columns = {'index':'Tissue'})
        temp_xx['Tissue'] = temp_xx.index
        print(pca_label.head())
        temp_xx = pd.merge(left=temp_xx, right=pca_label, how='left', left_on='Tissue', right_on='Tissue')
        print(temp_xx.head())
        new_temp_xx = pd.DataFrame()
        new_temp_xx['Type'] = temp_xx['Type']
        new_temp_xx['Labs'] = temp_xx['Labs']
        temp_xx = temp_xx.drop(['Tissue', 'Type','Labs'], axis=1)
        print(temp_xx.head())
        # xx = new_data_df.values.T
        xx = temp_xx.values
        # scaler = StandardScaler()
        # xx = scaler.fit_transform(xx)
        pca = PCA(n_components=2)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        # plt.figure(figsize=(16, 10))
        img = io.BytesIO()

        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
        fig = px.scatter_matrix(yy, labels=labels, dimensions=range(2), color=new_temp_xx['Type'],symbol=new_temp_xx['Labs'])
        fig.update_layout(
            height=600, width=1200,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig


    @staticmethod
    def Gene_PCA_3D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        pca_label = pd.read_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/PCA_Label.txt",
                                sep='\t')
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "variants_effect_5_prime_UTR_premature_start_codon_gain_variant" in df.columns:
            df.rename(columns={'variants_effect_5_prime_UTR_premature_start_codon_gain_variant':'var_eff_5_UTR_SCG'}, inplace=True)
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        gene_expression = new_data_df.values.tolist()
        # xx = new_data_df.values.T
        temp_xx = new_data_df.T
        # temp_xx = temp_xx.rename(columns = {'index':'Tissue'})
        temp_xx['Tissue'] = temp_xx.index
        print(pca_label.head())
        temp_xx = pd.merge(left=temp_xx, right=pca_label, how='left', left_on='Tissue', right_on='Tissue')
        print(temp_xx.head())
        new_temp_xx = pd.DataFrame()
        new_temp_xx['Type'] = temp_xx['Type']
        new_temp_xx['Labs'] = temp_xx['Labs']
        temp_xx = temp_xx.drop(['Tissue', 'Type', 'Labs'], axis=1)
        print(temp_xx.head())
        # xx = new_data_df.values.T
        xx = temp_xx.values
        # scaler = StandardScaler()
        # xx = scaler.fit_transform(xx)
        pca = PCA(n_components=3)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        # plt.figure(figsize=(16, 10))
        img = io.BytesIO()
        fig = px.scatter_3d(yy, x=0, y=1, z=2,  color=new_temp_xx['Type'],symbol=new_temp_xx['Labs'], title=f'Total Explained Variance: {var}',
                            labels={'0': 'PC1', '1': 'PC2', '2': 'PC3'})
        # fig.show()
        fig.update_layout(
            height=800, width=1400,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=15
        )
        return fig

    @staticmethod
    def Protein_PCA_2D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        pca_label = pd.read_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Protein_PCA_Label.txt",
                                sep='\t')
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        print("Number of NaNs: {}".format(new_data_df.shape[0] - new_data_df.shape[0]))
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        gene_expression = new_data_df.values.tolist()
        temp_xx = new_data_df.T
        # temp_xx = temp_xx.rename(columns = {'index':'Tissue'})
        temp_xx['Tissue'] = temp_xx.index
        print(pca_label.head())
        temp_xx = pd.merge(left=temp_xx, right=pca_label, how='left', left_on='Tissue', right_on='Tissue')
        print(temp_xx.head())
        new_temp_xx = pd.DataFrame()
        new_temp_xx['Type'] = temp_xx['Type']
        new_temp_xx['Labs'] = temp_xx['Labs']
        temp_xx = temp_xx.drop(['Tissue', 'Type', 'Labs'], axis=1)
        print(temp_xx.shape)
        print(new_temp_xx.head())
        print(new_temp_xx.shape)
        # xx = new_data_df.values.T
        xx = temp_xx.values
        # scaler = StandardScaler()
        # xx = scaler.fit_transform(xx)
        pca = PCA(n_components=2)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }

        fig = px.scatter_matrix(yy, labels=labels, dimensions=range(2), color=list(new_temp_xx['Type']),
                                symbol=list(new_temp_xx['Labs']))
        fig.update_layout(
            height=600, width=1200,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    @staticmethod
    def Protein_PCA_3D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        pca_label = pd.read_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/Protein_PCA_Label.txt",
                                sep='\t')
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        gene_expression = new_data_df.values.tolist()
        # xx = new_data_df.values.T
        temp_xx = new_data_df.T
        # temp_xx = temp_xx.rename(columns = {'index':'Tissue'})
        temp_xx['Tissue'] = temp_xx.index
        print(pca_label.head())
        temp_xx = pd.merge(left=temp_xx, right=pca_label, how='left', left_on='Tissue', right_on='Tissue')
        print(temp_xx.head())
        new_temp_xx = pd.DataFrame()
        new_temp_xx['Type'] = temp_xx['Type']
        new_temp_xx['Labs'] = temp_xx['Labs']
        temp_xx = temp_xx.drop(['Tissue', 'Type', 'Labs'], axis=1)
        print(temp_xx.head())
        # xx = new_data_df.values.T
        xx = temp_xx.values
        # scaler = StandardScaler()
        # xx = scaler.fit_transform(xx)
        pca = PCA(n_components=3)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        fig = px.scatter_3d(yy, x=0, y=1, z=2, color=new_temp_xx['Type'], symbol=new_temp_xx['Labs'],
                            title=f'Total Explained Variance: {var}',
                            labels={'0': 'PC1', '1': 'PC2', '2': 'PC3'})
        fig.update_layout(
            height=800, width=1400,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=15
        )
        return fig

    @staticmethod
    def PCA_2D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        gene_expression = new_data_df.values.tolist()
        xx = new_data_df.values.T
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=2)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
        fig = px.scatter_matrix(yy, labels=labels, dimensions=range(2), color=cid)
        fig.update_layout(
            height=600, width=1200,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig


    @staticmethod
    def PCA_3D_samples_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        label_df = df
        print(label_df.head())
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "classical_label" in df.columns:
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        gene_expression = new_data_df.values.tolist()
        xx = new_data_df.values.T
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=3)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        fig = px.scatter_3d(yy, x=0, y=1, z=2, color=cid,
                            title=f'Total Explained Variance: {var}',
                            labels={'0': 'PC1', '1': 'PC2', '2': 'PC3'})
        fig.update_layout(
            height=800, width=1400,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=15
        )
        return fig

    def ENC_plot(df):
        if "ID" in df.columns:
            df = df.drop("ID", axis=1)
        img = io.BytesIO()
        plt.style.use('seaborn-whitegrid')
        if "classical_label" in df.columns:
            target = df["classical_label"]
            df = df.drop(["classical_label"], axis=1)
        elif "core_label" in df.columns:
            target = df["core_label"]
            df = df.drop(["core_label"], axis=1)
        elif "Origin" in df.columns:
            target = df["Origin"]
            df = df.drop(["Origin"], axis=1)
        df = df[["Nc","GC3s"]]
        print((2+ 0.416+(29/(math.pow( 0.416,2)+math.pow((1-0.416),2)))))
        df['curve'] = df['GC3s'].apply(lambda x : (2+x+(29/(math.pow(x,2)+(math.pow((1-x),2))))))
        print(df.head())
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        plt.figure(figsize=(8, 6))
        plt.scatter(df['GC3s'],df['Nc'],c=target, cmap='plasma', alpha=0.4, edgecolors='black', s=65);
        plt.scatter(df['GC3s'], df['curve'], color='red')
        plt.xlabel('GC3s',fontsize=15,fontweight='bold')
        plt.ylabel('Nc',fontsize=15,fontweight='bold')
        plt.suptitle('ENC-plot analysis', fontsize=15,fontweight='bold')
        plt.tight_layout()
        plt.savefig('C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/plot.png')
        plt.savefig(img, format='png')
        img.seek(0)
        plot_url = base64.b64encode(img.getvalue()).decode()
        return plot_url

    def dinucleic_box_plots(df):
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            classical_val = df.loc[df['classical_label'] == 'classical Genes']['value']
            other_val = df.loc[df['classical_label'] == 'Other Genes']['value']
            ttest, pval = ttest_ind(classical_val, other_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in both Classical and Other genes" \
                        " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % ("dinucleotide count", pval)
            else:
                p_val_label = "There is no significant difference of %s in both Classical and Other genes" \
                  " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % ("dinucleotide count", pval)

            fig = px.box(df, x="variable", y="value", color="classical_label",title= p_val_label)
            fig.update_traces(quartilemethod="exclusive")
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='count')
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')

            core_val = df.loc[df['core_label'] == 'Core Genes']['value']
            nearCore_val = df.loc[df['core_label'] == 'Near-Core Genes']['value']
            Dispensible_val = df.loc[df['core_label'] == 'Dispensable Genes']['value']
            Private_val = df.loc[df['core_label'] == 'Private Genes']['value']
            None_val = df.loc[df['core_label'] == 'None']['value']

            anovatest, pval = stats.f_oneway(core_val, nearCore_val,Dispensible_val,Private_val,None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                              "dinucleotide count", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                              "dinucleotide count", pval)

            fig = px.box(df, x="variable", y="value", color="core_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='count')
            # fig.show()
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')

            WGD_val = df.loc[df['Origin'] == 'WGD']['value']
            Tandem_val = df.loc[df['Origin'] == 'Tandem']['value']
            Both_val = df.loc[df['Origin'] == 'Both']['value']
            None_val = df.loc[df['Origin'] == 'None']['value']

            anovatest, pval = stats.f_oneway(WGD_val, Tandem_val,Both_val,None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                              "dinucleotide count", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                              "dinucleotide count", pval)

            fig = px.box(df, x="variable", y="value", color="Origin",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='count')
            # fig.show()
        if "None" in df.columns:
            df['None'] = df['None'].replace(1, 'All Genes')
            fig = px.box(df, x="variable", y="value")
            fig.update_traces(quartilemethod="exclusive")
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='count')
        fig.update_layout(height=600, width=1400,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    def frequency_dinucleic_box_plots(df, cols=3, width=20, height=800, hspace=0.45, wspace=0.5):
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            classical_val = df.loc[df['classical_label'] == 'classical Genes']['value']
            other_val = df.loc[df['classical_label'] == 'Other Genes']['value']
            ttest, pval = ttest_ind(classical_val, other_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % (
                              "dinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % (
                              "dinucleotide frequency", pval)
            fig = px.box(df, x="variable", y="value", color="classical_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='frequency')
            # fig.show()
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')

            core_val = df.loc[df['core_label'] == 'Core Genes']['value']
            nearCore_val = df.loc[df['core_label'] == 'Near-Core Genes']['value']
            Dispensible_val = df.loc[df['core_label'] == 'Dispensable Genes']['value']
            Private_val = df.loc[df['core_label'] == 'Private Genes']['value']
            None_val = df.loc[df['core_label'] == 'None']['value']

            anovatest, pval = stats.f_oneway(core_val, nearCore_val, Dispensible_val, Private_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                                  "dinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                                  "dinucleotide frequency", pval)
            fig = px.box(df, x="variable", y="value", color="core_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='frequency')
            # fig.show()
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')

            WGD_val = df.loc[df['Origin'] == 'WGD']['value']
            Tandem_val = df.loc[df['Origin'] == 'Tandem']['value']
            Both_val = df.loc[df['Origin'] == 'Both']['value']
            None_val = df.loc[df['Origin'] == 'None']['value']

            anovatest, pval = stats.f_oneway(WGD_val, Tandem_val, Both_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                                  "dinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                                  "dinucleotide frequency", pval)
            fig = px.box(df, x="variable", y="value", color="Origin",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='frequency')

        elif "None" in df.columns:
            df['None'] = df['None'].replace(1, 'All Genes')
            fig = px.box(df, x="variable", y="value")
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dinucleotide')
            fig.update_yaxes(title_text='frequency')

            # fig.show()
        fig.update_layout(height=600, width=1400,
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    def dipeptide_box_plots(df):
        df = df.replace(r'^\s*$', 0, regex=True)
        print(df.columns)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            classical_val = df.loc[df['classical_label'] == 'classical Genes']['value']
            other_val = df.loc[df['classical_label'] == 'Other Genes']['value']
            ttest, pval = ttest_ind(classical_val, other_val)
            print("p-value", pval)
            # pval = 0.04
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 dipeptides </b>" % (
                                  "dipeptides frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 dipeptides </b>" % (
                                  "dipeptides frequency", pval)

            fig = px.box(df, x="variable", y="value", color="classical_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dipeptide')
            fig.update_layout(height=700, width=1550)
            fig.update_yaxes(title_text='frequency')
            # fig.show()
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')

            core_val = df.loc[df['core_label'] == 'Core Genes']['value']
            nearCore_val = df.loc[df['core_label'] == 'Near-Core Genes']['value']
            Dispensible_val = df.loc[df['core_label'] == 'Dispensable Genes']['value']
            Private_val = df.loc[df['core_label'] == 'Private Genes']['value']
            None_val = df.loc[df['core_label'] == 'None']['value']

            anovatest, pval = stats.f_oneway(core_val, nearCore_val, Dispensible_val, Private_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 dipeptides </b>" % (
                                  "dipeptide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 dipeptides </b>" % (
                                  "dipeptide frequency", pval)

            fig = px.box(df, x="variable", y="value", color="core_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dipeptide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=700, width=1550)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')

            WGD_val = df.loc[df['Origin'] == 'WGD']['value']
            Tandem_val = df.loc[df['Origin'] == 'Tandem']['value']
            Both_val = df.loc[df['Origin'] == 'Both']['value']
            None_val = df.loc[df['Origin'] == 'None']['value']

            anovatest, pval = stats.f_oneway(WGD_val, Tandem_val, Both_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 dipeptides </b>" % (
                                  "dipeptide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 dipeptides </b>" % (
                                  "dipeptide frequency", pval)

            fig = px.box(df, x="variable", y="value", color="Origin", title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dipeptide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=700, width=1550)
            # fig.show()
        elif "None" in df.columns:
            df['None'] = df["None"].replace(1, 'All Genes')

            fig = px.box(df, x="variable", y="value")
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='dipeptide')
            fig.update_layout(height=700, width=1550)
            fig.update_yaxes(title_text='frequency')
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    def tripeptide_box_plots(df):
        df = df.replace(r'^\s*$', 0, regex=True)
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            classical_val = df.loc[df['classical_label'] == 'classical Genes']['value']
            other_val = df.loc[df['classical_label'] == 'Other Genes']['value']
            ttest, pval = ttest_ind(classical_val, other_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 tripeptides </b>" % (
                                  "tripeptides frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of first 20 tripeptides </b>" % (
                                  "tripeptides frequency", pval)

            fig = px.box(df, x="variable", y="value", color="classical_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='tripeptide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=700, width=2000)
            # fig.show()
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')

            core_val = df.loc[df['core_label'] == 'Core Genes']['value']
            nearCore_val = df.loc[df['core_label'] == 'Near-Core Genes']['value']
            Dispensible_val = df.loc[df['core_label'] == 'Dispensable Genes']['value']
            Private_val = df.loc[df['core_label'] == 'Private Genes']['value']
            None_val = df.loc[df['core_label'] == 'None']['value']

            anovatest, pval = stats.f_oneway(core_val, nearCore_val, Dispensible_val, Private_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots first 20 tripeptides </b>" % (
                                  "tripeptide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots first 20 tripeptides </b>" % (
                                  "tripeptide frequency", pval)

            fig = px.box(df, x="variable", y="value", color="core_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='tripeptide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=700, width=2000)

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')

            WGD_val = df.loc[df['Origin'] == 'WGD']['value']
            Tandem_val = df.loc[df['Origin'] == 'Tandem']['value']
            Both_val = df.loc[df['Origin'] == 'Both']['value']
            None_val = df.loc[df['Origin'] == 'None']['value']

            anovatest, pval = stats.f_oneway(WGD_val, Tandem_val, Both_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots first 20 tripeptides </b>" % (
                                  "tripeptide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots first 20 tripeptides </b>" % (
                                  "tripeptide frequency", pval)

            fig = px.box(df, x="variable", y="value", color="Origin", title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='tripeptide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=700, width=2000)
            # fig.show()
        if "None" in df.columns:
            df['None'] = df['None'].replace(1, 'All Genes')
            fig = px.box(df, x="variable", y="value")
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='tripeptide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=700, width=2000)
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    def trinucleic_box_plots(df):
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')

            classical_val = df.loc[df['classical_label'] == 'classical Genes']['value']
            other_val = df.loc[df['classical_label'] == 'Other Genes']['value']
            ttest, pval = ttest_ind(classical_val, other_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % (
                                  "trinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % (
                                  "trinucleotide frequency", pval)

            fig = px.box(df, x="variable", y="value", color="classical_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='count')
            fig.update_layout(height=650, width=2000)
            # fig.show()
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')

            core_val = df.loc[df['core_label'] == 'Core Genes']['value']
            nearCore_val = df.loc[df['core_label'] == 'Near-Core Genes']['value']
            Dispensible_val = df.loc[df['core_label'] == 'Dispensable Genes']['value']
            Private_val = df.loc[df['core_label'] == 'Private Genes']['value']
            None_val = df.loc[df['core_label'] == 'None']['value']

            anovatest, pval = stats.f_oneway(core_val, nearCore_val, Dispensible_val, Private_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                                  "trinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                                  "trinucleotide frequency", pval)

            fig = px.box(df, x="variable", y="value", color="core_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='count')
            fig.update_layout(height=650, width=3500)

        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')

            WGD_val = df.loc[df['Origin'] == 'WGD']['value']
            Tandem_val = df.loc[df['Origin'] == 'Tandem']['value']
            Both_val = df.loc[df['Origin'] == 'Both']['value']
            None_val = df.loc[df['Origin'] == 'None']['value']

            anovatest, pval = stats.f_oneway(WGD_val, Tandem_val, Both_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                                  "trinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                                  "trinucleotide frequency", pval)

            fig = px.box(df, x="variable", y="value", color="Origin", title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='count')
            fig.update_layout(height=650, width=3500)
            # fig.show()
        if "None" in df.columns:
            df['None'] = df['None'].replace(1, 'All Genes')
            fig = px.box(df, x="variable", y="value")
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='count')
            fig.update_layout(height=650, width=2000)
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    def frequency_trinucleic_box_plots(df):
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')

            classical_val = df.loc[df['classical_label'] == 'classical Genes']['value']
            other_val = df.loc[df['classical_label'] == 'Other Genes']['value']
            ttest, pval = ttest_ind(classical_val, other_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % (
                                  "trinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in both Classical and Other genes" \
                              " as p-value = %.2E <br><b>Box plots of Classical genes </b>" % (
                                  "trinucleotide frequency", pval)
            fig = px.box(df, x="variable", y="value", color="classical_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=650, width=2000)
            # fig.show()
        elif "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            df['core_label'] = df['core_label'].replace(np.nan, 'None')
            core_val = df.loc[df['core_label'] == 'Core Genes']['value']
            nearCore_val = df.loc[df['core_label'] == 'Near-Core Genes']['value']
            Dispensible_val = df.loc[df['core_label'] == 'Dispensable Genes']['value']
            Private_val = df.loc[df['core_label'] == 'Private Genes']['value']
            None_val = df.loc[df['core_label'] == 'None']['value']

            anovatest, pval = stats.f_oneway(core_val, nearCore_val, Dispensible_val, Private_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                                  "trinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of core genes" \
                              " as p-value = %.2E <br><b>Box plots of Core genes </b>" % (
                                  "trinucleotide frequency", pval)
            fig = px.box(df, x="variable", y="value", color="core_label",title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=650, width=3500)
        elif "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            df['Origin'] = df['Origin'].replace(np.nan, 'None')

            WGD_val = df.loc[df['Origin'] == 'WGD']['value']
            Tandem_val = df.loc[df['Origin'] == 'Tandem']['value']
            Both_val = df.loc[df['Origin'] == 'Both']['value']
            None_val = df.loc[df['Origin'] == 'None']['value']

            anovatest, pval = stats.f_oneway(WGD_val, Tandem_val, Both_val, None_val)
            print("p-value", pval)
            if pval < 0.05:
                print("we reject null hypothesis")
                p_val_label = "There is significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                                  "trinucleotide frequency", pval)
            else:
                p_val_label = "There is no significant difference of %s in the different category of Origin genes" \
                              " as p-value = %.2E <br><b>Box plots of Origin genes </b>" % (
                                  "trinucleotide frequency", pval)
            fig = px.box(df, x="variable", y="value", color="Origin", title=p_val_label)
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=650, width=3500)
            # fig.show()
        if "None" in df.columns:
            df['None'] = df['None'].replace(1, 'All Genes')
            fig = px.box(df, x="variable", y="value")
            fig.update_traces(quartilemethod="exclusive")  # or "inclusive", or "linear" by default
            fig.update_xaxes(title_text='trinucleotide')
            fig.update_yaxes(title_text='frequency')
            fig.update_layout(height=650, width=2000)
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig

    @staticmethod
    def PCA_2D_label_plots(df):
        if "core_label" in df.columns:
            df = df[df['core_label'].notna()]
        if "Origin" in df.columns:
            df = df[df['Origin'].notna()]
        shapes = df['ID']
        color_data = []
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            color_data.extend(df['classical_label'])
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            # df['core_label'] = df['core_label'].replace(np.nan, 'None')
            color_data.extend(df['core_label'])
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            # df['core_label'] = df['core_label'].replace(np.nan, 'None')
            color_data.extend(df['Origin'])
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        features = new_data_df.columns
        xx = new_data_df.values
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=2)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
        fig = px.scatter_matrix(yy, labels=labels, dimensions=range(2), color=color_data,symbol=shapes,opacity=0.7)
        fig.update_layout(height=650, width=1300)
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=18
        )
        return fig


    @staticmethod
    def PCA_2D_label_biplot_plots(df):
        if "core_label" in df.columns:
            df = df[df['core_label'].notna()]
        if "Origin" in df.columns:
            df = df[df['Origin'].notna()]
        shapes = df['ID']
        color_data = []
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            color_data.extend(df['classical_label'])
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            # df['core_label'] = df['core_label'].replace(np.nan, 'None')
            color_data.extend(df['core_label'])
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1,'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3,'Both')
            # df['core_label'] = df['core_label'].replace(np.nan, 'None')
            color_data.extend(df['Origin'])
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        features = new_data_df.columns
        X = new_data_df[features]
        scaler = StandardScaler()
        xx = scaler.fit_transform(X)
        pca = PCA(n_components=2)
        components = pca.fit_transform(xx)

        loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
        labels = {
            str(i): f"PC {i + 1} ({var:.1f}%)"
            for i, var in enumerate(pca.explained_variance_ratio_ * 100)
        }
        fig = px.scatter(components,labels=labels,x=0, y=1, color=color_data,opacity=0.7)

        for i, feature in enumerate(features):
            fig.add_shape(
                type='line',
                x0=0, y0=0,
                x1=loadings[i, 0],
                y1=loadings[i, 1]
            )
            fig.add_annotation(
                x=loadings[i, 0],
                y=loadings[i, 1],
                ax=0, ay=0,
                xanchor="center",
                yanchor="bottom",
                text=feature,
            )
        fig.update_layout(height=650, width=1300)
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=11
        )
        return fig

    @staticmethod
    def PCA_3D_label_plots(df):
        if "core_label" in df.columns:
            df = df[df['core_label'].notna()]
        if "Origin" in df.columns:
            df = df[df['Origin'].notna()]
        shapes = df['ID']
        color_data = []
        if "ID" in df.columns:
            df = df.set_index('ID')
        if "classical_label" in df.columns:
            df['classical_label'] = df['classical_label'].replace(1, 'classical Genes')
            df['classical_label'] = df['classical_label'].replace(0, 'Other Genes')
            color_data.extend(df['classical_label'])
            df = df.drop(["classical_label"], axis=1)
        if "core_label" in df.columns:
            df['core_label'] = df['core_label'].replace(1, 'Core Genes')
            df['core_label'] = df['core_label'].replace(2, 'Near-Core Genes')
            df['core_label'] = df['core_label'].replace(3, 'Dispensable Genes')
            df['core_label'] = df['core_label'].replace(4, 'Private Genes')
            # df['core_label'] = df['core_label'].replace(np.nan, 'None')
            color_data.extend(df['core_label'])
            df = df.drop(["core_label"], axis=1)
        if "Origin" in df.columns:
            df['Origin'] = df['Origin'].replace(1, 'WGD')
            df['Origin'] = df['Origin'].replace(2, 'Tandem')
            df['Origin'] = df['Origin'].replace(3, 'Both')
            # df['core_label'] = df['core_label'].replace(np.nan, 'None')
            color_data.extend(df['Origin'])
            df = df.drop(["Origin"], axis=1)
        new_data_df = df.fillna(0)
        rid = list(new_data_df.index.values.tolist())
        print(rid[:5])
        cid = list(new_data_df.columns.values.tolist())
        print(cid[:5])
        features = new_data_df.columns
        xx = new_data_df.values
        scaler = StandardScaler()
        xx = scaler.fit_transform(xx)
        pca = PCA(n_components=3)
        yy = pca.fit_transform(xx)
        print(pca.explained_variance_ratio_)
        var = pca.explained_variance_ratio_.sum()
        fig = px.scatter_3d(yy, x=0, y=1, z=2, color=color_data,symbol=shapes,
                            title=f'Total Explained Variance: {var}',
                            labels={'0': 'PC1', '1': 'PC2', '2': 'PC3'},opacity=0.7)
        fig.update_layout(height=800, width=1400)
        fig.update_layout(
            legend=dict(title_font_family="Arial Black", title_font_size=22, title_font_color="black",
                        font=dict(size=18, family='Arial Black', color='black')),
            title_font=dict(size=20, family='Arial Black', color='black'),
            font_family="Arial Black",
            font_color="black",
            font_size=15
        )
        return fig














