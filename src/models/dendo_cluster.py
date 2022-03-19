from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.colors import rgb2hex, colorConverter
from scipy.cluster.hierarchy import set_link_color_palette
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt

from pylab import rcParams
rcParams['figure.figsize'] = 12, 9

import seaborn as sns
sns.set_style("whitegrid")

file = open("C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/cluster_dendo.html",'w')
cluster_classes = {}

class Clusters():
    @staticmethod
    def html(dict):
        html = '<table style="border: 0;">'
        for c in dict:
            hx = rgb2hex(colorConverter.to_rgb(c))
            html += '<tr style="border: 0;">' \
                    '<td style="background-color: {0}; ' \
                    'border: 0;">'.format(hx)
            html += c + '</td>'
            html += '<td style="border: 0">'
            html += str(dict[c])
            html += '</td></tr>'

        html += '</table>'

        return html

def get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        # print(c)
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))

    # cluster_classes = Clusters()
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l
    # print(cluster_classes)
    # print(Clusters.html(cluster_classes))
    print("Writing to file")
    file.write(Clusters.html(cluster_classes))
    # file.close()
    return Clusters.html(cluster_classes)

def get_clust_graph(df: object, label: object, numclust: object, transpose: object = False, dataname: object = None, save: object = False, xticksize: object = 8) -> object:
    if transpose==True:
        aml=df.transpose()
        xl="x-axis"
    else:
        aml=df
        xl="y-axis"

    data_dist = pdist(aml.transpose()) # computing the distance
    data_link = linkage(data_dist,  metric='euclidean', method='ward')#method="complete") # computing the linkage
    B = dendrogram(data_link, labels=label, p=numclust, truncate_mode="lastp", get_leaves=True,
                   count_sort='ascending', show_contracted=True)
    # get_cluster_classes(B)
    return get_cluster_classes(B)

def give_cluster_assigns(df, numclust, tranpose=True):
    if tranpose==True:
        data_dist = pdist(df.transpose())
        data_link = linkage(data_dist,  metric='euclidean', method='ward')
        cluster_assigns=pd.Series(sch.fcluster(data_link, numclust, criterion='maxclust', monocrit=None), index=df.columns)
    else:
        data_dist = pdist(df)
        data_link = linkage(data_dist,  metric='euclidean', method='ward')
        cluster_assigns=pd.Series(sch.fcluster(data_link, numclust, criterion='maxclust', monocrit=None), index=df.index)
    for i in range(1,numclust+1):
        print("Cluster ",str(i),": ( N =",len(cluster_assigns[cluster_assigns==i].index),")", ", ".join(list(cluster_assigns[cluster_assigns==i].index)))