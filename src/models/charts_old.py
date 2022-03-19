# # # import pandas as pd
# # # # # # col_file = pd.read_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/count_col.txt",sep='\t',header=None)
# # # # SBP_name = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/SBP_name.txt",'w')
# # # count_df = pd.read_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/collection_Copy/count.txt",sep='\t')
# # # # for name in count_df.columns:
# # # #     #print(name)
# # # #     if name.startswith("SBP"):
# # # #         #print(name)
# # # #         SBP_name.write(name+ "\n")
# # #
# # # col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/ARF_name.txt","r")
# # # col_name = []
# # # col_name.append('ID')
# # # for line in col_file:
# # #     col_name.append(line.strip())
# # # print(col_name)
# # #
# # # print(count_df[col_name].head())
# # # count_df['ARF'] = count_df[col_name].sum(axis=1)
# # # col_name.append('ARF')
# # # print(count_df[col_name].head())
# # #
# # # print(type(count_df['ARF_10'][0]))
# #
# # # # count_df.to_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/collection_Copy/count_1.txt", sep='\t', index=False, header=True)
# # #
# # # # import seaborn as sns
# # # # print(sns.color_palette(n_colors=4))
# # #
# # from bs4 import BeautifulSoup
# #
# # # soup = BeautifulSoup("<b></b>",features="html.parser")
# # soup = BeautifulSoup(open('C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/dendogram.html'), 'html.parser')
# # # original_tag = soup.b
# #
# # new_tag = soup.new_tag("a",**{'class':'download','download':'w3logo'},href="static/dist_plot/plot.png")
# # button_tag = soup.new_tag("button",**{'class':'btn'})
# # # i_tag = soup.new_tag("i",**{'class':'fa fa-download'})
# # button_tag .string = "Download text."
# # new_tag.append(button_tag)
# # # new_tag.append(i_tag)
# # # print(original_tag)
# # # <b><a href="http://www.example.com"></a></b>
# #
# #
# # soup.body.append(new_tag)
# #
# # with open("C:/Users/Shatabdi/PycharmProjects/Maize_FeatureStore/src/static/dist_plot/dendogram.html", "w") as file:
# #     file.write(str(soup))
# # # print(original_tag)
# #
# #
#
#
#
# ######################
# import csv, urllib.request
#
# url = 'http://winterolympicsmedals.com/medals.csv'
# # response = urllib.request.urlopen(url)
# # lines = [l.decode('utf-8') for l in response.readlines()]
# # cr = csv.reader(lines)
# #
# # for row in cr:
# #     print(row)
#
# import pandas as pd
# import io
# import requests
# url="http://winterolympicsmedals.com/medals.csv"
# s=requests.get(url).content
# c=pd.read_csv(io.StringIO(s.decode('utf-8')))
# print(c.head())


##########
# import plotly.express as px
# from sklearn.decomposition import PCA
# from sklearn import datasets
# from sklearn.preprocessing import StandardScaler
# import numpy as np
# df = px.data.iris()
# features = ['sepal_length', 'sepal_width', 'petal_length', 'petal_width']
# X = df[features]
#
# pca = PCA(n_components=2)
# components = pca.fit_transform(X)
#
# loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
#
# fig = px.scatter(components, x=0, y=1, color=df['species'])
#
# for i, feature in enumerate(features):
#     fig.add_shape(
#         type='line',
#         x0=0, y0=0,
#         x1=loadings[i, 0],
#         y1=loadings[i, 1]
#     )
#     fig.add_annotation(
#         x=loadings[i, 0],
#         y=loadings[i, 1],
#         ax=0, ay=0,
#         xanchor="center",
#         yanchor="bottom",
#         text=feature,
#     )
# fig.show()

###########
# import os
# print(os.getcwd())

# class hello():
#
# 	@staticmethod
#     def myname():
# 		return "my name is shatabdi"



# from flask import Flask
#
# app = Flask(__name__)
#
#
# @app.route("/")
# def index():
#     return "Hello, from WSGI"

#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Import the required Python modules and Flask libraries
from sklearn.utils import resample
import numpy as np
from scipy.sparse import coo_matrix

X = np.array([[1., 0.], [2., 1.], [0., 0.]])
y = np.array([0, 1, 2])
X_sparse = coo_matrix(X)
X, X_sparse, y = resample(X, X_sparse, y, random_state=0)
print(X)


