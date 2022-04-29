import json
import plotly
from common.database import Database
from models.seqFeature import DNA
from models.seqplot import plot_seq
import numpy as np
import pandas as pd
from flask import Flask, request, jsonify, render_template
import pickle
import subprocess
from flask import Flask, render_template, request
from models.forms import ContactForm
import pandas as pd
import sys
import os

__author__ = "shatabdi"

from flask import Flask, render_template, request, send_from_directory

col_file = open("/var/www/html/MFS/src/col_names/colfile_new.txt", "r")
col_name=[]
for line in col_file:
    col_name.append(line.strip())
# print(col_name)


app = Flask(__name__) #'__main__'
app.secret_key = 'secretKey'

model = pickle.load(open('/var/www/html/MFS/src/models/model.pkl', 'rb'))
scaler = pickle.load(open('/var/www/html/MFS/src/models/scaler.pkl', 'rb'))

model_seq = pickle.load(open('/var/www/html/MFS/src/models/seq_model.pkl', 'rb'))
scaler_seq = pickle.load(open('/var/www/html/MFS/src/models/seq_scaler.pkl', 'rb'))

@app.route('/') #www.mysite.com/api/
def base_method():
    return render_template('base.html')

@app.route('/model_advanced')
def models():
    return render_template('model_advanced.html')

@app.route('/model_basic')
def model_features():
    return render_template('model_basic.html')

@app.route('/features') #www.mysite.com/api/
def features_method():
    return render_template('featurebase.html')

@app.route('/DNAseq') #www.mysite.com/api/
def DNAseq_method():
    return render_template('Dna_seq.html')

@app.route('/Proteinseq') #www.mysite.com/api/
def Proteinseq_method():
    return render_template('Proteinseq.html')

@app.route('/Codon') #www.mysite.com/api/
def Codon_method():
    return render_template('Codon.html')

@app.route('/ProteinStructure') #www.mysite.com/api/
def ProteinStructure_method():
    return render_template('ProteinStructure.html')

@app.route('/ProteinLocalization') #www.mysite.com/api/
def ProteinLocalization_method():
    return render_template('ProteinLocalization.html')

@app.route('/Structure') #www.mysite.com/api/
def structure_method():
    return render_template('Gene_Features.html')

@app.route('/distance') #www.mysite.com/api/
def distance_method():
    return render_template('distance.html')

@app.route('/GeneExp') #www.mysite.com/api/
def GeneExpression_method():
    return render_template('GeneExp.html')

@app.route('/ProteinExp') #www.mysite.com/api/
def ProteinExpression_method():
    return render_template('ProteinExp.html')

@app.route('/ChromStates') #www.mysite.com/api/
def ChromatinStates_method():
    return render_template('ChromatinStates.html')

@app.route('/ATACseq') #www.mysite.com/api/
def ATACseq_method():
    return render_template('ATACseq.html')

@app.route('/HistoneModifications') #www.mysite.com/api/
def HistoneModifications_method():
    return render_template('HistoneModifications.html')

@app.route('/DNAmethylation') #www.mysite.com/api/
def DNAmethylation_method():
    return render_template('DNAmethylation.html')

@app.route('/miRNA') #www.mysite.com/api/
def miRNA_method():
    return render_template('miRNA.html')

@app.route('/PFAM') #www.mysite.com/api/
def PFAM_method():
    return render_template('PFAM.html')

@app.route('/Insertions') #www.mysite.com/api/
def Insertions_method():
    return render_template('Insertions.html')

@app.route('/TFbindingSite') #www.mysite.com/api/
def TFbindingSite_method():
    return render_template('TFbindingSite.html')


@app.route('/chromatinAssescibility') #www.mysite.com/api/
def chromatinAssescibility_method():
    return render_template('chromatinAssescibility.html')

@app.route('/TSS') #www.mysite.com/api/
def TSS_method():
    return render_template('TSS.html')

@app.route('/Enhancers') #www.mysite.com/api/
def Enhancers_method():
    return render_template('Enhancers.html')

@app.route('/TE') #www.mysite.com/api/
def TE_method():
    return render_template('TE.html')

@app.route('/Quadruplex') #www.mysite.com/api/
def Quadruplex_method():
    return render_template('Quadruplex.html')

@app.route('/Correlation') #www.mysite.com/api/
def Correlation_method():
    return render_template('Correlation.html')


@app.route('/Varionomic') #www.mysite.com/api/
def Varionomic_method():
    return render_template('Varionomic.html')

@app.route('/Other') #www.mysite.com/api/
def Other_method():
    return render_template('Other.html')

@app.route('/Labels') #www.mysite.com/api/
def Labels_method():
    return render_template('Labels.html')

@app.route('/exploratory') #www.mysite.com/api/
def features_analysis():
    return render_template('exploratory.html')

@app.before_first_request
def initialize_databases():
    Database.initialize()

def Convert(lst):
    res_dct = {lst[i]: 1 for i in range(0, len(lst), 1)}
    return res_dct

def iterator2dataframes(iterator, chunk_size: int):

  records = []
  list_cur = []
  for i, record in enumerate(iterator):
    records.append(record)
    if i % chunk_size == chunk_size - 1:
      list_cur.append(records)
      records = []
  if records:
      list_cur.append(records)
  return list_cur

def Extract(lst):
    return [item[0] for item in lst]

@app.route('/predict', methods=['POST'])
def predict():
    features = [float(x) for x in request.form.values()]
    print(features)
    final_features = np.array(features).reshape(1,-1)
    final_features = scaler.transform(final_features)
    print(final_features)
    prediction = model.predict(final_features)
    print(prediction)
    y_probabilities_test = model.predict_proba(final_features)
    y_prob_success = y_probabilities_test[:, 1]
    print("final features", final_features)
    print("prediction:", prediction)
    output = round(prediction[0], 2)
    y_prob = round(y_prob_success[0], 3)
    print(output)

    if output == 0:
        return render_template('model_advanced.html',
                               prediction_text='THE GENE IS MORE LIKELY TO BE NON-CORE WITH PROBABILITY VALUE  {}'.format(
                                   y_prob))
    else:
        return render_template('model_advanced.html',
                               prediction_text='THE GENE IS MORE LIKELY TO BE CORE WITH PROBABILITY VALUE  {}'.format(
                                   y_prob))


@app.route('/predict_api', methods=['POST'])
def predict_api():
    data = request.get_json(force=True)
    prediction = model.predict([np.array(list(data.values()))])

    output = prediction[0]
    return jsonify(output)

@app.route('/predict_basic', methods=['POST'])
def predict_basic():
    features = [float(x) for x in request.form.values()]
    print(features)
    final_features = np.array(features).reshape(1,-1)
    final_features = scaler_seq.transform(final_features)
    print(final_features)
    prediction = model_seq.predict(final_features)
    print(prediction)
    y_probabilities_test = model_seq.predict_proba(final_features)
    y_prob_success = y_probabilities_test[:, 1]
    print("final features", final_features)
    print("prediction:", prediction)
    output = round(prediction[0], 2)
    y_prob = round(y_prob_success[0], 3)
    print(output)

    if output == 0:
        return render_template('model_basic.html',
                               prediction_text='THE GENE IS MORE LIKELY TO BE NON-CORE WITH PROBABILITY VALUE  {}'.format(
                                   y_prob))
    else:
        return render_template('model_basic.html',
                               prediction_text='THE GENE IS MORE LIKELY TO BE CORE WITH PROBABILITY VALUE  {}'.format(
                                   y_prob))


@app.route('/predict_basic_api', methods=['POST'])
def predict_basic_api():
    data = request.get_json(force=True)
    prediction = model_seq.predict([np.array(list(data.values()))])

    output = prediction[0]
    return jsonify(output)

@app.route('/handle_data', methods=['POST'])
def handle_data():
    projectpath = request.form['projectFilepath']
    print(projectpath)
    # Change accordingly to your Rscript.exe & R script path
    r_path = "/usr/lib64/R/bin/Rscript"
    script_path = "/var/www/html/MFS/src/models/max.R"
    # Used as input arguments to the R code
    # args = ["42", "15", "87", "62"]
    args = [str(projectpath)]
    # Execute command
    cmd = [r_path, script_path] + args
    print(cmd)
    result = subprocess.check_output(cmd, universal_newlines=True)
    # Display result
    print("The result is:", result)
    seq_feature = pd.read_csv("/var/www/html/MFS/src/trying_protr.txt",
                              sep = '\t')
    seq_feature.reset_index(inplace=True)
    seq_feature = seq_feature[["TC_SAL","TC_PAP","TC_LPL","TC_RSA","TC_VVA","TC_ARG"
        ,"CTriad_VS111","CTDD_prop7.G3.residue0","CTDD_prop1.G2.residue0","CTDD_prop5.G2.residue0",
                               "CTDD_prop2.G3.residue0","CTDD_prop7.G1.residue0","PAAC_Xc1.P","PAAC_Xc1.A"]]
    # return render_template('model_basic.html')
    value_TC_SAL = seq_feature.loc[0,"TC_SAL"]
    value_TC_PAP = seq_feature.loc[0, "TC_PAP"]
    value_TC_LPL = seq_feature.loc[0, "TC_LPL"]
    value_TC_RSA = seq_feature.loc[0, "TC_RSA"]
    value_TC_VVA = seq_feature.loc[0, "TC_VVA"]
    value_TC_ARG = seq_feature.loc[0, "TC_ARG"]
    value_CTriad_VS111 = seq_feature.loc[0, "CTriad_VS111"]
    value_CTDD_prop7_G3_residue0 = seq_feature.loc[0, "CTDD_prop7.G3.residue0"]
    value_CTDD_prop1_G2_residue0 = seq_feature.loc[0, "CTDD_prop1.G2.residue0"]
    value_CTDD_prop5_G2_residue0 = seq_feature.loc[0, "CTDD_prop5.G2.residue0"]
    value_CTDD_prop2_G3_residue0 = seq_feature.loc[0, "CTDD_prop2.G3.residue0"]
    value_CTDD_prop7_G1_residue0 = seq_feature.loc[0, "CTDD_prop7.G1.residue0"]
    value_PAAC_Xc1_P = seq_feature.loc[0, "PAAC_Xc1.P"]
    value_PAAC_Xc1_A = seq_feature.loc[0, "PAAC_Xc1.A"]
    DNApath = request.form['DNAcoding']
    print(DNApath)

    script_path = "/var/www/html/MFS/src/models/maxDNA.R"
    DNAargs = [str(DNApath)]
    cmd2 = [r_path, script_path] + DNAargs
    print(cmd2)
    result2 = subprocess.check_output(cmd2, universal_newlines=True)
    print("The result2 is:", result2)
    DNAseq_feature = pd.read_csv(
        "/var/www/html/MFS/src/DNA.txt",
        sep='\t')
    DNAseq_feature.reset_index(inplace=True)
    DNAseq_feature = DNAseq_feature[["PseDNC_Xc1.AT", "PseKNC_3_Xc1.ACG", "PseDNC_Xc1.GG", "PseKNC_3_Xc1.TAG",
                                     "PseKNC_3_Xc1.TGT"]]
    print(DNAseq_feature)
    value_PseDNC_Xc1_AT = DNAseq_feature.loc[0, "PseDNC_Xc1.AT"]
    value_PseKNC_3_Xc1_ACG = DNAseq_feature.loc[0, "PseKNC_3_Xc1.ACG"]
    value_PseDNC_Xc1_GG = DNAseq_feature.loc[0, "PseDNC_Xc1.GG"]
    value_PseKNC_3_Xc1_TAG = DNAseq_feature.loc[0, "PseKNC_3_Xc1.TAG"]
    value_PseKNC_3_Xc1_TGT = DNAseq_feature.loc[0, "PseKNC_3_Xc1.TGT"]
    return render_template('model_basic.html',TC_SAL = float(value_TC_SAL), TC_PAP = float( value_TC_PAP),
                           TC_LPL = float(value_TC_LPL),TC_RSA = float(value_TC_RSA),TC_VVA = float(value_TC_VVA)
                           ,TC_ARG = float(value_TC_ARG), CTriad_VS111 = float( value_CTriad_VS111),CTDD_prop7_G3_residue0 = float(value_CTDD_prop7_G3_residue0)
                           ,CTDD_prop1_G2_residue0 = float(value_CTDD_prop1_G2_residue0 ), CTDD_prop5_G2_residue0 = float(value_CTDD_prop5_G2_residue0),CTDD_prop2_G3_residue0 = float(value_CTDD_prop2_G3_residue0)
                           ,CTDD_prop7_G1_residue0 = float(value_CTDD_prop7_G1_residue0),PAAC_Xc1_P = float(value_PAAC_Xc1_P),PAAC_Xc1_A = float(value_PAAC_Xc1_A)
                           , PseDNC_Xc1_AT=float(value_PseDNC_Xc1_AT), PseKNC_3_Xc1_ACG=float(value_PseKNC_3_Xc1_ACG),
                           PseDNC_Xc1_GG=float(value_PseDNC_Xc1_GG), PseKNC_3_Xc1_TAG=float(value_PseKNC_3_Xc1_TAG),
                           PseKNC_3_Xc1_TGT=float(value_PseKNC_3_Xc1_TGT)
                           )



@app.route('/Structure',methods=['POST']) #www.mysite.com/api/
def features_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.structure(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/structure_plot',methods=['POST']) #www.mysite.com/api/
def structure_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.structure_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots"):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[0] == "Count_and_distribution"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_structure_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_structure_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_structure_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots"):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[0] == "Count_and_distribution"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)


@app.route('/distance',methods=['POST']) #www.mysite.com/api/
def disatnce_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.structure(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/distance_plot',methods=['POST']) #www.mysite.com/api/
def distance_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.structure_plot(select, analysis)
        if (analysis[0] == "downsampled_dendogram_plots"):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_distance_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_distance_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_structure_plot(select, analysis, cluster)
        if (analysis[0] == "downsampled_dendogram_plots"):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/DNAseq',methods=['POST']) #www.mysite.com/api/
def DNA_Seq():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.DNA_seq(select, analysis)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/DNAseq_plot', methods=['POST'])  # www.mysite.com/api/
def DNAseq_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        print(str(select))
        print(str(analysis))
        plot_sns = plot_seq.seq_plot(select,analysis)
        # return render_template('distance_plot.html', plot_graph=plot_sns)
        graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_DNAseq_plot', methods=['POST'])  # www.mysite.com/api/
def downsampled_DNAseq_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        print(str(select))
        print(str(analysis))
        plot_sns=plot_seq.downsampled_seq_plot(select, analysis)
        # return render_template('distance_plot.html', plot_graph=plot_sns)
        graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

# @app.route('/downsampled_DNAseq_plot', methods=['POST'])  # www.mysite.com/api/
# def downsampled_DNAseq_plot():
#     if request.method == 'POST':
#         select = request.form.getlist('comp_select')
#         analysis = request.form.getlist('comp_analysis')
#         print(str(select))
#         print(str(analysis))
#         plot_sns=plot_seq.same_len_downsampled_seq_plot(select, analysis)
#         return render_template('distance_plot.html', plot_graph=plot_sns)

@app.route('/Proteinseq',methods=['POST']) #www.mysite.com/api/
def Proteinseq():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.Protein_seq(select, analysis)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Proteinseq_plot', methods=['POST'])  # www.mysite.com/api/
def Proteinseq_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        print(str(select))
        print(str(analysis))
        plot_sns = plot_seq.pep_seq_plot(select,analysis)
        graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_Proteinseq_plot', methods=['POST'])  # www.mysite.com/api/
def downsampled_Proteinseq_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        print(str(select))
        print(str(analysis))
        plot_sns=plot_seq.pep_downsampled_seq_plot(select, analysis)
        graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)


@app.route('/Codon',methods=['POST']) #www.mysite.com/api/
def codon_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.codon(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Codon_plot',methods=['POST']) #www.mysite.com/api/
def codon_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.codon_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots"):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)



@app.route('/downsampled_Codon_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_codon_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_codon_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
        analysis[0] == 'Count_down_dendogram_plots'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)


@app.route('/ProteinStructure', methods=['POST'])  # www.mysite.com/api/
def ProteinStructure_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.ProteinStructure(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html', select=json_data)


@app.route('/ProteinStructure_plot', methods=['POST'])  # www.mysite.com/api/
def ProteinStructure_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.ProteinStructure_plot(select, analysis)
        # return render_template('distance_plot.html', plot_graph=plot_sns)
        if (analysis[0] == "downsampled_dendogram_plots"):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)


@app.route('/downsampled_ProteinStructure_plot', methods=['POST'])  # www.mysite.com/api/
def downsampled_ProteinStructure_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        print(cluster)
        plot_sns = plot_seq.downsampled_ProteinStructure_plot(select, analysis,cluster)
        # return render_template('distance_plot.html', plot_graph=plot_sns)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/ProteinLocalization', methods=['POST'])  # www.mysite.com/api/
def ProteinLocalization_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.ProteinLocalization(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html', select=json_data)


@app.route('/ProteinLocalization_plot', methods=['POST'])  # www.mysite.com/api/
def ProteinLocalization_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.ProteinLocalization_plot(select, analysis)
        # return render_template('distance_plot.html', plot_graph=plot_sns)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot" or analysis[0] == "Kmode_cluster_plots"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)


@app.route('/downsampled_ProteinLocalization_plot', methods=['POST'])  # www.mysite.com/api/
def downsampled_ProteinLocalization_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.downsampled_ProteinLocalization_plot(select, analysis)
        # return render_template('distance_plot.html', plot_graph=plot_sns)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"  or analysis[0] == "Kmode_cluster_plots"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/GeneExp',methods=['POST']) #www.mysite.com/api/
def GeneExp_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.GeneExp(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/GeneExp_plot',methods=['POST']) #www.mysite.com/api/
def GeneExp_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        plot_sns = plot_seq.GeneExp_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)



@app.route('/downsampled_GeneExp_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_GeneExp_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_GeneExp_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/ProteinExp',methods=['POST']) #www.mysite.com/api/
def ProteinExp_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.ProteinExp(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/ProteinExp_plot',methods=['POST']) #www.mysite.com/api/
def ProteinExp_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.ProteinExp_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)



@app.route('/downsampled_ProteinExp_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_ProteinExp_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_ProteinExp_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/ChromStates',methods=['POST']) #www.mysite.com/api/
def ChromStates_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.ChromStates(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/ChromStates_plot',methods=['POST']) #www.mysite.com/api/
def ChromStates_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        final = ""
        plot_sns = plot_seq.ChromStates_plot(select, analysis)
        if (analysis[0] == 'Leaf_chromatin_state_plot') or (analysis[0] == 'Leaf_fold_enrichment_plot') or  \
                (analysis[0] == 'Leaf_TES_fold_enrichment_plot') or (analysis[0] == 'Leaf_TSS_fold_enrichment_plot') :
            final = render_template('leaf_chromhmm_emissions_plot.html',user_image = plot_sns)
        elif (analysis[0] == 'ear_chromatin_state_plot') or (analysis[0] == 'ear_fold_enrichment_plot')or \
                (analysis[0] == 'ear_TES_fold_enrichment_plot') or (analysis[0] == 'ear_TSS_fold_enrichment_plot') :
            final = render_template('ear_chromhmm_plot.html',user_image = plot_sns)
        elif (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            final =  render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            final = render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            final= render_template('notdash.html', graphJSON=graphJSON)
        return final



@app.route('/downsampled_ChromStates_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_ChromStates_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        final = ""
        plot_sns = plot_seq.downsampled_ChromStates_plot(select, analysis,cluster)
        if (analysis[0] == 'Leaf_chromatin_state_plot') or (analysis[0] == 'ear_chromatin_state_plot') or \
                (analysis[0] == 'Leaf_fold_enrichment_plot') or (analysis[0] == 'ear_fold_enrichment_plot') or \
                (analysis[0] == 'Leaf_TES_fold_enrichment_plot') or (analysis[0] == 'Leaf_TSS_fold_enrichment_plot') or \
                (analysis[0] == 'ear_TES_fold_enrichment_plot') or (analysis[0] == 'ear_TSS_fold_enrichment_plot'):
            final = render_template('leaf_chromhmm_emissions_plot.html', user_image=plot_sns)
        elif (analysis[0] == 'Leaf_chromatin_state_plot') or (analysis[0] == 'ear_chromatin_state_plot') or \
                (analysis[0] == 'Leaf_fold_enrichment_plot') or (analysis[0] == 'ear_fold_enrichment_plot') or \
                (analysis[0] == 'Leaf_TES_fold_enrichment_plot') or (analysis[0] == 'Leaf_TSS_fold_enrichment_plot') or \
                (analysis[0] == 'ear_TES_fold_enrichment_plot') or (analysis[0] == 'ear_TSS_fold_enrichment_plot'):
            final = render_template('ear_chromhmm_plot.html', user_image=plot_sns)
        elif (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
              analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            final = render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            final = render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            final = render_template('notdash.html', graphJSON=graphJSON)
        return final

@app.route('/ATACseq',methods=['POST']) #www.mysite.com/api/
def ATACseq_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.ATACseq(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/ATACseq_plot',methods=['POST']) #www.mysite.com/api/
def ATACseq_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.ATACseq_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)



@app.route('/downsampled_ATACseq_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_ATACseq_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_ATACseq_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/HistoneModifications',methods=['POST']) #www.mysite.com/api/
def HistoneModifications_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.HistoneModifications(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/HistoneModifications_plot',methods=['POST']) #www.mysite.com/api/
def HistoneModifications_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.ATACseq_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)



@app.route('/downsampled_HistoneModifications_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_HistoneModifications_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_ATACseq_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/DNAmethylation',methods=['POST']) #www.mysite.com/api/
def DNAmethylation_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.DNAmethylation(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/DNAmethylation_plot',methods=['POST']) #www.mysite.com/api/
def DNAmethylation_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.ATACseq_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)



@app.route('/downsampled_DNAmethylation_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_DNAmethylation_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_ATACseq_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/miRNA',methods=['POST']) #www.mysite.com/api/
def miRNA_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.miRNA(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/miRNA_plot',methods=['POST']) #www.mysite.com/api/
def miRNA_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.count_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_miRNA_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_miRNA_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = [0]
        plot_sns = plot_seq.downsampled_count_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

# @app.route('/PFAM',methods=['POST']) #www.mysite.com/api/
# def PFAM_Table():
#     if request.method == 'POST':
#         select = request.form.getlist('comp_select')
#         select = [i for i in select if i != "none"]
#         analysis = request.form.getlist('comp_analysis')
#         analysis = [i for i in analysis if i != "none"]
#         print(str(select))
#         print(str(analysis))
#         json_data = DNA.PFAM(select, analysis)
#         # print(json_data)
#         print(type(json_data))
#         return render_template('table.html',select=json_data)
#
# @app.route('/PFAM_plot',methods=['POST']) #www.mysite.com/api/
# def PFAM_plot():
#     if request.method == 'POST':
#         select = request.form.getlist('comp_select')
#         analysis = request.form.getlist('comp_analysis')
#         plot_sns = plot_seq.count_plot(select,analysis)
#         if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
#                 analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
#             graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
#             return render_template('nodash_One.html', graphJSON=graphJSON)
#         elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
#             0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
#             return render_template('distance_plot.html', plot_graph=plot_sns)
#         else:
#             graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
#             return render_template('notdash.html', graphJSON=graphJSON)
#
# @app.route('/downsampled_PFAM_plot',methods=['POST']) #www.mysite.com/api/
# def downsampled_PFAM_plot():
#     if request.method == 'POST':
#         select = request.form.getlist('comp_select')
#         analysis = request.form.getlist('comp_analysis')
#         plot_sns = plot_seq.downsampled_count_plot(select, analysis)
#         if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
#                 analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
#             graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
#             return render_template('nodash_One.html', graphJSON=graphJSON)
#         elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
#             0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
#             return render_template('distance_plot.html', plot_graph=plot_sns)
#         else:
#             graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
#             return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/Insertions',methods=['POST']) #www.mysite.com/api/
def Insertions_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.Insertions(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Insertions_plot',methods=['POST']) #www.mysite.com/api/
def Insertions_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.count_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_Insertions_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_Insertions_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_count_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/TFbindingSite',methods=['POST']) #www.mysite.com/api/
def TFbindingSite_Table():
    if request.method == 'POST':
        if 'top_singleslct' in request.form:
            select = request.form.getlist('top_singleslct')
        else:
            select = request.form.getlist('comp_select')
        select.append(request.form.getlist('comp_select1'))
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.TFbindingSite(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/TFbindingSite_plot',methods=['POST']) #www.mysite.com/api/
def TFbindingSite_plot():
    if request.method == 'POST':
        if 'top_singleslct' in request.form:
            select = request.form.getlist('top_singleslct')
        else:
            select = request.form.getlist('comp_select')
        select.append(request.form.getlist('comp_select1'))
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.TF_count_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_TFbindingSite_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_TFbindingSite_plot():
    if request.method == 'POST':
        if 'top_singleslct' in request.form:
            select = request.form.getlist('top_singleslct')
        else:
            select = request.form.getlist('comp_select')
        select.append(request.form.getlist('comp_select1'))
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_TF_count_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/TSS',methods=['POST']) #www.mysite.com/api/
def TSS_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.TSS(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/TSS_plot',methods=['POST']) #www.mysite.com/api/
def TSS_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.count_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_TSS_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_TSS_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_count_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/Enhancers',methods=['POST']) #www.mysite.com/api/
def Enhancers_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.Enhancers(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Enhancers_plot',methods=['POST']) #www.mysite.com/api/
def Enhancers_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.count_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_Enhancers_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_Enhancers_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_count_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/TE',methods=['POST']) #www.mysite.com/api/
def TE_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.TE(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/TE_plot',methods=['POST']) #www.mysite.com/api/
def TE_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.count_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_TE_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_TE_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_count_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/Quadruplex',methods=['POST']) #www.mysite.com/api/
def Quadruplex_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.Quadruplex(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Quadruplex_plot',methods=['POST']) #www.mysite.com/api/
def Quadruplex_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.count_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_Quadruplex_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_Quadruplex_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_count_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/Correlation',methods=['POST']) #www.mysite.com/api/
def Correlation_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.Correlation(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Correlation_plot',methods=['POST']) #www.mysite.com/api/
def Correlation_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.ProteinLocalization_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_Correlation_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_Correlation_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.downsampled_ProteinLocalization_plot(select, analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/Varionomic',methods=['POST']) #www.mysite.com/api/
def Varionomic_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.Varionomic(select, analysis)
        # print(json_data)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Varionomic_plot',methods=['POST']) #www.mysite.com/api/
def Varionomic_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.Varionomic_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_Varionomic_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_Varionomic_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_Varionomic_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/Other',methods=['POST']) #www.mysite.com/api/
def Other_Table():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        select = [i for i in select if i != "none"]
        analysis = request.form.getlist('comp_analysis')
        analysis = [i for i in analysis if i != "none"]
        print(str(select))
        print(str(analysis))
        json_data = DNA.Other(select, analysis)
        print(type(json_data))
        return render_template('table.html',select=json_data)

@app.route('/Other_plot',methods=['POST']) #www.mysite.com/api/
def Other_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        plot_sns = plot_seq.others_plot(select,analysis)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/downsampled_Other_plot',methods=['POST']) #www.mysite.com/api/
def downsampled_Other_plot():
    if request.method == 'POST':
        select = request.form.getlist('comp_select')
        analysis = request.form.getlist('comp_analysis')
        cluster = request.form.getlist('cluster')
        plot_sns = plot_seq.downsampled_others_plot(select, analysis,cluster)
        if (analysis[0] == "downsampled_dendogram_plots" or analysis[0] == "down_dendogram_plots" or
                analysis[0] == 'Count_down_dendogram_plots' or analysis[0] == 'Count_downsampled_dendogram_plot'):
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
            return render_template('nodash_One.html', graphJSON=graphJSON)
        elif (analysis[0] == "Violin_plots" or analysis[0] == "Correlation_plots" or analysis[
            0] == "Count_and_distribution" or analysis[0] == "ENC_plot"):
            return render_template('distance_plot.html', plot_graph=plot_sns)
        else:
            graphJSON = json.dumps(plot_sns, cls=plotly.utils.PlotlyJSONEncoder)
        return render_template('notdash.html', graphJSON=graphJSON)

@app.route('/contactus', methods=["GET", "POST"])
def get_contact():
    form = ContactForm()
    if request.method == 'POST':
        name = request.form["name"]
        email = request.form["email"]
        subject = request.form["subject"]
        message = request.form["message"]
        res = pd.DataFrame({'name': name, 'email': email, 'subject': subject, 'message': message}, index=[0])
        res.to_csv('/var/www/html/MFS/src/feedback/contactusMessage.csv', mode='a'
                   ,header=False)
        return render_template('contact.html', form=form)
    else:
        return render_template('contact.html', form=form)

@app.route('/FAQ')
def FAQS():
    return render_template('FAQ.html')

@app.route('/data_sources')
def data_sources():
    return render_template('data_sources.html')

@app.route('/tool_sources')
def tool_sources():
    return render_template('tool_sources.html')

@app.route('/favicon.ico') 
def favicon(): 
    return send_from_directory(os.path.join(app.root_path, 'static'), 'favicon.ico')

if __name__ == '__main__':
    app.run(debug=True)


