import pandas as pd
import numpy as np
import json
from pandas import DataFrame
from bson.json_util import dumps
from src.common.database import Database
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


class DNA():

    @staticmethod
    def structure(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('struc_dis', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')
                for column in struc_merged.columns:
                    if column == 'distance':
                        struc_merged[column] = struc_merged[column].abs()
                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def codon(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def ProteinStructure(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and  "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                if 'SignalP' in struc_merged.columns:
                    struc_merged['SignalP'] = struc_merged['SignalP'].replace(1, 'SignalP')
                    struc_merged['SignalP'] = struc_merged['SignalP'].replace(0, 'No SignalP')
                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def ProteinLocalization(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def GeneExp(select, analysis):
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
        dict['ID'] = 1
        new_dict = {}
        new_dict.update(take(5, dict.items()))
        new_dict['ID'] = 1
        new_dict['_id'] = 0
        label_dict = {}
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('gene_exp', new_dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def ProteinExp(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('prot_exp', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def ChromStates(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def ATACseq(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def HistoneModifications(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def DNAmethylation(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def miRNA(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def PFAM(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def Insertions(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def TFbindingSite(select, analysis):
        dict = {}
        if "ARF" in str(select):
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/ARF_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif "BZIP" in str(select):
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/BZIP_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif "EREB" in str(select):
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/EREB_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif "LBD" in str(select):
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/LBD_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        elif "SBP" in str(select):
            col_file = open("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/col_names/SBP_name.txt", "r")
            col_name = []
            col_name.append('ID')
            for line in col_file:
                col_name.append(line.strip())
                dict = Convert(col_name)
        else:
            dict = Convert(select)
        label_dict = {}
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def TSS(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def Enhancers(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def TE(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def Quadruplex(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('count', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def Correlation(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def Varionomic(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('varionomic', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def Other(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['_id'] = 0
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            cursor = Database.find_by_structure('many_col_data', dict)
            list_cur = list(cursor)
            struc_df = DataFrame(list_cur)
            struc_df.drop_duplicates(subset="ID", keep=False, inplace=True)
            print(struc_df.head())

            if "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                struc_merged = pd.merge(left=struc_df, right=label_df, how='left', left_on='ID', right_on='ID')

                records = json.loads(struc_merged.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)
            else:
                records = json.loads(struc_df.T.to_json()).values()
                list_cur = list(records)
                json_data = dumps(list_cur, indent=2)

        return json_data

    @staticmethod
    def DNA_seq(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            if 'dinucleic' in str(select) and len((select)) == 1:
                del dict['dinucleic']
                print(dict)
                dict1 = {'kmer_2_AA': 1, 'kmer_2_AC': 1, 'kmer_2_AG': 1, 'kmer_2_AT': 1}
                dict.update(dict1)
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('di_nuc', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:
                print("Hi")
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                if 'dinucleic' in str(select):

                    dict1 = {'ID': 1, 'kmer_2_AA': 1, 'kmer_2_AC': 1, 'kmer_2_AG': 1, 'kmer_2_AT': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('di_nuc', dict1)
                    list_cur = list(cursor1)
                    di_df = DataFrame(list_cur)
                    di_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(di_df.head())
                    di_merged = pd.merge(left=di_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    # di_merged.drop_duplicates(subset="ID",
                    #                      keep=False, inplace=True)
                    #       di_merged.to_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/my.txt", sep='\t',
                    # index=False, header=True)
                    records = json.loads(di_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'trinucleic' in str(select):
                    dict1 = {'ID': 1, 'kmer_3_AAA': 1, 'kmer_3_AAC': 1, 'kmer_3_AAG': 1, 'kmer_3_AAT': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('tri_nuc', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'DCC' in str(select):
                    cols = []
                    for name in col_name:
                        if name.startswith("DCC"):
                            cols.append(name)
                            dict1 = Convert(cols)
                    dict1['_id'] = 0
                    dict1['ID'] = 1
                    print(dict1)
                    cursor1 = Database.find_by_structure('auto_nuc', dict1)
                    list_cur = list(cursor1)
                    nuc_DCC_df = DataFrame(list_cur)
                    nuc_DCC_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(nuc_DCC_df.head())
                    nuc_DCC_merged = pd.merge(left=nuc_DCC_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(nuc_DCC_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'DACC' in str(select):
                    cols = []
                    for name in col_name:
                        if name.startswith("DACC"):
                            cols.append(name)
                            dict1 = Convert(cols)
                    dict1['_id'] = 0
                    dict1['ID'] = 1
                    print(dict1)
                    cursor1 = Database.find_by_structure('auto_nuc', dict1)
                    list_cur = list(cursor1)
                    nuc_DACC_df = DataFrame(list_cur)
                    nuc_DACC_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(nuc_DACC_df.head())
                    nuc_DACC_merged = pd.merge(left=nuc_DACC_df, right=label_df, how='left', left_on='ID',
                                               right_on='ID')
                    records = json.loads(nuc_DACC_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'TCC' in str(select):
                    cols = []
                    for name in col_name:
                        if name.startswith("TCC"):
                            cols.append(name)
                            dict1 = Convert(cols)
                    dict1['_id'] = 0
                    dict1['ID'] = 1
                    print(dict1)
                    cursor1 = Database.find_by_structure('auto_nuc', dict1)
                    list_cur = list(cursor1)
                    nuc_TCC_df = DataFrame(list_cur)
                    nuc_TCC_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(nuc_TCC_df.head())
                    nuc_TCC_merged = pd.merge(left=nuc_TCC_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(nuc_TCC_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'TACC' in str(select):
                    cols = []
                    for name in col_name:
                        if name.startswith("TACC"):
                            cols.append(name)
                            dict1 = Convert(cols)
                    dict1['_id'] = 0
                    dict1['ID'] = 1
                    print(dict1)
                    cursor1 = Database.find_by_structure('auto_nuc', dict1)
                    list_cur = list(cursor1)
                    nuc_TACC_df = DataFrame(list_cur)
                    nuc_TACC_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(nuc_TACC_df.head())
                    nuc_TACC_merged = pd.merge(left=nuc_TACC_df, right=label_df, how='left', left_on='ID',
                                               right_on='ID')
                    records = json.loads(nuc_TACC_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'PseDNC' in str(select):
                    cols = []
                    for name in col_name:
                        if name.startswith("PseDNC"):
                            cols.append(name)
                            dict1 = Convert(cols)
                    dict1['_id'] = 0
                    dict1['ID'] = 1
                    print(dict1)
                    cursor1 = Database.find_by_structure('auto_nuc', dict1)
                    list_cur = list(cursor1)
                    nuc_PseDNC_df = DataFrame(list_cur)
                    nuc_PseDNC_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(nuc_PseDNC_df.head())
                    nuc_PseDNC_merged = pd.merge(left=nuc_PseDNC_df, right=label_df, how='left', left_on='ID',
                                                 right_on='ID')
                    records = json.loads(nuc_PseDNC_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'PseKNC' in str(select):
                    dict1 = {'ID': 1, 'PseKNC_3_Xc1_AAA': 1, 'PseKNC_3_Xc1_AAC': 1, 'PseKNC_3_Xc1_AAG': 1, 'PseKNC_3_Xc1_AAT': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('auto_nuc', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
        return json_data

    @staticmethod
    def Protein_seq(select, analysis):
        label_dict = {}
        dict = Convert(select)
        dict['ID'] = 1
        print(dict)
        if len((analysis)) == 0:
            if 'AAC' in str(select) and len((select)) == 1:
                del dict['AAC']
                print(dict)
                dict1 = {'AAC_A': 1, 'AAC_R': 1, 'AAC_N': 1, 'AAC_D': 1}
                dict.update(dict1)
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('mono_pep', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif "classical_label" in str(select) and "core_label" in str(select) and "Origin" in str(select) and len((select)) == 3:
                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)
            elif ("classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select)) and len((select)) == 1:

                dict['_id'] = 0
                print(dict)
                cursor = Database.find_by_structure('label', dict)
                list_cur = list(cursor)
                json_data = dumps(list_cur, indent=2)

            elif "classical_label" in str(select) or "core_label" in str(select) or "Origin" in str(select):
                label_dict['_id'] = 0
                label_dict['ID'] = 1
                label_dict['classical_label'] = 1
                label_dict['core_label'] = 1
                label_dict['Origin'] = 1
                print(label_dict)
                cursor_label = Database.find_by_structure('label', label_dict)
                list_cur = list(cursor_label)
                # print(list_cur)
                label_df = DataFrame(list_cur)
                label_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                print(label_df)

                if "classical_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'classical_label']]
                elif "core_label" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'core_label']]
                elif "Origin" in str(select) and len((select)) == 2:
                    label_df = label_df[['ID', 'Origin']]
                if 'AAC' in str(select):
                    dict1 = {'ID': 1,'AAC_A': 1, 'AAC_R': 1, 'AAC_N': 1, 'AAC_D': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('mono_pep', dict1)
                    list_cur = list(cursor1)
                    di_df = DataFrame(list_cur)
                    di_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(di_df.head())
                    di_merged = pd.merge(left=di_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(di_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'DC' in str(select):

                    dict1 = {'ID': 1, 'DC_AA': 1, 'DC_RA': 1, 'DC_NA': 1, 'DC_DA': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('di_pep', dict1)
                    list_cur = list(cursor1)
                    di_df = DataFrame(list_cur)
                    di_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(di_df.head())
                    di_merged = pd.merge(left=di_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    # di_merged.drop_duplicates(subset="ID",
                    #                      keep=False, inplace=True)
                    #       di_merged.to_csv("C:/Users/Shatabdi/Documents/B73_V5_data/final_data_for_Db/my.txt", sep='\t',
                    # index=False, header=True)
                    records = json.loads(di_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'TC' in str(select):
                    dict1 = {'ID': 1, 'TC_AAA': 1, 'TC_RAA': 1, 'TC_NAA': 1, 'TC_DAA': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('tri_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'Moreau' in str(select):
                    dict1 = {'ID': 1, 'Moreau_CIDH920105_lag1': 1, 'Moreau_CIDH920105_lag2': 1, 'Moreau_CIDH920105_lag3': 1, 'Moreau_CIDH920105_lag4': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('Moreau_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'Moran' in str(select):
                    dict1 = {'ID': 1, 'Moran_CIDH920105_lag1': 1, 'Moran_CIDH920105_lag2': 1,
                             'Moran_CIDH920105_lag3': 1, 'Moran_CIDH920105_lag4': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('Moran_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'Geary' in str(select):
                    dict1 = {'ID': 1, 'Geary_CIDH920105_lag1': 1, 'Geary_CIDH920105_lag2': 1,
                             'Geary_CIDH920105_lag3': 1, 'Geary_CIDH920105_lag4': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('Geary_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'Composition' in str(select):
                    dict1 = {'ID': 1, 'CTDC_hydrophobicity_Group1': 1, 'CTDC_normwaalsvolume_Group1': 1,
                             'CTDC_polarity_Group1': 1, 'CTDC_polarizability_Group1': 1,'CTDC_charge_Group1':1,
                             'CTDC_secondarystruct_Group1':1,'CTDC_solventaccess_Group1':1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('CTD_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
                if 'Transition' in str(select):
                    dict1 = {'ID': 1, 'CTDT_prop1_Tr1221': 1, 'CTDT_prop2_Tr1221': 1,
                             'CTDT_prop3_Tr1221': 1, 'CTDT_prop4_Tr1221': 1,'CTDT_prop5_Tr1221':1,
                             'CTDT_prop6_Tr1221':1,'CTDT_prop7_Tr1221':1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('CTD_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'Distribution' in str(select):
                    dict1 = {'ID': 1, 'CTDD_prop1_G1_residue0': 1, 'CTDD_prop2_G1_residue0': 1,
                             'CTDD_prop3_G1_residue0': 1, 'CTDD_prop4_G1_residue0': 1,'CTDD_prop7_G1_residue0':1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('CTD_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'CTriad' in str(select):
                    dict1 = {'ID': 1, 'CTriad_VS111': 1, 'CTriad_VS211': 1,
                             'CTriad_VS311': 1, 'CTriad_VS411': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('CTriad_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'SOCN' in str(select):
                    dict1 = {'ID': 1, 'SOCN_Schneider_lag1': 1, 'SOCN_Schneider_lag2': 1, 'SOCN_Grantham_lag1': 1,
                             'SOCN_Grantham_lag2': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('Quasi_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'QSO' in str(select):
                    dict1 = {'ID': 1, 'QSO_Schneider_Xr_A': 1, 'QSO_Schneider_Xr_R': 1, 'QSO_Schneider_Xr_N': 1,
                             'QSO_Schneider_Xr_D': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('Quasi_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'PAAC' in str(select):
                    dict1 = {'ID': 1, 'PAAC_Xc1_A': 1, 'PAAC_Xc1_R': 1, 'PAAC_Xc1_N': 1,
                             'QSO_Schneider_Xr_D': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('Pseudo_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)

                if 'QSO' in str(select):
                    dict1 = {'ID': 1, 'APAAC_Pc2_Hydrophobicity_2': 1, 'APAAC_Pc2_Hydrophilicity_2': 1, 'APAAC_Pc2_Hydrophobicity_3': 1,
                             'APAAC_Pc2_Hydrophilicity_3': 1}
                    dict1['_id'] = 0
                    print(dict1)
                    cursor1 = Database.find_by_structure('Pseudo_pep', dict1)
                    list_cur = list(cursor1)
                    tri_df = DataFrame(list_cur)
                    tri_df.drop_duplicates(subset="ID", keep=False, inplace=True)
                    print(tri_df.head())
                    tri_merged = pd.merge(left=tri_df, right=label_df, how='left', left_on='ID', right_on='ID')
                    records = json.loads(tri_merged.T.to_json()).values()
                    list_cur = list(records)
                    json_data = dumps(list_cur, indent=2)
        return json_data