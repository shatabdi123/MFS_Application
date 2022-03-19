from mongoengine import *
import json

__author__ = 'shatabdi'

# connect(host="mongodb://127.0.0.1:27017/featureStore")


class features(Document):
    ID =  StringField(required=True,unique=True)
    Gene_length = IntField(null=True)
    Exon_number = IntField(null=True)
    p3UTR_length = IntField(null=True)
    p5UTR_length = IntField(null=True)
    CDS_length = IntField(null=True)
    ContentA = FloatField(null=True)
    ContentG = FloatField(null=True)
    ContentGC = FloatField(null=True)
    CAI = FloatField(null=True)

    @staticmethod
    def initialize():
        connect(host="mongodb://127.0.0.1:27017/featureStore")

    @staticmethod
    def create_structures_fromFile(filename):

        with open(filename) as file:
            file_data = json.load(file)

        if isinstance(file_data, list):

            feature_instances = [features(**data) for data in file_data]

            features.objects.insert(feature_instances, load_bulk=False)



    @staticmethod
    def insert_one(data):
        feature_instances = features(**data)
        features.objects.insert(feature_instances)

    @staticmethod
    def update(query, data):
        return features.objects.update(query, data, upsert=True)


