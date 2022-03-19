import pymongo

__author__ = 'shatabdi'


class Database(object):
    URI = "mongodb://127.0.0.1:27017"
    DATABASE = None

    @staticmethod
    def initialize():
        client = pymongo.MongoClient(Database.URI)
        # Database.DATABASE = client['bigdb']
        Database.DATABASE = client['bigfeatureStore']

    @staticmethod
    def find(collection, query):
        return Database.DATABASE[collection].find(query)

    @staticmethod
    def find_one(collection, query):
        return Database.DATABASE[collection].find_one(query)

    @staticmethod
    def get_collection_size():
        return Database.DATABASE.dataSize()

    @staticmethod
    def remove(collection, query):
        return Database.DATABASE[collection].remove(query)

    @staticmethod
    def find_by_structure(collection, structureName):
        return Database.DATABASE[collection].find({}, structureName)
        # return Database.DATABASE[collection].find({}, structureName).limit(10000)

    @staticmethod
    def find_all(collection,dict):
        return Database.DATABASE[collection].find({},dict)

    # @staticmethod
    # def left_label():
    #     res = Database.DATABASE['di_nuc'].aggregate([
    #     {
    #
    #         "$lookup":{
    #                 "from": "label",
    #                 "localField": "ID",
    #                 "foreignField": "ID",
    #                 "as": "common"
    #         }
    #     },{
    #         "$project": {'ID':1,'kmer_2_AA': 1, 'kmer_2_AC': 1, 'kmer_2_AG': 1, 'kmer_2_AT': 1,'classical_label':1,'core_label':1,
    #                      'common':1, '_id':0}
    #         }
    #     ])
    #     return res


