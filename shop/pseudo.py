'''
This object is used for pseudo list
'''

class PseudoSet:

    def __init__(self):
        self.filename = "info/pseudo_gene_names"
        self.data =set()
        self.get_pseudogene_list()

    def get_pseudogene_list(self):
        lines = open(self.filename).readlines()
        for line in lines:
            data = line.strip()
            self.data.add(data)

    def is_pseudo_gene(self, name):
        return name in self.data
