import constants
import sys

class ShopInfo:
    def __init__(self):
        self.cuffcompare_folder = sys.argv[1]
        self.jeweler_folder = sys.argv[2]
        self.bracelet_folder = sys.argv[3]
        self.mismatch_analyzer_folder = sys.argv[4]
        self.result_folder = sys.argv[5]
        if (not os.path.exists(self.result_folder)):
            os.makedirs(self.result_folder)

    @property
    def cuffcompare_file(self):
        return self.cuffcompare_folder + "/" + "cuffcompare.tracking"

    @property
    def mismatch_analyzer_file(self):
        return self.mismatch_analyzer_folder + "/" + "result.consistent.locations"

    @property
    def bracelet_file(self):
        return self.bracelet_folder + "/" + "result.bracelet"

    @property
    def gene_relationship_file(self):
        return self.result_folder + "/gene_relationships.obj"

    @property
    def gene_meta_file(self):
        return self.result_folder + "/gene_meta.obj"

    @property
    def blacklist_file(self):
        return self.result_folder + "/blacklist.obj"

    def load_num_reads(self, gene_id):
        ret = 0
        single = 0
        multiple = 0
        filename  = self.jeweler_folder + "/"+ gene_id + "/" + gene_id + ".mamf.multiple.reads"
        lines = open( filename).readlines()
        ret += len(lines)
        filename  = self.jeweler_folder + "/"+ gene_id +"/"+gene_id+ ".mamf.single.reads"
        lines = open( filename).readlines()
        ret += len(lines)
        return ret
