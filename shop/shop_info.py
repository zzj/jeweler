import sys
import os

class ShopInfo:
    def __init__(self):
        self.cuffcompare_folder = sys.argv[1]
        self.jeweler_folder = sys.argv[2]
        self.bracelet_folder = sys.argv[3]
        self.mismatch_analyzer_folder = sys.argv[4]
        self.result_folder = sys.argv[5]
        self.cufflinks_folder = sys.argv[6]
        self.sample_id = sys.argv[7]
        try:
            if sys.argv[8] == "--simulation":
                self.is_simulation = True
                self.simulation_profile = sys.args[9]
        except Exception as e:
            self.is_simulation = False
        if (not os.path.exists(self.result_folder)):
            os.makedirs(self.result_folder)

    @property
    def cuffcompare_file(self):
        return self.cuffcompare_folder + "/" + "cuffcompare.tracking"

    @property
    def mismatch_analyzer_file(self):
        return self.mismatch_analyzer_folder + "/" + "result"

    @property
    def bracelet_file(self):
        return self.bracelet_folder + "/" + "result.bracelet"

    @property
    def test_data_file(self):
        return self.result_folder + "/test_data.obj"

    @property
    def blacklist_file(self):
        return self.result_folder + "/blacklist.obj"

    @property
    def jeweler_db(self):
        return self.jeweler_folder