import mock
import jeweler_pb2
import cuffcompare
import os
from nose.tools import *

class TestCuffcompareResult:
    def cleanup(self):
        if os.path.exists(self.proto_filename):
            os.remove(self.proto_filename)
        
    def setUp(self):
        self.filename = "test_data/cuffcompare.result"
        self.proto_filename = "test_data/cuffcompare.result.proto.data"
        self.cleanup()
        self.cuffcompare_result = cuffcompare.CuffcompareResult(self.filename)
        
    def tearDown(self):
        self.cleanup()

    def test_construct(self):
        assert self.proto_filename == self.cuffcompare_result.protobuf_file
        data = self.cuffcompare_result.data.gene_data
        assert 6 == len(data)
        assert_equal('CUFF.1044', data[0].gene_id)
        assert_equal('CUFF.1044.1', data[1].transcript_id)
        assert_equal('', data[2].gene_name)
        assert_equal('', data[2].transcript_name)
        assert_equal(False, data[2].is_good)
        assert_equal('Fam72a', data[3].gene_name)
        assert_equal('ENSMUST00000068613', data[3].transcript_name)
        assert_equal(True, data[3].is_good)
        assert_equal(False, data[3].is_pseudo)
        assert_equal(True, data[5].is_pseudo)
        assert_equal(None, self.cuffcompare_result.get('CUFF.10123'))
        assert_not_equal(None, self.cuffcompare_result.get('CUFF.1044'))
