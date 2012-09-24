import mock
import zleveldb
import jeweler_pb2

def test_basic_file_op():
    filename = "test_data/test_python"
    fd = open(filename, "wb")
    bd = jeweler_pb2.BraceletData()
    bd.name = "TEST"
    zleveldb.write_protobuf_data(fd, bd)
    fd.close()
    fd = open(filename, "rb")
    bd = jeweler_pb2.BraceletData()
    zleveldb.load_protobuf_data(fd, bd)
    assert("TEST" == bd.name)


def test_loop_file_op():
    filename = "test_data/test_python"
    fd = open(filename, "wb")
    bd = jeweler_pb2.BraceletData()
    bd.name = "TEST"
    zleveldb.write_protobuf_data(fd, bd)
    fd.close()
    fd = open(filename, "rb")
    bd = jeweler_pb2.BraceletData()
    while zleveldb.load_protobuf_data(fd, bd):
        assert("TEST" == bd.name)
    fd.close()


def test_single_protobuf_data():
    filename = "test_data/test_single"
    bd = jeweler_pb2.BraceletData()
    bd.name = "TEST"
    zleveldb.write_single_protobuf_data(filename, bd)
    bd = jeweler_pb2.BraceletData()
    zleveldb.load_single_protobuf_data(filename, bd)
    assert("TEST" == bd.name)
