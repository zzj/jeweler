import struct

def load_protobuf_data(fd, data):
    k = fd.read(4)
    if (len(k) != 4):
        return False
    m = struct.unpack("i", k)[0]
    buf = fd.read(m)
    if (len(buf) != m):
        return False
    data.ParseFromString(buf)
    return True

def write_protobuf_data(fd, data):
    buf = data.SerializeToString()
    fd.write(struct.pack("i",len(buf)))
    fd.write(buf)


def load_single_protobuf_data(filename, data):
    fd = open(filename, "rb")
    data.ParseFromString(fd.read())
    fd.close()

def write_single_protobuf_data(filename, data):
    fd = open(filename, "wb")
    fd.write(data.SerializeToString())
    fd.close()
