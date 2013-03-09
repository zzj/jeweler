
wget http://csbio.unc.edu/gene_scissors/learning_model_all

wget http://downloads.sourceforge.net/project/boost/boost/1.52.0/boost_1_52_0.tar.gz
tar zxf boost_1_52_0.tar.gz -C lib/
cd lib/boost_1_52_0
./bootstrap.sh
./b2 install --prefix=`pwd`/../../third_party/bin/ -j8

cd ../..
wget http://leveldb.googlecode.com/files/leveldb-1.9.0.tar.gz
tar zxf leveldb-1.9.0.tar.gz -C lib/
cd lib/leveldb-1.9.0
make -j8

cd ../glog-0.3.2/
./configure --prefix=`pwd`/../../third_party/bin/
make -j8

cd ../protobuf-2.4.1
./configure --prefix=`pwd`/../../third_party/bin/
make -j8
make install
cd ..

cd bamtools/
mkdir build
cd build
cmake ..
make
cd ../../..
