mkdir $HOME/bin/
wget http://downloads.sourceforge.net/project/boost/boost/1.52.0/boost_1_52_0.tar.gz
tar zxf boost_1_52_0.tar.gz -C lib/
cd lib/boost_1_52_0
./bootstrap.sh
./b2 install --prefix=$HOME/bin/

cd ../leveldb-1.5.0
make

cd ../glog-0.3.2/
./configure --prefix=$HOME/bin/
make

cd ../protobuf-2.4.1
./configure --prefix=$HOME/bin/
make
make install

cd ../..
