#!/bin/bash
module load cmake

pushd thirdparty
wget https://github.com/fmtlib/fmt/archive/refs/tags/9.1.0.tar.gz
tar -zxf ./9.1.0.tar.gz
cd fmt-9.1.0
mkdir build
cd build
cmake -DBUILD_SHARED_LIBS=TRUE -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE ..
make -j4
# Install to <project_directory>/thirdparty/usr/...
make install DESTDIR=../../

popd
pushd thirdparty
wget https://github.com/leethomason/tinyxml2/archive/refs/tags/9.0.0.tar.gz
tar -zxf ./9.0.0.tar.gz
cd tinyxml2-9.0.0/
mkdir build
cd build
cmake ..
make -j4
# Install to <project_directory>/thirdparty/usr/...
make install DESTDIR=../../

popd
pushd thirdparty
# Install libuuid
yumdownloader --resolve libuuid-devel
rpm2cpio libuuid-devel-2.23.2-65.el7_9.1.x86_64.rpm | cpio -idmv
# The project depends on UUID_STR_LEN which is not defined in the centos7 release. So let's add that now.
patch -p0 < 0001_uuid_length.patch

# Install libmd
wget https://dl.fedoraproject.org/pub/epel/7/x86_64/Packages/l/libmd-devel-1.0.4-2.el7.x86_64.rpm
rpm2cpio libmd-devel-1.0.4-2.el7.x86_64.rpm | cpio -idmv  
wget https://dl.fedoraproject.org/pub/epel/7/x86_64/Packages/l/libmd-1.0.4-2.el7.x86_64.rpm
rpm2cpio libmd-1.0.4-2.el7.x86_64.rpm | cpio -idmv  
