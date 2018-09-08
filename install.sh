#!/bin/bash

#libbsc
rm -rf src/tools/libbsc
git clone https://github.com/shubhamchandak94/libbsc.git src/tools/libbsc
(cd src/tools/libbsc && make)
cp src/tools/libbsc/bsc ./

#p7zip
rm -rf src/tools/p7zip_16.02
wget https://sourceforge.net/projects/p7zip/files/p7zip/16.02/p7zip_16.02_x86_linux_bin.tar.bz2 -P src/tools/
tar -jxvf src/tools/p7zip_16.02_x86_linux_bin.tar.bz2 -C src/tools/
rm -rf src/tools/p7zip_16.02_x86_linux_bin.tar.bz2
cp src/tools/p7zip_16.02/bin/7z ./
cp src/tools/p7zip_16.02/bin/7z.so ./

# a pseudo file ‘config.h’ must be provided
cd src
make clean
cp config_pseudo.h config.h
make minicomsg
# echo -e "\033[40;37m Compress single reads is OK. \033[0m"

echo "#define ORDER" >> config.h
echo "int cmpcluster3(const void *a_, const void *b_);" >> config.h
make minicompe
# echo -e "\033[40;37m Compress paired-end reads is OK. \033[0m"

make decompress
# echo -e "\033[40;37m Decompress is OK. \033[0m"

cd ..
cp src/decompress ./

chmod +x minicom
