#!/bin/bash
OPTION="-lboost_serialization -lboost_iostreams -lboost_program_options -lboost_system -lboost_filesystem"
g++ hash.cpp -o hash.o ${OPTION} 
g++ similarity.cpp -o similarity.o ${OPTION} 

# anarci symlink
ln -sf anarci-1.3/lib/python/anarci anarci

# test
bash test.sh
