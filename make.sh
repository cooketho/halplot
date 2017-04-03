# Change PGPATH to match the location of the progressiveCactus installation on your machine.
PGPATH=/usr/local/bin/progressiveCactus

CPATH=${PGPATH}/submodules/hal/lib/
CPATH=${PGPATH}/submodules/sonLib/lib/:$CPATH
CPATH=/usr/local/lib/:$CPATH
export CPATH
LIBRARY_PATH=${PGPATH}/submodules/hal/lib
LIBRARY_PATH=${PGPATH}/submodules/hdf5/lib/:$PATH
LIBRARY_PATH=${PGPATH}/submodules/sonLib/lib/:$PATH
export LIBRARY_PATH
g++ -std=c++0x -c -o haltraverse.o haltraverse.cpp
g++ -std=c++0x haltraverse.o -l:halLib.a -L${PGPATH}/submodules/hdf5/lib -lhdf5 -lhdf5_cpp -l:sonLib.a -o haltraverse
