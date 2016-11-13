echo "Removing old build"
rm -rf CMakeFiles cmake_install.cmake CMakeCache.txt CodedAperture Makefile

echo "Building..."
cmake3 -DGeant4_DIR=/Users/Hellfeld/GEANT4.10.2.p1/geant4-install/lib/Geant4-10.2.1/ .

