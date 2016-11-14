echo "\nRemoving old build"
rm -rf CMakeFiles cmake_install.cmake CMakeCache.txt PRISM Makefile

echo "\nBuilding...\n"
cmake3 -DGeant4_DIR=/Users/Hellfeld/GEANT4.10.2.p1/geant4-install/lib/Geant4-10.2.1/ .

echo "\nDone! Now run $ make to compile.\n"
