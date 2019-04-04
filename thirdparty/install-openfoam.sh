mkdir vtk-build
cd vtk-build

cmake ../vtk
 -DBUILD_SHARED_LIBS:BOOL=OFF \
 -DCMAKE_BUILD_TYPE=Release \
 -DVTK_USE_SYSTEM_ZLIB:BOOL=ON 

make -j4
make install
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib' >> ~/.bashrc
