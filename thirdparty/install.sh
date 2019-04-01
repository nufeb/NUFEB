mkdir vtk-build
cd vtk-build

cmake ../VTK-6.3.0
 -DBUILD_SHARED_LIBS:BOOL=OFF \
 -DCMAKE_BUILD_TYPE=Release \
 -DVTK_USE_SYSTEM_ZLIB:BOOL=ON 

make -j4
make install
