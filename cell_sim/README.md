# CAPMD

## Build instructions

Dependencies:

* [CMake](https://cmake.org)
* [VTK](https://vtk.org)
* [PROJ](https://proj.org)
* [PDAL](https://pdal.io)
* [unixODBC](http://www.unixodbc.org)

```
git clone https://github.com/nrstillman/CAPMD.git
cd CAPMD
mkdir build
cd build
cmake ..
make
```


## Use instructions

There are two options to use capmd, as a c++ program or as a python module. 

For c++ program, build as capmd-debug:
	
	`cmake --build ./capmd/cmake-build-debug --target capmd-debug -- `

The current working version can then be run using main.cpp. 


For python module, build c++ as pycapmd:

	`cmake --build ./capmd/cmake-build-debug --target pycapmd -- `

Following succesful build, a python module (pycapmd.so) will be created. Move this to `python-run` and run using either `testing_pycamd.py` or `testing_pycamd.ipynb`. 

## Visualisation instructions

All visualisation occurs using Paraview. To view animation of simulation output, follow these steps:

* Open Paraview and load files from simulation output
* Select: Properties - Apply on loaded files. Then select Features-Common-Glyph
* This will create a new object called Glyph. Under properties of this object, select Glyph source - Sphere, Radius - Scalars, Scale mode - Scalars, Scale factor - 2, Glyph mode - all points. 
* Press Apply. Particles should now be visible in the rendering window. 




