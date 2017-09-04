
## Library compilation

Create a __build__ folder (could be here, for instance), create a folder to 
compile the library and proceed to use CMake there, i.e., type in a terminal:

    mkdir build
    cd build
    mkdir lib
    cd lib
    cmake <PATH_TO_LIBRARY>
    make
    make install
    
If the build folder is located here, then ```<PATH_TO_LIBRARY>``` is simply 
`../../Library`.



## Examples compilation

**Currently there is only a test inside the folder examples, but I plan to write 
examples as soon as the implementation is done.**

At the moment I am not using global variables in the CMake file and the following 
structure is assumed:

``` .
    |--build
       |-- lib
       |-- ex
```

Go to your existing __build__ folder (same as you used for the Library, see 
directory tree above), create a folder to compile the examples and proceed to 
use CMake there, i.e., type in a terminal:

    cd build
    mkdir ex
    cd ex
    cmake <PATH_TO_EXAMPLES>
    make
    
If the build folder is located here, then ```<PATH_TO_EXAMPLES>``` is simply 
`../../Examples`.

Note that you can create the build folder for Examples in a different location, 
but then you will need to change the CMake file in order to link it to the 
library.

