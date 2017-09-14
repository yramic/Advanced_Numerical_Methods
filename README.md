# RW/CSE MSc Course "Advanced Numerical Methods for CSE" -- Code Repository

## Download and Compile Everything

	$ git clone https://gitlab.math.ethz.ch/AdvNumCSE/Code
	$ cd Code
	$ mkdir build && cd build
	$ cmake ..
	$ make

All binaries can be found in the `build` folder.

## Download and Compile One Code

The path for each problem / lecture code is:

 - `<Chapter>/<ProblemName>/templates_nolabels` or `<Chapter>/<LectureCode>/templates_nolabels` for templates.
 - `<Chapter>/<ProblemName>/solutions_nolabels` or `<Chapter>/<LectureCode>/solutions_nolabels` for solutions.

Each of these folders contains an independent `CMake` file. Either within the subfolder of the cloned repository or using the `Download zip` button, you should be able to compile and execute the codes using:

```
$ cmake .
$ make
$ ./executable_name
```

Alternatively, after running `cmake ..` in the `build` folder, you can recover the list of all possible targets of `make` by typing:

    $ cmake --build . --target help

or

    $ make help

Then you can choose a target and run:

    $ make <target>

Targets will be in the format `<chapter>_<problem_name>_<code_type>`, where `<code_type>` is either `template` or `solution`.

The corresponding executable will be located in:

        $ ./build/<chapter>/<problem-name>/

**TIP**

Using:

    $ make -j<number_ob_processes>
    
may (or may not) speed up the compilation time.

## Extra

### Third-Party Libraries

**Required**

- C++ compiler (C++11 support required): only tested `gcc` and `clang`
- Git
- CMake
- Eigen (automatically downloaded when running `cmake` for the first time)

### F.A.Q.

Some package may be missing on your machine.

- On a fresh install of *Ubuntu*:

        sudo apt-get install git cmake
    
- On *Mac OS X*:

If you are missing CMake on Mac OS X, the easiest way to obtain those packages is via [Homebrew](http://brew.sh/). After Homebrew has been installed, you can install CMake:

    brew install cmake

### Additional links

- [Course webpage](https://moodle-app2.let.ethz.ch/course/view.php?id=3643)
- [Course VVZ](http://www.vvz.ethz.ch/Vorlesungsverzeichnis/lerneinheitPre.do?semkez=2017W&ansicht=KATALOGDATEN&lerneinheitId=117918&lang=en)
- [Git tutorial](https://gitlab.math.ethz.ch/tille/gitlab-introduction/blob/master/git/README.md)
- [Debugger tutorial](https://gitlab.math.ethz.ch/tille/debugging-cpp-code-with-lldb)

