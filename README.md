# RW/CSE MSc Course "Advanced Numerical Methods for CSE" -- Code Repository

This repository will host all the codes used in the lecture notes and assignments.

The [lecture homepage](https://moodle-app2.let.ethz.ch/course/view.php?id=3643) contains more information regarding the course, 
and will be used to distribute lecture notes and assignment handouts.

__Note:__ For submitting your code solutions, please use the online submission interface in the [lecture homepage](https://moodle-app2.let.ethz.ch/course/view.php?id=3643).

If you are new to gitLab, you can use [this tutorial](https://gitlab.math.ethz.ch/tille/gitlab-introduction) on Git and GitLab 
to make yourself familiar with the interface.

**Additional links**

- [Lecture homepage](https://moodle-app2.let.ethz.ch/course/view.php?id=4245)
- [Course VVZ](http://www.vvz.ethz.ch/Vorlesungsverzeichnis/lerneinheit.view?lerneinheitId=117918&semkez=2017W&ansicht=KATALOGDATEN&lang=en)
- [Git tutorial](https://gitlab.math.ethz.ch/tille/gitlab-introduction/blob/master/git/README.md)
- [Eigen documentation](http://eigen.tuxfamily.org/dox/)
- [Debugger tutorial](https://gitlab.math.ethz.ch/tille/debugging-cpp-code-with-lldb)

## Download and Compile Everything

        $ git clone --recurse-submodules https://gitlab.math.ethz.ch/AdvNumCSE/Code
        $ cd Code
        $ mkdir build && cd build
        $ cmake ..
        $ make

All binaries can be found at the `build/bin` path.

## Download and Compile One Code

The path for each problem / lecture code is:

 - `<Chapter>/<ProblemName>/templates_nolabels` or `<Chapter>/<LectureCode>/templates_nolabels` for templates.
 - `<Chapter>/<ProblemName>/solutions_nolabels` or `<Chapter>/<LectureCode>/solutions_nolabels` for solutions.

Each of these folders contains an independent `CMake` file.
You should be able to compile and execute a code using:

```
$ cmake .
$ make
$ cd bin
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

        $ ./build/bin/<chapter>/<problem-name>/

**TIP**

Using:

        $ make -j<number_of_processes>
    
may (or may not) speed up the compilation time.

## Extra

### Dependencies / Requirements

**Required**

- C++ compiler (C++11 support required): only tested `gcc` and `clang`
- Git
- CMake
- Eigen (automatically downloaded when running `cmake` for the first time)

__Debian__ (e.g., Ubuntu, Linux Mint)

        sudo apt-get install build-essential git cmake

__RedHat__ (e.g., Fedora, CentOS)

        sudo yum install make gcc gcc-c++ git cmake

__Mac OS X__

Dependencies can be installed with [brew](http://brew.sh/):

        brew tap homebrew/science
        brew install git cmake
	
__Windows__

If you are using Windows, it is important to note that the exam will be in Linux.
Therefore, it is a good idea to familiarize yourself with the Linux tools.

We recommend you to use a virtual machine through [VirtualBox](https://www.virtualbox.org/).
You can install [LUbuntu](http://lubuntu.net/) in VirtualBox,
as it is relatively light on resources.

### Repository Policy

__All relevant material for students can be found in the `master` branch.__

The assistants will temporary open other branches to develop particular features or assignments.
These branches will be merged into the `master` branch when their contents are ready for the students.

We kindly ask students to please __ignore the additional branches__ and to only __keep track of the 'master' branch__.

