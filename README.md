# RW/CSE MSc Course "Advanced Numerical Methods for CSE" (R. Hiptmair) -- Code Repository

This repository will host all the codes used in the lecture notes and assignments.

The [lecture homepage](https://metaphor.ethz.ch/x/2023/hs/401-4671-00L/) and the
introduction to the [lecture
document](https://people.math.ethz.ch/~grsam/ADVNCSE/ADVNCSE.pdf) provide more information
regarding the course. 

In particular, descriptions for all homework projects are contained in the [AdvNumCSE
Homework
Collection](https://people.math.ethz.ch/~grsam/ADVNCSE/HOMEWORK/ADVNCSEProblems.pdf). The
weekly homework assignments will be published on the [lecture homepage](https://metaphor.ethz.ch/x/2023/hs/401-4671-00L/).

If you are new to gitLab, you can use [this tutorial](https://gitlab.math.ethz.ch/tille/gitlab-introduction) on Git and GitLab 
to make yourself familiar with the interface.

**Additional links**

- [Lecture homepage]()https://metaphor.ethz.ch/x/2023/hs/401-4671-00L/
- [Course VVZ entry](https://www.vorlesungen.ethz.ch/Vorlesungsverzeichnis/lerneinheit.view?semkez=2023W&ansicht=ALLE&lerneinheitId=174538&lang=de)
- [Git tutorial](https://gitlab.math.ethz.ch/tille/gitlab-introduction/blob/master/git/README.md)
- [Eigen documentation](http://eigen.tuxfamily.org/dox/)
- [Debugger tutorial](https://gitlab.math.ethz.ch/tille/debugging-cpp-code-with-lldb)

## Download and Compile Everything

        $ git clone --recurse-submodules https://gitlab.math.ethz.ch/AdvNumCSE/Code
        $ cd Code
        $ mkdir build && cd build
        $ cmake ..
        $ make

All binaries can be found inside the `build` directory in respective subdirectories. 

## Download and Compile a Single Code

The path for each problem / lecture code is:

 - `<Chapter>/<ProblemName>/templates_nolabels` or `<Chapter>/<LectureCode>/templates_nolabels` for templates.
 - `<Chapter>/<ProblemName>/solutions_nolabels` or `<Chapter>/<LectureCode>/solutions_nolabels` for solutions.

Each of these codes correspond to a `target`, which is a `make` command that can be used to compile only that one code and the related dependencies.
Targets are in the format `a<assignment_number>problem<problem_number>` (e.g., `a1problem1`) or `a<assignment_number>problem<problem_number>_sol`, depending on whether the code is a template or a solution.

After running `cmake ..` in the `build` folder, you can recover the list of all possible targets by typing:

        $ cmake --build . --target help

or

        $ make help

Then you can choose a target and run:

        $ make <target>

The corresponding executable will be located in:

        $ ./build/bin/

**TIP**

Using:

        $ make -j<number_of_processes>
    
may (or may not) speed up the compilation time.


## Homework Projects Code Upload

After cloning the AdvNumCSE code repository you are supposed to set up _your own working
branch_ named after **your NETHZ login ID**

	$ cd <repository root directory>
	$ git branch <NETHZ ID>
	$ git checkout <NETHZ ID>

After you have edited files, created new ones (tell Git aboout them with `git add`), have
written and tested your codes, `git push` your project codes to the server. 

Your solution codes must reside in `<Chapter>/<ProblemName>/mysolution` and it should be
possble to generate an executable from them through issuing a suitable `make` command. 

Upload of your homework is important, because the **code review** will be based on it:
during the code review you will be shown some of your C++ source codes and you will be asked to
explain them. 

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

        brew tap brewsci/science
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

