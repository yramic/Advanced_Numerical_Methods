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

## Download and Compile

        $ git clone https://gitlab.math.ethz.ch/AdvNumCSE/Code
        $ cd Code
        $ mkdir build && cd build
        $ cmake ..

To compile the desired Assignment you should change to the corresponding path and compile it.

    $ cd <Chapter>/<ProblemName>
    $ make

Will create multiple executables.
```
Code
├── build (was created by you)
│   ├── <Chapter>
│   :   ├── <ProblemName>
│       :   ├── <ProblemName>_mastersolution      (executable)
│           ├── <ProblemName>_mysolution          (executable)
│           ├── <ProblemName>_test_mastersolution (executable)
│           ├── <ProblemName>_test_mysolution     (executable)
│           :
│
├── <Chapter>
:   ├── <ProblemName>
    :   ├── all            (developer folder which creates the other folders)
        ├── mastersolution (folder containing mastersolution)
        ├── mysolution     (folder containing source files, to be modified by you)
        ├── templates      (folder containing source files)
        :
```
    
    $ ./<ProblemName>_mysolution 
Will execute the `main()` function found inside `ProblemName_main.cpp`

    $ ./<ProblemName>_test_mysolution 
Will run a series of tests found inside `mysolution/test/ProblemName_test.cpp`

By replacing `mysolution` with `mastersolution` in the later steps you can execute the provided mastersolution.

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

### Explanations for  Developers

When creating a new Problem we recommend that you follow these guidelines.
-  Create a new Branch for your Problem, do not work on master!

-  In your new branch change the directory tree to the following structure:
```
Code
├── <Chapter>
    ├── CMakeLists.txt              (add your <ProblemName> to the list)
    ├── <ProblemName>               
       ├── CMakeLists.txt           (contains only a few lines)
        ├── all 
           ├── problemname.h        (contains the declaration and templated functions)
           ├── problemname.cpp      (contains the definitions)
           ├── problemname_main.cpp (contains the main function that will be executed)
           ├── dependencies.cmake   (contains the path to the source files and libraries you used)         
           ├── test
               ├── problemname_test.cpp (contains a set of tests created with the gtest framework)
               ├── dependencies.cmake   (contains the path to the source files and libraries you used)         
├── scripts 
:  ├── assignment_list.json (Add <Chapter/ProblemName> to the Problems list)       
```

- Run the command `find . -iname '*.h' -o -iname '*.cpp' | xargs clang-format -style=Google -i` inside `<Chapter/ProblemName/all>` 
to reformat the code using the Google style recommendations. To except a line of code from these changes, put the comment `// NOLINT` behind the line or 
`//NOLINTNEXTLINE` above the line that should not be affected by the reformatting. More information on this topic can be found [here](https://clang.llvm.org/extra/clang-tidy/).


- After writing the mastersolution inside `<Chapter/ProblemName/all>` and putting the SOLUTION macros execute:
`python3 deploy.py`, which will genereate the directories `<Chapter/ProblemName/mastersolution>` , `<Chapter/ProblemName/mysolution>` and `<Chapter/ProblemName/templates>`

#### Troubleshooting

- Your `<ProblemName>` does not appear in the build folder after executing CMake?

  - Did you add `<ProblemName>` to the `<Chapter/CMakeLists.txt>`?
  - Did you execute `cmake ..` from build?


- The folders `<Chapter/ProblemName/mastersolution>` etc. are not being created (correctly) by `python3 deploy.py`?

    - Did you add `<ProblemName>` to the `<scripts/assignment_list.json>`?
    - Did you put the SOLUTION macro?
