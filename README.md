# galactic-streaming
A new code to model Cosmic-Ray streaming in the Milky Way

## Prerequesites
To install the code you will currently need to install two
external libraries on your systems, both of which are freely
available. These are at this moment:

* gsl (some necessary computations)
* hdf5 (access and storage of data)

This list might grow in the future.
Apart from that a `C++` compiler and `cmake` are required to
build the code, where a `C++` compiler is needed that allows
at least `C++-14`

## Building the code
To build the code we rely on `cmake`. Correspondingly, you 
first need to creat a sub-directory `build` (or
any other name you like) in the base directory of the project. 
After switching to the new directory, call

    cmake ..

and then compile via

    make -j

This will produce all current executables. They can be found 
in the `bin`
sub-folder of the base directory.

## Using the code

Currently, the code is still in a testing phase. You can find
all tests in the `apps` directory with the corresponding executables
in the `bin` directory. 
