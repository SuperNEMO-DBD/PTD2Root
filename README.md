PTD2Root
========
A simple Falaise plugin to convert BRIO format particle tracking data to a flat ROOT TTree

Installing
==========
PTD2Root requires the following software:

- CMake 3.3 or higher
- Falaise 2.0 or higher
  - Falaise will provide the needed Bayeux and ROOT dependencies

These packages can easily be obtained through the cadfaelbrew package manager.
PTD2Root is also supplied through this system, so if you have this installed,
simply do

```console
$ brew install ptd2root
```

To build and install from source, simply do:

```console
$ git clone https://github.com/SuperNEMO-DBD/PTD2Root.git PTD2Root.git
$ cd PTD2Root.git
$ ls
CMakeLists.txt flptd2root.py  p2r.conf       ptd2root.cpp   ptd2root.h     readme.pdf
$ mkdir build
$ cd build
$ cmake ..
$ make && make install
```

Using the Converter
===================
Both a plugin module for Falaise's `flreconstruct` program is supplied, together with
a convenience command-line converter program.

If you already have BRIO output files from `flreconstruct`, these can be converted to
ROOT TTree format using the `flptd2root.py` program. This may be passed the name of the
input BRIO file and desired output ROOT file, e.g.

```console
$ flptd2root.py -i mybriofile.brio -o converted.root
```

Note that `flptd2root.py` requires the `flreconstruct` program to be in your `PATH`.

To use the plugin module in your own pipeline, simply add it as the last step. The
output file can be set with the `filename_out` string key. See the file [`p2r.conf`(p2r.conf)
for an example setup.

Format of the Output TFile/TTree
================================
The format of the output file, including Tree and Branch naming are describing in the
[accompanying document](readme.pdf)


