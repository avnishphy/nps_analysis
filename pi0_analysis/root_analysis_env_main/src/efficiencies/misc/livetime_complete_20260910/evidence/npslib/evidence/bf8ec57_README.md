# NPSlib

This is the reconstruction code for the Hall C Neutral Particle Spectrometer [(NPS)](https://wiki.jlab.org/cuawiki/index.php/Main_Page), implemented as an add-on library for the Hall C analyzer [hcana](https://github.com/JeffersonLab/hcana).

### Building and loading the library

Building the library currently requires a CMake build of hcana (version >= 0.97) and installation of hcana under a to-level directory.

* Set up hcana as usual. It should be sufficient to have `hcana` in your `PATH`. 
* Download or clone this repository and change into its top-level directory.
* Follow the usual CMake workflow:

```
% cmake -B BUILD -S . -DCMAKE_INSTALL_PREFIX=<install location>
% cmake --build BUILD -j4
```
* Troubleshooting: If `hcana` is not found, try setting `HCANA` (or `CMAKE_PREFIX_PATH`) to the hcana installation location, i.e. the top-level directory containing `bin`, `lib(64)`, and `include`. Note: Building against a hcana CMake build directory, or against hcana that was built with `scons` or with the plain Makefiles is not supported.
* Load the library from hcana and check availability of its classes. This example assumes that you are in the top-level directory of NPSlib after building with CMake as indicated above:

```
% hcana -l
hcana [0] gSystem->Load("BUILD/src/libNPS")
(int) 0
hcana [1] .class THcNPSApparatus
==========================================
class THcNPSApparatus
SIZE: 240 FILE: THcNPSApparatus.h LINE: 13
Base classes: ----------------------------
0x0        public THaApparatus
   0x0        public THaAnalysisObject
      0x0        public TNamed
         0x0        public TObject
...
```
* Install the library (optional but recommended):

```
% cmake --install BUILD
```
This will put the library in \<install location\>/lib(64). For convenient use from within `hcana`, add this location to `(DY)LD_LIBRARY_PATH`.

### Class description

The library currently contains the following classes: (Details TODO)

* `THcNPSApparatus`
* `THcNPSArray`
* `THcNPSCalorimeter`
* `THcNPSShowerHit`
* `THcNPSAnalyzer`
* `THcNPSCluster`
* `THcNPSSecondaryKine`
* `THcNPSCoinTime`
* `Decoder::VTPModule`: Driver for VTP front-end electronics.
* `Decoder::VLDModule`: Driver for VLD front-end electronics.

