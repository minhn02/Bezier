Bezier is a C++ library for constructing Bezier curves, specifically for trajectories. The library is designed to be used as a single header file which contains classes for single curve and spline construction. The implementation of the Spline is mostly tailored towards robotic trajectories from *Trajectory Planning for Automatic Machines and Robots*.

## Functionality

### Curve
The curve class defines a [Bezier Curve](https://en.wikipedia.org/wiki/B%C3%A9zier_curve) which can can hold ```n>0``` points in ```d``` dimensional space, defined over ```t0``` to ```t0+T``` (which defaults to being defined over ```[0, 1]```).

The curve can be evaluated in time and its n'th derivative.

### Spline
The spline class defines a sequence of cubic bezier curves that have continuity in position. The spline can be instantiated with a sequence of curves (although continuity is not checked), a sequence of points, or a function of time. The Spline constructor will interplate the latter two to create the sequence of cubic bezier curves.

The curve can be evaluated in time and its n'th derivative.

Cubic bezier curves make up the splines as of now, but it could be generalized to any order curves, although the required level of continuity would need to be specified somehow.

## Installation
The library depends on [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page).

It can be installed by putting the ```bezier.h``` file into a local project and adding ```#include "bezier.h"``` to the file. There is also a CMake install target to ```/usr/local/lib```.

## Testing
Basic tests are defined in the Cmake file and cover simple curve construction, derivative evaluation, and spline interpolation. More comprehensive coverage is yet to be implented.

The library has only been tested on Ubuntu.

## TODO
1. Add argument checking for safety 
2. add more spline interpolation methods, such as implementing different ways of time spacing along the spline.

## Author
[Minh Nguyen](https://github.com/minhn02) - minh@minh02.com