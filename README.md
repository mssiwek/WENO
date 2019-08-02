# WENO in python3

This repository contains the python library weno.py which contains WENO functions
that take a 5 point stencil of cell-centered hydrodynamical quantities as input and
return a 5-th order non-linearly weighted polynomial which numerically approximates
solutions of hyperbolic conservation laws to 5-th order convergence in smooth regions.
This fit becomes 1st order accurate around discontinuities, but avoids
oscillatory behaviour by essentially discarding stencils that contain discontinuities.
The repository further contains a 1-D advection code that demonstrates the WENO method: advect_weno.py

### Prerequisites

Python3, numpy, matplotlib.pyplot

## Installation

Clone the repository and ensure all required libraries are installed.
Then simply run the 1-D advection code by opening a terminal in the WENO directory and typing:
```
>> python advect_weno.py
```

The code asks for several input arguments:

```
>> Which IC function do you want to advect? Examples include: gaussian, simple_step, square_well, quadratic_well, trig_disc. See test_functions.py for more info.
```

Respond by choosing one of the functions you want to advect. Second, you will be asked for the resolution:

```
>> At which resolution?
```

Type an integer number of grid vertices. Then, decide if you want to calculate, store and plot the L2 error for a convergence test.

```
>> Do you want to plot the L2 error? This may take some time. [y/n]
```

Type y or n. Calculating the L2 error will take some time for high resolutions.

```
>> End time:
```

Choose for how long the wave pulse should be advected. The default wavespeed is 1, so t_final = 2 is equivalent to one cycle.
Finally, in addition to the 5th order convergent WENO method you can show a 1st order convergent in the same figure for comparison:

```
>> Do you also want to plot the first order solution for comparison? [y/n]
```

Type y or n.

## Authors

* Magdalena Siwek

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Jonathan Zrake, Andrew MacFadyen
