<!--Last edited: 2019-02-22-->
<!--Milos Atz-->

# NwPy: Nuclear Waste analysis in Python 

`nwpy` is a Python package for analyzing characteristics of nuclear wastes and disposal performance for different nuclear energy fuel cycles. `nwpy` is made up of methods that allow for the characterization of fuel cycle waste streams. In addition, `nwpy` currently contains two subpackages:

* `nwpy.attractiveness`: a package for calculating nonproliferation properties of nuclear material streams, including a figure of merit for material attractiveness.
* `nwpy.repository_area`: a package to determine the footprint of a geological repository for the disposal of heat-generating nuclear wastes.

This code repository is developed and maintained by Milos Atz of the UC Berkeley Nuclear Waste Management Group.

#### Characterization of fuel cycle wastes in `nwpy`

In `nwpy`, prescribed mass flow operations are applied to material streams after discharge from a nuclear reactor, producing objects that hold waste canister properties for each HLW stream. The code builds and runs decay cases in ORIGEN-S, a program distributed as a part of the [SCALE](https://www.ornl.gov/scale) software distributed by Oak Ridge National Laboratory through [RSICC](https://rsicc.ornl.gov).

The data that underlies this analysis comes from the DOE [Fuel Cycle Evaluation and Screening report](https://fuelcycleevaluation.inl.gov/SitePages/Home.aspx) (2014). Access to the data can be obtained through [Sandia Connect](https://connect.sandia.gov/). The code is distributed with sample data  from a [2011 Sandia report](http://energy.sandia.gov/wp-content/gallery/uploads/116202.pdf).

The waste streams produced by `nwpy` can be used to evaluate some of the waste management metrics utilized in the Fuel Cycle Evaluation and Screening Study, to characterize waste streams in new ways, or as input to the Python subpackages to analyze other aspects of nuclear waste management.

#### `nwpy.attractiveness`

The package `nwpy.attractiveness` evaluates the attractiveness for proliferation of special nuclear material, such as recycled or waste streams from nuclear fuel cycles. The primary output of the package is a [figure of merit](https://www.tandfonline.com/doi/abs/10.13182/NT10-203) describing this attractiveness. The figure of merit is calculated by first determining three characteristics of the nuclear material:

* the bare sphere critical mass (kg)
* heat content (W/kg)
* dose rate at 1-m from a sphere of 20% of the critical mass (rad/hr)




#### `nwpy.repository_area`

In `nwpy.repository_area`, a repository heat transfer model developed at SNL/LLNL is used to evaluate heat generation and dissipation in close-contact geological repositories for nuclear wastes. For different fuel cycles, in three generic repository concepts, the package can be used to calculate:

1. The repository footprint required to store the waste at a fixed emplacement time.
2. The required interim storage time before waste can be emplaced without violating thermal constraints
2. The amount of interim storage time required to emplace the waste within a given fooprint.

This package accepts data characterizing waste streams from `nwpy` to calculate the repository footprint for each waste stream from a given fuel cycle. The underlying data for the package is from a [2011 Sandia report](http://energy.sandia.gov/wp-content/gallery/uploads/116202.pdf).


## Dependencies

`nwpy` is written in Python 2.7 and relies on a few native packages, detailed below. Additionally, `nwpy` relies on ORIGEN-S, a depletion and decay program distributed with [SCALE](https://www.ornl.gov/scale). The `attractiveness` subpackage requires MCNP. 

To use `nwpy`, first obtain (and install) SCALE by applying for a software license through the [Radiation Safety Information Computational Center](https://rsicc.ornl.gov) (RSICC).

| Package      | Minimum Version    |
|--------------|--------------------|
| `SCALE`      | 6.2                |
| `MCNP`       | 6                  |
| `Python`     | 2.7 (not yet on 3) |
| `numpy`      | 1.11.3             |
| `scipy`      | 0.19.0             |
| `pandas`     | 0.20.1             |
| `matplotlib` | 1.5.1              |
| `pytest`     | 3.2.5              |

Note: the package requires you to add the SCALE `bin/` directory to your shell `PATH`. On Unix systems and/or for bash users, this requires modifying your `.bash_profile`. For example, if you installed SCALE-6.2.2 in the "Applications" subdirectory from the root folder, you might add:

```
export PATH=$PATH:/Applications/SCALE-6.2.2/Contents/Resource/bin/
```

In addition, you must have Python added to your `PATH`. If using `pip` for the installation, it is useful to have that in your `PATH` too.

## Installing nwpy

After SCALE is installed, you can install `nwpy` in any directory by entering the following commands:

```
git clone https://github.com/MilosAtz/nwpy
cd nwpy-master # how it downloads from github
pip install .
```

## Testing nwpy

To run the tests built into `nwpy`, navigate to the directory in which you installed `nwpy` and type:

```
pytest nwpy
```

When you run the tests, you'll notice an "output" directory is produced in the folder from which you've run them. This directory holds the output data from the calculations used in the tests. It is recommended that you delete this directory before running the full set of tests so as to let the code make a new one.

Tests are currently under development. The fuel cycle waste stream characterization methods in `nwpy` have approximately 70\% coverage. Tests for `nwpy.repository_area` are in progress. Tests for `nwpy.attractiveness` are not yet implemented.

## Using nwpy

`nwpy` is a Python package that can be called interactively. From within a Python environment, simply `import` the package.

```
import nwpy
# depending on which package or module you need:
from nwpy import fuelcycle
# or
from nwpy.fuelcycle import stage # less common
# or
from nwpy.repository import repository
# etc
```

Examples about using `nwpy` can be found in the Notebooks subdirectory of this repository. Command-line functionality may be introduced in the future.

## Documentation

`nwpy` is significantly commented and all functions contain descriptive docstrings. Eventually, that documentation will be built in a more formal way.

## Acknowledgements

This research and development is advised by [Professor Massimiliano Fratoni](https://www.nuc.berkeley.edu/people/massimiliano_fratoni) at the University of California Berkeley.

This material is based upon work supported by the National Science Foundation Graduate Research Fellowship under Grant No. 1752814.





<!--## Introduction



. In this report, 40 different fuel cycles were analyzed according to multiple metrics. One of these metrics was "Waste Management" and included characteristics such as mass/volume of waste produced, and radioactivity and toxicity of the waste at future times. Because the goal of the report was to compare the fuel cycles, the authors chose these metrics so as to not skew the results by considering characteristics intrinsic to specific technologies.

However, total radioactivity and toxicity do not represent disposal risk, and the study of disposal performance requires  consideration of additional parameters. My goal is to assess and compare various aspects of the waste management issue for different prospective fuel cycles to build upon the work and findings of the 2014 DOE FCE&S report.

## ORIGEN Decay Calculator
For each fuel cycle included in the DOE FCE&S report, an Analysis Example was used to provide quantitative data that could inform evaluation of performance according to the metrics. A product of the Analysis Examples is isotopic data for each reactor waste stream in the proposed fuel cycle. This data will serve as the source term for all calculations in this study.

However, the isotopic data is a snapshot (usually) at reactor discharge. Where applicable, the appropriate mass balances must be performed to reflect the separations processes included in the fuel cycle. To understand how the isotopic inventory (and associated characteristics, like heat generation) behave over time, this data requires further evaluation in ORIGEN-S. ORIGEN-S can quickly and reliably perform irradiation and decay calculations by the point depletion equation, but the input files are arcane and tedious.

This code imports waste isotopic data from each fuel cycle, along with other fuel cycle data required to obtain the waste composition, and produces ORIGEN-S input files for decay calculations. 

## Installation

More description to come-->