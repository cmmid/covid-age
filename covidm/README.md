# covidm
Dynamic model of SARS-nCoV-2 transmission


## Quick start guide

### Installing dependencies for Mac OS

You will need to install gfortran binaries from here: https://github.com/fxcoudert/gfortran-for-macOS/releases

Once installed, run `gcc --version` in terminal to get your current version, e.g. `Target: x86_64-apple-darwin18.8.2.0`. Then run below in terminal to add library path for R:

`cd ~
mkdir .R
cd .R
echo FLIBS=-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin18/8.2.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm >> Makevars
`

Finally, install nlopt: `brew install nlopt`

### Running example

The basic model framework is in the `examples` folder. This includes the following R files:

> `examples/0-libraries.R` - Load required R libraries

> `examples/1-getting-started.R` - Compile code and run simulation

> `examples/2-interventions.R` - Run intervention scenarios

> `examples/3-processes.R` - Set up an observation model

Model parameters are defined in `parameters_ref.txt`.
