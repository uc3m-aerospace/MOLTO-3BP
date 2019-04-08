# Spice folder
Download the Matlab interface of the [Spice Toolkit](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html) and include it in this directory.

The kernels required to run the default molto-3bp are already included. In case you need a different ephemerides file, include them in the /kernels folder and modify the function *load_spice_kernels.m* to load the new kernels you are going to use.
