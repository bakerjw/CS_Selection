# Conditional Spectrum ground motion selection
Software for selecting earthquake ground motions to match a target conditional spectrum

This software can be used to select ground motions with a target response spectrum mean and covariance matrix. It improves upon an earlier algorithm by Jayaram et al. (2011).

Further documentation can be found in the following paper:

Lee, C., and Baker, J. W. (2016). “An Improved Algorithm for Selecting Ground Motions to Match a Conditional Spectrum.” (in preparation).



Some potential features to be included in the future are:
* Generalize the get_target_spectrum.m function to compute targets for cases with multiple contributing earthquake sources and/or multiple ground motion prediction equations (see Lin et al. 2013)
* Built a compiled C function to perform the optimization, as this is the computationally expensive portion of the procedure
* Add metadata for additional ground motion databases as they become available
* Modify the compute_scale_factor function to optionally use the approach of Ha and Han (2016), for better computational performance

##References

Ha, S. J., and Han, S. W. (2016). “An efficient method for selecting and scaling ground motions matching target response spectrum mean and variance.” Earthquake Engineering & Structural Dynamics, (in press).

Jayaram, N., Lin, T., and Baker, J. W. (2011). "A computationally efficient ground-motion selection algorithm for matching a target response spectrum mean and variance." Earthquake Spectra, 27(3), 797?815.

Lin, T., Harmsen, S. C., Baker, J. W., and Luco, N. (2013). "Conditional Spectrum Computation Incorporating Multiple Causal Earthquakes and Ground Motion Prediction Models." Bulletin of the Seismological Society of America, 103(2A), 1103?1116.

