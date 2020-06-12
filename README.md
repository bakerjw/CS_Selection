# Conditional Spectrum ground motion selection
Software for selecting earthquake ground motions to match a target conditional spectrum

This software can be used to select ground motions with a target response spectrum mean and covariance matrix. It improves upon an earlier algorithm by Jayaram et al. (2011).

Further documentation can be found in the following paper:

Baker, J. W., and Lee, C. (2018). “An Improved Algorithm for Selecting Ground Motions to Match a Conditional Spectrum.” Journal of Earthquake Engineering, 22(4), 708–723.

Vertical component selection was added in 2020, and is documented in

Kwong, N. S., Jaiswal, K. S., Luco, N., and Baker, J. W. (2020). “Selecting three components of ground motions from conditional spectra for multiple stripe analyses.” 17th World Conference on Earthquake Engineering, Sendai, Japan.


## Potential features to be included in the future:
* Generalize the get_target_spectrum.m function to compute targets for cases with multiple contributing earthquake sources and/or multiple ground motion prediction equations (see Lin et al. 2013)
* Add metadata for additional ground motion databases as they become available
* Modify the compute_scale_factor function to optionally use the approach of Ha and Han (2016), for better computational performance

## Known issues:
* As of July 2016, PEER has stopped providing the time series for the NGA-West1 database, so the URLs provided from this software are no longer usable.
* The NGA-West2 database includes spectra for some ground motions that are not provided by PEER's online tool (e.g., ground motions from the Wenchuan earthquake).
* As of 2020, the BroadBand Platform time series are no longer available online, so the URLs listed in the metadata do not work.

## References

Ha, S. J., and Han, S. W. (2016). “An efficient method for selecting and scaling ground motions matching target response spectrum mean and variance.” Earthquake Engineering & Structural Dynamics, 45(8), 1381–1387.

Jayaram, N., Lin, T., and Baker, J. W. (2011). "A computationally efficient ground-motion selection algorithm for matching a target response spectrum mean and variance." Earthquake Spectra, 27(3), 797-815.

Lin, T., Harmsen, S. C., Baker, J. W., and Luco, N. (2013). "Conditional Spectrum Computation Incorporating Multiple Causal Earthquakes and Ground Motion Prediction Models." Bulletin of the Seismological Society of America, 103(2A), 1103-1116.

