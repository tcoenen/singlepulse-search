singlepulse-search
==================

Some library functions and scripts that can be used with the output of the
single-pulse search implementation of PRESTO [1]. All plots for ssps need the
brp plotting library [2].

Contents
========

* ssps_condense.py Plot a 2d-histogram of single-pulse candidates on the
  time-DM plane (useful for quickly assesing data quality and checking for the
  presence of bright pulsars).

* ssps_grab.py Extract pulse trains from the data and plot the results.

* Library functions (under ssps/) that can read / write / manipulate PRESTO
  single-pulse search output.


[1] https://github.com/scottransom/presto
[2] https://github.com/tcoenen/brp
