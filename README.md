
> This tool is no longer maintained please use [PMOIRED](https://github.com/amerand/PMOIRED) instead.

# [VLTI/GRAVITY](https://www.eso.org/sci/facilities/paranal/instruments/gravity.html) reduced and calibrated data Quick Look

A Python 3 Tkinter GUI aimed at visualizing VLTI/GRAVITY reduced and calibrated data by the [GRAVITY pipeline](https://www.eso.org/sci/software/pipelines/gravity/gravity-pipe-recipes.html) version 0.9.7 and above.

MacOS Users: This script will crash OS X graphical interface if you use Anaconda. Please use a framework python version (e.g. https://www.python.org/downloads/mac-osx/).

Note that the python 2.7 version is here for legacy ([`gravi_quick_look.py`](gravi_quick_look.py)) by not updated anymore!

## Overview

Run [`gravi_quick_look_3.py`](gravi_quick_look_3.py) as a Python 3 script and select a directory (alternatively, a directory can be given as an argument). The GUI will look like this: ![Figure 1](graviql.png)

* The top rows allow to start show the data plots (in Matplotlib). Note the spectra are roughly corrected for telluric features using a synthetic atmospheric model for 2.0mm of water vapor. Spectral ranges can be selected using buttons or by manually entering the range as the boundaries (in um) separated by a space. "SC/FT" drop down menu allows to select between the science spectrograph (default) of fringe tracker data.

* The "visibility model" box allows to overplot a simple V2 / differential phase / closure phase model. Note that a list of models is kept in file ".graviqlmod" file: if a plot is called without a model (empty string), and you used a model previously, the past model will be called up. The syntax is a coma separated list of "param=value". Possible parameters are:
  - 'ud=1.0' uniform disk diameter, in mas
  - 'f=1.0' continuum flux of the primary (default is 1.0)
  - 'f_2.166_0.002=0.1' gaussian emission (centered 2.166um, fwhm=0.002um)
  - 'fres=0.1' continuum resolved flux
  - 'fres_2.166_0.002=0.1' gaussian spectral emission (centered 2.166um, fwhm=0.002um) of resolved flux
  - 'udc=1.0, xc=1.0, yc=1.2, fc=0.1' uniform disk diameter of a companion, at -1.0mas towards east and 1.2 mas towards north. companion flux is 0.1
  - 'fc_2.166_0.002=0.1': gaussian line (centered 2.166um, fwhm=0.002um) for the companion flux
  - gaussian feature: 'dg' (size) 'xg', 'yg' (position), 'fg' (flux) and 'fg_...' for emission / absorption lines.

* The bottom part of the GUI contains the list of files which can be selected for display, either individually, or averaged. Object' names followed by a ''\*'' indicate data which are reduced and calibrated (as opposed to simply reduced data). Columns description:
  - Object: name of the SC object (same as FT in single feed)
  - Prog.ID and container: ESO specific
  - Mode: spectral and polarisation model
  - DIT: exposure time of SC
  - Baseline: actual quadruplet used for observation
  - See/T0 @V: Seeing (in arcseconds) and Coherence time (in milliseconds), as measured in visible by the observatory's DIMM
  - FT(T0@K): Fraction of reduced frame for the Fringe Tracker, as well as measured coherence time in K band.
  - SC(nB): Fraction of reduced frames for the Science Channel, as well as number of baseline reduced (can be <6 in bad conditions)
  - Correction: the type of correction applied to the SC data. Can be "NONE" or "V-FACTOR".
  - Date-Obs: UT time and date of observation
  - LST: local sidereal time

## Limitations
* Tested on MEDIUM and HIGH dispersion modes only.
* Telluric correction is very approximative and only indicative.
* Tested on [Anaconda3](https://www.continuum.io/downloads)/MacOS and [Anaconda3](https://www.continuum.io/downloads)/Linux.

## Dependencies

* Numpy, Matplotlib and Astropy. GUI in Tkinter (standard library).
