![Speed spectrum of surface‐quasigeostrophic turbulence](https://github.com/houssamyassin/pySQG/blob/master/docs/_static/speed.png?raw=true)

pySQG: a Python Surface Quasigeostrophic Model
================================================

A fork of [pyqg](https://github.com/pyqg/pyqg) that generalizes the model's surface quasigeostrophic module to allow for arbitrary vertical density stratification. 

Surface quasigeostrophy describes the evolution of buoyancy anomalies at the surface of a rapidly rotating, strongly stratified fluid. The theory was originally developed for the atmospheric tropopause ([Blumen 1978](https://doi.org/10.1175/1520-0469(1978)035<0774:UPVFPI>2.0.CO;2), [Held et al. 1995](https://doi.org/10.1017/S0022112095000012)), where the assumption of uniform stratification is appropriate.

Recent ocean observations, show that surface buoyancy anomalies dominate wintertime dynamics in major currents, suggesting that surface quasigeostrophy is a useful model for upper ocean turbulence ([Lapyere \& Klein 2006](https://doi.org/10.1175/JPO2840.1)). However, the surface quasigeostrophic model predicts kinetic energy spectra that are too shallow to be consistent with ocean measurements [(Callies & Ferrari 2013)](https://doi.org/10.1175/JPO-D-13-063.1). 

[Yassin & Griffies (2022)](https://doi.org/10.1175/JPO-D-22-0040.1) demonstrate that this inconsistency arises from the uniform vertical stratification assumption in the original model. The dynamics of surface buoyancy anomalies are highly sensitive to the underlying stratification. The figure above, which shows surface speeds and the underlying stratification profile, shows how vertical stratification controls the dynamical regime of surface buoyancy anomalies. In particular, the ocean is characterized by a weakly stratified mixed layer overlying a stratified interior, which extends the range of mixed-layer eddies and steepens their spectral slopes.

This repository extends the original [pyqg](https://github.com/pyqg/pyqg) surface quasigeostrophy module so it supports arbitrary vertical stratification. Given any vertical buoyancy frequency profile, the model calculates the spectral space inversion relation between the surface buoyancy and the streamfunction by solving a boundary-value problem for the vertical mode structure at each horizontal wavenumber. This fork also adds  spectrally localized stochastic forcing, which is convenient when measuring spectral slopes. The core solvers and documentation remain those of the original model. 

For more details, see:
- Yassin & Griffies (2022), ["Surface Quasigeostrophic Turbulence in Variable Stratification"](https://doi.org/10.1175/JPO-D-22-0040.1), *Journal of Physical Oceanography*, 52, 2995–3013.
- Yassin (2023),  ["The Buoyancy Staircase Limit in Surface Quasigeostrophic Turbulence"](https://doi.org/10.1017/jfm.2023.318), *Journal of Fluid Mechanics*, 962, A35.

