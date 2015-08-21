# IEEE-EMC2015
This repository contains the [presentation slides](./Presentation_IEEE_EMC15.pdf "Slides"), the source codes used to produce the figures in the article:
### *Source Stirring Analysis in a Reverberation Chamber Based on Modal Expansion of the Electric Field*
- Emmanuel Amador, EDF Lab, LME, emmanuel.amador@edf.fr
- Philippe Besnier IETR-UMR CNRS 6164, INSA Rennes, philippe.besnier@insa-rennes.fr

presented at the *Joint IEEE International Symposium on Electromagnetic Compatibility and EMC Europe, Dresden 2015.*

This repository allows to reproduce the results presented in the article.
Basically, the numerical model based on a modal expansion of the electric field in a rectangular cavity is given. This model allows to compute the electric field in the volume of the cavity for a given transmitter position.

The programs included produce the different figures of the article and some extra programs gives more latitude to anyone who would like to experiment with the model developed in the article.
Comments and suggestions are welcomed.

###Packages needed:
- Julia 0.3+ and the follwing packages: PyPlot, NPZ, Grid.
- Python 2.7+ and Matplotlib, Numpy.



Source stirring animation: `./code/animation_sourcestirring.jl`

[![ScreenShot](http://img.youtube.com/vi/84I4xAoDW-8/0.jpg)](http://youtu.be/84I4xAoDW-8)


