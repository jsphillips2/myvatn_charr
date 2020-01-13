Mývatn Arctic Charr
========

Analyses associated with the manuscript "Estimating time-varying demographic rates from age-structured abundance: application to a historically variable fishery."
-------

Joe Phillips, Hólar Universtiy and University of Wisconsin-Madison

## Description

This repository contains code for inferring temporal variation in demographic rates from age-strucutred observations of abundance, applied to the Arctic Charr (*Salvelinus alpinus*) population of Lake Mývatn, Iceland. The method models time-varying survival and recruitment as random walks, allowing them to varying smoothly through time while taking advantage of the full set of data for estimating all parameters simultaneously. The model is fit in a Bayesian framework using Stan, run from R using the 'rstan' package. This repository contains all of the raw data and code for reproducing the analyses. 

## Contents

* `data`: Raw data files and code for processing data.

* `analyses`: Analyses of the model fits, including code for model comparison and demographic analysis.

* `model`: Files for specifying and fitting the model in Stan. Different input and output folders (`model_full`,`model_full_bias`,`model_rho`, `model_phi`, and `model_fixed`) correspond to various model fits, differing either due to the input data add/or the Stan file used to specify the model. All of the model fits can be run from the file `fit_model.R`. The folder `stan` contains files specifying the structure of different versions of the model. The file `model_full.stan` specifies the full model allowing both survival and recruitment to vary through time. The files `model_rho.stan`, `model_phi.stan`, and `model_fixed.stan` fit reduced versions of the model with either survival, recruitment, and both fixed through time (respectively). 

* `supplement`: Rmd file for generating supplementary materials associated with manuscript.   