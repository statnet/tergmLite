# tergmLite

[![Version](http://img.shields.io/badge/Version-2.2.1-orange.svg?style=flat)](https://github.com/statnet/tergmLite/releases)
[![Build Status](https://travis-ci.org/statnet/tergmLite.svg?branch=master)](https://travis-ci.org/statnet/tergmLite)

*Fast Simulation of Simple Temporal Exponential Random Graph Models*

#### Overview

This package provides functions for the computationally efficient simulation and 
resimulation of dynamic networks estimated with the statistical framework of 
temporal exponential random graph models (TERGMs), implemented in the tergm package 
within the Statnet suite of R software. Networks are represented within an edgelist 
format only, with nodal attributes stored separately. Also includes efficient 
functions for the deletion and addition of nodes within that nework representation.

The statistical framework of temporal exponential random graph models (TERGMs)
provides a rigorous, flexible approach to estimating generative models for
dynamic networks and simulating from them for the purposes of modeling infectious
disease transmission dynamics. TERGMs are used within the `EpiModel` software
package to do just that. While estimation of these models is relatively fast,
the resimulation of them using the tools of the `tergm` package
is computationally burdensome, requiring hours to days to iteratively resimulate
networks with co-evolving demographic and epidemiological dynamics. The
primary reason for the computational burden is the use of the `network`
class of object (designed within the package of the same name); these objects
have tremendous flexibility in the types of networks they represent but at the
expense of object size. Continually reading and writing larger-than-necessary
data objects has the effect of slowing the iterative dynamic simulations.

The `tergmLite` package reduces that computational burden by representing
networks less flexibly, but much more efficiently. For epidemic models, the only
types of networks that we typically estimate and simulate from are undirected,
binary edge networks with no missing data (as it is simulated). Furthermore,
the network history (edges or node attributes) does not need to be stored for
research-level applications in which summary epidemiological statistics (e.g.,
disease prevalence, incidence, and variations on those) at the population-level
are the standard output metrics for epidemic models. Therefore, the network
may be stored as a cross-sectional edgelist, which is a two-column matrix
of current edges between one node (in column one) and another node (in column two).
Attributes of the edges that are called within ERGMs may be stored separately in
vector format, as they are in `EpiModel`. With this approach, the simulation
time is sped up by a factor of 25-50 fold, depending on the specific research
application.

#### Version 2.0 Installation Notes

Versions >= 2.0 implement the new `networkLite` API that is implemented across the
`tergmLite`, `network`, and `EpiModel` packages. If you would like to use the 
version before the implementation of this new API, you should install version 
`1.2.0` with:

```r
remotes::install_github("statnet/tergmLite@v1.2.0")
```
