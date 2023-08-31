iFraudSimulator: Generating synthetic insurance fraud network data.
====
  Package to generate synthetic insurance fraud network data.

<p align="left">
  <img src="inst/figures/iFraudSimulator.png" width="25%">
</p>

## Installation

### On current R (>= 3.0.0)
* Development version from Github:

```
library("devtools"); install_github("BavoDC/iFraudSimulator", dependencies = TRUE, build_vignettes = TRUE)
```

(This requires `devtools` >= 1.6.1, and installs the "master" (development) branch.)
This approach builds the package from source, i.e. `make` and compilers must be installed on your system -- see the R FAQ for your operating system; you may also need to install dependencies manually.

## Documentation
The basic functionality of the package is explained and demonstrated in the vignette, which you can access using
```
vignette("iFraudSimulator")
```

The working paper can be found on [arXiv](https://arxiv.org/abs/2308.11659).


## Contact
If you have questions, remarks or suggestions regarding the package, you can contact me at [bavo.campo@kuleuven.be](mailto:bavo.campo@kuleuven.be).

## Citation
If you use this package, please cite:

- Campo, Bavo D.C., and Antonio, Katrien (2023). An engine to simulate insurance fraud network data. arXiv:2308.11659, available at https://arxiv.org/abs/2308.11659.
