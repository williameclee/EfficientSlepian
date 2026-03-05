# EfficientSlepian

An open-source MATLAB re-implementation of the efficient way to compute Slepian functions, as described in [Bates, Alice P et al. (2017)](https://doi.org/10.1109/TSP.2017.2712122).

## Usage

The main entry function is `glmalpha_eff.m`, as in the following example:

```matlab
[G, V, N] = glmalpha_eff(domain, L)
```

## Dependencies

* [slepian_alpha](https://github.com/csdms-contrib/slepian_alpha.git) is the main codebase that is doing all the heavy lifting, including the integration and rotation of the Slepian functions.
* [ULMO](https://github.com/williameclee/ulmo.git) contains some utilities functions, but should be easily replaceable.

---

Author: [En-Chi Lee](https://github.com/williameclee) ([williameclee@arizona.edu](mailto:williameclee@arizona.edu))
