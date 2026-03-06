# EfficientSlepian

An open-source MATLAB re-implementation of the efficient way to compute Slepian functions, as described in [Bates, Alice P et al. (2017)](https://doi.org/10.1109/TSP.2017.2712122).

## Theoretical background

In this efficient formulation, there are essentially two stages to compute the Slepian basis:

1. Compute the Slepian functions for the polar cap that encloses the domain, which forms a much smaller basis.
2. Compute, again, the Slepian functions for the polar cap Slepian basis over the rotated domain.

The 'Slepian functions of the Slepian functions' are then projected back to the spherical harmonics and rotated back to the original domain to get the final projection matrix G.
Because the first stage can exploit the axisymmetry of the polar cap, and the second stage is only computed over a much smaller number of basis functions, this efficient formulation can be much faster than the direct formulation (`glmalpha.m` in [slepian_alpha](https://github.com/csdms-contrib/slepian_alpha.git)) for large bandwidths and small domains.

## Usage

The main entry function is `glmalpha_eff.m`, as in the following example:

```matlab
[G, V, N] = glmalpha_eff(domain, L)
```

## Dependencies

* [slepian_alpha](https://github.com/csdms-contrib/slepian_alpha.git) is the main codebase that is doing all the heavy lifting, including the integration and rotation of the Slepian functions.
* [ULMO](https://github.com/williameclee/ulmo.git) contains some utilities functions, but should be easily replaceable.

---

Whilst this repository does not contain the original, unpublished code by Alice Bates, some design choices are still inspired by the original implementation.

Author: 2026/03/05, [En-Chi Lee](https://github.com/williameclee) ([williameclee@arizona.edu](mailto:williameclee@arizona.edu))
