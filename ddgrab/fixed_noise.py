"""Expose a frozen C = diag(sigma**2) + basis @ diag(phi) @ basis.T to PINT.

Arrays use PINT TOA order, not enterprise's sorted order. This component is
runtime-only: save/reload the companion NPZ together with the result par file.
"""

import hashlib

import astropy.units as u
import numpy as np
from pint.models.noise_model import CorrelatedNoiseComponent, NoiseComponent


def toa_fingerprint(toas):
    """Identify TOAs, row order, frequencies, uncertainties and flags."""
    digest = hashlib.sha256()
    for values in (toas.get_mjds().value, toas.get_freqs().to_value(u.MHz),
                   toas.get_errors().to_value(u.s)):
        digest.update(np.asarray(values, dtype="<f8").tobytes())
    digest.update(repr(list(map(str, toas.get_obss()))).encode())
    digest.update(repr([sorted(dict(flags).items()) for flags in toas.table["flags"]]).encode())
    return digest.hexdigest()


class FixedTOANoise(CorrelatedNoiseComponent):
    register = False
    category = "fixed_toa_noise"
    is_time_correlated = True

    def __init__(self, toas, sigma, basis, phi):
        super().__init__()
        self.sigma = np.asarray(sigma, dtype=float).copy()
        self.basis = np.asarray(basis, dtype=float).copy()
        self.phi = np.asarray(phi, dtype=float).copy()
        n = len(toas)
        if self.sigma.shape != (n,) or self.basis.ndim != 2 or self.basis.shape[0] != n:
            raise ValueError("Fixed covariance dimensions do not match the TOAs")
        if self.phi.shape != (self.basis.shape[1],) or self.phi.size == 0:
            raise ValueError("Expected at least one correlated-noise basis column")
        if not all(np.all(np.isfinite(x)) for x in (self.sigma, self.basis, self.phi)):
            raise ValueError("Fixed covariance contains non-finite values")
        if np.any(self.sigma <= 0) or np.any(self.phi <= 0):
            raise ValueError("Fixed covariance variances must be positive")
        self.fingerprint = toa_fingerprint(toas)
        self.scaled_toa_sigma_funcs = [self.scaled_sigma]
        self.basis_funcs = [self.basis_weight]

    def check_toas(self, toas):
        if toa_fingerprint(toas) != self.fingerprint:
            raise ValueError("Fixed covariance requires the same TOAs in the same order")

    def scaled_sigma(self, toas):
        self.check_toas(toas)
        return self.sigma * u.s

    def basis_weight(self, toas):
        self.check_toas(toas)
        return self.basis, self.phi

    def get_noise_basis(self, toas):
        return self.basis_weight(toas)[0]

    def get_noise_weights(self, toas):
        return self.basis_weight(toas)[1]

    def save(self, path):
        np.savez_compressed(path, sigma=self.sigma, basis=self.basis, phi=self.phi,
                            toa_fingerprint=self.fingerprint, row_order="pint")


def attach_fixed_noise(model, toas, sigma, basis, phi):
    """Replace existing noise components to avoid double-counting covariance."""
    for name, component in list(model.components.items()):
        if isinstance(component, NoiseComponent):
            model.remove_component(name)
    component = FixedTOANoise(toas, sigma, basis, phi)
    model.add_component(component, validate=False)
    model.setup()
    model.validate()
    return component
