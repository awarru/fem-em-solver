"""Phantom material models for MRI-oriented simulations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..utils.constants import EPSILON_0


@dataclass(frozen=True)
class GelledSalinePhantomMaterial:
    """Minimal gelled-saline phantom material container for frequency-domain solves.

    Parameters
    ----------
    sigma
        Conductivity in S/m.
    epsilon_r
        Relative permittivity (dimensionless).
    frequency_hz
        Operating frequency in Hz.
    mu_r
        Relative permeability. Defaults to 1.0.
    """

    sigma: float
    epsilon_r: float
    frequency_hz: float
    mu_r: float = 1.0

    def validate(self) -> None:
        """Raise ValueError when properties are outside the MVP-valid range."""
        if self.sigma < 0:
            raise ValueError("phantom sigma must be non-negative")
        if self.epsilon_r <= 0:
            raise ValueError("phantom epsilon_r must be positive")
        if self.frequency_hz <= 0:
            raise ValueError("phantom frequency_hz must be positive")
        if self.mu_r <= 0:
            raise ValueError("phantom mu_r must be positive")

    @property
    def omega(self) -> float:
        return float(2.0 * np.pi * self.frequency_hz)

    @property
    def displacement_term(self) -> float:
        """Return omega*epsilon for the phantom material [S/m equivalent]."""
        return float(self.omega * EPSILON_0 * self.epsilon_r)

    @property
    def conduction_plus_displacement(self) -> float:
        """Return sigma + omega*epsilon (used in MVP diagnostics/pipeline wiring)."""
        return float(self.sigma + self.displacement_term)
