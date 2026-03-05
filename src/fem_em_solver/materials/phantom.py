"""Phantom material models for MRI-oriented simulations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from ..utils.constants import EPSILON_0


_PHANTOM_PRESET_TABLE: Mapping[str, tuple[float, float]] = {
    # (sigma [S/m], epsilon_r [-])
    "low": (0.45, 80.0),
    "mid": (0.75, 78.0),
    "high": (1.05, 76.0),
}


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

    @classmethod
    def from_preset(
        cls,
        preset: str,
        *,
        frequency_hz: float,
        mu_r: float = 1.0,
    ) -> "GelledSalinePhantomMaterial":
        """Build a gelled-saline material from a named low/mid/high preset.

        Notes
        -----
        Presets are intentionally lightweight MVP defaults for MRI-like workflows.
        """

        key = preset.strip().lower()
        if key not in _PHANTOM_PRESET_TABLE:
            raise ValueError(
                f"Unknown gelled saline preset {preset!r}. "
                f"Expected one of: {sorted(_PHANTOM_PRESET_TABLE)}"
            )

        sigma, epsilon_r = _PHANTOM_PRESET_TABLE[key]
        material = cls(sigma=sigma, epsilon_r=epsilon_r, frequency_hz=frequency_hz, mu_r=mu_r)
        material.validate()
        return material

    @staticmethod
    def preset_names() -> tuple[str, ...]:
        """Return available preset names in deterministic order."""
        return tuple(sorted(_PHANTOM_PRESET_TABLE))

    def with_frequency_adjustment(
        self,
        frequency_hz: float,
        *,
        conductivity_exponent: float = 0.0,
        permittivity_exponent: float = 0.0,
    ) -> "GelledSalinePhantomMaterial":
        """Return a frequency-adjusted material via lightweight power-law hooks.

        This is an MVP hook for frequency dependence. Setting exponents to 0 keeps
        constant properties. Non-zero exponents scale each property by:

            (f_new / f_ref) ** exponent

        where ``f_ref`` is this material's ``frequency_hz``.
        """

        if not np.isfinite(frequency_hz) or frequency_hz <= 0:
            raise ValueError(f"frequency_hz must be finite and positive, got {frequency_hz!r}")
        if not np.isfinite(conductivity_exponent):
            raise ValueError("conductivity_exponent must be finite")
        if not np.isfinite(permittivity_exponent):
            raise ValueError("permittivity_exponent must be finite")

        ratio = float(frequency_hz / self.frequency_hz)
        sigma = float(self.sigma * ratio ** conductivity_exponent)
        epsilon_r = float(self.epsilon_r * ratio ** permittivity_exponent)
        adjusted = GelledSalinePhantomMaterial(
            sigma=sigma,
            epsilon_r=epsilon_r,
            frequency_hz=frequency_hz,
            mu_r=self.mu_r,
        )
        adjusted.validate()
        return adjusted

    def validate(self) -> None:
        """Raise ValueError when properties are outside the MVP-valid range."""
        if not np.isfinite(self.sigma) or self.sigma < 0:
            raise ValueError("phantom sigma must be finite and non-negative")
        if not np.isfinite(self.epsilon_r) or self.epsilon_r <= 0:
            raise ValueError("phantom epsilon_r must be finite and positive")
        if not np.isfinite(self.frequency_hz) or self.frequency_hz <= 0:
            raise ValueError("phantom frequency_hz must be finite and positive")
        if not np.isfinite(self.mu_r) or self.mu_r <= 0:
            raise ValueError("phantom mu_r must be finite and positive")

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

    @property
    def complex_admittivity(self) -> complex:
        """Return sigma + j*omega*epsilon as a compact frequency-domain helper."""
        return complex(self.sigma, self.displacement_term)
