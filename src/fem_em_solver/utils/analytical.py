"""Analytical solutions for electromagnetic validation cases."""

import numpy as np
from typing import Union, Callable


class AnalyticalSolutions:
    """Collection of analytical solutions for EM validation."""
    
    @staticmethod
    def straight_wire_magnetic_field(
        points: np.ndarray,
        current: float,
        wire_position: np.ndarray = None
    ) -> np.ndarray:
        """Analytical B-field for infinite straight wire.
        
        For a wire along the z-axis carrying current I:
            B_φ = μ₀I / (2πr)
        
        In Cartesian coordinates:
            B_x = -μ₀I * y / (2πr²)
            B_y =  μ₀I * x / (2πr²)
            B_z = 0
        
        Parameters
        ----------
        points : np.ndarray
            Array of shape (n, 3) with evaluation points
        current : float
            Current in wire [A]
        wire_position : np.ndarray, optional
            (x, y) position of wire. Default is (0, 0)
            
        Returns
        -------
        np.ndarray
            B-field at points, shape (n, 3)
        """
        if wire_position is None:
            wire_position = np.array([0.0, 0.0])
        
        mu_0 = 4 * np.pi * 1e-7  # H/m
        
        # Extract coordinates relative to wire
        x = points[:, 0] - wire_position[0]
        y = points[:, 1] - wire_position[1]
        
        # Cylindrical radius
        r = np.sqrt(x**2 + y**2)
        
        # Avoid division by zero at wire location
        r = np.maximum(r, 1e-10)
        
        # B-field magnitude (azimuthal direction)
        B_phi = mu_0 * current / (2 * np.pi * r)
        
        # Convert to Cartesian
        # B_φ direction is (-y/r, x/r, 0)
        B = np.zeros_like(points)
        B[:, 0] = -B_phi * y / r  # B_x
        B[:, 1] =  B_phi * x / r  # B_y
        B[:, 2] = 0.0              # B_z
        
        return B
    
    @staticmethod
    def straight_wire_vector_potential(
        points: np.ndarray,
        current: float,
        wire_position: np.ndarray = None,
        A_ref: float = 0.0
    ) -> np.ndarray:
        """Analytical A-field for infinite straight wire.
        
        The magnetic vector potential for a wire along z-axis:
            A_z = -μ₀I / (2π) * ln(r/r_ref)
        
        Choosing r_ref such that A=0 at some reference point.
        
        Parameters
        ----------
        points : np.ndarray
            Array of shape (n, 3) with evaluation points
        current : float
            Current in wire [A]
        wire_position : np.ndarray, optional
            (x, y) position of wire
        A_ref : float
            Reference potential (gauge choice)
            
        Returns
        -------
        np.ndarray
            A-field at points, shape (n, 3)
        """
        if wire_position is None:
            wire_position = np.array([0.0, 0.0])
        
        mu_0 = 4 * np.pi * 1e-7
        
        x = points[:, 0] - wire_position[0]
        y = points[:, 1] - wire_position[1]
        r = np.sqrt(x**2 + y**2)
        
        # Avoid log(0)
        r = np.maximum(r, 1e-10)
        
        # A has only z-component
        A = np.zeros_like(points)
        A[:, 2] = -mu_0 * current / (2 * np.pi) * np.log(r) + A_ref
        
        return A
    
    @staticmethod
    def circular_loop_magnetic_field_on_axis(
        z: np.ndarray,
        current: float,
        radius: float,
        loop_center: float = 0.0
    ) -> np.ndarray:
        """B-field on axis of circular current loop.
        
        For a loop of radius a in the xy-plane, centered at z=z0:
            B_z(z) = μ₀Ia² / (2(a² + (z-z0)²)^(3/2))
        
        Parameters
        ----------
        z : np.ndarray
            Positions along z-axis [m]
        current : float
            Current in loop [A]
        radius : float
            Loop radius [m]
        loop_center : float
            z-position of loop center [m]
            
        Returns
        -------
        np.ndarray
            B_z at positions z
        """
        mu_0 = 4 * np.pi * 1e-7
        
        dz = z - loop_center
        denom = 2 * (radius**2 + dz**2)**(3/2)
        
        B_z = mu_0 * current * radius**2 / denom
        
        return B_z
    
    @staticmethod
    def helmholtz_coil_field_on_axis(
        z: np.ndarray,
        current: float,
        radius: float,
        separation: float = None
    ) -> np.ndarray:
        """B-field on axis of Helmholtz coil.
        
        Helmholtz coil: two identical loops separated by distance = radius,
        carrying current in same direction. This configuration gives
        maximally uniform field in the center.
        
        Parameters
        ----------
        z : np.ndarray
            Positions along z-axis [m]
        current : float
            Current in each loop [A]
        radius : float
            Loop radius [m]
        separation : float, optional
            Distance between loops. Default is radius (Helmholtz condition)
            
        Returns
        -------
        np.ndarray
            B_z at positions z
        """
        if separation is None:
            separation = radius  # Helmholtz condition
        
        # Center the configuration at z=0
        z1 = -separation / 2
        z2 = separation / 2
        
        B1 = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(
            z, current, radius, z1
        )
        B2 = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(
            z, current, radius, z2
        )
        
        return B1 + B2


class ErrorMetrics:
    """Error metrics for comparing numerical and analytical solutions."""
    
    @staticmethod
    def l2_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute L2 norm of error.
        
        ||u_num - u_ana||₂ = sqrt(∫ |u_num - u_ana|² dx)
        
        For discrete points:
            L2_error = sqrt(sum(|u_num - u_ana|²))
        """
        diff = numerical - analytical
        return np.sqrt(np.sum(diff**2))
    
    @staticmethod
    def l2_relative_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute relative L2 error."""
        l2_err = ErrorMetrics.l2_error(numerical, analytical)
        l2_ana = np.sqrt(np.sum(analytical**2))
        
        if l2_ana < 1e-15:
            return l2_err
        
        return l2_err / l2_ana
    
    @staticmethod
    def max_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute maximum absolute error (L∞ norm)."""
        return np.max(np.abs(numerical - analytical))
    
    @staticmethod
    def max_relative_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute maximum relative error."""
        abs_err = np.abs(numerical - analytical)
        max_abs_ana = np.max(np.abs(analytical))
        
        if max_abs_ana < 1e-15:
            return np.max(abs_err)
        
        return np.max(abs_err / (np.abs(analytical) + 1e-15))
