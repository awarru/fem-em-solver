"""Physical constants and unit conversions for EM simulations."""

import numpy as np

# Physical constants (SI units)
MU_0 = 4 * np.pi * 1e-7          # Permeability of free space [H/m]
EPSILON_0 = 8.854187817e-12      # Permittivity of free space [F/m]
C_0 = 299792458                  # Speed of light [m/s]
ETA_0 = 376.730313461            # Impedance of free space [Ohm]

# Common material properties
class MaterialProperties:
    """Material properties at common MRI frequencies."""
    
    # Vacuum / Air
    VACUUM = {
        'mu_r': 1.0,
        'epsilon_r': 1.0,
        'sigma': 0.0,
    }
    
    # Copper (conductors)
    COPPER = {
        'mu_r': 0.999991,
        'sigma': 5.8e7,  # S/m
    }
    
    # Gelled saline phantom (typical at 128 MHz / 3T)
    # Conductivity tunable with salt concentration
    GELLED_SALINE_3T = {
        'mu_r': 1.0,
        'epsilon_r': 78.0,
        'sigma': 0.7,  # S/m (tunable 0.6-0.9)
        'density': 1000,  # kg/m^3
    }
    
    # Brain tissue (gray matter, 128 MHz)
    BRAIN_GM_3T = {
        'mu_r': 1.0,
        'epsilon_r': 52.0,
        'sigma': 0.6,  # S/m
        'density': 1040,  # kg/m^3
    }
    
    # Brain tissue (white matter, 128 MHz)
    BRAIN_WM_3T = {
        'mu_r': 1.0,
        'epsilon_r': 42.0,
        'sigma': 0.4,  # S/m
        'density': 1040,  # kg/m^3
    }
    
    # CSF (cerebrospinal fluid, 128 MHz)
    CSF_3T = {
        'mu_r': 1.0,
        'epsilon_r': 74.0,
        'sigma': 2.0,  # S/m
        'density': 1007,  # kg/m^3
    }


# Frequency-dependent properties
class FrequencyDependent:
    """Frequency-dependent material properties."""
    
    @staticmethod
    def gelled_saline_conductivity(frequency_hz, base_sigma=0.7):
        """
        Estimate conductivity of gelled saline at different frequencies.
        
        Parameters
        ----------
        frequency_hz : float
            Frequency in Hz
        base_sigma : float
            Conductivity at 128 MHz (S/m)
            
        Returns
        -------
        float
            Conductivity at specified frequency (S/m)
        """
        # Conductivity increases slightly with frequency due to ionic mobility
        # This is a simplified model
        f_ref = 128e6  # 128 MHz reference
        return base_sigma * (1 + 0.1 * np.log10(frequency_hz / f_ref))
    
    @staticmethod
    def skin_depth(material_sigma, frequency_hz, mu_r=1.0):
        """
        Calculate skin depth.
        
        Parameters
        ----------
        material_sigma : float
            Conductivity (S/m)
        frequency_hz : float
            Frequency (Hz)
        mu_r : float
            Relative permeability
            
        Returns
        -------
        float
            Skin depth (m)
        """
        omega = 2 * np.pi * frequency_hz
        mu = mu_r * MU_0
        return np.sqrt(2 / (omega * mu * material_sigma))


# MRI-specific constants
class MRIConstants:
    """MRI system constants."""
    
    # Field strengths and corresponding Larmor frequencies (proton)
    B0_1_5T = 1.5  # Tesla
    F0_1_5T = 63.87e6  # Hz (63.87 MHz)
    
    B0_3T = 3.0  # Tesla
    F0_3T = 127.74e6  # Hz (127.74 MHz)
    
    B0_7T = 7.0  # Tesla
    F0_7T = 297.2e6  # Hz (297.2 MHz)
    
    B0_9_4T = 9.4  # Tesla
    F0_9_4T = 400.0e6  # Hz (400 MHz)
    
    # SAR limits (IEC 60601-2-33)
    SAR_WHOLE_BODY_LIMIT = 4.0  # W/kg (normal mode)
    SAR_PARTIAL_BODY_LIMIT = 8.0  # W/kg (normal mode)
    SAR_HEAD_LIMIT = 3.2  # W/kg (10g average, normal mode)
    
    # Phantom specifications
    HEAD_PHANTOM_DIAMETER = 0.16  # m (16 cm)
    BODY_PHANTOM_DIAMETER = 0.34  # m (34 cm)
    PHANTOM_HEIGHT = 0.30  # m (30 cm)
