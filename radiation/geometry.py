"""Detector solid-angle geometry."""
import numpy as np


def solid_angle_disk(r_cm, a_cm):
    """On-axis solid angle (sr) of a circular disk of radius a_cm at distance r_cm."""
    r_cm = np.asarray(r_cm, dtype=float)
    return 2 * np.pi * (1 - r_cm / np.sqrt(r_cm**2 + a_cm**2))
