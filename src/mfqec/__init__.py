"""Measurement-free versus feed-forward error correction for a trapped-ion memory."""
from .schemes import (Params, bare_fail, measurement_based_fail, measurement_free_fail,
                      qubit_count, total_time_mb, total_time_mf)

__all__ = ["Params", "bare_fail", "measurement_based_fail", "measurement_free_fail",
           "qubit_count", "total_time_mb", "total_time_mf"]
__version__ = "0.1.0"
