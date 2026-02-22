"""Post-processing and analysis."""

from .evaluation import evaluate_vector_field_parallel
from .phantom_fields import (
    compute_phantom_eb_metrics_and_export,
    compute_tagged_vector_magnitude_stats,
    export_tagged_field_samples_csv,
)

__all__ = [
    "evaluate_vector_field_parallel",
    "compute_tagged_vector_magnitude_stats",
    "export_tagged_field_samples_csv",
    "compute_phantom_eb_metrics_and_export",
]
