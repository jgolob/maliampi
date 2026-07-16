"""SV identity registry creation and validation for reference packages."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import cast

import pandas as pd

from .sv import SvValidationError, read_sv_h5ad

REGISTRY_COLUMNS = ("sv_id", "sequence_sha256", "sequence")


def write_sv_registry(sv_h5ad: str | Path, output_parquet: str | Path) -> None:
    """Write the refpkg SV identity registry from a validated source artifact."""
    artifact = read_sv_h5ad(sv_h5ad)
    registry = cast(pd.DataFrame, artifact.var).loc[:, REGISTRY_COLUMNS].copy()
    registry.to_parquet(output_parquet, index=False)


def validate_sv_registry(sv_h5ad: str | Path, registry_parquet: str | Path) -> dict[str, str]:
    """Reject any known sequence whose registry identity differs from the H5AD."""
    artifact = read_sv_h5ad(sv_h5ad)
    registry = pd.read_parquet(registry_parquet)
    if set(REGISTRY_COLUMNS) - set(registry.columns):
        raise SvValidationError("Refpkg SV registry lacks required identity columns")
    if registry["sequence_sha256"].duplicated().any() or registry["sv_id"].duplicated().any():
        raise SvValidationError("Refpkg SV registry contains duplicate identities")
    known = registry.set_index("sequence_sha256")
    source = cast(pd.DataFrame, artifact.var).loc[:, REGISTRY_COLUMNS]
    for sv_id, sequence_sha256, sequence in source.itertuples(index=False, name=None):
        if sequence_sha256 in known.index:
            registry_row = known.loc[sequence_sha256]
            if registry_row.sv_id != sv_id or registry_row.sequence != sequence:
                raise SvValidationError(
                    "Refpkg registry composite SV ID disagrees with full sequence hash"
                )
    digest = sha256(Path(registry_parquet).read_bytes()).hexdigest()
    return {"registry_sha256": digest, "known_sv_count": str(len(registry))}
