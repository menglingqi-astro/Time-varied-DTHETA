"""In-memory DDH conversion and TOA-dependent reference epoch selection."""

from io import StringIO
from pathlib import Path

import astropy.units as u
import numpy as np
from pint.models import get_model_and_toas as pint_get_model_and_toas


def get_model_and_toas(parfile, timfile, *, convert_ddh=True, **kwargs):
    """Load DDHAB and resolve a missing DTHEPOCH from the retained TOA span.

    Accepts paths or file-like inputs, and forwards PINT loading options.
    With convert_ddh=True, BINARY DDH is changed to DDHAB in memory.
    Existing DTHEPOCH and coefficient values/fit flags are preserved.
    No par/tim file is written. PINT's TOA caching is disabled by default.
    The midpoint is in barycentric TDB, matching the binary input time.
    Save model.as_parfile() to persist the selected epoch and coefficients.
    """
    text = parfile.read() if hasattr(parfile, "read") else Path(parfile).read_text()
    lines = text.splitlines()
    binary_rows = [i for i, line in enumerate(lines) if line.split() and line.split()[0].upper() == "BINARY"]
    if len(binary_rows) != 1:
        raise ValueError("Expected exactly one active BINARY line")
    index = binary_rows[0]
    binary = lines[index].split()[1].upper()
    if binary in ("DDH", "DDFWHE") and convert_ddh:
        lines[index] = "BINARY DDHAB"
    elif binary != "DDHAB":
        raise ValueError("Expected BINARY DDHAB (or DDH with convert_ddh=True)")
    if kwargs.get("force_binary_model") not in (None, "DDHAB"):
        raise ValueError("The DDHAB loader cannot force a different binary model")

    has_epoch = any(line.split() and line.split()[0].upper() == "DTHEPOCH" for line in lines)
    has_epoch = has_epoch or kwargs.get("DTHEPOCH") is not None
    if not has_epoch:
        # Temporary in-memory epoch permits PINT to build the model before TOAs
        # are loaded. No residual calculation or fitting uses this placeholder.
        lines.append("DTHEPOCH 54000 0")
    add_tzr = kwargs.pop("add_tzr_to_model", True)
    kwargs.setdefault("usepickle", False)
    model, toas = pint_get_model_and_toas(
        StringIO("\n".join(lines) + "\n"), timfile, add_tzr_to_model=False, **kwargs,
    )
    if not has_epoch:
        times = np.asarray(model.get_barycentric_toas(toas).to_value(u.day), dtype=np.longdouble)
        if times.size == 0 or not np.all(np.isfinite(times)):
            raise ValueError("Cannot select DTHEPOCH from empty or non-finite TOAs")
        model.DTHEPOCH.value = (times.min() + times.max()) / 2
        model.DTHEPOCH.frozen = True
    model.meta["ddhab_epoch_source"] = "input" if has_epoch else "barycentric_toa_span_midpoint"
    model.validate()
    model.components["BinaryDDHAB"].update_binary_object(toas)
    if add_tzr and "AbsPhase" not in model.components:
        model.add_tzr_toa(toas)
    return model, toas
