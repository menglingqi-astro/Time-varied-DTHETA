# ddab

This package provides a separately registered PINT binary model. See the main
[README](../README.md#models-and-parameter-conventions) for its mathematical
definition, parameter units, epoch convention, fit flags, and scientific scope.

From the repository root:

```python
from ddab import get_model_and_toas
model, toas = get_model_and_toas(
    "/path/to/source.par", "/path/to/source.tim",
    planets=False, usepickle=False,
)
```

The loader modifies the model in memory and preserves the input files. Supply
a fixed DTHEPOCH or let the TOA-aware loader select the barycentric span midpoint.
Import this extension before loading its model through PINT's native par-only API.
The installed PINT package and original model registrations remain intact.
