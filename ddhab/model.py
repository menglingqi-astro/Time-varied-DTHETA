"""DDH Shapiro parameterization with DDAB's total DTH Taylor series."""
from pint.models.binary_dd import BinaryDDH
from pint.models.stand_alone_psr_binaries.DDH_model import DDHmodel
from pint.exceptions import MissingParameter
from ddab.model import BinaryDDAB, DDABmodel


class DDHABmodel(DDABmodel, DDHmodel):
    """Combine DDAB's deformation/chain rules with native DDH Shapiro delay.

    The cooperative initialization first constructs DDH, then adds DDAB's
    DTHEPOCH/DTH1/DTH2. The shared DD base occurs only once in the MRO.
    H3 and STIGMA determine M2 and SINI as in native DDH.
    """
    def __init__(self, t=None, input_params=None):
        super().__init__(t=t, input_params=input_params)
        self.binary_name = "DDHAB"


class BinaryDDHAB(BinaryDDAB, BinaryDDH):
    """Registered BINARY DDHAB: independent PK parameters and H3/STIGMA."""
    register = True

    def __init__(self):
        super().__init__()
        self.binary_model_name = "DDHAB"
        self.binary_model_class = DDHABmodel

    def validate(self):
        if self.DTHEPOCH.value is None:
            raise MissingParameter("DDHAB", "DTHEPOCH",
                                   "Use ddhab.get_model_and_toas() to choose the TOA-span midpoint.")
        super().validate()
