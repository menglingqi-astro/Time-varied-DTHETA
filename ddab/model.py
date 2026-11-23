"""DD with a Taylor series for the total DTH; no installed PINT code is changed."""

import astropy.units as u
import numpy as np
from pint.exceptions import InvalidModelParameters, MissingParameter
from pint.models.binary_dd import BinaryDD
from pint.models.parameter import MJDParameter, floatParameter
from pint.models.stand_alone_psr_binaries.DD_model import DDmodel


class DDABmodel(DDmodel):
    """DTH(t) = DTH + DTH1*dt + DTH2*dt**2/2, dt in seconds.

    DTH remains the scalar coefficient at DTHEPOCH. dth() evaluates the
    total observable deformation, without separating GR and aberration.
    All other forward calculations are inherited from DD.
    """

    def __init__(self, t=None, input_params=None):
        super().__init__()
        self.binary_name = "DDAB"
        for name, value in {
            "DTHEPOCH": np.longdouble(54000)*u.day,
            "DTH1": np.longdouble(0)/u.s,
            "DTH2": np.longdouble(0)/u.s**2,
        }.items():
            self.add_binary_params(name, value)
        self.add_inter_vars(["dth"])
        if t is not None:
            self.t = np.atleast_1d(t)*u.day if not hasattr(t, "unit") else np.atleast_1d(t)
        if input_params is not None:
            self.update_input(**input_params)

    def dth_dt(self):
        if not hasattr(self, "_t"):
            return np.longdouble(0)*u.s
        return (self.t-self.DTHEPOCH).to(u.s)

    def dth(self):
        dt = self.dth_dt()
        return (self.DTH+self.DTH1*dt+self.DTH2*dt**2/2).decompose()

    def d_dth_d_par(self, par):
        if par not in self.binary_params:
            raise ValueError(f"{par} is not in binary parameter list")
        dt = self.dth_dt()
        if par == "DTH":
            return np.ones_like(dt.value, dtype=np.longdouble)*u.dimensionless_unscaled
        if par == "DTH1":
            return dt
        if par == "DTH2":
            return dt**2/2
        if par == "DTHEPOCH":
            return -(self.DTH1+self.DTH2*dt)*(u.day.to(u.s)*u.s/u.day)
        return np.zeros_like(dt.value, dtype=np.longdouble)/getattr(self, par).unit

    def eTheta(self):
        theta = self.ecc()*(1+self.dth())
        if np.any(~np.isfinite(theta)) or np.any(np.abs(theta) >= 1):
            raise InvalidModelParameters("DDAB requires finite |ecc*(1+DTH(t))| < 1 at every TOA")
        return theta

    def d_eTheta_d_par(self, par):
        return self.prtl_der("ecc", par)*(1+self.dth()) + self.ecc()*self.d_dth_d_par(par)

    def d_er_d_par(self, par):
        return self.prtl_der("ecc", par)*(1+self.DR) + self.ecc()*self.prtl_der("DR", par)

    def d_E_d_par(self, par):
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            return ((self.d_M_d_par(par)+np.sin(self.E())*self.prtl_der("ecc", par))
                    / (1-self.ecc()*np.cos(self.E()))).to(u.rad/getattr(self, par).unit)

    def d_nu_d_par(self, par):
        # The inherited sin(E)/sin(nu) expression is 0/0 at periastron.
        # Use the nonsingular chain rule also used in DDGRAB.
        e, anomaly = self.ecc(), self.E()
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            return (np.sqrt(1-e**2)/(1-e*np.cos(anomaly))*self.d_E_d_par(par)
                    + np.sin(anomaly)/(np.sqrt(1-e**2)*(1-e*np.cos(anomaly)))
                    * self.prtl_der("ecc", par)).to(u.rad/getattr(self, par).unit)

    def prtl_der(self, y, x):
        # Named DD shortcuts for beta/DTH omit e. Others also omit factors
        # in eccentricity derivatives; use the complete product/chain rule,
        # as DDGRAB does. dth's epoch is independent of T0.
        if y in {"eTheta", "beta", "er", "dth", "E", "nu", "omega"} and x in self.binary_params:
            result = getattr(self, f"d_{y}_d_par")(x)
            return result.to(getattr(self, y)().unit/getattr(self, x).unit,
                             equivalencies=u.dimensionless_angles())
        return super().prtl_der(y, x)

    def d_beta_d_DTH(self):
        # Also correct direct callers, beyond prtl_der's dispatch path.
        return self.d_beta_d_par("DTH")


class BinaryDDAB(BinaryDD):
    """Registered BINARY DDAB component retaining DD's independent PK parameters."""

    register = True

    def __init__(self):
        super().__init__()
        self.binary_model_name = "DDAB"
        self.binary_model_class = DDABmodel
        self.DTH.frozen = False
        self.DTH.long_double = True
        self.DTH.description = "Total DTH at DTHEPOCH"
        self.add_param(MJDParameter(
            name="DTHEPOCH", time_scale="tdb", frozen=True,
            description="Fixed reference epoch for the total DTH Taylor series",
            tcb2tdb_scale_factor=u.Quantity(1),
        ))
        for name, unit, description in (
            ("DTH1", "1/s", "First time derivative of total DTH at DTHEPOCH"),
            ("DTH2", "1/s2", "Second time derivative of total DTH at DTHEPOCH"),
        ):
            self.add_param(floatParameter(
                name=name, units=unit, value=0, frozen=False, long_double=True,
                description=description, tcb2tdb_scale_factor=u.Quantity(1),
            ))

    def validate(self):
        super().validate()
        if self.DTHEPOCH.value is None:
            raise MissingParameter("DDAB", "DTHEPOCH", "Use ddab.get_model_and_toas() to choose the TOA-span midpoint.")
        if not self.DTHEPOCH.frozen:
            raise ValueError("DTHEPOCH is a fixed reference epoch; its fit flag must be 0.")
        for name in ("DTHEPOCH", "DTH", "DTH1", "DTH2"):
            if not np.isfinite(getattr(self, name).value):
                raise ValueError(f"{name} must be finite")
