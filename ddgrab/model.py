"""Additive, time-dependent DTH extension of PINT's DDGR binary model.

Importing this module registers BinaryDDGRAB; no existing PINT class is changed.
"""

import astropy.constants as c
import astropy.units as u
import numpy as np
import pint
from pint.exceptions import InvalidModelParameters, MissingParameter
from pint.models.binary_dd import BinaryDDGR
from pint.models.parameter import MJDParameter, floatParameter
from pint.models.stand_alone_psr_binaries.DDGR_model import DDGRmodel


class DDGRABmodel(DDGRmodel):
    """DDGR with DTH = DTH_GR + DTHAB + DTHAB1*dt + DTHAB2*dt**2/2.

    DTH_GR is the scalar computed by the inherited DDGR _updatePK().
    dt uses the barycentric binary input time minus DTHEPOCH, in seconds.
    DTHAB1/2 are actual derivatives in s^-1/s^-2, with no implicit scaling.
    """

    def __init__(self, t=None, input_params=None):
        super().__init__()
        self.binary_name = "DDGRAB"
        for name, value in {
            "DTHEPOCH": np.longdouble(54000) * u.day,
            "DTHAB": np.longdouble(0) * u.dimensionless_unscaled,
            "DTHAB1": np.longdouble(0) / u.s,
            "DTHAB2": np.longdouble(0) / u.s**2,
        }.items():
            self.add_binary_params(name, value)
        if t is not None:
            self.t = np.atleast_1d(t) * u.day if not hasattr(t, "unit") else np.atleast_1d(t)
        if input_params is not None:
            self.update_input(**input_params)
        self._updatePK()

    @property
    def DTH_GR(self):
        """Constant GR contribution for the current masses and PB."""
        return self._DTH

    def dth_dt(self):
        if not hasattr(self, "_t"):
            # DDGR also initializes PK quantities before any TOAs are supplied.
            return np.longdouble(0) * u.s
        return (self.t - self.DTHEPOCH).to(u.s)

    def dth_ab(self):
        dt = self.dth_dt()
        return (self.DTHAB + self.DTHAB1 * dt + self.DTHAB2 * dt**2 / 2).decompose()

    def _updatePK(self, ARTOL=1e-10):
        if getattr(self, "_updating_input", False):
            return
        # PINT sends the same complete binary input once per design-matrix
        # column. Avoid recalculating the PK arrays for unchanged input.
        fields = (("PB", u.day), ("ECC", u.dimensionless_unscaled),
                  ("M2", u.Msun), ("MTOT", u.Msun), ("A1", u.lsec),
                  ("EDOT", 1/u.s), ("A1DOT", u.lsec/u.s), ("T0", u.day))
        ready = all(hasattr(self, name) for name, _ in fields)
        key = tuple(getattr(self, name).to_value(unit) for name, unit in fields)+(ARTOL,) if ready else None
        time = self.t.to_value(u.day) if hasattr(self, "_t") else None
        previous_time = getattr(self, "_pk_time", None)
        same_time = (time is None and previous_time is None) or (
            time is not None and previous_time is not None and np.array_equal(time, previous_time))
        if key is not None and key == getattr(self, "_pk_key", None) and same_time:
            return
        super()._updatePK(ARTOL=ARTOL)
        self._gr_jac = {}
        self._pk_key = key
        self._pk_time = None if time is None else np.array(time, copy=True)

    def update_input(self, **updates):
        # Inherited PB/mass/eccentricity setters normally update PK quantities
        # individually. Defer those calls until the whole input is consistent.
        self._updating_input = True
        try:
            super().update_input(**updates)
        finally:
            self._updating_input = False
        self._gr_jac = {}
        self._updatePK()

    def _identity(self, name, par):
        return (np.longdouble(1 if name == par else 0) * getattr(self, name).unit / getattr(self, par).unit)

    def _gr_derivatives(self, par):
        """Differentiate the inherited PK values, including all indirect paths.

        In particular GAMMA and k use the Newtonian arr0, while DR and DTH
        use the relativistic arr. Differentiating the implicit Kepler equation
        avoids approximating either by the instantaneous pb(). No GR forward
        calculation is replaced; all values still come from DDGR._updatePK.
        """
        if par in self._gr_jac:
            return self._gr_jac[par]
        mass, m2, m1 = self.MTOT, self.M2, self.MTOT - self.M2
        dm = self._identity("MTOT", par)
        dm2 = self._identity("M2", par)
        dm1 = dm - dm2
        dpb = self._identity("PB", par)
        de = super().prtl_der("ecc", par)
        dx = self.d_a1_d_par(par)
        e = self.ecc()
        dln_a0 = (dm / mass + 2 * dpb / self.PB) / 3
        q = m1 * m2 / mass**2 - 9
        dq = (dm1*m2 + m1*dm2) / mass**2 - 2*m1*m2/mass**2 * dm/mass
        factor = pint.GMsun * self.mtot / (2*self.arr*c.c**2)
        z = q * factor
        dln_a = (dln_a0 + (2/3)*factor*(dq + q*dm/mass)/(1+z)) / (1 + (2/3)*z/(1+z))
        d_sini = (self.SINI * (dm/mass - dm2/m2 - dln_a) + mass/(self.arr*m2)*dx)
        # Express de terms without division by eccentricity, so e=0 works.
        arr0 = (pint.GMsun*self.mtot/self._n**2)**(1/3)
        gamma_per_e = (pint.GMsun*self.m2*(self.mtot+self.m2)/(self._n*c.c**2*arr0*self.mtot)).to(u.s)
        d_gamma = gamma_per_e*de + self.GAMMA*(dm2/m2 + (dm+dm2)/(mass+m2) - dm/mass + dpb/self.PB - dln_a0)
        poly = 1 + (73/24)*e**2 + (37/96)*e**4
        dln_fe = (((73/12)*e + (37/24)*e**3)/poly + 7*e/(1-e**2))*de
        d_pbdot = self.PBDOT*(dm1/m1 + dm2/m2 - dm/(3*mass) - (5/3)*dpb/self.PB + dln_fe)
        d_k = self.k*(dm/mass - dln_a0 + 2*e*de/(1-e**2))
        nr = 3*m1**2 + 6*m1*m2 + 2*m2**2
        nt = 3.5*m1**2 + 6*m1*m2 + 2*m2**2
        dnr = (6*m1+6*m2)*dm1 + (6*m1+4*m2)*dm2
        dnt = (7*m1+6*m2)*dm1 + (6*m1+4*m2)*dm2
        d_dr = self.DR*(dnr/nr - dm/mass - dln_a)
        d_dth = self.DTH_GR*(dnt/nt - dm/mass - dln_a)
        if par == "DTHAB":
            d_dth = self.d_DTH_d_DTHAB()
        elif par == "DTHAB1":
            d_dth = self.d_DTH_d_DTHAB1()
        elif par == "DTHAB2":
            d_dth = self.d_DTH_d_DTHAB2()
        elif par == "DTHEPOCH":
            d_dth = self.d_DTH_d_DTHEPOCH()
        result = {"SINI": d_sini, "GAMMA": d_gamma, "PBDOT": d_pbdot,
                  "k": d_k, "DR": d_dr, "DTH": d_dth}
        self._gr_jac[par] = result
        return result

    @property
    def DTH(self):
        return self.DTH_GR + self.dth_ab()

    @DTH.setter
    def DTH(self, value):
        # Preserve DDGR's treatment of supplied derived DTH values.
        DDGRmodel.DTH.fset(self, value)

    def eTheta(self):
        # Do not use DDGR's cached _eth, which contains only the GR term.
        theta = self.ecc() * (1 + self.DTH)
        if np.any(np.abs(theta) >= 1):
            raise InvalidModelParameters("DDGRAB requires |ecc*(1+DTH_obs)| < 1 at every TOA")
        return theta

    def d_DTH_d_DTHAB(self):
        return np.ones_like(self.dth_dt().value, dtype=np.longdouble) * u.dimensionless_unscaled

    def d_DTH_d_DTHAB1(self):
        return self.dth_dt()

    def d_DTH_d_DTHAB2(self):
        return self.dth_dt()**2 / 2

    def d_DTH_d_DTHEPOCH(self):
        return -(self.DTHAB1 + self.DTHAB2 * self.dth_dt()) * (u.day.to(u.s) * u.s / u.day)

    # DDGR's formulas refer to self.DTH. Here they must use DTH_GR only:
    # the additive term has no mass or PB dependence at fixed input time.
    def d_DTH_d_MTOT(self):
        return self._gr_derivatives("MTOT")["DTH"]

    def d_DTH_d_M2(self):
        return self._gr_derivatives("M2")["DTH"]

    def d_DTH_d_PB(self):
        return self._gr_derivatives("PB")["DTH"]

    def d_DTH_d_par(self, par):
        return self._gr_derivatives(par)["DTH"]

    def d_er_d_par(self, par):
        return super().prtl_der("ecc", par)*(1+self.DR) + self.ecc()*self._gr_derivatives(par)["DR"]

    def d_eTheta_d_par(self, par):
        return self.prtl_der("ecc", par) * (1 + self.DTH) + self.ecc() * self.d_DTH_d_par(par)

    def d_M_d_par(self, par):
        if self.orbits_cls.name != "orbitPB":
            return super().d_M_d_par(par)
        pb = self.PB.to(u.s)
        dt = self.tt0
        ddt = -self._identity("T0", par).to(u.s/getattr(self, par).unit)
        dpb = self._identity("PB", par).to(u.s/getattr(self, par).unit)
        dratio = ddt/pb - dt*dpb/pb**2
        ratio = dt/pb
        ddecay = self._gr_derivatives(par)["PBDOT"] + self._identity("XPBDOT", par)
        return 2*np.pi*u.rad*(dratio*(1-(self.PBDOT+self.XPBDOT)*ratio) - ddecay*ratio**2/2)

    def d_pb_d_par(self, par):
        if self.orbits_cls.name != "orbitPB":
            return super().d_pb_d_par(par)
        return (self._identity("PB", par) + self._gr_derivatives(par)["PBDOT"]*self.tt0
                - self.PBDOT*self._identity("T0", par))

    def d_E_d_par(self, par):
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            return ((self.d_M_d_par(par) + np.sin(self.E())*super().prtl_der("ecc", par))
                    / (1-self.ecc()*np.cos(self.E()))).to(u.rad/getattr(self, par).unit)

    def d_nu_d_par(self, par):
        e, anomaly = self.ecc(), self.E()
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            return (np.sqrt(1-e**2)/(1-e*np.cos(anomaly))*self.d_E_d_par(par)
                    + np.sin(anomaly)/(np.sqrt(1-e**2)*(1-e*np.cos(anomaly)))
                    * super().prtl_der("ecc", par)).to(u.rad/getattr(self, par).unit)

    def d_omega_d_par(self, par):
        excess = self.XOMDOT/(self._n*u.rad)
        dexcess = self._identity("XOMDOT", par)/(self._n*u.rad) + excess*self._identity("PB", par)/self.PB
        return (self._identity("OM", par) + self.d_nu_d_par(par)*(self.k+excess)
                + self.nu()*(self._gr_derivatives(par)["k"]+dexcess))

    def d_nhat_d_par(self, par):
        e, anomaly = self.ecc(), self.E()
        with u.set_enabled_equivalencies(u.dimensionless_angles()):
            return self.nhat()*(-self.d_pb_d_par(par)/self.pb()
                + (super().prtl_der("ecc", par)*np.cos(anomaly)
                   - e*np.sin(anomaly)*self.d_E_d_par(par))/(1-e*np.cos(anomaly)))

    def prtl_der(self, y, x):
        # Bypass inherited named shortcuts that omit DDGR's indirect paths.
        if y in {"SINI", "GAMMA", "PBDOT", "DR", "DTH", "k"}:
            result = self._gr_derivatives(x)[y]
        elif y in {"er", "eTheta", "E", "nu", "omega", "beta"}:
            result = getattr(self, f"d_{y}_d_par")(x)
        else:
            return super().prtl_der(y, x)
        value = getattr(self, y)
        unit = value().unit if callable(value) else value.unit
        return result.to(unit/getattr(self, x).unit, equivalencies=u.dimensionless_angles())

    # prtl_der prefers named methods over the general dispatcher.
    def d_eTheta_d_MTOT(self):
        return self.d_eTheta_d_par("MTOT")

    def d_eTheta_d_M2(self):
        return self.d_eTheta_d_par("M2")

    def d_beta_d_MTOT(self):
        return self.d_beta_d_par("MTOT")

    def d_beta_d_M2(self):
        return self.d_beta_d_par("M2")


class BinaryDDGRAB(BinaryDDGR):
    """PINT component for BINARY DDGRAB, retaining DDGR's parameter set."""

    register = True

    def __init__(self):
        super().__init__()
        self.binary_model_name = "DDGRAB"
        self.binary_model_class = DDGRABmodel
        self.add_param(MJDParameter(
            name="DTHEPOCH", time_scale="tdb", frozen=True,
            description="Fixed reference epoch for the DTHAB Taylor series",
            tcb2tdb_scale_factor=u.Quantity(1),
        ))
        for name, unit, description in (
            ("DTHAB", "", "Additive DTH at DTHEPOCH"),
            ("DTHAB1", "1/s", "First time derivative of additive DTH at DTHEPOCH"),
            ("DTHAB2", "1/s2", "Second time derivative of additive DTH at DTHEPOCH"),
        ):
            self.add_param(floatParameter(
                name=name, units=unit, value=0, frozen=False, long_double=True,
                description=description, tcb2tdb_scale_factor=u.Quantity(1),
            ))

    def validate(self):
        self.check_required_params(["PB", "MTOT", "M2"])
        if not 0 < self.M2.value < self.MTOT.value:
            raise ValueError("DDGRAB requires 0 < M2 < MTOT")
        super().validate()
        if self.DTHEPOCH.value is None:
            raise MissingParameter("DDGRAB", "DTHEPOCH", "Use ddgrab.get_model_and_toas() to choose the TOA-span midpoint.")
        if not self.DTHEPOCH.frozen:
            raise ValueError("DTHEPOCH is a fixed reference epoch; its fit flag must be 0.")
        for name in ("DTHEPOCH", "DTHAB", "DTHAB1", "DTHAB2"):
            if not np.isfinite(getattr(self, name).value):
                raise ValueError(f"{name} must be finite")
