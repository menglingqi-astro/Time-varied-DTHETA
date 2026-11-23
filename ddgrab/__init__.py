"""Import to register DDGRAB with PINT, without altering DDGR."""

from .model import BinaryDDGRAB, DDGRABmodel
from .loader import get_model_and_toas

__all__ = ["BinaryDDGRAB", "DDGRABmodel", "get_model_and_toas"]
