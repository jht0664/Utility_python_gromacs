#!/usr/bin/env python

# =============================================================================
# MODULE DOCSTRING
# =============================================================================

"""
OEMC integrator for hybrid MC/MD simulations

DESCRIPTION

This module provides OEMC integrators for OpenMM.

EXAMPLES

COPYRIGHT

@author Hyuntae Jung <hjung52@wisc.edu>

All code in this repository is released under the MIT License.

This program is free software: you can redistribute it and/or modify it under
the terms of the MIT License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the MIT License for more details.

You should have received a copy of the MIT License along with this program.

"""

# ============================================================================================
# GLOBAL IMPORTS
# ============================================================================================

import logging

import simtk.unit as unit
import simtk.openmm as mm

from openmm_oemc.constants import 

logger = logging.getLogger(__name__)

# Energy unit used by OpenMM unit system
_OPENMM_ENERGY_UNIT = unit.kilojoules_per_mole

# ============================================================================================
# BASE CLASSES
# ============================================================================================



if __name__ == '__main__':
    import doctest
    doctest.testmod()
