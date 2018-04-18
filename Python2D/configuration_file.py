from __future__ import division
from enum import Enum

import numpy as np
import sympy as sp


n_elem = 7  		# number of elements in x and ydirection
init_temp = 20 		# initial temp
jump_temp = 0		# spike temp
bound_temp = 0		# boundary temp
steps = 100			# total steps
dt = 0.01			# dt (stepsize)
scalingFactor = 1	# scaling factor
# style = "explicit"
style = "implicit"
