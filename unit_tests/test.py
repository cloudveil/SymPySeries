from sympy import *
from utilities import NumericalSeries

E4 = NumericalSeries.E(4) 
E6 = NumericalSeries.E(6)
E8 = NumericalSeries.E(8)
E10 = NumericalSeries.E(10)

delta = Rational(1, 1728) * (E4 * E4 * E4 - E6 * E6)

f_0 = NumericalSeries(Rational(1, 12), 50)
