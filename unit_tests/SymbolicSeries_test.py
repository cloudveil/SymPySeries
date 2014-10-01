from sympy import *
from sympy.abc import n
from utilities import SymbolicSeries

sym_series = SymbolicSeries(Rational(1, 6))

print "symbolic q-series, general form:"
print sym_series.to_str(3) + "\n"

print "symbolic q-series, normalized:"
normalized_series = SymbolicSeries(Rational(1, 6))
normalized_series.coef = Piecewise((1, n < 1), (sym_series.coef, True))
print normalized_series.to_str(3) + "\n"

print "symbolic q-series, recurrence for differential equation f'' + E6 * f = 0:"
print (sym_series.diff(2) + SymbolicSeries.E(6) * sym_series).coef

print "\n\n"

# coefficient list courtesy of Mathematica
sine = SymbolicSeries(Rational(1, 1))
sine.coef = I * (I**(n+1)) * (-1 + (-1)**(n+1)) / (2 * factorial(n+1))
print "sine series:"
print sine.to_str(5) + "\n"

cosine = sine.diff(1)
print "cosine series (by series differentiation):"
print cosine.to_str(5) + "\n"

print "series addition, sin(q) + cos(q)"
print (sine + cosine).to_str(5) + "\n"

print "series multiplication, sin(q) * cos(q):"
print (sine * cosine).to_str(5) + "\n"

# coefficient list check for (sine * cosine)(x) also courtesy of Mathematica
sine_cosine_accepted = SymbolicSeries(Rational(1, 1))
sine_cosine_accepted.coef = (I * 2**(-2 + (n+1)) * ((-I)**(n+1) - I**(n+1)))/factorial(n+1)

print "accepted value for series multiplication, sin(q) * cos(q):"
print sine_cosine_accepted.to_str(5) + "\n"

print "combined series arithmetic, sin(q)**2 + cos(q)**2:"
print (sine*sine + cosine*cosine).to_str(5) + "\n"
