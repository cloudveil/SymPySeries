from sympy import *
from sympy.abc import n
from utilities import NumericalSeries

numerical_sine = NumericalSeries(Rational(1, 1))
numerical_sine.init(lambda k: I * (I**(k+1)) * (-1 + (-1)**(k+1)) / (2 * factorial(k+1))) 
numerical_cosine = numerical_sine.diff() 

print "Numerical sine and cosine series: note how the first"
print "coefficient in each list corresponds to a different power of q"
print "sine: ", numerical_sine.a[0:9]
print "cosine: ", numerical_cosine.a[0:9]
print

print "these NumericalSeries objects keep track of their powers of q, though,"
print "so adding them is not just elementwise addition:"
print "sine + cosine: ", (numerical_sine + numerical_cosine).a[0:9], "\n"

print "they also check for special cases, and give warnings for them"
print "sin(x) + (sin(x))'': ", (numerical_sine + numerical_sine.diff(2)).a[0:9]
print

print "like the SymbolicSeries, these NumericalSeries also have"
print "well-defined differentiation and multiplication operators"
print "(sin(x)): ", numerical_sine.diff(0).a[0:9]
print "(sin(x))': ", numerical_sine.diff(1).a[0:9]
print "(sin(x))'': ", numerical_sine.diff(2).a[0:9]
print "(sin(x))''': ", numerical_sine.diff(3).a[0:9]
print "(sin(x))'''': ", numerical_sine.diff(4).a[0:9]
print
print "sin(x) * sin(x): ", (numerical_sine * numerical_sine).a[0:9]
print "sin(x) * cos(x): ", (numerical_sine * numerical_cosine).a[0:9]
print "cos(x) * cos(x): ", (numerical_cosine * numerical_cosine).a[0:9]


numerical_tan = NumericalSeries(Rational(1, 1))
numerical_tan.init(lambda k: (I*(2*I)**(1+k)*(-1+(-1)**(1+k))*
                   (-1+2**(2+k))*bernoulli(2+k))/factorial(2+k))
print 
print "Unlike the symbolic series, they also have a division operator"
print "sin(x) / cos(x): ", (numerical_sine / numerical_cosine).a[0:9]  
print "accepted tan(x): ", numerical_tan.a[0:9]

numerical_exp = NumericalSeries(Rational(0, 1))
numerical_exp.init(lambda k: 1/factorial(k))

numerical_test = NumericalSeries(Rational(1, 1))
numerical_test.init(lambda k: ((-1-I)**k+(-1+I)**k+I*(-1-2*I)**k-I*(-1+2*I)**k)/(2*factorial(k)))

f = (((numerical_sine * numerical_sine).diff() + numerical_cosine) / numerical_exp)

print 
print "stress test"
print "((sin(x)*sin(x))' + cos(x)) / exp(x):"
print "accepted values: ", [simplify(a_n) for a_n in numerical_test.a[0:9]]
print "series arithmetic values: ", f.a[0:9]


print
print "mixed length operations test"

sine_10 = NumericalSeries(Rational(1, 1), nmax=10)
sine_10.init(lambda k: im(I**(k+1) / factorial(k+1)))

sine_5 = NumericalSeries(Rational(1, 1), nmax=5)
sine_5.init(lambda k: im(I**(k+1) / factorial(k+1)))

cosine_5 = NumericalSeries(Rational(0, 1), nmax=5)
cosine_5.init(lambda k: re(I**(k) / factorial(k)))

print "sine (10 terms) / cosine (5 terms): ", (sine_10 / cosine_5).a
print 
print "sine (10 terms) - sine (5 terms): ", (sine_10 - sine_5).a
