from sympy import *
from sympy.abc import n
import time

class Timer:    
  def __enter__(self):
    self.start = time.clock()
    return self

  def __exit__(self, *args):
    self.end = time.clock()
    self.interval = self.end - self.start

def sigma(n, k):
  return sum([d**k for d in divisors(n)])

class SymbolicSeries:

  # for a series of the form
  #
  #                   nmax
  #                  _____
  #                  \    `
  #         r         \          (r + n)
  #  f  =  q     +    /    a(n) q
  #                  /____,
  #                  n = 1
  #
  # where r and a(n) are rational numbers 

  def __init__(self, r, coef="a"):
    self.r = sympify(r)
    self.rank = 0
    self.coef_char = coef
    self.coef = Function(coef)(n)

  def diff(self, m=1):
    f = SymbolicSeries(self.r - m, self.coef_char)
    f.coef = self.coef * ff(self.r + n, m)
    f.rank = self.rank
    f.check_for_leading_zeroes()
    return f
  
  # Assuming the following definition of
  # modular derivative:
  #
  #           df     k
  # Dk f = q ---- - ---- E2 f
  #           dq     12
  #
  def mdiff(self, k, n):

    i = 0
    f = SymbolicSeries(self.r, self.coef_char)
    f.coef = self.coef

    E2 = SymbolicSeries.E(2)

    while(i < n):
      q_dfdq = f.diff() 
      q_dfdq.r += 1
      f = q_dfdq - Rational(k + 2*i, 12) * E2 * f
      i += 1

    return f

  def __add__(self, other):

    # first make sure the other
    # operand is compatible
    if isinstance(other, SymbolicSeries):

      # also check that the power series
      # "overlap" to some extent
      if(sympify(other.r - self.r).is_integer):

        f = SymbolicSeries(self.r if other.r >= self.r else other.r) 
        f.coef_char = self.coef_char
        delta = abs(other.r - self.r)
  
        if(other.r < self.r):
          f.coef = Piecewise((other.coef, n < delta),
                             (self.coef.subs(n, n-delta) + other.coef, True))
        elif(other.r > self.r):
          f.coef = Piecewise((self.coef, n < delta),
                             (self.coef + other.coef.subs(n, n-delta), True))
        else:
          f.coef = self.coef + other.coef
          
        f.check_for_leading_zeroes()       
        f.rank = max(self.rank, other.rank)
        return f 

      else:
        print "base exponent incompatible"

    else:
      print "objects incompatible"

  def __sub__(self, other):

    # first make sure the other
    # operand is compatible
    if isinstance(other, SymbolicSeries):

      # also check that the power series
      # "overlap" to some extent
      if(sympify(other.r - self.r).is_integer):

        f = SymbolicSeries(self.r if other.r >= self.r else other.r) 
        f.coef_char = self.coef_char
        delta = abs(other.r - self.r)
  
        if(other.r < self.r):
          f.coef = Piecewise((-other.coef, n < delta),
                             (self.coef.subs(n, n-delta) - other.coef, True))
        elif(other.r > self.r):
          f.coef = Piecewise((self.coef, n < delta),
                             (self.coef - other.coef.subs(n, n-delta), True))
        else:
          f.coef = self.coef - other.coef
          
        f.check_for_leading_zeroes()
        f.rank = max(self.rank, other.rank)
        return f 

      else:
        print "base exponent incompatible"

    else:
      print "objects incompatible"


  def __rmul__(self, other):
    
    if(sympify(other).is_rational):
      f = SymbolicSeries(self.r, self.coef_char)
      f.coef = self.coef * other
      f.rank = self.rank
      return f

  def __mul__(self, other):

    # check and make sure that the second
    # operand is compatible
    if isinstance(other, SymbolicSeries):
      rank = max(self.rank, other.rank) + 1
      j = Symbol("j_" + str(rank))
      f = SymbolicSeries(self.r + other.r, self.coef_char)
      f.coef = Sum(self.coef.subs(n, j) * 
                   other.coef.subs(n, n-j), (j, 0, n))
      f.rank = rank
      return f
   
  def check_for_leading_zeroes(self, i_max = 100):

    # in the event that the leading terms
    # vanish identically, look for the
    # smallest nonzero power of q to make
    # the new leading term

    i = 0
    while(self.coef.subs(n, i).doit() == 0):
      i = i + 1
      if(i > i_max):
        print "warning: first " + str(i_max) + " terms are zero" 
        break

    if(i <= i_max and i != 0):
      self.r += i 
      self.coef = self.coef.subs(n, n+i).doit() 

  def offset(self, r):
    if(sympify(r).is_rational):
      f = SymbolicSeries(self.r + r, self.coef_char)
      f.rank = self.rank
      f.coef = self.coef 
      return f

  def to_str(self, num_terms):

    my_str = ""
    for i in range(0, num_terms):
      if(self.coef.subs(n, i).doit() != 0):
        my_str += "("
        my_str += str(self.coef.subs(n, i).doit())
        my_str += ")*(q**("
        my_str += str(self.r + i)
        my_str += "))"
        my_str += " + "

    my_str += "..."

    return my_str

  @staticmethod
  def E(k):
    summation_coef = 2 * k / bernoulli(k)
    Ek = SymbolicSeries(Rational(0, 1))
    Ek.coef = Piecewise((1, n == 0), 
                        (2 * k * Function("sigma")(n, k-1) / bernoulli(k), True)) 
    return Ek  

class NumericalSeries:

  # for a series of the form
  #
  #                   nmax
  #                  _____
  #                  \    `
  #         r         \          (r + n)
  #  f  =  q     +    /    a(n) q
  #                  /____,
  #                  n = 1
  #
  # where r and a(n) are rational numbers 

  def __init__(self, r, nmax=100):
    self.r = r 
    self.nmax = nmax 
    self.a = [0] * nmax
    self.a[0] = 1
   
  def init(self, f):
    self.a = [f(n) for n in range(0, self.nmax)]  
    self.check_for_leading_zeroes()
 
  def __mul__(self, other):

    # check and make sure that the second
    # operand is compatible
    if isinstance(other, NumericalSeries):
      nmax = min(self.nmax, other.nmax)
      f = NumericalSeries(self.r + other.r, nmax)
      for i in range(0, nmax):
        f.a[i] = 0
        for j in range(0, i+1):
          f.a[i] += self.a[j] * other.a[i-j]
      return f
    elif(sympify(other).is_rational):
      f = NumericalSeries(self.r, self.nmax)
      f.a = [other * a_n for a_n in self.a]
      return f

  def __add__(self, other):

    # check and make sure that the second
    # operand is compatible
    if isinstance(other, NumericalSeries):

      # also check that the power series
      # "overlap" to some extent
      if(sympify(other.r - self.r).is_integer):

        # offset will be used to describe
        # which indices from each series 
        # share a common power of q
        r = min(self.r, other.r)
        offset = abs(other.r - self.r)
        nmax = min(self.r + self.nmax, other.r + other.nmax) - r

        f = NumericalSeries(r, nmax)

        for i in range(0, nmax):
          if(self.r < other.r):
            f.a[i] = self.a[i] if i < offset else self.a[i] + other.a[i-offset]       
          elif(self.r > other.r):
            f.a[i] = other.a[i] if i < offset else other.a[i] + self.a[i-offset]       
          else:
            f.a[i] = self.a[i] + other.a[i]       

        f.check_for_leading_zeroes() 
        return f

  def __sub__(self, other):
    
    # check and make sure that the second
    # operand is compatible
    if isinstance(other, NumericalSeries):
      nmax = self.nmax

      # also check that the power series
      # "overlap" to some extent
      if(sympify(other.r - self.r).is_integer):

        # offset will be used to describe
        # which indices from each series 
        # share a common power of q
        r = min(self.r, other.r)
        offset = abs(other.r - self.r)
        nmax = min(self.r + self.nmax, other.r + other.nmax) - r

        f = NumericalSeries(r, nmax)

        for i in range(0, nmax):
          if(self.r < other.r):
            f.a[i] = self.a[i] if i < offset else self.a[i] - other.a[i-offset]
          elif(self.r > other.r):
            f.a[i] = -other.a[i] if i < offset else -other.a[i] + self.a[i-offset]
          else:
            f.a[i] = self.a[i] - other.a[i]

        f.check_for_leading_zeroes() 
        return f

  def __rmul__(self, other):
    if(sympify(other).is_rational):
      f = NumericalSeries(self.r, self.nmax)
      f.a = [other * a_n for a_n in self.a]
      return f

  def __div__(self, other):

    # check and make sure that the second
    # operand is compatible
    if isinstance(other, NumericalSeries):

      nmax = min(self.nmax, other.nmax)

      f = NumericalSeries(self.r - other.r, nmax)

      # make a copy of the series coefficients
      # to modify during long division
      tmp = self.a[0:nmax] 

      for i in range(0, nmax):
        f.a[i] = tmp[i] / other.a[0]        

        for j in range(i, nmax):
          tmp[j] -= f.a[i] * other.a[j-i] 
      
      return f

  def mdiff(self, k, n=1):

    i = 0
    f = NumericalSeries(self.r, self.nmax)
    f.a = self.a[:]

    E2 = NumericalSeries.E(2, self.nmax)

    while(i < n):
      q_dfdq = f.diff() 
      q_dfdq.r += 1
      f = q_dfdq - Rational(k + 2*i, 12) * E2 * f
      i += 1
    
    return f

  def diff(self, n=1):
    f = NumericalSeries(self.r - n, self.nmax) # is nmax appropriate?
    f.a = [ff(m + self.r, n) * a_m for (m, a_m) in enumerate(self.a)]
    f.check_for_leading_zeroes()
    return f 

  def check_for_leading_zeroes(self):

    # in the event that the leading term
    # vanishes identically, look for the
    # smallest nonzero power of q to make
    # the new leading term
    error = True

    for i in range(0, self.nmax):
      if(self.a[0] == 0):
        self.a.pop(0) 
        self.nmax -= 1
        self.r += 1 
      else:
        error = False
        break

    if(error):
      print "warning: All terms have vanished identically"

  def nb_output(self, num_terms=5):
    q = symbols("q")
    return sum([self.a[i] * q**(i + self.r) for i in range(0, num_terms)])
 
  @staticmethod
  def E(k, nmax=100):
    summation_coef = - 2 * k / bernoulli(k)
    Ek = NumericalSeries(Rational(0, 1), nmax)

    Ek.a[0] = 1
    for i in range(1, nmax):
      Ek.a[i] = summation_coef * sigma(i, k-1) 

    return Ek
