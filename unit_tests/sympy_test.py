from sympy import *

class SigmaTable:

  def __init__(self, max_entries=100):
    self.values = []
    self.max_entries = 0
    self.calculate_up_to(max_entries)

  def calculate_up_to(self, n):
    for i in range(self.max_entries, n):
      self.values.append(sum(divisors(i)))       

  def __getitem__ (self, index):
    return self.values[index]  

class qSeries
 
s = SigmaTable(100) 

for i in range(1, 25):
  print(s[i])
