class ns_rng:
   def int_uniform(self, low, high):
      raise("int_uniform not implemented")
   def float_uniform(self, low, high, size=None):
      raise("float_uniform not implemented")
   def normal(self, std_dev):
      raise("normal not implemented")
   def shuffle_in_place(self, list):
      raise("shuffle_in_place not implemented")

import numpy as np

class ns_rng_numpy(ns_rng):
   def int_uniform(self, low, high):
      return np.random.randint(low,high)
   def float_uniform(self, low, high, size=None):
      return np.random.uniform(low, high, size)
   def normal(self, std_dev):
      return np.random.normal(0.0,std_dev)
   def shuffle_in_place(self, list):
      np.random.shuffle(list)

class ns_rng_julia(ns_rng):
   def __init__(self, j):
      self.j=j
   def int_uniform(self, low, high):
      return self.j.rand(xrange(high-low))+low
   def float_uniform(self, low, high, size=None):
      if size is None:
	 return self.j.rand()*(high-low)+low
      else:
	 return self.j.rand(size)*(high-low)+low
   # def normal(std_dev):
   #
   def shuffle_in_place(self, list):
      self.j.shuffle_b(list)

class ns_rng_internal(ns_rng):
   def int_uniform(self, low, high):
      return random.randint(low,high-1)
   def float_uniform(self, low,high, size=None):
      if size is not None:
	 sys.stderr.write("ns_rng_internal size not implemented for float_uniform\n")
	 sys.exit(1)
      return random.uniform(low,high)
   def normal(self, std_dev):
      return random.normalvariate(0.0, std_dev)
   def shuffle_in_place(self, list):
      random.shuffle(list)
