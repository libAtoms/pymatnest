class NsRng:
   def int_uniform(self, low, high):
      raise("int_uniform not implemented")
   def float_uniform(self, low, high, size=None):
      raise("float_uniform not implemented")
   def normal(self, std_dev):
      raise("normal not implemented")
   def shuffle_in_place(self, list):
      raise("shuffle_in_place not implemented")
   def seed_init(self):
      raise("seed_init not implemented")
   def seed_common(self):
      raise("seed_common not implemented")
   def seed_local(self):
      raise("seed_local not implemented")

import numpy as np
class NsRngNumpy(NsRng):
   def int_uniform(self, low, high):
      return np.random.randint(low,high)
   def float_uniform(self, low, high, size=None):
      return np.random.uniform(low, high, size)
   def normal(self, std_dev):
      return np.random.normal(0.0,std_dev)
   def shuffle_in_place(self, list):
      np.random.shuffle(list)
   def init_seeds(self, delta_seed, comm):
      if comm is None:
	 self.common_random_state = None
      else:
	 if rank == 0:
	    if delta_seed > 0:
	       np.random.seed(comm.size + 1 + delta_seed)
	    else:
	       np.random.seed()
	    self.common_random_state = np.random.get_state()
	 else:
	    self.common_random_state = None
	 self.common_random_state = comm.bcast(self.common_random_state, root=0)
      if delta_seed > 0:
	 np.random.seed(comm.rank + 1 + delta_seed)
      else:
	 np.random.seed()
   def switch_to_common(self):
      if self.common_random_state is not None:
	 self.local_random_state = np.random.get_state()
	 np.random.set_state(self.common_random_state)
   def switch_to_local(self):
      if self.common_random_state is not None:
	 self.common_random_state = np.random.get_state()
	 np.random.set_state(self.local_random_state)

class NsRngJulia(NsRng):
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
      list[:] = self.j.shuffle(list)

import random, sys
class NsRngInternal(NsRng):
   def int_uniform(self, low, high):
      return random.randint(low,high-1)
   def float_uniform(self, low,high, size=None):
      if size is not None:
	 out = np.zeros (size)
	 for x in np.nditer(out, op_flags=['readwrite']):
	    x[...] = random.uniform(low,high)
	 return out
      else:
	 return random.uniform(low,high)
   def normal(self, std_dev):
      return random.normalvariate(0.0, std_dev)
   def shuffle_in_place(self, list):
      random.shuffle(list)
   def init_seeds(self, delta_seed, comm):
      if comm is None:
	 self.common_random_state = None
      else:
	 if rank == 0:
	    if delta_seed > 0:
	       random.seed(comm.size + 1 + delta_seed)
	    else:
	       random.seed()
	    self.common_random_state = random.get_state()
	 else:
	    self.common_random_state = None
	 self.common_random_state = comm.bcast(self.common_random_state, root=0)
      if delta_seed > 0:
	 random.seed(comm.rank + 1 + delta_seed)
      else:
	 random.seed()
   def switch_to_common(self):
      if self.common_random_state is not None:
	 self.local_random_state = random.get_state()
	 random.set_state(self.common_random_state)
   def switch_to_local(self):
      if self.common_random_state is not None:
	 self.common_random_state = random.get_state()
	 random.set_state(self.local_random_state)
