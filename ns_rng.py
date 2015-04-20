class NsRng:
   def int_uniform(self, low, high):
      raise("int_uniform not implemented")
   def float_uniform(self, low, high, size=None):
      raise("float_uniform not implemented")
   def normal(self, std_dev):
      raise("normal not implemented")
   def shuffle_in_place(self, list):
      raise("shuffle_in_place not implemented")
   def seed_init(self, delta_seed=-1, comm=None):
      raise("seed_init not implemented")
   def seed_common(self):
      raise("seed_common not implemented")
   def seed_local(self):
      raise("seed_local not implemented")

import numpy as np
class NsRngNumpy(NsRng):
   def __init__(self, delta_seed=-1, comm=None):
      if comm is None:
	 self.common_random_state = None
      else:
	 if comm.rank == 0:
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
   def int_uniform(self, low, high):
      return np.random.randint(low,high)
   def float_uniform(self, low, high, size=None):
      return np.random.uniform(low, high, size)
   def normal(self, std_dev):
      return np.random.normal(0.0,std_dev)
   def shuffle_in_place(self, list):
      np.random.shuffle(list)

   def switch_to_common(self):
      if self.common_random_state is not None:
	 self.local_random_state = np.random.get_state()
	 np.random.set_state(self.common_random_state)
   def switch_to_local(self):
      if self.common_random_state is not None:
	 self.common_random_state = np.random.get_state()
	 np.random.set_state(self.local_random_state)

#class NsRngJulia(NsRng):
#   def __init__(self, j):
#      self.j=j
#   def int_uniform(self, low, high):
#      return self.j.rand(xrange(high-low))+low
#   def float_uniform(self, low, high, size=None):
#      if size is None:
#	 return self.j.rand()*(high-low)+low
#      else:
#	 return self.j.rand(size)*(high-low)+low
#   # def normal(std_dev):
#   #
#   def shuffle_in_place(self, list):
#      list[:] = self.j.shuffle(list)

import random, sys
class NsRngInternal(NsRng):
   def __init__(self, delta_seed=-1, comm=None):
      if comm is None:
	 self.common_random_state = None
      else:
	 if comm.rank == 0:
	    if delta_seed > 0:
	       random.seed(comm.size + 1 + delta_seed)
	    else:
	       random.seed()
	    self.common_random_state = random.getstate()
	 else:
	    self.common_random_state = None
	 self.common_random_state = comm.bcast(self.common_random_state, root=0)
      if delta_seed > 0:
	 random.seed(comm.rank + 1 + delta_seed)
      else:
	 random.seed()
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
   def switch_to_common(self):
      if self.common_random_state is not None:
	 self.local_random_state = random.getstate()
	 random.setstate(self.common_random_state)
   def switch_to_local(self):
      if self.common_random_state is not None:
	 self.common_random_state = random.getstate()
	 random.setstate(self.local_random_state)

def unpack_uint32(bytes):
   tup = struct.unpack("I", bytes)
   return tup[0]

# import rngstream
import rngstream, ctypes, os, struct
class NsRngStream(NsRng):

   def set_package_seed(self, delta_seed):
      seed = []
      for i in range(6):
	 if delta_seed > 0:
	    seed.append(delta_seed)
	 else:
	    rand_bytes = os.urandom(4)
	    seed.append(unpack_uint32(rand_bytes))
      self.r.set_package_seed((ctypes.c_ulong*6)(*seed))

   def __init__(self, delta_seed=-1, comm=None):
      self.r = rngstream.RngStream()

      self.set_package_seed(delta_seed)
      rng = self.r.create_stream()
      if comm is None:
	 self.rng = rng
	 self.l_rng = None
	 self.g_rng = None
      else:
	 seed = (ctypes.c_ulong * 6)()
	 seed_new = (ctypes.c_ulong * 6)()
	 if comm.rank == 0:
	    self.g_rng = rng
	    self.r.get_state(self.g_rng,seed)
	    seed_comm = list(seed)
	    seed_comm = comm.bcast(seed_comm, root=0)
	 else:
	    seed_comm = None
	    seed_comm = comm.bcast(seed_comm, root=0)
	    seed = (ctypes.c_ulong * 6)(*seed_comm)
	    self.r.set_package_seed(seed)
	    self.g_rng = self.r.create_stream()
	 if comm.rank == 0:
	    for i in range(comm.size):
	       t_rng = self.r.create_stream()
	       self.r.get_state(t_rng,seed_new)
	       if i == 0:
		  self.r.set_package_seed(seed_new)
		  self.l_rng = self.r.create_stream()
	       else:
		  seed_new_comm = list(seed_new)
		  comm.send(seed_new_comm, i, tag=1)
	 else:
	    seed_new_comm = comm.recv(source=0, tag=1)
	    seed_new = (ctypes.c_ulong * 6)(*seed_new_comm)
	    self.r.set_package_seed(seed_new)
	    self.l_rng = self.r.create_stream()
	 self.rng = self.l_rng
      self.r.get_state(self.l_rng,seed)
      print comm.rank, "end of init l_rng ", list(seed)
      self.r.get_state(self.g_rng,seed)
      print comm.rank, "end of init g_rng ", list(seed)

   def switch_to_common(self):
      if self.g_rng is not None:
	 self.rng = self.g_rng
   def switch_to_local(self):
      if self.l_rng is not None:
	 self.rng = self.l_rng

   def int_uniform(self, low, high):
      return self.r.int_uniform(self.rng, ctypes.c_int(low), ctypes.c_int(high-1) )
   def float_uniform(self, low, high, size=None):
      if size is None:
	 return self.r.float_uniform_01(self.rng) * (high-low)+low
      else:
	 out = np.zeros(size)
	 for x in np.nditer(out, op_flags=['readwrite']):
	    x[...] = self.r.float_uniform_01(self.rng) * (high-low)+low
	 return out
   # def normal(self, std_dev):
   def shuffle_in_place(self, list):
      random.shuffle(list, lambda : self.r.float_uniform_01(self.rng) )
