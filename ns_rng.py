class NsRng:
   def int_uniform(self, low, high):
      raise("int_uniform not implemented")
   def float_uniform(self, low, high, size=None):
      raise("float_uniform not implemented")
   def normal(self, std_dev, size=None):
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
   def normal(self, std_dev, size=None):
      return np.random.normal(0.0,std_dev, size)
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
      if size is None:
	 return random.uniform(low,high)
      else:
	 out = np.zeros (size)
	 for x in np.nditer(out, op_flags=['readwrite']):
	    x[...] = random.uniform(low,high)
	 return out
   def normal(self, std_dev, size=None):
      if size is None:
	 return random.normalvariate(0.0, std_dev)
      else:
	 out = np.zeros (size)
	 for x in np.nditer(out, op_flags=['readwrite']):
	    x[...] = random.normalvariate(0.0, std_dev)
	 return out
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
import rngstream, ctypes, os, struct, math
class NsRngStream(NsRng):

   def set_package_seed(self, delta_seed):
      seed = []
      for i in range(6):
	 if delta_seed > 0:
	    seed.append(delta_seed)
	 else:
	    i_rv = unpack_uint32(os.urandom(4))
	    while i_rv == 0 or (i <= 2 and i_rv >= 4294967087) or (i >= 3 and i_rv >= 4294944443):
	       i_rv = unpack_uint32(os.urandom(4))
	    seed.append(i_rv)
      self.r.set_package_seed((ctypes.c_ulong*6)(*seed))

   def __init__(self, delta_seed=-1, comm=None):
      self.r = rngstream.RngStream()

      self.set_package_seed(delta_seed)
      seed = (ctypes.c_ulong * 6)()
      if comm is None:
	 self.rng = self.r.create_stream()
	 self.l_rng = None
	 self.g_rng = None
      else:
	 if comm.rank == 0:
	    self.g_rng = self.r.create_stream()
	    self.r.get_state(self.g_rng,seed)
	    seed_comm = comm.bcast(list(seed), root=0)
	 else:
	    seed_comm = None
	    seed_comm = comm.bcast(seed_comm, root=0)
	    seed = (ctypes.c_ulong * 6)(*seed_comm)
	    self.r.set_package_seed(seed)
	    self.g_rng = self.r.create_stream()
	 if comm.rank == 0:
	    self.l_rng = self.r.create_stream()
	    seed_new = (ctypes.c_ulong * 6)()
	    for i in range(1,comm.size):
	       t_rng = self.r.create_stream()
	       self.r.get_state(t_rng,seed_new)
	       comm.send(list(seed_new), i, tag=i)
	 else:
	    seed_new_comm = comm.recv(source=0, tag=comm.rank)
	    self.r.set_package_seed((ctypes.c_ulong * 6)(*seed_new_comm))
	    self.l_rng = self.r.create_stream()
	 self.rng = self.l_rng
      if comm is None:
	 self.r.get_state(self.rng,seed)
	 print "Initialized RngStream rng to ", list(seed)
      else:
	 if comm.rank == 0:
	    self.r.get_state(self.g_rng,seed)
	    print comm.rank, "Initialized RngStream g_rng to ", list(seed)
	 self.r.get_state(self.l_rng,seed)
	 print comm.rank, "Initialized RngStream l_rng to ", list(seed)
      self.saved_normal_rv = None

   def switch_to_common(self):
      if self.g_rng is not None:
	 self.rng = self.g_rng
	 self.saved_normal_rv = None
   def switch_to_local(self):
      if self.l_rng is not None:
	 self.rng = self.l_rng
	 self.saved_normal_rv = None

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
   def normal_variate(self, std_dev):
      if self.saved_normal_rv is None:
	 rsq = 0.0
	 while rsq > 1.0 or rsq == 0.0:
	    v1 = 2.0*self.r.float_uniform_01(self.rng)-1.0
	    v2 = 2.0*self.r.float_uniform_01(self.rng)-1.0
	    rsq = v1*v1+v2*v2
	 fac = math.sqrt(-2.0*math.log(rsq)/rsq)
	 self.saved_normal_rv = v1*fac
	 return v2*fac
      else:
	 rv = self.saved_normal_rv
	 self.saved_normal_rv = None
	 return rv
   def normal(self, std_dev, size=None):
      if size is None:
	 return self.normal_variate(std_dev)
      else:
	 out = np.zeros(size)
	 for x in np.nditer(out, op_flags=['readwrite']):
	    x[...] = self.normal_variate(std_dev)
	 return out
   def shuffle_in_place(self, list):
      random.shuffle(list, lambda : self.r.float_uniform_01(self.rng) )
