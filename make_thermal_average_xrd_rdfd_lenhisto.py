import misc_calc_lib
import quippy
import numpy as np
import argparse


# This is a script to calculate thermally average powder spectra and rdfs as well as 
# thermally weighted histograms of cell vector length.
# Please consult the README for further explanations.


# IMPORTANT: THIS SCRIPT ONLY WORKS FOR ONE VERSION OF THIS SCRIPT RUNNING IN A FOLDER AT A GIVEN TIME

# Referencces:

# original formula for weights/remaing phase space volume gamma:
# S. Martiniani, J. D. Stevenson, D. J. Wales, D. Frenkel,
# Superposition enhanced nested sampling, Physical Review X 4 (3) (2014) 031034.
# (I approximated the recorded points to have the same distance in the logarithm of the phase space volume (for one iteration))

# principle of integration of the partition function and the average quantities:
# L. B. Partay, A. P. Bartok, G. Csanyi, 
# Efficient sampling of atomic configurational spaces, The Journal of Physical Chemistry B 114 (32) (2010) 10502-10512.


# Path to the QUIP build. Needs to be adjusted on each system.
# (In my case it's "/home/lsc23/QUIP_git_with_GAP/build/linux_x86_64_gfortran_openmp")
QUIP_path = "$QUIP_path"

if QUIP_path == "$QUIP_path":
   print("Error! QUIP_path needs to be set in make_thermal_average_xrd_rdfd_lenhisto.py. Aborting.")
   quit()



k_B = 8.6173303*10.0**(-5) # [eV/K] https://physics.nist.gov/ (accessed 2017/10/04 16:50)

do_rdfd = False # RDFs in QUIP are not using periodic cells. This makes it very hard to compare different cells of the same structure. Hence, it is turned off!
               # If set to "True" the script uses a 6x6x6 supercell for the comparison structures.


parser = argparse.ArgumentParser()

parser.add_argument("-fn", "--filename", help="Name of '.extxyz'/'.xyz' file to analyse")
parser.add_argument("-Ts", "--T_array", help='array of T in format "T_1 T_2 ... T_N". Converts to integers at the moment.')
parser.add_argument("-nc", "--n_cull", help="n_cull of nested sampling run")
parser.add_argument("-nw", "--n_walkers", help="n_walkers of nested sampling run")
parser.add_argument("-s", "--comp_strucs", help="Structures for xrd spectrum identification in format \"struc1 struc2 struc3\"")

args = parser.parse_args()

filename = args.filename
T_range = []  # Temperatures to be weighted at
for el in args.T_array.split():
   T_range.append(int(el))
# These are the comparison structures whose rdfds (if on) and xrds get automatically calculated. They must be appropriately defined in create_at_accord_struc (see misc_calc_lib.py). 
# E.g. ["hcp_Hennig_MEAM", "omega_Hennig_MEAM", "bcc", "fcc"]
comparison_structures = []
for el in args.comp_strucs.split():
   comparison_structures.append(el)

n_cull = int(args.n_cull)
n_walker = int(args.n_walkers)

if filename[len(filename) - len(".extxyz") ::].find(".extxyz") == 0:
   extens = ".extxyz"
elif filename[len(filename) - len(".xyz") ::].find(".xyz") == 0:
   extens = ".xyz"
else:
   print("ERROR: Filename neither '.extxyz' nor '.xyz' file! Aborting!")
   quit()


if filename.find("/") < 0:
   reduced_filename = filename
else:
   reduced_filename = filename[-filename[::-1].find("/"):]
raw_reduced_filename = reduced_filename[:len(reduced_filename)-len(extens)]




print(reduced_filename)

#quit()

iter_nr = []
enthalpy = []
gamma_log = []
box_volume = []

# rdf range parameters:

a_0 = 0.0
a_end = 10.0
n_a = 100
r_range = [a_0, a_end]
# xrd parameters

two_theta_range = '"0.0 180.0"'
n_two_theta = 361

do_xrd = True

#This defines the percentage (according to probabilities of each structure) which we define siginficant enough to calculate xrds on. We only calculate for 'siginficant_part' most likely structures.
significant_part = 0.95#1 - 10.0**(-16)

threshold = (1 - significant_part)/2.0


inputs = quippy.AtomsReader(filename)


for at in inputs:
   iter_nr.append(at.info["iter"])
   enthalpy.append(at.info["ns_energy"])
   box_volume.append(at.get_volume())
    

print("len(iter_nr) = " + str(len(iter_nr)) + " len(enthalpy) = " + str(len(enthalpy)))



volume_one_iter = 1.0
for i in range(0,n_cull):
   volume_one_iter = volume_one_iter * (n_walker - i) / (n_walker + 1 - i)

volume_one_iter_log = np.log(volume_one_iter)
#print "volume_one_iter = " + str(volume_one_iter)

# we assume that we start counting at iteration 0
nr_last_iter_change = 0

last_iters = [0,iter_nr[0]]

# Creates gamma_log for weigth calculation including if the trajectory has changing number of iterations between samples

for i in range(0,len(iter_nr)):
   # detects if iteration number changes
   if iter_nr[i] > iter_nr[i-1]:
      for j in range(nr_last_iter_change,i):
         if j==0:
            gamma_log.append( ( last_iters[1] - last_iters[0] +  1.0 / (i - nr_last_iter_change) ) * volume_one_iter_log )
         else:
            if j==nr_last_iter_change:
               gamma_log.append(gamma_log[j-1] + ( last_iters[1] - last_iters[0] -1 +  1.0 / (i - nr_last_iter_change) ) * volume_one_iter_log)
            else:
               gamma_log.append(gamma_log[j-1] + ( 1.0 / (i - nr_last_iter_change) ) * volume_one_iter_log)
      nr_last_iter_change = i
      last_iters[0] = last_iters[1]
      last_iters[1] = iter_nr[i]

   # special case for last recoreded iteration
   if (i == len(iter_nr) -1):
      for j in range(nr_last_iter_change,len(iter_nr)):
         if j==nr_last_iter_change:
            gamma_log.append(gamma_log[j-1] + ( last_iters[1] - last_iters[0] -1 +  1.0 / (len(iter_nr) - nr_last_iter_change) ) * volume_one_iter_log)
         else:
            gamma_log.append(gamma_log[j-1] + ( 1.0 / (len(iter_nr) - nr_last_iter_change) ) * volume_one_iter_log)
      nr_last_iter_change = i
      last_iters[0] = last_iters[1]
      last_iters[1] = iter_nr[i]


xrd_matrix = []
rdf_matrix = []


at = inputs[len(inputs)-1]
z = at.get_atomic_numbers()[0]


rdfd_results = misc_calc_lib.rdfd_QUIP(QUIP_path,at,n_a,r_range)


xrd_results = misc_calc_lib.xrd_QUIP(QUIP_path,at,n_two_theta,two_theta_range)
angle = xrd_results[0]
xrd_temp = xrd_results[1]


# Making the thermal average
for T in T_range:
   beta = 1/(k_B * T)

   weight = []
   gamma_log_1_beta = []
   gamma_log_2_beta = []

   for i in range(0,len(iter_nr)-1):
      gamma_log_1_beta.append(gamma_log[i] - beta * enthalpy[i])
      gamma_log_2_beta.append(gamma_log[i+1] - beta * enthalpy[i])

   shift = max(max(gamma_log_1_beta),max(gamma_log_2_beta))

   for i in range(0,len(iter_nr)-1):
      #print str(gamma[i+1])
      weight.append( np.exp(gamma_log[i] - beta * enthalpy[i] - shift) - np.exp(gamma_log[i+1] - beta * enthalpy[i] - shift) )
      #print weight

   partion_fct = sum(weight)


   rdf_null = 0.0 * rdfd_results[1]
   xrd_null = 0.0 * xrd_temp
   V_array = []
   a_lat_array = []
   b_lat_array = []
   c_lat_array = []

   xrd_matrix = []
   rdf_matrix = []

   cumulative_weights = 0.0

   part_fct_red = 0.0

   for i_at,at in enumerate(inputs):
      V_array.append(at.get_volume())
      a_b_c = at.get_cell_lengths_and_angles()[0:3]
      a_b_c.sort()
      a_lat_array.append(a_b_c[0])
      b_lat_array.append(a_b_c[1])
      c_lat_array.append(a_b_c[2])

      if i_at < len(weight):
         cumulative_weights += weight[i_at]/partion_fct


      if ((i_at < len(weight)) and (cumulative_weights >= threshold)) and (cumulative_weights <= 1 - threshold): #at.get_volume() < at.get_number_of_atoms()*50:

#         print("cumu_weight = " + str(cumulative_weights))

         if do_xrd == True:
            rdfd_results = misc_calc_lib.rdfd_QUIP(QUIP_path,at,n_a, r_range)

            xrd_results = misc_calc_lib.xrd_QUIP(QUIP_path,at,n_two_theta,two_theta_range)
#            angle = xrd_results[0]
            xrd_temp = xrd_results[1]
#      misc_calc_lib.submit_commands(["rm temp_atom.xrd.raw temp_atom.xrd temp_atom.xyz"])

            xrd_matrix.append(xrd_temp)
            rdf_matrix.append(rdfd_results[1])
            part_fct_red = part_fct_red + weight[i_at] 
            r = rdfd_results[0]

#            print("xrd_temp is:")
#            print(xrd_temp)

#            print(xrd_matrix)
#            quit()
         else:
            xrd_matrix.append(xrd_null)
            rdf_matrix.append(rdf_null)
      else:
         xrd_matrix.append(xrd_null)
         rdf_matrix.append(rdf_null)

#   print(xrd_matrix)


   rdf = rdf_matrix[0]*0.0
   xrd = xrd_matrix[0]*0.0
   V = V_array[0]*0.0

   print("len(weight) = " + str(len(weight)) + " len(rdf_matrix) = " + str(len(rdf_matrix)) + " len(xrd_matrix) = " + str(len(xrd_matrix)))
   cumulative_weights = 0
   for i in range(0,len(xrd_matrix) - 1):
#      cumulative_weights += weight[i]/partion_fct

#      if (cumulative_weights >= threshold) and (cumulative_weights <= 1 - threshold):
      print("iteration for calculating averages = " + str(i))
      print("len(weight) = " + str(len(weight)) + " len(rdf_matrix) = " + str(len(rdf_matrix)) + " len(xrd_matrix) = " + str(len(xrd_matrix)))
      rdf = rdf + rdf_matrix[i]*weight[i]
      xrd = xrd + xrd_matrix[i]*weight[i]
      V = V + V_array[i]*weight[i]

   rdf = rdf/part_fct_red
   xrd = xrd/part_fct_red
   V = V/partion_fct


   a_histo, bin_limits = np.histogram(a_lat_array[:-1],bins=n_a,range=(a_0,a_end),weights=weight)
   b_histo, bin_limits = np.histogram(b_lat_array[:-1],bins=n_a,range=(a_0,a_end),weights=weight)
   c_histo, bin_limits = np.histogram(c_lat_array[:-1],bins=n_a,range=(a_0,a_end),weights=weight)

   lat_vec_mean_histo = (a_histo + b_histo + c_histo)/3.0

   print(bin_limits)

   if do_rdfd:
#      print("we got this far!")
   #   quit()
      with open(raw_reduced_filename + "_signifpart_" + str(significant_part) + ".custom_T_" + str(T) +"_rdfd","w") as rdf_f:
         for i in range(0,len(r)):
            rdf_f.write(str(r[i]) + "   " + str(rdf[i]) + "\n")

   if do_xrd:
      with open(raw_reduced_filename + "_signifpart_" + str(significant_part) + ".custom_T_" + str(T) + "_xrd", "w") as xrd_f:
         for i in range(0,len(angle)):
            xrd_f.write(str(angle[i]) + "   " +str(xrd[i]) + "\n") 

   print "partion_fct at " + str(T) + " K = " + str(partion_fct)


   with open(raw_reduced_filename + "_signifpart_" + str(significant_part) + ".custom_T_" + str(T) + "_lattice_len_histo", "w") as histo_f:
      histo_f.write("#   len lat vec [Angstrom]   a      b     c      mean(a,b,c)\n")
      for i in range(0,len(a_histo)):
         d_a = (a_end - a_0)/n_a
         bin_middle = a_0 + 0.5 * d_a + i * d_a
         histo_f.write(str(bin_middle) + "   " + str(a_histo[i])  + "   " + str(b_histo[i]) + "   " + str(c_histo[i]) + "   " + str(lat_vec_mean_histo[i]) + "\n")


   z = at.get_atomic_numbers()[0]
   V_aver_per_at = V/len(at)

   for struc in comparison_structures:

      at_average = misc_calc_lib.create_at_accord_struc(V_aver_per_at ,z , struc)
      comparison_struc_name_raw = struc + "_V_mean_of_" + raw_reduced_filename + ".T_" + str(T)
      comparison_struc_name = comparison_struc_name_raw + ".xyz"

      comparison_struc_name_rdfd = comparison_struc_name_raw + "_rdfd"
      comparison_struc_name_xrd = comparison_struc_name_raw + "_xrd"

      at_average.write(comparison_struc_name)
      rdfd_results = misc_calc_lib.rdfd_QUIP(QUIP_path,quippy.supercell(at_average,6,6,6),n_a,r_range)  # I'm using a supercell to get at least the positions right. 
      xrd_results = misc_calc_lib.xrd_QUIP(QUIP_path,at_average,n_two_theta,two_theta_range)
      if do_rdfd:
         with open(comparison_struc_name_rdfd, "w") as rdf_f:
            for i,el in enumerate(rdfd_results[0]):
               rdf_f.write(str(rdfd_results[0][i]) + "   " + str(rdfd_results[1][i]) + "\n")

      if do_xrd:
         with open(comparison_struc_name_xrd, "w") as xrd_f:
            for i,el in enumerate(xrd_results[0]):
               xrd_f.write(str(xrd_results[0][i]) + "   " + str(xrd_results[1][i]) + "\n")
