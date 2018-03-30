#These creations are in part different from the ones used for the MEAM quasiharmonic calculations.
import subprocess
import numpy as np

def create_hcp_custom(c_vs_a, V_per_at, z):
   import quippy
   import numpy as np

   V = 2.0 * V_per_at
   a = (2 * V / (3.0**(0.5) * c_vs_a))**(1.0/3.0)
   c = c_vs_a * a

   lattice = []
   lattice.append([3.0**(0.5) /2.0 * a,-a/2.0,0])
   lattice.append([3.0**(0.5) /2.0 * a, a/2.0,0])
   lattice.append([0,0,c])
   lattice = np.transpose(lattice)
   unitcell = quippy.Atoms(n=0, lattice=lattice)

   pos = []
   pos.append([3.0**(0.5) /6.0 * a,0,0.0])
   pos.append([3.0**(0.5) /2.0 * a,0,c/2.0])

   for i in range(0,len(pos)):
      unitcell.add_atoms(pos[i],z)

   return unitcell


def create_omega_custom(c_vs_a, V_per_at, z):
   import quippy

   V = 3.0 * V_per_at
   a = (2 * V / (3.0**(0.5) * c_vs_a))**(1.0/3.0)
   c = c_vs_a * a

   lattice = []
   lattice.append([3.0**(0.5) /2.0 * a,-a/2.0,0])
   lattice.append([3.0**(0.5) /2.0 * a, a/2.0,0])
   lattice.append([0,0,c])
   lattice = np.transpose(lattice)
   unitcell = quippy.Atoms(n=0, lattice=lattice)

   pos = []
   pos.append([3.0**(0.5) /6.0 * a,0,0.0])
   pos.append([3.0**(0.5) /2.0 * a,0,c/2.0])
   pos.append([3.0**(0.5) * 5.0/6.0 * a,0,c/2.0])

   for i in range(0,len(pos)):
      unitcell.add_atoms(pos[i],z)

   return unitcell


def create_beta_V(V_per_at, z):
   import quippy

   a = (2.0 * V_per_at)**(1.0/3.0)

   unitcell = quippy.structures.bcc1(a, z)

   return unitcell

def create_fcc_V(V_per_at, z):
   import quippy

   a = (4.0 * V_per_at)**(1.0/3.0)

   unitcell = quippy.fcc(a, z)

   return unitcell

def create_at_accord_struc(V ,z , struc):
#  c/a ratios from (see structure names): Hennig et al., PHYSICAL REVIEW B 78, 054121 (2008)
#                                         Trinkle et al., PHYSICAL REVIEW B 73, 094123 (2006)
   allowed_struc = ["bcc", "fcc", "hcp_Hennig_MEAM", "omega_Hennig_MEAM","hcp_Hennig_DFT", "omega_Hennig_DFT","bcc_Trinkle_TB", "hcp_Trinkle_TB", "omega_Trinkle_TB","bcc_Trinkle_fitting", "hcp_Trinkle_fitting", "omega_Trinkle_fitting"]

   if struc in allowed_struc:
      print("\nAccepted structure name!\n")


#   print("IMPORTANT: WARNING! c_vs_a set to 1. Only use bcc!")

   if struc == "bcc":
      at = create_beta_V(V, z)
   elif struc == "fcc":
      at = create_fcc_V(V,z)
   elif struc == "hcp_Hennig_MEAM":
      c_vs_a = 1.596
      at = create_hcp_custom(c_vs_a, V, z)
   elif struc == "omega_Hennig_MEAM":
      c_vs_a = 0.611
      at = create_omega_custom(c_vs_a, V, z)
   elif struc == "hcp_Hennig_DFT":
      c_vs_a = 1.583
      at = create_hcp_custom(c_vs_a, V, z)
   elif struc == "omega_Hennig_DFT":
      c_vs_a = 0.619
      at = create_omega_custom(c_vs_a, V, z)
   elif struc == "bcc_Trinkle_TB":
      at = create_beta_custom(V, z)
   elif struc == "hcp_Trinkle_TB":
      c_vs_a = 4.71/2.94
      at = create_hcp_custom(c_vs_a, V, z)
   elif struc == "omega_Trinkle_TB":
      c_vs_a = 2.84/4.58
      at = create_omega_custom(c_vs_a, V, z)
   elif struc == "bcc_Trinkle_fitting":
      at = create_beta_custom(V, z)
   elif struc == "hcp_Trinkle_fitting":
      c_vs_a = 1.588
      at = create_hcp_custom(c_vs_a, V, z)
   elif struc == "omega_Trinkle_fitting":
      c_vs_a = 0.613
      at = create_omega_custom(c_vs_a, V, z)
   else:
      print("\nERROR: Structure name '" + struc + "' not known!\n\nAllowed names are:")

      for dummy_struc in allowed_struc:
         print(dummy_struc)
      print("")

      quit()

   return at

# Writes .cell as well as .param files based on templates
def write_files(at, name_raw, template_name_cell):

   write_file_cell = name_raw + ".cell"

   write_cell(at, write_file_cell, template_name_cell)
   write_param(write_file_cell, template_name_cell)

# Write a .cell file based on a template
def write_cell(at, write_file_cell, template_name_cell):

   with open(write_file_cell, "w") as write_lines:
      with open(template_name_cell, "r") as template_lines:
         for line in template_lines:
            write_lines.write(line)
            if line.find("%block lattice_cart") >= 0:
               cell = at.get_cell()
               for i in xrange(0,3):
                  write_lines.write("   ".join(map(str, cell[i,:])) + "\n")

            if line.find("%block positions_frac") >= 0:
               scale_pos = at.get_scaled_positions()
               print("Writing positions!")
               for i in xrange(0,at.get_number_of_atoms()):
                  write_lines.write("Ti " + " ".join(map(str, scale_pos[i,:])) + "\n")

def write_param(write_file_cell, template_name_cell):

   with open(write_file_cell[0:-len(".cell")] + ".param", "w") as write_lines:
      with open(template_name_cell[0:-len(".cell")] + ".param", "r") as template_lines:
         for line in template_lines:
#            print("test")
            write_lines.write(line)


def write_sub_file(submis_com, nr_node_tot, sub_file):

   sub_template = "CASTEP_TEMPLATE.sub"

   with open(sub_file, "w") as write_lines:
      with open(sub_template, "r") as template_lines:
         for line in template_lines:
            write_lines.write(line)
      write_lines.write(submis_com + "\n")
      write_lines.write("") 

   command = []
   command.append("sed -i 's/\$nr_node_tot/" + str(nr_node_tot) + "/g' " + sub_file)
   submit_commands(command)
   command = []


def submit_commands(command):

   for i in range(0,len(command)):
      print(command[i])
      subprocess.call(command[i],shell=True)


def calc_aspec_ratio(at):

   cell = at.get_cell()

   V = at.get_volume()

   aspec_ratio_array = []

   for i in range(0,3):

      # cell vector not in plane parallel to the cell surfaces whose distance are to measure
      vec_0 = cell[i,:]

      # cell vectors defining the plane
      vec_1 = cell[(i+1)%3,:]
      vec_2 = cell[(i+2)%3,:]

      # vector orthogonal to plane 
      cross = np.cross(vec_1,vec_2)
      # normalized orthogonal vector
      cross_norm =  cross/np.sqrt(np.dot(cross,cross))
      # distance between cell surfaces normalized by cell volume:
      aspec_ratio_array.append(abs(np.dot(cross_norm,vec_0))/V**(1.0/3))



   return aspec_ratio_array



# Calcualtes powders spectrum via QUIP
def xrd_QUIP(QUIP_path, at, n_two_theta, two_theta_range):

   xyz_fime = "temp_atom.xyz"

   at.write(xyz_fime)

   xrd_raw_fime = xyz_fime[:len(xyz_fime)-len(".xyz")] + ".xrd.raw"
   xrd_final_fime = xyz_fime[:len(xyz_fime)-len(".xyz")] + ".xrd"

   command = []
   command.append( QUIP_path + '/structure_analysis_traj type=xrd xrd_2theta_range=' + two_theta_range + " xrd_n_2theta=" + str(n_two_theta)+ " infile=" + xyz_fime + " outfile=" + xrd_raw_fime)
   command.append( QUIP_path + "/mean_var_correl infile=" + xrd_raw_fime + " outfile=" + xrd_final_fime + " mean")

   submit_commands(command)
   command = []

   angle, xrd_temp = np.loadtxt( "temp_atom.xrd", skiprows=1, unpack = True)

   submit_commands(["rm " + xyz_fime + "*" ])
   submit_commands(["rm " + xyz_fime[:len(xyz_fime)-len(".xyz")] + ".xrd*"])

   return [angle,xrd_temp]



#Calculates radial distribution function via QUIP
def rdfd_QUIP(QUIP_path, at, n_a, r_range):


   rdfd_n_bins = int(n_a)
   rdfd_bin_width = float(r_range[1] - r_range[0])/rdfd_n_bins

   xyz_fime = "temp_atom.xyz"

   at.write(xyz_fime)

   rdfd_raw_fime = xyz_fime[:len(xyz_fime)-len(".xyz")] + ".rdfd.raw"
   rdfd_final_fime = xyz_fime[:len(xyz_fime)-len(".xyz")] + ".rdfd"

   command = []
   command.append( QUIP_path + '/structure_analysis_traj type=rdfd rdfd_n_bins=' + str(rdfd_n_bins) + " rdfd_bin_width=" + str(rdfd_bin_width)+ " infile=" + xyz_fime + " outfile=" + rdfd_raw_fime)
   command.append( QUIP_path + "/mean_var_correl infile=" + rdfd_raw_fime + " outfile=" + rdfd_final_fime + " mean")

   submit_commands(command)
   command = []

   angle, rdfd_temp = np.loadtxt( "temp_atom.rdfd", skiprows=1, unpack = True)

   submit_commands(["rm " + xyz_fime + "*" ])
   submit_commands(["rm " + xyz_fime[:len(xyz_fime)-len(".xyz")] + ".rdfd*"])

   return [angle,rdfd_temp]



"""
# Returns the x value (usually temperature) of the maximum y of a two column file over x_range = [x_start, x_end].
# Comments in the file are signfied by "#".
def find_x_of_max_y_of_file(filename, x_range):

   vals_of_interests = []
   with open(filename, "r") as flines:
      for line in flines:
         if line[0] != "#":
            line_float_split = map(float, line.split())
            # If temperature is in the temeprature range of interest, append
            if line_float_split[0] >= x_range[0] and line_float_split[0] <= x_range[1]:
               vals_of_interests.append(line_float_split)

   vals_trans = zip(*vals_of_interests)
   print(vals_trans)

   max_index = find_max_index(vals_trans[1])

   if len(max_index) > 1:
      print("Error! There should only be one maximum. Either you got more two data points which are the same values and the maximum, which is very unlikely, or there's something wrong with the input file. Aborting!")
      quit()

   return


#start for recursive function to find maximum values
def find_max_index(array):

   return find_max_index_recursive(array, [])


# Yields a list of the index of the maximum values.
def find_max_index_recursive(array, index_list):

   max_array = max(array)
   max_index = array.index(max_array)
   index_list.append(max_index)

   # array with maximum replaced(hence "repl_" array):
   repl_array = array[0:max_index] + [max_array - 999] + array[max_index + 1:]
   if max_array == max(repl_array):

      index_list = find_max_index_recursive(repl_array, index_list)

   result_list = list(index_list)
   return result_list


# Automatic calculation of average and standard deviation of maximum position (presumably of the C_P curve) for a number of results of different runs
# This assumes that the results are however already calculated.
def find_avg_std_max_pos(filename_raw, run_nr_start, run_nr_end, x_range):

   pos_list = []
   for i in xrange(run_nr_start, run_nr_end, 1):

      filename = filename_raw + str(i)
      pos_list.append(find_x_of_max_y_of_file(filename, x_range))

   mean = np.mean(pos_list)
   std = np.std(pos_list, ddof=1)




   return [mean, std]
"""
