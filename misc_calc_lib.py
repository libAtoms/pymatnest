from __future__ import print_function

#These creations are in part different from the ones used for the MEAM quasiharmonic calculations.
import subprocess
import numpy as np
import matplotlib.pyplot as plt

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
   allowed_struc = ["bcc", "fcc", "hcp", "hcp_Hennig_MEAM", "omega_Hennig_MEAM","hcp_Hennig_DFT", "omega_Hennig_DFT","bcc_Trinkle_TB", "hcp_Trinkle_TB", "omega_Trinkle_TB","bcc_Trinkle_fitting", "hcp_Trinkle_fitting", "omega_Trinkle_fitting"]

   if struc in allowed_struc:
      print("\nAccepted structure name!\n")


#   print("IMPORTANT: WARNING! c_vs_a set to 1. Only use bcc!")

   if struc == "bcc":
      at = create_beta_V(V, z)
   elif struc == "fcc":
      at = create_fcc_V(V,z)
   elif struc == "hcp":
      c_vs_a = (8.0/3.0)**(0.5)
      at = create_hcp_custom(c_vs_a, V, z)
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
def write_cell_and_param(at, name_raw, template_name_cell):

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
               for i in xrange(0,len(at)):
                  write_lines.write(str(at.get_atomic_numbers()[i]) + " " + " ".join(map(str, scale_pos[i,:])) + "\n")

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
      subprocess.check_output(command[i],shell=True)


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


def raw_filename_xyz_extxyz(filepath):


   if filepath[len(filepath) - len(".extxyz") ::].find(".extxyz") == 0:
      extens = ".extxyz"
   elif filepath[len(filepath) - len(".xyz") ::].find(".xyz") == 0:
      extens = ".xyz"
   else:
      print("ERROR: Filename neither '.extxyz' nor '.xyz' file! Aborting!")
      quit()
   
   
   if filepath.find("/") < 0:
      reduced_filename = filepath
   else:
      reduced_filename = filepath[-filepath[::-1].find("/"):]
   raw_reduced_filename = reduced_filename[:len(reduced_filename)-len(extens)]

   return raw_reduced_filename

# Returns the x value (usually temperature) of the maximum y of a two column file over x_range = [x_start, x_end].
# Comments in the file are signfied by "#".
def find_x_of_max_y_of_file(filename, x_range, x_col_index, y_col_index):

   vals_of_interests = []
   with open(filename, "r") as flines:
      for i, line in enumerate(flines):
         if line[0] != "#" and i > 0: #Ommitting first lines and comments
            line_float_split = map(float, line.split())
            # If temperature is in the temeprature range of interest, append
            if line_float_split[0] >= x_range[0] and line_float_split[0] <= x_range[1]:
               vals_of_interests.append(line_float_split)

   vals_trans = []
   for el in zip(*vals_of_interests):
      vals_trans.append(list(el))

   max_index = find_max_index(vals_trans[y_col_index])

   if len(max_index) > 1:
      print("Error! There should only be one maximum. Either you got more two data points which are the same values and the maximum, which is very unlikely, or there's something wrong with the input file. Aborting!")
      quit()

   return vals_trans[x_col_index][max_index[0]]


#start for recursive function to find maximum values
def find_max_index(array):

   return find_max_index_recursive(array, [])



# Yields a list of the index of the maximum values (In case there are more than one).
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
# x_col_index and y_col_index are the indices of the columns of interest. In our case they should be T and C_P which are 0 and 3.
def find_avg_std_max_pos(filename_beginning, filename_ending, run_nr_start, run_nr_end, x_range, x_col_index, y_col_index):

   pos_list = []
   for i in xrange(run_nr_start, run_nr_end + 1, 1):

      filename = filename_beginning + str(i) + filename_ending
      pos_list.append(find_x_of_max_y_of_file(filename, x_range, x_col_index, y_col_index))

   mean = np.mean(pos_list)
   std = np.std(pos_list, ddof=1)
   mean_error = std/np.sqrt(run_nr_end - run_nr_start + 1)




   return [mean, std, mean_error]


# gets data from a number of files and returns the maximum y_value:
# returns an array with 0th entry being the list of data from the files and the 1st entry being
# the y_max value over all the data
# y_index denotes the index of the column we use as y in an x-y plot in case there are more than two.
def get_data_for_plot(result_name_list, y_index):

   print(result_name_list)

   max_vals = []
   result_array = []
   for filename in result_name_list:
   
      data = np.loadtxt(filename, skiprows = 1, unpack = True)
      max_vals.append(max(data[y_index]))
      result_array.append(data)
   

   max_over_runs = max(max_vals)

   return [result_array, max_over_runs]


def make_name_list(name_raw, prefix, suffix, start_no, end_no):

   name_list = []

   for i in xrange(start_no, end_no + 1):

      name_list.append(prefix + name_raw + str(i) + suffix)

   return name_list

# plots C_P curves. If ask_for_check == True, it asks the user whether the C_P curve is acceptable
def plot_C_P(result_name_list, labels, x_index, y_index, T_mean, T_aim, ask_for_check_bool):


   get_data_results = get_data_for_plot(result_name_list, y_index)

   result_array = get_data_results[0]

   max_over_runs = get_data_results[1]

   if labels == []:
      labels = result_name_list

   for result in result_array:
      plt.plot(result[x_index], result[y_index], '-')
   plt.plot((T_mean, T_mean), (0, 1.2 * max_over_runs), '-')
   plt.plot((T_aim, T_aim), (0, 1.2 * max_over_runs), '--')
   plt.legend(labels)
   plt.show()

   if ask_for_check_bool:
      ask_for_check()


# plots the xrds of result_name_list and comp_name_list (the latter adjusted to fit on screen). It also asks the user whether the plots are okay and quits if the plots are deemed bad.
def plot_xrd(result_name_list, comp_name_list, result_labels, comp_labels, x_index, y_index, n_per_average_batch, ask_for_check_bool, shift_bool):

   if comp_labels == []:
      comp_labels = comp_name_list
   if result_labels == []:
      result_labels = result_name_list

   labels = result_labels + comp_labels
      

# Data for results
   get_data_results = get_data_for_plot(result_name_list, y_index)


# Data for comparison
   get_data_comp = get_data_for_plot(comp_name_list, y_index)

   result_array = get_data_results[0]

   max_over_runs = get_data_results[1]

   if shift_bool == True:
      shift = max_over_runs * 0.1
   elif shift_bool == False:
      shift = 0.0
   else:
      print("Error! shift_bool must be boolean! Aborting")
      quit()

#  This averages batches of runs defined by n_per_average_batch
   intens_batch_np_array_list = []
   angles_array = []
   for i in range(0, len(result_array)/n_per_average_batch):
      batch = []
      for j in range(0, n_per_average_batch):
         batch.append(np.array(result_array[i*n_per_average_batch + j][y_index]))
      intens_batch_np_array_list.append(np.array(batch))

      angles_array.append(result_array[i * n_per_average_batch][x_index])


   intens_mod = [np.ndarray.mean(batch, axis = 0) + shift * (i + len(comp_name_list)) for i, batch in enumerate(intens_batch_np_array_list)] 

   comp_array = get_data_comp[0]


   for i, intens in enumerate(intens_mod):
      plt.plot(angles_array[i], intens, '-')
   
#  Adding the plots of the comparison structures, each scaled so that its respective maximum is the same as the overall maximum of the runs. This enables better comparison as we care more about the peak position than absolute height.
   for i, comp in enumerate(comp_array):
      plt.plot(comp[x_index], comp[y_index]/max(comp[y_index]) * max_over_runs + shift * i, '--')

   plt.legend(labels)

   plt.show()

   if ask_for_check_bool:
      ask_for_check()


# Asks user whether plot is accetable. If user replies with "no", it quits the program.
def ask_for_check():

   got_answer = False

   while got_answer == False:

      go_on = raw_input("Are these results acceptable? Reply with 'yes' or 'no':")

      if go_on == "yes":
         go_on_bool = True
         got_answer = True
      elif go_on == "no":
         go_on_bool = False
         print("ALOGRITHM ERROR: The results have been deemed unsatisfactory by the user. Is something not converged? Did we get the wrong phase? Aborting!")
         quit()
      else:
         print("Error! Expecting 'yes' or 'no' but got '" + go_on + "'. Asking again!")


#returns a list of postions of occurance of a string (word) in another string (orig_string)
def return_pos_in_string(orig_string, word):

   string_old = ""
   string = orig_string
   word_pos_list = []
#  make a list of the positions of the word
   while string.find(word) >= 0:
      word_pos_list.append(string.find(word) + len(string_old))
      string_old = string_old + string[:string.find(word) + len(word)]
      string = string[string.find(word) + len(word):]

   return word_pos_list

#returns the last part of a path
def last_part_path(path):

   return path[return_pos_in_string(path, "/")[-1] + 1:]
