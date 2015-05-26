set term pdf

set key outside

set output "cluster.Cv.pdf"
set title "LJ6 cluster, sigma=3"
set xlabel "T / epsilon"
set ylabel "Cv per atom"
set autoscale y
plot "analysis.cluster.MC.Livia" u ($1/0.1):($4/6) w l lt 1 color black ti "Livia reference", \
     "analysis.run_test.cluster.MC" u ($1/0.25):($4/6) w l lt 1 color red ti "MC", \
     "analysis.run_test.cluster.MD" u ($1/0.25):($4/6) w l lt 1 color cyan ti "MD"
unset log x
unset log y
set autoscale x
set autoscale y

set output "cluster.energies.pdf"
set xlabel "n iter / n walkers"
set ylabel "E / epsilon"
set yrange [:1]
plot "< tail -n +2 run_test.cluster.MC/*energies" u ($1/1024):($2/6) w l lt 1 color red ti "MC", \
     "< tail -n +2 run_test.cluster.MD/*energies" u ($1/1024):($2/6) w l lt 1 color cyan ti "MD"
unset log x
unset log y
set autoscale x
set autoscale y

####################################################################################################

set output "periodic.Cp.pdf"
set title "LJ64 periodic, sigma=3"
set xlabel "T / epsilon"
set ylabel "Cp per atom"
set log y
set yrange [1:]
set xrange [0:1.5]
plot "analysis.periodic.MC.Livia" u ($1/0.1):($4/64) w l lt 1 color black ti "Livia reference", \
     "analysis.run_test.periodic.MC" u ($1/0.25):($4/64) w l lt 1 color red ti "MC", \
     "analysis.run_test.periodic.MC_velo" u ($1/0.25):($4/64) w l lt 1 color magenta ti "MC velo", \
     "analysis.run_test.periodic.MC_velo_microcanonical" u ($1/0.25):($4/64) w l lt 1 color grey ti "MC velo mc", \
     "analysis.run_test.periodic.MD_orig" u ($1/0.25):($4/64) w l lt 1 lw .5 color blue ti "MD orig", \
     "analysis.run_test.periodic.MD_velo_rando" u ($1/0.25):($4/64) w l lt 1 lw .5 color cyan ti "MD velo rando", \
     "analysis.run_test.periodic.MD_velo_random_walk" u ($1/0.25):($4/64) w l lt 2 lw .5 color cyan ti "MD velo random walk", \
     "analysis.run_test.periodic.MD_velo_rando_long" u ($1/0.25):($4/64) w l lt 1 lw .5 color orange ti "MD velo rando long", \
     "analysis.run_test.periodic.MD_velo_rando_Efuzz" u ($1/0.25):($4/64) w l lt 1 lw .5 color yellow ti "MD velo rando Efuzz", \
     "analysis.run_test.periodic.MD_velo_rev" u ($1/0.25):($4/64) w l lt 1 lw .5 color green ti "MD velo rev"
unset log x
unset log y
set autoscale x
set autoscale y

set output "periodic.energies.pdf"
set xlabel "n iter / n walkers"
set ylabel "H / epsilon"
set yrange [:15]
plot "< cat run_test.periodic.MC/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 color red ti "MC", \
     "< cat run_test.periodic.MC_velo/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 color magenta ti "MC velo", \
     "< cat run_test.periodic.MC_velo_microcanonical/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 color grey ti "MC velo mc", \
     "< cat run_test.periodic.MD_orig/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color blue ti "MD orig", \
     "< cat run_test.periodic.MD_velo_rando/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color cyan ti "MD velo rando", \
     "< cat run_test.periodic.MD_velo_random_walk/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 2 lw .5 color cyan ti "MD velo random walk", \
     "< cat run_test.periodic.MD_velo_rando_long/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color orange ti "MD velo rando long", \
     "< cat run_test.periodic.MD_velo_rando_Efuzz/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color yellow ti "MD velo rando Efuzz", \
     "< cat run_test.periodic.MD_velo_rev/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color green ti "MD velo rev"
unset log x
unset log y
set autoscale x
set autoscale y

set output "periodic.volumes.pdf"
set xlabel "n iter / n walkers"
set ylabel "V / sigma cubed"
set log y
# set yrange [:15]
plot "proc_traj.run_test.periodic.MC" every 100 u ($1/128):($3/64) w l lt 1 color red ti "MC", \
     "proc_traj.run_test.periodic.MC_velo" every 100 u ($1/128):($3/64) w l lt 1 color magenta ti "MC velo", \
     "proc_traj.run_test.periodic.MC_velo_microcanonical" every 100 u ($1/128):($3/64) w l lt 1 color grey ti "MC velo mc", \
     "proc_traj.run_test.periodic.MD_orig" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color blue ti "MD orig", \
     "proc_traj.run_test.periodic.MD_velo_rando" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color cyan ti "MD velo rando", \
     "proc_traj.run_test.periodic.MD_velo_random_walk" every 100 u ($1/128):($3/64) w l lt 2 lw .5 color cyan ti "MD velo random walk", \
     "proc_traj.run_test.periodic.MD_velo_rando_long" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color orange ti "MD velo rando long", \
     "proc_traj.run_test.periodic.MD_velo_rando_Efuzz" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color yellow ti "MD velo rando Efuzz", \
     "proc_traj.run_test.periodic.MD_velo_rev" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color green ti "MD velo rev"
unset log x
unset log y
set autoscale x
set autoscale y




# ####################################################################################################
# 
# set output "periodic.Cp.pdf"
# set title "LJ64 periodic, sigma=3"
# set xlabel "T / epsilon"
# set ylabel "Cp per atom"
# set log y
# set yrange [1:]
# set xrange [0:1.5]
# plot "analysis.periodic.MC.Livia" u ($1/0.1):($4/64) w l lt 1 color black ti "Livia reference", \
#      "analysis.run_test.periodic.MC" u ($1/0.25):($4/64) w l lt 1 color red ti "MC", \
#      "analysis.run_test.periodic.MC_velo" u ($1/0.25):($4/64) w l lt 1 color magenta ti "MC velo", \
#      "analysis.run_test.periodic.MC_velo_alt1" u ($1/0.25):($4/64) w l lt 1 color yellow ti "MC velo alt1", \
#      "analysis.run_test.periodic.MC_velo_microcanonical" u ($1/0.25):($4/64) w l lt 1 color purple ti "MC velo microcanonical", \
#      "analysis.run_test.periodic.MD_orig" u ($1/0.25):($4/64) w l lt 1 lw .5 color blue ti "MD orig", \
#      "analysis.run_test.periodic.MD_alt1" u ($1/0.25):($4/64) w l lt 1 lw .5 color green ti "MD alt1", \
#      "analysis.run_test.periodic.MD_alt2" u ($1/0.25):($4/64) w l lt 1 lw .5 color orange ti "MD alt2", \
#      "analysis.run_test.periodic.MD_alt3" u ($1/0.25):($4/64) w l lt 1 lw .5 color grey ti "MD alt3", \
#      "analysis.run_test.periodic.MD_long" u ($1/0.25):($4/64) w l lt 1 lw .5 color cyan ti "MD long", \
#      "analysis.run_test.periodic.MD_bigsteps" u ($1/0.25):($4/64) w l lt 2 lw .5 color red ti "MD big steps", \
#      "analysis.run_test.periodic.MD_bigsteps2" u ($1/0.25):($4/64) w l lt 2 lw .5 color green ti "MD big steps 2", \
#      "analysis.run_test.periodic.MC_velo_MD" u ($1/0.25):($4/64) w l lt 2 lw .5 color purple ti "MC velo MD"
# unset log x
# unset log y
# set autoscale x
# set autoscale y
# 
# set output "periodic.energies.pdf"
# set xlabel "n iter / n walkers"
# set ylabel "H / epsilon"
# unset log y
# set yrange [:15]
# plot "< cat run_test.periodic.MC/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 color red ti "MC", \
#      "< cat run_test.periodic.MC_velo/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 color magenta ti "MC velo", \
#      "< cat run_test.periodic.MC_velo_alt1/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 color yellow ti "MC velo alt1", \
#      "< cat run_test.periodic.MC_velo_microcanonical/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 color purple ti "MC velo microcanonical", \
#      "< cat run_test.periodic.MD_orig/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color blue ti "MD orig", \
#      "< cat run_test.periodic.MD_alt1/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color green ti "MD alt1", \
#      "< cat run_test.periodic.MD_alt2/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color orange ti "MD alt2", \
#      "< cat run_test.periodic.MD_alt3/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw .5 color grey ti "MD alt3", \
#      "< cat run_test.periodic.MD_long/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 1 lw 0.5 color cyan ti "MD long", \
#      "< cat run_test.periodic.MD_bigsteps/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 2 lw 0.5 color red ti "MD big steps", \
#      "< cat run_test.periodic.MD_bigsteps2/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 2 lw 0.5 color green ti "MD big steps 2", \
#      "< cat run_test.periodic.MC_velo_MD/*energies | tail -n +2" every 100 u ($1/128):($2/64) w l lt 2 lw .5 color purple ti "MC velo MD", \
#      -7.75/64 ti "MC velo E at 90000 iters"
# unset log x
# unset log y
# set autoscale x
# set autoscale y
# 
# set output "periodic.volumes.pdf"
# set xlabel "n iter / n walkers"
# set ylabel "V / sigma cubed"
# unset log y
# # set yrange [:15]
# plot "proc_traj.run_test.periodic.MC" every 100 u ($1/128):($3/64) w l lt 1 color red ti "MC", \
#      "proc_traj.run_test.periodic.MC_velo" every 100 u ($1/128):($3/64) w l lt 1 color magenta ti "MC velo", \
#      "proc_traj.run_test.periodic.MC_velo_alt1" every 100 u ($1/128):($3/64) w l lt 1 color yellow ti "MC velo alt1", \
#      "proc_traj.run_test.periodic.MC_velo_microcanonical" every 100 u ($1/128):($3/64) w l lt 1 color purple ti "MC velo microcanonical", \
#      "proc_traj.run_test.periodic.MD_orig" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color blue ti "MD orig", \
#      "proc_traj.run_test.periodic.MD_alt1" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color green ti "MD alt1", \
#      "proc_traj.run_test.periodic.MD_alt2" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color orange ti "MD alt2", \
#      "proc_traj.run_test.periodic.MD_alt3" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color grey ti "MD alt3", \
#      "proc_traj.run_test.periodic.MD_long" every 100 u ($1/128):($3/64) w l lt 1 lw .5 color cyan ti "MD long", \
#      "proc_traj.run_test.periodic.MD_bigsteps" every 100 u ($1/128):($3/64) w l lt 2 lw .5 color red ti "MD big steps", \
#      "proc_traj.run_test.periodic.MD_bigsteps2" every 100 u ($1/128):($3/64) w l lt 2 lw .5 color green ti "MD big steps 2", \
#      "proc_traj.run_test.periodic.MC_velo_MD" every 100 u ($1/128):($3/64) w l lt 2 lw .5 color purple ti "MC velo MD"
# unset log x
# unset log y
# set autoscale x
# set autoscale y
# 
# set output "periodic.E_vs_V.pdf"
# set xlabel "V / sigma cubed"
# set ylabel "H / epsilon"
# set yrange [:1000]
# plot "proc_traj.run_test.periodic.MC" every 100 u 3:2 w l lt 1 color red ti "MC", \
#      "proc_traj.run_test.periodic.MC_velo" every 100 u 3:2 w l lt 1 color magenta ti "MC velo", \
#      "proc_traj.run_test.periodic.MD_orig" every 100 u 3:2 w l lt 1 color cyan ti "MD orig"
# unset log x
# unset log y
# set autoscale x
# set autoscale y
# 
# 
# set output "periodic.E_vs_V.log.pdf"
# set xlabel "V / sigma cubed"
# set ylabel "H / epsilon"
# set yrange [:1000]
# set log x
# plot "proc_traj.run_test.periodic.MC" every 100 u 3:2 w l lt 1 color red ti "MC", \
#      "proc_traj.run_test.periodic.MC_velo" every 100 u 3:2 w l lt 1 color magenta ti "MC velo", \
#      "proc_traj.run_test.periodic.MD_orig" every 100 u 3:2 w l lt 1 color cyan ti "MD orig"
