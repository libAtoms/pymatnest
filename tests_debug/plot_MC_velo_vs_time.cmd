# os.system('for n in `seq 1 20`; do head -$(( 10000 * $n )) < run_test.periodic.MC_velo/*.energies > t; ../ns_analyse -M 0.002 -D 0.0015625 -n 400 -k 1.0 t; echo ""; echo ""; done > MC_velo_analysis_vs_time')

set term pdf
set output "MC_velo_analysis_vs_time.pdf"
set xlabel "$T/\epsilon$"
set ylabel "$C_p$/atom"
set log y
plot "MC_velo_analysis_vs_time" i 0 u ($1/0.25):($4/64) w l ti "10", \
     "MC_velo_analysis_vs_time" i 1 u ($1/0.25):($4/64) w l ti "20", \
     "MC_velo_analysis_vs_time" i 2 u ($1/0.25):($4/64) w l ti "30", \
     "MC_velo_analysis_vs_time" i 3 u ($1/0.25):($4/64) w l ti "40", \
     "MC_velo_analysis_vs_time" i 4 u ($1/0.25):($4/64) w l ti "50", \
     "MC_velo_analysis_vs_time" i 5 u ($1/0.25):($4/64) w l ti "60", \
     "MC_velo_analysis_vs_time" i 6 u ($1/0.25):($4/64) w l ti "70", \
     "MC_velo_analysis_vs_time" i 7 u ($1/0.25):($4/64) w l ti "80", \
     "MC_velo_analysis_vs_time" i 8 u ($1/0.25):($4/64) w l ti "90", \
     "MC_velo_analysis_vs_time" i 9 u ($1/0.25):($4/64) w l lt 2 color red ti "100", \
     "MC_velo_analysis_vs_time" i 11 u ($1/0.25):($4/64) w l lt 2 color blue ti "120", \
     "MC_velo_analysis_vs_time" i 13 u ($1/0.25):($4/64) w l lt 2 color magenta ti "140", \
     "MC_velo_analysis_vs_time" i 15 u ($1/0.25):($4/64) w l lt 2 color cyan ti "160"
set nolog y

####################################################################################################

set output "MC_velo_vs_time.pdf"

set xlabel "$n_\mathrm{iter}$"

set ylabel "$H/\epsilon$"

set y2label "$V/\sigma^3$"
set log y2

set y3label "step size"
set y3range [0:6]

set arrow 0 from first 0, axis3 5 to first 81000, axis3 5 with nohead color blue
set arrow 1 from first 81000, axis3 5 to first 81000,axis3 2.5   with nohead color blue 
plot "proc_traj.run_test.periodic.MC_velo" every 100 u 1:2 axes x1y1 w l ti "$H$", \
     "proc_traj.run_test.periodic.MC_velo" every 100 u 1:3 w l axes x1y2 ti "$V$", \
     "< egrep 'adju.*MC_atom ' run_test.periodic.MC_velo/*out" u 2:8 axes x1y3 w l ti "step size"
