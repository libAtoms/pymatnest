set term pdf
set output "MD_vs_time.pdf"
set xlabel "$n_\mathrm{iter}$"
set ylabel "v"
set y2label "time step"
set log y2
plot "proc_traj.run_test.periodic.MD_orig" every 100 u 1:(sqrt($4/64)) w l ti "sqrt(KE/N)", \
     "< fgrep adjusted run_test.periodic.MD_orig/*out | grep MD_atom" u 2:8 axes x1y2 w l ti "time step"
