set term pdf

set output "cluster.Cv.pdf"
set title "LJ6 cluster, sigma=3"
set xlabel "T / epsilon"
set ylabel "Cv per atom"
set autoscale y
plot "analysis.cluster.MC.Livia" u ($1/0.1):($4/6) w l lt 1 color black ti "Livia reference", \
     "test.cluster.MD.fortran.energies.analysis" u ($1/0.25):($4/6) w l lt 1 color red ti "MD fortran", \
     "test.cluster.MD.quip.energies.analysis" u ($1/0.25):($4/6) w l lt 1 color cyan ti "MD quip", \
     "test.cluster.MC.fortran.energies.analysis" u ($1/0.25):($4/6) w l lt 1 color blue ti "MC fortran", \
     "test.cluster.MC.quip.energies.analysis" u ($1/0.25):($4/6) w l lt 1 color green ti "MC quip"
unset log x
unset log y
set autoscale x
set autoscale y


set output "periodic.Cp.pdf"
set title "LJ64 periodic, sigma=3"
set xlabel "T / epsilon"
set ylabel "Cp per atom"
set log y
set yrange [1:]
set key left top
set xrange [0:0.40]
plot "analysis.periodic.MC.Baldock" u ($1*8.6173324e-5/0.4):($2/8.6e-5/64) w l lt 1 color black ti "Rob reference", \
     "test.periodic.MD.fortran.energies.analysis" u ($1/0.4):($4/64) w l lt 1 color red ti "MD fortran", \
     "test.periodic.MD.quip.energies.analysis" u ($1/0.4):($4/64) w l lt 1 color cyan ti "MD quip"
unset log x
unset log y
set autoscale x
set autoscale y
