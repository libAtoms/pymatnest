set term pdf

set output "cluster.Cv.pdf"
set title "LJ6 cluster, sigma=3"
set xlabel "T / epsilon"
set ylabel "Cv per atom"
set autoscale y
plot "analysis.cluster.MC.Livia" u ($1/0.1):($4/6) w l lt 1 color black ti "Livia reference", \
     "analysis.test.cluster.MC.energies" u ($1/0.25):($4/6) w l lt 1 color red ti "MC", \
     "analysis.test.cluster.MD.energies" u ($1/0.25):($4/6) w l lt 1 color cyan ti "MD"
unset log x
unset log y
set autoscale x
set autoscale y

set output "cluster.energies.pdf"
set xlabel "n iter / n walkers"
set ylabel "E / atom epsilon"
set yrange [:4]
plot "< tail -n +2 test.cluster.MC.energies" u ($1/1024):($2/6/0.25) w l lt 1 color red ti "MC", \
     "< tail -n +2 test.cluster.MD.energies" u ($1/1024):($2/6/0.25) w l lt 1 color cyan ti "MD"
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
set xrange [0:1.5]
plot "analysis.periodic.MC.Livia" u ($1/0.1):($4/64) w l lt 1 color black ti "Livia reference", \
     "analysis.test.periodic.MC.energies" u ($1/0.25):($4/64) w l lt 1 color red ti "MC", \
     "analysis.test.periodic.MD.energies" u ($1/0.25):($4/64) w l lt 1 color cyan ti "MD"
unset log x
unset log y
set autoscale x
set autoscale y

set output "periodic.energies.pdf"
set xlabel "n iter / n walkers"
set ylabel "H / atom epsilon"
unset log y
set yrange [:60]
plot "< tail -n +2 test.periodic.MC.energies" u ($1/128):($2/64/0.25) w l lt 1 color red ti "MC", \
     "< tail -n +2 test.periodic.MD.energies" u ($1/128):($2/64/0.25) w l lt 1 color cyan ti "MD"
unset log x
unset log y
set autoscale x
set autoscale y

set output "periodic.volumes.pdf"
set xlabel "n iter / n walkers"
set ylabel "H / atom epsilon"
unset log y
set yrange [:60]
plot "< tail -n +2 test.periodic.MC.energies" u ($1/128):($2/64/0.25) w l lt 1 color red ti "MC", \
     "< tail -n +2 test.periodic.MD.energies" u ($1/128):($2/64/0.25) w l lt 1 color cyan ti "MD"
unset log x
unset log y
set autoscale x
set autoscale y
