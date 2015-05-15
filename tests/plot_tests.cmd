set term pdf

set output "cluster.Cv.pdf"
set title "LJ6 cluster, sigma=3"
set xlabel "T / epsilon"
set ylabel "Cv per atom"
set autoscale y
plot "cluster.MC.Livia.analysis" u ($1/0.1):($4/6) w l lt 1 color black ti "Livia reference", \
     "test.cluster.MC.energies.analysis" u ($1/0.25):($4/6) w l lt 1 color red ti "MC", \
     "test.cluster.MD.energies.analysis" u ($1/0.25):($4/6) w l lt 1 color cyan ti "MD"
set autoscale x
set autoscale y

set output "cluster.energies.pdf"
set xlabel "n iter / n walkers"
set ylabel "E / epsilon"
set yrange [:1]
plot "< tail -n +2 test.cluster.MC.energies" u ($1/1024):($2/6) w l lt 1 color red ti "MC", \
     "< tail -n +2 test.cluster.MD.energies" u ($1/1024):($2/6) w l lt 1 color cyan ti "MD"
set autoscale x
set autoscale y


set output "periodic.Cv.pdf"
set title "LJ64 periodic, sigma=3"
set xlabel "T / epsilon"
set ylabel "Cv per atom"
set log y
set yrange [1:]
set xrange [0:1.5]
plot "periodic.MC.Livia.analysis" u ($1/0.1):($4/64) w l lt 1 color black ti "Livia reference", \
     "test.periodic.MC.energies.analysis" u ($1/0.25):($4/64) w l lt 1 color red ti "MC", \
     "test.periodic.MC_velo.energies.analysis" u ($1/0.25):($4/64) w l lt 1 color magenta ti "MC velo"
set autoscale x
set autoscale y

set output "periodic.energies.pdf"
set xlabel "n iter / n walkers"
set ylabel "E / epsilon"
unset log y
set yrange [:15]
plot "< tail -n +2 test.periodic.MC.energies" u ($1/128):($2/64) w l lt 1 color red ti "MC", \
     "< tail -n +2 test.periodic.MC_velo.energies" u ($1/128):($2/64) w l lt 1 color magenta ti "MC velo"
