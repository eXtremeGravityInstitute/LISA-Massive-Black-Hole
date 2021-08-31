set term png enhanced truecolor crop font Helvetica 26  size 1600,800
#set term png enhanced truecolor crop font Helvetica 18  size 2400,800
set output "Qscan.png"
set pm3d map corners2color c1
set ylabel "f (Hz)"
set xlabel "t (s)"
set logscale y
set xrange [1e5:2571440]
set yrange [1e-5:0.1]
set ytics (1e-5,1e-4,1e-3,1e-2,0.1)
set cbrange [0:9]
# colorbrewer palette colors
set palette defined (0 '#fff7ec', 1'#fee8c8', 2 '#fdd49e', 3 '#fdbb84', 4 '#fc8d59', 5 '#ef6548', 6 '#d7301f',7 '#990000')
splot "Qtransform.dat" using 1:2:3 notitle




