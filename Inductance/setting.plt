set datafile separator ","
set terminal gif animate delay 2 size 400,800
unset key
set ticslevel 0
set xrange[-0.02:0.02]
set yrange[-0.02:0.02]
set zrange[0:0.08]
set view equal xyz
set output 'rotate.gif'
do for [i=0:360]{
    set view 70,i,1
    splot 'output.csv' using 1:2:3:4:5:6 with vector
}
set output
