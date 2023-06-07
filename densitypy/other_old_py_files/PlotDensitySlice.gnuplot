# 
#   Possible commands one may want to use to generate the input files
#   and run gnuplot recursively
#   find ChDen*   -exec  Split3DGrid -i {} -o Z0{} -n 4 \;
#   find ChDen*Z0 -exec gnuplot -e file=\'{}\' PlotDensitySlice.gnuplot \;
#   ls ChDen*-*.png | sort -r > list.txt
#   ls ChDen[0-9,_]*.png >> list.txt   <<< not OK
#   mencoder mf://@listfiles.txt -mf w=600:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o CO2ChargeMigration.avi


#Instructions (Right Now)
#1. Go inside simulation directory
#2. Run $ ./testing 
#	This will make a file containing the list of the names for timestep files in order, Convert the 3d grid to 2d by summing over the chosen orientation, and finally call gnuplot script to create png files.
#3. Run $ mencoder mf://@listfiles.txt -mf w=600:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o CO2ChargeMigration.avi

unset key
set pm3d 
set pm3d interpolate 8,8
set view map
set hidden3d
#file='SliceZ0'
coord1=1
coord2=2
coord3=4
c1min =-4.0
c1max = 4.0
c2min =-4.0
c2max = 4.0
cbmax =2e-5

set palette defined   ( \
0.00  0.50  1.00  1.00, \
0.10  0.00  1.00  1.00, \
0.25  0.00  0.00  1.00, \
0.45  0.00  0.00  0.50, \
0.50  0.00  0.00  0.00, \
0.55  0.50  0.00  0.00, \
0.75  1.00  0.00  0.00, \
0.90  1.00  1.00  0.00, \
1.00  1.00  1.00  0.50 )

set palette defined   ( \
0.00  0.00  0.00  0.50, \
0.15  0.00  0.00  1.00, \
0.30  0.00  0.50  1.00, \
0.45  0.00  1.00  1.00, \
0.50  1.00  1.00  1.00, \
0.55  1.00  1.00  0.00, \
0.70  1.00  0.50  0.00, \
0.85  1.00  0.00  0.00, \
1.00  0.50  0.00  0.00 )

set xrange  [ c1min:c1max]
set yrange  [ c2min:c2max]
set cbrange [-cbmax:cbmax]

unset colorbox
unset tics
set border lw 0.1
ratyx=(c2max-c2min)/(c1max-c1min)
set size ratio ratyx
nxpts=600
nypts=nxpts*ratyx
set terminal png truecolor size nxpts,nypts crop
set out file.'.png'

splot file u (column(coord1)):(column(coord2)):(column(coord3)) w l
