# set palette defined   ( \
# 0.00  0.50  1.00  1.00, \
# 0.10  0.00  1.00  1.00, \
# 0.25  0.00  0.00  1.00, \
# 0.45  0.00  0.00  0.50, \
# 0.50  0.00  0.00  0.00, \
# 0.55  0.50  0.00  0.00, \
# 0.75  1.00  0.00  0.00, \
# 0.90  1.00  1.00  0.00, \
# 1.00  1.00  1.00  0.50 )

# set palette defined   ( \
# 0.00  0.00  0.00  0.50, \
# 0.15  0.00  0.00  1.00, \
# 0.30  0.00  0.50  1.00, \
# 0.45  0.00  1.00  1.00, \
# 0.50  1.00  1.00  1.00, \
# 0.55  1.00  1.00  0.00, \
# 0.70  1.00  0.50  0.00, \
# 0.85  1.00  0.00  0.00, \
# 1.00  0.50  0.00  0.00 )

set palette defined   ( \
0.00  0.00  0.00  0.50, \
0.10  0.00  0.50  1.00, \
0.17  0.00  1.00  1.00, \
0.25  0.00  1.00  0.50, \
0.34  0.50  1.00  0.00, \
0.46  1.00  1.00  0.00, \
0.61  1.00  0.50  0.00, \
1.00  0.50  0.00  0.00 )

# average_ALL = (($18**2+$19**2+$20**2+$21**2+$22**2+$23**2)/3)**(0.5)
# average_ww = ((($3+$4+$5)-($6+$7+$8))/3)**(0.5)

# set terminal png size 3840,2880
#set term eps
#set term svg
###Setting terminal type

set term png transparent truecolor size 3840,2304
set lmargin 30
set border lw 18
unset colorbox





###General
set pm3d
set pm3d interpolate 8,8
set view map
set hidden3d
# set cbrange[0:1e-1]
# set xtics 0.02
# set ytics 0.02



CBTicsFont="Arial,72"
CBLabelFont="Arial,96"

#DipoleFT_ALL
# set output "NMA_DipoleFT_ALL.png"
set xtics 500
set ytics 0.04
set cbrange[0:7e-3]
# set ytics font "Arial, 120" offset -1,0,0
# set xtics font "Arial, 120" offset 1,-7,0
# set xlabel "Time delay (a.u)"      font "Arial, 156" offset  0,-20.0,0
# set ylabel "Emission energy (a.u)" font "Arial, 156" offset  -34.0,0,0
# Average
# splot [][0.2:0.48] 'sim1/Dipole/DipoleFT_ALL' u 9:17:((($18**2+$19**2+$20**2+$21**2+$22**2+$23**2)/3)**(0.5)) w l notitle

#DipoleFT_ww
# set output "NMA_DipoleFT_ww.png"
# set xtics 0.04
# set ytics 0.04
# set cbrange[0:1]
# set ytics font "Arial, 120" offset -1,0,0
# set xtics font "Arial, 120" offset 1,-7,0
# set xlabel "Excitation energy (a.u)"      font "Arial, 26" offset  0,-0.5,0
# set ylabel "Emission energy (a.u)" font "Arial, 26" offset  -4.5,0,0
#Single Axis
#Average
# splot [0.2:0.48][0.2:0.48] 'sim/Dipole/DipoleFT_ww' u ($1+$2):1:((($3**2+$4**2+$5**2+$6**2+$7**2+$8**2)/3)**(0.5)) w l notitle



###Dipole Response
##DipoleFT_ALL
# set output "NMA_DipoleFT_ALL.png"
# set cbrange[0:7e-3]
# set ytics font "Arial, 20" offset 1,0,0
# set xtics font "Arial, 20" offset 1,0,0
# set xlabel "Time delay (a.u)"      font "Arial, 26" offset  0,-0.5,0
# set ylabel "Emission energy (a.u)" font "Arial, 26" offset  -4.5,0,0
#splot [-600:1600][0.28:0.62] 'sim/Dipole/DipoleFT_ALL' u 9:17:(($20**2)+($21**2)) w l notitle

# splot [][0.2:0.48] 'sim/Dipole/DipoleFT_ALL' u 9:17:((($18**2+$19**2+$20**2+$21**2+$22**2+$23**2)/3)**(0.5)) w l notitle


# set cbrange[-2.5e-3:2.5e-3]
#set output "NMA_Reconstructed_Dipole_ww.png"
#set xtics font "Arial, 20"
#set ytics font "Arial, 20" offset 1,0,0
#set xlabel "Excitation energy (a.u)" font "Arial, 26" offset  0,-0.2,0
#set ylabel "Emission energy (a.u)"   font "Arial, 26" offset  -0.4,0,0
#splot [0.35:0.55][0.28:0.62]'difference_sim' u ($2+$1):1:($7) w l notitle

##DipoleFT_ww
#set cbrange [-5:5]
#set output "NMA_DipoleFT_ww.png"
#set xtics font "Arial, 20"
#set ytics font "Arial, 20" offset 1,0,0
#set xlabel "Excitation energy (a.u)" font "Arial, 26" offset  0,-0.2,0
#set ylabel "Emission energy (a.u)"   font "Arial, 26" offset  -0.4,0,0
#splot [0.35:0.55][0.28:0.62] 'sim/Dipole/DipoleFT_ww' u ($2+$1):1:(($5**2)+($6**2)) w l notitle


##Atomic Charge
#AtomicChargeFT_ww
#set cbrange [-5:5]
#set xtics font "Arial, 20"
#set ytics font "Arial, 20" offset 1,0,0
#set xlabel "Excitation energy (a.u)" font "Arial, 26" offset  0,-0.2,0
#set ylabel "Emission energy (a.u)"   font "Arial, 26" offset  -0.4,0,0
#set output "NMA_charges_ww.png"

###########
# set cbrange [-5:5]
#set logscale
#splot [0.35:0.55][0.28:0.62] 'difference' u ($2+$1):1:($7) w l notitle
#splot [0.35:0.55][0.28:0.62] 'difference' u ($2+$1):1:($4) w l notitle
#splot [0.35:0.55][0.28:0.62] 'difference' u ($2+$1):1:($4 - 2.8 *$7) w l notitle #y pol
#splot [0.35:0.55][0.28:0.62] 'difference' u ($2+$1):1:(($3+$4+$5)/3 - 2.4 *($6+$7+$8)/3) w l notitle #average

#splot [0.35:0.55][0.28:0.62] 'sim/Dipole/DipoleFT_ww_reconstructed' u ($1+$2):1:((($5**2+$6**2)>0.01?1:0)*atan2($6,$5)) w l notitle
#splot [0.35:0.55][0.28:0.62] 'sim/Dipole/DipoleFT_ww' u ($1+$2):1:((($5**2+$6**2)>0.01?1:0)*atan2($6,$5)) w l notitle
##########


#set output "NMA_O_ww.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($2+$1):1:(($3**2)+($4**2)) w l notitle

#set output "NMA_N_ww.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($2+$1):1:(($5**2)+($6**2)) w l notitle

#set output "NMA_C_ww.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($2+$1):1:(($7**2)+($8**2)) w l notitle

#set output "NMA_RC_ww.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($2+$1):1:(($9**2)+($10**2)) w l notitle

#set output "NMA_LC_ww.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($2+$1):1:(($11**2)+($12**2)) w l notitle



###1d
#set output "Dipoleresponse.png"
#plot [-230:1700] [-0.014:0.014] '../Dipole/DipolePP1000.0' u 2:5 w l

#set output "Pulsesequenceblue.svg"
#plot [-300:1800] [-8e-17:8e-17] 'pulsePP1000.0' u 1:2 w l lc "blue"
#set output "Pulsesequencered.svg"
#plot [-300:1800] [-8e-17:8e-17] 'pulsePP1000.0' u 1:2 w l lc "red"


#Phase peak A
# set cbrange [-5:5]

#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/ Dipoleresponse.png"
#splot  'sim_test/Dipole/DipoleFT_ww' u ($1+$2):1:($5**2 + $6**2) w l

#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/ Dipole.png"
#splot [0.36:0.40][0.42:0.50] 'sim_test/Dipole/DipoleFT_ww' u ($1+$2):1:((($5**2+$6**2)>0.7?1:0)*atan2($6,$5))

#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/ O.png"
#splot [0.36:0.40][0.42:0.50] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($3**2+$4**2)>0.01?1:0)*atan2($4,$3)) w l

#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/S N.png"
#splot [0.36:0.40][0.42:0.50] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($5**2+$6**2)>0.00015?1:0)*atan2($6,$5)) w l

#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/ C.png"
#splot [0.36:0.40][0.42:0.50] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($7**2+$8**2)>0.01?1:0)*atan2($8,$7)) w l

#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/ RC.png"
#splot [0.36:0.40][0.42:0.50] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($9**2+$10**2)>0.000185?1:0)*atan2($10,$9)) w l

#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/ LC.png"
#splot [0.36:0.40][0.42:0.50] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($11**2+$12**2)>0.00006?1:0)*atan2($12,$11)) w l


##Phase peak ALL
#set cbrange [-5:5]
#
##set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/A/ Dipoleresponse.png"
#splot  'sim_test/Dipole/DipoleFT_ww' u ($1+$2):1:($5**2 + $6**2) w l
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ Dipole.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/Dipole.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/Dipole/DipoleFT_ww' u ($1+$2):1:((($5**2+$6**2)>0.7?1:0)*atan2($6,$5))
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ O.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/Oresclaed.png"
##splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($3**2+$4**2)>0.01?1:0)*atan2($4,$3)) w l notitle
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($3**2+$4**2)>0.92?1:0)*atan2($4,$3)) w l notitle
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/O.png
#
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ N.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/N.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($5**2+$6**2)>0.00015?1:0)*atan2($6,$5)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ C.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/C.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($7**2+$8**2)>0.01?1:0)*atan2($8,$7)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ RC.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/RC.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($9**2+$10**2)>0.000185?1:0)*atan2($10,$9)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ LC.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/LC.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($11**2+$12**2)>0.00006?1:0)*atan2($12,$11)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ H1.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/H1.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($13**2+$14**2)>0.00074?1:0)*atan2($14,$13)) w l notitle #still off
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ H2.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/H2.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($15**2+$16**2)>0.00074?1:0)*atan2($16,$15)) w l notitle #still off
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ H3.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/H3.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($17**2+$18**2)>0.000024?1:0)*atan2($18,$17)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ H4.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/H4.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($19**2+$20**2)>0.00018?1:0)*atan2($20,$19)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ H5.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/H5.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($21**2+$22**2)>0.000025?1:0)*atan2($22,$21)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ H6.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/H6.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($23**2+$24**2)>0.000025?1:0)*atan2($24,$23)) w l notitle
#
#set output "../../../../../../../Desktop/Conferences/SURE/NMA/Phase\ comparison/F/ H7.png"
##set output "../../../../../home/../../../Desktop/Picture\ collage\ c\,/phasespectrum/H7.png"
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($25**2+$26**2)>0.000025?1:0)*atan2($26,$25)) w l notitle
#

#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($25**2+$26**2)>0.7?1:0)*atan2($26,$25)+(($25**2+$26**2)>0.7?1:0)*atan2($26,$25)+(($23**2+$24**2)>0.7?1:0)*atan2($24,$23)+(($21**2+$22**2)>0.7?1:0)*atan2($22,$21)+(($19**2+$20**2)>0.7?1:0)*atan2($20,$19)+(($17**2+$18**2)>0.7?1:0)*atan2($18,$17)+(($15**2+$16**2)>0.7?1:0)*atan2($16,$15)+(($13**2+$14**2)>0.7?1:0)*atan2($14,$13)+(($11**2+$12**2)>0.7?1:0)*atan2($12,$11)+(($9**2+$10**2)>0.7?1:0)*atan2($10,$9)+(($7**2+$8**2)>0.7?1:0)*atan2($8,$7)+(($5**2+$6**2)>0.7?1:0)*atan2($6,$5)+(($3**2+$4**2)>0.7?1:0)*atan2($4,$3)) w l notitle
#
#
#
#
#
#splot [0.35:0.55][0.28:0.62] 'sim_test/AtomicCharge/AtomicChargeFT_ww' u ($1+$2):1:((($25**2+$26**2+$23**2+$24**2+$21**2+$22**2+$19**2+$20**2+$17**2+$18**2+$15**2+$16**2+$13**2+$14**2+$11**2+$12**2+$9**2+$10**2+$7**2+$8**2+$5**2+$6**2+$3**2+$4**2)>0.7?1:0)*atan2($26+$24+$22+$20+$18+$16+$14+$12+$10+$8+$6+$4,$25+$3+$5+$7+$9+$11+$13+$15+$17+$19+$21+$23+$25)) w l notitle
#+(()>0.7?1:0)*atan2($26,$25)
#+(()>0.7?1:0)*atan2($24,$23)
#+(()>0.7?1:0)*atan2($22,$21)
#+(()>0.7?1:0)*atan2($20,$19)
#+(()>0.7?1:0)*atan2($18,$17)
#+(()>0.7?1:0)*atan2($16,$15)
#+(()>0.7?1:0)*atan2($14,$13)
#+(()>0.7?1:0)*atan2($12,$11)
#+(()>0.7?1:0)*atan2($10,$9)
#+(()>0.7?1:0)*atan2($8,$7)
#+(()>0.7?1:0)*atan2($6,$5)
#+(()>0.7?1:0)*atan2($4,$3))
#
#
#
#
#splot [0.35:0.55][0.28:0.62] 'sim_test/Dipole/DipoleFT_ww' u ($1+$2):1:((($5**2+$6**2+$3**2+$4**2)>0.7?1:0)*atan2($6+$4,$5+$3))
#
#
#set output "../../../../../../../Desktop/Research/2DMAPS/X_Dipole_"
#splot 'sim_test/Dipole/DipoleFT_ww' u ($2+$1):1:($3**2+$4**2) w l
#set output "../../../../../../../Desktop/Research/2DMAPS/X_Dipole_reconstructed"
#splot 'sim_test/Dipole/DipoleFT_ww_reconstructed' u ($2+$1):1:($3**2+$4**2) w l
#set output "../../../../../../../Desktop/Research/2DMAPS/Y_Dipole_"
#splot 'sim_test/Dipole/DipoleFT_ww' u ($2+$1):1:($5**2+$6**2) w l
#set output "../../../../../../../Desktop/Research/2DMAPS/Y_Dipole_reconstructed"
#splot 'sim_test/Dipole/DipoleFT_ww_reconstructed' u ($2+$1):1:($5**2+$6**2) w l


