
  �0  r   k820309    �          2021.10.0   䝶d                                                                                                          
       src/Module_CD_IO.f90 MODULE_CD_IO              SET_CD_IO_VERBOUS LOADGEOMETRY LOADENERGIES LOADGRID LOADORBITALS LOADDIPOLEME LOADTDMS LOADDIPOLEMO Write_Dipole Write_Q_Charge Write_Dipole1 Write_Q_Charge1 Write_Charge_Density WRITE_WEIGHTS READ_WEIGHTS WRITE_SUMMARY            0                                                
                            @                                      
                      �   @                                       '                    #TEAM_PTR    #VERSION_NUMBER    #CHECK_VALID                � D                                                                #C_PTR                                                                                    X#C_PTR       n                                      0                                               �   @                                      '                    #PTR                 � D                                                                     � D                                                                     � D                                                                      @   @                                 	     'P              
      #ENVELOPE 
   #T    #O    #F 
   #D    #I    #A    #P    #THETA    #PHI                 �                                        
                                �                                                      
                �                                                      
                �                                       
               
                �                                                       
                �                                            (          
                �                                            0          
                �                                            8          
                �                                            @       	   
                �                                            H       
   
   #         @                                                                #LOGI              
                                                      #         @                                                                #NATOMS    #ATCOORD    #FILENAME    #ATOMNAME              D                                                              D                                                            
               &                   &                                                     
  @                                                          1 ,        D                                                                           &                                                   #         @                                                                #FILENAME    #NSTATES    #EVEC              
                                                             1           D                                                              D                                                            
               &                                           #         @                                                                #FILENAME     #NPTS !   #GRIDV "             
                                                              1           D                                          !                    D                                         "                   
               &                   &                                           #         @                                            #                    #DIR $   #NORB %   #IVORB &   #NPTS '   #ORBTAB (             
                                         $                    1           
                                           %                     
                                           &                                 &                                                     
                                           '                   D                                         (                   
 	              &                   &                                           #         @                                            )                    #DMAT *   #INPDIR +   #NSTATES ,           D @                                       *                   
 
              &                   &                   &                                                     
                                         +                    1           
  @                                        ,           #         @                                            -                    #FILENAMEDM .   #FILENAMETDM /   #NSTATES 0   #NORB 1   #TDM 2             
                                         .                    1           
                                         /                    1           
  @                                        0                     D                                          1                    D                                         2                   
 
              &                   &                   &                   &                                           #         @                                            3                    #INPDIR 4   #NORB 5   #IVORB 6   #MUORB 7             
                                         4                    1           
                                           5                     
                                           6                                 &                                                   D                                         7                   
               &                   &                                           #         @                                            8                    #FILENAME 9   #DIPOLE :   #NTIMES ;   #TMIN <   #DT =             
                                         9                    1           
                                          :                   
              &                   &                                                     
                                           ;                     
                                          <     
                
                                          =     
      #         @                                            >                    #FILENAME ?   #CHARGE @   #NTIMES A   #TMIN B   #DT C   #NATOMS D             
                                         ?                    1           
                                          @                   
              &                   &                                                     
                                           A                     
                                          B     
                
                                          C     
                
                                           D           #         @                                            E                    #FILENAME F   #DIPOLE G   #NTIMES H   #TMIN I   #DT J             
                                         F                    1           
                                          G                   
              &                   &                                                     
                                           H                     
                                          I     
                
                                          J     
      #         @                                            K                    #FILENAME L   #CHARGE M   #NTIMES N   #TMIN O   #DT P   #NATOMS Q             
                                         L                    1           
                                          M                   
              &                   &                   &                                                     
                                           N                     
                                          O     
                
                                          P     
                
                                           Q           #         @                                            R                    #FILENAME S   #NPTS T   #GRIDV U   #CHDEN V   #WEIGHTV W   #NATOMS X             
                                         S                    1           
                                           T                   
                                          U                   
              &                   &                                                   
                                          V                   
              &                                                   
                                          W                   
              &                   &                                                     
                                           X           #         @                                            Y                    #FILENAME Z   #WEIGHTV [   #GRIDV \   #NATOMS ]   #NPTS ^             
                                         Z                    1           
                                          [                   
              &                   &                                                     
                                          \                   
              &                   &                                                     
                                           ]                     
                                           ^           #         @                                            _                    #FILENAME `   #WEIGHTV a   #NATOMS b   #NPTS c             
                                         `                    1         D                                         a                   
               &                   &                                                     
                                           b                     
                                           c           #         @                                            d                    #FILENAME e   #NPTS f   #NATOMS g   #VOLUME h   #COMPUTED_VOLUME i   #NTIMES j   #TMIN k   #TMAX l   #ATOMNAME m   #RADIUS_BS n   #NORB o   #ORBTAB p             
                                         e                    1           
                                           f                     
                                           g                     
                                          h     
                
                                          i     
                
                                           j                     
                                          k     
                
                                          l     
      ,          
                                          m                                 &                                                             
                                          n                   
              &                                                     
                                           o                     
                                          p                   
              &                   &                                              �   *      fn#fn "   �   �   b   uapp(MODULE_CD_IO     �  H   J  ISO_FORTRAN_ENV     �  H   J  MODULEPULSES_3D *   8  �      TEAM_TYPE+ISO_FORTRAN_ENV <   �    %   TEAM_TYPE%TEAM_PTR+ISO_FORTRAN_ENV=TEAM_PTR $   �  a      C_PTR+ISO_C_BINDING ,   ?  P   %   C_PTR%PTR+ISO_C_BINDING=PTR H   �  P   %   TEAM_TYPE%VERSION_NUMBER+ISO_FORTRAN_ENV=VERSION_NUMBER B   �  P   %   TEAM_TYPE%CHECK_VALID+ISO_FORTRAN_ENV=CHECK_VALID &   /  �      PULSE+MODULEPULSES_3D /   �  P   a   PULSE%ENVELOPE+MODULEPULSES_3D (   *  P   a   PULSE%T+MODULEPULSES_3D (   z  P   a   PULSE%O+MODULEPULSES_3D (   �  P   a   PULSE%F+MODULEPULSES_3D (     P   a   PULSE%D+MODULEPULSES_3D (   j  P   a   PULSE%I+MODULEPULSES_3D (   �  P   a   PULSE%A+MODULEPULSES_3D (   
  P   a   PULSE%P+MODULEPULSES_3D ,   Z  P   a   PULSE%THETA+MODULEPULSES_3D *   �  P   a   PULSE%PHI+MODULEPULSES_3D "   �  Z       SET_CD_IO_VERBOUS '   T	  H   a   SET_CD_IO_VERBOUS%LOGI    �	  �       LOADGEOMETRY $   !
  H   a   LOADGEOMETRY%NATOMS %   i
  �   a   LOADGEOMETRY%ATCOORD &     T   a   LOADGEOMETRY%FILENAME &   i  �   a   LOADGEOMETRY%ATOMNAME      u       LOADENERGIES &   z  T   a   LOADENERGIES%FILENAME %   �  H   a   LOADENERGIES%NSTATES "   
  �   a   LOADENERGIES%EVEC    �
  s       LOADGRID "     T   a   LOADGRID%FILENAME    q  H   a   LOADGRID%NPTS    �  �   a   LOADGRID%GRIDV    e  �       LOADORBITALS !   �  T   a   LOADORBITALS%DIR "   =  H   a   LOADORBITALS%NORB #   �  �   a   LOADORBITALS%IVORB "     H   a   LOADORBITALS%NPTS $   a  �   a   LOADORBITALS%ORBTAB    
  s       LOADDIPOLEME "   �  �   a   LOADDIPOLEME%DMAT $   D  T   a   LOADDIPOLEME%INPDIR %   �  H   a   LOADDIPOLEME%NSTATES    �  �       LOADTDMS $   q  T   a   LOADTDMS%FILENAMEDM %   �  T   a   LOADTDMS%FILENAMETDM !     H   a   LOADTDMS%NSTATES    a  H   a   LOADTDMS%NORB    �  �   a   LOADTDMS%TDM    �  |       LOADDIPOLEMO $     T   a   LOADDIPOLEMO%INPDIR "   U  H   a   LOADDIPOLEMO%NORB #   �  �   a   LOADDIPOLEMO%IVORB #   1  �   a   LOADDIPOLEMO%MUORB    �  �       Write_Dipole %   e  T   a   Write_Dipole%FILENAME #   �  �   a   Write_Dipole%DIPOLE #   e  H   a   Write_Dipole%NTIMES !   �  H   a   Write_Dipole%TMIN    �  H   a   Write_Dipole%DT    =  �       Write_Q_Charge '   �  T   a   Write_Q_Charge%FILENAME %   %  �   a   Write_Q_Charge%CHARGE %   �  H   a   Write_Q_Charge%NTIMES #     H   a   Write_Q_Charge%TMIN !   a  H   a   Write_Q_Charge%DT %   �  H   a   Write_Q_Charge%NATOMS    �  �       Write_Dipole1 &   y  T   a   Write_Dipole1%FILENAME $   �  �   a   Write_Dipole1%DIPOLE $   y  H   a   Write_Dipole1%NTIMES "   �  H   a   Write_Dipole1%TMIN     	   H   a   Write_Dipole1%DT    Q   �       Write_Q_Charge1 (   �   T   a   Write_Q_Charge1%FILENAME &   9!  �   a   Write_Q_Charge1%CHARGE &   �!  H   a   Write_Q_Charge1%NTIMES $   E"  H   a   Write_Q_Charge1%TMIN "   �"  H   a   Write_Q_Charge1%DT &   �"  H   a   Write_Q_Charge1%NATOMS    #  �       Write_Charge_Density #   �#  T   a   Write_Charge_Density%FILENAME    $  H   a   Write_Charge_Density%NPTS     P$  �   a   Write_Charge_Density%GRIDV     �$  �   a   Write_Charge_Density%CHDEN "   �%  �   a   Write_Charge_Density%WEIGHTV !   <&  H   a   Write_Charge_Density%NATOMS    �&  �       WRITE_WEIGHTS '   '  T   a   WRITE_WEIGHTS%FILENAME &   d'  �   a   WRITE_WEIGHTS%WEIGHTV $   (  �   a   WRITE_WEIGHTS%GRIDV %   �(  H   a   WRITE_WEIGHTS%NATOMS #   )  H   a   WRITE_WEIGHTS%NPTS    L)  �       READ_WEIGHTS &   �)  T   a   READ_WEIGHTS%FILENAME %   !*  �   a   READ_WEIGHTS%WEIGHTV $   �*  H   a   READ_WEIGHTS%NATOMS "   +  H   a   READ_WEIGHTS%NPTS    ]+  �       WRITE_SUMMARY '   E,  T   a   WRITE_SUMMARY%FILENAME #   �,  H   a   WRITE_SUMMARY%NPTS %   �,  H   a   WRITE_SUMMARY%NATOMS %   )-  H   a   WRITE_SUMMARY%VOLUME .   q-  H   a   WRITE_SUMMARY%COMPUTED_VOLUME %   �-  H   a   WRITE_SUMMARY%NTIMES #   .  H   a   WRITE_SUMMARY%TMIN #   I.  H   a   WRITE_SUMMARY%TMAX '   �.  �   a   WRITE_SUMMARY%ATOMNAME (   -/  �   a   WRITE_SUMMARY%RADIUS_BS #   �/  H   a   WRITE_SUMMARY%NORB %   	0  �   a   WRITE_SUMMARY%ORBTAB 