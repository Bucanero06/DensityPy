  ýV  Ë   k820309    é          2021.10.0   LsÖd                                                                                                          
       src/Module_CD_IO.f90 MODULE_CD_IO              SET_CD_IO_VERBOUS LOADGEOMETRY LOADENERGIES WRITEDIPOLEFTFILE LOAD_Q_CHARGE_AND_WRITE2ALL WRITEATOMICCHARGEFT APPENDDIPOLE2FTALLFILE WRITEALLATOMICCHARGEFTTOSINGLEFILE LOAD_DIPOLE LOADATOMICCHARGES LOADFTDIPOLE_ASFUNCITONOF_TIMEDELAY SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM LOADFTOFCHARGEASFUNCOFTIMEDELAY WRITE_BIDIMENTIONALCHARGEFTWW LOAD_XUVATOMICCHARGE            0                                                
                            @                                      
                                                                   
                            @                                      
                     @   @                                        '                   #N    #P    #PRINTPULSES    #WRITE    #WRITEFTA    #GETMAXFIELDSTRENGTH    #FREE    #ADD    #A     #E "   #FTA $   #FTE &   #GETFREQUENCYSPACE (   #PULSE_TRAINWRITE    #PULSE_TRAINWRITEFTA    #PULSE_TRAINGETMAXFIELDSTRENGTH    #FREEPULSETRAIN    #ADDPULSESTRN    #ADDPULSEPARAMS    #AFIELD !   #EFIELD #   #FTAFIELD %   #FTEFIELD '   #PULSE_TRAINGETFREQUENCYSPACE )   #PULSE_TRAIN_PRINTPULSES                  $                                                                        $                                            @              P             #PULSE    p          p @           p @                                       @   @                                      'P              
      #ENVELOPE 	   #T 
   #O    #F    #D    #I    #A    #P    #THETA    #PHI                                                         	                                                                       
               
                                                                      
                                                                      
                                                                       
                                                            (          
                                                            0          
                                                            8          
                                                            @       	   
                                                            H       
   
   4             $                          @                                3             $                          @                     u #PULSE_TRAIN    #PULSE_TRAIN_PRINTPULSES    4             $                          @                                3             $                          @                     u #PULSE_TRAIN    #PULSE_TRAINWRITE    4             $                          @                                3             $                          @                     u #PULSE_TRAIN    #PULSE_TRAINWRITEFTA    4             $                          @                                3             $                          @                     u #PULSE_TRAIN    #PULSE_TRAINGETMAXFIELDSTRENGTH    4             $                          @                                3             $                          @                     u #PULSE_TRAIN    #FREEPULSETRAIN    4             $                          @                                3             $                          @                     u #PULSE_TRAIN    #ADDPULSESTRN    #ADDPULSEPARAMS    4             $                          @                          	       3             $                          @                     u #PULSE_TRAIN    #AFIELD !   4             $                          @            "             
       3             $                          @                     u #PULSE_TRAIN    #EFIELD #   4             $                         @            $                    3             $                         @                     u #PULSE_TRAIN    #FTAFIELD %   4             $                          @            &                    3             $                          @                     u #PULSE_TRAIN    #FTEFIELD '   4             $                          @            (                    3             $                          @                     u #PULSE_TRAIN    #PULSE_TRAINGETFREQUENCYSPACE )   1         À    D                                                      #PULSE_TRAINWRITE *   #         @     @                                     *                    #TRAIN +   #FILENAME ,   #TMIN -   #TMAX .   #DT /             
  @                                       +                  #PULSE_TRAIN              
  @                                      ,                    1           
                                          -     
                
                                          .     
                
                                          /     
      1         À    D                                                      #PULSE_TRAINWRITEFTA 0   #         @     @                                     0                    #SELF 1   #FILENAME 2             
 @                                       1                   #PULSE_TRAIN              
                                         2                    1 1         À    D                                                     #PULSE_TRAINGETMAXFIELDSTRENGTH 3   %         @   @                                   3                    
       #TRAIN 4   #TMIN 5   #TMAX 6   #DT 7             
  @                                       4                  #PULSE_TRAIN              
                                          5     
                
                                          6     
                
                                          7     
      1         À    D                                                      #FREEPULSETRAIN 8   #         @     @                                     8                    #TRAIN 9             
                                         9                   #PULSE_TRAIN    1         À    D                                                      #ADDPULSESTRN :   #         @     @                                     :                    #TRAIN ;   #STRN <             
                                         ;                   #PULSE_TRAIN              
  @                                      <                    1 1         À    D                                                      #ADDPULSEPARAMS =   #         @     @                                     =                    #TRAIN >   #ENVELOPE ?   #T0 @   #OMEGA0 A   #NT B   #CEP C   #I0 D             
                                         >                   #PULSE_TRAIN              
                                           ?                                     
                                          @     
                
                                          A     
                
                                          B     
                
                                          C     
                
                                          D     
      1         À    D                                   !                  #AFIELD E   %         @   @                                   E                    
       #TRAIN F   #TIME G             
 @                                       F                   #PULSE_TRAIN              
  @                                       G     
      1         À    D                                   #                  #EFIELD H   %         @   @                                   H                    
       #TRAIN I   #TIME J             
 @                                       I                   #PULSE_TRAIN              
  @                                       J     
      1         À    D                                   %              	    #FTAFIELD K   %         @   @                                   K                           #TRAIN L   #OMEGA M   #MU N             
 @                                       L                   #PULSE_TRAIN              
  @                                       M     
                
 @                                        N           1         À    D                                   '              
    #FTEFIELD O   %         @   @                                   O                           #TRAIN P   #OMEGA Q   #MU R             
 @                                       P                   #PULSE_TRAIN              
  @                                       Q     
                
 @                                        R           1         À    D                                    )                  #PULSE_TRAINGETFREQUENCYSPACE S   #         @     @                                     S                    #SELF T   #W U             
                                          T                  #PULSE_TRAIN             @                                       U                   
               &                                           1         À    D                                                      #PULSE_TRAIN_PRINTPULSES V   #         @     @                                     V                    #TRAIN W   #UID X             
                                          W                  #PULSE_TRAIN              
                                           X                          À   @                                  Y     '                    #TEAM_PTR Z   #VERSION_NUMBER ]   #CHECK_VALID ^                D                                      Z                          #C_PTR [                                                                                   X#C_PTR [      n                                      0                                               À   @                                 [     '                    #PTR \                 D                                     \                                 D                                      ]                                D                                      ^                  #         @                                            _                    #LOGI `             
                                           `           #         @                                            a                    #NATOMS b   #ATCOORD c   #FILENAME d   #ATOMNAME e             D                                          b                    D                                         c                   
               &                   &                                                     
  @                                      d                    1 ,        D                                         e                                  &                                                   #         @                                            f                    #FILENAME g   #NSTATES h   #EVEC i             
                                         g                    1           D                                          h                    D                                         i                   
               &                                           #         @                                            j                    #FILENAME k   #DIPOLEFTTOTAL l   #OMEGAVEC m   #NOMEGAS n             
                                         k                    1           
                                          l                                 &                   &                                                                                            m                   
               &                                                     
                                           n           #         @                                            o                 	   #FILENAME p   #CHARGE q   #NTIMES r   #TMIN s   #DT t   #NATOMS u   #ISIM v   #ATOM_NAMES w   #UID_ATOMICCHARGEALL x             
                                         p                    1         D @                                       q                   
               &                   &                   &                                                     
                                           r                     
                                          s     
                
                                          t     
                
                                           u                     
  @                                        v           ,          
                                          w                    	             &                                                             
                                           x           #         @                                            y                    #FILENAME z   #ATOMICCHARGEFT {   #OMEGAVEC |   #NOMEGAS }   #NATOMS ~   #ATOMNAME              
                                         z                    1           
                                          {                                 &                   &                                                                                              |                   
               &                                                     
                                           }                     
                                           ~           ,          
                                                                           &                                                   #         @                                                                #FILENAME    #DIPOLEFTTOTAL    #OMEGAVEC    #NOMEGAS    #ISIM    #TRAIN              
                                                             1           
                                                                           &                   &                                                                                                                 
               &                                                     
                                                                
                                                               
  P                                                                         &                                           #PULSE_TRAIN    #         @                                                                #FILENAME    #ATOMICCHARGEFT_NEW    #OMEGAVEC    #NOMEGAS    #NATOMS    #ISIM    #TRAIN    #ATOMNAME              
                                                             1                                                                                      &                   &                   &                                                                                                                 
               &                                                     
                                                                
                                                                
                                                               
  P                                                                         &                                           #PULSE_TRAIN    ,          
                                                                           &                                                   #         @                                                                #FILENAME    #DIPOLE    #NTIMES              
                                                             1         D                                                                           &                   &                                                     
                                                      #         @                                                                #FILENAME    #ATOMICCHARGEEVOLUTION    #NTIMES    #NATOMS    #ISIM    #TMIN    #DT    #UID_ATOMICCHARGEALL              
                                                             1         D                                                            
 
              &                   &                                                     
                                                                
                                                                
  @                                                             
                                               
                
                                               
                
                                                      #         @                                                                #FILENAME    #N_SIMULATIONS    #NOMEGAS     #DIPOLEFTWT ¡   #TVEC ¢             
                                                             1           
                                                                
                                                               D                                         ¡                                  &                   &                   &                                                   D                                         ¢                   
               &                                           #         @                                            £                    #FILENAME ¤   #DIPOLEFTWW ¥   #TAUOMEGAVEC ¦   #OMEGAVEC §   #NTAUOMEGAS ¨   #NOMEGAS ©             
                                         ¤                    1           
                                          ¥                                 &                   &                   &                                                     
                                          ¦                   
              &                                                     
                                          §                   
              &                                                     
                                           ¨                     
                                           ©           #         @                                            ª                    #FILENAME «   #N_SIMULATIONS ¬   #NOMEGAS ­   #NATOMS ®   #TVEC ¯   #CHARGEFTWT_NEW °             
                                         «                    1           
                                           ¬                     
                                           ­                     
                                           ®                   D                                         ¯                   
               &                                                   D                                         °                                  &                   &                   &                   &                                           #         @                                            ±                    #FILENAME ²   #CHARGEFTWW_NEW ³   #TAUOMEGAVEC ´   #OMEGAVEC µ   #NTAUOMEGAS ¶   #NOMEGAS ·   #NATOMS ¸   #ATOMNAME ¹             
                                         ²                    1           
                                          ³                                  &                   &                   &                   &                                                     
                                          ´                   
 !             &                                                     
                                          µ                   
 "             &                                                     
                                           ¶                     
                                           ·                     
                                           ¸           ,          
                                          ¹                    #             &                                                   #         @                                            º                    #FILENAME »   #CHARGE ¼   #NTIMES ½   #NATOMS ¾             
                                         »                    1         D                                         ¼                                  &                   &                   &                                                     
                                           ½                     
                                           ¾                  *      fn#fn "   Ê   n  b   uapp(MODULE_CD_IO     8  H   J  ISO_FORTRAN_ENV       H   J  MODULEPULSES_3D     È  H   J  MODULECONSTANTS $     H   J  MODULEERRORHANDLING ,   X  ó     PULSE_TRAIN+MODULEPULSES_3D .   K  P   a   PULSE_TRAIN%N+MODULEPULSES_3D .     ¯   a   PULSE_TRAIN%P+MODULEPULSES_3D &   J  «      PULSE+MODULEPULSES_3D /   õ  P   a   PULSE%ENVELOPE+MODULEPULSES_3D (   E  P   a   PULSE%T+MODULEPULSES_3D (     P   a   PULSE%O+MODULEPULSES_3D (   å  P   a   PULSE%F+MODULEPULSES_3D (   5  P   a   PULSE%D+MODULEPULSES_3D (     P   a   PULSE%I+MODULEPULSES_3D (   Õ  P   a   PULSE%A+MODULEPULSES_3D (   %	  P   a   PULSE%P+MODULEPULSES_3D ,   u	  P   a   PULSE%THETA+MODULEPULSES_3D *   Å	  P   a   PULSE%PHI+MODULEPULSES_3D 8   
  P   a   PULSE_TRAIN%PRINTPULSES+MODULEPULSES_3D 0   e
  v   `   gen@PRINTPULSES+MODULEPULSES_3D 2   Û
  P   a   PULSE_TRAIN%WRITE+MODULEPULSES_3D *   +  o   `   gen@WRITE+MODULEPULSES_3D 5     P   a   PULSE_TRAIN%WRITEFTA+MODULEPULSES_3D -   ê  r   `   gen@WRITEFTA+MODULEPULSES_3D @   \  P   a   PULSE_TRAIN%GETMAXFIELDSTRENGTH+MODULEPULSES_3D 8   ¬  }   `   gen@GETMAXFIELDSTRENGTH+MODULEPULSES_3D 1   )  P   a   PULSE_TRAIN%FREE+MODULEPULSES_3D )   y  m   `   gen@FREE+MODULEPULSES_3D 0   æ  P   a   PULSE_TRAIN%ADD+MODULEPULSES_3D (   6     `   gen@ADD+MODULEPULSES_3D .   µ  P   a   PULSE_TRAIN%A+MODULEPULSES_3D &     e   `   gen@A+MODULEPULSES_3D .   j  P   a   PULSE_TRAIN%E+MODULEPULSES_3D &   º  e   `   gen@E+MODULEPULSES_3D 0     P   a   PULSE_TRAIN%FTA+MODULEPULSES_3D (   o  g   `   gen@FTA+MODULEPULSES_3D 0   Ö  P   a   PULSE_TRAIN%FTE+MODULEPULSES_3D (   &  g   `   gen@FTE+MODULEPULSES_3D >     P   a   PULSE_TRAIN%GETFREQUENCYSPACE+MODULEPULSES_3D 6   Ý  {   `   gen@GETFREQUENCYSPACE+MODULEPULSES_3D N   X  f   %   PULSE_TRAIN%PULSE_TRAINWRITE+MODULEPULSES_3D=PULSE_TRAINWRITE 1   ¾        PULSE_TRAINWRITE+MODULEPULSES_3D 7   C  a   a   PULSE_TRAINWRITE%TRAIN+MODULEPULSES_3D :   ¤  T   a   PULSE_TRAINWRITE%FILENAME+MODULEPULSES_3D 6   ø  H   a   PULSE_TRAINWRITE%TMIN+MODULEPULSES_3D 6   @  H   a   PULSE_TRAINWRITE%TMAX+MODULEPULSES_3D 4     H   a   PULSE_TRAINWRITE%DT+MODULEPULSES_3D T   Ð  i   %   PULSE_TRAIN%PULSE_TRAINWRITEFTA+MODULEPULSES_3D=PULSE_TRAINWRITEFTA 4   9  h      PULSE_TRAINWRITEFTA+MODULEPULSES_3D 9   ¡  a   a   PULSE_TRAINWRITEFTA%SELF+MODULEPULSES_3D =     T   a   PULSE_TRAINWRITEFTA%FILENAME+MODULEPULSES_3D j   V  t   %   PULSE_TRAIN%PULSE_TRAINGETMAXFIELDSTRENGTH+MODULEPULSES_3D=PULSE_TRAINGETMAXFIELDSTRENGTH ?   Ê        PULSE_TRAINGETMAXFIELDSTRENGTH+MODULEPULSES_3D E   I  a   a   PULSE_TRAINGETMAXFIELDSTRENGTH%TRAIN+MODULEPULSES_3D D   ª  H   a   PULSE_TRAINGETMAXFIELDSTRENGTH%TMIN+MODULEPULSES_3D D   ò  H   a   PULSE_TRAINGETMAXFIELDSTRENGTH%TMAX+MODULEPULSES_3D B   :  H   a   PULSE_TRAINGETMAXFIELDSTRENGTH%DT+MODULEPULSES_3D J     d   %   PULSE_TRAIN%FREEPULSETRAIN+MODULEPULSES_3D=FREEPULSETRAIN /   æ  [      FREEPULSETRAIN+MODULEPULSES_3D 5   A  a   a   FREEPULSETRAIN%TRAIN+MODULEPULSES_3D F   ¢  b   %   PULSE_TRAIN%ADDPULSESTRN+MODULEPULSES_3D=ADDPULSESTRN -     e      ADDPULSESTRN+MODULEPULSES_3D 3   i  a   a   ADDPULSESTRN%TRAIN+MODULEPULSES_3D 2   Ê  T   a   ADDPULSESTRN%STRN+MODULEPULSES_3D J     d   %   PULSE_TRAIN%ADDPULSEPARAMS+MODULEPULSES_3D=ADDPULSEPARAMS /           ADDPULSEPARAMS+MODULEPULSES_3D 5     a   a   ADDPULSEPARAMS%TRAIN+MODULEPULSES_3D 8   y  X   a   ADDPULSEPARAMS%ENVELOPE+MODULEPULSES_3D 2   Ñ  H   a   ADDPULSEPARAMS%T0+MODULEPULSES_3D 6     H   a   ADDPULSEPARAMS%OMEGA0+MODULEPULSES_3D 2   a  H   a   ADDPULSEPARAMS%NT+MODULEPULSES_3D 3   ©  H   a   ADDPULSEPARAMS%CEP+MODULEPULSES_3D 2   ñ  H   a   ADDPULSEPARAMS%I0+MODULEPULSES_3D :   9  \   %   PULSE_TRAIN%AFIELD+MODULEPULSES_3D=AFIELD '     m      AFIELD+MODULEPULSES_3D -     a   a   AFIELD%TRAIN+MODULEPULSES_3D ,   c  H   a   AFIELD%TIME+MODULEPULSES_3D :   «  \   %   PULSE_TRAIN%EFIELD+MODULEPULSES_3D=EFIELD '      m      EFIELD+MODULEPULSES_3D -   t   a   a   EFIELD%TRAIN+MODULEPULSES_3D ,   Õ   H   a   EFIELD%TIME+MODULEPULSES_3D >   !  ^   %   PULSE_TRAIN%FTAFIELD+MODULEPULSES_3D=FTAFIELD )   {!  v      FTAFIELD+MODULEPULSES_3D /   ñ!  a   a   FTAFIELD%TRAIN+MODULEPULSES_3D /   R"  H   a   FTAFIELD%OMEGA+MODULEPULSES_3D ,   "  H   a   FTAFIELD%MU+MODULEPULSES_3D >   â"  ^   %   PULSE_TRAIN%FTEFIELD+MODULEPULSES_3D=FTEFIELD )   @#  v      FTEFIELD+MODULEPULSES_3D /   ¶#  a   a   FTEFIELD%TRAIN+MODULEPULSES_3D /   $  H   a   FTEFIELD%OMEGA+MODULEPULSES_3D ,   _$  H   a   FTEFIELD%MU+MODULEPULSES_3D f   §$  r   %   PULSE_TRAIN%PULSE_TRAINGETFREQUENCYSPACE+MODULEPULSES_3D=PULSE_TRAINGETFREQUENCYSPACE =   %  a      PULSE_TRAINGETFREQUENCYSPACE+MODULEPULSES_3D B   z%  a   a   PULSE_TRAINGETFREQUENCYSPACE%SELF+MODULEPULSES_3D ?   Û%     a   PULSE_TRAINGETFREQUENCYSPACE%W+MODULEPULSES_3D \   o&  m   %   PULSE_TRAIN%PULSE_TRAIN_PRINTPULSES+MODULEPULSES_3D=PULSE_TRAIN_PRINTPULSES 8   Ü&  d      PULSE_TRAIN_PRINTPULSES+MODULEPULSES_3D >   @'  a   a   PULSE_TRAIN_PRINTPULSES%TRAIN+MODULEPULSES_3D <   ¡'  H   a   PULSE_TRAIN_PRINTPULSES%UID+MODULEPULSES_3D *   é'        TEAM_TYPE+ISO_FORTRAN_ENV <   t(    %   TEAM_TYPE%TEAM_PTR+ISO_FORTRAN_ENV=TEAM_PTR $   )  a      C_PTR+ISO_C_BINDING ,   ð)  P   %   C_PTR%PTR+ISO_C_BINDING=PTR H   @*  P   %   TEAM_TYPE%VERSION_NUMBER+ISO_FORTRAN_ENV=VERSION_NUMBER B   *  P   %   TEAM_TYPE%CHECK_VALID+ISO_FORTRAN_ENV=CHECK_VALID "   à*  Z       SET_CD_IO_VERBOUS '   :+  H   a   SET_CD_IO_VERBOUS%LOGI    +         LOADGEOMETRY $   ,  H   a   LOADGEOMETRY%NATOMS %   O,  ¬   a   LOADGEOMETRY%ATCOORD &   û,  T   a   LOADGEOMETRY%FILENAME &   O-     a   LOADGEOMETRY%ATOMNAME    ë-  u       LOADENERGIES &   `.  T   a   LOADENERGIES%FILENAME %   ´.  H   a   LOADENERGIES%NSTATES "   ü.     a   LOADENERGIES%EVEC "   /         WRITEDIPOLEFTFILE +   0  T   a   WRITEDIPOLEFTFILE%FILENAME 0   p0  ¬   a   WRITEDIPOLEFTFILE%DIPOLEFTTOTAL +   1     a   WRITEDIPOLEFTFILE%OMEGAVEC *   °1  H   a   WRITEDIPOLEFTFILE%NOMEGAS ,   ø1  Ç       LOAD_Q_CHARGE_AND_WRITE2ALL 5   ¿2  T   a   LOAD_Q_CHARGE_AND_WRITE2ALL%FILENAME 3   3  Ä   a   LOAD_Q_CHARGE_AND_WRITE2ALL%CHARGE 3   ×3  H   a   LOAD_Q_CHARGE_AND_WRITE2ALL%NTIMES 1   4  H   a   LOAD_Q_CHARGE_AND_WRITE2ALL%TMIN /   g4  H   a   LOAD_Q_CHARGE_AND_WRITE2ALL%DT 3   ¯4  H   a   LOAD_Q_CHARGE_AND_WRITE2ALL%NATOMS 1   ÷4  H   a   LOAD_Q_CHARGE_AND_WRITE2ALL%ISIM 7   ?5     a   LOAD_Q_CHARGE_AND_WRITE2ALL%ATOM_NAMES @   Û5  H   a   LOAD_Q_CHARGE_AND_WRITE2ALL%UID_ATOMICCHARGEALL $   #6  §       WRITEATOMICCHARGEFT -   Ê6  T   a   WRITEATOMICCHARGEFT%FILENAME 3   7  ¬   a   WRITEATOMICCHARGEFT%ATOMICCHARGEFT -   Ê7     a   WRITEATOMICCHARGEFT%OMEGAVEC ,   ^8  H   a   WRITEATOMICCHARGEFT%NOMEGAS +   ¦8  H   a   WRITEATOMICCHARGEFT%NATOMS -   î8     a   WRITEATOMICCHARGEFT%ATOMNAME '   9  ¡       APPENDDIPOLE2FTALLFILE 0   +:  T   a   APPENDDIPOLE2FTALLFILE%FILENAME 5   :  ¬   a   APPENDDIPOLE2FTALLFILE%DIPOLEFTTOTAL 0   +;     a   APPENDDIPOLE2FTALLFILE%OMEGAVEC /   ¿;  H   a   APPENDDIPOLE2FTALLFILE%NOMEGAS ,   <  H   a   APPENDDIPOLE2FTALLFILE%ISIM -   O<  ¥   a   APPENDDIPOLE2FTALLFILE%TRAIN 3   ô<  À       WRITEALLATOMICCHARGEFTTOSINGLEFILE <   ´=  T   a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%FILENAME F   >  Ä   a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%ATOMICCHARGEFT_NEW <   Ì>     a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%OMEGAVEC ;   `?  H   a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%NOMEGAS :   ¨?  H   a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%NATOMS 8   ð?  H   a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%ISIM 9   8@  ¥   a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%TRAIN <   Ý@     a   WRITEALLATOMICCHARGEFTTOSINGLEFILE%ATOMNAME    yA  v       LOAD_DIPOLE %   ïA  T   a   LOAD_DIPOLE%FILENAME #   CB  ¬   a   LOAD_DIPOLE%DIPOLE #   ïB  H   a   LOAD_DIPOLE%NTIMES "   7C  Æ       LOADATOMICCHARGES +   ýC  T   a   LOADATOMICCHARGES%FILENAME 8   QD  ¬   a   LOADATOMICCHARGES%ATOMICCHARGEEVOLUTION )   ýD  H   a   LOADATOMICCHARGES%NTIMES )   EE  H   a   LOADATOMICCHARGES%NATOMS '   E  H   a   LOADATOMICCHARGES%ISIM '   ÕE  H   a   LOADATOMICCHARGES%TMIN %   F  H   a   LOADATOMICCHARGES%DT 6   eF  H   a   LOADATOMICCHARGES%UID_ATOMICCHARGEALL 4   ­F         LOADFTDIPOLE_ASFUNCITONOF_TIMEDELAY =   EG  T   a   LOADFTDIPOLE_ASFUNCITONOF_TIMEDELAY%FILENAME B   G  H   a   LOADFTDIPOLE_ASFUNCITONOF_TIMEDELAY%N_SIMULATIONS <   áG  H   a   LOADFTDIPOLE_ASFUNCITONOF_TIMEDELAY%NOMEGAS ?   )H  Ä   a   LOADFTDIPOLE_ASFUNCITONOF_TIMEDELAY%DIPOLEFTWT 9   íH     a   LOADFTDIPOLE_ASFUNCITONOF_TIMEDELAY%TVEC 1   I  ª       SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM :   +J  T   a   SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM%FILENAME <   J  Ä   a   SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM%DIPOLEFTWW =   CK     a   SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM%TAUOMEGAVEC :   ×K     a   SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM%OMEGAVEC <   kL  H   a   SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM%NTAUOMEGAS 9   ³L  H   a   SAVEBIDIMENTIOAL_DIPOLE_SPECTRUM%NOMEGAS 0   ûL  ¨       LOADFTOFCHARGEASFUNCOFTIMEDELAY 9   £M  T   a   LOADFTOFCHARGEASFUNCOFTIMEDELAY%FILENAME >   ÷M  H   a   LOADFTOFCHARGEASFUNCOFTIMEDELAY%N_SIMULATIONS 8   ?N  H   a   LOADFTOFCHARGEASFUNCOFTIMEDELAY%NOMEGAS 7   N  H   a   LOADFTOFCHARGEASFUNCOFTIMEDELAY%NATOMS 5   ÏN     a   LOADFTOFCHARGEASFUNCOFTIMEDELAY%TVEC ?   cO  Ü   a   LOADFTOFCHARGEASFUNCOFTIMEDELAY%CHARGEFTWT_NEW .   ?P  È       WRITE_BIDIMENTIONALCHARGEFTWW 7   Q  T   a   WRITE_BIDIMENTIONALCHARGEFTWW%FILENAME =   [Q  Ü   a   WRITE_BIDIMENTIONALCHARGEFTWW%CHARGEFTWW_NEW :   7R     a   WRITE_BIDIMENTIONALCHARGEFTWW%TAUOMEGAVEC 7   ËR     a   WRITE_BIDIMENTIONALCHARGEFTWW%OMEGAVEC 9   _S  H   a   WRITE_BIDIMENTIONALCHARGEFTWW%NTAUOMEGAS 6   §S  H   a   WRITE_BIDIMENTIONALCHARGEFTWW%NOMEGAS 5   ïS  H   a   WRITE_BIDIMENTIONALCHARGEFTWW%NATOMS 7   7T     a   WRITE_BIDIMENTIONALCHARGEFTWW%ATOMNAME %   ÓT         LOAD_XUVATOMICCHARGE .   UU  T   a   LOAD_XUVATOMICCHARGE%FILENAME ,   ©U  Ä   a   LOAD_XUVATOMICCHARGE%CHARGE ,   mV  H   a   LOAD_XUVATOMICCHARGE%NTIMES ,   µV  H   a   LOAD_XUVATOMICCHARGE%NATOMS 