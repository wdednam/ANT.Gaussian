!**********************************************************
!*********************  ANT.G-2.5.0  **********************
!**********************************************************
!*                                                        *
!*  Copyright (c) by                                      *
!*                                                        *
!*  Juan Jose Palacios (1)                                *
!*  David Jacob (2)                                       *
!*  Maria Soriano (1)                                     *
!*                                                        *
!* (1) Departamento de Fisica de la Materia Condensada    *
!*     Universidad Autonoma de Madrid                     *      
!*     28049 Madrid (SPAIN)                               *
!* (2) Theory Department                                  *
!*     Max-Planck-Institute for Microstructure Physics    *
!*     Halle, 06120 (GERMANY)                             *
!*                                                        *
!**********************************************************
  MODULE parameters
!**********************************************************
!* Module for calculation / evaluation parameters         *
!**********************************************************
  USE preproc, ONLY: MaxAtm
  IMPLICIT NONE
  SAVE

  ! ****************************
  ! Basic calculation parameters
  ! ****************************

  ! Mixing parameter (0.0 < alpha < 0.1)
  REAL*8 :: alpha = 0.03d0         
  CHARACTER(len=10), PARAMETER ::  alpha_keyw = "ALPHA"     

  ! Pulay parameter (1 < NPulay < 1000)
  INTEGER :: NPulay = 4
  CHARACTER(len=10), PARAMETER ::  NPulay_keyw = "NPULAY"     

  ! Error in total charge
  REAL*8 :: ChargeAcc = 1.0d-2
  REAL*8 :: ChargeA = 1.0d-5    
  CHARACTER(len=10), PARAMETER ::  ChargeAcc_keyw = "CHARGEACC"       

  ! Precision for Fermi level (eV)
  REAL*8 :: FermiAcc = 1.0d-2
  REAL*8 :: FermiA = 1.0d-5
  CHARACTER(len=10), PARAMETER :: FermiAcc_keyw = "FERMIACC"  

  ! Error in numerical integration of density matrix
  REAL*8 :: PAcc
  REAL*8 :: PA = 1.0d-7
  CHARACTER(len=10), PARAMETER ::  PAcc_keyw = "PACC"       

  ! Work at full accuracy from the begining
  LOGICAL :: FullAcc = .FALSE. ; CHARACTER(len=10), PARAMETER :: FullAcc_keyw = "FULLACC"

  ! Maximum number of steps when searching for the Fermi level
  INTEGER :: Max = 10
  CHARACTER(len=10), PARAMETER :: Max_keyw = "MAX"  

  ! Self-energy accuracy
  REAL*8 :: SelfAcc = 1.0d-6     
  CHARACTER(len=10), PARAMETER :: SelfAcc_keyw = "SELFACC"  

  ! Parameters for identifying Bethe lattice directions
  REAL*8 :: Small = 0.1
  REAL*8 :: SmallD = 0.1
  CHARACTER(len=10), PARAMETER :: Small_keyw = "SMALL"  
  CHARACTER(len=10), PARAMETER :: SmallD_keyw = "SMALLD"  

  ! Imaginary part
  REAL*8 :: eta = 1.0d-10         
  CHARACTER(len=10), PARAMETER :: eta_keyw = "ETA"       

  ! Bethe lattice glue parameter
  REAL*8 :: Glue = 1.0     
  CHARACTER(len=10), PARAMETER :: Glue_keyw = "GLUE"       

  ! Bias voltage (V)
  REAL*8 :: BiasVoltage = 0.00d0
  CHARACTER(len=10), PARAMETER :: Bias_keyw = "BIASVOLT"    

  ! Excess charge
  REAL*8 :: QExcess = 0.00d0     
  CHARACTER(len=10), PARAMETER :: QExcess_keyw = "QEXCESS" 

  ! Convergence before switching on leads
  REAL*8 :: SL = 1.0d-2           
  CHARACTER(len=10), PARAMETER :: SL_keyw = "SL"        

  ! Starting value for Fermi level search (eV)
  REAL*8 :: FermiStart = -5.00d0   
  CHARACTER(len=10), PARAMETER :: FermiStart_keyw = "FERMISTART"

  ! Type of left and right electrode: "BETHE" or "1DLEAD"
  CHARACTER(LEN=10), DIMENSION(2) :: ElType = (/"BETHE","BETHE"/)
  CHARACTER(len=10), DIMENSION(2), PARAMETER :: ElType_keyw = (/"TYPE1","TYPE2"/)

  ! Bethe lattice parameters
  CHARACTER(LEN=20), DIMENSION(2) :: BLPar = (/"Papacon","Papacon"/)
  CHARACTER(len=10), DIMENSION(2), PARAMETER :: BLPar_keyw = (/"BLPAR1","BLPAR2"/)

  ! Keyword to use supplementary density matrix
  LOGICAL :: PFix = .FALSE. 
  CHARACTER(len=10), PARAMETER :: PFix_keyw = "PFIX"
  CHARACTER(len=50) :: densitymatrixx 

  ! Keyword to keep atoms unchanged when using PFIX
  INTEGER :: NFix = 0
  INTEGER, DIMENSION( MaxAtm ) :: IFix = 0
  CHARACTER(len=10), PARAMETER :: NFix_keyw = "NFIX"
  
  ! Keyword to determine the type of overlap in the Bethe lattice basis set
  REAL*8 :: Overlap = 0.0d0           
  CHARACTER(len=10), PARAMETER :: Overlap_keyw = "OVERLAP"

  ! Number of atoms to be connected to Bethe lattice
  CHARACTER(LEN=10), DIMENSION(2), PARAMETER :: NEmbed_keyw = (/"NEMBED1","NEMBED2"/)
  INTEGER, DIMENSION(2) :: NEmbed = (/0,0/)

  ! Manual selection of number of atoms in electrodes for single element systems
  CHARACTER(LEN=10), DIMENSION(2), PARAMETER :: NAtomEl_keyw = (/"NATOMEL1","NATOMEL2"/)
  INTEGER, DIMENSION(2) :: NAtomEL = (/0,0/)

  ! Switch on Fock matrix mixing (instead of density matrix) for damping
  LOGICAL :: FMixing = .FALSE. ; CHARACTER(len=10), PARAMETER :: FMixing_keyw = "FMIXING"
   
  ! Whether to compute the density matrix by integration along imaginary axis instead of 
  ! complex contour (default)
  LOGICAL :: DMImag = .false.; CHARACTER(len=10), PARAMETER :: DMImag_keyw = "DMIMAG"

  ! *****************
  ! Correlated Blocks
  ! *****************
  CHARACTER(LEN=10), PARAMETER :: CorrBl_keyw = "CORRBLOCKS"
  INTEGER :: NCorrBl = 0
  INTEGER, dimension(MaxAtm) :: CorrBeg=0, CorrEnd=0
  REAL*8, dimension(MaxAtm) :: UCoul=0, JHund=0

  ! Whether to perform DFT+U calculation
  LOGICAL :: DFTU = .false.; CHARACTER(len=10), PARAMETER :: DFTU_keyw = "DFT+U"
  
  ! Whether to fix density in correlated subspace (for charge self-consistency)
  Logical :: FixCorrDens = .false.; CHARACTER(len=10), PARAMETER :: FixCorrDens_keyw = "FIXCORRDE"

  ! Whether to calculate hybridization function
  Logical :: HybFunc = .false.; CHARACTER(len=10), PARAMETER :: HybFunc_keyw = "HYBFUNC"
  
  ! Non-orthogonal Bethe lattice: Deorthogonalize Bethe Lattice Self-energy
  LOGICAL :: Portho = .FALSE. ; CHARACTER(len=10), PARAMETER :: Portho_keyw = "PORTHO"

  ! Whether to diagonalize the correlated blocks
  LOGICAL :: DiagCorrBl = .FALSE. ; CHARACTER(len=10), PARAMETER :: DiagCorrBl_keyw = "DIAGCORRBL"

  ! *************************
  ! Spin transport parameters
  ! *************************

  ! Magnetic boundary conditions: U=up, D=down
  LOGICAL :: UD    = .FALSE. ; CHARACTER(len=10), PARAMETER :: UD_keyw = "UD"   
  LOGICAL :: DU    = .FALSE. ; CHARACTER(len=10), PARAMETER :: DU_keyw = "DU"   
  LOGICAL :: DD    = .FALSE. ; CHARACTER(len=10), PARAMETER :: DD_keyw = "DD"   

  ! Magnetization reversal of initial guess for all atoms starting from specified value
  INTEGER :: MRStart = 0 
  CHARACTER(len=10), PARAMETER :: MRStart_keyw = "MRSTART"
  
  ! Erase spin density for all atoms 
  LOGICAL :: SpinDel = .false. 
  CHARACTER(len=10), PARAMETER :: SpinDel_keyw = "SPINDEL"

  ! Atom spin manipulation of initial guess
  INTEGER :: NSpinEdit = 0
  INTEGER, DIMENSION( MaxAtm ) :: SpinEdit = 1
  CHARACTER(LEN=10), PARAMETER :: SpinEdit_keyw = "SPINEDIT"

  ! Number of steps with spinlock
  INTEGER :: NSpinLock = 0
  CHARACTER(len=10), PARAMETER :: NSpinLock_keyw = "NSPINLOCK"

  ! Number of alpha and beta electrons    
  INTEGER :: Nalpha = -1
  INTEGER :: Nbeta = -1
  CHARACTER(len=10), PARAMETER :: Nalpha_keyw = "NALPHA"
  CHARACTER(len=10), PARAMETER :: Nbeta_keyw = "NBETA"

  ! Convergence before switching off spin lock, (  SwOffSpL < SL  )
  REAL*8 :: SwOffSpL = 1.0d-3            
  CHARACTER(len=10), PARAMETER ::  SwOffSpL_keyw = "SWOFFSPL"        

  ! SOC
  LOGICAL :: soc = .FALSE. 
  CHARACTER(len=10), PARAMETER :: SOC_keyw = "SOC"  
  ! REAL*8 :: soc = 0.0d0           
  ! CHARACTER(len=10), PARAMETER :: SOC_keyw = "SOC"
  
  ! SOC multiplicative factor due to lack of nodal structure in basis set
  INTEGER :: socfac = 1
  CHARACTER(LEN=10), PARAMETER :: SOCFAC_keyw = "SOCFAC"  

  ! SOC_COEFF_P
  REAL*8 :: soc_cff_p = 0.0d0                           
  CHARACTER(len=10), PARAMETER :: SOC_CFF_P_keyw = "SOC_CFF_P"
  
  ! SOC_COEFF_D                                             
  REAL*8 :: soc_cff_d = 0.0d0                             
  CHARACTER(len=10), PARAMETER :: SOC_CFF_D_keyw = "SOC_CFF_D"    

  ! SOC_COEFF_F                                             
  REAL*8 :: soc_cff_f = 0.0d0                             
  CHARACTER(len=10), PARAMETER :: SOC_CFF_F_keyw = "SOC_CFF_F"   
    
  ! THETA
  REAL*8 :: theta = 0.0d0                           
  CHARACTER(len=10), PARAMETER :: THETA_keyw = "THETA"
  
  ! PHI                                             
  REAL*8 :: phi = 0.0d0                             
  CHARACTER(len=10), PARAMETER :: PHI_keyw = "PHI"

  ! Atom SOC definition
  INTEGER :: NSOCEdit = 0
  INTEGER, DIMENSION( MaxAtm ) :: SOCEdit = 1
  CHARACTER(LEN=10), PARAMETER :: SOCEdit_keyw = "SOCEDIT"
  
  ! Atom spin orientation manipulation of initial guess
  INTEGER :: NSpinRot = 0
  REAL*8, DIMENSION( MaxAtm) :: SpinRotTheta = 0.0d0, SpinRotPhi = 0.0d0
  CHARACTER(LEN=10), PARAMETER :: SpinRot_keyw = "SPINROT"       
  

  ! *********************
  ! Output parameters
  ! *********************

  ! Print out Hamiltonian
  LOGICAL :: Hamilton = .FALSE. 
  CHARACTER(len=10), PARAMETER :: Hamilton_keyw = "HAMILTON"
  
  ! Print out Mulliken population analysis
  LOGICAL :: Mulliken = .FALSE. 
  CHARACTER(len=10), PARAMETER :: Mulliken_keyw = "MULLIKEN"
  
  
  ! Energy at which to evaluate the Projected DOS
  REAL*8:: DOSEnergy= 0.0d0       
  CHARACTER(len=10), PARAMETER :: DOSEnergy_keyw = "DOSENERGY"   

  ! Energy step 
  REAL*8:: EStep = 0.01d0       
  CHARACTER(len=10), PARAMETER :: EStep_keyw = "ESTEP"   
  
  ! Energy window from EW1 to EW2
  REAL*8:: EW1 = -3.00d0     
  CHARACTER(len=10), PARAMETER :: EW1_keyw = "EW1"  
  REAL*8:: EW2 =  3.00d0     
  CHARACTER(len=10), PARAMETER :: EW2_keyw = "EW2"  

  ! calculate LDOS for atom numbers from LDOS_BEG up to LDOS_END
  INTEGER :: LDOS_Beg = 1      
  CHARACTER(len=10), PARAMETER :: LDOS_Beg_keyw = "LDOS_BEG" 
  INTEGER :: LDOS_End = 0      
  CHARACTER(len=10), PARAMETER :: LDOS_End_keyw = "LDOS_END" 

  ! Use hermitian expression for transmission matrix.
  LOGICAL :: HTransm = .FALSE. 
  CHARACTER(len=10), PARAMETER :: HTransm_keyw = "HTRANSM"
  
  ! Number of eigen channels to print out ( 0 = No diagonalization of T )
  INTEGER :: NChannels = 0     
  CHARACTER(len=10), PARAMETER :: NChannels_keyw = "NCHANNELS"
  
  ! Print out Bulk DOS
  LOGICAL :: LeadDOS = .FALSE. 
  CHARACTER(len=10), PARAMETER :: LeadDOS_keyw = "LEADDOS"  

  ! Eigen channel analysis with reduced transmission matrix
  INTEGER :: RedTransmB = 1      
  CHARACTER(len=10), PARAMETER :: RedTransmB_keyw = "RTM_BEG" 
  INTEGER :: RedTransmE = 0  
  CHARACTER(len=10), PARAMETER :: RedTransmE_keyw = "RTM_END" 
  
  ! Whether to generate input for ANT.1D program
  LOGICAL :: ANT1DInp = .FALSE.
  CHARACTER(len=10), PARAMETER :: ANT1DInp_keyw = "ANT1D"

  ! Number of atom whose Hamiltonian to print out
  INTEGER :: PrtHatom = 1     
  CHARACTER(len=10), PARAMETER :: PrtHatom_keyw = "PRTHATOM"

  PRIVATE :: read_line

CONTAINS
  ! **********************
  ! Read single input line
  ! **********************
  SUBROUTINE read_line( logfile, inpfile, ios )
    use util, only: upcase
    IMPLICIT NONE
    
    INTEGER, INTENT(in)  :: logfile, inpfile
    INTEGER, INTENT(out) :: ios
    
    CHARACTER          :: eqsign, comment
    CHARACTER(len=10)  :: keyword
    CHARACTER(LEN=100) :: strval
    REAL*8             :: rval,rvall
    INTEGER            :: ival, index, i, imax
    
    ! Jump comment lines
    READ (unit=inpfile,fmt=*,iostat=ios), comment
    IF( ios /= 0 .OR. trim(comment) == '!' ) RETURN
    BACKSPACE inpfile

    READ (unit=inpfile,fmt=*,iostat=ios), keyword
    IF( ios /= 0 ) RETURN

    call upcase(keyword)

    ! Evaluating keys 
    SELECT CASE ( keyword )
    CASE ( alpha_keyw     ,&
         & ChargeAcc_keyw ,&
         & PAcc_keyw ,&
         & FermiAcc_keyw  ,&
         & SelfAcc_keyw   ,&
         & Small_keyw   ,&
         & SmallD_keyw   ,&
         & Bias_keyw      ,&
         & QExcess_keyw   ,&
         & eta_keyw       ,&
         & glue_keyw      ,&
         & SL_keyw        ,&
         & SwOffSpL_keyw  ,&
         & FermiStart_keyw,&
         & EStep_keyw     ,&
         & DOSEnergy_keyw     ,&
         & Overlap_keyw   ,&
!         & SOC_keyw   ,&
!         & SOCFAC_keyw  ,&  
         & SOC_CFF_P_keyw   ,&
         & SOC_CFF_D_keyw   ,&
         & SOC_CFF_F_keyw   ,&         
         & THETA_keyw   ,&
         & PHI_keyw   ,&
         & EW1_keyw       ,&
         & EW2_keyw   )
       !
       ! 1. looking for real variables
       !
       BACKSPACE inpfile
       READ (unit=inpfile,fmt=*,iostat=ios), keyword, eqsign, rval      
       IF( ios /= 0 .OR. eqsign /= '=' ) RETURN 

       call upcase(keyword)
       
       SELECT CASE ( keyword )
       CASE ( alpha_keyw )
          alpha = rval
       CASE ( ChargeAcc_keyw )
          ChargeA = rval
       CASE ( FermiAcc_keyw )
          FermiA = rval
       CASE ( PAcc_keyw )
          PA = rval
       CASE ( SelfAcc_keyw )
          SelfAcc = rval
       CASE ( Small_keyw )
          Small = rval
       CASE ( SmallD_keyw )
          SmallD = rval
       CASE ( Bias_keyw )
          BiasVoltage = rval
       CASE ( QExcess_keyw )
          QExcess = rval 
       CASE ( eta_keyw )
          eta = rval
       CASE ( glue_keyw )
          glue = rval
       CASE ( SL_keyw )
          SL = rval
       CASE ( SwOffSpL_keyw )
          SwOffSpL = rval
       CASE ( FermiStart_keyw ) 
          FermiStart = rval
       CASE ( DOSEnergy_keyw ) 
          DOSEnergy = rval
       CASE ( EStep_keyw ) 
          EStep = rval
       CASE ( Overlap_keyw ) 
          Overlap = rval
       CASE ( SOC_CFF_P_keyw ) 
          soc_cff_p = rval
       CASE ( SOC_CFF_D_keyw ) 
          soc_cff_d = rval   
       CASE ( SOC_CFF_F_keyw ) 
          soc_cff_f = rval                     
       CASE ( THETA_keyw )
          theta = rval    
       CASE ( PHI_keyw )   
          phi = rval     
       CASE ( EW1_keyw ) 
          EW1 = rval
       CASE ( EW2_keyw ) 
          EW2 = rval
       END SELECT
       
    CASE ( LDOS_Beg_keyw, LDOS_End_keyw, NChannels_keyw, RedTransmB_keyw, RedTransmE_keyw, &
         MRStart_keyw, NSpinLock_keyw, NEmbed_keyw(1), NEmbed_keyw(2), NAtomEl_keyw(1), NAtomEl_keyw(2), &
         NPulay_keyw, Nalpha_keyw, Nbeta_keyw, Max_keyw, SOCFAC_keyw, PrtHatom_keyw  )
       !
       ! 2. looking for integer variables
       !
       BACKSPACE inpfile
       READ (unit=inpfile,fmt=*,iostat=ios), keyword, eqsign, ival      
       IF( ios /= 0 .OR. eqsign /= '=' ) RETURN 
       
       call upcase(keyword)

       SELECT CASE ( keyword )
       CASE ( LDOS_Beg_keyw )
          LDOS_Beg = ival
       CASE ( LDOS_End_keyw )
          LDOS_End = ival
       CASE( NChannels_keyw )
          NChannels = ival
       CASE( RedTransmB_keyw )
          RedTransmB = ival
       CASE( RedTransmE_keyw )
          RedTransmE = ival
       CASE( MRStart_keyw )
           MRStart = ival
       CASE( NSpinLock_keyw )
          NSpinLock = ival
       CASE( Nalpha_keyw )
          Nalpha = ival
       CASE( Nbeta_keyw )
          Nbeta = ival
       CASE( NEmbed_keyw(1) )
          NEmbed(1) = ival
       CASE( NEmbed_keyw(2) )
          NEmbed(2) = ival
       CASE( NAtomEl_keyw(1) )
          NAtomEl(1) = ival
       CASE( NAtomEl_keyw(2) )
          NAtomEl(2) = ival
       CASE( NPulay_keyw )
          NPulay = ival
       CASE( Max_keyw )
          Max = ival
       CASE( SOCFAC_keyw )
          socfac = ival     
       CASE( PrtHatom_keyw )
          PrtHatom = ival     

       END SELECT
       
    CASE ( ElType_keyw(1), ElType_keyw(2), BLPar_keyw(1), BLPar_keyw(2) )
       !
       ! 3. looking for string variables
       !
       BACKSPACE inpfile
       READ (unit=inpfile,fmt=*,iostat=ios), keyword, eqsign, strval
       IF( ios /= 0 .OR. eqsign /= '=' ) RETURN 

       call upcase(keyword)

       SELECT CASE ( keyword )
       CASE( ElType_keyw(1) ) 
          ElType(1) = strval
       CASE( ElType_keyw(2) )
          ElType(2) = strval
       CASE( BLPar_keyw(1) ) 
          BLPar(1) = strval
       CASE( BLPar_keyw(2) )
          BLPar(2) = strval
       END SELECT
       !
       ! 4. Logical switches
       !
    CASE ( DU_keyw )
       DU = .TRUE.
    CASE ( UD_keyw )
       UD = .TRUE.
    CASE ( DD_keyw )
       DD = .TRUE.
    CASE ( LeadDOS_keyw )
       LeadDOS = .TRUE.
    CASE ( HTransm_keyw )
       HTransm = .TRUE.
    CASE ( Hamilton_keyw )
       Hamilton = .TRUE.
    CASE ( Mulliken_keyw )
       Mulliken = .TRUE.
    CASE ( FullAcc_keyw )
       FullAcc = .TRUE.
    CASE ( ANT1DInp_keyw )
       ANT1DInp = .TRUE.
    CASE ( DFTU_keyw )
       DFTU = .true.
    CASE ( HybFunc_keyw )
       HybFunc = .true.
    CASE ( FixCorrDens_keyw )
       FixCorrDens = .true.
    CASE ( Portho_keyw )
       Portho = .true.
    CASE ( DiagCorrBl_keyw )
       DiagCorrBl = .true.
    CASE ( SOC_keyw )
       soc = .true.        
    CASE ( SpinDel_keyw )
       SpinDel = .true.
    CASE ( FMixing_keyw )
       FMixing = .true.
    CASE ( DMImag_keyw )
       DMImag = .true.
       !
       ! 5. Integer arrays
       !
    CASE ( SPINEDIT_keyw )
       READ (unit=inpfile,fmt=*,iostat=ios), NSpinEdit
       IF( ios /= 0 ) RETURN 
       DO i=1,NSpinEdit
          READ (unit=inpfile,fmt=*,iostat=ios), index, ival
          IF( ios /= 0 ) RETURN 
          IF( abs(ival) > 1 )THEN
             WRITE( unit=logfile, fmt=* ) "ERROR - Illegal spin value in SPINEDIT field"
             WRITE( unit=logfile, fmt=* ) "Allowed values: -1, 0, 1"
             ios = 1
             RETURN
          END IF
          IF( index > MaxAtm .OR. index < 1 )THEN
             WRITE( unit=logfile, fmt=* ) "ERROR - Illegal atom number in SPINEDIT field"
             WRITE( unit=logfile, fmt=* ) "Allowed values: 1 ... 10000"
             ios = 1
             RETURN
          END IF
          SpinEdit( index ) = ival
       END DO
       
     CASE ( SOCEDIT_keyw )
       READ (unit=inpfile,fmt=*,iostat=ios), NSOCEdit
       IF( ios /= 0 ) RETURN 
       DO i=1,NSOCEdit
          READ (unit=inpfile,fmt=*,iostat=ios), index, ival
          IF( ios /= 0 ) RETURN 
          IF( abs(ival) > 1 )THEN
             WRITE( unit=logfile, fmt=* ) "ERROR - Illegal soc value in SOCEDIT field"
             WRITE( unit=logfile, fmt=* ) "Allowed values: 0, 1"
             ios = 1
             RETURN
          END IF
          IF( index > MaxAtm .OR. index < 1 )THEN
             WRITE( unit=logfile, fmt=* ) "ERROR - Illegal atom number in SOCEDIT field"
             WRITE( unit=logfile, fmt=* ) "Allowed values: 1 ... 10000"
             ios = 1
             RETURN
          END IF
          SOCEdit( index ) = ival
       END DO

    CASE ( CorrBl_keyw )
       READ (unit=inpfile,fmt=*,iostat=ios), NCorrBl
       print *, "Correlated Blocks: NCorrBlocks = ", NCorrBl
       IF( ios /= 0 ) RETURN 
       DO i=1,NCorrBl
          READ (unit=inpfile,fmt=*,iostat=ios), CorrBeg(i), CorrEnd(i), UCoul(i), JHund(i)
          print *, CorrBeg(i), CorrEnd(i)
          print *, UCoul(i), JHund(i)
          IF( ios /= 0 ) RETURN 
          IF( CorrBeg(i) < 1 .OR. CorrEnd(i) < 1 )THEN
             WRITE( unit=logfile, fmt=* ) "ERROR - Negative number for begin or end of correlated block in CORRBLOCK Field"
             ios = 1
             RETURN
          END IF
       END DO

    CASE ( PFix_keyw )
       PFix = .TRUE.
       READ (unit=inpfile,fmt=*,iostat=ios), densitymatrixx
       if (ios /= 0) RETURN

    CASE ( NFix_keyw )
       BACKSPACE inpfile
       READ (unit=inpfile,fmt=*,iostat=ios), keyword, eqsign, NFix      
       call upcase(keyword)
       IF( ios /= 0 .OR. eqsign /= '=' ) RETURN 
       READ (unit=inpfile,fmt=*,iostat=ios), (IFix(i),i=1,NFix)
       if (ios /= 0) RETURN
       DO i=1,NFix
          IF( IFix(i) > MaxAtm .OR. IFix(i) < 1 )THEN
             WRITE( unit=logfile, fmt=* ) "ERROR - Illegal atom number in NFix field"
             WRITE( unit=logfile, fmt=* ) "Allowed values: 1 ... MaxAtm"
             ios = 1
             RETURN
          END IF
       END DO      
       !                  
       ! 5. Real arrays
       !                  
    CASE ( SPINROT_keyw )
           READ (unit=inpfile,fmt=*,iostat=ios), NSpinRot
           IF( ios /= 0 ) RETURN 
           DO i=1,NSpinRot
              READ (unit=inpfile,fmt=*,iostat=ios), index, rval, rvall
              IF( ios /= 0 ) RETURN 
              IF( abs(rval) > 180.0d0 )THEN
                 WRITE( unit=logfile, fmt=* ) "ERROR - Illegal theta value in SPINROT field"
                 WRITE( unit=logfile, fmt=* ) "Allowed values: -180\BA ... 180\BA"
                 ios = 1
                 RETURN
              END IF
              IF( abs(rvall) > 360.0d0 )THEN                                               
                 WRITE( unit=logfile, fmt=* ) "ERROR - Illegal phi value in SPINROT field"
                 WRITE( unit=logfile, fmt=* ) "Allowed values: -360\BA ... 360\BA"                 
                 ios = 1                                                                  
                 RETURN                                                                   
              END IF                                                                      
              IF( index > MaxAtm .OR. index < 1 )THEN
                 WRITE( unit=logfile, fmt=* ) "ERROR - Illegal atom number in SPINROT field"
                 WRITE( unit=logfile, fmt=* ) "Allowed values: 1 ... 10000"
                 ios = 1
                 RETURN
              END IF
              SpinRotTheta( index ) = rval
              SpinRotPhi( index ) = rvall
           END DO              

    CASE default
       WRITE( unit=logfile, fmt=* ) "ERROR - Undefined keyword: ", keyword
       ios = 1 ! Abort reading 
    END SELECT    

    RETURN
  END SUBROUTINE read_line
  
  ! ***************
  ! Read input file
  ! ***************
  ! Rem: Does not give the correct line number because 
  !      empty lines are simply overread in subroutine readline
  !      What could we do?
  !
  SUBROUTINE read_parameters( logfile, inpfile, ios, xxx )
    IMPLICIT NONE
    
    INTEGER, INTENT(in)  :: logfile, inpfile
    INTEGER, INTENT(out) :: ios
    INTEGER :: nline = 0
    CHARACTER(LEN=50) :: xxx
    
    WRITE( unit=logfile, fmt=* ) 'Scanning ',trim(xxx)//'.ant',' file for parameter specifications...'
    
    DO 
       CALL read_line( logfile, inpfile, ios )
       nline = nline+1
       !print *,nline
       !IF( ios == 1 ) WRITE ( unit=logfile, fmt=* ) "  in line No. " , nline
       IF( ios /= 0 ) EXIT
    END DO
    
    IF( ios > 0  ) THEN
       WRITE ( unit=logfile, fmt=* ) "Error reading ANT input file"
       STOP
    END IF
    IF( ios < 0 ) WRITE ( unit=logfile, fmt=* ) "Done."    
  END SUBROUTINE read_parameters
  

  ! **************************************
  ! write actual parameters to output file
  ! **************************************
  SUBROUTINE write_parameters( logfile )
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: logfile

    INTEGER :: i
    
    WRITE(unit=logfile,fmt=*) "****************************"
    WRITE(unit=logfile,fmt=*) "Basic calculation parameters"
    WRITE(unit=logfile,fmt=*) "****************************"
    WRITE(unit=logfile,fmt=*) alpha_keyw, " = ", alpha
    WRITE(unit=logfile,fmt=*) NPulay_keyw, " = ", NPulay
    WRITE(unit=logfile,fmt=*) ChargeAcc_keyw, " = ", ChargeA, " %"
    WRITE(unit=logfile,fmt=*) FermiAcc_keyw, " = ", FermiA, " %"
    WRITE(unit=logfile,fmt=*) PAcc_keyw, " = ", PA, " %"
    WRITE(unit=logfile,fmt=*) FullAcc_keyw, " = ", FullAcc
    WRITE(unit=logfile,fmt=*) Max_keyw, " = ", Max
    WRITE(unit=logfile,fmt=*) SelfAcc_keyw, " = ", SelfAcc
    WRITE(unit=logfile,fmt=*) Small_keyw, " = ", Small
    WRITE(unit=logfile,fmt=*) SmallD_keyw, " = ", SmallD, " A"
    WRITE(unit=logfile,fmt=*) Bias_keyw, " = ", BiasVoltage, " V"
    WRITE(unit=logfile,fmt=*) QExcess_keyw, " = ", QExcess, "electrons"
    WRITE(unit=logfile,fmt=*) eta_keyw, " = ", eta
    WRITE(unit=logfile,fmt=*) glue_keyw, " = ", glue
    WRITE(unit=logfile,fmt=*) FermiStart_keyw, " = ", FermiStart, " eV"
    WRITE(unit=logfile,fmt=*) SOC_keyw, " = ", soc
    WRITE(unit=logfile,fmt=*) SOCFAC_keyw, " = ", socfac    
    WRITE(unit=logfile,fmt=*) SOC_CFF_P_keyw, " = ", soc_cff_p, " eV"
    WRITE(unit=logfile,fmt=*) SOC_CFF_D_keyw, " = ", soc_cff_d, " eV"    
    WRITE(unit=logfile,fmt=*) SOC_CFF_F_keyw, " = ", soc_cff_f, " eV"        
    WRITE(unit=logfile,fmt=*) THETA_keyw, " = ", theta, " degrees"
    WRITE(unit=logfile,fmt=*) PHI_keyw, " = ", phi, " degrees"   
    WRITE(unit=logfile,fmt=*) SL_keyw, " = ", SL
    WRITE(unit=logfile,fmt=*) DMImag_keyw, " = ", DMImag
    WRITE(unit=logfile,fmt=*) FMixing_keyw, " = ", FMixing
    WRITE(unit=logfile,fmt=*) "************************"
    WRITE(unit=logfile,fmt=*) "Bethe lattice parameters"
    WRITE(unit=logfile,fmt=*) "************************"
    WRITE(unit=logfile,fmt=*) ElType_keyw(1), " = ", ElType(1)
    WRITE(unit=logfile,fmt=*) ElType_keyw(2), " = ", ElType(2)
    WRITE(unit=logfile,fmt=*) BLPar_keyw(1), " = ", BLPar(1)
    WRITE(unit=logfile,fmt=*) BLPar_keyw(2), " = ", BLPar(2)
    if (NEmbed(1) == 0) then
    WRITE(unit=logfile,fmt=*) NEmbed_keyw(1), " = ? "
    else
    WRITE(unit=logfile,fmt=*) NEmbed_keyw(1), " = ", NEmbed(1)
    end if
    if (NEmbed(2) == 0) then
    WRITE(unit=logfile,fmt=*) NEmbed_keyw(2), " = ? "
    else
    WRITE(unit=logfile,fmt=*) NEmbed_keyw(2), " = ", NEmbed(2)
    end if
    if (NAtomEl(1) == 0) then
    WRITE(unit=logfile,fmt=*) NAtomEl_keyw(1), " = ? "
    else
    WRITE(unit=logfile,fmt=*) NAtomEl_keyw(1), " = ", NAtomEl(1)
    end if
    if (NAtomEl(2) == 0) then
    WRITE(unit=logfile,fmt=*) NAtomEl_keyw(2), " = ? "
    else
    WRITE(unit=logfile,fmt=*) NAtomEl_keyw(2), " = ", NAtomEl(2)
    end if
    if (Overlap < 0.0) then
    WRITE(unit=logfile,fmt=*) Overlap_keyw, " = System Overlap"
    else
    WRITE(unit=logfile,fmt=*) Overlap_keyw, " = ", Overlap
    end if
    if (PFix) then
       WRITE(unit=logfile,fmt=*) "****************************"
       WRITE(unit=logfile,fmt=*) "Supplementary density matrix"
       WRITE(unit=logfile,fmt=*) "****************************"
       WRITE(unit=logfile,fmt=*) "Supplementary Density Matrix  ", densitymatrixx
       WRITE(unit=logfile,fmt=*) NFix_keyw, " = ", (IFix(i), i=1,NFix)
    end if
    WRITE(unit=logfile,fmt=*) "*************************"
    WRITE(unit=logfile,fmt=*) "Spin transport parameters"
    WRITE(unit=logfile,fmt=*) "*************************"
    WRITE(unit=logfile,fmt=*) SwOffSpL_keyw, " = ", SwOffSpL
    WRITE(unit=logfile,fmt=*) NSpinLock_keyw, " = ", NSpinLock
    WRITE(unit=logfile,fmt=*) Nalpha_keyw, " = ", Nalpha
    WRITE(unit=logfile,fmt=*) Nbeta_keyw, " = ", Nbeta
    WRITE(unit=logfile,fmt=*) DU_keyw, " = ", DU
    WRITE(unit=logfile,fmt=*) UD_keyw, " = ", UD
    WRITE(unit=logfile,fmt=*) DD_keyw, " = ", DD
    WRITE(unit=logfile,fmt=*) MRStart_keyw, " = ", MRStart
    WRITE(unit=logfile,fmt=*) SpinDel_keyw, " = ", SpinDel
    WRITE(unit=logfile,fmt=*) SpinEdit_keyw, " = ", NSpinEdit
    DO i=1,MaxAtm
       IF( SpinEdit(i) .NE. 1 ) WRITE(unit=logfile,fmt=*) i, SpinEdit(i)
    END DO
    WRITE(unit=logfile,fmt=*) SpinRot_keyw, " = ", NSpinRot
    DO i=1,NSpinRot
       WRITE(unit=logfile,fmt='(I4,F11.4,F11.4)') i, SpinRotTheta(i), SpinRotPhi(i)
    END DO                                                              
    WRITE(unit=logfile,fmt=*) "***********************"
    WRITE(unit=logfile,fmt=*) "Correlations parameters"
    WRITE(unit=logfile,fmt=*) "***********************"
    WRITE(unit=logfile,fmt=*) CorrBl_keyw, " = ", NCorrBl
    DO i=1,NCorrBl
       WRITE(unit=logfile,fmt='(A,I2,A,I4,A,I4,A,F5.3,A,F5.3)'), " CorrBl. #", i, ":  ", CorrBeg(i), " - ", CorrEnd(i), " U=", UCoul(i), " J=", JHund(i)
    END DO
    WRITE(unit=logfile,fmt=*) POrtho_keyw, " = ", POrtho
    WRITE(unit=logfile,fmt=*) HybFunc_keyw, " = ", HybFunc
    WRITE(unit=logfile,fmt=*) DFTU_keyw, " = ", DFTU
    WRITE(unit=logfile,fmt=*) DiagCorrBl_keyw, " = ", DiagCorrBl
    WRITE(unit=logfile,fmt=*) "*****************"
    WRITE(unit=logfile,fmt=*) "Output parameters"
    WRITE(unit=logfile,fmt=*) "*****************"
    WRITE(unit=logfile,fmt=*) Hamilton_keyw, " = ", Hamilton
    WRITE(unit=logfile,fmt=*) Mulliken_keyw, " = ", Mulliken
    WRITE(unit=logfile,fmt=*) HTransm_keyw, " = ", HTransm
    WRITE(unit=logfile,fmt=*) EStep_keyw, " = ", EStep
    WRITE(unit=logfile,fmt=*) DOSEnergy_keyw, " = ", DOSEnergy
    WRITE(unit=logfile,fmt=*) EW1_keyw, " = ", EW1
    WRITE(unit=logfile,fmt=*) EW2_keyw, " = ", EW2
    WRITE(unit=logfile,fmt=*) LDOS_Beg_keyw, " = ", LDOS_Beg
    WRITE(unit=logfile,fmt=*) LDOS_End_keyw, " = ", LDOS_End
    WRITE(unit=logfile,fmt=*) RedTransmB_keyw, " = ",RedTransmB
    WRITE(unit=logfile,fmt=*) RedTransmE_keyw, " = ",RedTransmE
    WRITE(unit=logfile,fmt=*) NChannels_keyw, " = ", NChannels
    WRITE(unit=logfile,fmt=*) LeadDOS_keyw, " = ", LeadDOS
    WRITE(unit=logfile,fmt=*) ANT1DInp_keyw, " = ", ANT1DInp
    WRITE(unit=logfile,fmt=*) PrtHatom_keyw, " = ", PrtHatom

    IF (NAtomEl(1)/=0 .and. (NAtomEl(2)==0.and.ElType(2)/='GHOST')) THEN
       print *, 'Bad specification of number of electrode atoms'
       stop
    ELSE IF (NAtomEl(2)/=0 .and. (NAtomEl(1)==0.and.ElType(1)/='GHOST')) THEN
       print *, 'Bad specification of number of electrode atoms'
       stop
    ELSE IF ((NAtomEl(1)/=0 .and. NEmbed(1)/=0 .and. NAtomEl(1) < NEmbed(1)) .or. (NAtomEl(2)/=0 .and. NEmbed(2)/=0 .and. NAtomEl(2) < NEmbed(2))) THEN
       print *, 'Bad specification of number of electrode atoms'
       stop
    ELSE
       continue
    END IF

  END SUBROUTINE write_parameters

END MODULE parameters

