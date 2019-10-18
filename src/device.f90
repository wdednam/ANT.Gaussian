!*********************************************************!
!*********************  ANT.G-2.5.0  *********************!
!*********************************************************!
!                                                         !
!  Copyright (c) by                                       !
!                                                         !
!  Juan Jose Palacios (1)                                 !
!  David Jacob (2)                                        !
!  Maria Soriano (1)                                      !
!  Angel J. Perez-Jimenez (3)                             !
!  Emilio SanFabian (3)                                   !
!  Jose Antonio Antonio Verges (4)                        !
!  Enrique Louis (5)                                      !
!  Wynand Dednam (5)                                      !
!                                                         !
! (1) Departamento de Fisica de la Materia Condensada     !
!     Universidad Autonoma de Madrid                      !    
!     28049 Madrid (SPAIN)                                !
! (2) Theory Department                                   !
!     Max-Planck-Institute for Microstructure Physics     !
!     Halle, 06120 (GERMANY)                              !
! (3) Departamento de Quimica Fisica                      !
!     Universidad de Alicante                             !
!     03690 Alicante (SPAIN)                              !
! (4) Insto. de Ciencias de Materiales de Madrid (ICMM)   !
!     Consejo Superior de Investigacion y Ciencia (CSIC)  !
!     28049 Madrid (SPAIN)                                !
! (5) Departamento de Fisica Aplicada                     !
!     Universidad de Alicante                             !    
!     03690 Alicante (SPAIN)                              !
!                                                         !
!*********************************************************!
  MODULE device
!*********************************************************!
!  Main module for computation of device Green's function !
!  and transport                                          !
!*********************************************************!
  use ANTCommon
  implicit none
  save

  private

  public :: DevNAOrbs, DevNSpin, DevShift, DevFockMat, DevDensMat, SetDevDensMat
  public :: LeadsOn, SwitchOnLeads, SecantOn, SwitchOnSecant, SwitchOffSecant
  public :: EvaluationOn, SwitchOnEvaluation
  public :: SwitchOnChargeCntr, SwitchOffChargeCntr, SwitchOffSpinLock, SwitchOnSpinLock
  public :: InitDevice, ReadDensMat, ReadFockMat, CleanUpDevice, InitElectrodes, Transport


  !*****************************
  ! Module's internal variables
  !*****************************

  ! *** Number of atomic orbitals, number of non-degenerate spin-channels ***
  integer :: NAOrbs, NSpin, DNAOrbs

  ! *** Total number of electrons in central device region
  integer :: NCDEl, NCDAO1, NCDAO2

  ! *** Actual electron charge for certain Fermi energy ***
  real*8 :: QAlpha, QBeta, Q_SOC

  ! *** Overlap matrix S of device ***
  real*8, dimension(:,:),allocatable :: SD, InvSD
  real*8, dimension(:,:), allocatable :: S_SOC

  ! *** Complex S^+1/2 matrix ***
  complex*16, dimension(:,:),allocatable :: SPH
  
  ! *** Hamiltonian and density matrix of device
  real*8, dimension(:,:,:),allocatable :: HD
  real*8, dimension(:,:,:),allocatable :: PD
  complex*16, dimension(:,:,:),allocatable :: PDOUT
  complex*16, dimension(:,:),allocatable :: H_SOC, PD_SOC
  complex*16, dimension(:,:),allocatable :: PDOUT_SOC

  ! *** Orthognalization matrix for device ***
  real*8, dimension(:,:),allocatable :: OD

  ! *** Energy shifts ***
  real*8 :: shift, ShiftUp, ShiftDown

  ! *** spin excess charge ***
  real*8 :: exspin

  ! *** internal spin variable 1=up,2=down ***
  integer  :: ispin
  
  ! *** Lowest and highest eigen value of Fock matrix ***
  real*8 :: LEV,HEV
 
  ! *** lower and upper band edge of electrode ***
  real*8 :: EMinEc, EMaxEc

  ! *** lower and upper energy bound ***
  real*8 :: EMin, EMax

  ! *** Density of states projected on atoms at the Fermi energy
  real*8, dimension(:,:), allocatable :: AtomDOSEF

  ! *** Control Switches ***
  logical :: Leads       = .false.
  logical :: Secant      = .false.
  logical :: Evaluation  = .false.
  logical :: UDTrans     = .false.
  logical :: ChargeCntr  = .false.
  logical :: SpinLock    = .true.

  ! Whether device Hamiltonain has been orthogonalized
  logical :: HDOrtho = .false.
  
!!$OMP THREADPRIVATE(shift)
  contains

  !********************************
  !*** Public module procedures ***
  !********************************

  !**************************************
  ! Access functions to private entities
  !**************************************


  ! *** Total number of atomic orbitals in device Hilbert space ***
  integer function DevNAOrbs()
    implicit none
    DevNAOrbs=NAOrbs
  end function DevNAOrbs

  ! *** Number of non-degenerate spin bands ***
  integer function DevNSpin() 
    implicit none
    DevNSpin=NSpin
  end function DevNSpin

  ! *** Spin excess charge to equilibrate Fermi level up and down ***
  real*8 function DevShift()
    implicit none
    DevShift=-shift
  end function DevShift

  ! *** Get matrix Element of Density matrix ***
  real*8 function DevFockMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevFockMat = HD(is, i, j)
  end function DevFockMat

  ! *** Get matrix Element of Density matrix ***
  real*8 function DevDensMat( is, i, j )
    implicit none
    integer, intent(in) :: is, i, j
    DevDensMat = PD(is, i, j)
  end function DevDensMat

  ! *** Set matrix Element of Density matrix ***
  subroutine SetDevDensMat( is, i, j, pij )
    implicit none
    integer, intent(in) :: is, i, j
    real*8, intent(in) :: pij
    PD(is, i, j) = pij
  end subroutine SetDevDensMat

  ! ***
  logical function LeadsOn()
    implicit none
    LeadsOn = Leads
  end function LeadsOn

  ! *** 
  subroutine SwitchOnLeads
    implicit none
    Leads = .true.
  end subroutine SwitchOnLeads

  ! ***
  logical function SecantOn()
    implicit none
    SecantOn = Secant
  end function SecantOn

  ! ***
  subroutine SwitchOnSecant()
    implicit none
    Secant = .true.
  end subroutine SwitchOnSecant

  ! ***
  subroutine SwitchOffSecant()
    implicit none
    Secant = .false.
  end subroutine SwitchOffSecant

  ! ***
  logical function EvaluationOn()
    implicit none
    EvaluationOn = Evaluation
  end function EvaluationOn

  ! ***
  subroutine SwitchOnEvaluation()
    implicit none
    Evaluation = .true.
  end subroutine SwitchOnEvaluation

  subroutine SwitchOnChargeCntr()
    implicit none
    ChargeCntr = .true.
  end subroutine SwitchOnChargeCntr

  ! ***
  subroutine SwitchOffChargeCntr()
    implicit none
    ChargeCntr = .false.
  end subroutine SwitchOffChargeCntr

  ! ***
  subroutine SwitchOffSpinLock()
    implicit none
    SpinLock = .false.
  end subroutine SwitchOffSpinLock

  ! ***
  subroutine SwitchOnSpinLock()
    implicit none
    SpinLock = .true.
  end subroutine SwitchOnSpinLock
  

  !***********************************************
  !* Initialize device for transport calculation *
  !***********************************************
  subroutine InitDevice( NBasis, UHF, S )
    use constants, only: d_zero
    use numeric, only: RMatPow
    use parameters, only: ElType, FermiStart, Overlap, HybFunc, SOC, biasvoltage
    use cluster, only: AnalyseCluster, AnalyseClusterElectrodeOne, AnalyseClusterElectrodeTwo, NAOAtom, NEmbedBL
    use g09Common, only: GetNAtoms, GetAtmChg
    use correlation
    use orthogonalization
    use ANTCommon
    
    implicit none

    integer, intent(in) :: NBasis
    logical, intent(in) :: UHF
    real*8, dimension(NBasis,NBasis),intent(in) :: S

    integer :: AllocErr, ios, iatom, NEmbed1, NEmbed2

    real*8, dimension(NBasis,NBasis) :: RSPH 

    write(ifu_log,*) "-------------------"
    write(ifu_log,*) "Initializing device"
    write(ifu_log,*) "-------------------"

    NAOrbs = NBasis
    DNAOrbs = 2*NBasis
    NSpin=1
    if( UHF ) NSpin=2
    exspin=d_zero

    !if (biasvoltage /= 0.0) allocate(PDOUT(NSpin,NAOrbs,NAOrbs), STAT=AllocErr) 
    allocate(PDOUT(NSpin,NAOrbs,NAOrbs), STAT=AllocErr) 
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Allocation error for PDOUT"
       stop
    end if
    ! Dynamic arrays 
    allocate( SD(NAOrbs,NAOrbs), InvSD(NAOrbs,NAOrbs),  HD(NSpin,NAOrbs,NAOrbs),  PD(NSpin,NAOrbs,NAOrbs),  STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Allocation error for SD, InvSD, SMH, SPH, H, P"
       stop
    end if

    SD = S
    call RMatPow( SD, -1.0d0, InvSD )

    if( HybFunc ) call InitCorrelation(NAOrbs,NSpin)

    allocate( SPH(NAorbs,NAOrbs), STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/InitDevice: Allocation error for SPH(:,:)"
       stop
    end if
    ! Computing transformation matrix S^+1/2
    call RMatPow( SD,  0.5d0, RSPH )
    SPH = RSPH

    shift = -FermiStart
    shiftup = shift
    shiftdown = shift

    IF( ElType(1) == "BETHE" .and. ElType(2) == "BETHE" ) THEN 
      call AnalyseCluster
    ELSE IF  (ElType(1) == "BETHE" .and. ElType(2) == "GHOST" ) THEN
      call AnalyseClusterElectrodeOne
    ELSE IF  (ElType(1) == "GHOST" .and. ElType(2) == "BETHE" ) THEN
      call AnalyseClusterElectrodeTwo
    ELSE IF  (ElType(1) == "GHOST" .and. ElType(2) == "GHOST" ) THEN
      continue                           
    ELSE 
      print *, 'These electrodes are not implemented yet !!!'
      stop
    END IF

    call InitElectrodes

    EMin = 0.0d0
    EMax = 0.0d0

    LEV = 0.0d0
    HEV = 0.0d0

    ! Compute number of electrons 
    ! in central device region
    

    IF (Overlap < 0.01) THEN
       NEmbed1=0
       NEmbed2=0
    ELSE
       NEmbed1=NEmbedBL(1)
       NEmbed2=NEmbedBL(2)
    END IF

   
    print *, "---------------------------------------------------"
    print *, " Details on device and contacting atoms -----------"
    print *, "---------------------------------------------------"

    print *, "NEmbed(1) =", NEmbedBL(1)
    print *, "NEmbed(2) =", NEmbedBL(2)

    NCDEl = 0
    do iatom=NEmbed1+1,GetNAtoms()-NEmbed2
      NCDEl = NCDEl + GetAtmChg(iatom)
    end do
    
    print *, "Number of electrons in neutral reduced device"
    print *, "NCDEl = ", NCDEl
    
    print *, "First and last orbital in reduced device"
    IF (Overlap < 0.01) THEN
       NCDAO1 = 1
    ELSE
       NCDAO1 = 0
    END IF
    do iatom=1,NEmbed1
       NCDAO1 = NCDAO1 + NAOAtom(iatom) 
    end do

    print *, "NCDAO1 = ", NCDAO1

    NCDAO2 = NCDAO1-1
    do iatom = NEmbed1+1,GetNAtoms()-NEmbed2
       NCDAO2 = NCDAO2 + NAOAtom(iatom) 
    end do

    IF  (ElType(1) == "GHOST" .and. ElType(2) == "GHOST" ) NCDAO2 = NAOrbs
    print *, "NCDAO2 = ", NCDAO2
    print *, "---------------------------------------------------"

  end subroutine InitDevice
  
  !***********************************************
  !* Read initial density matrix from file P.dat *
  !***********************************************
  subroutine ReadDensMat(densitymatrix)
    use parameters, only: NSpinEdit, SpinEdit, MRStart, SpinDel, PFix, NFix, IFix, densitymatrixx
    use constants, only: d_zero
    use numeric, only: RMatPow
    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
    use g09Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
    use ANTCommon
    implicit none
    
    integer :: norb, ni, nj, isp, n, iatom, jatom, is, i, j, ii, jj, AOStart, AO_BEG, AO_END, iorb, jorb, iato, jato
    real*8 :: density, TrP, xxx
    real*8, dimension(NSpin,NAOrbs,NAOrbs) :: PDMod
    real*8, dimension(NAOrbs,NAOrbs) :: PDAux
    character (len=50) :: densitymatrix
  
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"
    write(ifu_log,*) "Starting calculation with density matrix from file  ", densitymatrix
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"

    PD=0.0d0
    PDMod=0.0d0

    open(ifu_dm,file=densitymatrix,status='old')
    read(ifu_dm,*,end=21) shift
    shift = -shift
    do i=1,2*NAOrbs*NAOrbs
       read(ifu_dm,*,end=21)is,ii,jj,density
       PDMod(is,ii,jj)=density
       PDMod(is,jj,ii)=density
    end do
 21 close (ifu_dm)

    if (PFix) then
       write(ifu_log,*) "-------------------------------------------------------------------------------------------"
       write(ifu_log,*) "  ... and using supplementary density matrix from file  ", densitymatrixx
       write(ifu_log,*) "-------------------------------------------------------------------------------------------"
       if (NFix == 0) print *,'Warning ... NFix = 0'
       open(ifu_dmx,file=densitymatrixx,status='old')
       read(ifu_dmx,*) xxx   
       do norb=1,8*NAOrbs*NAOrbs
          read(ifu_dmx,*,end=22)is,ni,nj,density,iato,iorb,jato,jorb
          do isp=1,NSpin
             i=0
             do iAtom=1,GetNAtoms()
                do n=1,NFix
                   if (iAtom == IFix(n)) then 
                      i=i+NAOAtom(iAtom)
                      goto 11
                   end if
                end do
                do ii=1,NAOAtom(iAtom)
                   i=i+1
                   j=0
                   do jAtom=1,GetNAtoms()
                      do n=1,NFix
                         if (jAtom == IFix(n)) then 
                            j=j+NAOAtom(jAtom)
                            goto 12
                         end if
                      end do
                      do jj=1,NAOAtom(jAtom)
                         j=j+1
                         if (is == isp .and. iato == iAtom .and. jato == jAtom .and. iorb == ii .and. jorb == jj) then
                            PDMod(isp,i,j)=density
                            PDMod(isp,j,i)=density
                         end if
                      end do
12                 end do
                end do
11           end do
          end do
       end do
22     close(ifu_dmx)
       
    end if
  
    !open(111,file='readdensitymatrix',status='unknown')
    !do is=1,NSpin
    !i=0
    !do iAtom=1,GetNAtoms()
    ! do ii=1,NAOAtom(iAtom)
    !    i=i+1
    !    j=0
    !    do jAtom=1,GetNAtoms()
    !       do jj=1,NAOAtom(jAtom)
    !          j=j+1
    !          write(111,'(i3,4i5,e18.10)')is,iAtom,ii,jAtom,jj,PDMod(is,i,j)
    !end do
    !end do
    !end do
    !end do
    !end do
    !close(111)

    !
    ! Manipulate atomic spins of initial guess 
    ! if SpinEdit or MRStart set
    !
    if( NSpinEdit > 0 .or. MRStart > 0 .or. SpinDel )then

       if( MRStart > 0 )then
          do iAtom=MRStart,GetNAtoms()
             SpinEdit( iAtom ) = -1
          end do
       end if

       if( SpinDel )then
          do iAtom=1,GetNAtoms()
             SpinEdit( iAtom ) = 0
          end do
       end if

       PD = PDMod
       PDMod = d_zero

       do iAtom=1,GetNAtoms()
          do jAtom=1,GetNAtoms()
             if( SpinEdit(iAtom) ==  1 .and. SpinEdit(jAtom) ==  1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)
                      PDMod(1,i,j) = PD(1,i,j)
                      PDMod(2,i,j) = PD(2,i,j)
                   end do
                end do
             else if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      PDMod(2,i,j) = PD(1,i,j)
                      PDMod(1,i,j) = PD(2,i,j)
                   end do
                end do
             else
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      PDMod(1,i,j) = 0.5d0*(PD(1,i,j)+PD(2,i,j))            
                      PDMod(2,i,j) = 0.5d0*(PD(1,i,j)+PD(2,i,j))            
                   end do
                end do
             end if
          end do
       end do

    end if

    ! Normalize to correct number of electrons when recommended
     
    TrP = d_zero
    do is=1,NSpin
       PDAux=MATMUL(PDMod(is,:,:),SD)
       do i=1,NAOrbs
          TrP = TrP + PDAux(i,i)
       end do
    end do
    if (NSpin ==1) then
       PRINT *, "Tr[P*S] of initial guess =", TrP*2.0
       if( NSpinEdit > 0 .or. MRStart > 0 .or. PFIX) then
         PD = PDMod*GetNE()/(TrP*2.0)
       else
         PD = PDMod
       end if
    else
       PRINT *, "Tr[P*S] of initial guess =", TrP
       if( NSpinEdit > 0 .or. MRStart > 0 .or. PFIX) then
          PD = PDMod*GetNE()/(TrP)
       else
          PD = PDMod
       end if
    end if
    PRINT *, "--------------------------------"
       
  end subroutine ReadDensMat

  !********************************************
  !* Read initial Fock matrix from file F.dat *
  !********************************************
  subroutine ReadFockMat(fockmatrix)
    use parameters, only: NSpinEdit, SpinEdit, MRStart, SpinDel !!!, PFix, NFix, IFix, densitymatrixx
    use constants, only: d_zero
    use numeric, only: RMatPow
    use cluster, only: LoAOrbNo, HiAOrbNo,NAOAtom
    use g09Common, only: GetNE, GetNAtoms, GetNAE, GetNBE
    use ANTCommon
    implicit none
    
    integer :: norb, ni, nj, isp, n, iatom, jatom, is, i, j, ii, jj, AOStart, AO_BEG, AO_END, iorb, jorb, iato, jato
    real*8 :: fock !, TrP, xxx

    real*8, dimension(NSpin,NAOrbs,NAOrbs) :: HDMod
    character (len=50) :: fockmatrix
  
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"
    write(ifu_log,*) "Starting calculation with fock matrix from file  ", fockmatrix
    write(ifu_log,*) "-------------------------------------------------------------------------------------------"

    HD=0.0d0
    HDMod=0.0d0

    open(ifu_fm,file=fockmatrix,status='old')
    read(ifu_fm,*,end=21) shift
    shift = -shift
    do i=1,2*NAOrbs*NAOrbs
       read(ifu_fm,*,end=21)is,ii,jj,fock
       HD(is,ii,jj)=fock
       HD(is,jj,ii)=fock
    end do
21 close (ifu_fm)

    !
    ! Manipulate atomic spins of initial guess 
    ! if SpinEdit or MRStart set
    !
    if( NSpinEdit > 0 .or. MRStart > 0 .or. SpinDel )then

       if( MRStart > 0 )then
          do iAtom=MRStart,GetNAtoms()
             SpinEdit( iAtom ) = -1
          end do
       end if

       if( SpinDel )then
          do iAtom=1,GetNAtoms()
             SpinEdit( iAtom ) = 0
          end do
       end if

       HDMod = d_zero

       do iAtom=1,GetNAtoms()
          do jAtom=1,GetNAtoms()
             if( SpinEdit(iAtom) ==  1 .and. SpinEdit(jAtom) ==  1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)
                      HDMod(1,i,j) = HD(1,i,j)
                      HDMod(2,i,j) = HD(2,i,j)
                   end do
                end do
             else if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      HDMod(2,i,j) = HD(1,i,j)
                      HDMod(1,i,j) = HD(2,i,j)
                   end do
                end do                
             else
                do i=LoAOrbNo(iAtom),HiAOrbNo(iAtom)
                   do j=LoAOrbNo(jAtom),HiAOrbNo(jAtom)                
                      HDMod(1,i,j) = 0.5d0*(HD(1,i,j)+HD(2,i,j))            
                      HDMod(2,i,j) = 0.5d0*(HD(1,i,j)+HD(2,i,j))            
                   end do
                end do
             end if
          end do
       end do

       HD = HDMod

    end if

  end subroutine ReadFockMat


  !**************************
  ! Deallocate dynamic arrays
  !**************************
  subroutine CleanUpDevice
    use BetheLattice, only: CleanUpBL, LeadBL
    use parameters, only: ElType
    integer :: AllocErr, LeadNo

    deallocate( SD, InvSD, HD, PD, PDOUT, SPH, STAT=AllocErr )
    if( AllocErr /= 0 ) then
       print *, "DEVICE/Deallocation error for SD, InvSD, HD, PD, PDOUT, SPH"
       stop
    end if
    do LeadNo=1,2
       select case( ElType(LeadNo) )
       case( "BETHE" )
          call CleanUpBL( LeadBL(LeadNo) ) 
       end select
    end do
  end subroutine CleanUpDevice


  !*************************
  !* Initialize electrodes *
  !*************************
  subroutine InitElectrodes
    use BetheLattice, only: InitBetheLattice, LeadBL, BL_EMin, BL_EMax
    use OneDLead, only: Init1DLead, Lead1d, L1D_EMin, L1D_EMax
    use parameters, only: ElType
    implicit none 
    integer :: LeadNo
    real*8, dimension(2) :: EMin, EMax

    do LeadNo=1,2
       select case( ElType(LeadNo) )
       case( "BETHE" )
          call InitBetheLattice( LeadBL(LeadNo), LeadNo )
          EMin(LeadNo) = BL_EMin( LeadBL(LeadNo) )
          EMax(LeadNo) = BL_EMax( LeadBL(LeadNo) )
       case( "1DLEAD" )
          call Init1DLead( Lead1d(LeadNo), LeadNo )
          EMin(LeadNo) = L1D_EMin( Lead1d(LeadNo) )
          EMax(LeadNo) = L1D_EMax( Lead1d(LeadNo) )
       case( "GHOST" )
          EMin(LeadNo) = -100.0                       
          EMax(LeadNo) =  100.0                       
       case DEFAULT
          print *, "Unknown option for ElType:", ElType(LeadNo)
          stop
       end select
    end do
    
    EMinEc = min( Emin(1),EMin(2) )
    EMaxEc = max( EMax(1),EMax(2) )
  end subroutine InitElectrodes

  
  !***************************
  !* Solve transport problem *
  !***************************
  subroutine Transport(F,ADDP) 
    use parameters, only: RedTransmB, RedTransmE, ANT1DInp, ElType, HybFunc, POrtho, DFTU, DiagCorrBl, DMImag, LDOS_Beg, LDOS_End, NSpinEdit, SpinEdit, SOC, PrtHatom
    use numeric, only: RMatPow, RSDiag
    use cluster, only: LoAOrbNo, HiAOrbNo
    use correlation
    use orthogonalization
    use g09Common, only: GetNAtoms
    implicit none

    logical,intent(out) :: ADDP
    real*8, dimension(NSpin,NAOrbs,NAOrbs),intent(in) :: F

    real*8, dimension(NAOrbs) :: evals
    real*8, dimension(:,:),allocatable :: SPM

    real*8 :: diff !!,TrP,QD
    integer :: i,j,is, info, AllocErr, iatom, jatom, Atom

    HD = F
    !
    ! Estimate upper bound for maximal eigenvalue of HD and use it for upper and lower energy boundaries
    !
    if( NSpin == 1 ) EMax = maxval(sum(abs(HD(1,:,:)),1))
    if( NSpin == 2 ) EMax = max(maxval(sum(abs(HD(1,:,:)),1)),maxval(sum(abs(HD(2,:,:)),1)))
    EMin = -EMax
    print *, "EMin=", EMin
    print *, "EMax=", EMax

    if( .not. DMImag .and. ChargeCntr )then
       print *, "--------------"
       print *, "Charge Control"
       print *, "--------------"
       ! Find upper and lower energy bound 
       call FindEnergyBounds
    endif

    if( DFTU ) call Add_DFT_plus_U_Pot( PD, HD )

    if(.not.DMImag) call CompDensMat(ADDP)
    if(DMImag) call CompDensMat2(ADDP)

    if( Evaluation )then
       print *
       print *, "****************************************** "
       print *, "*                                        * "
       print *, "*        ANT.G09 final analysis          * "
       print *, "*                                        * "       
       print *, "****************************************** "
       print *

       IF( ANT1DInp ) call WriteANT1DInput

       if( POrtho )then
          allocate( OD(NAorbs,NAOrbs), STAT=AllocErr )
          if( AllocErr /= 0 ) then
             print *, "DEVICE/InitDevice: Allocation error for OD(:,:)"
             stop
          end if
          do ispin=1,NSpin
             PD(ispin,:,:) = matmul( SD, matmul(PD(ispin,:,:), SD ) )
          end do
          HDOrtho = .true.
          call ProjOrtho(cix, SD, OD )
          ! Othogonalize density matrix density matrix and Hamiltonian
          do ispin=1,NSpin
             HD(ispin,:,:) = matmul( transpose(OD), matmul(F(ispin,:,:), OD) )
             PD(ispin,:,:) = matmul( OD, matmul(PD(ispin,:,:), transpose(OD) ) )
          end do
          HDOrtho = .true.
       end if
       if( DiagCorrbl ) call DiagCorrBlocks( HD, SD )
       call Hamiltonian
       IF ( HybFunc ) call CompHybFunc
       IF ((ElType(1) == "GHOST" .or. ElType(2) == "GHOST") .and. LDOS_Beg <= LDOS_End) CALL LDOS
       IF (ElType(1) /= "GHOST" .and. ElType(2) /= "GHOST") THEN            
          IF( RedTransmE >= RedTransmB  ) call EigenChannelAnalysis
          call transmission
       END IF
       
       !
       ! Print hamiltonian of atoms 
       ! with SpinEdit set
       !
       if( NSpinEdit > 0 )then                         
       
       Atom = PrtHatom
 
          do iAtom=1,GetNAtoms()
             do jAtom=1,GetNAtoms()
                !if( SpinEdit(iAtom) == -1 .and. SpinEdit(jAtom) == -1 )then
                if( iAtom == Atom .and. jAtom == Atom )then
                   PRINT *, " Hamiltonian of edited spin (atom) ",iAtom," is: "
                   PRINT *, " Up-Up "  
                   do i=LoAOrbNo(Atom),HiAOrbNo(Atom)                                                               
                      PRINT '(1000(F11.5))', ( (HD(1,i,j)), j=LoAOrbNo(Atom),HiAOrbNo(Atom) )               
                   end do                                                                                                       
                   PRINT *, " Down-Down "                                                                              
                   do i=LoAOrbNo(Atom),HiAOrbNo(Atom)                                                            
                      PRINT '(1000(F11.5))', ( (HD(2,i,j)), j=LoAOrbNo(Atom),HiAOrbNo(Atom) )             
                   end do                                                                              
                end if   
             end do
          end do
       end if   
        
       if (SOC) then 
          call MullPop_SOC
       else 
          call MullPop
       end if   
       return
    end if
  end subroutine transport
  
  subroutine WriteANT1DInput
    use parameters, only: eta
    use AntCommon
    implicit none

    real*8 :: dsmall
    integer :: i,j

    CHARACTER(len=55) :: fname

    dsmall = eta
    
    fname='dev.'//trim(ant1dname)//'.dat'
    open(ifu_ant,file=fname,status='unknown')

    write(ifu_ant,'(A)'),       "&DevParams"
    write(ifu_ant,'(A,I1)'),    "NDSpin = ", NSpin
    write(ifu_ant,'(A,I4)'),    "NDAO = ", NAOrbs
    write(ifu_ant,'(A,I4)'),    "NDEl = ", NCDEl
    write(ifu_ant,'(A,F12.8)'), "EFermi = ", -shift
    write(ifu_ant,'(A)'),       "sparse = .true."
    write(ifu_ant,'(A)'),       "/"
    write(ifu_ant,*)
    write(ifu_ant,'(A)'),       "! Hamiltonian"
    do ispin=1,NSpin
       if( NSpin == 2 .and. ispin == 1 ) write(ifu_ant,'(A)'),       "! Spin-up"
       if( NSpin == 2 .and. ispin == 2 ) write(ifu_ant,'(A)'),       "! Spin-down"
       do i=1,NAOrbs
          do j=1,NAOrbs
             if(abs(HD(ispin,i,j))>=dsmall) write(ifu_ant,'(I6,I6,ES20.8)'), i, j, HD(ispin,i,j)
          enddo
       enddo
       write(ifu_ant,'(I6,I6,ES20.8)'), 0, 0, 0.0d0
       write(ifu_ant,*)
    end do
    write(ifu_ant,'(A)'),       "! Overlap"
    do i=1,NAOrbs
       do j=1,NAOrbs
          if(abs(SD(i,j))>=dsmall) write(ifu_ant,'(I6,I6,ES20.8)'), i, j, SD(i,j)
       enddo
    enddo
    write(ifu_ant,'(I6,I6,ES20.8)'), 0, 0, 0.0d0
    write(ifu_ant,*)
    close(ifu_ant)
  end subroutine WriteANT1DInput
  
  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
  subroutine CompDensMat(ADDP)
    use numeric, only: SECANT, MULLER, BISEC
    use g09Common, only: GetNAE, GetNBE
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real*8 :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail
    
    Z=10.0d0*FermiAcc
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    if( NSpin == 2 .and. SPINLOCK )then
       print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
       do ispin=1,NSpin
          root_fail = .true.
          if (ispin.eq.1) E1=shiftup
          if (ispin.eq.2) E1=shiftdown
          E0=E1-Z
          E2=E1+Z
          if( root_fail )then
             print*,'MULLER method'
             call MULLER(F,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
                print *, 'Warning: MULLER method failed to find root. Using SECANT.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if( root_fail )then
             print*,'SECANT method'
             call SECANT(F,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax) then
                print *, 'Warning: SECANT method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if (root_fail) then
             print *, 'BISEC method'
             if (ispin.eq.1) shiftup = BISEC(F,EMin,EMax,Delta,5*Max,K)
             if (ispin.eq.2) shiftdown = BISEC(F,EMin,EMax,Delta,5*Max,K)
             DE=Delta
             if(k.lt.5*Max) root_fail = .false.
             if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
          end if
          write(ifu_log,*)'--------------------------------------------------------'
          if (ispin.eq.1) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -shiftup,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
          end if
          if (ispin.eq.2) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -shiftdown,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
          end if
          write(ifu_log,*)'--------------------------------------------------------'
       end do
    else
       root_fail = .true.
       E0=shift-Z 
       E1=shift
       E2=shift+Z
       if (root_fail) then
          print*,'MULLER method'
          call MULLER(QXTot,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
             print *, 'Warning: MULLER method failed to find root. Using SECANT.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print*,'SECANT method'
          call SECANT(QXTot,E0,E2,Delta,Epsilon,Max,E3,DE,Cond,K)
          if(k .eq. Max .or. E3<EMin .or. E3>EMax) then
             print *, 'Warning: SECANT method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if (root_fail) then
          print *, 'BISEC method'
          shift = BISEC(QXTot,EMin,EMax,Delta,5*Max,K)
          DE=Delta
          if(k.lt.5*Max) root_fail = .false.
          if(k.ge.5*Max) print *, 'Warning: BISECT method failed to find root. Skipping this cycle.'
       end if

       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
       write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
       write(ifu_log,'(A,F10.5)') ' Total number of electrons:  ', QAlpha+QBeta
       write(ifu_log,*)'-----------------------------------------------'
    end if
    ADDP = .not. root_fail
  end subroutine CompDensMat

  !**************************************************************
  !* Subroutine for determining Fermi energy and density matrix *
  !* for some total charge                                      *
  !**************************************************************
  !* Pre-condition:                                             *
  !*   HC: Fock-matrix                                          *
  !*   shiftup,shiftdown: starting values for Fermi-energies    *
  !* Results:                                                   *
  !*   PC: Density-matrix                                       *
  !*   shiftup,shiftdown: Fermi-energies                        *
  !**************************************************************
  subroutine CompDensMat2(ADDP)
    use numeric, only: SECANT, MULLER, BISEC
    use g09Common, only: GetNAE, GetNBE
    use parameters, only: FermiAcc,ChargeAcc,Max,QExcess
    !use ieee_arithmetic
    implicit none

    logical,intent(out) :: ADDP
    integer :: i,j, k,cond
    real*8 :: E0,E1,E2,E3,DE,Z, Delta, Epsilon
    
    logical :: root_fail
    
    Z=10.0d0*FermiAcc
    Delta=FermiAcc
    Epsilon=ChargeAcc*(NCDEl+QExcess)

    root_fail = .true.
    if( NSpin == 2 .and. SPINLOCK )then
       print'(A,f10.1)',' SPINLOCK: NAlpha-NBeta = ',dble((GetNAE()-GetNBE())*(NCDEl+QExcess))/dble(GetNAE()+GetNBE())
       do ispin=1,NSpin
          if (ispin.eq.1) E0=shiftup
          if (ispin.eq.2) E0=shiftdown
          E1=E0-Z
          E2=E0+Z 
          if( root_fail ) then
             print*,'Secant method'
             call SECANT(CompSpinPD,E0,E1,Delta,Epsilon,Max,E2,DE,Cond,K)
             if(k.eq.Max .or. E2<EMin .or. E2>EMax)then
                print *, 'Warning: Secant method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E2
                if (ispin.eq.2)shiftdown=E2
                root_fail = .false.
             end if
          end if
          if( root_fail ) then
             print*,'Muller method'
             call MULLER(CompSpinPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
             if(k.eq.Max .or. E3<EMin .or. E3>EMax)then
                print *, 'Warning: Muller method failed to find root. Using BISEC.'
                root_fail = .true.
             else
                if (ispin.eq.1)shiftup=E3
                if (ispin.eq.2)shiftdown=E3
                root_fail = .false.
             end if
          end if
          if( root_fail )then
             print *, 'BISEC method'
             if( ispin.eq.1) shiftup = BISEC(CompSpinPD,EMin,EMax,Delta,Max,K)
             if( ispin.eq.2) shiftdown = BISEC(CompSpinPD,EMin,EMax,Delta,Max,K)
             DE=Delta
             if(k.lt.Max) root_fail = .false.
          end if
          write(ifu_log,*)'--------------------------------------------------------'
          if (ispin.eq.1) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for alpha electrons= ', -shiftup,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
          end if
          if (ispin.eq.2) then
             write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy for beta electrons=  ', -shiftdown,'  +/-',dabs(DE)
             write(ifu_log,*)
             write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
          end if
          write(ifu_log,*)'--------------------------------------------------------'
       end do
    else
       E0=shift
       E1=E0-Z 
       E2=E0+Z 
       if( root_fail )then
          print*,'Secant method'
          call SECANT(CompPD,E0,E1,Delta,Epsilon,Max,E2,DE,Cond,K)
          if(k.eq.Max .or. E2<EMin .or. E2>EMax)then
             print *, 'Warning: Secant method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E2
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print*,'Muller method'
          call MULLER(CompPD,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
          if(k.eq.Max .or. E3<EMin .or. E3>EMax)then
             print *, 'Warning: Muller method failed to find root. Using BISEC.'
             root_fail = .true.
          else
             shift = E3
             root_fail = .false.
          end if
       end if
       if( root_fail )then
          print *, 'BISEC method'
          print *, EMin, EMax
          shift = BISEC(CompPD,EMin,EMax,Delta,Max,K)
          DE=Delta
          if(k.lt.Max) root_fail = .false.
       end if
       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of alpha electrons: ', QAlpha
       write(ifu_log,'(A,F10.5)') ' Number of beta electrons:  ', QBeta
       write(ifu_log,*)'-----------------------------------------------'
    end if
    ADDP = .not. root_fail
  end subroutine CompDensMat2
  
  ! 
  ! Computes the density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real*8 function CompPD( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
!    USE IFLPORT
     use omp_lib
    use g09Common, only: GetNAE, GetNBE
     implicit none

     ! chemical potential
     real*8, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
     real*8, dimension(2*npmax) :: x, w
     real*8 :: Ei, dEdx, Q, QQ, DPD, E0 !, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     
     real*8, parameter :: x0 = 0.5d0
     real*8 :: aa, bb, cc, EM

     EM = EMax
     aa = EM/x0
     bb = aa*(1.0d0-x0)**2
     cc = EM - bb/(1-x0)

     shift = mu

     Q = d_zero

     do n=nmin,nmax

        np=2**n-1
        
        QQ = Q
        PD=d_zero
        QAlpha = d_zero; QBeta = d_zero
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

        do ispin=1,NSpin
!$OMP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMP DO
           do i=1,np
              Ei = 2.0d0*EMax*x(i)
              if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
              dEdx = 2.0d0*EMax
              if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
              call GPlus0( ui*Ei, GD )
!$OMP CRITICAL
              do k=1,NAOrbs
                 do l=1,NAOrbs
                    DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                    PD(ispin,k,l) = PD(ispin,k,l) + DPD
                    if(ispin.eq.1) QAlpha = QAlpha + DPD*SD(l,k)
                    if(ispin.eq.2) QBeta = QBeta + DPD*SD(l,k)
                 end do
              end do
!$OMP END CRITICAL
           end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        end do
        if(NSpin.eq.1) QBeta=QAlpha
        Q = QAlpha + QBeta        
        if( n > nmin .and. abs(Q-QQ) < PAcc*NCDEl ) exit
     end do
        
     print '(A,I4,A)', ' Integration of P has needed ', np, ' points.'
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     CompPD = Q - dble(NCDEl)

  end function CompPD

  
  ! 
  ! Computes the spin-resolved density matrix for a fixed chemical potential mu
  ! by integrating Greens function on matsubara axis. No lower energy
  ! bound required anymore!
  ! - Returns number of electrons in device region
  !
  ! - replaces old code in functions F(x) and QXTot(x)
  !
  real*8 function CompSpinPD( mu )
     use constants
     use util
     use numeric, only: CHDiag, gauleg, RTrace
     use parameters, only: eta, PAcc
!    USE IFLPORT
     use omp_lib
    use g09Common, only: GetNAE, GetNBE
     implicit none

     ! chemical potential
     real*8, intent(in) :: mu

     complex*16, dimension(NAOrbs,NAOrbs) :: GD

     integer, parameter :: nmin=1, nmax=8, npmax=2**nmax-1
     real*8, dimension(2*npmax) :: x, w
     real*8 :: Ei, dEdx, Q, QQ, DPD !, E0, Qi
     integer :: n, np, i, j, k, l !, info, ierr
     real*8 :: NCDAB

     ! Number of alpha OR beta electrons in central region C of device
     if(ispin.eq.1) NCDAB = dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE())
     if(ispin.eq.2) NCDAB = dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE())
     
     shift = mu

     Q=0.0d0
     do n=nmin,nmax
        np=2**n-1
        
        QQ = Q
        Q = 0.0d0
        PD(ispin,:,:) = 0.0d0
        
        ! Compute Gauss-Legendre abcsissas and weights
        call gauleg(0.0d0,2.0d0,x(1:2*np),w(1:2*np),2*np)

!$OMP PARALLEL PRIVATE(Ei,dEdx,GD,DPD)
!$OMP DO
        do i=1,np
           Ei = 2.0d0*EMax*x(i)
           if( x(i) > 0.5d0 ) Ei = 0.5d0*EMax/(1.0d0-x(i))
           dEdx = 2.0d0*EMax
           if( x(i) > 0.5d0 ) dEdx = 0.5d0*EMax/(1.0d0-x(i))**2   
           call GPlus0( ui*Ei, GD )
!$OMP CRITICAL
           do k=1,NAOrbs
              do l=1,NAOrbs
                 DPD = w(i)*(dEdx*real(GD(k,l))/d_pi + 0.5d0*InvSD(k,l))
                 PD(ispin,k,l) = PD(ispin,k,l) + DPD
                 Q = Q + DPD*SD(l,k)
              end do
           end do
!$OMP END CRITICAL
        end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
        if( n > nmin .and. abs(Q-QQ) < PAcc*NCDAB ) exit
     end do
        
     print '(A,I4,A)', ' Integration of P has needed ', np, ' points.'
     print '(A,F10.5,A,F10.5,A,F8.5)', ' mu =', -mu, '  Num. of electrons =', Q, ' +/-', abs(Q-QQ) 

     if(ispin.eq.1) QAlpha = Q
     if(ispin.eq.2) QBeta = Q

     CompSpinPD = Q - NCDAB

   end function CompSpinPD

  
  !*************************************
  !* Compute retarded Green's function *
  !*************************************
  subroutine gplus0(z,green)
    use PARAMETERS, only: eta,glue
    use constants, only: c_zero, ui
     use lapack_blas, only: zgetri, zgetrf
!   use lapack95, only: zgetri, zgetrf
    implicit none
    external zgetrf, zgetri    

    integer :: i, j, info, omp_get_thread_num
    integer, dimension(NAOrbs) :: ipiv
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr, temp 
    complex*16 :: work(4*NAOrbs) 

    complex*16, intent(in) :: z 

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: green
    

    ! Initilization 
    green=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

  end subroutine gplus0


  !********************************************************
  !* Compute retarded Green's function and gamma matrices *
  !********************************************************
  subroutine gplus(z,green,gammar,gammal)
    use PARAMETERS, only: eta, glue
    use constants, only: c_zero, ui
     use lapack_blas, only: zgetri, zgetrf
!   use lapack95, only: zgetri, zgetrf

    implicit none
    external zgetrf, zgetri 
    
    integer :: i, j, info, omp_get_thread_num
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl,sigr, temp 
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: green
    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gammar,gammal
    
    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************

    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

    !write(ifu_log,*)omp_get_thread_num(),green(1,1)
  end subroutine gplus

  !*************************************
  !* Compute lesser Green's function *
  !*************************************
  subroutine glesser(z,gless)
    use parameters, only: eta, biasvoltage, glue
    use constants, only: c_zero, ui, c_one
     use lapack_blas, only: zgetri, zgetrf, zgemm
!   use lapack95, only: zgetri, zgetrf, zgemm

    implicit none
    external zgetrf, zgetri, zgemm     

    integer :: i, j, info, error
    integer, dimension(NAOrbs) :: ipiv
    complex*16 :: work(4*NAOrbs)

    complex*16, intent(in) :: z

    complex*16, dimension(NAOrbs,NAOrbs), intent(out) :: gless
    complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL

     allocate(glessR(NAOrbs,NAOrbs),glessL(NAOrbs,NAOrbs),green(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then 
        print*,"Problems allocating" 
        stop
     end if
     allocate(sigl(NAOrbs,NAOrbs),sigr(NAOrbs,NAOrbs),gammar(NAOrbs,NAOrbs),gammal(NAOrbs,NAOrbs),stat=error)
     if (error /= 0 ) then 
        print*,"Problems allocating" 
        stop
     end if

    gammar=c_zero
    gammal=c_zero
    sigr=-ui*eta*SD 
    sigl=-ui*eta*SD 

    call CompSelfEnergies( ispin, z, sigl, sigr )
    sigr=glue*sigr
    sigl=glue*sigl
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    gammar=ui*(sigr-conjg(transpose(sigr)))
    gammal=ui*(sigl-conjg(transpose(sigl)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************
    do i=1,NAOrbs
       do j=1,NAOrbs
          green(i,j)=(z-shift)*SD(i,j)-HD(ispin,i,j)-sigl(i,j)-sigr(i,j)
       enddo
    enddo

    call zgetrf(NAOrbs,NAOrbs,green,NAOrbs,ipiv,info)
    call zgetri(NAOrbs,green,NAOrbs,ipiv,work,4*NAOrbs,info)

    !************************************************************************
    !* G< (Gless)
    !************************************************************************

     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammal,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessL,NAOrbs)
     call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,gammar,NAOrbs,green,NAOrbs, &
          &           c_zero,gless,NAOrbs)
     call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,green,NAOrbs,gless,NAOrbs, &
          &           c_zero,glessR,NAOrbs)

     gless=ui*(glessL*theta(dble(z)-biasvoltage/2)+glessR*theta(dble(z)+biasvoltage/2))

     deallocate(glessR,glessL,sigl,sigr,gammar,gammal,stat=error)
     if (error /= 0 ) then 
        print*,"Problems deallocating" 
        stop
     end if

  end subroutine glesser

  !*********************************************************************
  !* Compute lesser Green's function  with SOC (doubling the matrices) *
  !*********************************************************************
  subroutine glesser_SOC(z,gless)
    use parameters, only: eta, glue, biasvoltage
    use constants, only: c_zero, ui, c_one
     use lapack_blas, only: zgetri, zgetrf, zgemm
!   use lapack95, only: zgetri, zgetrf

    implicit none
    external zgetrf, zgetri, zgemm     

    integer :: i, j, info, omp_get_thread_num, error
    integer, dimension(DNAOrbs) :: ipiv
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: gless
    complex*16, dimension(DNAOrbs,DNAOrbs) :: gammar,gammal
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    complex*16, dimension(DNAOrbs,DNAOrbs) :: green,glessL,glessR

!   complex*16, dimension(:,:), allocatable :: green,gammar,gammal,sigl,sigr,glessR,glessL

!    allocate(glessR(NAOrbs,NAOrbs),glessL(NAOrbs,NAOrbs),green(NAOrbs,NAOrbs),stat=error)
!    if (error /= 0 ) then 
!       print*,"Problems allocating" 
!       stop
!    end if
!    allocate(sigl(NAOrbs,NAOrbs),sigr(NAOrbs,NAOrbs),gammar(NAOrbs,NAOrbs),gammal(NAOrbs,NAOrbs),stat=error)
!    if (error /= 0 ) then 
!       print*,"Problems allocating" 
!       stop
!    end if

   
    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigmar=c_zero
    sigmal=c_zero


    sigr1=-ui*eta*SD 
    sigl1=-ui*eta*SD 
    sigr2=-ui*eta*SD 
    sigl2=-ui*eta*SD 

    call CompSelfEnergies( 1, z, sigl1, sigr1 )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2 )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

    gammar=ui*(sigmar-conjg(transpose(sigmar)))
    gammal=ui*(sigmal-conjg(transpose(sigmal)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)


    !************************************************************************
    !* G< (Gless)
    !************************************************************************

     call zgemm('N','C',DNAOrbs,DNAOrbs,DNAOrbs,c_one,gammal,DNAOrbs,green,DNAOrbs, &
          &           c_zero,gless,DNAOrbs)
     call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one,green,DNAOrbs,gless,DNAOrbs, &
          &           c_zero,glessL,DNAOrbs)
     call zgemm('N','C',DNAOrbs,DNAOrbs,DNAOrbs,c_one,gammar,DNAOrbs,green,DNAOrbs, &
          &           c_zero,gless,DNAOrbs)
     call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one,green,DNAOrbs,gless,DNAOrbs, &
          &           c_zero,glessR,DNAOrbs)

     gless=ui*(glessL*theta(dble(z)-biasvoltage/2)+glessR*theta(dble(z)+biasvoltage/2))

!    deallocate(glessR,glessL,sigl,sigr,gammar,gammal,stat=error)
!    if (error /= 0 ) then 
!       print*,"Problems deallocating" 
!       stop
!    end if

  end subroutine glesser_SOC



  ! *************
  ! Step function
  ! *************
  function theta(x)
    implicit none
    real*8,intent(in) :: x
    real*8 :: theta
    theta=0.0d0
    if(x.lt.0.0d0) theta=1.0
    return
  end function theta

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  !****************************************
  double precision function f(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real*8, intent(in) :: x

    real*8 :: chargeup, chargedown, rrr, a, b, Q
    integer :: i,j,M

    shift=x
    
    ! Radius of complex contour integration
    ! add 10eV just in case 
    rrr = 0.5*abs(EMin)+10.0d0;

    !c c Integral limits ... (a,b)
    a = 0.d0
    b = d_pi
    M=1000
    call IntCompPlane(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))
      
    ! Density matrix out of equilibirum
    if (biasvoltage /= 0.0) then
       M=1000
       call IntRealAxis(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
       PD = PD + PDOUT
    end if
    ! Density matrix out of equilibirum
      
    Q=d_zero
    !do i=1,NAOrbs
    do i=NCDAO1, NCDAO2
       do j=1,NAOrbs
          Q=Q+PD(ispin,i,j)*SD(j,i)
       end do
    end do
    if (ispin.eq.1) QAlpha = Q
    if (ispin.eq.2) QBeta  = Q
 !if (ispin.eq.1) f = QAlpha - dble(GetNAE())
 !if (ispin.eq.2) f = QBeta  - dble(GetNBE())
    if (ispin.eq.1) f = QAlpha - dble(GetNAE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
    if (ispin.eq.2) f = QBeta  - dble(GetNBE()*NCDEl)/dble(GetNAE()+GetNBE()) -QExcess/2.0
    return
  end function f


  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x
  !****************************************
  double precision function QXTot(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real*8, intent(in) :: x

    real*8 :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

    do ispin=1,NSpin
       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+10.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000

       call IntCompPlane(rrr,a,b,M,d_zero-dabs(biasvoltage/2.0))

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.0) then
          M=1000
          call IntRealAxis(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          PD=PD+PDOUT
       end if
       ! Density matrix out of equilibirum

       !do i=1,NAOrbs
       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+PD(ispin,i,j)*SD(j,i)
          end do
       end do
       if( ispin == 1 ) QAlpha = Q
       if( ispin == 2 ) QBeta  = Q
    end do

    if( NSpin == 1 ) QBeta = QAlpha
    !QXTot = QAlpha + QBeta -dble(GetNAE()) -dble(GetNBE())
    QXTot = QAlpha + QBeta - dble(NCDEl) - QExcess
    return
  end function QXTot

  !****************************************
  ! Function gives excess charge of device
  ! in dependence of Fermi energy x with SOC
  !****************************************
  double precision function QXTot_SOC(x)
    use parameters, only: biasvoltage,QExcess
    use constants, only: d_pi, d_zero
    use g09Common, only: GetNAE, GetNBE
    implicit none

    real*8, intent(in) :: x

    real*8 :: rrr, a, b, Q
    integer :: i,j,M,omp_get_thread_num

       shift=x
    !write(ifu_log,*)omp_get_thread_num(),'in qxtot',shift
       Q = d_zero

       ! Radius of complex contour integration
       ! add 10eV just in case 
       rrr = 0.5*abs(EMin)+20.0d0;

       !c c Integral limits ... (a,b)
       a = 0.d0
       b = d_pi
       M=1000

       call IntRealAxis_SOC(Emin,-dabs(biasvoltage/2.0),M)
       PD_SOC = PDOUT_SOC

       ! Density matrix out of equilibirum
       if (biasvoltage /= 0.0) then
          M=1000
          call IntRealAxis_SOC(-dabs(biasvoltage/2.0),dabs(biasvoltage/2.0),M)
          PD_SOC=PD_SOC+PDOUT_SOC
       end if
       ! Density matrix out of equilibirum

       do i=NCDAO1, NCDAO2
          do j=1,NAOrbs
             Q=Q+REAL(PD_SOC(i,j)*S_SOC(j,i))
          end do
       end do
       do i=NCDAO1+NAOrbs, NCDAO2+NAOrbs
          do j=NAOrbs+1,DNAOrbs
             Q=Q+REAL(PD_SOC(i,j)*S_SOC(j,i))
          end do
       end do

    Q_SOC = Q
    QXTot_SOC = Q - dble(NCDEl) - QExcess
    return
  end function QXTot_SOC


  !ccccccccccccccccccccccccccccccc
  !c                                                                              c
  !c     Change of variable no. 3 for the numerical integration:                  c
  !c                                                                              c
  !c     int_{-infty}^{Eq} DOS(E)dE = int_{-1}^{1} DOS(E(x))*(dE/dx)dx            c
  !c                                                                              c
  !c     E = Em*(1-bx)/(1+x)                                                      c
  !c                                                                              c
  !ccccccccccccccccccccccccccccccc
  real*8 function edex3(Em,b,x)
    implicit none
    real*8 :: Em, b ,x
    edex3 = 0.5d0*((Em-b)*x + (Em+b))
    return
  end function edex3
  
  !******************************!
  ! Writing the Hamiltonian      !
  !******************************!
  subroutine Hamiltonian
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol, NEmbedBL
    USE parameters, only: Hamilton
    use g09Common, only: GetAtmCo
    use constants, only: Bohr
    implicit none

    integer :: i,j, I1, is ,n
    real*8 :: ro_a, ro_b

    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Self-consistent Hamiltonian  ---'
    write(ifu_log,*)'-------------------------------------'

    
    !open(ifu_ham,file='V.'//trim(xxx)//'.dat',status='unknown')
    write(ifu_log,*)'---------'
    write(ifu_log,*)'Left lead'
    write(ifu_log,*)'---------'
    I1=0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_b=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+hd(1,i,i)
          if(NSpin==2) ro_b=ro_b+hd(2,i,i)
       end do
       ro_a=ro_a/NAOAtom(j)
       ro_b=ro_b/NAOAtom(j)
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j, ro_a
       if(NSpin ==2 ) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
       IF (Hamilton) THEN
       if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
       if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
       END IF
       I1 = I1 + NAOAtom(j)
    end do
    
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'--------'
       write(ifu_log,*)'Molecule'
       write(ifu_log,*)'--------'
       do j = NALead(1)+1,NALead(1)+NAMol()
          ro_a = 0.0d0
          ro_b = 0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+hd(1,i,i)
             if(NSpin==2) ro_b=ro_b+hd(2,i,i)
          end do
          ro_a=ro_a/NAOAtom(j)
          ro_b=ro_b/NAOAtom(j)
          if(NSpin==1) write(ifu_log,1011)'Atom:',j, ro_a
          if(NSpin==2) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
          IF (Hamilton) THEN
          if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
          if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
          END IF
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
    write(ifu_log,*)'----------'
    write(ifu_log,*)'Right lead'
    write(ifu_log,*)'----------'
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_b = 0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+hd(1,i,i)
          if(NSpin==2) ro_b=ro_b+hd(2,i,i)
       end do
       ro_a=ro_a/NAOAtom(j)
       ro_b=ro_b/NAOAtom(j)
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,1011)'Atom:',j, ro_a
       if(NSpin==2) write(ifu_log,1011)'Atom:',j, (ro_a+ro_b)/2.0
       IF (Hamilton) THEN
       if(NSpin ==1 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a
       if(NSpin ==2 ) write(ifu_ham,1012)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b)/2.0
       END IF
    end do
1012 format(4f10.5)
1011 format(a6,i4,f12.6)
     !close(ifu_ham)
  end subroutine Hamiltonian 

  !******************************!
  ! Mulliken population analysis !
  !******************************!
  subroutine MullPop
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol
    USE parameters, only: Mulliken, LDOS_Beg, LDOS_End
    use g09Common, only: GetAtmCo
    use constants, only: Bohr
    implicit none

    integer :: i,j, I1, is ,n, l
    real*8 :: sdeg, ro_a, ro_b, chargemol, chargelead1, chargelead2, spinlead1, spinlead2, spinmol
    real*8, dimension(NAOrbs,NAOrbs) :: rho_a, rho_b, tmp  
 
    write(ifu_log,*)'-------------------------------------'
    write(ifu_log,*)'---  Mulliken population analysis ---'
    write(ifu_log,*)'-------------------------------------'
    
    if (NSpin.eq.2) sdeg=1.0d0
    if (NSpin.eq.1) sdeg=2.0d0

    rho_a = matmul( PD(1,:,:), SD )
    if( NSpin == 2 ) rho_b = matmul( PD(2,:,:), SD )

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 1'
    write(ifu_log,*)'----------------------'
    I1=0
    chargelead1=0.0
    spinlead1=0.0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_b=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       if(NSpin ==1 ) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin ==2 ) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead1=chargelead1+(ro_a+ro_b)
       if (NSpin == 2) spinlead1=spinlead1+(ro_a-ro_b)
       I1 = I1 + NAOAtom(j)
    end do
    
   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 1:',chargelead1*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 1:',spinlead1
   write(ifu_log,*)'---------------------------------------------------------'

    chargemol=0.0d0
    spinmol=0.0d0
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'-------------------'
       write(ifu_log,*)'Charges in molecule'
       write(ifu_log,*)'-------------------'
       do j = NALead(1)+1,NALead(1)+NAMol()
          ro_a = 0.0d0
          ro_b = 0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+rho_a(i,i)
             if(NSpin==2) ro_b=ro_b+rho_b(i,i)
          end do
          if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',sdeg*ro_a
          if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
          IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
          END IF
          IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
            if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
          END IF
          chargemol=chargemol+ro_a+ro_b
          if (NSpin == 2) spinmol=spinmol+(ro_a-ro_b)
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in molecule:',chargemol*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in molecule:',spinmol
   write(ifu_log,*)'----------------------------------------------------'

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 2'
    write(ifu_log,*)'----------------------'
    chargelead2=0.0
    spinlead2=0.0
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_b = 0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          if(NSpin==2) ro_b=ro_b+rho_b(i,i)
       end do
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,1011)'Atom:',j,' El.dens:',ro_a*sdeg
       if(NSpin==2) write(ifu_log,1012)'Atom:',j,' El.dens:', (ro_a+ro_b),' Sp.dens:', (ro_a-ro_b)
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,1013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_a-ro_b)
       END IF
       chargelead2=chargelead2+ro_a+ro_b
       if (NSpin == 2) spinlead2=spinlead2+(ro_a-ro_b)
    end do

   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 2:',chargelead2*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 2:',spinlead2
   write(ifu_log,*)'---------------------------------------------------------'

1011 format(a6,i4,a10,f8.4)
1012 format(a6,i4,a10,f8.4,a10,2f9.4)
1013 format(7f9.4)
  end subroutine MullPop
  
  !*******************************************************************!
  ! Mulliken population analysis with non-collinear or SOC hamiltonian!
  !*******************************************************************!
  subroutine MullPop_SOC
    use cluster, only: NALead, NAMol, NAOAtom, NAOMol
    USE parameters, only: Mulliken, LDOS_Beg, LDOS_End
    use g09Common, only: GetAtmCo
    use constants, only: Bohr, d_zero, c_zero
    implicit none

    integer :: i,j, I1, is ,n, l
    real*8 :: sdeg, ro_a, ro_ab, ro_ba, ro_ab_I, ro_ba_I, ro_b, spindens, chargemol, chargelead1, chargelead2, spinlead1, spinlead2, spinmol
    real*8, dimension(NAOrbs,NAOrbs) :: rho_a, rho_ab, rho_ba, rho_ab_I, rho_ba_I, rho_b, tmp, PD_SOC_UU, S_SOC_UU, PD_SOC_UD, PD_SOC_UD_I, S_SOC_UD, PD_SOC_DU, PD_SOC_DU_I, S_SOC_DU, PD_SOC_DD, S_SOC_DD  

    write(ifu_log,*)'-----------------------------------------------'
    write(ifu_log,*)'--- Mulliken population analysis (with SOC) ---'
    write(ifu_log,*)'-----------------------------------------------'

    S_SOC_UU=d_zero
    S_SOC_UD=d_zero
    S_SOC_DU=d_zero
    S_SOC_DD=d_zero
    PD_SOC_UU=d_zero
    PD_SOC_UD=d_zero
    PD_SOC_UD_I=d_zero
    PD_SOC_DU=d_zero
    PD_SOC_DU_I=d_zero
    PD_SOC_DD=d_zero
    rho_a = d_zero
    rho_ab = d_zero
    rho_ab_I = d_zero
    rho_ba = d_zero
    rho_ba_I = d_zero
    rho_b = d_zero

    
    do i=1,NAOrbs
    do j=1,NAOrbs
       S_SOC_UU(i,j)=S_SOC(i,j)
       S_SOC_UD(i,j)=S_SOC(i+NAOrbs,j)
       S_SOC_UD(i,j)=S_SOC(i,j+NAOrbs)
       S_SOC_DD(i,j)=S_SOC(i+NAOrbs,j+NAOrbs)
       PD_SOC_UU(i,j)=REAL(PD_SOC(i,j))
       PD_SOC_UD(i,j)=REAL(PD_SOC(i,j+NAOrbs))
       PD_SOC_UD_I(i,j)=DIMAG(PD_SOC(i,j+NAOrbs))
       PD_SOC_DU(i,j)=REAL(PD_SOC(i+NAOrbs,j))    
       PD_SOC_DU_I(i,j)=DIMAG(PD_SOC(i+NAOrbs,j))                  
       PD_SOC_DD(i,j)=REAL(PD_SOC(i+NAOrbs,j+NAOrbs))
    end do
    end do         
    
    if (NSpin.eq.2) sdeg=1.0d0
    if (NSpin.eq.1) sdeg=2.0d0

    rho_a = matmul( PD_SOC_UU, S_SOC_UU )+matmul( PD_SOC_UD,S_SOC_DU )
    IF( NSpin == 2 ) THEN 
      rho_ab = matmul( PD_SOC_UU, S_SOC_UD )+matmul( PD_SOC_UD, S_SOC_DD )
      rho_ab_I = matmul( PD_SOC_UD_I, S_SOC_DD ) !matmul( PD_SOC_UU, S_SOC_UD )+
      rho_ba = matmul( PD_SOC_DU, S_SOC_UU )+matmul( PD_SOC_DD,S_SOC_DU )
      rho_ba_I = matmul( PD_SOC_DU_I, S_SOC_UU )!+matmul( PD_SOC_DD,S_SOC_DU )
      rho_b = matmul( PD_SOC_DU, S_SOC_UD )+matmul( PD_SOC_DD, S_SOC_DD )
    END IF

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 1'
    write(ifu_log,*)'----------------------'
    I1=0
    chargelead1=0.0
    spinlead1=0.0
    do j=1,NALead(1)
       ro_a=0.0d0
       ro_ab=0.0d0
       ro_ab_I=0.0d0
       ro_ba=0.0d0
       ro_ba_I=0.0d0
       ro_b=0.0d0
       spindens=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          IF(NSpin==2) THEN 
            ro_ab=ro_ab+rho_ab(i,i)
            ro_ab_I=ro_ab_I+rho_ab_I(i,i)
            ro_ba=ro_ba+rho_ba(i,i)  
            ro_ba_I=ro_ba_I+rho_ba_I(i,i)
            ro_b=ro_b+rho_b(i,i)
          END IF 
       end do
       if(NSpin ==1 ) write(ifu_log,2011)'Atom:',j,' El.dens:',ro_a*sdeg
       IF(NSpin ==2 ) THEN     
         spindens = sqrt((ro_ab+ro_ba)**2+(ro_ab_I-ro_ba_I)**2+(ro_a-ro_b)**2)
         write(ifu_log,2012)'Atom:',j,' El.dens:',(ro_a+ro_b),' Sp.dens.x:',(ro_ab+ro_ba),' Sp.dens.y:',(ro_ab_I-ro_ba_I),' Sp.dens.z:',(ro_a-ro_b),' Coll. sp.dens:',spindens
       END IF
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
       END IF
       chargelead1=chargelead1+(ro_a+ro_b)
       if (NSpin == 2) spinlead1=spinlead1+spindens
       I1 = I1 + NAOAtom(j)
    end do
    
   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 1:',chargelead1*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 1:',spinlead1
   write(ifu_log,*)'---------------------------------------------------------'

    chargemol=0.0d0
    spinmol=0.0d0
    if (NAOMol().ge.1) then 
       write(ifu_log,*)'-------------------'
       write(ifu_log,*)'Charges in molecule'
       write(ifu_log,*)'-------------------'
       do j = NALead(1)+1,NALead(1)+NAMol()
       ro_a=0.0d0
       ro_ab=0.0d0
       ro_ab_I=0.0d0
       ro_ba=0.0d0
       ro_ba_I=0.0d0
       ro_b=0.0d0
       spindens=0.0d0
          do i=I1+1,I1+NAOAtom(j)
             ro_a=ro_a+rho_a(i,i)
             IF(NSpin==2) THEN 
                ro_ab=ro_ab+rho_ab(i,i)
                ro_ab_I=ro_ab_I+rho_ab_I(i,i)
                ro_ba=ro_ba+rho_ba(i,i)  
                ro_ba_I=ro_ba_I+rho_ba_I(i,i)
                ro_b=ro_b+rho_b(i,i)
             END IF 
          end do
          if(NSpin==1) write(ifu_log,2011)'Atom:',j,' El.dens:',sdeg*ro_a
          IF(NSpin ==2 ) THEN 
            spindens = sqrt((ro_ab+ro_ba)**2+(ro_ab_I-ro_ba_I)**2+(ro_a-ro_b)**2)
            write(ifu_log,2012)'Atom:',j,' El.dens:',(ro_a+ro_b),' Sp.dens.x:',(ro_ab+ro_ba),' Sp.dens.y:',(ro_ab_I-ro_ba_I),' Sp.dens.z:',(ro_a-ro_b),' Coll. sp.dens:',spindens
          END IF
          IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
            if(NSpin ==2 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
          END IF
          IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
            if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
            if(NSpin ==2 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
          END IF
          chargemol=chargemol+ro_a+ro_b
          if (NSpin == 2) spinmol=spinmol+spindens
          I1 = I1 + NAOAtom(j)
       end do
    end if
    
   write(ifu_log,*)'----------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in molecule:',chargemol*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in molecule:',spinmol
   write(ifu_log,*)'----------------------------------------------------'

    write(ifu_log,*)'----------------------'
    write(ifu_log,*)'Charges in electrode 2'
    write(ifu_log,*)'----------------------'
    chargelead2=0.0
    spinlead2=0.0
    do j=NALead(1)+NAMol()+1,NALead(1)+NAMol()+NALead(2)
       ro_a=0.0d0
       ro_ab=0.0d0
       ro_ab_I=0.0d0
       ro_ba=0.0d0
       ro_ba_I=0.0d0
       ro_b=0.0d0
       spindens=0.0d0
       do i=I1+1,I1+NAOAtom(j)
          ro_a=ro_a+rho_a(i,i)
          IF(NSpin==2) THEN 
            ro_ab=ro_ab+rho_ab(i,i)
            ro_ab_I=ro_ab_I+rho_ab_I(i,i)
            ro_ba=ro_ba+rho_ba(i,i)  
            ro_ba_I=ro_ba_I+rho_ba_I(i,i)
            ro_b=ro_b+rho_b(i,i)
          END IF 
       end do
       I1 = I1 + NAOAtom(j)
       if(NSpin==1) write(ifu_log,2011)'Atom:',j,' El.dens:',ro_a*sdeg
       IF(NSpin ==2 ) THEN    
         spindens = sqrt((ro_ab+ro_ba)**2+(ro_ab_I-ro_ba_I)**2+(ro_a-ro_b)**2)
         write(ifu_log,2012)'Atom:',j,' El.dens:',(ro_a+ro_b),' Sp.dens.x:',(ro_ab+ro_ba),' Sp.dens.y:',(ro_ab_I-ro_ba_I),' Sp.dens.z:',(ro_a-ro_b),' Coll. sp.dens:',spindens
       END IF
       IF(Mulliken .and. LDOS_Beg <= LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg,AtomDOSEF(1,j)
         if(NSpin ==2 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b),(AtomDOSEF(l,j)*(-1)**(l+1),l=1,2)
       END IF
       IF(Mulliken .and. LDOS_Beg > LDOS_End) THEN
         if(NSpin ==1 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),ro_a*sdeg
         if(NSpin ==2 ) write(ifu_mul,2013)(GetAtmCo(n,j)*Bohr,n=1,3),(ro_a+ro_b),(ro_ab+ro_ba),(ro_ab_I-ro_ba_I),(ro_a-ro_b)
       END IF
       chargelead2=chargelead2+ro_a+ro_b
       if (NSpin == 2) spinlead2=spinlead2+spindens
    end do

   write(ifu_log,*)'---------------------------------------------------------'
   write(ifu_log,*)'Total num. of electrons in electrode 2:',chargelead2*sdeg
   if (NSpin == 2) write(ifu_log,*)'Total spin in electrode 2:',spinlead2
   write(ifu_log,*)'---------------------------------------------------------'    

2011 format(a6,i4,a10,f8.4)
2012 format(a6,i4,a10,f8.4,3(a12,f9.4),a16,f9.4)
2013 format(10f9.4)
  end subroutine MullPop_SOC  

!*******************************************************************************
!* Subroutine to compute LDOS(E) when it is not computed by Transmission        *
!*******************************************************************************
  SUBROUTINE LDOS
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: c_one, c_zero, d_zero, d_pi
    use parameters, only: LDOS_Beg, LDOS_End, EW1, EW2, EStep, DOSEnergy
    use preproc, only: MaxAtm
!   USE IFLPORT
    use omp_lib
    use g09Common, only: GetNAtoms

    real*8 :: energy, trans, DOS, energ
    real*8, dimension(10001) :: xxx
    real*8, dimension(MaxAtm) :: AtomDOS
    complex*16 :: cenergy,ctrans
    integer :: n, nsteps, i, imin, imax, info ,j
    
    complex*16, dimension(NAOrbs,NAOrbs) :: GammaL, GammaR, Green, T, temp, SG

    print *
    print *, "-------------------------"
    print *, "--- Calculating  LDOS ---"
    print *, "-------------------------"
    print *

    allocate ( AtomDOSEF(2,MaxAtm) )

    nsteps = (EW2-EW1)/EStep + 1
    do ispin=1,NSpin

       open(333,file='tempDOS',status='unknown')
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,green,gammar,gammal) 
!$OMP DO SCHEDULE(STATIC,10)
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          call gplus(cenergy,Green,GammaR,GammaL)

          ! Mulliken DOS 
!$OMP CRITICAL
          SG = matmul( SD, green )
          ! computing total DOS
          DOS=d_zero
          AtomDOS=d_zero
          do j=1,GetNAtoms()
          do i=LoAOrbNo(j),HiAOrbNo(j)
             AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
             DOS=DOS-dimag(SG(i,i))/d_pi
          end do
          end do

          if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

          ! print out DOS and atomic orbital resolved DOS ***
          imin = LoAOrbNo(LDOS_Beg)
          if( imin < 1 ) imin = 1
          imax = HiAOrbNo(LDOS_End)
          if( imax > NAOrbs ) imax = NAOrbs
          call flush(333)
          write(333,3333) energy,DOS*(-1)**(ispin+1),(AtomDOS(j)*(-1)**(ispin+1),j=LDOS_Beg,LDOS_End),(-dimag(SG(i,i))*(-1)**(ispin+1)/d_pi,i=imin,imax)

!$OMP END CRITICAL

       end do ! End of energy loop
      
!$OMP END DO
!$OMP END PARALLEL

  ! Reordering in energy for nice output
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(333)
          do i=1,10000000000
          read(333,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(333)
             read(333,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             write(ifu_dos,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             exit
          end if
          end do
       end do
 
      close(333,status='delete')

        write(ifu_dos,*) '    '                    
      end do ! End of spin loop


3333 format(f10.5,10000E14.5)

  END SUBROUTINE LDOS

!*******************************
!* Subroutine to evaluate T(E) *
!*******************************
  subroutine transmission
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: c_one, c_zero, d_zero, d_pi
    use parameters, only: NChannels,HTransm,EW1,EW2,EStep,LDOS_Beg,LDOS_End, DOSEnergy, SOC, FermiAcc, QExcess, ChargeAcc
    use numeric, only: CMatPow, CHDiag, CDiag, sort, MULLER
     use lapack_blas, only: zgemm
!   use lapack95, only: zgemm
    use preproc, only: MaxAtm
!   USE IFLPORT
    use omp_lib
    use g09Common, only: GetNAtoms
    implicit none
    external zgemm    

    real*8 :: energy, trans, DOS, energ, E0, E1, E2, E3, Delta, Epsilon, DE
    real*8, dimension(MaxAtm) :: AtomDOS
    real*8, dimension(10001) :: xxx
    complex*16 :: cenergy,ctrans
    integer :: n, nsteps, i, imin, imax, info, j, AllocErr, cond, k
    integer :: Max = 20
    complex*16, dimension(:,:), allocatable :: GammaL, GammaR, Green, T, temp, SG
    complex*16, dimension(:,:), allocatable :: DGammaL, DGammaR, DGreen, DT, Dtemp, DSG
    complex*16, dimension(:,:),allocatable :: dummy
    real*8, dimension(:), allocatable :: tn,Dtn
    complex*16, dimension(:), allocatable :: ctn,Dctn
    real*8, dimension(:),allocatable   :: tchan1,tchan2

    print *
    print *, "--------------------------------"
    print *, "--- Calculating Transmission ---"
    print *, "---  (and DOS if required)   ---"
    print *, "--------------------------------"
    print *

    if (SOC) then
       write(ifu_log,*)' Adding spin-orbit coupling ...'
       write(ifu_log,*)'... and finding new Fermi level'
       allocate(DGammaL(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DGammaR(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DGreen(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DT(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Dtemp(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(DSG(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(S_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(H_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(PD_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(PDOUT_SOC(DNAOrbs,DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Dtn(DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Dctn(DNAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       DSG=c_zero
       DT=c_zero
       Dtemp=c_zero
       S_SOC=d_zero

       call spin_orbit
       
    else
       allocate(GammaL(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(GammaR(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(Green(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(T(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(temp(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(SG(NAOrbs,NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(tn(NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
       allocate(ctn(NAOrbs), STAT=AllocErr)
       if( AllocErr /= 0 ) stop
    end if

    if( NChannels > 0 ) then
       print *, "Number of eigen-channels to print out: ", NChannels 
       allocate( tchan1(NChannels), tchan2(NChannels), dummy(NChannels,NChannels), STAT = AllocErr )
       if( AllocErr /= 0 ) stop
    end if

    allocate( AtomDOSEF(2,MaxAtm), STAT=AllocErr)
    if( AllocErr /= 0 ) stop

    if (SOC) then
      
    ! finding new Fermi energy with SOC

       E0=shift-10.0d0*FermiAcc
       E1=shift
       E2=shift+10.0d0*FermiAcc
       Delta=FermiAcc
       Epsilon=ChargeAcc*(NCDEl+QExcess)
       print*,'MULLER method'
       call MULLER(QXTot_SOC,E0,E1,E2,Delta,Epsilon,Max,E3,DE,K,Cond)
       if(k .eq. Max .or. E2<EMin .or. E2>EMax) then
         write(ifu_log,*) 'I could not accurately find the new Fermi level ...'
         write(ifu_log,*) ' ...using the best approximation'
       end if
       shift = E3
       write(ifu_log,*)'-----------------------------------------------'
       write(ifu_log,'(A,F9.5,A,f9.5)') ' Fermi energy= ',-shift,'  +/-',dabs(DE)
       write(ifu_log,*)
       write(ifu_log,'(A,F10.5)') ' Number of electrons:  ', Q_SOC
       write(ifu_log,*)'-----------------------------------------------'
    end if
         
    nsteps = (EW2-EW1)/EStep + 1

    do ispin=1,NSpin

      if (ispin == 2 .and. SOC) exit

      if (LDOS_Beg <= LDOS_End ) open(333,file='tempDOS',status='unknown')
      open(334,file='tempT',status='unknown')

      if (.not. SOC) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,green,gammar,gammal,T,temp) FIRSTPRIVATE(nsteps)
!$OMP DO SCHEDULE(STATIC,10)
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
          call gplus(cenergy,Green,GammaR,GammaL)

          if( .not. HTransm )then
             !*************************************************************
             !* Here we use the following non-Hermitian expression  for T *
             !* [Gamma_L G^a Gamma_R G^r]                                 *
             !* It works better for large clusters                        *
             !*************************************************************
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one, GammaL,NAorbs, Green,  NAOrbs, c_zero, T,    NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, T,     NAOrbs, GammaR, NAOrbs, c_zero, temp, NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, temp,  NAOrbs, Green,  NAOrbs, c_zero, T,    NAOrbs)
          else
             !********************************************************
             !* Here we use the following Hermitian expression for T *
             !* [Gamma_L^1/2 G^a Gamma_R G^r Gamma_L^1/2]            *
             !********************************************************
             call CMatPow(GammaL,0.5d0,temp)
             GammaL=temp
             call zgemm('N','C',NAOrbs,NAOrbs,NAOrbs,c_one,GammaL,NAOrbs,Green, NAOrbs,c_zero,temp,NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,temp,  NAOrbs,GammaR,NAOrbs,c_zero,T,   NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,T,     NAOrbs,Green, NAOrbs,c_zero,temp,NAOrbs)
             call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one,temp,  NAOrbs,GammaL,NAOrbs,c_zero,T,   NAOrbs)
          end if

!$OMP CRITICAL
          ! Mulliken DOS 
          if (LDOS_Beg <= LDOS_End ) then
            SG = matmul( SD, green )
            ! computing total DOS
            DOS=d_zero
            AtomDOS=d_zero
            do j=1,GetNAtoms()
            do i=LoAOrbNo(j),HiAOrbNo(j)
               AtomDOS(j)=AtomDOS(j)-dimag(SG(i,i))/d_pi
               DOS=DOS-dimag(SG(i,i))/d_pi
            end do
            end do

            if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

            ! print out DOS and atomic orbital resolved DOS ***
            imin = LoAOrbNo(LDOS_Beg)
            if( imin < 1 ) imin = 1
            imax = HiAOrbNo(LDOS_End)
            if( imax > NAOrbs ) imax = NAOrbs
            call flush(333)
            write(333,3333) energy,DOS*(-1)**(ispin+1),(AtomDOS(j)*(-1)**(ispin+1),j=LDOS_Beg,LDOS_End),(-dimag(SG(i,i))*(-1)**(ispin+1)/d_pi,i=imin,imax)
          end if

          ! computing transmission T
          ctrans=c_zero
          do i=1,NAOrbs
             ctrans=ctrans + T(i,i)
          end do

          if (dimag(ctrans).gt.1.0d-5) then
             write(ifu_log,*)'Transmission not real !!!'
             stop
          end if
          trans=ctrans

          ! Diagonalize the T matrix 
          ! to get eigen channels
          if( NChannels > 0 )then
             if( HTransm ) then 
                call CHDiag( T, tn, info )
             else
                call CDiag( T, ctn, info )
                do i=1,NAOrbs
                  tn(i) = dble( ctn(i) )
                end do
                ! sort eigenvalues smallest to biggest
                call sort(NAOrbs,tn)
             end if
             if( n > 3 ) call SeparateSpaghettis( tchan1, tchan2, tn(NAOrbs-NChannels+1:NAOrbs), dummy, NChannels)
             tchan1=tchan2
             tchan2=tn(NAOrbs-NChannels+1:NAOrbs)
          end if
          call flush(334)
          write(334,1002)energy,trans,(tn(i),i=NAOrbs,NAOrbs-NChannels+1,-1)
          
!$OMP END CRITICAL
       end do ! End of energy loop
!$OMP END DO
!$OMP END PARALLEL

       else !SOC case

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,cenergy,energy,Dgreen,Dgammar,Dgammal,DT,Dtemp) 
!$OMP DO SCHEDULE(STATIC,10)
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          cenergy=dcmplx(energy)

          !*********************************************************************
          !* Evaluation of the retarded "Green" function and coupling matrices *
          !*********************************************************************
             call gplus_SOC(cenergy,DGreen,DGammaR,DGammaL)
             call zgemm('N','C',DNAOrbs,DNAOrbs,DNAOrbs,c_one, DGammaL,DNAOrbs, DGreen,  DNAOrbs, c_zero, DT,    DNAOrbs)
             call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one, DT,     DNAOrbs, DGammaR, DNAOrbs, c_zero, Dtemp, DNAOrbs)
             call zgemm('N','N',DNAOrbs,DNAOrbs,DNAOrbs,c_one, Dtemp,  DNAOrbs, DGreen,  DNAOrbs, c_zero, DT,    DNAOrbs)

!$OMP CRITICAL
      ! Mulliken DOS 
           if (LDOS_Beg <= LDOS_End ) then
             DSG = matmul( S_SOC, DGreen )
             DOS=d_zero
             AtomDOS=d_zero
             do j=1,GetNAtoms()
             do i=LoAOrbNo(j),HiAOrbNo(j)
                AtomDOS(j)=AtomDOS(j)-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs))/(2*d_pi)
                DOS=DOS-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs))/(2*d_pi)
             end do
             end do

             if (dabs(energy-DOSEnergy) < EStep/2) AtomDOSEF(ispin,:)=AtomDOS

     ! print out DOS and atomic orbital resolved DOS ***
             imin = LoAOrbNo(LDOS_Beg)
             if( imin < 1 ) imin = 1
             imax = HiAOrbNo(LDOS_End)
             if( imax > NAOrbs ) imax = NAOrbs
             call flush(333)
             write(333,3333) energy,DOS,(AtomDOS(j),j=LDOS_Beg,LDOS_End),((-dimag(DSG(i,i))/(2*d_pi)-dimag(DSG(i+NAOrbs,i+NAOrbs)))/(2*d_pi),i=imin,imax)
           end if

          ! computing transmission T
          ctrans=c_zero
          do i=1,DNAOrbs
             ctrans=ctrans + DT(i,i)
          end do

          if (dimag(ctrans).gt.1.0d-5) then
             write(ifu_log,*)'Transmission not real !!!'
             stop
          end if
          trans=ctrans/2.0   ! in units of 2e^2/h

          ! Diagonalize the T matrix 
          ! to get eigen channels
          if( NChannels > 0 )then
             if( HTransm ) then 
                call CHDiag( DT, Dtn, info )
             else
                call CDiag( DT, Dctn, info )
                do i=1,DNAOrbs
                  Dtn(i) = dble( Dctn(i) )
                end do
                ! sort eigenvalues smallest to biggest
                call sort(DNAOrbs,Dtn)
             end if
             if( n > 3 ) call SeparateSpaghettis( tchan1, tchan2, Dtn(DNAOrbs-NChannels+1:DNAOrbs), dummy, NChannels)
             tchan1=tchan2
             tchan2=tn(DNAOrbs-NChannels+1:DNAOrbs)
          end if

          call flush(334)
          !write(334,1002)energy,trans,(Dtn(i),i=DNAOrbs,DNAOrbs-NChannels+1,-1)
          write(334,1002)energy,trans

!$OMP END CRITICAL
       end do ! End of energy loop
!$OMP END DO
!$OMP END PARALLEL

       end if

       if (LDOS_Beg <= LDOS_End ) then
  ! Reordering in energy for nice DOS output
       do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(333)
          do i=1,10000000000
          read(333,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(333)
             read(333,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             write(ifu_dos,3333) (xxx(j),j=1,2+(LDOS_End-LDOS_Beg+1)+(imax-imin+1))
             exit
          end if
          end do
       end do
      write(ifu_dos,*)'   '
      close(333,status='delete')
      end if

  ! Reordering in energy for nice T output
      do n=1,nsteps
          energy=EW1+EStep*(n-1)
          rewind(334)
          do i=1,10000000000
          read(334,*)energ
          if (dabs(energy-energ) < 0.000001) then
             backspace(334)
             read(334,1002) (xxx(j),j=1,2+NChannels)
             write(ifu_tra,1002) (xxx(j),j=1,2+NChannels)
             exit
          end if
          end do
       end do
      close(334,status='delete')

      write(ifu_tra,*)'   '

      end do ! End of spin loop

      if (SOC) then
         deallocate(DGammaL)
         deallocate(DGammaR)
         deallocate(DGreen)
         deallocate(DT)
         deallocate(Dtemp)
         deallocate(DSG)
         deallocate(Dtn)
         deallocate(Dctn)
         !deallocate(S_SOC)   ! DO NOT deallocate when calculating Mulliken population analysis with S_SOC!!!
         deallocate(H_SOC)
      else 
         deallocate(GammaL)
         deallocate(GammaR)
         deallocate(Green)
         deallocate(T)
         deallocate(temp)
         deallocate(SG)
         deallocate(tn)
         deallocate(ctn)
      end if

    if( NChannels > 0 ) then
       deallocate( tchan1, tchan2, dummy )
    end if

1002 format(f10.5,10000E14.5)
3333 format(f10.5,10000E14.5)

  end subroutine transmission

  !
  ! Computes Hybridization functions for all correlated subspaces
  ! for all energy points defined in mesh.dat and writes them to 
  ! file
  !
  subroutine CompHybFunc
    use constants
    use parameters, only: NCorrBl, CorrBeg, CorrEnd, eta
    use util
    use correlation
    use antcommon
!   USE IFLPORT
    use omp_lib
    implicit none

    integer :: ios, n, nmax, iblock, i, j, iao,jao
    real*8 :: En, Ep
    complex*16 :: zn
    character(len=3) :: istr
    character(len=100) :: fname

    real*8, dimension(:), allocatable :: EMesh
    
    ! Device Green's function
    complex*16, dimension(NAorbs,NAOrbs) :: GD

    complex*16, dimension(:,:,:,:,:), allocatable :: delta

    print *, "-----------------------------------------"
    print *, "--- Computing Hybridization functions ---"
    print *, "-----------------------------------------"

    call SetHamOvl( HD, SD )
    !
    ! Read mesh file mesh.dat
    ! 
    print *, "Reading energy mesh from file mesh.dat" 
    ! open mesh file
    open(unit=ifu_msh,file='mesh.dat',status='old',iostat=ios)
    if( ios /= 0 )then
       print *, "Device/CompHybFunc/Error: Could not open energy mesh file mesh.dat. Abort."
       STOP
    end if
    ios = 0; nmax=0; Ep=-1e+10
    do while( ios == 0 )
       read(unit=ifu_msh,fmt=*,iostat=ios), En
       if( ios /= 0 ) exit
       if( En <= Ep ) exit
       Ep = En
       nmax=nmax+1
       !print *, nmax, En
    end do
    print *, "Mesh file has", nmax, " data points."
    allocate( EMesh( nmax ) )
    rewind(ifu_msh)
    do n=1,nmax
       read(unit=ifu_msh,fmt=*,iostat=ios), EMesh(n)
       if( ios /= 0 ) exit
    end do
    close(ifu_msh)
    !
    ! calculate hybridization function for all mesh points
    ! 
    allocate( delta(nmax,NSpin,NCorrBl,NMaxCorr,NMaxCorr) ) 
    delta = c_zero

    do ispin=1,NSpin
!$OMP PARALLEL PRIVATE(En,zn,GD)
!$OMP DO
       do n=1,nmax
          En = EMesh(n)
          zn = En+ui*eta
          call gplus0(zn,GD)
!$OMP CRITICAL
          call CompDelta( ispin, En, -shift, GD, delta(n,ispin,:,:,:) )
!$OMP END CRITICAL
       end do
!$OMP END DO
!$OMP END PARALLEL
    end do
    !
    ! write hybridization functions to files
    !
    do iblock=1,NCorrBl
       print *, "Correlated block", iblock
       call int2str( iblock, istr )
       fname='delta.out.'//istr
       print *, "Output file for Hybridization:", fname
       open(unit=ifu_hyb,file=fname,status='unknown',iostat=ios)
       fname='Ac.out.'//istr 
       print *, "Output file for Bath Sepctral function:", fname
       open(unit=ifu_ac,file=fname,status='unknown',iostat=ios)
       do n=1,nmax
          En=EMesh(n)
          !
          ! write hybridization function Delta
          !
          write(ifu_hyb,fmt='(E20.10,1000E14.6)'),& 
               En, ( (delta(n,ispin,iblock,i,i),ispin=1,NDSpin),i=1,ncorrao(iblock) )
          call flush(ifu_hyb)
          !
          ! write bath spectral function = -Im(Delta)/pi
          !
          write(ifu_ac,fmt='(E20.10,1000E14.6)'), &
               En, ( (-AIMAG(delta(n,ispin,iblock,i,i))/d_pi,ispin=1,NDSpin),i=1,ncorrao(iblock) )
          call flush(ifu_ac)
       end do
       close(ifu_hyb)
       close(ifu_ac)
    end do
    print *, "done." 

    deallocate( EMesh, delta )
    
  end subroutine CompHybFunc


  !*******************************************!
  ! Routine for orbital eigen-channel         !
  ! analysis with reduced transmission matrix !
  !*******************************************!
  subroutine EigenChannelAnalysis
    use parameters, only:  RedTransmB,RedTransmE, eta, EW1, EW2,EStep
    use Cluster, only : hiaorbno, loaorbno
    use constants, only: ui, d_pi
    use numeric, only: CInv, CMatPow, CHDiag
    implicit none
    
    real*8 :: rho, phi,DE,energy,delta,mindelta,tsave
    real*8, parameter :: dsmall = 1.0d-10
    integer :: NMaxEV = 10 ! At how many points to print out eign vectors: 2*NMaxEV+1

    integer :: NSD, NLD, NRD, NS1, NS2, N,NMax, info

    complex*16,dimension(NAOrbs,NAOrbs) :: SigmaL, SigmaR
    complex*16,dimension(:,:),allocatable :: GSD, SigmaLD, SigmaRD, GammaLD, GammaRD, GammaLDph, TS, TS0
    complex*16,dimension(:,:),allocatable :: gLD, gRD
    complex*16,dimension(:,:),allocatable :: VLS, VRS

    complex*16 :: cenergy, ci
    integer :: is, i, nchan, j, jmin, k

    real*8,dimension(:),allocatable :: tchan,tchan1,tchan2

    ! Dimension of scattering region inside device
    NSD = HiAOrbNo(RedTransmE)-LoAOrbNo(RedTransmB)+1
    ! Dimension of left region of device
    NLD = LoAOrbNo(RedTransmB)-1
    ! Dimension of right region of device
    NRD = NAOrbs-HiAOrbNo(RedTransmE)
    
    NS1 = LoAOrbNo(RedTransmB)
    NS2 = HiAOrbNo(RedTransmE)

    allocate( GSD(NSD,NSD), SigmaLD(NSD,NSD), SigmaRD(NSD,NSD), &
         GammaLD(NSD,NSD), GammaRD(NSD,NSD), GammaLDph(NSD,NSD), TS(NSD,NSD), TS0(NSD,NSD), &
         gLD(NLD,NLD), gRD(NRD,NRD), VLS(NLD,NSD), VRS(NRD,NSD), tchan(NSD), tchan1(NSD), tchan2(NSD) )

    !open(ifu_red,file='t.dat',status='unknown')

    print *
    print *, "-----------------------------------------------------------------------"
    print *, "--- Orbital Eigen-channel Analysis with reduced Transmission matrix ---"
    print *, "-----------------------------------------------------------------------"
    print *
    print *, "Begin of scattering region: AO ", NS1
    print *, "End of scattering region:   AO ", NS2
    print *

    do is=1,NSpin
       if(NSpin ==2 .and. is==1)print*,"Spin UP"
       if(NSpin ==2 .and. is==2)print*,"Spin DOWN"       

       NMax = int(EW2/EStep)

       !DE = d_zero
       !IF(NMaxE>0) DE = EWindow/DBLE(NMaxE)

       do N = -NMax,NMax!-NMaxE,NMaxE  
          energy = N*EStep
          cenergy = energy

          call CompSelfEnergies(is,cenergy,SigmaL,SigmaR)
       
          do i=1,NLD
             do j=1,NLD
                gLD(i,j) = (energy-shift+ui*eta)*SD(i,j) - HD(is,i,j) - SigmaL(i,j)
             end do
          end do

          do i=1,NRD
             do j=1,NRD
                gRD(i,j) = (energy-shift+ui*eta)*SD(NS2+i,NS2+j) - HD(is,NS2+i,NS2+j) - SigmaR(NS2+i,NS2+j)
             end do
          end do

          info = CInv( gLD )
          info = CInv( gRD )
          
          VLS = HD(is,1:NLD, NS1:NS2) - (energy-shift)*SD(1:NLD, NS1:NS2)
          VRS = HD(is,NS2+1:NAOrbs,NS1:NS2) - (energy-shift)*SD(NS2+1:NAOrbs,NS1:NS2)
          
          SigmaLD = matmul( conjg(transpose(VLS)), matmul( gLD, VLS ) )
          SigmaRD = matmul( conjg(transpose(VRS)), matmul( gRD, VRS ) )
          
          GSD = (energy-shift+ui*eta)*SD(NS1:NS2,NS1:NS2) & 
               - HD(is,NS1:NS2,NS1:NS2) - SigmaLD - SigmaRD

          info = CInv( GSD )

          GammaLD = ui*(SigmaLD-conjg(transpose(SigmaLD)))
          GammaRD = ui*(SigmaRD-conjg(transpose(SigmaRD)))

          call CMatPow( GammaLD, 0.5d0, GammaLDph )
          
          ![Gamma_L^1/2 G^a Gamma_R G^r Gamma_L^1/2] 
          TS = matmul( GammaLDph, conjg(transpose(GSD)) )
          TS = matmul( TS, GammaRD )
          TS = matmul( TS, GSD )
          TS = matmul( TS, GammaLDph )

          call CHDiag( TS, tchan, info )
          if( info /= 0 )then
             print*, "Error diagonalizing Reduced transmission matrix: info=", info
             return
          end if

          ! *** Ordering of eigenchannels ***
          if( N >= -NMax+2 ) call SeparateSpaghettis( tchan1, tchan2, tchan, TS, NSD )
          tchan1=tchan2
          tchan2=tchan

          write(ifu_red,'(F10.5,100E14.5)'), energy,(tchan(i),i=NSD,1,-1)

          if( N ==0 )then 
             print *, "Eigenchannel composition at Fermi level:"
             print *, "----------------------------------------"
             do nchan=NSD,1,-1
                print '(A,I2,A,F9.3)',"Channel ", NSD-nchan+1, ": Transmission = ", tchan(nchan)
                do i=1,NSD
                   ci = TS(i,nchan)
                   rho = abs(ci)
                   if( abs(real(ci)) < dsmall * abs(aimag(ci)) )then
                      phi = sign(0.5*d_pi,aimag(ci))
                   else
                      phi = atan( aimag(ci)/real(ci) )
                      if( real(ci) < 0 .and. aimag(ci) > 0 ) phi = phi + 0.5*d_pi
                      if( real(ci) < 0 .and. aimag(ci) < 0 ) phi = phi - 0.5*d_pi
                   end if
                   print '(I4,A,F8.4,F8.4)', i, " : ", rho, phi
                end do
             end do
          end if

       end do
       write(ifu_red,*)
       print *, " "
    end do

    deallocate( GSD, SigmaLD, SigmaRD, GammaLD, GammaRD, GammaLDph, TS, TS0, gLD, gRD, VLS, VRS , tchan, tchan1, tchan2 )

    !close(ifu_red)

  end subroutine EIGENCHANNELANALYSIS


  !
  ! *** Subroutine to separate individual eigen channel     ***
  ! *** transmissions (= spaghettis) using first derivative ***
  !
  subroutine SeparateSpaghettis( evals1, evals2, evals, evecs, N )
    implicit none
    
    ! eigen values at before last, last and actual energy
    real*8, dimension(N) :: evals1, evals2, evals
    ! eigenvectors at actual energy
    complex*16, dimension(N,N) :: evecs
    ! number of eigenvectors
    integer :: N
    
    integer :: i, j, jmin
    real*8 :: yex, delta, mindelta, valsave

    complex*16, dimension(N) :: vecsave

    do i=1,N

       ! extrapolate actual eigenvalue from last 2 eigenvalues
       yex = 2.0d0*evals2(i)-evals1(i)

       ! Find actual eigenvalue which deviates minimally 
       ! from extrapolated value 
       mindelta = abs(yex-evals(i))
       jmin = i
       do j=i+1,N
          delta =  abs(yex-evals(j))
          if( delta < mindelta )then 
             jmin = j 
             mindelta = delta
          end if
       end do

       ! change eigenvector i with eigenvector jmin
       if( jmin /= i )then 
          vecsave = evecs(:,jmin)
          evecs(:,jmin)=evecs(:,i)
          evecs(:,i)=vecsave
          valsave = evals(jmin)
          evals(jmin) = evals(i)
          evals(i) = valsave
       end if
    end do
    
  end subroutine SeparateSpaghettis
  

  ! ***********************************
  ! Compute Self energies of electrodes
  ! ***********************************
  subroutine CompSelfEnergies( spin, cenergy, Sigma1, Sigma2 )
    use parameters, only: ElType, DD, UD, DU, Overlap
    use BetheLattice, only: CompSelfEnergyBL, LeadBL 
    use OneDLead, only: CompSelfEnergy1D, Lead1D
    use constants
     use lapack_blas, only: zgemm
!   use lapack95, only: zgemm
    implicit none
    external zgemm
    
    integer, intent(in) :: spin
    complex*16, intent(in) :: cenergy
    complex*16, dimension(:,:),intent(inout) :: Sigma1
    complex*16, dimension(:,:),intent(inout) :: Sigma2
    complex*16, dimension(NAOrbs,NAOrbs) :: temp
    integer :: is,omp_get_thread_num
    
    !Sigma1=(0.0,0.0)
    !Sigma2=(0.0,0.0)

    is = spin
    if( DD .and. spin == 1 ) is=2
    if( DD .and. spin == 2 ) is=1
    if( UD .and. spin == 1 ) is=1
    if( UD .and. spin == 2 ) is=2
    if( DU .and. spin == 1 ) is=2
    if( DU .and. spin == 2 ) is=1
    !write(ifu_log,*)omp_get_thread_num(),'in CompSelfEnergies',shift,cenergy,'lead 1'
    select case( ElType(1) )
    case( "BETHE" ) 
       call CompSelfEnergyBL( LeadBL(1), is, cenergy, Sigma1 )
       ! transform to non-orthogonal basis
       if( Overlap < -0.01 .and. .not. HDOrtho )then
          ! Sigma1 -> S^1/2 * Sigma1 * S^1/2
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, Sigma1,NAorbs, SPH,  NAOrbs, c_zero, temp,   NAOrbs)
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, SPH,   NAorbs, temp, NAOrbs, c_zero, Sigma1, NAOrbs)
       endif
    case( "1DLEAD" )
       call CompSelfEnergy1D( Lead1D(1), is, cenergy, Sigma1 )
    case( "GHOST" )
        continue
    end select

    if( DD .and. spin == 1 ) is=2
    if( DD .and. spin == 2 ) is=1
    if( UD .and. spin == 1 ) is=2
    if( UD .and. spin == 2 ) is=1
    if( DU .and. spin == 1 ) is=1
    if( DU .and. spin == 2 ) is=2
    
    !write(ifu_log,*)omp_get_thread_num(),'in CompSelfEnergies',shift,cenergy,'lead 2'
    select case( ElType(2) )
    case( "BETHE" )
       call CompSelfEnergyBL( LeadBL(2), is, cenergy, Sigma2 )
       ! transform to non-orthogonal basis 
       if( Overlap < -0.01 .and. .not. HDOrtho )then
          ! Sigma2 -> S^1/2 * Sigma2 * S^1/2
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, Sigma2,NAorbs, SPH,  NAOrbs, c_zero, temp,   NAOrbs)
          call zgemm('N','N',NAOrbs,NAOrbs,NAOrbs,c_one, SPH,   NAorbs, temp, NAOrbs, c_zero, Sigma2, NAOrbs)
       endif
    case( "1DLEAD" )
       call CompSelfEnergy1D( Lead1D(2), is, cenergy, Sigma2 )
    case( "GHOST" )
        continue
    end select
  end subroutine CompSelfEnergies


  !
  ! *** Integrand for charge integration on complex contour ***
  !
  real*8 function DDOS( E0, R, phi )
    use constants, only: c_zero, ui, d_pi

    real*8, intent(in) :: phi, E0, R
    integer :: i, j !!, ispin 
    complex*16,dimension(NAOrbs,NAOrbs) :: green,gammar,gammal
    complex*16 :: TrGS, z

    z = E0 - R*(cos(phi) - ui*sin(phi)) 

    TrGS=c_zero
    do ispin=1,NSpin
       call gplus(z,green,gammar,gammal)
       !do i=1,NAOrbs
       do i=NCDAO1,NCDAO2
          do j=1,NAOrbs
             TrGS = TrGS + green(i,j)*SD(j,i)
          end do
       end do
    end do
    ! Account for spin degeneracy
    if(NSpin==1) TrGS = TrGS * 2.0d0
    DDOS = -DIMAG(R*(sin(phi)+ui*cos(phi))*TrGS)/d_pi 
  end function DDOS

  ! 
  ! *** Total charge up to energy E, lower bound is EMin ***
  ! 
  real*8 function TotCharge( E )
    use constants, only: d_zero, d_pi
    use parameters, only: ChargeAcc
    use numeric, only: gauleg 
    implicit none
    
    real*8, intent(in) :: E

    integer, parameter :: nmax = 2047
    
    real*8 :: q,qq, E0, R, w_j, phi_j
    integer :: n, n1, n2, j, i
    
    real*8, dimension(nmax) :: x, w
    
    ! Integration contour parameters:
    E0 = 0.5*(E + EMin)
    R  = 0.5*(E - EMin)
    ! Computing integral of DOS over 
    ! complex contour using Gauss-Legendre 
    ! quadrature
    n=1
    do  
       q = d_zero
       do j=1,n
          call gauleg(d_zero,d_pi,x(1:n),w(1:n),n)
          q = q + w(j)*DDOS( E0, R, x(j) )
       end do
       if( n > 1 .and. (q == d_zero .or. abs(q-qq) < ChargeAcc*NCDEl ) ) exit  
       n=2*n+1
       if( n > nmax )then
          print *, "TotCharge/gaussian quadrature has not converged after", nmax, " steps."
          TotCharge = 2.0d0*(NCDAO2-NCDAO1+1) - 10.0d0*ChargeAcc*NCDEl
          return
       end if
       qq = q
    end do
    !print *, "gaussian quadrature converged after", n, " steps. Error:", abs(q-qq)

    TotCharge = q
  end function TotCharge

  ! *************************************
  !  Estimates upper/lower energy 
  !  boundary EMin/EMax,
  !  above/below which DOS is gauranteed 
  !  to be zero
  ! *************************************
  subroutine FindEnergyBounds
    use parameters, only: ChargeAcc,eta

    integer :: i, cond,k 
    real*8 :: EStep, Q
    
    print *, "Searching energy boundaries [EMin, EMax] such that"
    print '(A,I4)', " Total Integrated Spectral Weight (TISW)=", 2*(NCDAO2-NCDAO1+1)

    EStep = 10.0d0 +10000.0 *eta

    do
       Q = TotCharge( EMax )
       print '(A,F15.3,A,F15.3,A,F12.5)', " EMin=", EMin, "  EMax=", EMax , "  TISW=", Q
       if( abs(Q - 2.0d0*(NCDAO2-NCDAO1+1)) < ChargeAcc*NCDEl*10.0 ) then
          exit
       end if
       EMin = EMin - EStep
       EMax = EMax + EStep
    end do
    print *, "--------------------------------------------------"

  end subroutine FindEnergyBounds

!ccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccc
  subroutine IntRealAxis(Er,El,M)

    use constants, only: ui,d_pi,d_zero
    use parameters, only: PAcc 
!   USE IFLPORT
    use omp_lib
   
    real*8,intent(in) :: Er,El
    integer,intent(inout) :: M
    real*8, dimension(M) :: xs,xcc
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(NAOrbs,NAOrbs) :: greenp,greenm 
    complex*16, dimension(NAOrbs,NAOrbs) :: p,q
    complex*16, dimension(NAOrbs,NAOrbs) :: PDP
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real*8 :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

      pi=d_pi

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
      call glesser(E0,green)
      do i=1,NAOrbs
       do j=1,NAOrbs
        PDOUT(ispin,i,j) = -ui*green(i,j)/(2*pi)
        p(i,j) = PDOUT(ispin,i,j)
       enddo
      enddo
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      do i=1,NAOrbs
       do j=1,NAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*PDOUT(ispin,i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       call glesser(Em,greenm)
!!$OMP  SECTION
       call glesser(Ep,greenp)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          do i=1,NAOrbs
           do j=1,NAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
      enddo
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDOUT(ispin,i,j)=PDOUT(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      do i=1,NAOrbs
       do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl.and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      do i=1,NAOrbs
      do j=1,NAOrbs
        rCHp=dble(PDOUT(ispin,i,j)-p(i,j))
        aCHp=dimag(PDOUT(ispin,i,j)-p(i,j))
        rCHq=dble(PDOUT(ispin,i,j)-q(i,j))
        aCHq=dimag(PDOUT(ispin,i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl) M = 1
        PDOUT(ispin,i,j) = 16*PDOUT(ispin,i,j)/(3*(n+1))
        PDOUT(ispin,i,j) = PDOUT(ispin,i,j)*(El-Er)/2
      enddo
      enddo
      write(ifu_log,'(A51,i4)')' Integration of PD_OUT has needed ',(((n-1)/2)+1)/2, ' integration points along real axis'

      return
    end subroutine IntRealAxis

!ccccccccccccccccccccccccccccccc
!    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
!    second kind                                                               c
!        eps: Tolerance                                                        c
!        M:   On input, maximum number of points allowed                       c
!             On output, 0 for an alleged successfull calculation, 1 otherwise c
!        F(): External function to be integrated.                              c
!        CH:  The value of the integral. Interval [-1,1]                       c
!ccccccccccccccccccccccccccccccc
  subroutine IntRealAxis_SOC(Er,El,M)

    use constants, only: ui,d_pi,d_zero,c_zero
    use parameters, only: PAcc 
!   USE IFLPORT
    use omp_lib
   
    real*8,intent(in) :: Er,El
    integer,intent(inout) :: M
    real*8, dimension(M) :: xs,xcc
    complex*16, dimension(DNAOrbs,DNAOrbs) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs) :: greenp,greenm 
    complex*16, dimension(DNAOrbs,DNAOrbs) :: p,q
    complex*16, dimension(DNAOrbs,DNAOrbs) :: PDP
    complex*16 :: E0,Em,Ep
    integer :: n,i,j,l,k,k1
    real*8 :: pi,S0,c0,rchp,rchq,xp,c1,s1,s,cc,x,xx,achp,achq

      pi=d_pi
       
      PDOUT_SOC=c_zero

! Initializing M, n, S0, C0, CH and p

      M = (M-1)*0.5d0
      n = 1
      S0=1
      C0=0
      E0=edex3(El,Er,d_zero)
      call glesser_SOC(E0,green)
      do i=1,DNAOrbs
       do j=1,DNAOrbs
        PDOUT_SOC(i,j) = -ui*green(i,j)/(2*pi)
        p(i,j) = PDOUT_SOC(i,j)
       enddo
      enddo
! Computing the (2n+1) points quadrature formula ...
! ... updating q, p, C1, S1, C0, S0, s and c
1     continue
      do i=1,DNAOrbs
       do j=1,DNAOrbs
         q(i,j) = 2*p(i,j)
         p(i,j) = 2*PDOUT_SOC(i,j)
       enddo
      enddo
      C1 = C0
      S1 = S0
      C0 = sqrt((1+C1)*0.5d0)
      S0 = S1/(2*C0)
      !s = S0
      !cc = C0
      xs(1) = S0
      xcc(1) = C0
      do l=1,n,2
         xs(l+2)=xs(l)*C1+xcc(l)*S1
         xcc(l+2)=xcc(l)*C1-xs(l)*S1
      end do
! ... computing F() at the new points
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,xx,Em,Ep,greenp,greenm,i,j,pdp)
       PDP=d_zero
!$OMP DO SCHEDULE(STATIC,1)
      do l=1,n,2
        xx = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
        Em=edex3(El,Er,-xx)
        Ep=edex3(El,Er,xx)
!!$OMP  PARALLEL DEFAULT(SHARED)
!!$OMP  SECTIONS
!!$OMP  SECTION
       call glesser_SOC(Em,greenm)
!!$OMP  SECTION
       call glesser_SOC(Ep,greenp)
!!$OMP  END SECTIONS
!!$OMP  END PARALLEL
          do i=1,DNAOrbs
           do j=1,DNAOrbs
            pdp(i,j) = pdp(i,j)-ui*(greenm(i,j)+greenp(i,j))*xs(l)**4/(2*pi)
           enddo
          enddo
      enddo
!$OMP END DO
!$OMP CRITICAL
       do i = 1,DNAOrbs
          do j = 1,DNAOrbs
             PDOUT_SOC(i,j)=PDOUT_SOC(i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL

! ... replacing n by 2n+1
         n = n + n + 1
! Stopping?
      do i=1,DNAOrbs
       do j=1,DNAOrbs
        rCHp=dble(PDOUT_SOC(i,j)-p(i,j))
        aCHp=dimag(PDOUT_SOC(i,j)-p(i,j))
        rCHq=dble(PDOUT_SOC(i,j)-q(i,j))
        aCHq=dimag(PDOUT_SOC(i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl.and.n.le.M) goto 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl.and.n.le.M) goto 1
       enddo
      enddo
! Test for successfullness and integral final value
      M = 0
      do i=1,DNAOrbs
      do j=1,DNAOrbs
        rCHp=dble(PDOUT_SOC(i,j)-p(i,j))
        aCHp=dimag(PDOUT_SOC(i,j)-p(i,j))
        rCHq=dble(PDOUT_SOC(i,j)-q(i,j))
        aCHq=dimag(PDOUT_SOC(i,j)-q(i,j))
        if (rCHp*rCHp*16.gt.3*(n+1)*abs(rCHq)*PAcc*NCDEl) M = 1
        if (aCHp*aCHp*16.gt.3*(n+1)*abs(aCHq)*PAcc*NCDEl) M = 1
        PDOUT_SOC(i,j) = 16*PDOUT_SOC(i,j)/(3*(n+1))
        PDOUT_SOC(i,j) = PDOUT_SOC(i,j)*(El-Er)/2
      enddo
      enddo
      write(ifu_log,'(A51,i4)')' Integration of PD_SOC has needed ',(((n-1)/2)+1)/2, ' integration points along real axis'

      return
    end subroutine IntRealAxis_SOC


  !ccccccccccccccccccccccccccccccc
  !c    Numerical integration with the GAUSS-CHEBYSHEV quadrature formula of the  c
  !c second kind                                                                  c
  !c        eps: Tolerance                                                        c
  !c        b: parameter of the change of variable                                c
  !c        Em: maximum value of the energy range                                 c
  !c        M:   On input, maximum number of points allowed                       c
  !c             On output, 0 for an alleged successfull calculation, 1 otherwise c
  !c        dn:  On output, the density matrix                                    c
  !c        CH:  On output, the value of the integral (charge density).           c
  !c             Interval [-1,1]                                                  c
  !c        Eq:  On output, the value of the upper bound of the integral.         c
  !c        The rest of arguments are neeed by the subrtn. gplus                  c
  !ccccccccccccccccccccccccccccccc
  subroutine IntCompPlane(rrr,bi,Emi,M,Eq)
    use parameters, only: PAcc 
    use constants, only: d_pi, d_zero, ui
!   USE IFLPORT
    use omp_lib

    implicit none

    real*8,intent(in) :: rrr, bi, Emi, Eq
    integer,intent(inout) :: M
    real*8, dimension(NAOrbs,NAOrbs) :: PDP

    real*8 :: a,b,Em,S0,c0,x0,er0,der0,ch,xp,q,c1,s1,s,cc,x,erp,erm
    integer :: n,i,j,l,k,k1,chunk!,omp_get_thread_num,omp_get_num_threads
    real*8, dimension(M) :: xs,xcc

    complex*16 :: E0,E
    complex*16, dimension(2) :: EE
    complex*16, dimension(NAOrbs,NAOrbs) :: green
    complex*16, dimension(2,NAOrbs,NAOrbs) :: greenn 

   !logical :: omp_get_nested
 
    a = 1.d0/d_pi
    b = bi
    Em = Emi
    PD(ispin,:,:) = d_zero
    M = (M-1)*0.5
    n = 1
    S0 = 1
    C0 = 0
    x0 = 0.d0
    er0 = edex3(Em,b,x0)
    der0 = 0.5d0*(Em-b)
    E0 = rrr*exp(ui*er0)-rrr+Eq
    call gplus0(E0,green)
    CH = 0.d0
    do i = 1,NAOrbs
       do j =1,NAOrbs
          PD(ispin,i,j)= a*dimag(ui*rrr*exp(ui*er0)*green(i,j))*der0
          CH = CH + PD(ispin,i,j)*SD(j,i)
       enddo
    enddo

    xp = CH
1   q = xp + xp
    xp = CH + CH
    C1 = C0
    S1 = S0
    C0 = sqrt((1+C1)*0.5d0)
    S0 = S1/(C0+C0)
    xs(1) = S0
    xcc(1) = C0
    do l=1,n,2
       xs(l+2)=xs(l)*C1+xcc(l)*S1
       xcc(l+2)=xcc(l)*C1-xs(l)*S1
    end do
    !call omp_set_nested(.true.)
    !call omp_set_num_threads(2)
    !print *, omp_get_nested()
    !print *,'--------------------------' 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l,x,erp,erm,EE,greenn,i,j,pdp)
    PDP=d_zero
     chunk=max(((n+1)/2)/omp_get_num_threads(),1)
    !chunk=1
!$OMP DO SCHEDULE(STATIC,chunk)
    do l=1,n,2
       !write(ifu_log,*)'thread',omp_get_thread_num(),'l=',l
       x = 1+0.21220659078919378103*xs(l)*xcc(l)*(3+2*xs(l)*xs(l))-dble(l)/(n+1)
       erp = 0.5d0*((Em-b)*x + (Em+b))
       erm = 0.5d0*(-(Em-b)*x + (Em+b))
       EE(1) = rrr*exp(ui*erp)-rrr+Eq
       EE(2) = rrr*exp(ui*erm)-rrr+Eq
! Useful in case of nesting is allowed
!$OMP  PARALLEL SHARED(l) PRIVATE(k)
!$OMP  DO SCHEDULE(DYNAMIC,1)
       do k=1,2
          !print *,'l',l,'k',k,omp_get_thread_num()
          call gplus0(EE(k),greenn(k,:,:))
       end do
!$OMP  END DO
!$OMP  END PARALLEL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PDP(i,j) = PDP(i,j)+ a*(dimag(ui*rrr*exp(ui*erp)*greenn(1,i,j))*der0 &
                  &   +dimag(ui*rrr*exp(ui*erm)*greenn(2,i,j))*der0)*xs(l)**4
          end do
       end do
    end do
!$OMP END DO
!$OMP CRITICAL
       do i = 1,NAOrbs
          do j = 1,NAOrbs
             PD(ispin,i,j)=PD(ispin,i,j)+PDP(i,j)
          end do
       end do
!$OMP END CRITICAL
!$OMP END PARALLEL
    CH = 0.d0
    do k=1,NAOrbs
       ! Non-orthogonal basis: Ch = Tr[P*SD]   
       do k1=1,NAOrbs
          CH = CH + PD(ispin,k,k1)*SD(k1,k)
       end do
    enddo
    ! ... replacing n by 2n+1
    n = n + n + 1
    ! Stopping?
    !if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl.and.n.le.M) goto 1
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*10.0d-10.and.n.le.M) goto 1
    ! Test for successfullness and integral final value
    M = 0
    !if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*PAcc*NCDEl) M = 1
    if ((CH-xp)*(CH-xp)*16.gt.3*(n+1)*abs(CH-q)*10.0d-10) M = 1
    CH = 16*CH/(3*(n+1))
    do k=1,NAOrbs
       do l=1,NAOrbs
          PD(ispin,k,l) = 16*PD(ispin,k,l)/(3*(n+1))
       enddo
    enddo
    write(ifu_log,'(A47,I4)')' Integration of P has needed a max no. of loops=',(((n-1)/2)+1)/2

    return
  end subroutine IntCompPlane


  !*********************************************************************
  !* Compute retarded Green's function with SOC (doubling the matrices)*
  !*********************************************************************
  subroutine gplus_SOC(z,green,gammar,gammal)
    use parameters, only: eta, glue
    use constants, only: c_zero, ui
     use lapack_blas, only: zgetri, zgetrf
!   use lapack95, only: zgetri, zgetrf

    implicit none
    external zgetri, zgetrf     

    integer :: i, j, info, omp_get_thread_num
    integer, dimension(DNAOrbs) :: ipiv
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: gammar,gammal
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    
    ! Initilization 
    green=c_zero
    gammar=c_zero
    gammal=c_zero
    sigmar=c_zero
    sigmal=c_zero


    sigr1=-ui*eta*SD 
    sigl1=-ui*eta*SD 
    sigr2=-ui*eta*SD 
    sigl2=-ui*eta*SD 

    call CompSelfEnergies( 1, z, sigl1, sigr1 )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2 )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

    gammar=ui*(sigmar-conjg(transpose(sigmar)))
    gammal=ui*(sigmal-conjg(transpose(sigmal)))

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)

    !write(ifu_log,*)omp_get_thread_num(),green(1,1)
  end subroutine gplus_SOC

  !****************************************************
  !* Compute retarded Green's function with SOC*
  !****************************************************
  subroutine gplus0_SOC(z,green)
    use parameters, only: eta, glue
    use constants, only: c_zero, ui
     use lapack_blas, only: zgetri, zgetrf
!   use lapack95, only: zgetri, zgetrf

    implicit none
    external zgetri, zgetrf 
    
    integer :: i, j, info, omp_get_thread_num
    integer, dimension(DNAOrbs) :: ipiv
    complex*16, dimension(4*DNAOrbs) :: work
    complex*16, intent(in) :: z 
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl1,sigr1
    complex*16, dimension(NAOrbs,NAOrbs) :: sigl2,sigr2
    complex*16, dimension(DNAOrbs,DNAOrbs), intent(out) :: green
    complex*16, dimension(DNAOrbs,DNAOrbs) :: sigmar,sigmal
    
    ! Initilization 
    green=c_zero
    sigmar=c_zero
    sigmal=c_zero
    sigr1=-ui*eta*SD 
    sigl1=-ui*eta*SD 
    sigr2=-ui*eta*SD 
    sigl2=-ui*eta*SD 

    call CompSelfEnergies( 1, z, sigl1, sigr1 )
   
    if (NSpin == 2) call CompSelfEnergies( 2, z, sigl2, sigr2 )
    sigr1=glue*sigr1
    sigl1=glue*sigl1
    sigr2=glue*sigr2
    sigl2=glue*sigl2
    
    !************************************************************************
    !c Coupling matrices                                   
    !************************************************************************
    if (NSpin == 2) then
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr2(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl2(i,j)
       end do
       end do
    else
       do i=1,NAOrbs
       do j=1,NAOrbs
          sigmar(i,j)=sigr1(i,j)
          sigmar(i+NAOrbs,j+NAOrbs)=sigr1(i,j)
          sigmal(i,j)=sigl1(i,j)
          sigmal(i+NAOrbs,j+NAOrbs)=sigl1(i,j)
       end do
       end do
    end if

    !************************************************************************
    !c Retarded "Green" function
    !************************************************************************

    do i=1,DNAOrbs
       do j=1,DNAOrbs
          green(i,j)=(z-shift)*S_SOC(i,j)-H_SOC(i,j)-sigmal(i,j)-sigmar(i,j)
       enddo
    enddo

    call zgetrf(DNAOrbs,DNAOrbs,green,DNAOrbs,ipiv,info)
    call zgetri(DNAOrbs,green,DNAOrbs,ipiv,work,4*DNAOrbs,info)

    !write(ifu_log,*)omp_get_thread_num(),green(1,1)
  end subroutine gplus0_SOC


!ccccccccccccccccccccccccccccccc
! SOC subroutine
!ccccccccccccccccccccccccccccccc

 subroutine spin_orbit

    use SpinOrbit, only: CompHSO
    use cluster, only: LoAOrbNo, HiAOrbNo
    use parameters, only: PrtHatom
    use g09Common, only: GetNShell, GetNAtoms
      
    complex*16, dimension(DNAOrbs,DNAOrbs) :: hamil, hamil_SO, overlap_SO
    integer :: i,j,totdim,nshell,Atom
 
 Atom = PrtHatom   
 totdim=NAOrbs   
 hamil = dcmplx(0.0d0,0.0d0)
 hamil_SO = dcmplx(0.0d0,0.0d0)
 overlap_SO = dcmplx(0.0d0,0.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Duplicate the size of the Hamiltonian and Overlap matrix to include up and down
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (NSpin == 2) then
    do i=1,NAOrbs
    do j=1,NAOrbs
       hamil(i,j)=dcmplx(HD(1,i,j),0.0d0)
       hamil(i+NAOrbs,j+NAOrbs)=dcmplx(HD(2,i,j),0.0d0)
       overlap_SO(i,j)=dcmplx(SD(i,j),0.0d0)              
       overlap_SO(i+NAOrbs,j+NAOrbs)=dcmplx(SD(i,j),0.0d0)
    end do
    end do
 else 
    do i=1,NAOrbs
    do j=1,NAOrbs
       hamil(i,j)=dcmplx(HD(1,i,j),0.0d0)
       hamil(i+NAOrbs,j+NAOrbs)=dcmplx(HD(1,i,j),0.0d0)
       overlap_SO(i,j)=dcmplx(SD(i,j),0.0d0)              
       overlap_SO(i+NAOrbs,j+NAOrbs)=dcmplx(SD(i,j),0.0d0)
    end do
    end do
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Return the input matrix with double size plus the soc interaction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 nshell = GetNShell()
 
 CALL CompHSO(hamil_SO,HD,hamil,SD,overlap_SO,NAOrbs,nshell)
 
!PRINT *, "hamil_SO matrix real part:"
!do i=1, totdim*2
!    PRINT '(1000(ES14.4))',  ( REAL(hamil_SO( i, j )), j=1,totdim*2 ) 
!end do 
!
!PRINT *, "hamil_SO matrix imaginary part:"
!do i=1, totdim*2
!    PRINT '(1000(ES14.4))',  ( AIMAG(hamil_SO( i, j )), j=1,totdim*2 ) 
!end do 

!PRINT *, "HD matrix for atom ", PrtHatom
!PRINT *, "Up-Up" 
!do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)
!    PRINT '(1000(F11.5))',  ( (HD(1, i, j )), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom) ) 
!end do  
!PRINT *, "Down-Down"
!do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)                                                     
!    PRINT '(1000(F11.5))',  ( (HD(2, i, j )), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)) 
!end do                                                            
!
!PRINT *, "Real part of hamil_so matrix for atom ", PrtHatom
!PRINT *, "Up-Up" 
!do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)
!    PRINT '(1000(F11.5))',  ( REAL(hamil_so( i, j )), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom) ) 
!end do  
!PRINT *, "Up-Down" 
!do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)
!    PRINT '(1000(F11.5))',  ( REAL(hamil_so( i, j )), j=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom) ) 
!end do  
!PRINT *, "Down-Up"
!do i=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom)                                                     
!    PRINT '(1000(F11.5))',  ( REAL(hamil_so( i, j )), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom) ) 
!end do                                   
!PRINT *, "Down-Down"
!do i=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom)                                                     
!    PRINT '(1000(F11.5))',  ( REAL(hamil_so( i, j )), j=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom) ) 
!end do                                   
!
!PRINT *, "Imaginary part of hamil_so matrix for atom ", PrtHatom
!PRINT *, "Up-Up" 
!do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)
!    PRINT '(1000(F11.5))',  ( AIMAG(hamil_so( i, j )), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom) ) 
!end do  
!PRINT *, "Up-Down" 
!do i=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom)
!    PRINT '(1000(F11.5))',  ( AIMAG(hamil_so( i, j )), j=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom) ) 
!end do  
!PRINT *, "Down-Up"
!do i=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom)                                                     
!    PRINT '(1000(F11.5))',  ( AIMAG(hamil_so( i, j )), j=LoAOrbNo(PrtHatom),HiAOrbNo(PrtHatom) ) 
!end do                                   
!PRINT *, "Down-Down"
!do i=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom)                                                     
!    PRINT '(1000(F11.5))',  ( AIMAG(hamil_so( i, j )), j=totdim+LoAOrbNo(PrtHatom),totdim+HiAOrbNo(PrtHatom) ) 
!end do  


PRINT *, "Hamil matrix for atom ",Atom," : "
PRINT *, "Up-Up" 
do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
end do  
PRINT *, "Up-Down" 
do i=LoAOrbNo(Atom),HiAOrbNo(Atom)
    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
end do  
PRINT *, "Down-Up"
do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                     
    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=LoAOrbNo(Atom),HiAOrbNo(Atom) ) 
end do                                   
PRINT *, "Down-Down"
do i=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom)                                                     
    PRINT '(1000(F11.5))',  ( REAL(hamil( i, j )), j=totdim+LoAOrbNo(Atom),totdim+HiAOrbNo(Atom) ) 
end do                                                                                             

 do i=1, totdim*2
    do j=1, totdim*2
       H_SOC(i,j)=hamil(i,j)+hamil_SO(i,j)
       S_SOC(i,j)=REAL(overlap_SO(i,j))
    end do
 end do

 return
 end subroutine spin_orbit

END MODULE device
