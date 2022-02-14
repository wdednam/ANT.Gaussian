!**********************************************************
!*********************  ANT.G-2.5.2  **********************
!**********************************************************
!*                                                        *
!*  Copyright (c) by                                      *
!*                                                        *
!*  Juan Jose Palacios (1)                                *
!*  David Jacob (2)                                       *
!*  Maria Soriano (1)                                     *
!*  Angel J. Perez-Jimenez (3)                            *
!*  Emilio SanFabian (3)                                  *
!*  Jose Antonio Antonio Verges (4)                       *
!*  Enrique Louis (5)                                     *
!*                                                        *
!* (1) Departamento de Fisica de la Materia Condensada    *
!*     Universidad Autonoma de Madrid                     *      
!*     28049 Madrid (SPAIN)                               *
!* (2) Theory Department                                  *
!*     Max-Planck-Institute for Microstructure Physics    *
!*     Halle, 06120 (GERMANY)                             *
!* (3) Departamento de Quimica Fisica                     *
!*     Universidad de Alicante                            *
!*     03690 Alicante (SPAIN)                             *
!* (4) Insto. de Ciencias de Materiales de Madrid (ICMM)  *
!*     Consejo Superior de Investigacion y Ciencia (CSIC) *
!*     28049 Madrid (SPAIN)                               *
!* (5) Departamento de Fisica Aplicada                    *
!*     Universidad de Alicante                            *      
!*     03690 Alicante (SPAIN)                             *
!*                                                        *
!**********************************************************
  SUBROUTINE ANT (UHF,JCycle,IRwGen,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig,denerrj,Crit,ANTOn,NBasis)
!**********************************************************************************************************************
!* Interface subroutine with Gaussian                                                                                 *
!**********************************************************************************************************************
  USE Parameters, ONLY: SL, SwOffSPL, alpha, Read_Parameters, Write_Parameters, NSpinLock, npulay
  USE Parameters, ONLY: ChargeAcc,ChargeA,FermiAcc,FermiA,PAcc,PA,FullAcc,RedTransmB,RedTransmE,ElType,LDOS_Beg,LDOS_End
  USE Parameters, ONLY: Mulliken, Hamilton, PFix, DFTU, FMixing, SOC, ROT, HFEnergy
  USE constants, ONLY: Hart, d_zero
  USE Numeric, ONLY: RTrace  
  USE preproc
  USE OneDLead, only: CleanUp1DLead 
  USE device, ONLY: InitDevice, DevFockMat, DevDensMat, ReadDensMat, LeadsOn, DevShift, SwitchOnLeads, &
       SwitchOnEvaluation, SwitchOnSecant, SwitchOffSecant, SwitchOnSpinLock, SwitchOffSpinLock, SwitchOn1DElectrodes, &
       SwitchOff1DElectrodes, SwitchOnChargeCntr, SwitchOffChargeCntr, transport, CleanUpDevice, SetDevDensMat, ReadFockMat
#ifdef G03ROOT
  USE g03Common, ONLY: GetNShell, GetAtm4Sh, Get1stAO4Sh, GetNBasis, GetAN, GetAtmChg, GetAtmCo, GetNAtoms
#endif
#ifdef G09ROOT
  USE g09Common, ONLY: GetNShell, GetAtm4Sh, Get1stAO4Sh, GetNBasis, GetAN, GetAtmChg, GetAtmCo, GetNAtoms, GetAtTyp, GetNE, GetIChg, GetMultip  
#endif
  use ANTCommon
  use util
  IMPLICIT NONE

  ! dummy arguments
  LOGICAL, INTENT(in)      :: UHF          
  REAL*8,INTENT(in)        :: denerrj, Crit
  INTEGER, INTENT(in)    :: NBasis,IRwGen,IRwH,IRwPA,IRwPB,IRwFA,IRwFB,IRwS1,IRwEig
  INTEGER, INTENT(inout) :: JCycle
  LOGICAL,INTENT(inout)    :: ANTOn
  LOGICAL :: ADDP
  LOGICAL, SAVE :: FInit

  ! local variables
  REAL*8    :: val, density, fock, exspin, norma, ENR, E1A, E1B, Ex, Ec, EJ, ETot
  INTEGER :: ic, ndim, i, j, info, ii, jj, ispin, acount, AllocErr, ios ,k, iatom, jatom, ntot, NAtoms, NBas6D, mdv, isym, jsym, ISCF, IROHF, &
             ipa, ipb, icja, icjb, iv, mdv1, IPrint, IGrid, IHMeth, lseall, FMFlag, FMFlg1, IPFlag, NFxFlg, IOpCl, ICntrl, IExCor, IExchn, ICorr, &
             AccDes, IPureD,IPureF, NOpUse , NBTI, ISym2E, IRadAn, IRanWt, IRanGd, NMtPBC, IJunk
  INTEGER, DIMENSION(MaxAtm) :: NAO                                                                                    

  CHARACTER(len=50), SAVE :: densitymatrix, fockmatrix

  INTEGER, SAVE :: isw, NCycLeadsOn, NSpin
  INTEGER, SAVE :: ntt, ntt6d, len
  INTEGER, PARAMETER :: zero=0, two=2, one=1

  REAL*8, DIMENSION(:),ALLOCATABLE   :: TS, VV, PDA, PDB, CJA, CJB, IIAn, IIAtTyp, AAtmChg
  LOGICAL, DIMENSION(50) ::  AllowP
  LOGICAL ::  FMM, HaveDB, OKSym, AtBlock
  REAL*8, DIMENSION(1), SAVE :: XX = 0.0d0
  INTEGER, DIMENSION(1), SAVE :: JJJ = 0
  LOGICAL, DIMENSION(1), SAVE :: LLJ = .False. 
  REAL*8, DIMENSION(4) :: ScaDFX
  REAL*8 :: ScaHFX, RJunk
  REAL*8, DIMENSION(5,6) :: Omega 
  INTEGER, DIMENSION(200) :: IOpArr
  REAL*8, DIMENSION(:,:),ALLOCATABLE :: S, CC

  !
  ! Matrices in lower trinagular form (alpha and beta) 
  ! for communication with gaussian
  !
  REAL*8, DIMENSION(:,:),ALLOCATABLE,SAVE ::  MT, CT
  !
  ! Fock and density matrices in regular form
  !
  REAL*8, DIMENSION(:,:,:),ALLOCATABLE,SAVE :: F, D, HC
  REAL*8, DIMENSION(:,:,:,:),ALLOCATABLE,SAVE ::  DD
  REAL*8, DIMENSION(:,:,:),ALLOCATABLE,SAVE ::  DDD
  REAL*8, DIMENSION(npulay+1,npulay+1) ::  b
  REAL*8, DIMENSION(npulay+1) ::  c
  INTEGER, DIMENSION(npulay+1) ::  ipiv

  real*8 :: fmix
  
  !*******************************************************************
  ! Before first cycle: Initialize module device, allocate memory etc 
  !*******************************************************************
  IF(JCycle.EQ.0) THEN

     NSpin = 1
     IF(UHF) NSpin = 2
     ntt = (NBasis*(NBasis+1))/2
     Call GetNB6(NBas6D)
     ntt6d = (NBas6D*(NBas6D+1))/2
     
     ! Read IOp common block for later use
     CALL FileIO(two,-996,200,IOpArr,zero)    
     
     !WRITE(ifu_log,'(200(I17))') IOpArr     
     
     ANTOn = .FALSE.
          
     PRINT *
     PRINT *, " ****************************************************************** "
     PRINT *, " ***                                                            *** "
     PRINT *, " ***                      A l i c a n t e                       *** "
     PRINT *, " ***                      N a n o                               *** "
     PRINT *, " ***                      T r a n s p o r t                     *** "
#ifdef G03ROOT
     PRINT *, " ***                      G 0 3                                 *** "
#endif
#ifdef G09ROOT
     PRINT *, " ***                      G 0 9                                 *** "
#endif
     PRINT *, " ***                                                            *** "
     PRINT *, " ****************************************************************** "
     PRINT *, " ***                     Version: 2.5.2                         *** "
     PRINT *, " ****************************************************************** "
     PRINT *, " *  Copyright (c) by                                              * "
     PRINT *, " *                                                                * "
     PRINT *, " *  Juan Jose Palacios (1)                                        * "
     PRINT *, " *  David Jacob (2)                                               * "
     PRINT *, " *  Maria Soriano (1,5)                                           * "
     PRINT *, " *  Angel J. Perez-Jimenez (3)                                    * "
     PRINT *, " *  Emilio SanFabian (3)                                          * "
     PRINT *, " *  Jose Antonio Antonio Verges (4)                               * "
     PRINT *, " *  Enrique Louis (5)                                             * "
     PRINT *, " *  Wynand Dednam (5)                                             * "
     PRINT *, " *                                                                * "
     PRINT *, " * (1) Departamento de Fisica de la Materia Condensada            * "
     PRINT *, " *     Universidad Autonoma de Madrid                             * "
     PRINT *, " *     28049 Madrid (SPAIN)                                       * "
     PRINT *, " * (2) Theory Department                                          * "
     PRINT *, " *     Max-Planck-Institute for Microstructure Physics            * "
     PRINT *, " *     Halle, 06120 (GERMANY)                                     * "
     PRINT *, " * (3) Departamento de Quimica Fisica                             * " 
     PRINT *, " *     Universidad de Alicante                                    * " 
     PRINT *, " *     03690 Alicante (SPAIN)                                     * "
     PRINT *, " * (4) Insto. de Ciencias de Materiales de Madrid (ICMM)          * "
     PRINT *, " *     Consejo Superior de Investigacion y Ciencia (CSIC)         * "
     PRINT *, " *     28049 Madrid (SPAIN)                                       * "
     PRINT *, " * (5) Departamento de Fisica Aplicada                            * "
     PRINT *, " *     Universidad de Alicante                                    * "
     PRINT *, " *     03690 Alicante (SPAIN)                                     * "
     PRINT *, " *                                                                * "
     PRINT *, " ****************************************************************** "
     PRINT *
     PRINT *, " ****************************************************************** "
     PRINT *, " *                                                                * "
#ifdef G03ROOT
     PRINT *, " *              Initializing ANT.G03 modules                      * "
#endif
#ifdef G09ROOT
     PRINT *, " *              Initializing ANT.G09 modules                      * "
#endif
     PRINT *, " *                                                                * "
     PRINT *, " ****************************************************************** "
     PRINT *

     ! Reading parameter file
     OPEN(UNIT=ifu_nam,FILE="name",IOSTAT=ios,STATUS='OLD')
     IF( ios == 0 ) THEN
        READ(ifu_nam,*) jobname
        OPEN(UNIT=ifu_ini,FILE=trim(jobname)//'.ant',IOSTAT=ios,STATUS='OLD')
        IF( ios == 0 ) THEN
          CALL read_parameters( ifu_log, ifu_ini, ios, jobname )
        ELSE
          WRITE(ifu_log,*) " "
          WRITE(ifu_log,*) "*** No parameter file found. Using default values."
          WRITE(ifu_log,*) " "
        END IF
        CLOSE(ifu_ini)
     ELSE
        WRITE(ifu_log,*) " "
        WRITE(ifu_log,*) "*** No name file found. Using default values."
        WRITE(ifu_log,*) " "
     END IF
     CLOSE(UNIT=ifu_nam)
     ant1dname = jobname
     call locase(ant1dname)

     CALL write_parameters( ifu_log )

     !Opening writting files
     OPEN(ifu_xyz,file=trim(jobname)//'.xyz',status='unknown')
     IF (ElType(1) /= 'GHOST' .and. ElType(2) /= 'GHOST' .and. (SOC .or. ROT)) OPEN(ifu_tra,file='T.'//trim(jobname)//'.SOC.dat',status='unknown')
     IF (ElType(1) /= 'GHOST' .and. ElType(2) /= 'GHOST' .and. (.not. SOC) .and. (.not. ROT)) OPEN(ifu_tra,file='T.'//trim(jobname)//'.dat',status='unknown')
     !IF (ElType(1) /= 'GHOST' .and. ElType(2) /= 'GHOST') OPEN(ifu_tra,file='T.'//trim(jobname)//'.dat',status='unknown')
     IF (RedTransmB < RedTransmE) OPEN(ifu_red,file='t.'//trim(jobname)//'.dat',status='unknown')
     IF (Hamilton) OPEN(ifu_ham,file='V.'//trim(jobname)//'.dat',status='unknown')
     IF (Mulliken) OPEN(ifu_mul,file='Q.'//trim(jobname)//'.dat',status='unknown')
     IF (LDOS_Beg <= LDOS_End) OPEN(ifu_dos,file='DOS.'//trim(jobname)//'.dat',status='unknown')

     ! Creating density matrix file name
     densitymatrix='P.' // trim(jobname) // '.dat'
     ! Creating Fock matrix file name
     fockmatrix='F.' // trim(jobname) // '.dat'
     !
     ! Allocate memory for dynamic arrays
     !
     ALLOCATE( MT( NSpin, ntt ), STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not allocate memory for MT(:)"
        STOP
     END IF
     ALLOCATE( F(NSpin,NBasis,NBasis), STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not allocate memory for F(:,:,:)"
        STOP
     END IF
     if(.not. FMixing )then
        ALLOCATE( D(NSpin,NBasis,NBasis), DD(npulay+1,NSpin,NBasis,NBasis),DDD(npulay,NSpin,NBasis*NBasis), STAT=AllocErr )
        IF( AllocErr /= 0 ) THEN
           PRINT *, "ESF/Error: could not allocate memory for D(:,:,:), DD(:,:,:,:),DDD(:,:,:)"
           STOP
        END IF
        DD=0.0
        DDD=0.0
     end if
     if( HFEnergy )then  
        ALLOCATE( CT( NSpin, ntt ), STAT=AllocErr )
        IF( AllocErr /= 0 ) THEN
           PRINT *, "ESF/Error: could not allocate memory for CT(:,:)"
           STOP
        END IF
        ALLOCATE( HC(NSpin,NBasis,NBasis), STAT=AllocErr )
        IF( AllocErr /= 0 ) THEN
           PRINT *, "ESF/Error: could not allocate memory for HC(:,:,:)"
           STOP
        END IF
        HC=0.0        
     end if     

     !
     ! Obtain overlap matrix in lower triangular form
     ! and transform to regular form matrix
     !
     ALLOCATE( TS(ntt), S(NBasis,NBasis), STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not allocate memory for S(:,:), TS(:)"
        STOP
     END IF

     CALL FileIO(two,-IRwS1,ntt,TS,zero)
     acount = 1
     DO i=1,NBasis
        DO j=1,i
           S(i,j) = TS(acount)
           S(j,i) = TS(acount)
           acount = acount+1
        END DO
     END DO  

     ! 
     ! Initialize module device 
     !
     CALL InitDevice( NBasis, UHF, S )           

     CLOSE(ifu_xyz)

     if(FMixing)OPEN(ifu_fm,file=fockmatrix,iostat=ios,status='old')
     if(.not.FMixing)OPEN(ifu_dm,file=densitymatrix,iostat=ios,status='old')
     IF (ios /= 0) THEN
        FInit=.false.
     ELSE
        FInit=.true.
        if(FMixing)CLOSE(ifu_fm)
        if(.not.FMixing)CLOSE(ifu_dm)
     END IF

     IF( FInit )THEN
        !
        ! Initialization from file F.dat or P.dat
        !
        if(FMixing)then
           CALL ReadFockMat(fockmatrix)
           DO i=1,NBasis
              DO j=1,NBasis
                 F(1,i,j) = DevFockMat(1,i,j)
                 if(UHF) F(2,i,j) = DevFockMat(2,i,j)
              end DO
           end DO
        else
           CALL ReadDensMat(densitymatrix)        
           !
           ! Transform to triangular form
           !
           acount = 1
           DO i=1,NBasis
              DO j=1,i
                 D(1,i,j)=DevDensMat(1,i,j)
                 MT(1,acount)=D(1,i,j)
                 IF(UHF) D(2,i,j) = DevDensMat(2,i,j)
                 IF(UHF) MT(2,acount) = D(2,i,j)
                 acount = acount +1
              ENDDO
           ENDDO
           !
           ! PASS DENSITY MATRIX TO GAUSSIAN VIA RWF
           !
           CALL FileIO(one,-IRwPA,ntt,MT(1,:),zero)
           IF (UHF) CALL FileIO(one,-IRwPB,ntt,MT(2,:),zero)
        end if
     ENDIF

     DEALLOCATE( S, TS, STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Error: could not DEAllocate memory for S(:,:), TS(:)"
        STOP
     END IF
     
     RETURN
  END IF

  IF( .NOT. LeadsOn() )THEN
     !
     ! Connect Leads in first cycle
     ! when reinitialized
     !
     IF( FInit )THEN
        ANTOn = .true.
        CALL SwitchOnLeads
        NCycLeadsOn = 0
        !
        ! Or if SL criterion is met:
        !
     ELSEIF( denerrj .LE. SL ) THEN
        ANTOn = .true.
        Call SwitchOnLeads
        NCycLeadsOn = 0
        if(.not.FMixing)then
           !
           ! Obtain density matrix from gaussian only in first step 
           ! with leads switched on and if not initialized from file
           !
           CALL FileIO(two,-IRwPA,ntt,MT(1,:),zero)
           IF(UHF) CALL FileIO(two,-IRwPB,ntt,MT(2,:),zero)
           acount = 1
           DO i=1,NBasis
              DO j=1,i
                 D(1,i,j)=MT(1,acount)
                 D(1,j,i)=D(1,i,j)
                 call SetDevDensMat( 1, i, j, D(1,i,j) )
                 call SetDevDensMat( 1, j, i, D(1,j,i) )
                 IF (UHF) THEN
                    D(2,i,j)=MT(2,acount)
                    D(2,j,i)=D(2,i,j)
                    call SetDevDensMat( 2, i, j, D(2,i,j) )
                    call SetDevDensMat( 2, j, i, D(2,j,i) )
                 ENDIF
                 acount = acount+1
              ENDDO
           ENDDO
        end if
     END IF
  END IF

  IF( jcycle.EQ.1000 )THEN
     NCycLeadsOn = 999
     CALL SwitchOnEvaluation()
     CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. .not. FullAcc .and. jcycle > 2) THEN
     IF( denerrj .LT. Crit*10) THEN
        ChargeAcc = ChargeA
        FermiAcc = FermiA
        PAcc = PA
     ELSE IF( denerrj .LT. Crit*100) THEN
        ChargeAcc = ChargeA*10.0
        FermiAcc = FermiA*10.0
        PAcc = PA*10.0
     ELSE IF( denerrj .LT. Crit*1000) THEN
        ChargeAcc = ChargeA*100.0
        FermiAcc = FermiA*100.0
        PAcc = PA*100.0
     ELSE
        ChargeAcc = ChargeA*1000.0
        FermiAcc = FermiA*1000.0
        PAcc = PA*1000.0
     END IF
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. FullAcc) THEN
     ChargeAcc = ChargeA
     FermiAcc = FermiA
     PAcc = PA
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE IF( LeadsOn() .and. .not. FullAcc .and. jcycle <= 2) THEN
     ChargeAcc = ChargeA*10.0
     FermiAcc = FermiA*10.0
     PAcc = PA*10.0
     IF( NSpinLock .LT. 0 .AND. denerrj .LT. SwOffSPL  ) CALL SwitchOffSpinLock()
     IF( NSpinLock .GE. 0 .AND. NCycLeadsOn .GE. NSpinLock ) CALL SwitchOffSpinLock()
  ELSE
     !
     ! If leads are not switched on or if not in cycle 1000 return to Gaussian
     ! without calculation of density matrix with leads connected
     !
     RETURN
  END IF

  
  !*********************************************
  ! COMPUTE DENSITY MATRIX WITH LEADS CONNECTED
  !*********************************************

  WRITE(ifu_log,*) '-------------------------------------------------------------------------'
  WRITE(ifu_log,*) 'Computing the density matrix at cycle:', JCycle
  WRITE(ifu_log,*) '-------------------------------------------------------------------------'
  
  ! Calculate Total energy in open boundary Hartree-Fock calculations. See Example link program in Gaussian Programming reference
  IF(HFEnergy .and. jcycle.EQ.1000 ) THEN 
    
    CALL execute_command_line('grep " MDV=  " *.log | cut -b61-79 > mdv.txt' ) 
    OPEN(unit=10,file='mdv.txt') 
    READ(10,*) mdv 
    CLOSE(10)
    CALL execute_command_line('rm -r mdv.txt')
    
    !print *, "MDV = ",mdv
    
    ALLOCATE( VV( mdv ), PDA (ntt6d), PDB(ntt6d), CJA (ntt6d), CJB(ntt6d), STAT=AllocErr )
    IF( AllocErr /= 0 ) THEN
       PRINT *, "ESF/Error: could not allocate memory for VV(:), PDA(:), PDB(:), CJA(:), CJB(:)"
       STOP
    END IF  
    
    NAtoms = GetNAtoms()    
    
    ALLOCATE( IIAn( NAtoms ), IIAtTyp(NAtoms), AAtmChg(NAtoms), CC(3,NAtoms), STAT=AllocErr )
    IF( AllocErr /= 0 ) THEN
       PRINT *, "ESF/Error: could not allocate memory for IIAn(:),IIAtTyp(:),AAtmChg(:),CC(:,:)"
       STOP
    END IF
    
    DO i=1,NAtoms
         IIAn(i) = GetAN(i)
         !Print *, "Atomic number of Atom ",i," is ", IIAn(i)
         IIAtTyp(i) = GetAtTyp(i)
         !Print *, "Atom type of Atom ",i," is ", IIAtTyp(i)
         AAtmChg(i) = GetAtmChg(i)
         !Print *, "Atomic charge of Atom ",i," is ", AAtmChg(i)
         CC(1,i) = GetAtmCo(1,i)
         CC(2,i) = GetAtmCo(2,i)
         CC(3,i) = GetAtmCo(3,i)      
         !Print *, "Coordinates of Atom ",i," are ", CC(1,i), CC(2,i), CC(3,i)
    ENDDO    
             
    CALL GetNOp(isym,jsym)  
    
    CALL ILSW(two,one,IOpCl)
    
    ipa = 1
    ipb = ipa + ntt6d
    
    !icja = ipb + ntt6d
    !icjb = icja + ntt6d    
    
    !iv = icjb + 2*ntt6d
    
    iv = ipb + 2*ntt6d
    
    CALL TstCor(iv,mdv, 'ANT')
    
    mdv1 = mdv - iv + 1    

    CALL ILSW(2,1,ISCF) 
    CALL ILSW(2,22,IROHF)
    
    CALL FileIO(two,-IRwFA,ntt,VV(ipa),zero)
    
    IF (UHF) Call FileIO(two,-IRwFB,ntt,VV(ipb),zero) 
       
    
    PDA = VV(ipa)
    IF(UHF) PDB = VV(ipb)   
    
    !CALL ILSW(2,1,ISCF) 
    !CALL ILSW(2,22,IROHF)
    !
    !CALL FileIO(two,-528,ntt,VV(ipa),zero)
    !
    !IF(ISCF.eq.0 .and. IROHF.eq.0) THEN
    !   Call AScale(ntt,0.5d0,VV(ipa),VV(ipa)) 
    !   Call AMove(ntt,VV(ipa),VV(ipa))
    !ELSE
    !   Call FileIO(two,-530,ntt,VV(ipb),zero) 
    !ENDIF    
    !
    !PDA = VV(ipa)
    !IF(UHF) PDB = VV(ipb)   
    !
    CALL FileIO(two,-501,one,ENR,40)   
    !    
    !CALL FileIO(two,-515,ntt,CJA,zero)
    !
    !IF(UHF) THEN
    !  CALL FileIO(two,515,ntt,CJB,zero)
    !END IF
    !
    !E1A = 0.0d0
    !acount = 1
    !DO i=1,NBasis
    !   DO j=1,i
    !      E1A=E1A+PDA(acount)*CJA(acount)
    !      IF (UHF) THEN
    !         E1B=E1B+PDB(acount)*CJB(acount)
    !      ENDIF
    !      acount = acount +1
    !   ENDDO
    !ENDDO    
    !
    !print *,"E1A = ", E1A
    !IF (UHF) THEN 
    !   print *,"E1B = ", E1B
    !ENDIF   
        
    IPrint = IOpArr(33)
    IGrid = IOpArr(90)
    WRITE(ifu_log,'(A,I10,A)') ' Evaluate the HF functional using Xc grid ',IGrid,'.'    
        
    !CALL SetPFL(6,Iprint,zero,IPFlag,AllowP,FMM,FMFlag,FMFlg1,NFxFlg,IHMeth,Omega,one,lseall,iv,VV,mdv) 
    
    !WRITE(ifu_log,'(A,I10)') ' IPFLag is ', IPFLag
    !WRITE(ifu_log,'(A,I10)') ' FMM is ', FMM
    !WRITE(ifu_log,'(A,I10)') ' FMFlag is ', FMFlag
    !WRITE(ifu_log,'(A,I10)') ' FMFlg1 is ', FMFlg1
    !WRITE(ifu_log,'(A,I10)') ' NFxFlg is ', NFxFlg
    !WRITE(ifu_log,'(A,I10)') ' Hamiltonian type is ', IHMeth
    !WRITE(ifu_log,'(A,I10)') ' LSEALL is ', lseall
    !WRITE(ifu_log,*) ' AllowP is ', AllowP        
    !WRITE(ifu_log,*) ' Omega is ', Omega
    
    !mdv1 = mdv - iv + 1
    
    !HaveDB=IOpCl.eq.1
    !CALL PurCar(6,IPrint, .False., .True., .False., .False., .False.,.False., .False. ,HaveDB, .True., &
    !     NBasis,NBas6D,one,one,one,zero,zero,IPureD,IPureF,PDA,PDB,XX,XX,XX,XX,XX,XX,XX,VV(iv),mdv1)  
         
    
    CALL execute_command_line('grep " ISym2E= " *.log | cut -b68-69 > isym2e.txt' ) 
    OPEN(unit=10,file='isym2e.txt') 
    READ(10,*) ISym2E
    CLOSE(10)
    CALL execute_command_line('rm -r isym2e.txt')  
    
    !print *, "ISym2E = ",ISym2E  
    
    OKSym = ISym2E.ne.0
    AtBlock = .not. OKSym
	
    CALL XLSW(two, 11, 1, IExCor,RJunk)
    CALL XLSW(two, 13, 1, IRadAn,RJunk)
    CALL XLSW(two, 35, 1, IRanWt,RJunk)
    CALL XLSW(two, 36, 1, IRanGd,RJunk)
    CALL XLSW(two, 46, 1, IHMeth,RJunk)
    CALL XLSW(two, 51, 1, NMtPBC,RJunk)  
    CALL XLSW(two,1,2,IJunk,ScaDFX)
    CALL XLSW(two,2,2,IJunk,ScaHFX)
    CALL XLSW(two,3,2,IJunk,Omega)

    WRITE(ifu_log,'(A,I10)') ' IExCor is ', IExCor
    WRITE(ifu_log,'(A,I10)') ' IRadAn is ', IRadAn
    WRITE(ifu_log,'(A,I10)') ' IRanWt is ', IRanWt
    WRITE(ifu_log,'(A,I10)') ' IRanGd is ', IRanGd
    WRITE(ifu_log,'(A,I10)') ' NMtPBC is ', NMtPBC
    WRITE(ifu_log,'(A,I10)') ' Hamiltonian type is ', IHMeth
    WRITE(ifu_log,*) ' ScaDFX is ', ScaDFX
    WRITE(ifu_log,*) ' ScaHFX is ', ScaHFX
    WRITE(ifu_log,*) ' Omega is ', Omega
        
    AccDes = 0.0d0
    
    CALL HarFok(6,IPrint,one,.true.,.false.,OKSym,AtBlock,&
       .false.,IExCor,IRadAn,IRanWt,IRanGd,AccDes,IOpCl,.true.,NAtoms,&
       NBasis,NBas6D,NBTI,GetIChg(),GetMultip(),GetNE(),IIAn,IIAtTyp,AAtmChg,CC,ScaDFX,&
       NMtPBC,ETot,PDA,PDB,JJJ,XX,zero,JJJ,JJJ,XX,&
       Omega,VV(iv),mdv1)    
           
    
!    CALL ASet(4,1.0d0,ScaDFX)              
!    
!    IExChn = 2
!    ICorr = 0
!    ICntrl = 1
!    AccDes = 0.0d0
!    CALL CalDSu(6,IPrint,IPFlag,one,one,one,NAtoms, .False. ,Ex,Ec,JJ,XX,&
!                XX,XX,XX,XX,XX,XX,XX,ICntrl,IExchn,ICorr,zero,IGrid,zero,zero,zero,ScaDFX,&
!                CC,IIAn,IIAtTyp,AAtmChg,NAtoms,NAtoms,JJ,PDA,XX,XX,PDB,XX,XX,zero,XX,zero,&
!                XX,NBas6D,IOpCl,AccDes, .False., .False. ,one,zero,one,one,one,JJJ,JJJ,JJJ,XX,JJJ, &
!                 XX,GetNE(),zero,JJJ,JJJ,XX,XX,Omega,VV(iv),mdv1)                                  
!                 
!    ICntrl = 500
!    NOpUse = one
!    ScaHFX = 0.0d0
!    NBTI = zero
!    If(IOpCl.ne.zero) Call ASASB(ntt6d,0.5d0,PDA,0.5d0,PDB) 
!    Call FoFJK(6,IPrint,IHMeth,zero,zero,ICntrl,Omega,zero,FMM,FMFlag, &
!               FMFlg1,NFxFlg,IPFlag,AllowP,lseall,VV(lseall),one,NAtoms,.True., &
!               .True., .False., .False., .False. ,zero,AccDes,ScaHFX,one,one,one,zero,one,zero,zero, &
!                NBasis,zero,NOpUse,one,JJJ,JJJ,JJJ,XX,JJJ,NBTI,JJJ,PDA,PDB,XX,XX,XX,XX,XX, &
!                XX,XX,JJJ,CJA,CJB,JJJ,XX,XX,NAtoms,IIAn,AAtmChg,CC,IIAtTyp,one,zero,XX,XX, &
!                XX,JJJ,XX,RJunk,1.0d0,VV(iv),mdv1)
!    
!    EJ = 0.0d0
!    ! Transform from lower trinagular form to regular form 
!    acount = 1
!    DO i=1,NBas6D
!       DO j=1,i
!          EJ=EJ+PDA(acount)*CJA(acount)
!          acount = acount +1
!       ENDDO
!    ENDDO    
        
    !write(ifu_log,'(A)') ' Hartree-Fock energy breakdown: '
    !write(ifu_log,'(A,F24.12)') ' ENR =     ', ENR
    !write(ifu_log,'(A,F24.12)') ' E1A =     ', E1A
    !write(ifu_log,'(A,F24.12)') ' E1B =     ', E1B
    !write(ifu_log,'(A,F24.12)') ' Ex =      ', Ex
    !write(ifu_log,'(A,F24.12)') ' Ec =      ', Ec    
    !write(ifu_log,'(A,F24.12)') ' EJ =      ', EJ
    !write(ifu_log,'(A,F24.12)') ' ETotal =  ', ENR+E1A+E1B+EJ
    write(ifu_log,'(A,F24.12)') ' ETotal =  ', Etot      
    
    DEALLOCATE(VV, PDA, PDB, CJA, CJB, STAT=AllocErr )
    IF( AllocErr /= 0 ) THEN
       PRINT *, "ESF/Deallocation error for VV(:),PDA(:),PDB(:),CJA(:),CJB(:)"
       STOP
    END IF        
    DEALLOCATE(IIAn, IIAtTyp, AAtmChg, CC, STAT=AllocErr )
    IF( AllocErr /= 0 ) THEN
       PRINT *, "ESF/Deallocation error for IIAn(:),IIAtTyp(:),AAtmChg(:),CC(:,:)"
       STOP
    END IF     
    
    CALL FileIO(two,-IRwH,ntt,CT(1,:),zero)
    IF(UHF) CALL FileIO(two,IRwH,ntt,CT(2,:),zero)   
         
  END IF    

  ! Obtain Fock matrix from Gaussian RWF
  CALL FileIO(two,-IRwFA,ntt,MT(1,:),zero)
  IF(UHF) CALL FileIO(two,-IRwFB,ntt,MT(2,:),zero)
  fmix=1.0d0
  ! For Fock matrix mixing 
  if(FMixing.and.NCycLeadsOn>0) fmix=alpha
  !if(NCycLeadsOn>0) fmix=alpha
  if(FMixing.and.FINIT.and.NCycLeadsOn.eq.0) fmix=0.0d0
  ! Transform from lower trinagular form to regular form
  !print *, "NCycLeadsOn=", NCycLeadsOn
  !print *, "fmix =", fmix
  acount = 1
  DO i=1,NBasis
     DO j=1,i
        ! Mixing with old Fock matrix if fmix<1
        F(1,i,j)=(1.0d0-fmix)*F(1,i,j)+fmix*Hart*MT(1,acount)
        F(1,j,i)=F(1,i,j)
        IF(HFEnergy  .and. jcycle.EQ.1000 ) THEN 
          HC(1,i,j)=CT(1,acount)
          HC(1,j,i)=HC(1,i,j)                       
        ENDIF  
        IF (UHF) THEN
           F(2,i,j)=(1.0d0-fmix)*F(2,i,j)+fmix*Hart*MT(2,acount)
           F(2,j,i)=F(2,i,j)
           IF(HFEnergy  .and. jcycle.EQ.1000 ) THEN 
             HC(2,i,j)=CT(2,acount)
             HC(2,j,i)=HC(2,i,j)                   
           ENDIF             
        ENDIF
        acount = acount +1
     ENDDO
  ENDDO
  
  if(FMixing)then
     ntot=0
     DO i=1,GetNShell()-1
        IF (GetAtm4Sh(i).NE.GetAtm4Sh(i+1)) THEN
           NAO(GetAtm4Sh(i))=Get1stAO4Sh(i+1)-(ntot+1)
           ntot=ntot+NAO(GetAtm4Sh(i))
        ENDIF
     ENDDO
     NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot
     ! Write density matrix to file F.' // trim(jobname) // '.dat
     if (NCycLeadsOn  < 1000) then
        OPEN(ifu_fm,file=fockmatrix,status='unknown')
        WRITE(ifu_fm,*)DevShift()
        DO ispin=1,NSpin
           i=0
           do iAtom=1,GetNAtoms()
              do ii=1,NAO(iAtom)
                 i=i+1
                 j=0
                 do jAtom=1,GetNAtoms() 
                    do jj=1,NAO(jAtom)
                       j=j+1
                       if (i >= j) WRITE(ifu_fm,'(i2,2i6,e18.10,i5,i6,i5,i6)')ispin,i,j,F(ispin,i,j),iAtom,ii,jAtom,jj
                    end do
                 end do
              end do
           end do
        end do
        close(ifu_fm)
     end if
  end if
  
  ! Count cycles with leads connected
  NCycLeadsOn = NCycLeadsOn + 1
  !print*,NCycLeadsOn
  
  ! Turn on charge control every 5 steps in the first cycles
  IF(MOD(NCycLeadsOn-1,10) == 0 .and. NCycLeadsOn <= 21) CALL SwitchOnChargeCntr
  IF(MOD(NCycLeadsOn-1,20) == 0 .and. NCycLeadsOn > 21) CALL SwitchOnChargeCntr
    
  ! Call subroutine that solves transport problem
  CALL Transport(F,HC,ENR,JCYCLE,ADDP)
  
   if (ElType(1) == "1DLEAD" .and. ElType(2) == "1DLEAD") THEN
       IF(MOD(NCycLeadsOn-1,10) == 0 .and. NCycLeadsOn >= 6) THEN
          CALL CleanUp1DLead(1)
          CALL CleanUp1DLead(2)
          CALL SwitchOff1DElectrodes
       END IF 
   end if   

  CALL SwitchOffChargeCntr
  
  IF( SL <= 0.0d0 ) alpha = 1.0d0
 
  !
  ! Pulay accelaration for density matrix mixing (FMixing == false)
  !
  if(.not.FMixing)then 
  
     DD(1,:,:,:)=D(:,:,:)
  
     DO k=npulay,1,-1
        DD(k+1,:,:,:)=DD(k,:,:,:)
     END DO
  
     if (ADDP) then
        DO ispin=1,NSpin
           DO i=1,NBasis
              DO j=1,NBasis
                 DD(1,ispin,i,j)=(1.0d0-alpha)*DD(2,ispin,i,j)+alpha*DevDensMat(ispin,i,j)           
              END DO
           END DO
        END DO
     else
        DD(1,:,:,:)=0.5*DD(3,:,:,:)+0.5*DD(2,:,:,:)
     end if
  
     DO k=1,npulay
        DO ispin=1,NSpin
           ic=0
           DO i=1,NBasis
              DO j=1,NBasis
                 ic=ic+1
                 DDD(k,ispin,ic)=DD(k,ispin,i,j)-DD(k+1,ispin,i,j)   
              END DO
           END DO
        END DO
     END DO

     ! Pulay acceleration procedure
     ndim=npulay+1
     if (mod(NCycLeadsOn,ndim)==npulay) then
        print *, "Pulay kick ............."
     
        D=0.0
        DO ispin=1,NSpin
           b=0.0
           c=0.0
           c(ndim)=-1.0
           !Matrix to obtain the coefficients
           DO k=1,ndim
              DO j=1,ndim
                 if (k.lt.ndim.and.j.lt.ndim) then
                    DO i=1,NBasis*NBasis
                       b(k,j)=b(k,j)+DDD(k,ispin,i)*DDD(j,ispin,i)
                    END DO
                 else if (k==ndim.and.j.lt.ndim) then
                    b(k,j)=-1.0
                 else if (j==ndim.and.k.lt.ndim) then
                    b(k,j)=-1.0
                 end if
              END DO
           END DO
        
           Call DGETRF( ndim, ndim, b, ndim, IPIV, INFO) 
           Call DGETRS( 'N', ndim, 1, b, ndim, IPIV, c, ndim, INFO)
        
           DO k=1,npulay
              DO i=1,NBasis
                 DO j=1,NBasis
                    D(ispin,i,j)=D(ispin,i,j)+c(k)*DD(k,ispin,i,j)
                 END DO
              END DO
           END DO
        END DO
  
     ELSE
        D=DD(1,:,:,:)
     END IF
     ! End Pulay accelaration

     !     
     ! For DFT+U calculations the damped density matrix has
     ! to be fed back to device module in order to calculate
     ! the correct DFT+U potential in the next step 
     !  
     if( DFTU )then
        do ispin=1,NSpin
           do i=1,NBasis
              do j=1,NBasis
                 call SetDevDensMat( ispin, i, j, D(ispin,i,j) )
              end do
           end do
        end do
     end if
  
     ntot=0
     DO i=1,GetNShell()-1
        IF (GetAtm4Sh(i).NE.GetAtm4Sh(i+1)) THEN
           NAO(GetAtm4Sh(i))=Get1stAO4Sh(i+1)-(ntot+1)
           ntot=ntot+NAO(GetAtm4Sh(i))
        ENDIF
     ENDDO
     NAO(GetAtm4Sh(GetNShell()))=GetNBasis()-ntot

     ! Write density matrix to file P.' // trim(jobname) // '.dat
     if (NCycLeadsOn  < 1000) then
        OPEN(ifu_dm,file=densitymatrix,status='unknown')
        WRITE(ifu_dm,*)DevShift()
        DO ispin=1,NSpin
           i=0
           do iAtom=1,GetNAtoms()
              do ii=1,NAO(iAtom)
                 i=i+1
                 j=0
                 do jAtom=1,GetNAtoms() 
                    do jj=1,NAO(jAtom)
                       j=j+1
                       if (i >= j) WRITE(ifu_dm,'(i2,2i6,e18.10,i5,i6,i5,i6)')ispin,i,j,D(ispin,i,j),iAtom,ii,jAtom,jj
                    end do
                 end do
              end do
           end do
        end do
        close(ifu_dm)
     end if
     
     if (NCycLeadsOn == 1 .and. PFix .and. .not. FInit) then
        CALL ReadDensMat(densitymatrix)        
        do ispin=1,NSpin
           do i=1,NBasis
              do j=1,NBasis
                 D(ispin,i,j)=DevDensMat(ispin,i,j)
              end do
           end do
        end do
     end if

  end if

  !*******************************************
  ! RETURN DENSITY MATRIX TO GAUSSIAN VIA RWF
  !*******************************************
  
  ! Transform density matrix D 
  ! to lower triangular form
  ! and return to Gaussian via RWF
  acount = 1
  DO i=1,NBasis
     DO j=1,i
        if(FMixing)then
           MT(1,acount)=DevDensMat(1,i,j)
           IF(UHF) MT(2,acount) = DevDensMat(2,i,j)
        else
           MT(1,acount)=D(1,i,j)
           IF(UHF) MT(2,acount) = D(2,i,j)
        end if
        acount = acount +1
     ENDDO
  ENDDO
  CALL FileIO(one,-IRwPA,ntt,MT(1,:),zero)
  IF (UHF) CALL FileIO(one,-IRwPB,ntt,MT(2,:),zero)
  
  ! In last cycle (1000): clean up and return to Gaussian 
  ! without returning density matrix
  IF (jcycle.EQ.1000) THEN
     CLOSE(ifu_ham)
     CLOSE(ifu_tra)
     CLOSE(ifu_dos)
     CLOSE(ifu_mul)
     CLOSE(ifu_red)
     
#ifdef G03ROOT
     print *, " ----------------- End of ANT.G03 ------------------- "
#endif
#ifdef G09ROOT
     print *, " ----------------- End of ANT.G09 ------------------- "
#endif
     DEALLOCATE( MT, F, STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Deallocation error for MT(:,:),F(:,:,:)"
        STOP
     END IF
     IF(HFEnergy) THEN
       DEALLOCATE( CT, HC, STAT=AllocErr )
       IF( AllocErr /= 0 ) THEN
          PRINT *, "ESF/Deallocation error for CT(:,:),HC(:,:,:)"
          STOP
       END IF     
     ENDIF
     if(.not.FMixing)then
        DEALLOCATE( D, DD, DDD, STAT=AllocErr )
     IF( AllocErr /= 0 ) THEN
        PRINT *, "ESF/Deallocation error for D(:,:,:), DD(:,:,:,:), DDD(:,:,:)"
        STOP
     END IF
     end if
     CALL CleanUpDevice
     RETURN
  END IF
  
  !*********
  ! BYE-BYE
  !*********
  RETURN
END SUBROUTINE ANT
