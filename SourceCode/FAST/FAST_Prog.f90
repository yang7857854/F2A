!=======================================================================
SUBROUTINE FAST (INPUTFILE,iSts,iStps,PtfmMotIn,PtfmVelIn,PtfmAcceIn)


   ! This program models 2- or 3-bladed turbines of a standard configuration.


USE                             ADAMSInput
USE                             DOFs
USE                             General
USE                             InitCond
USE                             Linear
USE                             SimCont
USE                             NWTC_Library
USE                             AeroDyn
USE                             FAST_IO_Subs       ! FAST_Begin(), FAST_Input(), PrintSum(), RunTimes()
USE                             FASTsubs           ! FAST_Initialize(), TimeMarch()
USE                             FAST2ADAMSSubs     ! MakeAdm(), MakeACF(), MakeACF_Lin
USE                             FAST_Lin_Subs      ! CalcSteady(), Linearize()
USE                             HydroDyn
USE                             Noise

! Yang: USE AQWA & Output modules for the variables
USE                             AQWAinFAST
USE                             Output
USE                             TurbConf
! Yang: End of the change

IMPLICIT                        NONE


   ! Local variables:

!REAL(ReKi)                   :: GBoxTrq                                         ! Unused gearbox torque on the LSS side in N-m.

CHARACTER (LEN=1024)         :: INPUTFILE
!INTEGER(4)                   :: L                                              ! Generic index.
INTEGER(4)                   :: iSts                                           ! status of calling FAST. 0 = first,1 = the rest, 2 = the end.
INTEGER(4)                   :: iStps                                          ! number of steps for time marching solver.= dT_AQWA/dT_FAST. In the first round, shall plus 1.

INTEGER                       :: ErrStat

! Yang: define the platform motion and velocity matrix
REAL (Reki)           :: PtfmMotIn(6),PtfmVelIn(6),PtfmAcceIn(6)                    ! Platfrom motion, velocity and acceleration in the six DoFs
! End of the change

! Set Prifile
if (iSts == 0 .and. PrtFASTRslts) then      ! adding the condition (PrtFASTRslts) in this statement is just to avoid the second initilization when stage = 2 and time = 0
    Prifile = INPUTFILE
! Get the current time.
    
    CALL DATE_AND_TIME ( Values=StrtTime )                                           ! Let's time the whole simulation
    CALL CPU_TIME ( UsrTime1 )                                                       ! Initial time (this zeros the start time when used as a MATLAB function)
    
! Set version & initialize NWTC Library (open console, set pi constants)
    ! Yang: Not output FAST notice when couping with AQWA
    CALL SetVersion
    CALL NWTC_Init()                                                                 ! sets the pi constants
    
       ! Tell our nice users what they're running.
    ! Yang: Not output FAST notice when couping with AQWA
    CALL DispNVD()
    
       ! Open and read input files, initialize global parameters.
    
    CALL FAST_Begin()
    CALL FAST_Input()

   ! Set up initial values for all degrees of freedom.
    
    CALL FAST_Initialize()

   ! Print summary information to "*.fsm"?

    IF ( SumPrint )  CALL PrintSum
endif
! Yang: Set values for Q and QD with regards to the platform DoFs
IF (iSts == 0 .or. iSts == 1)  then
   PtfmMotAQWA   = PtfmMotIn
   PtfmVelAQWA   = PtfmVelIn
   PtfmAcceAQWA   = PtfmAcceIn
  ! write(*,*) 'Running TimeMarch'
   CALL TimeMarch (iSts,iStps)
! Rotor, Nacelle, Tower mass center

ENDIF

if (iSts == 2) then
   CALL RunTimes
   CALL FAST_Terminate( ErrStat )
   CALL AD_Terminate(   ErrStat )
   CALL Hydro_Terminate( )
   CALL Noise_Terminate( )
   IF ( BEEP ) CALL UsrAlarm
   CALL NormStop( )
endif

END SUBROUTINE FAST
!=======================================================================


