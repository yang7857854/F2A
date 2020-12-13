SUBROUTINE USER_FORCE(MODE,I_CONTROL,R_CONTROL,NSTRUC,TIME,TIMESTEP,STAGE, &
                       POSITION,VELOCITY,COG, &
                       FORCE,ADDMASS,ERRORFLAG)

!DECLARATION TO MAKE USER_FORCE PUBLIC WITH UN-MANGLED NAME 

!DEC$ attributes dllexport , STDCALL , ALIAS : "USER_FORCE" :: user_force

!DEC$ ATTRIBUTES REFERENCE :: I_CONTROL, R_CONTROL
!DEC$ ATTRIBUTES REFERENCE :: POSITION, VELOCITY, COG, FORCE, ADDMASS
!DEC$ ATTRIBUTES REFERENCE :: MODE, NSTRUC, TIME, TIMESTEP, STAGE
!DEC$ ATTRIBUTES REFERENCE :: ERRORFLAG
!***************************************************************
! This DLL is called by AQWA for external force calculation
! It is developed by Yang Yang(a PDRA at LJMU) on 18/03/2019
!***************************************************************
! This DLL is revised to employ Aerodyn v15 for aerodynamic force calculation.
! Revision date: 01-05-2019.
! By: Yang Yang, PhD from Unversity of Shanghai for Science and Technology.
! All copyrights reserved. For any kind of use please contact Y.Yang for authorization.
! Email: 15216702797@163.com .or. y.yang@ljmu.ac.uk
!***************************************************************
! Revision for considering overhang of rotor
! Date: 22-10-2019
! By: Yang Yang
!***************************************************************
! Revision for coupling with FAST
! Date: 22-10-2019
! By: Yang Yang
!***************************************************************
! Revision: The DLL is modified to incoprate FAST with AQWA
! By: Yang Yang
! Date: 10-March-2020
!***************************************************************

USE Precision
USE AQWAinFAST
USE TurbConf, ONLY: rZT0zt,PtfmRef     ! distance between tower-base and platform reference point
USE MassInert,ONLY: TwrMass,RotMass,NacMass,NacYIner,TurbMass
USE FASTSubs

IMPLICIT NONE


INTEGER MODE, NSTRUC, STAGE, ERRORFLAG
REAL TIME, TIMESTEP
INTEGER, DIMENSION (100) :: I_CONTROL
REAL, DIMENSION (100) :: R_CONTROL
REAL, DIMENSION (3,NSTRUC) :: COG
REAL, DIMENSION (6,NSTRUC) :: POSITION, VELOCITY, FORCE
REAL, DIMENSION (6,6,NSTRUC) :: ADDMASS
!
INTEGER, DIMENSION (2) :: INITP

!
! Input Parameter Description:
!
! MODE(Int)     - 0 = Initialisation. This routine is called once with mode 0
!                     before the simulation. All parameters are as described
!                     below except for STAGE, which is undefined. FORCES and
!                     ADDMASS are assumed undefined on exit.
!                     IERR if set to > 0 on exit will cause
!                     the simulation to stop.
!
!                 1 = Called during the simulation. FORCE/ADDMASS output expected.
! 
!                99 = Termination. This routine is called once with mode 99
!                     at the end of the simulation.
!
! I_CONTROL(100)- User-defined integer control parameters input in .DAT file.
! R_CONTROL(100)- User-defined real control parameters input in .DAT file.
!
! NSTRUC(Int)   - Number of structures in the simulation
!
! TIME          - The current time (see STAGE below)
!
! TIMESTEP      - The timestep size (DT, see STAGE below)
!
! STAGE(Int)    - The stage of the integration scheme. AQWA time integration is
!                 based on a 2-stage predictor corrector method. This routine is
!                 therefore called twice at each timestep, once with STAGE=1 and
!                 once with STAGE=2. On stage 2 the position and velocity are
!                 predictions of the position and velocity at TIME+DT. 
!                 e.g. if the initial time is 0.0 and the step 1.0 seconds then
!                 calls are as follows for the 1st 3 integration steps:
!
!                 CALL USER_FORCE(.....,TIME=0.0,TIMESTEP=1.0,STAGE=1 ...)
!                 CALL USER_FORCE(.....,TIME=0.0,TIMESTEP=1.0,STAGE=2 ...)
!                 CALL USER_FORCE(.....,TIME=1.0,TIMESTEP=1.0,STAGE=1 ...)
!                 CALL USER_FORCE(.....,TIME=1.0,TIMESTEP=1.0,STAGE=2 ...)
!                 CALL USER_FORCE(.....,TIME=2.0,TIMESTEP=1.0,STAGE=1 ...)
!                 CALL USER_FORCE(.....,TIME=2.0,TIMESTEP=1.0,STAGE=2 ...)
!
! COG(3,NSTRUC) - Position of the Centre of Gravity in the Definition axes. 
!
! POSITION(6,NSTRUC) - Position of the structure in the FRA - angles in radians
!
! VELOCITY(6,NSTRUC) - Velocity of the structure in the FRA 
!                      angular velocity in rad/s
!
!
! Output Parameter Description:
!
! FORCE(6,NSTRUC) - Force on the Centre of gravity of the structure. NB: these
!                   forces are applied in the Fixed Reference axis e.g.
!                   the surge(X) force is ALWAYS IN THE SAME DIRECTION i.e. in
!                   the direction of the X fixed reference axis.
!
! ADDMASS(6,6,NSTRUC)
!                 - Added mass matrix for each structure. As the value of the
!                   acceleration is dependent on FORCES, this matrix may be used
!                   to apply inertia type forces to the structure. This mass
!                   will be added to the total added mass of the structure at
!                   each timestep at each stage.
!
! ERRORFLAG       - Error flag. The program will abort at any time if this
!                   error flag is non-zero. The values of the error flag will
!                   be output in the abort message.

!------------------------------------------------------------------------
! MODE=0 - Initialise any summing variables/open/create files.
!          This mode is executed once before the simulation begins.
!------------------------------------------------------------------------
! Define variables for the calculation
! UnIO variables
  integer,save :: UnPIn   = 10001                           ! I/O unit of the primary input file (InputFileForFAST2AQWA.txt).
  integer,save :: UnFstIn = 10002                           ! I/O unit of the FAST primary input file (xxx.fst).
  integer,save :: UnOuTF  = 10003                           ! I/O unit of the output file containing the platform loads
  integer,save :: UnOuTM  = 10004                           ! I/O unit of the output file containing the platform motions
  integer,save :: UnOuTD  = 10005                           ! I/O unit of the output file containing the relative displacement bwtween orignial point and CoG
  integer,save :: UnOuTAM  = 10006                           ! I/O unit of the output file containing the added mass due to the turbine
! Constants
!  real,parameter :: PI      = 3.1415926                     ! Parameter PI
!  real,parameter :: R2D     = 57.295780                     ! Factor to convert radians to degrees.
!  real,parameter :: RPS2RPM = 9.5492966                     ! Factor to convert radians per second(rps) to revolutions per minute (rpm).
! Input variables defined in InputFileForFAST2AQWA.txt
  character (len=1024) :: Prifile       ! - Primary input file for FAST
  logical,save         :: CouplingFlag  ! - Flag of the coupling interface
  logical,save         :: OPtfmFrc      ! - Flag of whether output the platform forces
  logical,save         :: OPtfmMot      ! - Flag of whether output the platform motions (CoG w.s.t. O)
  logical,save         :: OPtfmRld      ! - Flag of whether output the relative displacement at the reference point((0,0,0)) due to the rotations
  logical,save         :: OPtfmAdm      ! - Flag of whether output the added mass due to the wind turbine
  integer,save         :: IndexTwrStr    ! - Index of the structure connecting to the rotor directly.
! Local variables
  integer,save         :: RatioStep     ! Ratio of timestep in AQWA to time step in FAST
  integer,save         :: iSts          ! Status of the simulation. 0: Initilize FAST. 1: Run time marching in FAST. 2: Terminate FAST 
  integer,save         :: iStps         ! Number of steps for time marching solver = dT_AQWA/dT_FAST. In the first round, shall plus 1.
  real                 :: Mot(6)        ! Current motion of the platform
  real                 :: Vel(6)        ! Current velocity of the platform
  real,save            :: Acce(6)       ! Current Acceleration of the platform
  real,save            :: AddedPtfmLds(6)   ! Platform loads due to the turbine
  real,save            :: DT_FAST       ! - Time step in FAST (Unit: s)
  real,save            :: CoGCorr(3)    ! Corrected CoG
  real,save            :: VelLast(6)      ! Acceleration at last time step
  real                 :: deltaD(3)     ! Motion difference between the local and the orignial frames
  real                 :: TransMat(3,3) ! Transformation matrix
  integer              :: i,j           ! Indicator for calculation
  real                 :: AddedMass (6) ! Added mass due to rotor, tower and nacelle
  real                 :: Ixx,Iyy,Izz   ! Temporary inertia
  real                 :: TempVec       ! Temporary vector
  real                 :: z1(3),z2(3),z3(3)      ! Coordinate vector of platform reference frame
  real                 :: a1(3),a2(3),a3(3)      ! Coordinate vector of inertial frame
  real                 :: PosVec(3)        ! Position vector from tower base to the CoG
  real                 :: PosVecTraned(3)  ! Position vector transformed
  real                 :: VelVec(3)        ! Velocity vector
  real                 :: MomVec(3)        ! Moment vector
  
! Time and date information that will be output in the result file.
character(len=8) :: Currdate
character(len=10) :: CurrTime
character(len=50) :: Time_Date


!-------------------------
! Body of the subroutine
!-------------------------
IF (MODE.EQ.0) THEN
! Time and date of generating the current file.
   call date_and_time(DATE=Currdate,Time=CurrTime)
   Time_Date = 'at '//CurrTime(1:2)//':'//CurrTime(3:4)//':'//CurrTime(5:6)
   Time_Date = trim(Time_Date)//' on '//Currdate(7:8)//'-'//Currdate(5:6)//'-'//Currdate(1:4)
   
   INITP(1)=R_CONTROL(3)
   INITP(2)=R_CONTROL(4)

! First run of the DLL, open the input file
   open (UnPIn,file='InputFileForFAST2AQWA.txt')
    read (UnPIn,*) ! Input file of the interface FAST2AQWA that is developed by Y.Yang (PDRA in LJMU) on 10-March-2020 for performing fully-coupled analysis of floating offshore wind turbines (FOWTs).
    read (UnPIn,*) ! The FAST2AQWA is implemented through the user_force DLL. FAST v7 is used to examine the aero-elastic effects of the wind turbine.
    read (UnPIn,*) ! This file is specified for DTU 10 MW wind turbine supported by a multi-body platform (Do not remove anyline below)
    read (UnPIn,*) ! --------------- FAST Configuration ------------------
    read (UnPIn,*) Prifile       ! - Primary input file for FAST
    read (UnPIn,*) CouplingFlag  ! - Flag of the coupling interface
    read (UnPIn,*) ! --------------- AQWA structure properties -----------
    read (UnPIn,*) IndexTwrStr   ! - Index of the structure connecting to the rotor directly.
    read (UnPIn,*) OPtfmFrc      ! - Flag of whether output the platform forces
    !read (UnPIn,*) OPtfmMot      ! - Flag of whether output the platform motions (CoG w.s.t. O)
    !read (UnPIn,*) OPtfmRld      ! - Flag of whether output the relative displacement at the reference point((0,0,0)) due to the rotations
    !read (UnPIn,*) OPtfmAdm      ! - Flag of whether output the added mass due to the wind turbine
   close(UnPIn)
   open (UnFstIn, file = trim(Prifile))
    do i=1,10
       read(UnFstIn,*)
    enddo
    read(UnFstIn,*) DT_FAST
   close(UnFstIn)
   
! Write notice to the user:
    write (*,*) '---------------------------------------------------------------------------'
    write (*,*) '------------------       FAST2AQWA Interface        -----------------------'
    write (*,*) '------------------   (v1.01.01d-yy, 24-Nov-2020)    -----------------------'
    write (*,*) '---------------------------------------------------------------------------',&
                ' FAST2AQWA is developed by Dr Y.Yang ( University of Shanghai for Science  ',&
                ' and Technology) to incorporate FAST with AQWA through the user_force DLL  ',&
                ' for performing fully-coupled analysis of floating offshore wind turbines. ',&
                '---------------------------------------------------------------------------'
    write (*,*) 
    write (*,*) ' Using the FAST input file: ' ,trim(Prifile)

! Determine the variables that need to be passed into FAST
      RatioStep = NINT (TIMESTEP/(DT_FAST*2.0))    ! loop times in FAST
! Done reading.
  FAST2AQWA = CouplingFlag
! Open the file to store the forces that will be applied on the platform.
   if (OPtfmFrc) then
      open (UnOuTF,file = trim(Prifile(1:LEN(TRIM(Prifile))-4)//'_PtfmLoads.dat')) ! To store the forces will be applied on the platform.
       write(UnOuTF,'(A)') 'Forces calcuated using the DLL developed by Y.Yang (PDRA in LJMU) for AQWA earternal force calculation. Generated '//trim(Time_Date)
       write(UnOuTF,'(A10,6(A15))') 'Time','PtfmFrcX','PtfmFrcY','PtfmFrcZ','PtfmMmtX','PtfmMmtY','PtfmMmtZ'
       write(UnOuTF,'(A10,6(A15))') '(s)','(N)','(N)','(N)','(N-m)','(N-m)','(N-m)' 
   endif


continue

!------------------------------------------------------------------------
! MODE=1 - On-going - calculation of forces/mass
!------------------------------------------------------------------------
ELSEIF (MODE.EQ.1) THEN

!ERRORFLAG = 0

! Initialize force and added mass.
      FORCE = 0.0
      AddedMass = 0.0
! Set Motion
      Mot = POSITION (:,IndexTwrStr)
      Vel = VELOCITY (:,IndexTwrStr)
! Correct position vector from (Origin to CoG) TO (Origin to Reference point)
      if (FAST2AQWA) then
         PosVec = COG(:,IndexTwrStr) 
         PosVec(3) = PosVec(3) - PtfmRef
         call getTransMatEuler(TransMat, Mot(4:6))
         PosVecTraned(1) = TransMat(1,1) * PosVec(1) + TransMat(1,2) * PosVec(2) + TransMat(1,3) * PosVec(3)
         PosVecTraned(2) = TransMat(2,1) * PosVec(1) + TransMat(2,2) * PosVec(2) + TransMat(2,3) * PosVec(3)
         PosVecTraned(3) = TransMat(3,1) * PosVec(1) + TransMat(3,2) * PosVec(2) + TransMat(3,3) * PosVec(3)
         Mot(1:3) = Mot(1:3) - PosVecTraned
         call CrossProd(VelVec,Vel(4:6),PosVecTraned) ! Velocity due to rotations
         Vel(1:3) = Vel(1:3) - VelVec 
      else
         Mot  = 0.0
         Vel = 0.0
      endif
    
! Detertime status imported to FAST  
       if (TIME == 0 ) then
          iSts = 0              ! initialize FAST
          if (STAGE == 1) then
              iStps = 1
              PrtFASTRslts = .TRUE.
              VelLast = Vel
              Acce = 0.0
          else                 ! Time marching
              Acce = (Vel - VelLast)/TIMESTEP
              VelLast = Vel
              iStps = RatioStep
              PrtFASTRslts = .FALSE.
          endif
      else     ! time marching
          Acce = (Vel - VelLast)/TIMESTEP
          VelLast = Vel
          iSts = 1
          iStps = RatioStep
      endif
            
! Call FAST to march 1/2 TIMESTEP for STAGE 1 and STAGE 2, respectively
      call FAST(PriFile,iSts,iStps,Mot,Vel,Acce)
      
! Yang: convert tower-base loads to platform loads with respect to the platform reference frame
      PosVec = COG(:,IndexTwrStr)
! Position vector from tower-base to CoG      
      PosVec(3) = PosVec(3) - rZT0zt 
      AddedPtfmLds(1:3) =  PtfmFrcAQWA(1:3)
      call CrossProd (MomVec,AddedPtfmLds(1:3),PosVec)
      AddedPtfmLds(4:6) = PtfmMomAQWA + MomVec
     
! Convert the platform loads from platform reference frame to the inertial frame coordinate system
      call getTransMat (Mot(4),Mot(5),Mot(6), TransMat)
      TransMat = transpose(TransMat)
! Loads
      FORCE(1,IndexTwrStr) = DOT_PRODUCT(AddedPtfmLds(1:3),TransMat(1,:))   ! Surge force
      FORCE(2,IndexTwrStr) = DOT_PRODUCT(AddedPtfmLds(1:3),TransMat(2,:))   ! Sway  force
      FORCE(3,IndexTwrStr) = DOT_PRODUCT(AddedPtfmLds(1:3),TransMat(3,:))   ! Heave force
      FORCE(4,IndexTwrStr) = DOT_PRODUCT(AddedPtfmLds(4:6),TransMat(1,:))   ! Roll  moment
      FORCE(5,IndexTwrStr) = DOT_PRODUCT(AddedPtfmLds(4:6),TransMat(2,:))   ! Pitch moment
      FORCE(6,IndexTwrStr) = DOT_PRODUCT(AddedPtfmLds(4:6),TransMat(3,:))   ! Yaw   moment      
        

! Write the results in output
      if (STAGE == 1 .and. OPtfmFrc ) then
         write(UnOuTF,'(F10.3,6(ES15.6))') TIME,(FORCE(j,IndexTwrStr),j=1,6)       
      endif  

!------------------------------------------------------------------------
! MODE=99 - Termination - Output/print any summaries required/Close Files
!           This mode is executed once at the end of the simulation
!------------------------------------------------------------------------

ELSEIF (MODE.EQ.99) THEN
! Terminate FAST
iSts = 2
call FAST(PriFile,iSts,iStps,Mot,Vel,Acce)
! Close the files
if (OPtfmFrc) then
   close(UnOuTF)
endif

!------------------------------------------------------------------------
! MODE# ERROR - OUTPUT ERROR MESSAGE
!------------------------------------------------------------------------

ELSE	

ENDIF
RETURN

END SUBROUTINE USER_FORCE


subroutine getTransMat (Theta1,Theta2,Theta3, TransMat)
implicit none
real  :: TransMat(3,3)
real  :: Theta1,Theta2,Theta3,Theta11,Theta22,Theta33,sumSq,SqrtSumSq,ComDenom,Theta12S,Theta13S,Theta23S

Theta11 = Theta1*Theta1
Theta22 = Theta2*Theta2
Theta33 = Theta3*Theta3
sumSq = Theta11 + Theta22 + Theta33
SqrtSumSq = sqrt(sumSq+1.0)
ComDenom = sumSq*SqrtSumSq
Theta12S = Theta1*Theta2*(SqrtSumSq - 1.0)
Theta13S = Theta1*Theta3*(SqrtSumSq - 1.0)
Theta23S = Theta2*Theta3*(SqrtSumSq - 1.0)

if (ComDenom == 0) then
    TransMat(1,:) = (/ 1.0, 0.0, 0.0 /)
    TransMat(2,:) = (/ 0.0, 1.0, 0.0 /)
    TransMat(3,:) = (/ 0.0, 0.0, 1.0 /)

else
    TransMat(1,1) = (Theta11*SqrtSumSq+Theta22+Theta33)/ComDenom
    TransMat(2,2) = (Theta22*SqrtSumSq+Theta11+Theta33)/ComDenom
    TransMat(3,3) = (Theta33*SqrtSumSq+Theta11+Theta22)/ComDenom
    TransMat(1,2) = ( Theta3*sumSq + Theta12S)/ComDenom
    TransMat(2,1) = (-Theta3*sumSq + Theta12S)/ComDenom
    TransMat(1,3) = (-Theta2*sumSq + Theta13S)/ComDenom
    TransMat(3,1) = ( Theta2*sumSq + Theta13S)/ComDenom
    TransMat(2,3) = ( Theta1*sumSq + Theta23S)/ComDenom
    TransMat(3,2) = (-Theta1*sumSq + Theta23S)/ComDenom
endif
end subroutine getTransMat
!===============================================
SUBROUTINE   getTransMatEuler(TransMat, Rots)

implicit none

  real :: TransMat (3,3)       ! Transformation matrix
  real :: Rots (3)             ! Rotations
  real :: cx,cy,cz             ! Cosines
  real :: sx,sy,sz             ! Sines


  cx = cos(Rots(1));
  sx = sin(Rots(1));
  cy = cos(Rots(2));
  sy = sin(Rots(2));
  cz = cos(Rots(3));
  sz = sin(Rots(3));
  
  TransMat(1,:) = (/ cz*cy, -sz*cx+cz*sy*sx, sz*sx+cz*sy*cx /)
  TransMat(2,:) = (/sz*cy,  cz*cx+sz*sy*sx, -cz*sx+sz*sy*cx /)
  TransMat(3,:) = (/ -sy  ,   cy*sx ,     cy*cx       /)
  
  END SUBROUTINE getTransMatEuler
! ========================================================================================