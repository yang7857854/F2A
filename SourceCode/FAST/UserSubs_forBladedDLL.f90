   ! NOTE: This source file contains dummy placeholders for ALL of the
   !       user-specified routines available in FAST.  These routines
   !       are as follows:
   !          Routine       Description
   !          ------------  ---------------------------------------------------
   !          PitchCntrl()  User-specified blade pitch control (either
   !                        independent or rotor-collective) model.
   !          UserGen()     User-specified generator torque and power model.
   !          UserHSSBr()   User-specified high-speed shaft brake model.
   !          UserPtfmLd()  User-specified platform loading model.
   !          UserRFrl()    User-specified rotor-furl spring/damper model.
   !          UserTeet()    User-specified rotor-teeter spring/damper model.
   !          UserTFin()    User-specified tail fin aerodynamics model.
   !          UserTFrl()    User-specified tail-furl spring/damper model.
   !          UserVSCont()  User-specified variable-speed torque and power
   !                        control model.
   !          UserYawCont() User-specified nacelle-yaw control model.
   !       In order to interface FAST with your own user-specified routines,
   !       you can develop your own logic within these dummy placeholders and
   !       recompile FAST; OR comment out the appropriate dummy placeholders,
   !       create your own routines in their own source files, and recompile
   !       FAST while linking in these additional source files.  For example,
   !       the executable version of FAST that is distributed with the FAST
   !       archive is linked with the example PitchCntrl() routine contained in
   !       source file PitchCntrl_ACH.f90 and the example UserGen() and
   !       UserVSCont() routines contained in source file UserVSCont_KP.f90;
   !       thus, the dummy placeholders for routines PitchCntrl(), UserGen(),
   !       and UserVSCont() are commented out within this source file.  The
   !       example pitch controller was written by Craig Hansen (ACH) and the
   !       example generator and variable speed controllers were written by
   !       Kirk Pierce (KP).  Please see the aforementioned source files for
   !       additional information on these example user-specified routines.

   ! NOTE: If you (the user) wants to access the current value of ANY of the
   !       output parameters available as outputs from FAST from your
   !       user-defined routines, then do the following:
   !          (1) USE MODULE Output() in your routine.
   !          (2) Access the output parameter by typing "AllOuts(OutName)",
   !              where OutName is the PRIMARY name of the output parameter.
   !              For example, to access the current value of the in-plane
   !              bending moment at the root of blade 1 (in kN·m), type in
   !              "AllOuts(RootMxc1)", since RootMxc1 is the primary name of
   !              this output parameter--RootMIP1 will not work in place of
   !              RootMxc1, since it is a SECONDARY name.  Also, you CANNOT use
   !              the prefixes ("-", "_", "m", or "M") in front of OutName to
   !              reverse the sign of the selected output channel.
   !       Note that OutName DOES NOT have to be one of the output parameters
   !       you listed in OutList from your primary input file.  Also note that
   !       this technique WILL also work for user-defined routines written for
   !       ADAMS datasets extracted using the FAST-to-ADAMS preprocessor.

!=======================================================================
!SUBROUTINE PitchCntrl ( BlPitch, ElecPwr, HSS_Spd, GBRatio, TwrAccel, NumBl, ZTime, DT, DirRoot, BlPitchCom )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! blade pitch control model (either independent or rotor-collective).
!   ! Modify this code to create your own model.
!
!
!USE                             Precision
!
!
!IMPLICIT                        NONE
!
!
!   ! Passed variables:
!
!INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
!
!REAL(ReKi), INTENT(IN )      :: BlPitch   (NumBl)                               ! Current values of the blade pitch angles, rad.
!REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
!REAL(ReKi), INTENT(IN )      :: ElecPwr                                         ! Electrical power, watts.
!REAL(ReKi), INTENT(IN )      :: GBRatio                                         ! Gearbox ratio, (-).
!REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
!REAL(ReKi), INTENT(OUT)      :: BlPitchCom(NumBl)                               ! Commanded blade pitch angles (demand pitch angles), rad.
!REAL(ReKi), INTENT(IN )      :: TwrAccel                                        ! Tower Acceleration, m/s^2.
!REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!
!BlPitchCom = 0.0
!
!
!
!RETURN
!END SUBROUTINE PitchCntrl
!=======================================================================
!SUBROUTINE UserGen ( HSS_Spd, GBRatio, NumBl, ZTime, DT, GenEff, DelGenTrq, DirRoot, GenTrq, ElecPwr )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! generator torque and power model.  Modify this code to create your
!   ! own model.
!
!   ! NOTE: If you (the user) wants to switch on-or-off the generator DOF at
!   !       runtime from this user-defined routine, then do the following:
!   !          (1) USE MODULE DOFs().
!   !          (2) Type in "DOF_Flag(DOF_GeAz) = VALUE" where VALUE = .TRUE. or
!   !              .FALSE. depending on whether you want to turn-on or turn-off
!   !              the DOF, respectively.  Turning off the DOF forces the
!   !              current RATE to remain fixed.  If the rate is currently zero,
!   !              the current POSITION will remain fixed as well.
!   !       Note that this technique WILL NOT work for user-defined routines
!   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
!   !       preprocessor.
!
!
!USE                             Precision
!
!
!IMPLICIT                        NONE
!
!
!   ! Passed Variables:
!
!INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
!
!REAL(ReKi), INTENT(IN )      :: DelGenTrq                                       ! Pertubation in generator torque used during FAST linearization (zero otherwise), N-m.
!REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
!REAL(ReKi), INTENT(OUT)      :: ElecPwr                                         ! Electrical power (account for losses), watts.
!REAL(ReKi), INTENT(IN )      :: GBRatio                                         ! Gearbox ratio, (-).
!REAL(ReKi), INTENT(IN )      :: GenEff                                          ! Generator efficiency, (-).
!REAL(ReKi), INTENT(OUT)      :: GenTrq                                          ! Electrical generator torque, N-m.
!REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
!REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!
!GenTrq  = 0.0 + DelGenTrq  ! Make sure to add the pertubation on generator torque, DelGenTrq.  This is used only for FAST linearization (it is zero otherwise).
!
!
!   ! The generator efficiency is either additive for motoring,
!   !   or subtractive for generating power.
!
!IF ( GenTrq > 0.0 )  THEN
!   ElecPwr = GenTrq*HSS_Spd*GenEff
!ELSE
!   ElecPwr = GenTrq*HSS_Spd/GenEff
!ENDIF
!
!
!
!RETURN
!END SUBROUTINE UserGen
!=======================================================================

!=======================================================================
!SUBROUTINE UserPtfmLd ( X, XD, ZTime, DirRoot, PtfmAM, PtfmFt )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! platform loading model.  Modify this code to create your own model.
!   ! The local variables and associated calculations below provide a
!   ! template for making this user-specified platform loading model
!   ! include linear 6x6 damping and stiffness matrices.  These are
!   ! provided as an example only and can be modified or deleted as
!   ! desired by the user without detriment to the interface (i.e., they
!   ! are not necessary for the interface).
!
!   ! The platform loads returned by this routine should contain contributions
!   !   from any external load acting on the platform other than loads
!   !   transmitted from the wind turbine.  For example, these loads should
!   !   contain contributions from foundation stiffness and damping [not
!   !   floating] or mooring line restoring and damping [floating], as well as
!   !   hydrostatic and hydrodynamic contributions [offshore].  The platform
!   !   loads will be applied on the platform at the instantaneous platform
!   !   reference position within FAST and ADAMS.
!
!   ! This routine assumes that the platform loads are transmitted through a
!   !   medium like soil [foundation] and/or water [offshore], so that added
!   !   mass effects are important.  Consequently, the routine assumes that the
!   !   total platform load can be written as:
!   !
!   ! PtfmF(i) = SUM( -PtfmAM(i,j)*XDD(j), j=1,2,..,6) + PtfmFt(i) for i=1,2,...,6
!   !
!   ! where,
!   !   PtfmF(i)    = the i'th component of the total load applied on the
!   !                 platform; positive in the direction of positive motion of
!   !                 the i'th DOF of the platform
!   !   PtfmAM(i,j) = the (i,j) component of the platform added mass matrix
!   !                 (output by this routine)
!   !   XDD(j)      = the j'th component of the platform acceleration vector
!   !   PtfmFt(i)   = the i'th component of the portion of the platform load
!   !                 associated with everything but the added mass effects;
!   !                 positive in the direction of positive motion of the i'th
!   !                 DOF of the platform (output by this routine)
!
!   ! The order of indices in all arrays passed to and from this routine is as
!   !   follows:
!   !      1 = Platform surge / xi-component of platform translation (internal DOF index = DOF_Sg)
!   !      3 = Platform sway  / yi-component of platform translation (internal DOF index = DOF_Sw)
!   !      3 = Platform heave / zi-component of platform translation (internal DOF index = DOF_Hv)
!   !      4 = Platform roll  / xi-component of platform rotation    (internal DOF index = DOF_R )
!   !      5 = Platform pitch / yi-component of platform rotation    (internal DOF index = DOF_P )
!   !      6 = Platform yaw   / zi-component of platform rotation    (internal DOF index = DOF_Y )
!
!   ! NOTE: The added mass matrix returned by this routine, PtfmAM, must be
!   !       symmetric.  FAST and ADAMS will abort otherwise.
!   !
!   !       Please also note that the hydrostatic restoring contribution to the
!   !       hydrodynamic force returned by this routine should not contain the
!   !       effects of body weight, as is often done in classical marine
!   !       hydrodynamics.  The effects of body weight are included within FAST
!   !       and ADAMS.
!
!
!USE                             Precision
!
!
!IMPLICIT                        NONE
!
!
!   ! Passed Variables:
!
!REAL(ReKi), INTENT(OUT)      :: PtfmAM (6,6)                                    ! Platform added mass matrix, kg, kg-m, kg-m^2.
!REAL(ReKi), INTENT(OUT)      :: PtfmFt   (6)                                    ! The 3 components of the portion of the platform force (in N  ) acting at the platform reference and the 3 components of the portion of the platform moment (in N-m  ) acting at the platform reference associated with everything but the added-mass effects; positive forces are in the direction of motion.
!REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement    (in m  )        of the platform reference and the 3 components of the rotational displacement        (in rad  )        of the platform relative to the inertial frame.
!REAL(ReKi), INTENT(IN )      :: XD       (6)                                    ! The 3 components of the translational velocity        (in m/s)        of the platform reference and the 3 components of the rotational (angular) velocity  (in rad/s)        of the platform relative to the inertial frame.
!REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!   ! Local Variables:
!
!REAL(ReKi)                   :: Damp   (6,6)                                    ! Damping matrix.
!REAL(ReKi)                   :: Stff   (6,6)                                    ! Stiffness/restoring matrix.
!
!INTEGER(4)                   :: I                                               ! Generic index.
!INTEGER(4)                   :: J                                               ! Generic index.
!
!
!
!Damp  (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Damp  (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Damp  (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Damp  (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Damp  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Damp  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!
!Stff  (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Stff  (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Stff  (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Stff  (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Stff  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!Stff  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!
!PtfmAM(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!PtfmAM(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!PtfmAM(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!PtfmAM(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!PtfmAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!PtfmAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!
!PtfmFt(1)   = 0.0
!PtfmFt(2)   = 0.0
!PtfmFt(3)   = 0.0
!PtfmFt(4)   = 0.0
!PtfmFt(5)   = 0.0
!PtfmFt(6)   = 0.0
!
!DO J = 1,6
!   DO I = 1,6
!      PtfmFt(I) = PtfmFt(I) - Damp(I,J)*XD(J) - Stff(I,J)*X(J)
!   ENDDO
!ENDDO
!
!
!
!RETURN
!END SUBROUTINE UserPtfmLd
!=======================================================================
!=======================================================================
! Below is a replacement subroutine for UserPtfmLd that provides the 
! ability to conduct simulations of earthquake loading in FAST.
!=======================================================================

!=======================================================================
SUBROUTINE UserPtfmLd ( X, XD, ZTime, DirRoot, PtfmAM, PtfmFt )
! Yang: Start of modifying the UserPtfmLd subroutine.   
! Yang: This subroutine has been modified to consider soil effects.
! Date: 2018-05-21
      USE                             Precision
      USE                             Constants
      !USE                             SSI                                         ! Use SSI module to access the SSI parameters
      USE                             NWTC_IO
      !
      IMPLICIT                        NONE
      
      
      ! Passed Variables:
      
      REAL(ReKi), INTENT(OUT)      :: PtfmAM (6,6)                                    ! Platform added mass matrix, kg, kg-m, kg-m^2.
      REAL(ReKi), INTENT(OUT)      :: PtfmFt   (6)                                    ! The 3 components of the portion of the platform force (in N  ) acting at the platform reference and the 3 components of the portion of the platform moment (in N-m  ) acting at the platform reference associated with everything but the added-mass effects; positive forces are in the direction of motion.
      REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement    (in m  )        of the platform reference and the 3 components of the rotational displacement        (in rad  )        of the platform relative to the inertial frame.
      REAL(ReKi), INTENT(IN )      :: XD       (6)                                    ! The 3 components of the translational velocity        (in m/s)        of the platform reference and the 3 components of the rotational (angular) velocity  (in rad/s)        of the platform relative to the inertial frame.
      REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
      CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
      !
      !! Local Variables:
      !REAL(ReKi),dimension(6)      :: kTrans = 0.0                                    ! The translational stiffness of soil effect
      !REAL(ReKi),dimension(6)      :: cTrans = 0.0                                    ! The translational damping of soil effect
      !REAL(ReKi),dimension(6)      :: mTrans = 0.0                                    ! The mass properties of the soil effect
      !! Yang: modification
      !
      !REAL(ReKi)                   :: Damp   (6,6)                                    ! Damping matrix.
      !REAL(ReKi)                   :: Stff   (6,6)                                    ! Stiffness/restoring matrix.
      !REAL(ReKi)                   :: CsWaveVel                                       ! Wave velocity
      !REAL(ReKi),dimension(6)      :: DimlessGama                                     ! Dimensionless dashpot gama
      !REAL(ReKi),dimension(6)      :: DimlessMue                                      ! Coefficient of mass Mue
      !INTEGER(4)                   :: i,j,k                                         ! Generic index.
      !LOGICAL (1)                  :: ex
      !LOGICAL (1)                  :: PrintInfo = .TRUE.
      !
! Init!ializing the PtfmFt matrix
      !
      !PtfmFt(1)   = 0.0
      !PtfmFt(2)   = 0.0
      !PtfmFt(3)   = 0.0
      !PtfmFt(4)   = 0.0
      !PtfmFt(5)   = 0.0
      !PtfmFt(6)   = 0.0
      !
! Soil! effect
      !  IF (SSIMode == 1)    THEN          !  CS model£¨Wolf method) 
      !     IF (PrintInfo) THEN
      !        CALL WrScr1 ( 'The SSI is modelled using Wolf method.')
      !        PrintInfo = .FALSE.
      !     ENDIF
! Calc!ulate relative stiffness and damping properties.
      !     CsWaveVel = SQRT(ShrModuSoil*1.0e6/DenSoil)
      !     kTrans(1) = 8.0*ShrModuSoil*1.0e6*RadiusPtfm/(2-PoissRatio)
      !     kTrans(2) = kTrans(1)
      !     kTrans(3) = 4.0*ShrModuSoil*1.0e6*RadiusPtfm/(1-PoissRatio)
      !     kTrans(4) = 8.0*ShrModuSoil*1.0e6*RadiusPtfm**3/(3*(1-PoissRatio)) 
      !     kTrans(5) = kTrans(4)
      !     kTrans(6) = 16.0/3*ShrModuSoil*1.0e6*RadiusPtfm**3
      !     DimlessGama(:) = (/ 0.58, 0.58, 0.85, 0.15, 0.15, 0.21 /)
      !     DimlessMue(:) = (/ 0.095, 0.095, 0.27, 0.24, 0.24, 0.045 /)
      !     DO i=1,6
      !        cTrans(i) = RadiusPtfm/CsWaveVel*DimlessGama(i)*kTrans(i)
      !        mTrans(i) = (RadiusPtfm/CsWaveVel)**2*DimlessMue(i)*kTrans(i)
      !     ENDDO
      !     Damp  (1,:) = (/ cTrans(1), 0.0, 0.0, 0.0, 0.0, 0.0 /)
      !     Damp  (2,:) = (/ 0.0, cTrans(2), 0.0, 0.0, 0.0, 0.0 /)
      !     Damp  (3,:) = (/ 0.0, 0.0, cTrans(3), 0.0, 0.0, 0.0 /)
      !     Damp  (4,:) = (/ 0.0, 0.0, 0.0, cTrans(4), 0.0, 0.0 /)
      !     Damp  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, cTrans(5), 0.0 /)
      !     Damp  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, cTrans(6) /)
      !     
      !     Stff  (1,:) = (/ kTrans(1), 0.0, 0.0, 0.0, 0.0, 0.0 /)
      !     Stff  (2,:) = (/ 0.0, kTrans(2), 0.0, 0.0, 0.0, 0.0 /)
      !     Stff  (3,:) = (/ 0.0, 0.0, kTrans(3), 0.0, 0.0, 0.0 /)
      !     Stff  (4,:) = (/ 0.0, 0.0, 0.0, kTrans(4), 0.0, 0.0 /)
      !     Stff  (5,:) = (/ 0.0, 0.0, 0.0, 0.0, kTrans(5), 0.0 /)
      !     Stff  (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, kTrans(6) /)
      !     
      !     PtfmAM(1,:) = (/ mTrans(1), 0.0, 0.0, 0.0, 0.0, 0.0 /)
      !     PtfmAM(2,:) = (/ 0.0, mTrans(2), 0.0, 0.0, 0.0, 0.0 /)
      !     PtfmAM(3,:) = (/ 0.0, 0.0, mTrans(3), 0.0, 0.0, 0.0 /)
      !     PtfmAM(4,:) = (/ 0.0, 0.0, 0.0, mTrans(4), 0.0, 0.0 /)
      !     PtfmAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, mTrans(5), 0.0 /)
      !     PtfmAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, mTrans(6) /)
      !     
      ! ELSEIF (SSIMode == 2)       THEN  ! Specified matrices for the stiffness, damping and added-mass properties
      !     IF (PrintInfo) THEN
      !        CALL WrScr1 ( 'The SSI is modelled using specified spring and dashpot properties.')
      !        PrintInfo = .FALSE.
      !     ENDIF
      !     INQUIRE (file=TRIM(InFileGKCM), exist=ex)
      !     IF (ex .eq. .FALSE.) THEN
      !         CALL ProgAbort( 'The File'// TRIM(InFileGKCM)// ' which contains the specified matrices cannot be found' )
      !     ELSE
      !        OPEN (1001,file = TRIM(InFileGKCM))
      !         !skip the comments 
      !         DO j=1,4 
      !            READ (1001,*)
      !         ENDDO
      !         DO j=1,6
      !            READ (1001,*) (Stff(j,k),k=1,6)
      !         ENDDO
      !         READ (1001,*)
      !         DO j=1,6
      !            READ (1001,*) (Damp(j,k),k=1,6)
      !         enddo
      !         READ (1001,*)
      !         DO j=1,6
      !            READ (1001,*) (PtfmAM(j,k),k=1,6)
      !         ENDDO
      !        CLOSE (1001)
      !     ENDIF
      !  ELSE   ! SSI mode is none or DS, platform will be fixed.
      !     Damp = 0.0 
      !     Stff = 0.0
      !     PtfmAM = 0.0
      !  ENDIF 
  ! ca!lculate the force 
      !  DO j = 1,6
      !     DO i = 1,6
      !        PtfmFt(i) = PtfmFt(i) - Damp(i,j)*XD(j) - Stff(i,j)*X(j)
      !     ENDDO
      !  ENDDO
RETURN
      ! Yang: End of modifying the UserTwrLd subroutine
 END SUBROUTINE UserPtfmLd

!=======================================================================
! This is the end of the seismic UserPtfmLd
!=======================================================================
!JASON: MOVE THIS USER-DEFINED ROUTINE (UserTwrLd) TO THE UserSubs.f90 OF FAST WHEN THE MONOPILE LOADING FUNCTIONALITY HAS BEEN DOCUMENTED!!!!!
SUBROUTINE UserTwrLd ( JNode, X, XD, ZTime, DirRoot, TwrAM, TwrFt )


   ! This is a dummy routine for holding the place of a user-specified
   ! tower loading model.  Modify this code to create your own model.
   ! The local variables and associated calculations below provide a
   ! template for making this user-specified tower loading model
   ! include linear 6x6 damping and stiffness matrices.  These are
   ! provided as an example only and can be modified or deleted as
   ! desired by the user without detriment to the interface (i.e., they
   ! are not necessary for the interface).

   ! The tower loads returned by this routine should contain contributions from
   !   any external load acting on the current tower element (indicated by
   !   JNode) other than loads transmitted from tower aerodynamics.  For
   !   example, these tower forces should contain contributions from foundation
   !   stiffness and damping [not floating] or mooring line/guy wire restoring
   !   and damping, as well as hydrostatic and hydrodynamic contributions
   !   [offshore].

   ! This routine assumes that the tower loads are transmitted through a medium
   !   like soil [foundation] and/or water [offshore], so that added mass
   !   effects are important.  Consequently, the routine assumes that the total
   !   load per unit length on the current tower element can be written as:
   !
   ! TwrF(i) = SUM( -TwrAM(i,j)*XDD(j), j=1,2,..,6) + TwrFt(i) for i=1,2,...,6
   !
   ! where,
   !   TwrF(i)    = the i'th component of the total load per unit length
   !                applied on the current tower element; positive in the
   !                direction of positive motion of the i'th DOF of the current
   !                tower element
   !   TwrAM(i,j) = the (i,j) component of the tower added mass matrix per unit
   !                length (output by this routine)
   !   XDD(j)     = the j'th component of the current tower element
   !                acceleration vector
   !   TwrFt(i)   = the i'th component of the portion of the current tower
   !                element load per unit length associated with everything but
   !                the added mass effects; positive in the direction of
   !                positive motion of the i'th DOF of the current tower
   !                element (output by this routine)

   ! The order of indices in all arrays passed to and from this routine is as
   !   follows:
   !      1 = Current tower element surge / xi-component of translation
   !      3 = Current tower element sway  / yi-component of translation
   !      3 = Current tower element heave / zi-component of translation
   !      4 = Current tower element roll  / xi-component of rotation
   !      5 = Current tower element pitch / yi-component of rotation
   !      6 = Current tower element yaw   / zi-component of rotation

   ! NOTE: The added mass matrix returned by this routine, TwrAM, must be
   !       symmetric.  FAST and ADAMS will abort otherwise.
   !
   !       Please also note that the hydrostatic restoring contribution to the
   !       hydrodynamic force returned by this routine should not contain the
   !       effects of body weight, as is often done in classical marine
   !       hydrodynamics.  The effects of body weight are included within FAST
   !       and ADAMS.

USE                             NWTC_Library

IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(OUT)      :: TwrAM  (6,6)                                    ! Added mass matrix per unit length of current tower element, kg/m, kg-m/m, kg-m^2/m.
REAL(ReKi), INTENT(OUT)      :: TwrFt    (6)                                    ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force per unit length (in N/m) at the current tower element and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the tower moment per unit length (in N-m/m) acting at the current tower element associated with everything but the added-mass effects; positive forces are in the direction of motion.
REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement (in m  ) of the current tower node and the 3 components of the rotational displacement       (in rad  ) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
REAL(ReKi), INTENT(IN )      :: XD       (6)                                    ! The 3 components of the translational velocity     (in m/s) of the current tower node and the 3 components of the rotational (angular) velocity (in rad/s) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

INTEGER(4), INTENT(IN )      :: JNode                                           ! The number of the current tower node / element, (-). [1 to TwrNodes]
!
CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.


   ! Local Variables:

REAL(ReKi)                   :: Damp   (6,6)                                    ! Damping matrix.
REAL(ReKi)                   :: Stff   (6,6)                                    ! Stiffness/restoring matrix.

INTEGER(4)                   :: I                                               ! Generic index.
INTEGER(4)                   :: J                                               ! Generic index.



Damp (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Damp (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

Stff (1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
Stff (6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

TwrAM(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
TwrAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

TwrFt(1)   = 0.0
TwrFt(2)   = 0.0
TwrFt(3)   = 0.0
TwrFt(4)   = 0.0
TwrFt(5)   = 0.0
TwrFt(6)   = 0.0

DO J = 1,6
   DO I = 1,6
      TwrFt(I) = TwrFt(I) - Damp(I,J)*XD(J) - Stff(I,J)*X(J)
   ENDDO
ENDDO



RETURN
END SUBROUTINE UserTwrLd
!=======================================================================
!=======================================================================
!SUBROUTINE UserTwrLd ( JNode, X, XD, ZTime, DirRoot, TwrAM, TwrFt )
!
!
!   ! This routine implements the distributed springs foundation model
!   !   for the NREL 5MW Offshore Baseline Wind Turbine.  A direct CALL
!   !   to routine MorisonTwrLd() is made in order to obtain the
!   !   hydrodynamic loads on the monopile between the mudline and the
!   !   free surface.  The distributed springs at each tower node are
!   !   determined by interpolating into the array of discrete springs
!   !   specified by Patrik Passon of the University of Stuttgart in
!   !   Germany as documented in:
!   !      "OC3-Derivation and Description of the Soil-Pile-Interaction Models.pdf"
!   !      "OC3-Soil-Pile_Interaction_Model.xls"
!
!   ! The order of indices in all arrays passed to and from this routine is as
!   !   follows:
!   !      1 = Current tower element surge / xi-component of translation
!   !      2 = Current tower element sway  / yi-component of translation
!   !      3 = Current tower element heave / zi-component of translation
!   !      4 = Current tower element roll  / xi-component of rotation
!   !      5 = Current tower element pitch / yi-component of rotation
!   !      6 = Current tower element yaw   / zi-component of rotation
!
!   ! This routine was modified by Yang Yang at March 23, 2018 for the earthquake loads calculation using distributed springs method.
!   ! Yang: Start of modifying the UserTwrLd subroutine
!   !USE                             FixedBottomSupportStructure
!   !USE                             NWTC_Num
!   !USE                             Precision
!   !USE                             Waves, ONLY:WtrDpth, NWaveKin0, WaveKinzi0, DZNodes
!   !USE                             Constants
!   !USE                             MassInert
!   !USE                             SimCont,ONLY:DT
!   !USE                             SSI
!   !
!   IMPLICIT                        NONE
!   !
!   !
!   !! Passed Variables:
!   !
!   !REAL(ReKi), INTENT(OUT)      :: TwrAM  (6,6)                                    ! Added mass matrix per unit length of current tower element, kg/m, kg-m/m, kg-m^2/m.
!   !REAL(ReKi), INTENT(OUT)      :: TwrFt    (6)                                    ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force per unit length (in N/m) at the current tower element and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the tower moment per unit length (in N-m/m) acting at the current tower element associated with everything but the added-mass effects; positive forces are in the direction of motion.
!   !REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement (in m  ) of the current tower node and the 3 components of the rotational displacement       (in rad  ) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
!   !REAL(ReKi), INTENT(IN )      :: XD       (6)                                    ! The 3 components of the translational velocity     (in m/s) of the current tower node and the 3 components of the rotational (angular) velocity (in rad/s) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
!   !REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!   !
!   !INTEGER(4), INTENT(IN )      :: JNode                                           ! The number of the current tower node / element, (-). [1 to TwrNodes]
!   !
!   !CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!   !
!   !
!   !! Local Variables:
!   !
!   !REAL(ReKi)                   :: DZFract                                         ! The fraction of the current tower element that is below the seabed (0.0 <= DZFract <= 1.0): 0.0 = the element is entirely above the seabed, 1.0 = element is entirely below the seabed (-)
!   !REAL(ReKi), ALLOCATABLE      :: Spring   (:)                                    ! Discrete springs as specified by Patrik Passon (N/m)
!   !REAL(ReKi), ALLOCATABLE      :: Dashpot  (:)                                    ! Discrete Dashpot as specified by Patrik Passon (Ns/m)
!   !REAL(ReKi), ALLOCATABLE,SAVE :: Stff     (:)                                    ! The values of Spring(:) interpolated to each of the undeflected tower node elevations (as stored in the WaveKinzi0(:) array) (N/m)
!   !REAL(ReKi), ALLOCATABLE,SAVE :: Damp     (:)                                    ! The values of Dashpot(:) interpolated to each of the undeflected tower node elevations (as stored in the WaveKinzi0(:) array) (Ns/m)
!   !REAL(ReKi), ALLOCATABLE      :: zi       (:)                                    ! Elevations (w.r.t. MSL) were the discrete springs are specified by Patrik Passon (meters)
!   !
!   !INTEGER(4),ALLOCATABLE       :: SpringID (:)
!   !INTEGER(4)                   :: i,j,k                                           ! Generic index
!   !INTEGER(4)                   :: LastInd  = 1                                    ! Index into the arrays saved from the last call as a starting point for this call
!   !INTEGER(4)                   :: NSpring                                         ! Number of discrete springs specifed by Patrik Passon (-) 
!   !INTEGER(4)                   :: Sttus                                           ! Status returned by an attempted allocation.
!   !INTEGER(4)                   :: iTimeStep                                       ! The current time step
!   !
!   !
!   !LOGICAL(1), SAVE             :: FirstPas  = .TRUE.                              ! When .TRUE., indicates we're on the first pass
!   !LOGICAL(1)                   :: DampProp                                        ! When .TRUE. Contains the damping properties
!   !LOGICAL (1)                  :: ex
!   ! 
!   !!TwrAM(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!   !!TwrAM(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!   !!TwrAM(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!   !!TwrAM(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!   !!TwrAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!   !!TwrAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!   !!
!   !!TwrFt(1)   = 0.0
!   !!TwrFt(2)   = 0.0
!   !!TwrFt(3)   = 0.0
!   !!TwrFt(4)   = 0.0
!   !!TwrFt(5)   = 0.0
!   !!TwrFt(6)   = 0.0
!   !!iTimeStep = ZTime/DT+1
!   !
!   !   ! Perform initialization calculations on first pass only:
!   !
!   !IF ( FirstPas )  THEN
!   !         select case (SSIMode)
!   !           case (3)   ! DS model 
!   !              INQUIRE (file=trim(InFileDS), exist=ex)
!   !              if (ex .eq. .FALSE.) then
!   !                  CALL ProgAbort( 'The File'// trim(InFileDS)// ' which contains the properties of distributed springs cannot be found' )
!   !              else
!   !                 open (1002,file = trim(InFileDS))
!   !                  do i=1,3
!   !                     read(1002,*)  ! skip the comments
!   !                  enddo
!   !                  read (1002,*) NSpring
!   !                  allocate ( zi    (NSpring  ) , STAT=Sttus )
!   !                  if ( Sttus /= 0 )  then
!   !                     CALL Abort(' Error allocating memory for the zi array.')
!   !                  endif
!   !                  
!   !                  allocate ( Spring(NSpring  ) , STAT=Sttus )
!   !                  if ( Sttus /= 0 )  then
!   !                     CALL Abort(' Error allocating memory for the Spring array.')
!   !                  endif
!   !                  allocate ( Dashpot(NSpring  ) , STAT=Sttus )
!   !                  if ( Sttus /= 0 )  then
!   !                     CALL Abort(' Error allocating memory for the Dashpot array.')
!   !                  endif
!   !                  allocate ( SpringID(NSpring  ) , STAT=Sttus )
!   !                  if ( Sttus /= 0 )  then
!   !                     CALL Abort(' Error allocating memory for the SpringID array.')
!   !                  endif
!   !                  read(1002,*) DampProp
!   !                  read(1002,*)
!   !                  read(1002,*)
!   !                  ! read the spring and dashpot properties
!   !                  do i=1,NSpring
!   !                     if (DampProp) then
!   !                        read(1002,*) SpringID(i),zi(i),Spring(i),Dashpot(i)
!   !                     else
!   !                        read(1002,*) SpringID(i),zi(i),Spring(i)
!   !                        Dashpot(i) = 0
!   !                     endif
!   !                  enddo
!   !                 close (1002)
!   !              endif
!   !       end select
!   !
!   !   ! Interpolate to find the foundation stiffness, Stff(:), at each of the
!   !   !   undeflected tower node elevations (as stored in the WaveKinzi0(:)
!   !   !   array):
!   !   
!   !   ALLOCATE ( Stff  (NWaveKin0) , STAT=Sttus )
!   !   IF ( Sttus /= 0 )  THEN
!   !      CALL Abort(' Error allocating memory for the Stff array.')
!   !   ENDIF
!   !   ALLOCATE ( Damp  (NWaveKin0) , STAT=Sttus )
!   !   IF ( Sttus /= 0 )  THEN
!   !      CALL Abort(' Error allocating memory for the Damp array.')
!   !   ENDIF
!   !
!   !   DO j = 1,NWaveKin0   ! Loop through the tower nodes / elements
!   !     ! write (*,*) 'NWaveKin0= ',NWaveKin0
!   !      IF ( WaveKinzi0(j) > -WtrDpth )  THEN  ! .TRUE. if the current node lies above the seabed.
!   !         Stff(j) = 0.0
!   !         Damp(j) = 0
!   !      ELSE                                   ! The current tower node must lie  below the seabed.
!   !         Stff(j) = InterpStp( WaveKinzi0(j), zi(:), Spring(:), LastInd, NSpring )
!   !         Damp(j) = InterpStp( WaveKinzi0(j), zi(:), Dashpot(:), LastInd, NSpring )
!   !      ENDIF
!   !   ENDDO                ! J - Tower nodes / elements
!   !   FirstPas   = .FALSE. ! Don't enter here again!
!   !ENDIF
!   !
!   !
!   !
!   !! Yang: We have called MorisonTwrLd in TwrLoading subroutine for the calculation of hydrodynamic loads of local tower element.
!   !!       CALL MorisonTwrLd ( JNode, TwrDiam, TwrCA, TwrCD, X, XD, ZTime, TwrAM, TwrFt )
!   !
!   !! Find the fraction of the current tower element that is below the seabed:
!   !
!   !if (     ( WaveKinzi0(JNode) - 0.5*DZNodes(JNode) ) >= -WtrDpth )  then ! .TRUE. if the current tower element lies entirely above the seabed.
!   !   DZFract = 0.0
!   !elseif ( ( WaveKinzi0(JNode) + 0.5*DZNodes(JNode) ) <= -WtrDpth )  then ! .TRUE. if the current tower element lies entirely below the seabed.
!   !   DZFract = 1.0
!   !else                                                                    ! The seabed must fall somewhere along the current tower element; thus, interpolate.
!   !   DZFract = ( -WtrDpth - ( WaveKinzi0(JNode) - 0.5*DZNodes(JNode) ) )/DZNodes(JNode)
!   !endif
!   !! Compute the foundation loads using uncoupled linear springs for the
!   !!   portion of the current tower element that lies below the seabed: 
!   !
!   !if ( DZFract > 0.0 )  then ! .TRUE. if a portion of the current tower element lies below the seabed.
!   !   do k = 1,2     ! Loop through the xi- (1) and yi- (2) directions
!   !
!   !      TwrFt(k) = TwrFt(k) - Stff(JNode)*DZFract*X(k)-Damp(JNode)*DZFract*XD(k)
!   !      ! Debug
!   !      !if (ZTime > 110) write(*,*) 'Ft=', TwrFt(k), 'PtfmDisp= ', PtfmDisp(k,iTimeStep), 'Stff=',Stff(JNode), 'JNode=',JNode 
!   !      
!   !      
!   !   enddo          ! k - The xi- (1) and yi- (2) directions(Currently, the vertical springs are ignored)
!   !
!   !endif
!   
!   
!   
!RETURN
!! Yang: End of modifying the UserTwrLd subroutine
!END SUBROUTINE UserTwrLd
!=======================================================================


SUBROUTINE UserRFrl ( RFrlDef, RFrlRate, ZTime, DirRoot, RFrlMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! rotor-furl spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the rotor-furl DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_RFrl) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the rotor-furl hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: RFrlDef                                         ! Rotor-furl angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: RFrlMom                                         ! Rotor-furl restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: RFrlRate                                        ! Rotor-furl angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



RFrlMom = 0.0



RETURN
END SUBROUTINE UserRFrl
!=======================================================================
SUBROUTINE UserTeet ( TeetDef, TeetRate, ZTime, DirRoot, TeetMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! teeter spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the teeter DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_Teet) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the teeter hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: TeetDef                                         ! Rotor-teeter angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: TeetMom                                         ! Rotor-teeter restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: TeetRate                                        ! Rotor-teeter angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TeetMom = 0.0



RETURN
END SUBROUTINE UserTeet
!=======================================================================
SUBROUTINE UserTFin ( TFrlDef , TFrlRate, ZTime   , DirRoot, &
                      TFinCPxi, TFinCPyi, TFinCPzi,          &
                      TFinCPVx, TFinCPVy, TFinCPVz,          &
                      TFinAOA , TFinQ   ,                    &
                      TFinCL  , TFinCD  ,                    &
                      TFinKFx , TFinKFy                        )


   ! This is a dummy routine for holding the place of user-specified
   ! computations for tail fin aerodynamic loads.  Modify this code to
   ! create your own logic.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(OUT)      :: TFinAOA                                         ! Angle-of-attack between the relative wind velocity and tail fin chordline, rad.
REAL(ReKi), INTENT(OUT)      :: TFinCD                                          ! Tail fin drag            coefficient resulting from current TFinAOA, (-).
REAL(ReKi), INTENT(OUT)      :: TFinCL                                          ! Tail fin lift            coefficient resulting from current TFinAOA, (-).
REAL(ReKi), INTENT(IN )      :: TFinCPVx                                        ! Absolute Velocity of the tail center-of-pressure along tail fin chordline pointing toward tail fin trailing edge, m/s.
REAL(ReKi), INTENT(IN )      :: TFinCPVy                                        ! Absolute Velocity of the tail center-of-pressure normal to plane of tail fin pointing towards suction surface   , m/s.
REAL(ReKi), INTENT(IN )      :: TFinCPVz                                        ! Absolute Velocity of the tail center-of-pressure in plane of tail fin normal to chordline and nominally upward  , m/s.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Improve the description of input arguments TFinCPxi, TFinCPyi, and
!jmj   TFinCPzi:
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPxi                                        ! Downwind distance from the inertial frame origin to the tail fin center-of-pressure, m.
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPyi                                        ! Lateral  distance from the inertial frame origin to the tail fin center-of-pressure, m.
!remove6.02aREAL(ReKi), INTENT(IN )      :: TFinCPzi                                        ! Vertical distance from the inertial frame origin to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPxi                                        ! Downwind distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPyi                                        ! Lateral  distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
REAL(ReKi), INTENT(IN )      :: TFinCPzi                                        ! Vertical distance from the inertial frame origin at ground level [onshore] or MSL [offshore] to the tail fin center-of-pressure, m.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi), INTENT(OUT)      :: TFinKFx                                         ! Aerodynamic force  at the tail fin center-of-pressure (point K) along tail fin chordline pointing toward tail fin trailing edge, N.
REAL(ReKi), INTENT(OUT)      :: TFinKFy                                         ! Aerodynamic force  at the tail fin center-of-pressure (point K) normal to plane of tail fin pointing towards suction surface   , N.
REAL(ReKi), INTENT(OUT)      :: TFinQ                                           ! Dynamic pressure of the relative wind velocity, Pa.
REAL(ReKi), INTENT(IN )      :: TFrlDef                                         ! Tail-furl angular deflection, rad.
REAL(ReKi), INTENT(IN )      :: TFrlRate                                        ! Tail-furl angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TFinAOA = 0.0
TFinCL  = 0.0
TFinCD  = 0.0
TFinQ   = 0.0
TFinKFx = 0.0
TFinKFy = 0.0



RETURN
END SUBROUTINE UserTFin
!=======================================================================
SUBROUTINE UserTFrl ( TFrlDef, TFrlRate, ZTime, DirRoot, TFrlMom )


   ! This is a dummy routine for holding the place of a user-specified
   ! tail-furl spring/damper.  Modify this code to create your own device.

   ! NOTE: If you (the user) wants to switch on-or-off the tail-furl DOF at
   !       runtime from this user-defined routine, then do the following:
   !          (1) USE MODULE DOFs().
   !          (2) Type in "DOF_Flag(DOF_TFrl) = VALUE" where VALUE = .TRUE. or
   !              .FALSE. depending on whether you want to turn-on or turn-off
   !              the DOF, respectively.  Turning off the DOF forces the
   !              current RATE to remain fixed.  If the rate is currently zero,
   !              the current POSITION will remain fixed as well.
   !       This technique is useful, for example, if the tail-furl hinge has
   !       an electromagnetic latch that will unlock and relock the hinge under
   !       certain specified conditions.
   !       Note that this technique WILL NOT work for user-defined routines
   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
   !       preprocessor.


USE                             Precision


IMPLICIT                        NONE


   ! Passed Variables:

REAL(ReKi), INTENT(IN )      :: TFrlDef                                         ! Tail-furl angular deflection, rad.
REAL(ReKi), INTENT(OUT)      :: TFrlMom                                         ! Tail-furl restoring moment, N-m.
REAL(ReKi), INTENT(IN )      :: TFrlRate                                        ! Tail-furl angular rate, rad/s
REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.

CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



TFrlMom = 0.0



RETURN
END SUBROUTINE UserTFrl
!=======================================================================
!SUBROUTINE UserVSCont ( HSS_Spd, GBRatio, NumBl, ZTime, DT, GenEff, DelGenTrq, DirRoot, GenTrq, ElecPwr )
!
!
!   ! This is a dummy routine for holding the place of a user-specified
!   ! variable-speed torque and power control model.  Modify this code to
!   ! create your own model.
!
!   ! NOTE: If you (the user) wants to switch on-or-off the generator DOF at
!   !       runtime from this user-defined routine, then do the following:
!   !          (1) USE MODULE DOFs().
!   !          (2) Type in "DOF_Flag(DOF_GeAz) = VALUE" where VALUE = .TRUE. or
!   !              .FALSE. depending on whether you want to turn-on or turn-off
!   !              the DOF, respectively.  Turning off the DOF forces the
!   !              current RATE to remain fixed.  If the rate is currently zero,
!   !              the current POSITION will remain fixed as well.
!   !       Note that this technique WILL NOT work for user-defined routines
!   !       written for ADAMS datasets extracted using the FAST-to-ADAMS
!   !       preprocessor.
!
!
!USE                             Precision
!
!
!IMPLICIT                        NONE
!
!
!   ! Passed Variables:
!
!INTEGER(4), INTENT(IN )      :: NumBl                                           ! Number of blades, (-).
!
!REAL(ReKi), INTENT(IN )      :: DelGenTrq                                       ! Pertubation in generator torque used during FAST linearization (zero otherwise), N-m.
!REAL(ReKi), INTENT(IN )      :: DT                                              ! Integration time step, sec.
!REAL(ReKi), INTENT(OUT)      :: ElecPwr                                         ! Electrical power (account for losses), watts.
!REAL(ReKi), INTENT(IN )      :: GBRatio                                         ! Gearbox ratio, (-).
!REAL(ReKi), INTENT(IN )      :: GenEff                                          ! Generator efficiency, (-).
!REAL(ReKi), INTENT(OUT)      :: GenTrq                                          ! Electrical generator torque, N-m.
!REAL(ReKi), INTENT(IN )      :: HSS_Spd                                         ! HSS speed, rad/s.
!REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!
!GenTrq  = 0.0 + DelGenTrq  ! Make sure to add the pertubation on generator torque, DelGenTrq.  This is used only for FAST linearization (it is zero otherwise).
!
!
!   ! The generator efficiency is either additive for motoring,
!   !   or subtractive for generating power.
!
!IF ( GenTrq > 0.0 )  THEN
!   ElecPwr = GenTrq*HSS_Spd*GenEff
!ELSE
!   ElecPwr = GenTrq*HSS_Spd/GenEff
!ENDIF
!
!
!
!RETURN
!END SUBROUTINE UserVSCont
!=======================================================================

