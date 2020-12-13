MODULE UserWind
!  The purpose of this module is to allow user-defined wind.  
!----------------------------------------------------------------------------------------------------

   USE                           NWTC_Library
   USE                           SharedInflowDefns

   IMPLICIT                      NONE
   PRIVATE
    
    
      ! define variables for UserWind here
      
   LOGICAL, SAVE              :: Initialized = .FALSE.         ! This variable indicates if the initialization routine has been run
   
   REAL(ReKi)                 :: UWmeanU                       ! Possibly instantaneous, disk-averaged wind speeds.
   REAL(ReKi)                 :: UWmeanV                       !
   REAL(ReKi)                 :: UWmeanW                       !   
   

      ! allow the initialization and termination routines to be public (called from outside)

   PUBLIC                     :: UsrWnd_Init
   PUBLIC                     :: UsrWnd_Terminate
   PUBLIC                     :: UsrWnd_GetValue
   PUBLIC                     :: UsrWnd_GetWindSpeed

CONTAINS
!====================================================================================================
SUBROUTINE UsrWnd_Init(ErrStat)
!  This subroutine is called at the beginning of
!----------------------------------------------------------------------------------------------------

   INTEGER,    INTENT(OUT)    :: ErrStat           ! return 0 if no errors; non-zero otherwise

   !-------------------------------------------------------------------------------------------------
   ! Check that the module hasn't already been initialized.
   !-------------------------------------------------------------------------------------------------
      
   IF ( Initialized ) THEN  
      CALL WrScr( ' UserWind has already been initialized.' )
      ErrStat = 1
      RETURN
   ELSE
      ErrStat = 0
      CALL NWTC_Init()
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Perform any initialization steps here (read input files, etc.)
   !-------------------------------------------------------------------------------------------------
   
   CALL WrScr( '***** NOTE: User-defined wind employed *****' )


      ! Set the disk-average wind vector.
   
   UWmeanU = 10.0
   UWmeanV =  0.0
   UWmeanW =  0.0

   
   !-------------------------------------------------------------------------------------------------
   ! Set the initialization flag
   !-------------------------------------------------------------------------------------------------
   
   Initialized = .TRUE.
   
   RETURN

END SUBROUTINE UsrWnd_Init
!====================================================================================================
FUNCTION UsrWnd_GetValue(VarName, ErrStat)
!  This function returns a real scalar value whose name is listed in the VarName input argument.
!  If the name is not recognized, an error is returned in ErrStat.
!----------------------------------------------------------------------------------------------------
   
   CHARACTER(*),   INTENT(IN)    :: VarName
   INTEGER,        INTENT(OUT)   :: ErrStat           ! return 0 if no errors; non-zero otherwise
   REAL(ReKi)                    :: UsrWnd_GetValue

   
   CHARACTER(20)                 :: VarNameUC         ! upper-case VarName
   
   
   !-------------------------------------------------------------------------------------------------
   ! Check that the module has been initialized.
   !-------------------------------------------------------------------------------------------------   

   IF ( .NOT. Initialized ) THEN   
      CALL WrScr( 'Initialize UserWind before calling its subroutines.' )
      ErrStat = 1
      RETURN
   ELSE
      ErrStat = 0
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Return the requested values.
   !-------------------------------------------------------------------------------------------------   

   VarNameUC = VarName
   CALL Conv2UC( VarNameUC )

   SELECT CASE ( TRIM(VarNameUC) )
   
      CASE ('MEANU' )
         UsrWnd_GetValue = UWmeanU
         
      CASE ('MEANV' )
         UsrWnd_GetValue = UWmeanV

      CASE ('MEANW' )
         UsrWnd_GetValue = UWmeanW
         
      CASE DEFAULT
         CALL WrScr( ' Invalid variable name in UsrWnd_GetValue().' )
         ErrStat = 1
         
   END SELECT
      
   

END FUNCTION UsrWnd_GetValue
!====================================================================================================
FUNCTION UsrWnd_GetWindSpeed(Time, InputPosition, ErrStat)
! This function receives time and position (in InputInfo) where (undisturbed) velocities are 
! requested. It returns the velocities at the specified time and space.
!----------------------------------------------------------------------------------------------------
   
   REAL(ReKi),        INTENT(IN) :: Time
   REAL(ReKi),        INTENT(IN) :: InputPosition(3)        ! X,Y,Z (z is 0 at ground level)
   INTEGER,           INTENT(OUT):: ErrStat                 ! return 0 if no errors; non-zero otherwise
   TYPE(InflIntrpOut)            :: UsrWnd_GetWindSpeed
   
   
   !-------------------------------------------------------------------------------------------------
   ! Check that the module has been initialized.
   !-------------------------------------------------------------------------------------------------   
   IF ( .NOT. Initialized ) THEN   
      CALL WrScr( 'Initialize UserWind before calling its subroutines.' )
      ErrStat = 1
      RETURN
   ELSE
      ErrStat = 0
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Calculate the wind speed at this time and position.
   !-------------------------------------------------------------------------------------------------   
   !     Time
   !     X = InputPosition(1)           ! relative to the undeflected tower centerline (positive downwind)
   !     Y = InputPosition(2)           ! relative to the undeflected tower centerline (positive left when looking downwind)
   !     Z = InputPosition(3)           ! relative to the ground (0 is ground level)
   !-------------------------------------------------------------------------------------------------

      ! We'll test this with steady winds for now.

   UsrWnd_GetWindSpeed%Velocity(1) = 10.0    ! U velocity (along positive X)
   UsrWnd_GetWindSpeed%Velocity(2) =  0.0    ! V velocity (along positive Y)
   UsrWnd_GetWindSpeed%Velocity(3) =  0.0    ! V velocity (along positive Z)
   

END FUNCTION UsrWnd_GetWindSpeed
!====================================================================================================
SUBROUTINE UsrWnd_Terminate(ErrStat)
!  This subroutine is called at the end of program execution (including after fatal errors occur).  
!  It should close any files that could be open and deallocate any arrays that have been allocated.
!----------------------------------------------------------------------------------------------------
      
   INTEGER,    INTENT(OUT)    :: ErrStat           ! return 0 if no errors; non-zero otherwise

   ErrStat = 0

   !-------------------------------------------------------------------------------------------------
   ! Close files
   !-------------------------------------------------------------------------------------------------
      
      
   !-------------------------------------------------------------------------------------------------
   ! Deallocate arrays
   !-------------------------------------------------------------------------------------------------

      
   !-------------------------------------------------------------------------------------------------
   ! Set the initialization flag
   !-------------------------------------------------------------------------------------------------
   
   Initialized = .FALSE.

END SUBROUTINE UsrWnd_Terminate
!====================================================================================================
END MODULE UserWind
