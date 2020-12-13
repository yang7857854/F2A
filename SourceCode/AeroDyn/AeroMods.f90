! AeroDyn Modules
!=======================================================================
MODULE AD_IOParams

!bjj: why aren't these parameters?

   ! Contains input/output parameters.

INTEGER(4)                   :: UnADin  = 90 ! ipt file
INTEGER(4)                   :: UnADopt = 92 ! opt file
INTEGER(4)                   :: UnAirfl = 93 ! Airfoil data file
INTEGER(4)                   :: UnWind  = 91 ! HH or FF wind file

LOGICAL                      :: WrOptFile  = .TRUE.   ! Write the .opt file?


END MODULE AD_IOParams
!=======================================================================
MODULE AeroTime


   ! Contains aero calc information.


USE                             Precision


REAL(ReKi),SAVE                   :: OLDTIME   ! The previous time AeroDyn's loads were calculated
REAL(DbKi),SAVE                   :: TIME      ! Current time simulation time

REAL(ReKi),SAVE                   :: DT        ! actual difference between Time and OldTime when loads are calculated
REAL(ReKi),SAVE                   :: DTAERO    ! desired time interval for aerodynamics calculations


END MODULE AeroTime
!=======================================================================
MODULE Airfoil


   ! Contains airfoil information.


USE                             Precision


REAL(ReKi), ALLOCATABLE,SAVE      :: AL    ( :, : )        ! Table of angles of attack
REAL(ReKi), ALLOCATABLE,SAVE      :: CD    ( :, :, : )     ! Table of drag coefficients
REAL(ReKi), ALLOCATABLE,SAVE      :: CL    ( :, :, : )     ! Table of lift coefficients
REAL(ReKi), ALLOCATABLE,SAVE      :: CM    ( :, :, : )     ! Table of pitching moment coefficients
REAL(ReKi)                   :: MulTabLoc = 0.0
REAL(ReKi), ALLOCATABLE,SAVE      :: MulTabMet ( :, :)
REAL(ReKi)             ,SAVE      :: PMC

INTEGER   , PARAMETER        :: MAXTABLE = 10 !bjj: pjm increased this to 20
INTEGER   , ALLOCATABLE,SAVE      :: NFOIL ( : )           ! indices of the airfoil data file used for each element
INTEGER   , ALLOCATABLE ,SAVE     :: NLIFT ( : )           ! Number of aerodata points in each airfoil file
INTEGER   , ALLOCATABLE,SAVE      :: NTables  ( : )        ! number of airfoil data tables
INTEGER                ,SAVE      :: NumCL                 ! maximum number of aerodata points in all airfoil files {=max(NFoil(:)}
INTEGER                ,SAVE      :: NumFoil               ! number of different airfoil files used

CHARACTER(1024), ALLOCATABLE,SAVE   :: FOILNM ( : )          ! names of the data files that contain airfoil data

END MODULE Airfoil
!=======================================================================
MODULE Bedoes


   ! Contains Beddoes dynamic stall info.

!bjj: some "constants" could probably be parameters instead of set in BedDat()

USE                             Precision


REAL(ReKi), ALLOCATABLE,SAVE      :: ADOT  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: ADOT1 ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: AFE   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: AFE1  ( :, : )
REAL(ReKi)             ,SAVE      :: AN
REAL(ReKi), ALLOCATABLE,SAVE      :: ANE   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: ANE1  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: AOD   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: AOL   ( :, : )
REAL(ReKi)             ,SAVE      :: AS              ! Speed of sound for Mach number calculation
REAL(ReKi)             ,SAVE      :: CC
REAL(ReKi), ALLOCATABLE,SAVE      :: CDO   ( :, : )
REAL(ReKi)             ,SAVE      :: CMI
REAL(ReKi)             ,SAVE      :: CMQ
REAL(ReKi)             ,SAVE      :: CN
REAL(ReKi), ALLOCATABLE,SAVE      :: CNA   ( :, : )
REAL(ReKi)             ,SAVE      :: CNCP
REAL(ReKi)             ,SAVE      :: CNIQ
REAL(ReKi), ALLOCATABLE,SAVE      :: CNP   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNP1  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNPD  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNPD1 ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNPOT ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNPOT1( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNS   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNSL  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNV   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CVN   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CVN1  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DF    ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DFAFE ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DFAFE1( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DFC   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DN    ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DPP   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DQ    ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DQP   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DQP1  ( :, : )
REAL(ReKi)             ,SAVE      :: DS
REAL(ReKi)             ,SAVE      :: FK
REAL(ReKi)             ,SAVE      :: FP
REAL(ReKi)             ,SAVE      :: FPC
REAL(ReKi), ALLOCATABLE,SAVE      :: FSP   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: FSP1  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: FSPC  ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: FSPC1 ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: FTB   ( :, :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: FTBC  ( :, :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDCNV( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDDF ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDDFC( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDDN ( :, : )
REAL(ReKi), ALLOCATABLE ,SAVE     :: OLDDPP( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDDQ ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDTAU( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDXN ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: OLDYN ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: QX    ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: QX1   ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: TAU   ( :, : )
REAL(ReKi)             ,SAVE      :: TF              ! Time constant applied to location of the separation point
REAL(ReKi)             ,SAVE      :: TP              ! Time constant for pressure lag
REAL(ReKi)             ,SAVE      :: TV              ! Time constant for strength of shed vortex
REAL(ReKi)             ,SAVE      :: TVL             ! Non-dimensional time of transit for the vortex moving across the airfoil surface
REAL(ReKi), ALLOCATABLE,SAVE      :: XN    ( :, : )
REAL(ReKi), ALLOCATABLE,SAVE      :: YN    ( :, : )

LOGICAL,    ALLOCATABLE,SAVE      :: BEDSEP ( :, : )
LOGICAL,    ALLOCATABLE,SAVE      :: OLDSEP ( :, : )
LOGICAL                ,SAVE      :: SHIFT
LOGICAL                ,SAVE      :: VOR


END MODULE Bedoes
!=======================================================================
MODULE Blade


   ! Contains blade information.


USE                             Precision


REAL(ReKi), ALLOCATABLE,SAVE      :: C       (:)     ! Chord of each blade element (FROM INPUT FILE)
REAL(ReKi), ALLOCATABLE,SAVE      :: DR      (:)     ! Span-wise width of the element (length of the element, centered at RELM(i)) (FROM INPUT FILE)
REAL(ReKi)             ,SAVE      :: R               ! rotor radius

INTEGER                ,SAVE      :: NB              ! number of blades


END MODULE Blade
!=======================================================================
MODULE DynInflow


   ! Contains dynamic inflow information.


USE                             Precision


INTEGER   , PARAMETER        :: MAXINFL  = 6
INTEGER   , PARAMETER        :: MAXINFL0 = 2
INTEGER             ,SAVE         :: MminR    ( maxInfl, maxInfl )
INTEGER             ,SAVE         :: MminusR  ( maxInfl, maxInfl )
INTEGER             ,SAVE         :: MplusR   ( maxInfl, maxInfl )
INTEGER             ,SAVE         :: MRvector ( maxInfl )
INTEGER             ,SAVE         :: NJvector ( maxInfl )

REAL(ReKi)          ,SAVE         :: dAlph_dt ( maxInfl, 4 )
REAL(ReKi)          ,SAVE         :: dBeta_dt ( maxInfl0+1 : maxInfl, 4 )
REAL(ReKi)          ,SAVE         :: DT0
REAL(ReKi)          ,SAVE         :: GAMMA    ( maxInfl, maxInfl )
REAL(ReKi)          ,SAVE         :: old_Alph (              maxInfl )
REAL(ReKi)          ,SAVE         :: old_Beta ( maxInfl0+1 : maxInfl )
REAL(ReKi)          ,SAVE         :: old_LmdM
REAL(ReKi)          ,SAVE         :: oldKai
REAL(ReKi)          ,SAVE         :: PhiLqC   (              maxInfl )
REAL(ReKi)          ,SAVE         :: PhiLqS   ( maxInfl0+1 : maxInfl )
REAL(ReKi)          ,SAVE         :: Pzero
REAL(ReKi), ALLOCATABLE      :: RMC_SAVE( : , :, : )  !Store element parameters for GDW
REAL(ReKi), ALLOCATABLE      :: RMS_SAVE( : , :, : )
REAL(ReKi)          ,SAVE         :: TipSpeed
REAL(ReKi)          ,SAVE         :: totalInf
REAL(ReKi)          ,SAVE         :: Vparam
REAL(ReKi)          ,SAVE         :: Vtotal
REAL(ReKi)          ,SAVE         :: xAlpha   (              maxInfl )
REAL(ReKi)          ,SAVE         :: xBeta    ( maxInfl0+1 : maxInfl )
REAL(ReKi)          ,SAVE         :: xKai
REAL(ReKi)          ,SAVE         :: XLAMBDA_M
REAL(ReKi)          ,SAVE         :: xLcos    ( maxInfl, maxInfl )
REAL(ReKi)          ,SAVE         :: xLsin    ( maxInfl0+1 : maxInfl , maxInfl0+1 : maxInfl )
REAL(ReKi)          ,SAVE         :: xMinv    ( maxInfl )


END MODULE DynInflow
!=======================================================================
MODULE Element


   ! Contains element specific information.


USE                             Precision


REAL(ReKi), ALLOCATABLE,SAVE      :: A       (:,:)       ! induction factor?
REAL(ReKi), ALLOCATABLE,SAVE      :: AP      (:,:)
REAL(ReKi), ALLOCATABLE,SAVE      :: HLCNST  (:)         ! Hub-loss constant at each element
REAL(ReKi)          ,SAVE         :: PITNOW
REAL(ReKi), ALLOCATABLE,SAVE      :: RELM    (:)         ! Location of the center of the element; measured from the blade root. (INPUT FILE) Supposedly ignored by ADAMS.
REAL(ReKi), ALLOCATABLE,SAVE      :: TLCNST  (:)         ! Tip-loss constant at each element
REAL(ReKi), ALLOCATABLE,SAVE      :: TWIST   (:)         ! Twist of each blade element  (INPUT FILE)

INTEGER             ,SAVE         :: NELM                ! Number of elements per blade (INPUT FILE)

END MODULE Element
!=======================================================================
MODULE ElemInflow


   ! Contains element specific information associated with the inflow.


USE                             Precision

REAL(ReKi), ALLOCATABLE,SAVE      :: ALPHA(:,:)                                      ! Angle of attack                       of the inflow for the current blade, element, and time step.
REAL(ReKi), ALLOCATABLE,SAVE      :: W2(:,:)                                         ! The square of the relative wind speed of the inflow for the current blade, element, and time step.

END MODULE ElemInflow
!=======================================================================
MODULE ElOutParams


   ! Contains element output information.


USE                             Precision


REAL(ReKi), ALLOCATABLE,SAVE      :: AAA     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: AAP     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: ALF     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CDD     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CLL     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CMM     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: CNN     ( : )
REAL(ReKi), ALLOCATABLE ,SAVE     :: CTT     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DFNSAV  ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DFTSAV  ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: DynPres ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: PITSAV  ( : )
REAL(ReKi), ALLOCATABLE ,SAVE     :: PMM     ( : )
REAL(ReKi), ALLOCATABLE,SAVE      :: ReyNum  ( : )
REAL(ReKi)          ,SAVE         :: VXSAV
REAL(ReKi)          ,SAVE         :: VYSAV
REAL(ReKi)          ,SAVE         :: VZSAV

REAL(ReKi), ALLOCATABLE ,SAVE     :: SaveVX  ( :,: )         ! The velocity in the x direction at requested element, on each blade
REAL(ReKi), ALLOCATABLE,SAVE      :: SaveVY  ( :,: )         ! The velocity in the y direction at requested element, on each blade
REAL(ReKi), ALLOCATABLE ,SAVE     :: SaveVZ  ( :,: )         ! The velocity in the z direction at requested element, on each blade

INTEGER             ,SAVE         :: UnWndOut = 96           ! The unit for wind output at each element on each blade
INTEGER             ,SAVE         :: NumWndElOut             ! Number of wind elements to print
INTEGER   , ALLOCATABLE  ,SAVE    :: WndElPrList (:)
INTEGER   , ALLOCATABLE,SAVE      :: WndElPrNum  (:)
INTEGER   , ALLOCATABLE  ,SAVE    :: ElPrList (:)
INTEGER   , ALLOCATABLE  ,SAVE    :: ElPrNum  (:)
INTEGER             ,SAVE         :: NumElOut
INTEGER             ,SAVE         :: UnElem = 94


END MODULE ElOutParams
!=======================================================================
MODULE ErrCount


   ! Contains error counters.


INTEGER             ,SAVE         :: NumErr
INTEGER             ,SAVE         :: NumWarn


END MODULE ErrCount
!=======================================================================
MODULE InducedVel


   ! Contains induced velocity information.


USE                             Precision


REAL(ReKi)          ,SAVE         :: ATOLER                                   ! Convergence tolerance for induction factor
REAL(ReKi)          ,SAVE         :: EqAIDmult                                ! Multiplier for the drag term in the axial-induction equation.
REAL(ReKi)                        :: SumInfl = 0.0                            ! Initialize this value here for the first pass


END MODULE InducedVel
!=======================================================================
MODULE Rotor


   ! Contains rotor configuration information.


USE                             Precision


REAL(ReKi)          ,SAVE         :: AVGINFL         ! average induduced velocity at the previous time
REAL(ReKi)          ,SAVE         :: CTILT
REAL(ReKi)          ,SAVE         :: CYaw
REAL(ReKi)          ,SAVE         :: HH
REAL(ReKi)          ,SAVE         :: REVS
REAL(ReKi)          ,SAVE         :: STILT
REAL(ReKi)          ,SAVE         :: SYaw
REAL(ReKi)          ,SAVE         :: TILT
REAL(ReKi)          ,SAVE         :: YawAng
REAL(ReKi)          ,SAVE         :: YAWVEL


END MODULE Rotor
!=======================================================================
MODULE Switch


   ! Defines variables to control program options.


LOGICAL             ,SAVE         :: DSTALL       ! Dynamic stall model: TRUE = BEDDOES; FALSE = STEADY
LOGICAL             ,SAVE         :: DYNINFL      ! Dynamic inflow: TRUE = DYNIN; FALSE = EQUIL
LOGICAL             ,SAVE         :: DYNINIT
LOGICAL             ,SAVE         :: ELEMPRN
LOGICAL             ,SAVE         :: EquilDA
LOGICAL             ,SAVE         :: EquilDT
LOGICAL             ,SAVE         :: GTECH
LOGICAL             ,SAVE         :: HLOSS        ! Hub loss: TRUE = PRAND; FALSE = NONE
LOGICAL             ,SAVE         :: MultiTab
LOGICAL             ,SAVE         :: PMOMENT      ! Pitching moment: TRUE = USE_CM; FALSE = NO_CM
LOGICAL             ,SAVE         :: Reynolds
LOGICAL             ,SAVE         :: SIUNIT       ! TRUE = scientific units; FALSE = english
LOGICAL             ,SAVE         :: SKEW
LOGICAL             ,SAVE         :: SWIRL
LOGICAL             ,SAVE         :: TLOSS
LOGICAL             ,SAVE         :: WAKE


END MODULE Switch
!=======================================================================
MODULE TwrProps


   ! Contains tower aero information.


USE                             Precision


REAL(ReKi), ALLOCATABLE ,SAVE     :: TwrHtFr ( : )
REAL(ReKi), ALLOCATABLE ,SAVE     :: TwrWid  ( : )
REAL(ReKi), ALLOCATABLE ,SAVE     :: TwrCD   ( :, : )
REAL(ReKi), ALLOCATABLE ,SAVE     :: TwrRe   ( : )

REAL(ReKi)          ,SAVE         :: VTwr(3)
REAL(ReKi)          ,SAVE         :: Tower_Wake_Constant   ! Constant for tower wake model = 0 full potential flow = 0.1 model of Bak et al.


INTEGER,    ALLOCATABLE  ,SAVE    :: NTwrCDCol (:)         ! The tower CD column that represents a particular tower height
INTEGER             ,SAVE         :: NTwrHt                ! The number of tower height rows in the table
INTEGER             ,SAVE         :: NTwrRe                ! The number of tower Re entry rows in the table
INTEGER             ,SAVE         :: NTwrCD                ! The number of tower CD columns in the table

LOGICAL             ,SAVE         :: TwrPotent             ! Tower potential flow calculation
LOGICAL             ,SAVE         :: TwrShadow             ! Tower Shadow calculation


REAL(ReKi)          ,SAVE         :: SHADHWID                                !
REAL(ReKi)          ,SAVE         :: TSHADC1                                 !
REAL(ReKi)          ,SAVE         :: TSHADC2                                 !
REAL(ReKi)          ,SAVE         :: TWRSHAD                                 !

REAL(ReKi)          ,SAVE         :: T_Shad_Refpt !This was a local variable -- with new tower influence, it should be removed


LOGICAL             ,SAVE         :: PJM_Version = .FALSE.

CHARACTER(1024)     ,SAVE         :: TwrFile               ! Name of the tower properties input file

END MODULE TwrProps
!=======================================================================

MODULE Wind


   ! Module Wind is used for wind variables.


USE                             Precision


REAL(ReKi)          ,SAVE         :: ANGFLW                                  !
REAL(ReKi)          ,SAVE         :: CDEL                                    !
REAL(ReKi)          ,SAVE         :: KinVisc                                 ! KINEMATIC VISCOSITY   Units^2/SEC
REAL(ReKi)          ,SAVE         :: RHO                                     ! Ambient Air Density
REAL(ReKi)          ,SAVE         :: SDEL                                    !
REAL(ReKi)          ,SAVE         :: VROTORX                                 !
REAL(ReKi)          ,SAVE         :: VROTORY                                 !
REAL(ReKi)          ,SAVE         :: VROTORZ                                 !

END MODULE Wind
!=======================================================================

