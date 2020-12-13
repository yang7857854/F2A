!=======================================================================
MODULE ADAMSInput


   ! This MODULE stores FAST-to-ADAMS, ADAMS-specifc input parameters.


USE                             Precision


REAL(ReKi)                  ,SAVE :: BoomRad                                         ! Radius of the tail boom used for tail boom GRAPHICS.
REAL(ReKi)                  ,SAVE :: BPActrDmp                                       ! Blade pitch actuator damping          constant, (N-m/rad/s).
REAL(ReKi)                  ,SAVE :: BPActrSpr                                       ! Blade pitch actuator spring stiffness constant, (N-m/rad).
REAL(ReKi)                  ,SAVE :: CRatioBEA                                       ! The ratio of CMatrix to KMatrix for the blade extensional deflection.
REAL(ReKi)                  ,SAVE :: CRatioBGJ                                       ! The ratio of CMatrix to KMatrix for the blade torsion     deflection.
REAL(ReKi)                  ,SAVE :: CRatioTEA                                       ! The ratio of CMatrix to KMatrix for the tower extensional deflection.
REAL(ReKi)                  ,SAVE :: CRatioTGJ                                       ! The ratio of CMatrix to KMatrix for the tower torsion     deflection.
REAL(ReKi), PARAMETER             :: FrSrfcSpc =  5.0                                ! Distance between points on the still water level plane along the incident wave propogation heading direction for depicting the free surface where the elevation of the incident waves will be computed used for free surface GRAPHICS. (meters)  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
REAL(ReKi)                  ,SAVE :: GBoxLength                                      ! Length, width, height of the gearbox for gearbox GRAPHICS.
REAL(ReKi)                  ,SAVE :: GenLength                                       ! Length of the generator used for gen. GRAPHICS.
REAL(ReKi)                  ,SAVE :: GenRad                                          ! Radius of the generator used for gen. GRAPHICS.
REAL(ReKi)                  ,SAVE :: HubCylRad                                       ! Radius of hub cylincder used for hub GRAPHICS.
REAL(ReKi)                  ,SAVE :: HSSLength                                       ! Length of high-speed shaft for HSS GRAPHICS.
REAL(ReKi)                  ,SAVE :: HSSRad                                          ! Radius of the high-speed shaft used for HSS GRAPHICS.
REAL(ReKi)                  ,SAVE :: LSSLength                                       ! Length of low-speed shaft for LSS GRAPHICS.
REAL(ReKi)                  ,SAVE :: LSSRad                                          ! Radius of the low-speed shaft used for LSS GRAPHICS.
REAL(ReKi)                  ,SAVE :: NacLength                                       ! Length of nacelle used for the nacelle GRAPHICS.
REAL(ReKi)                  ,SAVE :: NacRadBot                                       ! Bottom radius of nacelle FRUSTUM used for the nacelle GRAPHICS.
REAL(ReKi)                  ,SAVE :: NacRadTop                                       ! Top    radius of nacelle FRUSTUM used for the nacelle GRAPHICS.
REAL(ReKi)                  ,SAVE :: ThkOvrChrd                                      ! Ratio of blade thickness to blade chord used for blade element GRAPHICS.
REAL(ReKi)                  ,SAVE :: TwrBaseRad                                      ! Tower base radius used for linearly tapered tower GRAPHICS.
REAL(ReKi)                  ,SAVE :: TwrTopRad                                       ! Tower top  radius used for linearly tapered tower GRAPHICS.

INTEGER(4)                  ,SAVE :: NFreeSrfc = -1                                  ! Number of points on free surface (not including the zero'th point) where the elevation of the incident waves will be computed (computed every FrSrfcSpc meters along the incident wave propogation heading direction for a length of the rotor diameter).
INTEGER(4)                  ,SAVE :: NLnNodes  = 10                                  ! Number of nodes per line for mooring line GRAPHICS.  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
INTEGER(4)                  ,SAVE :: NSides                                          ! The number of sides used in GRAPHICS CYLINDER and FRUSTUM statements.

LOGICAL                     ,SAVE :: MakeLINacf                                      ! Switch for making an ADAMS/LINEAR control command file.  To prevent an ADAMS/LINEAR control command file to be made, and to not include the RESULTS statement in the ADAMS dataset, set to .FALSE.
LOGICAL                     ,SAVE :: SaveGrphcs                                      ! Switch to determine whether or note GRAPHICS output is saved in an ADAMS analysis.



END MODULE ADAMSInput
!=======================================================================
MODULE AeroElem


   ! This MODULE stores FAST/AeroDyn interface variables.


USE                             Precision
USE                             AeroDyn  ! for type;  Precision is also included so the previous line could be removed, too.


TYPE(AllAeroMarkers)         ,SAVE :: ADAeroMarkers
TYPE(AeroLoadsOptions)       ,SAVE :: ADIntrfaceOptions
TYPE(AllAeroLoads)           ,SAVE :: ADAeroLoads
TYPE(AeroConfig)             ,SAVE :: ADInterfaceComponents                        ! The configuration markers that make up the bodies where aerodynamic calculations will be needed

INTEGER                      ,SAVE :: NumADBldNodes = 0                               ! Number of blade nodes in AeroDyn


END MODULE AeroElem
!=======================================================================
MODULE Blades


   ! This MODULE stores input variables for the blades.


USE                             Precision


REAL(ReKi)                  ,SAVE :: AdjBlMs                                         ! Factor to adjust blade mass density.
REAL(ReKi)                  ,SAVE :: AdjEdSt                                         ! Factor to adjust edge stiffness.
REAL(ReKi)                  ,SAVE :: AdjFlSt                                         ! Factor to adjust flap stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: AerCen    (:)                                   ! Aerodynamic center for distributed input data.
REAL(ReKi), ALLOCATABLE     ,SAVE :: AeroCent  (:,:)                                 ! Aerodynamic center for analysis nodes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: AeroTwst  (:)                                   ! Aerodynamic twist of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: Alpha     (:)                                   ! Blade coupling coefficient between flap and twist for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: AxRedBld  (:,:,:,:)                             ! The axial-reduction terms of the blade shape function.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BAlpha    (:,:)                                 ! Interpolated blade coupling coefficient between flap and twist.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BldEDamp  (:,:)                                 ! Blade edgewise damping coefficients.
REAL(ReKi)                  ,SAVE :: BldEdDmp  (1)                                   ! Blade structural damping ratios in edgewise direction.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BldFDamp  (:,:)                                 ! Blade flapwise damping coefficients.
REAL(ReKi)                  ,SAVE :: BldFlDmp  (2)                                   ! Blade structural damping ratios in flapwise direction.
REAL(ReKi)                  ,SAVE :: BldFlexL                                        ! Flexible blade length.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BlFract   (:)                                   ! Blade fractional radius for distributed input data.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BMassDen  (:)                                   ! Blade mass density for distributed input data.
REAL(ReKi), ALLOCATABLE     ,SAVE :: CAeroTwst (:)                                   ! Cosine of the aerodynamic twist of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: CBE       (:,:,:)                               ! Generalized edgewise damping of the blades.
REAL(ReKi), ALLOCATABLE     ,SAVE :: CBF       (:,:,:)                               ! Generalized flapwise damping of the blades.
REAL(ReKi), ALLOCATABLE     ,SAVE :: cgOffBEdg (:,:)                                 ! Interpolated blade edge (along local aerodynamic yb-axis) mass cg offset.
REAL(ReKi), ALLOCATABLE     ,SAVE :: cgOffBFlp (:,:)                                 ! Interpolated blade flap (along local aerodynamic xb-axis) mass cg offset.
REAL(ReKi), ALLOCATABLE     ,SAVE :: Chord     (:)                                   ! Chord of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: CThetaS   (:,:)                                 ! COS( ThetaS )
REAL(ReKi), ALLOCATABLE     ,SAVE :: DRNodes   (:)                                   ! Length of variable-spaced blade elements.
REAL(ReKi), ALLOCATABLE     ,SAVE :: EAOffBEdg (:,:)                                 ! Interpolated blade edge (along local aerodynamic yb-axis) elastic axis offset.
REAL(ReKi), ALLOCATABLE     ,SAVE :: EAOffBFlp (:,:)                                 ! Interpolated blade flap (along local aerodynamic xb-axis) elastic axis offset.
REAL(ReKi), ALLOCATABLE     ,SAVE :: EAStff    (:)                                   ! Blade extensional stiffness for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: EdgcgOf   (:)                                   ! Blade edge (along local aerodynamic yb-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: EdgEAOf   (:)                                   ! Blade edge (along local aerodynamic yb-axis) elastic axis offset for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: EdgIner   (:)                                   ! Blade edge (about local structural xb-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: EdgStff   (:)                                   ! Blade edge stiffness for distributed input data.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FlpcgOf   (:)                                   ! Blade flap (along local aerodynamic xb-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FlpEAOf   (:)                                   ! Blade flap (along local aerodynamic xb-axis) elastic axis offset for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FlpIner   (:)                                   ! Blade flap (about local structural yb-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FlpStff   (:)                                   ! Blade flap stiffness for distributed input data.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FStTunr   (:,:)                                 ! Blade flapwise modal stiffness tuners (stored for all blades).
REAL(ReKi)                  ,SAVE :: FlStTunr  (2)                                   ! Blade flapwise modal stiffness tuners (input).
REAL(ReKi), ALLOCATABLE     ,SAVE :: GJStff    (:)                                   ! Blade torsional stiffness for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: InerBEdg  (:,:)                                 ! Interpolated blade edge (about local structural xb-axis) mass inertia per unit length.
REAL(ReKi), ALLOCATABLE     ,SAVE :: InerBFlp  (:,:)                                 ! Interpolated blade flap (about local structural yb-axis) mass inertia per unit length.
REAL(ReKi), ALLOCATABLE     ,SAVE :: KBE       (:,:,:)                               ! Generalized edgewise stiffness of the blades.
REAL(ReKi), ALLOCATABLE     ,SAVE :: KBF       (:,:,:)                               ! Generalized flapwise stiffness of the blades.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MassB     (:,:)                                 ! Interpolated lineal blade mass density.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PrecrvRef (:)                                   ! Offset for defining the reference axis from the pitch axis for precurved blades at a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PreswpRef (:)                                   ! Offset for defining the reference axis from the pitch axis for preswept  blades at a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: RefAxisxb (:,:)                                 ! Interpolated Offset for defining the reference axis from the pitch axis for precurved blades at a given input station (along xb-axis).
REAL(ReKi), ALLOCATABLE     ,SAVE :: RefAxisyb (:,:)                                 ! Interpolated Offset for defining the reference axis from the pitch axis for preswept  blades at a given input station (along yb-axis).
REAL(ReKi), ALLOCATABLE     ,SAVE :: RNodes    (:)                                   ! Radius to analysis nodes relative to hub ( 0 < RNodes(:) < BldFlexL )
REAL(ReKi), ALLOCATABLE     ,SAVE :: RNodesNorm(:)                                   ! Normalized radius to analysis nodes relative to hub ( 0 < RNodesNorm(:) < 1 )
REAL(ReKi), ALLOCATABLE     ,SAVE :: rSAerCenn1(:,:)                                 ! Distance from point S on a blade to the aerodynamic center in the n1 direction (m).
REAL(ReKi), ALLOCATABLE     ,SAVE :: rSAerCenn2(:,:)                                 ! Distance from point S on a blade to the aerodynamic center in the n2 direction (m).
REAL(ReKi), ALLOCATABLE     ,SAVE :: SAeroTwst (:)                                   ! Sine of the aerodynamic twist of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffBE   (:,:)                                 ! Interpolated edgewise blade stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffBEA  (:,:)                                 ! Interpolated blade extensional stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffBF   (:,:)                                 ! Interpolated flapwise blade stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffBGJ  (:,:)                                 ! Interpolated blade torsional stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: SThetaS   (:,:)                                 ! SIN( ThetaS )
REAL(ReKi), ALLOCATABLE     ,SAVE :: StrcTwst  (:)                                   ! Structural twist for distributed input data.
REAL(ReKi), ALLOCATABLE     ,SAVE :: ThetaS    (:,:)                                 ! Structural twist for analysis nodes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwistedSF (:,:,:,:,:)                           ! Interpolated lineal blade mass density.

INTEGER(4)                  ,SAVE :: BldNodes                                        ! Number of blade nodes used in the analysis.
INTEGER(4)                  ,SAVE :: NBlInpSt                                        ! Number of blade input stations.
INTEGER(4)                  ,SAVE :: TipNode                                         ! Index of the additional node located at the blade tip = BldNodes + 1


END MODULE Blades
!=======================================================================
MODULE CoordSys


   ! This MODULE stores coordinate sytems used internally by FAST.  The 3
   !   components of each vector correspond to the z1, z2, and z3 components
   !   of the individual vectors.
   ! NOTE: the orientations of most of these coordinate systems will change
   !   every time step.


USE                             Precision


REAL(ReKi)                  ,SAVE :: a1       (3)                                    ! Vector / direction a1 (=  xt from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: a2       (3)                                    ! Vector / direction a2 (=  zt from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: a3       (3)                                    ! Vector / direction a3 (= -yt from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: b1       (3)                                    ! Vector / direction b1 (=  xp from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: b2       (3)                                    ! Vector / direction b2 (=  zp from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: b3       (3)                                    ! Vector / direction b3 (= -yp from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: c1       (3)                                    ! Vector / direction c1 (=  xs from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: c2       (3)                                    ! Vector / direction c2 (=  zs from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: c3       (3)                                    ! Vector / direction c3 (= -ys from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: d1       (3)                                    ! Vector / direction d1 (=  xn from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: d2       (3)                                    ! Vector / direction d2 (=  zn from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: d3       (3)                                    ! Vector / direction d3 (= -yn from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: e1       (3)                                    ! Vector / direction e1 (=  xa from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: e2       (3)                                    ! Vector / direction e2 (=  ya from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: e3       (3)                                    ! Vector / direction e3 (=  za from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: f1       (3)                                    ! Vector / direction f1.
REAL(ReKi)                  ,SAVE :: f2       (3)                                    ! Vector / direction f2.
REAL(ReKi)                  ,SAVE :: f3       (3)                                    ! Vector / direction f3.
REAL(ReKi)                  ,SAVE :: g1       (3)                                    ! Vector / direction g1 (=  xh from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: g2       (3)                                    ! Vector / direction g2 (=  yh from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: g3       (3)                                    ! Vector / direction g3 (=  zh from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: i1       (:,:)                                  ! i1(K,:) = vector / direction i1 for blade K (=  xcK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: i2       (:,:)                                  ! i2(K,:) = vector / direction i2 for blade K (=  ycK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: i3       (:,:)                                  ! i3(K,:) = vector / direction i3 for blade K (=  zcK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: j1       (:,:)                                  ! j1(K,:) = vector / direction j1 for blade K (=  xbK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: j2       (:,:)                                  ! j2(K,:) = vector / direction j2 for blade K (=  ybK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: j3       (:,:)                                  ! j3(K,:) = vector / direction j3 for blade K (=  zbK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: m1       (:,:,:)                                ! m1(K,J,:) = vector / direction m1 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
REAL(ReKi), ALLOCATABLE     ,SAVE :: m2       (:,:,:)                                ! m2(K,J,:) = vector / direction m2 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
REAL(ReKi), ALLOCATABLE     ,SAVE :: m3       (:,:,:)                                ! m3(K,J,:) = vector / direction m3 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
REAL(ReKi), ALLOCATABLE     ,SAVE :: n1       (:,:,:)                                ! n1(K,J,:) = vector / direction n1 for node J of blade K (= LxbK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: n2       (:,:,:)                                ! n2(K,J,:) = vector / direction n2 for node J of blade K (= LybK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: n3       (:,:,:)                                ! n3(K,J,:) = vector / direction n3 for node J of blade K (= LzbK from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: p1       (3)                                    ! Vector / direction p1 (used to calc. and return tail aerodynamic loads from AeroDyn).
REAL(ReKi)                  ,SAVE :: p2       (3)                                    ! Vector / direction p2 (used to calc. and return tail aerodynamic loads from AeroDyn).
REAL(ReKi)                  ,SAVE :: p3       (3)                                    ! Vector / direction p3 (used to calc. and return tail aerodynamic loads from AeroDyn).
REAL(ReKi)                  ,SAVE :: rf1      (3)                                    ! Vector / direction rf1 (rotor-furl coordinate system = d1 when rotor-furl angle = 0).
REAL(ReKi)                  ,SAVE :: rf2      (3)                                    ! Vector / direction rf2 (rotor-furl coordinate system = d2 when rotor-furl angle = 0).
REAL(ReKi)                  ,SAVE :: rf3      (3)                                    ! Vector / direction rf3 (rotor-furl coordinate system = d3 when rotor-furl angle = 0).
REAL(ReKi)                  ,SAVE :: rfa      (3)                                    ! Vector / direction of the rotor-furl axis.
REAL(ReKi), ALLOCATABLE     ,SAVE :: t1       (:,:)                                  ! Vector / direction t1 for tower node J (=  Lxt from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: t2       (:,:)                                  ! Vector / direction t2 for tower node J (=  Lzt from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: t3       (:,:)                                  ! Vector / direction t3 for tower node J (= -Lyt from the IEC coord. system).
REAL(ReKi), ALLOCATABLE     ,SAVE :: te1      (:,:,:)                                ! te1(K,J,:) = vector / direction te1 for node J of blade K (used to calc. noise).
REAL(ReKi), ALLOCATABLE     ,SAVE :: te2      (:,:,:)                                ! te2(K,J,:) = vector / direction te2 for node J of blade K (used to calc. noise).
REAL(ReKi), ALLOCATABLE     ,SAVE :: te3      (:,:,:)                                ! te3(K,J,:) = vector / direction te3 for node J of blade K (used to calc. noise).
REAL(ReKi)                  ,SAVE :: tf1      (3)                                    ! Vector / direction tf1 (tail-furl coordinate system = d1 when rotor-furl angle = 0).
REAL(ReKi)                  ,SAVE :: tf2      (3)                                    ! Vector / direction tf2 (tail-furl coordinate system = d2 when rotor-furl angle = 0).
REAL(ReKi)                  ,SAVE :: tf3      (3)                                    ! Vector / direction tf3 (tail-furl coordinate system = d3 when rotor-furl angle = 0).
REAL(ReKi)                  ,SAVE :: tfa      (3)                                    ! Vector / direction of the tail-furl axis.
REAL(ReKi)                  ,SAVE :: z1       (3)                                    ! Vector / direction z1 (=  xi from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: z2       (3)                                    ! Vector / direction z2 (=  zi from the IEC coord. system).
REAL(ReKi)                  ,SAVE :: z3       (3)                                    ! Vector / direction z3 (= -yi from the IEC coord. system).
! Yang: Start of adding TMD coordinates. 14-Aug-2018
REAL(ReKi)                  ,SAVE :: tmdnx1   (3)                                    ! Vector / direction tmdn1  (same as d coordinate system if TmdAngle = 0.)
REAL(ReKi)                  ,SAVE :: tmdnx2   (3)                                    ! Vector / direction tmdn2 
REAL(ReKi)                  ,SAVE :: tmdnx3   (3)                                    ! Vector / direction tmdn3 

REAL(ReKi)                  ,SAVE :: tmdpx1   (3)                                    ! Vector / direction tmdp1  (same as a coordinate system if TmdAngle = 0.)
REAL(ReKi)                  ,SAVE :: tmdpx2   (3)                                    ! Vector / direction tmdp2 
REAL(ReKi)                  ,SAVE :: tmdpx3   (3)                                    ! Vector / direction tmdp3 

REAL(ReKi)                  ,SAVE :: tmdny1   (3)                                    ! Vector / direction tmdn1  (same as d coordinate system if TmdAngle = 0.)
REAL(ReKi)                  ,SAVE :: tmdny2   (3)                                    ! Vector / direction tmdn2 
REAL(ReKi)                  ,SAVE :: tmdny3   (3)                                    ! Vector / direction tmdn3 

REAL(ReKi)                  ,SAVE :: tmdpy1   (3)                                    ! Vector / direction tmdp1  (same as a coordinate system if TmdAngle = 0.)
REAL(ReKi)                  ,SAVE :: tmdpy2   (3)                                    ! Vector / direction tmdp2 
REAL(ReKi)                  ,SAVE :: tmdpy3   (3)                                    ! Vector / direction tmdp3 
! Yang: End of adding TMD coordinates

END MODULE CoordSys
!=======================================================================
MODULE Constants


  ! This MODULE stores various constants.


USE                             Precision


REAL(ReKi), PARAMETER             :: Inv2Pi   =  0.15915494                          ! 0.5/Pi.
REAL(ReKi)                  ,SAVE :: TwoPiNB                                         ! 2*Pi/NumBl.  This constant is calculated in fast_io.f90/Inputs()

END MODULE Constants
!=======================================================================
MODULE DOFs


   ! This MODULE stores variables related to degrees of freedom.


USE                             Precision


REAL(ReKi), ALLOCATABLE          ,SAVE :: Q        (:,:)                                  ! Displacement matrix.
REAL(ReKi), ALLOCATABLE          ,SAVE :: QD       (:,:)                                  ! Velocity matrix.
REAL(ReKi), ALLOCATABLE     ,SAVE :: QD2      (:,:)                                  ! Acceleration matrix.

INTEGER(4), ALLOCATABLE     ,SAVE :: Diag     (:)                                    ! Array containing the indices of SrtPS() associated with each enabled DOF; that is, SrtPS(Diag(I)) = I.
INTEGER(4), ALLOCATABLE     ,SAVE :: DOF_BE   (:,:)                                  ! DOF indices for blade edge.
INTEGER(4), ALLOCATABLE     ,SAVE :: DOF_BF   (:,:)                                  ! DOF indices for blade flap.
INTEGER(4), PARAMETER             :: DOF_DrTr = 14                                   ! DOF index for drivetrain rotational-flexibility.
INTEGER(4), PARAMETER             :: DOF_GeAz = 13                                   ! DOF index for the generator azimuth.
INTEGER(4), PARAMETER             :: DOF_Hv   =  3                                   ! DOF index for platform heave.
INTEGER(4), PARAMETER             :: DOF_P    =  5                                   ! DOF index for platform pitch.
INTEGER(4), PARAMETER             :: DOF_R    =  4                                   ! DOF index for platform roll.
INTEGER(4), PARAMETER             :: DOF_RFrl = 12                                   ! DOF index for rotor-furl.
INTEGER(4), PARAMETER             :: DOF_Sg   =  1                                   ! DOF index for platform surge.
INTEGER(4), PARAMETER             :: DOF_Sw   =  2                                   ! DOF index for platform sway.
! Yang: Due to the TMDs, the DOF of Teet has to be changed to 24 from 22. 14-Aug-2018
! INTEGER(4), PARAMETER                  :: DOF_Teet = 22                                   ! DOF index for rotor-teeter.
INTEGER(4), PARAMETER             :: DOF_Teet = 24                                   ! DOF index for rotor-teeter.
! Yang: End of changing the value of DOF_Teet.
INTEGER(4), PARAMETER             :: DOF_TFA1 =  7                                   ! DOF index for 1st tower fore-aft mode.
INTEGER(4), PARAMETER             :: DOF_TFA2 =  9                                   ! DOF index for 2nd tower fore-aft mode.
INTEGER(4), PARAMETER             :: DOF_TFrl = 15                                   ! DOF index for tail-furl.
INTEGER(4), PARAMETER             :: DOF_TSS1 =  8                                   ! DOF index for 1st tower side-to-side mode.
INTEGER(4), PARAMETER             :: DOF_TSS2 = 10                                   ! DOF index for 2nd tower side-to-side mode.
INTEGER(4), PARAMETER             :: DOF_Y    =  6                                   ! DOF index for platform yaw.
INTEGER(4), PARAMETER             :: DOF_Yaw  = 11                                   ! DOF index for nacelle-yaw.
! Yang: Start of adding the DOFs indexes of the TMDs. 14-Aug-2018
INTEGER(4), PARAMETER             :: DOF_TmdX = 16                                   ! DOF index for axial TMD.
INTEGER(4), PARAMETER             :: DOF_TmdY = 17                                   ! DOF index for transverse TMD.
! Yang: End of adding the DOFs indexes of the TMDs.
INTEGER(4), ALLOCATABLE     ,SAVE :: IC       (:)                                    ! Array which stores pointers to predictor-corrector results.
INTEGER(4)                  ,SAVE :: NActvDOF                                        ! The number of active (enabled) DOFs in the model.
INTEGER(4)                  ,SAVE :: NAug                                            ! Dimension of augmented solution matrix.
INTEGER(4)                  ,SAVE :: NDOF                                            ! Number of total DOFs.
INTEGER(4), PARAMETER             :: NMX      =  9                                   ! Used in updating predictor-corrector values.
INTEGER(4)                  ,SAVE :: NPA                                             ! Number of DOFs                  that contribute to the angular velocity of the tail                                                      (body A) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPB                                             ! Number of DOFs                  that contribute to the angular velocity of the tower top / baseplate                                     (body B) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPCE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the hub center of mass                                                              (point C) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPDE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPF                                             ! Number of DOFs                  that contribute to the angular velocity of the tower elements                                            (body F) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPG                                             ! Number of DOFs                  that contribute to the angular velocity of the generator                                                 (body G) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPH                                             ! Number of DOFs                  that contribute to the angular velocity of the hub                                                       (body H) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPIE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the tail boom center of mass                                                        (point I) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPL                                             ! Number of DOFs                  that contribute to the angular velocity of the low-speed shaft                                           (body L) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPM                                             ! Number of DOFs                  that contribute to the angular velocity of the blade elements                                            (body M) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPN                                             ! Number of DOFs                  that contribute to the angular velocity of the nacelle                                                   (body N) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPTE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the tower nodes                                                                     (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPTTE                                           ! Number of tower DOFs            that contribute to the QD2T-related linear accelerations of the tower nodes                                                                     (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPR                                             ! Number of DOFs                  that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: NPSBE    (:)                                    ! Number of blade DOFs            that contribute to the QD2T-related linear accelerations of the blade nodes                                                                     (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: NPSE     (:)                                    ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the blade nodes                                                                     (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPUE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the nacelle center of mass                                                          (point U) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPX                                             ! Number of DOFs                  that contribute to the angular velocity of the platform                                                  (body X) in the inertia frame.
INTEGER(4)                  ,SAVE :: NPYE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the platform center of mass                                                         (point Y) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PA       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the tail                                                      (body A) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PB       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower top / baseplate                                     (body B) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PCE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the hub center of mass                                                              (point C) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PDE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PF       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower elements                                            (body F) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PG       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the generator                                                 (body G) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PH       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the hub                                                       (body H) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PIE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the tail boom center of mass                                                        (point I) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PL       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the low-speed shaft                                           (body L) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PM       (:,:)                                  ! Array of DOF indices (pointers) that contribute to the angular velocity of the blade elements                                            (body M) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PN       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the nacelle                                                   (body N) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PTE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the tower nodes                                                                     (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PTTE     (:)                                    ! Array of tower DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the tower nodes                                                               (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PR       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PS       (:)                                    ! Array of DOF indices (pointers) to the active (enabled) DOFs/states.
INTEGER(4), ALLOCATABLE     ,SAVE :: PSBE     (:,:)                                  ! Array of blade DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the blade nodes                                                               (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PSE      (:,:)                                  ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the blade nodes                                                                     (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PUE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the nacelle center of mass                                                          (point U) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PX       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the platform                                                  (body X) in the inertia frame.
INTEGER(4), ALLOCATABLE     ,SAVE :: PYE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the platform center of mass                                                         (point Y) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: SrtPS    (:)                                    ! Sorted (from smallest to largest DOF index) version of PS().
INTEGER(4), ALLOCATABLE     ,SAVE :: SrtPSNAUG(:)                                    ! SrtPS() with the additional value of NAUG.
! Yang: Start of creating number of DOF variables for the TMDs and pointer arrays. 14-Aug-2018
INTEGER(4)                  ,SAVE :: NPTmdXE                                         ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the TmdX center of mass                                                          (point U) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PTmdXE      (:)                                 ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the TmdX center of mass                                                          (point U) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                  ,SAVE :: NPTmdYE                                         ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the TmdY center of mass                                                          (point U) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE     ,SAVE :: PTmdYE      (:)                                 ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the TmdY center of mass   
! Yang: End of creating number of DOF variables for the TMDs and pointer arrays. 14-Aug-2018
LOGICAL,    ALLOCATABLE     ,SAVE :: DOF_Flag (:)                                    ! Array which stores values of the feature flags for each DOF.

CHARACTER(99), ALLOCATABLE  ,SAVE :: DOF_Desc (:)                                    ! Array which stores descriptions of each DOF.


END MODULE DOFs
!=======================================================================
MODULE DriveTrain


   ! This MODULE stores variables for the drivetrain.


USE                             Precision


REAL(ReKi)                  ,SAVE :: DTTorDmp                                        ! Drivetrain torsional damper
REAL(ReKi)                  ,SAVE :: DTTorSpr                                        ! Drivetrain torsional spring
REAL(ReKi)                  ,SAVE :: ElecPwr                                         ! Electrical power, W.
REAL(ReKi)                  ,SAVE :: GBRatio                                         ! Gearbox ratio
REAL(ReKi)                  ,SAVE :: GBoxEff                                         ! Gearbox efficiency.
REAL(ReKi)                  ,SAVE :: GenCTrq                                         ! Constant generator torque.
REAL(ReKi)                  ,SAVE :: GenEff                                          ! Generator efficiency
REAL(ReKi)                  ,SAVE :: GenSpRZT                                        ! Difference between rated and zero-torque generator speeds for SIG.
REAL(ReKi)                  ,SAVE :: GenSpRat                                        ! Rated generator speed.
REAL(ReKi)                  ,SAVE :: GenSpZT                                         ! Zero-torque generator speed.
REAL(ReKi)                  ,SAVE :: GenTrq                                          ! Electrical generator torque.
REAL(ReKi)                  ,SAVE :: HSSBrDT                                         ! Time it takes for HSS brake to reach full deployment once deployed.
REAL(ReKi)                  ,SAVE :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-). !bjj: used to be local variable in FAST.f90/Subroutine DrvTrTrq()
REAL(ReKi)                  ,SAVE :: HSSBrTqF                                        ! Fully deployed HSS brake torque
REAL(ReKi)                  ,SAVE :: HSSBrTrq                                        ! Instantaneous HSS brake torque
REAL(ReKi)                  ,SAVE :: HSSBrTrqC                                       ! A copy of the value of HSSBrTrq calculated in SUBROUTINE DrvTrTrq().
REAL(ReKi)                  ,SAVE :: SIG_PORt                                        ! Pull-out ratio (Tpullout/Trated).
REAL(ReKi)                  ,SAVE :: SIG_POSl                                        ! Pullout slip.
REAL(ReKi)                  ,SAVE :: SIG_POTq                                        ! Pullout torque.
REAL(ReKi)                  ,SAVE :: SIG_RtSp                                        ! Rated speed.
REAL(ReKi)                  ,SAVE :: SIG_RtTq                                        ! Rated torque.
REAL(ReKi)                  ,SAVE :: SIG_SlPc                                        ! Rated generator slip percentage.
REAL(ReKi)                  ,SAVE :: SIG_Slop                                        ! Torque/Speed slope for simple induction generator.
REAL(ReKi)                  ,SAVE :: SIG_SySp                                        ! Synchronous (zero-torque) generator speed.
REAL(ReKi)                  ,SAVE :: TEC_A0                                          ! A0 term for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_C0                                          ! C0 term for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_C1                                          ! C1 term for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_C2                                          ! C2 term for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_Freq                                        ! Line frequency for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_K1                                          ! K1 term for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_K2                                          ! K2 term for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_MR                                          ! Magnetizing reactance for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_Re1                                         ! Thevenin's equivalent stator resistance (ohms)
REAL(ReKi)                  ,SAVE :: TEC_RLR                                         ! Rotor leakage reactance for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_RRes                                        ! Rotor resistance for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_SLR                                         ! Stator leakage reactance for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_SRes                                        ! Stator resistance for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_SySp                                        ! Synchronous speed for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_V1a                                         ! Source voltage for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_VLL                                         ! Line-to-line RMS voltage for Thevenin-equivalent circuit.
REAL(ReKi)                  ,SAVE :: TEC_Xe1                                         ! Thevenin's equivalent stator leakage reactance (ohms)

INTEGER(4)                  ,SAVE :: GenDir                                          ! Direction of the generator = +/- 1 (+ 1 = same direction as LSS; -1 = opposite direction of LSS).
INTEGER(4)                  ,SAVE :: TEC_NPol                                        ! Number of poles for Thevenin-equivalent circuit.

LOGICAL                     ,SAVE :: GBRevers                                        ! Gearbox reversal flag.


END MODULE DriveTrain
!=======================================================================
MODULE EnvCond


   ! This MODULE stores input variables for environmental conditions.


USE                             Precision


REAL(ReKi)                  ,SAVE :: AirDens                                         ! Air density = RHO.
REAL(ReKi)                  ,SAVE :: CurrDIDir = 0.0                                 ! Depth-independent current heading direction.
REAL(ReKi)                  ,SAVE :: CurrDIV   = 0.0                                 ! Depth-independent current velocity.
REAL(ReKi)                  ,SAVE :: CurrNSDir = 0.0                                 ! Near-surface current heading direction.
REAL(ReKi)                  ,SAVE :: CurrNSRef = 0.0                                 ! Near-surface current reference depth.
REAL(ReKi)                  ,SAVE :: CurrNSV0  = 0.0                                 ! Near-surface current velocity at still water level.
REAL(ReKi)                  ,SAVE :: CurrSSDir = 0.0                                 ! Sub-surface current heading direction.
REAL(ReKi)                  ,SAVE :: CurrSSV0  = 0.0                                 ! Sub-surface current velocity at still water level.
REAL(ReKi)                  ,SAVE :: Gravity                                         ! Gravitational acceleration.

REAL(ReKi)                  ,SAVE :: WaveDir   = 0.0                                 ! Wave heading direction.
REAL(ReKi)                  ,SAVE :: WaveDT    = 0.0                                 ! Time step for incident wave calculations.
REAL(ReKi)                  ,SAVE :: WaveHs    = 0.0                                 ! Significant wave height.
REAL(ReKi)                  ,SAVE :: WavePkShp = 1.0                                 ! Peak shape parameter of incident wave spectrum.
REAL(ReKi)                  ,SAVE :: WaveTMax  = 0.0                                 ! Analysis time for incident wave calculations.
REAL(ReKi)                  ,SAVE :: WaveTp    = 0.0                                 ! Peak spectral period.
REAL(ReKi)                  ,SAVE :: WtrDens                                         ! Water density.
REAL(ReKi)                  ,SAVE :: WtrDpth                                         ! Water depth.

INTEGER(4)                  ,SAVE :: CurrMod                                         ! Current profile model switch.
INTEGER(4)                  ,SAVE :: WaveStMod = 0                                   ! Model switch for stretching incident wave kinematics to instantaneous free surface.
INTEGER(4)                  ,SAVE :: WaveMod   = 0                                   ! Incident wave kinematics model switch.
INTEGER(4)                  ,SAVE :: WaveSeed (2) = 0                                ! Random seeds of incident waves.

CHARACTER(1024)             ,SAVE :: GHWvFile  = ''                                  ! The root name of GH Bladed files containing wave data.


END MODULE EnvCond
!=======================================================================
MODULE Features


   ! This MODULE stores input variables for feature switches.

LOGICAL                  ,SAVE :: CompAero                                        ! Compute aerodynamic forces switch.
LOGICAL                  ,SAVE :: CompHydro = .FALSE.                             ! Compute hydrodynamic forces switch.

LOGICAL                  ,SAVE :: CompNoise                                       ! Compute aerodynamic noise  switch.
LOGICAL                  ,SAVE :: DrTrDOF                                         ! Drivetrain rotational-flexibility DOF.
LOGICAL                  ,SAVE :: EdgeDOF                                         ! Edgewise blade mode DOF.
LOGICAL                  ,SAVE :: FlapDOF1                                        ! First flapwise blade mode DOF.
LOGICAL                  ,SAVE :: FlapDOF2                                        ! Second flapwise blade mode DOF.
LOGICAL                  ,SAVE :: GenDOF                                          ! Generator DOF.
LOGICAL                  ,SAVE :: PtfmHvDOF = .FALSE.                             ! Platform vertical heave translation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                  ,SAVE :: PtfmPDOF  = .FALSE.                             ! Platform pitch tilt rotation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                  ,SAVE :: PtfmRDOF  = .FALSE.                             ! Platform roll tilt rotation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                  ,SAVE :: PtfmSgDOF = .FALSE.                             ! Platform horizontal surge translation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                  ,SAVE :: PtfmSwDOF = .FALSE.                             ! Platform horizontal sway translation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                  ,SAVE :: PtfmYDOF  = .FALSE.                             ! Platform yaw rotation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                  ,SAVE :: RFrlDOF   = .FALSE.                             ! Rotor-furl DOF. (Initialized to .FALSE. b/c not all models read in FurlFile)
LOGICAL                  ,SAVE :: TeetDOF   = .FALSE.                             ! Rotor-teeter DOF. (Initialized to .FALSE. b/c the 3-blader requires it to be .FALSE.)
LOGICAL                  ,SAVE :: TFrlDOF   = .FALSE.                             ! Tail-furl DOF. (Initialized to .FALSE. b/c not all models read in FurlFile)
LOGICAL                  ,SAVE :: TwFADOF1                                        ! First tower fore-aft bending-mode DOF.
LOGICAL                  ,SAVE :: TwFADOF2                                        ! Second tower fore-aft bending-mode DOF.
LOGICAL                  ,SAVE :: TwSSDOF1                                        ! First tower side-to-side bending-mode DOF.
LOGICAL                  ,SAVE :: TwSSDOF2                                        ! Second tower side-to-side bending-mode DOF.
LOGICAL                  ,SAVE :: YawDOF                                          ! Nacelle-yaw DOF.


END MODULE Features
!=======================================================================
MODULE General


   ! This MODULE stores input variables for general program control.


INTEGER(4)                  ,SAVE :: ADAMSPrep                                       ! ADAMS preprocessor mode {1: Run FAST, 2: use FAST as a preprocessor to create equivalent ADAMS model, 3: do both} (switch).
INTEGER(4)                  ,SAVE :: AnalMode                                        ! FAST analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch).
INTEGER(4)                  ,SAVE :: PtfmModel                                       ! Platform model {0: none, 1: onshore, 2: fixed bottom offshore, 3: floating offshore} (switch).
INTEGER(4)                  ,SAVE :: StrtTime (8)                                    ! Start time of simulation.
INTEGER(4)                  ,SAVE :: UnAC      = 24                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS SIMULATE analysis.
INTEGER(4)                  ,SAVE :: UnAD      = 23                                  ! I/O unit number for the ADAMS dataset output file (.adm).
INTEGER(4)                  ,SAVE :: UnAL      = 25                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS LINEAR analysis.
INTEGER(4)                  ,SAVE :: UnIn      = 20                                  ! I/O unit number for the input files.
INTEGER(4)                  ,SAVE :: UnLn      = 26                                  ! I/O unit number for the FAST linear output file (.lin).
INTEGER(4)                  ,SAVE :: UnNoSpec  = 27                                  ! I/O unit number for the noise spectr output file.
INTEGER(4)                  ,SAVE :: UnNoSPL   = 28                                  ! I/O unit number for the noise SPL output file.
INTEGER(4)                  ,SAVE :: UnOu      = 21                                  ! I/O unit number for the tabular output file.
INTEGER(4)                  ,SAVE :: UnOuBin   = 29                                  ! I/O unit number for the binary output file.
INTEGER(4)                  ,SAVE :: UnSu      = 22                                  ! I/O unit number for the summary output file.
! Yang: Add a Unit number for Structural control input file
INTEGER(4)                  ,SAVE :: UnSc      = 30                                  ! I/O unit number for the Structural control input file.

 
LOGICAL                     ,SAVE :: Cmpl4SFun  = .FALSE.                            ! Is FAST being compiled as an S-Function for Simulink?
LOGICAL                     ,SAVE :: Cmpl4LV    = .FALSE.                            ! Is FAST being compiled for Labview?
LOGICAL                     ,SAVE :: Furling                                         ! Read in additional model properties for furling turbine?
LOGICAL                     ,SAVE :: SumDisp                                         ! Display summary data on screen?
LOGICAL                     ,SAVE :: SumPrint                                        ! Print summary data to "*.fsm"?

CHARACTER(1024)             ,SAVE :: ADAMSFile                                       ! The name of the file containing ADAMS-specific data inputs.
CHARACTER(1024)             ,SAVE :: ADFile                                          ! The name of the AeroDyn input file.
CHARACTER(1024), ALLOCATABLE,SAVE :: BldFile  (:)                                    ! The names of the blade-data input files.
CHARACTER(1024)             ,SAVE :: DirRoot                                         ! The name of the root file including the full path to the current working directory.
CHARACTER(1024)             ,SAVE :: DynBrkFi                                        ! The name of the dynamic generator brake input file.
CHARACTER(1024)             ,SAVE :: FTitle                                          ! The title line from the primary input file.
CHARACTER(1024)             ,SAVE :: FurlFile                                        ! The name of the furling-data input file.
CHARACTER(1024)             ,SAVE :: LinFile                                         ! The name of the file containing FAST linearization control input parameters.
CHARACTER(1024)             ,SAVE :: NoiseFile                                       ! The name of the file containing aerodynamic noise input parameters.
CHARACTER(1024)             ,SAVE :: PriFile   = 'primary.fst'                       ! The name of the primary input file.  Can be overwritten on command line.
CHARACTER(1024)             ,SAVE :: PtfmFile                                        ! The name of the platform-data input file.
CHARACTER(1024)             ,SAVE :: RootName                                        ! The root name of the input and output files.
CHARACTER(1024)             ,SAVE :: TwrFile                                         ! The name of the tower-data input file.


! Yang Start of adding the defination of Structural control module in primary file of FAST. 14-Aug-2018
CHARACTER(1024)             ,SAVE :: StrcCtrlFile                                     ! The name of the seismic configuration file.
LOGICAL                     ,SAVE :: StrcCtrlMode                                   ! Consider earthquake or not.
! Yang End of adding the defination of Structural control module in primary file of FAST

END MODULE General
!=======================================================================
MODULE InitCond


   ! This MODULE stores input variables for initial conditions.


USE                             Precision


REAL(ReKi)                  ,SAVE :: Azimuth                                         ! Initial azimuth angle for blade 1.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BlPitchInit(:)                                  ! Initial blade pitch angles at the start of the simulation.
REAL(ReKi)                  ,SAVE :: IPDefl                                          ! Initial in-plane blade-tip deflection.
REAL(ReKi)                  ,SAVE :: NacYaw                                          ! Initial or fixed nacelle-yaw angle.
REAL(ReKi)                  ,SAVE :: OoPDefl                                         ! Initial out-of-plane blade-tip displacement.
REAL(ReKi)                  ,SAVE :: PtfmHeave = 0.0                                 ! Initial or fixed vertical heave translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmPitch = 0.0                                 ! Initial or fixed pitch tilt rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmRoll  = 0.0                                 ! Initial or fixed roll tilt rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmSurge = 0.0                                 ! Initial or fixed horizontal surge translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmSway  = 0.0                                 ! Initial or fixed horizontal sway translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmYaw   = 0.0                                 ! Initial or fixed yaw rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: QAzimInit                                       ! Initial value of the internal generator azimuth DOF (Q(DOF_GeAz)).
REAL(ReKi)                  ,SAVE :: RotFurl   = 0.0                                 ! Initial or fixed rotor-furl angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RotSpeed                                        ! Initial or fixed rotor speed.
REAL(ReKi)                  ,SAVE :: TailFurl  = 0.0                                 ! Initial or fixed tail-furl angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TTDspFA                                         ! Initial fore-aft tower-top displacement.
REAL(ReKi)                  ,SAVE :: TTDspSS                                         ! Initial side-to-side tower-top displacement.
REAL(ReKi)                  ,SAVE :: TeetDefl  = 0.0                                 ! Initial or fixed teeter angle. (Initialized to zero b/c the 3-blader requires it to be zero)


! Yang: Add the definations of initial variables of the TMDs
REAL(ReKi)                  ,SAVE :: TmdXDsp   = 0.0                                 ! Initial or fixed displacement of the axial TMD. (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdXDspInit = 0.0                               ! Initial or fixed displacement of the axial TMD. (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdXSprInit = 0.0                               ! Initial or fixed spring stiffness of the axial TMD. (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdXDampInit  = 0.0                             ! Initial or fixed damping of the axial TMD. (Initialized to zero b/c not all models use a TMD)

REAL(ReKi)                  ,SAVE :: TmdYDsp   = 0.0                                 ! Initial or fixed displacement of the transverse TMD. (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdYDspInit = 0.0                               ! Initial or fixed displacement of the axial TMD. (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdYSprInit = 0.0                               ! Initial or fixed spring stiffness of the transverse TMD. (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdYDampInit  = 0.0                             ! Initial or fixed damping of the transverse TMD. (Initialized to zero b/c not all models use a TMD)

REAL(ReKi)                  ,SAVE :: TmdXFextInit  = 0.0                             ! Initial or fixed external force of the axial TMD. (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdYFextInit  = 0.0                             ! Initial or fixed external force of the side-side TMD. (Initialized to zero b/c not all models use a TMD)

REAL(ReKi)                  ,SAVE :: TmdXAngleInit  = 0.0                             ! Initial or fixed TmdX rotation angle about vertical axis in degrees (Initialized to zero b/c not all models use a TMD)
REAL(ReKi)                  ,SAVE :: TmdYAngleInit  = 0.0                             ! Initial or fixed TmdY rotation angle about vertical axis in degrees (Initialized to zero b/c not all models use a TMD)
! Yang: End of the change 
LOGICAL,    ALLOCATABLE     ,SAVE :: DOF_FlagInit(:)                                 ! Array which stores initial values of the feature flags for each DOF (at the start of the simulation).


END MODULE InitCond
!=======================================================================
MODULE Linear


   ! This MODULE stores variables for a FAST linearization analysis.


USE                             Precision


REAL(ReKi)                  ,SAVE :: AbsQDNorm = 0.0                                 ! 2-norm of the absolute difference between the velocites     of two consecutive periods.
REAL(ReKi)                  ,SAVE :: AbsQNorm  = 0.0                                 ! 2-norm of the absolute difference between the displacements of two consecutive periods.
REAL(ReKi)                  ,SAVE :: DelGenTrq = 0.0                                 ! Pertubation in generator torque using during FAST linearization (zero otherwise).
REAL(ReKi)                  ,SAVE :: DispTol                                         ! Convergence tolerance for the 2-norm of the absolute difference between the displacements of two consecutive periods (rad).
REAL(ReKi)                  ,SAVE :: Period                                          ! Steady state period of solution.
REAL(ReKi), ALLOCATABLE     ,SAVE :: QD2op    (:,:)                                  ! Periodic steady state operating accelerations.
REAL(ReKi), ALLOCATABLE     ,SAVE :: QDop     (:,:)                                  ! Periodic steady state operating velocities.
REAL(ReKi), ALLOCATABLE     ,SAVE :: Qop      (:,:)                                  ! Periodic steady state operating displacements.
REAL(ReKi)                  ,SAVE :: VelTol                                          ! Convergence tolerance for the 2-norm of the absolute difference between the velocities    of two consecutive periods (rad/s).

INTEGER(4)                  ,SAVE :: CntrlInpt(7)                                    ! List   of control inputs [1 to NInputs] {1: nacelle yaw angle, 2: nacelle yaw rate, 3: generator torque, 4: collective blade pitch, 5: individual pitch of blade 1, 6: individual pitch of blade 2, 7: individual pitch of blade 3 [unavailable for 2-bladed turbines]} (-) [unused if NInputs=0]
INTEGER(4)                  ,SAVE :: Disturbnc(7)                                    ! List   of input wind disturbances [1 to NDisturbs] {1: horizontal hub-height wind speed, 2: horizontal wind direction, 3: vertical wind speed, 4: horizontal wind shear, 5: vertical power law wind shear, 6: linear vertical wind shear, 7: horizontal hub-height wind gust} (-) [unused if NDisturbs=0]
INTEGER(4)                  ,SAVE :: Iteration = 0                                   ! Current iteration (number of periods to convergence)
INTEGER(4)                  ,SAVE :: MdlOrder                                        ! Order of output linearized model (1: 1st order A, B, Bd; 2: 2nd order M, C, K, F, Fd) (switch)
INTEGER(4)                  ,SAVE :: NAzimStep                                       ! Number of azimuth steps in periodic linearized model (-).
INTEGER(4)                  ,SAVE :: NDisturbs                                       ! Number of wind disturbances [0 to 7] (-)
INTEGER(4)                  ,SAVE :: NInputs                                         ! Number of control inputs [0 (none) or 1 to 4+NumBl] (-)
INTEGER(4)                  ,SAVE :: NStep                                           ! Number of time steps in one Period.
INTEGER(4)                  ,SAVE :: TrimCase                                        ! Trim case {1: find nacelle yaw, 2: find generator torque, 3: find collective blade pitch} (switch) [used only when CalcStdy=True and GenDOF=True]

LOGICAL                     ,SAVE :: CalcStdy                                        ! Calculate periodic steady state condition (False: linearize about zero) (switch).
LOGICAL                     ,SAVE :: IgnoreMOD = .FALSE.                             ! Ignore the use of function MOD in SUBROUTINE CalcOuts()?


END MODULE Linear
!=======================================================================
MODULE MassInert


   ! This MODULE stores input variables for turbine mass and inertias.


USE                             Precision


REAL(ReKi)                  ,SAVE :: AtfaIner                                        ! Inertia of tail boom about the tail-furl axis whose origin is the tail boom center of mass.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BldCG     (:)                                   ! Blade center of mass wrt the blade root.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BldMass   (:)                                   ! Blade masses
REAL(ReKi)                  ,SAVE :: BoomMass  = 0.0                                 ! Tail boom mass. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi), ALLOCATABLE     ,SAVE :: FirstMom  (:)                                   ! First mass moment of inertia of blades wrt the root.
REAL(ReKi)                  ,SAVE :: GenIner                                         ! Generator inertia about HSS.
REAL(ReKi)                  ,SAVE :: Hubg1Iner                                       ! Inertia of hub about g1-axis (rotor centerline).
REAL(ReKi)                  ,SAVE :: Hubg2Iner                                       ! Inertia of hub about g2-axis (transverse to the cyclinder and passing through its c.g.).
REAL(ReKi)                  ,SAVE :: HubIner                                         ! Hub inertia about teeter axis (2-blader) or rotor axis (3-blader).
REAL(ReKi)                  ,SAVE :: HubMass                                         ! Hub mass.
REAL(ReKi)                  ,SAVE :: Nacd2Iner                                       ! Inertia of nacelle about the d2-axis whose origin is the nacelle center of mass.
REAL(ReKi)                  ,SAVE :: NacMass                                         ! Nacelle mass.
REAL(ReKi)                  ,SAVE :: NacYIner                                        ! Nacelle yaw inertia.
REAL(ReKi)                  ,SAVE :: PtfmMass  = 0.0                                 ! Platform mass. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmPIner = 0.0                                 ! Platform inertia for pitch tilt rotation about the platform CM. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmRIner = 0.0                                 ! Platform inertia for roll tilt rotation about the platform CM. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmYIner = 0.0                                 ! Platform inertia for yaw rotation about the platform CM. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                  ,SAVE :: RFrlIner  = 0.0                                 ! Rotor-furl inertia about rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlMass  = 0.0                                 ! Rotor-furl mass. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RotIner                                         ! Inertia of rotor about its centerline.
REAL(ReKi)                  ,SAVE :: RotMass                                         ! Rotor mass (blades, tips, and hub)
REAL(ReKi)                  ,SAVE :: RrfaIner                                        ! Inertia of structure that furls with the rotor (not including rotor) about the rotor-furl axis whose origin is the center of mass of the structure that furls with the rotor (not including rotor).
REAL(ReKi), ALLOCATABLE     ,SAVE :: SecondMom (:)                                   ! Second mass moment of inertia of blades wrt the root.
REAL(ReKi)                  ,SAVE :: TFinMass  = 0.0                                 ! Tail fin mass. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlIner  = 0.0                                 ! Tail boom inertia about tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi), ALLOCATABLE     ,SAVE :: TipMass   (:)                                   ! Tip-brake masses.
REAL(ReKi)                  ,SAVE :: TotalMass                                       ! Mass of turbine + platform.
REAL(ReKi)                  ,SAVE :: TurbMass                                        ! Mass of turbine (tower + rotor + nacelle).
REAL(ReKi)                  ,SAVE :: TwrMass                                         ! Mass of tower.
REAL(ReKi)                  ,SAVE :: TwrTpMass                                       ! Tower-top mass (rotor + nacelle).
REAL(ReKi)                  ,SAVE :: YawBrMass                                       ! Yaw bearing mass.
! Yang: add the  difinations of  tower modal mass
REAL(ReKi)                  ,SAVE :: MTFA      (2,2)                                 ! Generalized fore-aft mass of the tower.
REAL(ReKi)                  ,SAVE :: MTSS      (2,2)                                 ! Generalized side-to-side mass of the tower.
REAL(ReKi)                  ,SAVE :: MTFA_eq   (2,2) = 0.0                           ! Generalized fore-aft mass of the tower.
REAL(ReKi)                  ,SAVE :: MTSS_eq   (2,2) = 0.0                           ! Generalized side-to-side mass of the tower.


END MODULE MassInert
!=======================================================================
MODULE Modes


   ! This MODULE stores variables for mode shapes and CONTAINs FUNCTION SHP.


USE                             Precision


INTEGER(4), PARAMETER             :: PolyOrd = 6                                     ! Order of the polynomial describing the mode shape.

REAL(ReKi), ALLOCATABLE     ,SAVE :: BldEdgSh (:,:)                                  ! Blade-edge-mode shape coefficients.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BldFl1Sh (:,:)                                  ! Blade-flap-mode-1 shape coefficients.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BldFl2Sh (:,:)                                  ! Blade-flap-mode-2 shape coefficients.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FreqBE   (:,:,:)                                ! Blade edgewise natural frequencies (both w/ and w/o centrifugal stiffening)
REAL(ReKi), ALLOCATABLE     ,SAVE :: FreqBF   (:,:,:)                                ! Blade flapwise natural frequencies (both w/ and w/o centrifugal stiffening)
REAL(ReKi)                  ,SAVE :: FreqTFA  (2,2)                                  ! Computed fore-aft tower natural frequencies.
REAL(ReKi)                  ,SAVE :: FreqTSS  (2,2)                                  ! Computed side-to-side tower natural frequencies.
REAL(ReKi)                  ,SAVE :: TwFAM1Sh (2:PolyOrd)                            ! Tower fore-aft mode-1 shape coefficients.
REAL(ReKi)                  ,SAVE :: TwFAM2Sh (2:PolyOrd)                            ! Tower fore-aft mode-2 shape coefficients.    
REAL(ReKi)                  ,SAVE :: TwSSM1Sh (2:PolyOrd)                            ! Tower side-to-side mode-1 shape coefficients.
REAL(ReKi)                  ,SAVE :: TwSSM2Sh (2:PolyOrd)                            ! Tower side-to-side mode-2 shape coefficients.

LOGICAL                     ,SAVE :: CalcBMode                                       ! T: calculate blade mode shapes internally, F: use blade mode shapes from the blade file.
LOGICAL,    ALLOCATABLE     ,SAVE :: CalcBModes(:)                                   ! Holds CalcBMode for all of the blades.
LOGICAL                     ,SAVE :: CalcTMode                                       ! T: calculate tower mode shapes internally, F: use tower mode shapes from the tower file.


CONTAINS
!=======================================================================

   FUNCTION SHP(Fract, FlexL, ModShpAry, Deriv)

   ! SHP calculates the Derive-derivative of the shape function ModShpAry at Fract.
   ! NOTE: This function only works for Deriv = 0, 1, or 2.

      USE                           Precision


      IMPLICIT                      NONE


   ! Passed variables:

      REAL(ReKi), INTENT(IN )    :: FlexL          ! Length of flexible beam, (m)
      REAL(ReKi), INTENT(IN )    :: Fract          ! Fractional distance along flexible beam, 0<=Frac<=1
      REAL(ReKi), INTENT(IN )    :: ModShpAry(2:PolyOrd)                        ! Array holding mode shape coefficients
      REAL(ReKi)                 :: SHP            ! The shape function returned by this function.

      INTEGER(4), INTENT(IN )    :: Deriv          ! Which derivative to compute Deriv = 0 (regular function SHP), 1 (D(SHP)/DZ), 2 (D2(SHP)/DZ2)


   ! Lccal variables:

      INTEGER(4)                 :: CoefTmp        !Temporary coefficient
      INTEGER(4)                 :: I              !Counts through polynomial array.
      INTEGER(4)                 :: Swtch(0:2)     !Corresponds to which derivative to compute.  Sets all portions of the coefficient = 0 except those that are relevant.


      Swtch        = 0 !Initialize Swtch(:) to 0
      Swtch(Deriv) = 1
      SHP          = 0.0

      DO I = 2,PolyOrd
         CoefTmp = Swtch(0) + ( Swtch(1)*I ) + ( Swtch(2)*I*( I - 1 ) )

         IF ( (I == 2) .AND. (Deriv == 2) ) THEN
            SHP = ModShpAry(I)*CoefTmp/( FlexL**Deriv )
         ELSE
            SHP = SHP + ModShpAry(I)*CoefTmp*( Fract**( I - Deriv ) )/( FlexL**Deriv )
         ENDIF
      ENDDO !I

      RETURN
   END FUNCTION SHP
!=======================================================================
END MODULE Modes
!=======================================================================
MODULE NacelleYaw


   ! This MODULE stores variables for nacelle yaw.


USE                             Precision


REAL(ReKi)                  ,SAVE :: YawSpr                                          ! Nacelle-yaw spring constant.
REAL(ReKi)                  ,SAVE :: YawDamp                                         ! Nacelle-yaw constant.
REAL(ReKi)                  ,SAVE :: YawNeut                                         ! Neutral yaw position.
REAL(ReKi)                  ,SAVE :: YawRateNeut = 0.0                               ! Neutral yaw rate.


END MODULE NacelleYaw
!=======================================================================
MODULE Output


   ! This MODULE stores variables used for output.


USE                             NWTC_Library

! ==================================================================================================="
! NOTE: The following lines of code were generated by a Matlab script called "Write_ChckOutLst.m"
!      using the parameters listed in the "OutListParameters.xlsx" Excel file. Any changes to these 
!      lines should be modified in the Matlab script and/or Excel worksheet as necessary. 
! ==================================================================================================="
! This code was generated by Write_ChckOutLst.m at 30-Jan-2012 12:30:08.
   TYPE     OutParmType                                                          ! User-defined type for output parameters
      INTEGER                   :: Indx                                          ! Index into AllOuts output array
      CHARACTER(10)             :: Name                                          ! Name of the output parameter.
      CHARACTER(10)             :: Units                                         ! Units corresponding to the output parameter.
      INTEGER                   :: SignM                                         ! sign (output multiplier).
   END TYPE OutParmType


     ! Parameters:


     ! Indices for computing output channels:
     ! NOTES: 
     !    (1) These parameters are in the order stored in "OutListParameters.xlsx"
     !    (2) Array AllOuts() must be dimensioned to the value of the largest output parameter
     !    (3) If an index (MaxOutPts) ever becomes greater or equal to 1000, the logic to create ARRAY/1 in the FAST-to-ADAMS preprocessor will have to be changed.

     !  Time: 

   INTEGER, PARAMETER                  :: Time      =    0


  ! Wind Motions:

   INTEGER, PARAMETER                  :: WindVxi   =    1
   INTEGER, PARAMETER                  :: WindVyi   =    2
   INTEGER, PARAMETER                  :: WindVzi   =    3
   INTEGER, PARAMETER                  :: TotWindV  =    4
   INTEGER, PARAMETER                  :: HorWindV  =    5
   INTEGER, PARAMETER                  :: HorWndDir =    6
   INTEGER, PARAMETER                  :: VerWndDir =    7


  ! Blade 1 Tip Motions:

   INTEGER, PARAMETER                  :: TipDxc1   =    8
   INTEGER, PARAMETER                  :: TipDyc1   =    9
   INTEGER, PARAMETER                  :: TipDzc1   =   10
   INTEGER, PARAMETER                  :: TipDxb1   =   11
   INTEGER, PARAMETER                  :: TipDyb1   =   12
   INTEGER, PARAMETER                  :: TipALxb1  =   13
   INTEGER, PARAMETER                  :: TipALyb1  =   14
   INTEGER, PARAMETER                  :: TipALzb1  =   15
   INTEGER, PARAMETER                  :: TipRDxb1  =   16
   INTEGER, PARAMETER                  :: TipRDyb1  =   17
   INTEGER, PARAMETER                  :: TipRDzc1  =   18
   INTEGER, PARAMETER                  :: TipClrnc1 =   19


  ! Blade 2 Tip Motions:

   INTEGER, PARAMETER                  :: TipDxc2   =   20
   INTEGER, PARAMETER                  :: TipDyc2   =   21
   INTEGER, PARAMETER                  :: TipDzc2   =   22
   INTEGER, PARAMETER                  :: TipDxb2   =   23
   INTEGER, PARAMETER                  :: TipDyb2   =   24
   INTEGER, PARAMETER                  :: TipALxb2  =   25
   INTEGER, PARAMETER                  :: TipALyb2  =   26
   INTEGER, PARAMETER                  :: TipALzb2  =   27
   INTEGER, PARAMETER                  :: TipRDxb2  =   28
   INTEGER, PARAMETER                  :: TipRDyb2  =   29
   INTEGER, PARAMETER                  :: TipRDzc2  =   30
   INTEGER, PARAMETER                  :: TipClrnc2 =   31


  ! Blade 3 Tip Motions:

   INTEGER, PARAMETER                  :: TipDxc3   =   32
   INTEGER, PARAMETER                  :: TipDyc3   =   33
   INTEGER, PARAMETER                  :: TipDzc3   =   34
   INTEGER, PARAMETER                  :: TipDxb3   =   35
   INTEGER, PARAMETER                  :: TipDyb3   =   36
   INTEGER, PARAMETER                  :: TipALxb3  =   37
   INTEGER, PARAMETER                  :: TipALyb3  =   38
   INTEGER, PARAMETER                  :: TipALzb3  =   39
   INTEGER, PARAMETER                  :: TipRDxb3  =   40
   INTEGER, PARAMETER                  :: TipRDyb3  =   41
   INTEGER, PARAMETER                  :: TipRDzc3  =   42
   INTEGER, PARAMETER                  :: TipClrnc3 =   43


  ! Blade 1 Local Span Motions:

   INTEGER, PARAMETER                  :: Spn1ALxb1 =   44
   INTEGER, PARAMETER                  :: Spn1ALyb1 =   45
   INTEGER, PARAMETER                  :: Spn1ALzb1 =   46
   INTEGER, PARAMETER                  :: Spn2ALxb1 =   47
   INTEGER, PARAMETER                  :: Spn2ALyb1 =   48
   INTEGER, PARAMETER                  :: Spn2ALzb1 =   49
   INTEGER, PARAMETER                  :: Spn3ALxb1 =   50
   INTEGER, PARAMETER                  :: Spn3ALyb1 =   51
   INTEGER, PARAMETER                  :: Spn3ALzb1 =   52
   INTEGER, PARAMETER                  :: Spn4ALxb1 =   53
   INTEGER, PARAMETER                  :: Spn4ALyb1 =   54
   INTEGER, PARAMETER                  :: Spn4ALzb1 =   55
   INTEGER, PARAMETER                  :: Spn5ALxb1 =   56
   INTEGER, PARAMETER                  :: Spn5ALyb1 =   57
   INTEGER, PARAMETER                  :: Spn5ALzb1 =   58
   INTEGER, PARAMETER                  :: Spn6ALxb1 =   59
   INTEGER, PARAMETER                  :: Spn6ALyb1 =   60
   INTEGER, PARAMETER                  :: Spn6ALzb1 =   61
   INTEGER, PARAMETER                  :: Spn7ALxb1 =   62
   INTEGER, PARAMETER                  :: Spn7ALyb1 =   63
   INTEGER, PARAMETER                  :: Spn7ALzb1 =   64
   INTEGER, PARAMETER                  :: Spn8ALxb1 =   65
   INTEGER, PARAMETER                  :: Spn8ALyb1 =   66
   INTEGER, PARAMETER                  :: Spn8ALzb1 =   67
   INTEGER, PARAMETER                  :: Spn9ALxb1 =   68
   INTEGER, PARAMETER                  :: Spn9ALyb1 =   69
   INTEGER, PARAMETER                  :: Spn9ALzb1 =   70
   INTEGER, PARAMETER                  :: Spn1TDxb1 =   71
   INTEGER, PARAMETER                  :: Spn1TDyb1 =   72
   INTEGER, PARAMETER                  :: Spn1TDzb1 =   73
   INTEGER, PARAMETER                  :: Spn2TDxb1 =   74
   INTEGER, PARAMETER                  :: Spn2TDyb1 =   75
   INTEGER, PARAMETER                  :: Spn2TDzb1 =   76
   INTEGER, PARAMETER                  :: Spn3TDxb1 =   77
   INTEGER, PARAMETER                  :: Spn3TDyb1 =   78
   INTEGER, PARAMETER                  :: Spn3TDzb1 =   79
   INTEGER, PARAMETER                  :: Spn4TDxb1 =   80
   INTEGER, PARAMETER                  :: Spn4TDyb1 =   81
   INTEGER, PARAMETER                  :: Spn4TDzb1 =   82
   INTEGER, PARAMETER                  :: Spn5TDxb1 =   83
   INTEGER, PARAMETER                  :: Spn5TDyb1 =   84
   INTEGER, PARAMETER                  :: Spn5TDzb1 =   85
   INTEGER, PARAMETER                  :: Spn6TDxb1 =   86
   INTEGER, PARAMETER                  :: Spn6TDyb1 =   87
   INTEGER, PARAMETER                  :: Spn6TDzb1 =   88
   INTEGER, PARAMETER                  :: Spn7TDxb1 =   89
   INTEGER, PARAMETER                  :: Spn7TDyb1 =   90
   INTEGER, PARAMETER                  :: Spn7TDzb1 =   91
   INTEGER, PARAMETER                  :: Spn8TDxb1 =   92
   INTEGER, PARAMETER                  :: Spn8TDyb1 =   93
   INTEGER, PARAMETER                  :: Spn8TDzb1 =   94
   INTEGER, PARAMETER                  :: Spn9TDxb1 =   95
   INTEGER, PARAMETER                  :: Spn9TDyb1 =   96
   INTEGER, PARAMETER                  :: Spn9TDzb1 =   97
   INTEGER, PARAMETER                  :: Spn1RDxb1 =   98
   INTEGER, PARAMETER                  :: Spn1RDyb1 =   99
   INTEGER, PARAMETER                  :: Spn1RDzb1 =  100
   INTEGER, PARAMETER                  :: Spn2RDxb1 =  101
   INTEGER, PARAMETER                  :: Spn2RDyb1 =  102
   INTEGER, PARAMETER                  :: Spn2RDzb1 =  103
   INTEGER, PARAMETER                  :: Spn3RDxb1 =  104
   INTEGER, PARAMETER                  :: Spn3RDyb1 =  105
   INTEGER, PARAMETER                  :: Spn3RDzb1 =  106
   INTEGER, PARAMETER                  :: Spn4RDxb1 =  107
   INTEGER, PARAMETER                  :: Spn4RDyb1 =  108
   INTEGER, PARAMETER                  :: Spn4RDzb1 =  109
   INTEGER, PARAMETER                  :: Spn5RDxb1 =  110
   INTEGER, PARAMETER                  :: Spn5RDyb1 =  111
   INTEGER, PARAMETER                  :: Spn5RDzb1 =  112
   INTEGER, PARAMETER                  :: Spn6RDxb1 =  113
   INTEGER, PARAMETER                  :: Spn6RDyb1 =  114
   INTEGER, PARAMETER                  :: Spn6RDzb1 =  115
   INTEGER, PARAMETER                  :: Spn7RDxb1 =  116
   INTEGER, PARAMETER                  :: Spn7RDyb1 =  117
   INTEGER, PARAMETER                  :: Spn7RDzb1 =  118
   INTEGER, PARAMETER                  :: Spn8RDxb1 =  119
   INTEGER, PARAMETER                  :: Spn8RDyb1 =  120
   INTEGER, PARAMETER                  :: Spn8RDzb1 =  121
   INTEGER, PARAMETER                  :: Spn9RDxb1 =  122
   INTEGER, PARAMETER                  :: Spn9RDyb1 =  123
   INTEGER, PARAMETER                  :: Spn9RDzb1 =  124


  ! Blade 2 Local Span Motions:

   INTEGER, PARAMETER                  :: Spn1ALxb2 =  125
   INTEGER, PARAMETER                  :: Spn1ALyb2 =  126
   INTEGER, PARAMETER                  :: Spn1ALzb2 =  127
   INTEGER, PARAMETER                  :: Spn2ALxb2 =  128
   INTEGER, PARAMETER                  :: Spn2ALyb2 =  129
   INTEGER, PARAMETER                  :: Spn2ALzb2 =  130
   INTEGER, PARAMETER                  :: Spn3ALxb2 =  131
   INTEGER, PARAMETER                  :: Spn3ALyb2 =  132
   INTEGER, PARAMETER                  :: Spn3ALzb2 =  133
   INTEGER, PARAMETER                  :: Spn4ALxb2 =  134
   INTEGER, PARAMETER                  :: Spn4ALyb2 =  135
   INTEGER, PARAMETER                  :: Spn4ALzb2 =  136
   INTEGER, PARAMETER                  :: Spn5ALxb2 =  137
   INTEGER, PARAMETER                  :: Spn5ALyb2 =  138
   INTEGER, PARAMETER                  :: Spn5ALzb2 =  139
   INTEGER, PARAMETER                  :: Spn6ALxb2 =  140
   INTEGER, PARAMETER                  :: Spn6ALyb2 =  141
   INTEGER, PARAMETER                  :: Spn6ALzb2 =  142
   INTEGER, PARAMETER                  :: Spn7ALxb2 =  143
   INTEGER, PARAMETER                  :: Spn7ALyb2 =  144
   INTEGER, PARAMETER                  :: Spn7ALzb2 =  145
   INTEGER, PARAMETER                  :: Spn8ALxb2 =  146
   INTEGER, PARAMETER                  :: Spn8ALyb2 =  147
   INTEGER, PARAMETER                  :: Spn8ALzb2 =  148
   INTEGER, PARAMETER                  :: Spn9ALxb2 =  149
   INTEGER, PARAMETER                  :: Spn9ALyb2 =  150
   INTEGER, PARAMETER                  :: Spn9ALzb2 =  151
   INTEGER, PARAMETER                  :: Spn1TDxb2 =  152
   INTEGER, PARAMETER                  :: Spn1TDyb2 =  153
   INTEGER, PARAMETER                  :: Spn1TDzb2 =  154
   INTEGER, PARAMETER                  :: Spn2TDxb2 =  155
   INTEGER, PARAMETER                  :: Spn2TDyb2 =  156
   INTEGER, PARAMETER                  :: Spn2TDzb2 =  157
   INTEGER, PARAMETER                  :: Spn3TDxb2 =  158
   INTEGER, PARAMETER                  :: Spn3TDyb2 =  159
   INTEGER, PARAMETER                  :: Spn3TDzb2 =  160
   INTEGER, PARAMETER                  :: Spn4TDxb2 =  161
   INTEGER, PARAMETER                  :: Spn4TDyb2 =  162
   INTEGER, PARAMETER                  :: Spn4TDzb2 =  163
   INTEGER, PARAMETER                  :: Spn5TDxb2 =  164
   INTEGER, PARAMETER                  :: Spn5TDyb2 =  165
   INTEGER, PARAMETER                  :: Spn5TDzb2 =  166
   INTEGER, PARAMETER                  :: Spn6TDxb2 =  167
   INTEGER, PARAMETER                  :: Spn6TDyb2 =  168
   INTEGER, PARAMETER                  :: Spn6TDzb2 =  169
   INTEGER, PARAMETER                  :: Spn7TDxb2 =  170
   INTEGER, PARAMETER                  :: Spn7TDyb2 =  171
   INTEGER, PARAMETER                  :: Spn7TDzb2 =  172
   INTEGER, PARAMETER                  :: Spn8TDxb2 =  173
   INTEGER, PARAMETER                  :: Spn8TDyb2 =  174
   INTEGER, PARAMETER                  :: Spn8TDzb2 =  175
   INTEGER, PARAMETER                  :: Spn9TDxb2 =  176
   INTEGER, PARAMETER                  :: Spn9TDyb2 =  177
   INTEGER, PARAMETER                  :: Spn9TDzb2 =  178
   INTEGER, PARAMETER                  :: Spn1RDxb2 =  179
   INTEGER, PARAMETER                  :: Spn1RDyb2 =  180
   INTEGER, PARAMETER                  :: Spn1RDzb2 =  181
   INTEGER, PARAMETER                  :: Spn2RDxb2 =  182
   INTEGER, PARAMETER                  :: Spn2RDyb2 =  183
   INTEGER, PARAMETER                  :: Spn2RDzb2 =  184
   INTEGER, PARAMETER                  :: Spn3RDxb2 =  185
   INTEGER, PARAMETER                  :: Spn3RDyb2 =  186
   INTEGER, PARAMETER                  :: Spn3RDzb2 =  187
   INTEGER, PARAMETER                  :: Spn4RDxb2 =  188
   INTEGER, PARAMETER                  :: Spn4RDyb2 =  189
   INTEGER, PARAMETER                  :: Spn4RDzb2 =  190
   INTEGER, PARAMETER                  :: Spn5RDxb2 =  191
   INTEGER, PARAMETER                  :: Spn5RDyb2 =  192
   INTEGER, PARAMETER                  :: Spn5RDzb2 =  193
   INTEGER, PARAMETER                  :: Spn6RDxb2 =  194
   INTEGER, PARAMETER                  :: Spn6RDyb2 =  195
   INTEGER, PARAMETER                  :: Spn6RDzb2 =  196
   INTEGER, PARAMETER                  :: Spn7RDxb2 =  197
   INTEGER, PARAMETER                  :: Spn7RDyb2 =  198
   INTEGER, PARAMETER                  :: Spn7RDzb2 =  199
   INTEGER, PARAMETER                  :: Spn8RDxb2 =  200
   INTEGER, PARAMETER                  :: Spn8RDyb2 =  201
   INTEGER, PARAMETER                  :: Spn8RDzb2 =  202
   INTEGER, PARAMETER                  :: Spn9RDxb2 =  203
   INTEGER, PARAMETER                  :: Spn9RDyb2 =  204
   INTEGER, PARAMETER                  :: Spn9RDzb2 =  205


  ! Blade 3 Local Span Motions:

   INTEGER, PARAMETER                  :: Spn1ALxb3 =  206
   INTEGER, PARAMETER                  :: Spn1ALyb3 =  207
   INTEGER, PARAMETER                  :: Spn1ALzb3 =  208
   INTEGER, PARAMETER                  :: Spn2ALxb3 =  209
   INTEGER, PARAMETER                  :: Spn2ALyb3 =  210
   INTEGER, PARAMETER                  :: Spn2ALzb3 =  211
   INTEGER, PARAMETER                  :: Spn3ALxb3 =  212
   INTEGER, PARAMETER                  :: Spn3ALyb3 =  213
   INTEGER, PARAMETER                  :: Spn3ALzb3 =  214
   INTEGER, PARAMETER                  :: Spn4ALxb3 =  215
   INTEGER, PARAMETER                  :: Spn4ALyb3 =  216
   INTEGER, PARAMETER                  :: Spn4ALzb3 =  217
   INTEGER, PARAMETER                  :: Spn5ALxb3 =  218
   INTEGER, PARAMETER                  :: Spn5ALyb3 =  219
   INTEGER, PARAMETER                  :: Spn5ALzb3 =  220
   INTEGER, PARAMETER                  :: Spn6ALxb3 =  221
   INTEGER, PARAMETER                  :: Spn6ALyb3 =  222
   INTEGER, PARAMETER                  :: Spn6ALzb3 =  223
   INTEGER, PARAMETER                  :: Spn7ALxb3 =  224
   INTEGER, PARAMETER                  :: Spn7ALyb3 =  225
   INTEGER, PARAMETER                  :: Spn7ALzb3 =  226
   INTEGER, PARAMETER                  :: Spn8ALxb3 =  227
   INTEGER, PARAMETER                  :: Spn8ALyb3 =  228
   INTEGER, PARAMETER                  :: Spn8ALzb3 =  229
   INTEGER, PARAMETER                  :: Spn9ALxb3 =  230
   INTEGER, PARAMETER                  :: Spn9ALyb3 =  231
   INTEGER, PARAMETER                  :: Spn9ALzb3 =  232
   INTEGER, PARAMETER                  :: Spn1TDxb3 =  233
   INTEGER, PARAMETER                  :: Spn1TDyb3 =  234
   INTEGER, PARAMETER                  :: Spn1TDzb3 =  235
   INTEGER, PARAMETER                  :: Spn2TDxb3 =  236
   INTEGER, PARAMETER                  :: Spn2TDyb3 =  237
   INTEGER, PARAMETER                  :: Spn2TDzb3 =  238
   INTEGER, PARAMETER                  :: Spn3TDxb3 =  239
   INTEGER, PARAMETER                  :: Spn3TDyb3 =  240
   INTEGER, PARAMETER                  :: Spn3TDzb3 =  241
   INTEGER, PARAMETER                  :: Spn4TDxb3 =  242
   INTEGER, PARAMETER                  :: Spn4TDyb3 =  243
   INTEGER, PARAMETER                  :: Spn4TDzb3 =  244
   INTEGER, PARAMETER                  :: Spn5TDxb3 =  245
   INTEGER, PARAMETER                  :: Spn5TDyb3 =  246
   INTEGER, PARAMETER                  :: Spn5TDzb3 =  247
   INTEGER, PARAMETER                  :: Spn6TDxb3 =  248
   INTEGER, PARAMETER                  :: Spn6TDyb3 =  249
   INTEGER, PARAMETER                  :: Spn6TDzb3 =  250
   INTEGER, PARAMETER                  :: Spn7TDxb3 =  251
   INTEGER, PARAMETER                  :: Spn7TDyb3 =  252
   INTEGER, PARAMETER                  :: Spn7TDzb3 =  253
   INTEGER, PARAMETER                  :: Spn8TDxb3 =  254
   INTEGER, PARAMETER                  :: Spn8TDyb3 =  255
   INTEGER, PARAMETER                  :: Spn8TDzb3 =  256
   INTEGER, PARAMETER                  :: Spn9TDxb3 =  257
   INTEGER, PARAMETER                  :: Spn9TDyb3 =  258
   INTEGER, PARAMETER                  :: Spn9TDzb3 =  259
   INTEGER, PARAMETER                  :: Spn1RDxb3 =  260
   INTEGER, PARAMETER                  :: Spn1RDyb3 =  261
   INTEGER, PARAMETER                  :: Spn1RDzb3 =  262
   INTEGER, PARAMETER                  :: Spn2RDxb3 =  263
   INTEGER, PARAMETER                  :: Spn2RDyb3 =  264
   INTEGER, PARAMETER                  :: Spn2RDzb3 =  265
   INTEGER, PARAMETER                  :: Spn3RDxb3 =  266
   INTEGER, PARAMETER                  :: Spn3RDyb3 =  267
   INTEGER, PARAMETER                  :: Spn3RDzb3 =  268
   INTEGER, PARAMETER                  :: Spn4RDxb3 =  269
   INTEGER, PARAMETER                  :: Spn4RDyb3 =  270
   INTEGER, PARAMETER                  :: Spn4RDzb3 =  271
   INTEGER, PARAMETER                  :: Spn5RDxb3 =  272
   INTEGER, PARAMETER                  :: Spn5RDyb3 =  273
   INTEGER, PARAMETER                  :: Spn5RDzb3 =  274
   INTEGER, PARAMETER                  :: Spn6RDxb3 =  275
   INTEGER, PARAMETER                  :: Spn6RDyb3 =  276
   INTEGER, PARAMETER                  :: Spn6RDzb3 =  277
   INTEGER, PARAMETER                  :: Spn7RDxb3 =  278
   INTEGER, PARAMETER                  :: Spn7RDyb3 =  279
   INTEGER, PARAMETER                  :: Spn7RDzb3 =  280
   INTEGER, PARAMETER                  :: Spn8RDxb3 =  281
   INTEGER, PARAMETER                  :: Spn8RDyb3 =  282
   INTEGER, PARAMETER                  :: Spn8RDzb3 =  283
   INTEGER, PARAMETER                  :: Spn9RDxb3 =  284
   INTEGER, PARAMETER                  :: Spn9RDyb3 =  285
   INTEGER, PARAMETER                  :: Spn9RDzb3 =  286


  ! Blade Pitch Motions:

   INTEGER, PARAMETER                  :: PtchPMzc1 =  287
   INTEGER, PARAMETER                  :: PtchPMzc2 =  288
   INTEGER, PARAMETER                  :: PtchPMzc3 =  289


  ! Teeter Motions:

   INTEGER, PARAMETER                  :: TeetPya   =  290
   INTEGER, PARAMETER                  :: TeetVya   =  291
   INTEGER, PARAMETER                  :: TeetAya   =  292


  ! Shaft Motions:

   INTEGER, PARAMETER                  :: LSSTipPxa =  293
   INTEGER, PARAMETER                  :: LSSTipVxa =  294
   INTEGER, PARAMETER                  :: LSSTipAxa =  295
   INTEGER, PARAMETER                  :: LSSGagPxa =  296
   INTEGER, PARAMETER                  :: LSSGagVxa =  297
   INTEGER, PARAMETER                  :: LSSGagAxa =  298
   INTEGER, PARAMETER                  :: HSShftV   =  299
   INTEGER, PARAMETER                  :: HSShftA   =  300
   INTEGER, PARAMETER                  :: TipSpdRat =  301


  ! Nacelle IMU Motions:

   INTEGER, PARAMETER                  :: NcIMUTVxs =  302
   INTEGER, PARAMETER                  :: NcIMUTVys =  303
   INTEGER, PARAMETER                  :: NcIMUTVzs =  304
   INTEGER, PARAMETER                  :: NcIMUTAxs =  305
   INTEGER, PARAMETER                  :: NcIMUTAys =  306
   INTEGER, PARAMETER                  :: NcIMUTAzs =  307
   INTEGER, PARAMETER                  :: NcIMURVxs =  308
   INTEGER, PARAMETER                  :: NcIMURVys =  309
   INTEGER, PARAMETER                  :: NcIMURVzs =  310
   INTEGER, PARAMETER                  :: NcIMURAxs =  311
   INTEGER, PARAMETER                  :: NcIMURAys =  312
   INTEGER, PARAMETER                  :: NcIMURAzs =  313


  ! Rotor-Furl Motions:

   INTEGER, PARAMETER                  :: RotFurlP  =  314
   INTEGER, PARAMETER                  :: RotFurlV  =  315
   INTEGER, PARAMETER                  :: RotFurlA  =  316


  ! Tail-Furl Motions:

   INTEGER, PARAMETER                  :: TailFurlP =  317
   INTEGER, PARAMETER                  :: TailFurlV =  318
   INTEGER, PARAMETER                  :: TailFurlA =  319


  ! Nacelle Yaw Motions:

   INTEGER, PARAMETER                  :: YawPzn    =  320
   INTEGER, PARAMETER                  :: YawVzn    =  321
   INTEGER, PARAMETER                  :: YawAzn    =  322
   INTEGER, PARAMETER                  :: NacYawErr =  323


  ! Tower-Top / Yaw Bearing Motions:

   INTEGER, PARAMETER                  :: YawBrTDxp =  324
   INTEGER, PARAMETER                  :: YawBrTDyp =  325
   INTEGER, PARAMETER                  :: YawBrTDzp =  326
   INTEGER, PARAMETER                  :: YawBrTDxt =  327
   INTEGER, PARAMETER                  :: YawBrTDyt =  328
   INTEGER, PARAMETER                  :: YawBrTDzt =  329
   INTEGER, PARAMETER                  :: YawBrTAxp =  330
   INTEGER, PARAMETER                  :: YawBrTAyp =  331
   INTEGER, PARAMETER                  :: YawBrTAzp =  332
   INTEGER, PARAMETER                  :: YawBrRDxt =  333
   INTEGER, PARAMETER                  :: YawBrRDyt =  334
   INTEGER, PARAMETER                  :: YawBrRDzt =  335
   INTEGER, PARAMETER                  :: YawBrRVxp =  336
   INTEGER, PARAMETER                  :: YawBrRVyp =  337
   INTEGER, PARAMETER                  :: YawBrRVzp =  338
   INTEGER, PARAMETER                  :: YawBrRAxp =  339
   INTEGER, PARAMETER                  :: YawBrRAyp =  340
   INTEGER, PARAMETER                  :: YawBrRAzp =  341


  ! Local Tower Motions:

   INTEGER, PARAMETER                  :: TwHt1ALxt =  342
   INTEGER, PARAMETER                  :: TwHt1ALyt =  343
   INTEGER, PARAMETER                  :: TwHt1ALzt =  344
   INTEGER, PARAMETER                  :: TwHt2ALxt =  345
   INTEGER, PARAMETER                  :: TwHt2ALyt =  346
   INTEGER, PARAMETER                  :: TwHt2ALzt =  347
   INTEGER, PARAMETER                  :: TwHt3ALxt =  348
   INTEGER, PARAMETER                  :: TwHt3ALyt =  349
   INTEGER, PARAMETER                  :: TwHt3ALzt =  350
   INTEGER, PARAMETER                  :: TwHt4ALxt =  351
   INTEGER, PARAMETER                  :: TwHt4ALyt =  352
   INTEGER, PARAMETER                  :: TwHt4ALzt =  353
   INTEGER, PARAMETER                  :: TwHt5ALxt =  354
   INTEGER, PARAMETER                  :: TwHt5ALyt =  355
   INTEGER, PARAMETER                  :: TwHt5ALzt =  356
   INTEGER, PARAMETER                  :: TwHt6ALxt =  357
   INTEGER, PARAMETER                  :: TwHt6ALyt =  358
   INTEGER, PARAMETER                  :: TwHt6ALzt =  359
   INTEGER, PARAMETER                  :: TwHt7ALxt =  360
   INTEGER, PARAMETER                  :: TwHt7ALyt =  361
   INTEGER, PARAMETER                  :: TwHt7ALzt =  362
   INTEGER, PARAMETER                  :: TwHt8ALxt =  363
   INTEGER, PARAMETER                  :: TwHt8ALyt =  364
   INTEGER, PARAMETER                  :: TwHt8ALzt =  365
   INTEGER, PARAMETER                  :: TwHt9ALxt =  366
   INTEGER, PARAMETER                  :: TwHt9ALyt =  367
   INTEGER, PARAMETER                  :: TwHt9ALzt =  368
   INTEGER, PARAMETER                  :: TwHt1TDxt =  369
   INTEGER, PARAMETER                  :: TwHt1TDyt =  370
   INTEGER, PARAMETER                  :: TwHt1TDzt =  371
   INTEGER, PARAMETER                  :: TwHt2TDxt =  372
   INTEGER, PARAMETER                  :: TwHt2TDyt =  373
   INTEGER, PARAMETER                  :: TwHt2TDzt =  374
   INTEGER, PARAMETER                  :: TwHt3TDxt =  375
   INTEGER, PARAMETER                  :: TwHt3TDyt =  376
   INTEGER, PARAMETER                  :: TwHt3TDzt =  377
   INTEGER, PARAMETER                  :: TwHt4TDxt =  378
   INTEGER, PARAMETER                  :: TwHt4TDyt =  379
   INTEGER, PARAMETER                  :: TwHt4TDzt =  380
   INTEGER, PARAMETER                  :: TwHt5TDxt =  381
   INTEGER, PARAMETER                  :: TwHt5TDyt =  382
   INTEGER, PARAMETER                  :: TwHt5TDzt =  383
   INTEGER, PARAMETER                  :: TwHt6TDxt =  384
   INTEGER, PARAMETER                  :: TwHt6TDyt =  385
   INTEGER, PARAMETER                  :: TwHt6TDzt =  386
   INTEGER, PARAMETER                  :: TwHt7TDxt =  387
   INTEGER, PARAMETER                  :: TwHt7TDyt =  388
   INTEGER, PARAMETER                  :: TwHt7TDzt =  389
   INTEGER, PARAMETER                  :: TwHt8TDxt =  390
   INTEGER, PARAMETER                  :: TwHt8TDyt =  391
   INTEGER, PARAMETER                  :: TwHt8TDzt =  392
   INTEGER, PARAMETER                  :: TwHt9TDxt =  393
   INTEGER, PARAMETER                  :: TwHt9TDyt =  394
   INTEGER, PARAMETER                  :: TwHt9TDzt =  395
   INTEGER, PARAMETER                  :: TwHt1RDxt =  396
   INTEGER, PARAMETER                  :: TwHt1RDyt =  397
   INTEGER, PARAMETER                  :: TwHt1RDzt =  398
   INTEGER, PARAMETER                  :: TwHt2RDxt =  399
   INTEGER, PARAMETER                  :: TwHt2RDyt =  400
   INTEGER, PARAMETER                  :: TwHt2RDzt =  401
   INTEGER, PARAMETER                  :: TwHt3RDxt =  402
   INTEGER, PARAMETER                  :: TwHt3RDyt =  403
   INTEGER, PARAMETER                  :: TwHt3RDzt =  404
   INTEGER, PARAMETER                  :: TwHt4RDxt =  405
   INTEGER, PARAMETER                  :: TwHt4RDyt =  406
   INTEGER, PARAMETER                  :: TwHt4RDzt =  407
   INTEGER, PARAMETER                  :: TwHt5RDxt =  408
   INTEGER, PARAMETER                  :: TwHt5RDyt =  409
   INTEGER, PARAMETER                  :: TwHt5RDzt =  410
   INTEGER, PARAMETER                  :: TwHt6RDxt =  411
   INTEGER, PARAMETER                  :: TwHt6RDyt =  412
   INTEGER, PARAMETER                  :: TwHt6RDzt =  413
   INTEGER, PARAMETER                  :: TwHt7RDxt =  414
   INTEGER, PARAMETER                  :: TwHt7RDyt =  415
   INTEGER, PARAMETER                  :: TwHt7RDzt =  416
   INTEGER, PARAMETER                  :: TwHt8RDxt =  417
   INTEGER, PARAMETER                  :: TwHt8RDyt =  418
   INTEGER, PARAMETER                  :: TwHt8RDzt =  419
   INTEGER, PARAMETER                  :: TwHt9RDxt =  420
   INTEGER, PARAMETER                  :: TwHt9RDyt =  421
   INTEGER, PARAMETER                  :: TwHt9RDzt =  422
   INTEGER, PARAMETER                  :: TwHt1TPxi =  423
   INTEGER, PARAMETER                  :: TwHt1TPyi =  424
   INTEGER, PARAMETER                  :: TwHt1TPzi =  425
   INTEGER, PARAMETER                  :: TwHt2TPxi =  426
   INTEGER, PARAMETER                  :: TwHt2TPyi =  427
   INTEGER, PARAMETER                  :: TwHt2TPzi =  428
   INTEGER, PARAMETER                  :: TwHt3TPxi =  429
   INTEGER, PARAMETER                  :: TwHt3TPyi =  430
   INTEGER, PARAMETER                  :: TwHt3TPzi =  431
   INTEGER, PARAMETER                  :: TwHt4TPxi =  432
   INTEGER, PARAMETER                  :: TwHt4TPyi =  433
   INTEGER, PARAMETER                  :: TwHt4TPzi =  434
   INTEGER, PARAMETER                  :: TwHt5TPxi =  435
   INTEGER, PARAMETER                  :: TwHt5TPyi =  436
   INTEGER, PARAMETER                  :: TwHt5TPzi =  437
   INTEGER, PARAMETER                  :: TwHt6TPxi =  438
   INTEGER, PARAMETER                  :: TwHt6TPyi =  439
   INTEGER, PARAMETER                  :: TwHt6TPzi =  440
   INTEGER, PARAMETER                  :: TwHt7TPxi =  441
   INTEGER, PARAMETER                  :: TwHt7TPyi =  442
   INTEGER, PARAMETER                  :: TwHt7TPzi =  443
   INTEGER, PARAMETER                  :: TwHt8TPxi =  444
   INTEGER, PARAMETER                  :: TwHt8TPyi =  445
   INTEGER, PARAMETER                  :: TwHt8TPzi =  446
   INTEGER, PARAMETER                  :: TwHt9TPxi =  447
   INTEGER, PARAMETER                  :: TwHt9TPyi =  448
   INTEGER, PARAMETER                  :: TwHt9TPzi =  449
   INTEGER, PARAMETER                  :: TwHt1RPxi =  450
   INTEGER, PARAMETER                  :: TwHt1RPyi =  451
   INTEGER, PARAMETER                  :: TwHt1RPzi =  452
   INTEGER, PARAMETER                  :: TwHt2RPxi =  453
   INTEGER, PARAMETER                  :: TwHt2RPyi =  454
   INTEGER, PARAMETER                  :: TwHt2RPzi =  455
   INTEGER, PARAMETER                  :: TwHt3RPxi =  456
   INTEGER, PARAMETER                  :: TwHt3RPyi =  457
   INTEGER, PARAMETER                  :: TwHt3RPzi =  458
   INTEGER, PARAMETER                  :: TwHt4RPxi =  459
   INTEGER, PARAMETER                  :: TwHt4RPyi =  460
   INTEGER, PARAMETER                  :: TwHt4RPzi =  461
   INTEGER, PARAMETER                  :: TwHt5RPxi =  462
   INTEGER, PARAMETER                  :: TwHt5RPyi =  463
   INTEGER, PARAMETER                  :: TwHt5RPzi =  464
   INTEGER, PARAMETER                  :: TwHt6RPxi =  465
   INTEGER, PARAMETER                  :: TwHt6RPyi =  466
   INTEGER, PARAMETER                  :: TwHt6RPzi =  467
   INTEGER, PARAMETER                  :: TwHt7RPxi =  468
   INTEGER, PARAMETER                  :: TwHt7RPyi =  469
   INTEGER, PARAMETER                  :: TwHt7RPzi =  470
   INTEGER, PARAMETER                  :: TwHt8RPxi =  471
   INTEGER, PARAMETER                  :: TwHt8RPyi =  472
   INTEGER, PARAMETER                  :: TwHt8RPzi =  473
   INTEGER, PARAMETER                  :: TwHt9RPxi =  474
   INTEGER, PARAMETER                  :: TwHt9RPyi =  475
   INTEGER, PARAMETER                  :: TwHt9RPzi =  476


  ! Platform Motions:

   INTEGER, PARAMETER                  :: PtfmTDxt  =  477
   INTEGER, PARAMETER                  :: PtfmTDyt  =  478
   INTEGER, PARAMETER                  :: PtfmTDzt  =  479
   INTEGER, PARAMETER                  :: PtfmTDxi  =  480
   INTEGER, PARAMETER                  :: PtfmTDyi  =  481
   INTEGER, PARAMETER                  :: PtfmTDzi  =  482
   INTEGER, PARAMETER                  :: PtfmTVxt  =  483
   INTEGER, PARAMETER                  :: PtfmTVyt  =  484
   INTEGER, PARAMETER                  :: PtfmTVzt  =  485
   INTEGER, PARAMETER                  :: PtfmTVxi  =  486
   INTEGER, PARAMETER                  :: PtfmTVyi  =  487
   INTEGER, PARAMETER                  :: PtfmTVzi  =  488
   INTEGER, PARAMETER                  :: PtfmTAxt  =  489
   INTEGER, PARAMETER                  :: PtfmTAyt  =  490
   INTEGER, PARAMETER                  :: PtfmTAzt  =  491
   INTEGER, PARAMETER                  :: PtfmTAxi  =  492
   INTEGER, PARAMETER                  :: PtfmTAyi  =  493
   INTEGER, PARAMETER                  :: PtfmTAzi  =  494
   INTEGER, PARAMETER                  :: PtfmRDxi  =  495
   INTEGER, PARAMETER                  :: PtfmRDyi  =  496
   INTEGER, PARAMETER                  :: PtfmRDzi  =  497
   INTEGER, PARAMETER                  :: PtfmRVxt  =  498
   INTEGER, PARAMETER                  :: PtfmRVyt  =  499
   INTEGER, PARAMETER                  :: PtfmRVzt  =  500
   INTEGER, PARAMETER                  :: PtfmRVxi  =  501
   INTEGER, PARAMETER                  :: PtfmRVyi  =  502
   INTEGER, PARAMETER                  :: PtfmRVzi  =  503
   INTEGER, PARAMETER                  :: PtfmRAxt  =  504
   INTEGER, PARAMETER                  :: PtfmRAyt  =  505
   INTEGER, PARAMETER                  :: PtfmRAzt  =  506
   INTEGER, PARAMETER                  :: PtfmRAxi  =  507
   INTEGER, PARAMETER                  :: PtfmRAyi  =  508
   INTEGER, PARAMETER                  :: PtfmRAzi  =  509


  ! Blade 1 Root Loads:

   INTEGER, PARAMETER                  :: RootFxc1  =  510
   INTEGER, PARAMETER                  :: RootFyc1  =  511
   INTEGER, PARAMETER                  :: RootFzc1  =  512
   INTEGER, PARAMETER                  :: RootFxb1  =  513
   INTEGER, PARAMETER                  :: RootFyb1  =  514
   INTEGER, PARAMETER                  :: RootMxc1  =  515
   INTEGER, PARAMETER                  :: RootMyc1  =  516
   INTEGER, PARAMETER                  :: RootMzc1  =  517
   INTEGER, PARAMETER                  :: RootMxb1  =  518
   INTEGER, PARAMETER                  :: RootMyb1  =  519


  ! Blade 2 Root Loads:

   INTEGER, PARAMETER                  :: RootFxc2  =  520
   INTEGER, PARAMETER                  :: RootFyc2  =  521
   INTEGER, PARAMETER                  :: RootFzc2  =  522
   INTEGER, PARAMETER                  :: RootFxb2  =  523
   INTEGER, PARAMETER                  :: RootFyb2  =  524
   INTEGER, PARAMETER                  :: RootMxc2  =  525
   INTEGER, PARAMETER                  :: RootMyc2  =  526
   INTEGER, PARAMETER                  :: RootMzc2  =  527
   INTEGER, PARAMETER                  :: RootMxb2  =  528
   INTEGER, PARAMETER                  :: RootMyb2  =  529


  ! Blade 3 Root Loads:

   INTEGER, PARAMETER                  :: RootFxc3  =  530
   INTEGER, PARAMETER                  :: RootFyc3  =  531
   INTEGER, PARAMETER                  :: RootFzc3  =  532
   INTEGER, PARAMETER                  :: RootFxb3  =  533
   INTEGER, PARAMETER                  :: RootFyb3  =  534
   INTEGER, PARAMETER                  :: RootMxc3  =  535
   INTEGER, PARAMETER                  :: RootMyc3  =  536
   INTEGER, PARAMETER                  :: RootMzc3  =  537
   INTEGER, PARAMETER                  :: RootMxb3  =  538
   INTEGER, PARAMETER                  :: RootMyb3  =  539


  ! Blade 1 Local Span Loads:

   INTEGER, PARAMETER                  :: Spn1MLxb1 =  540
   INTEGER, PARAMETER                  :: Spn1MLyb1 =  541
   INTEGER, PARAMETER                  :: Spn1MLzb1 =  542
   INTEGER, PARAMETER                  :: Spn2MLxb1 =  543
   INTEGER, PARAMETER                  :: Spn2MLyb1 =  544
   INTEGER, PARAMETER                  :: Spn2MLzb1 =  545
   INTEGER, PARAMETER                  :: Spn3MLxb1 =  546
   INTEGER, PARAMETER                  :: Spn3MLyb1 =  547
   INTEGER, PARAMETER                  :: Spn3MLzb1 =  548
   INTEGER, PARAMETER                  :: Spn4MLxb1 =  549
   INTEGER, PARAMETER                  :: Spn4MLyb1 =  550
   INTEGER, PARAMETER                  :: Spn4MLzb1 =  551
   INTEGER, PARAMETER                  :: Spn5MLxb1 =  552
   INTEGER, PARAMETER                  :: Spn5MLyb1 =  553
   INTEGER, PARAMETER                  :: Spn5MLzb1 =  554
   INTEGER, PARAMETER                  :: Spn6MLxb1 =  555
   INTEGER, PARAMETER                  :: Spn6MLyb1 =  556
   INTEGER, PARAMETER                  :: Spn6MLzb1 =  557
   INTEGER, PARAMETER                  :: Spn7MLxb1 =  558
   INTEGER, PARAMETER                  :: Spn7MLyb1 =  559
   INTEGER, PARAMETER                  :: Spn7MLzb1 =  560
   INTEGER, PARAMETER                  :: Spn8MLxb1 =  561
   INTEGER, PARAMETER                  :: Spn8MLyb1 =  562
   INTEGER, PARAMETER                  :: Spn8MLzb1 =  563
   INTEGER, PARAMETER                  :: Spn9MLxb1 =  564
   INTEGER, PARAMETER                  :: Spn9MLyb1 =  565
   INTEGER, PARAMETER                  :: Spn9MLzb1 =  566
   INTEGER, PARAMETER                  :: Spn1FLxb1 =  567
   INTEGER, PARAMETER                  :: Spn1FLyb1 =  568
   INTEGER, PARAMETER                  :: Spn1FLzb1 =  569
   INTEGER, PARAMETER                  :: Spn2FLxb1 =  570
   INTEGER, PARAMETER                  :: Spn2FLyb1 =  571
   INTEGER, PARAMETER                  :: Spn2FLzb1 =  572
   INTEGER, PARAMETER                  :: Spn3FLxb1 =  573
   INTEGER, PARAMETER                  :: Spn3FLyb1 =  574
   INTEGER, PARAMETER                  :: Spn3FLzb1 =  575
   INTEGER, PARAMETER                  :: Spn4FLxb1 =  576
   INTEGER, PARAMETER                  :: Spn4FLyb1 =  577
   INTEGER, PARAMETER                  :: Spn4FLzb1 =  578
   INTEGER, PARAMETER                  :: Spn5FLxb1 =  579
   INTEGER, PARAMETER                  :: Spn5FLyb1 =  580
   INTEGER, PARAMETER                  :: Spn5FLzb1 =  581
   INTEGER, PARAMETER                  :: Spn6FLxb1 =  582
   INTEGER, PARAMETER                  :: Spn6FLyb1 =  583
   INTEGER, PARAMETER                  :: Spn6FLzb1 =  584
   INTEGER, PARAMETER                  :: Spn7FLxb1 =  585
   INTEGER, PARAMETER                  :: Spn7FLyb1 =  586
   INTEGER, PARAMETER                  :: Spn7FLzb1 =  587
   INTEGER, PARAMETER                  :: Spn8FLxb1 =  588
   INTEGER, PARAMETER                  :: Spn8FLyb1 =  589
   INTEGER, PARAMETER                  :: Spn8FLzb1 =  590
   INTEGER, PARAMETER                  :: Spn9FLxb1 =  591
   INTEGER, PARAMETER                  :: Spn9FLyb1 =  592
   INTEGER, PARAMETER                  :: Spn9FLzb1 =  593


  ! Blade 2 Local Span Loads:

   INTEGER, PARAMETER                  :: Spn1MLxb2 =  594
   INTEGER, PARAMETER                  :: Spn1MLyb2 =  595
   INTEGER, PARAMETER                  :: Spn1MLzb2 =  596
   INTEGER, PARAMETER                  :: Spn2MLxb2 =  597
   INTEGER, PARAMETER                  :: Spn2MLyb2 =  598
   INTEGER, PARAMETER                  :: Spn2MLzb2 =  599
   INTEGER, PARAMETER                  :: Spn3MLxb2 =  600
   INTEGER, PARAMETER                  :: Spn3MLyb2 =  601
   INTEGER, PARAMETER                  :: Spn3MLzb2 =  602
   INTEGER, PARAMETER                  :: Spn4MLxb2 =  603
   INTEGER, PARAMETER                  :: Spn4MLyb2 =  604
   INTEGER, PARAMETER                  :: Spn4MLzb2 =  605
   INTEGER, PARAMETER                  :: Spn5MLxb2 =  606
   INTEGER, PARAMETER                  :: Spn5MLyb2 =  607
   INTEGER, PARAMETER                  :: Spn5MLzb2 =  608
   INTEGER, PARAMETER                  :: Spn6MLxb2 =  609
   INTEGER, PARAMETER                  :: Spn6MLyb2 =  610
   INTEGER, PARAMETER                  :: Spn6MLzb2 =  611
   INTEGER, PARAMETER                  :: Spn7MLxb2 =  612
   INTEGER, PARAMETER                  :: Spn7MLyb2 =  613
   INTEGER, PARAMETER                  :: Spn7MLzb2 =  614
   INTEGER, PARAMETER                  :: Spn8MLxb2 =  615
   INTEGER, PARAMETER                  :: Spn8MLyb2 =  616
   INTEGER, PARAMETER                  :: Spn8MLzb2 =  617
   INTEGER, PARAMETER                  :: Spn9MLxb2 =  618
   INTEGER, PARAMETER                  :: Spn9MLyb2 =  619
   INTEGER, PARAMETER                  :: Spn9MLzb2 =  620
   INTEGER, PARAMETER                  :: Spn1FLxb2 =  621
   INTEGER, PARAMETER                  :: Spn1FLyb2 =  622
   INTEGER, PARAMETER                  :: Spn1FLzb2 =  623
   INTEGER, PARAMETER                  :: Spn2FLxb2 =  624
   INTEGER, PARAMETER                  :: Spn2FLyb2 =  625
   INTEGER, PARAMETER                  :: Spn2FLzb2 =  626
   INTEGER, PARAMETER                  :: Spn3FLxb2 =  627
   INTEGER, PARAMETER                  :: Spn3FLyb2 =  628
   INTEGER, PARAMETER                  :: Spn3FLzb2 =  629
   INTEGER, PARAMETER                  :: Spn4FLxb2 =  630
   INTEGER, PARAMETER                  :: Spn4FLyb2 =  631
   INTEGER, PARAMETER                  :: Spn4FLzb2 =  632
   INTEGER, PARAMETER                  :: Spn5FLxb2 =  633
   INTEGER, PARAMETER                  :: Spn5FLyb2 =  634
   INTEGER, PARAMETER                  :: Spn5FLzb2 =  635
   INTEGER, PARAMETER                  :: Spn6FLxb2 =  636
   INTEGER, PARAMETER                  :: Spn6FLyb2 =  637
   INTEGER, PARAMETER                  :: Spn6FLzb2 =  638
   INTEGER, PARAMETER                  :: Spn7FLxb2 =  639
   INTEGER, PARAMETER                  :: Spn7FLyb2 =  640
   INTEGER, PARAMETER                  :: Spn7FLzb2 =  641
   INTEGER, PARAMETER                  :: Spn8FLxb2 =  642
   INTEGER, PARAMETER                  :: Spn8FLyb2 =  643
   INTEGER, PARAMETER                  :: Spn8FLzb2 =  644
   INTEGER, PARAMETER                  :: Spn9FLxb2 =  645
   INTEGER, PARAMETER                  :: Spn9FLyb2 =  646
   INTEGER, PARAMETER                  :: Spn9FLzb2 =  647


  ! Blade 3 Local Span Loads:

   INTEGER, PARAMETER                  :: Spn1MLxb3 =  648
   INTEGER, PARAMETER                  :: Spn1MLyb3 =  649
   INTEGER, PARAMETER                  :: Spn1MLzb3 =  650
   INTEGER, PARAMETER                  :: Spn2MLxb3 =  651
   INTEGER, PARAMETER                  :: Spn2MLyb3 =  652
   INTEGER, PARAMETER                  :: Spn2MLzb3 =  653
   INTEGER, PARAMETER                  :: Spn3MLxb3 =  654
   INTEGER, PARAMETER                  :: Spn3MLyb3 =  655
   INTEGER, PARAMETER                  :: Spn3MLzb3 =  656
   INTEGER, PARAMETER                  :: Spn4MLxb3 =  657
   INTEGER, PARAMETER                  :: Spn4MLyb3 =  658
   INTEGER, PARAMETER                  :: Spn4MLzb3 =  659
   INTEGER, PARAMETER                  :: Spn5MLxb3 =  660
   INTEGER, PARAMETER                  :: Spn5MLyb3 =  661
   INTEGER, PARAMETER                  :: Spn5MLzb3 =  662
   INTEGER, PARAMETER                  :: Spn6MLxb3 =  663
   INTEGER, PARAMETER                  :: Spn6MLyb3 =  664
   INTEGER, PARAMETER                  :: Spn6MLzb3 =  665
   INTEGER, PARAMETER                  :: Spn7MLxb3 =  666
   INTEGER, PARAMETER                  :: Spn7MLyb3 =  667
   INTEGER, PARAMETER                  :: Spn7MLzb3 =  668
   INTEGER, PARAMETER                  :: Spn8MLxb3 =  669
   INTEGER, PARAMETER                  :: Spn8MLyb3 =  670
   INTEGER, PARAMETER                  :: Spn8MLzb3 =  671
   INTEGER, PARAMETER                  :: Spn9MLxb3 =  672
   INTEGER, PARAMETER                  :: Spn9MLyb3 =  673
   INTEGER, PARAMETER                  :: Spn9MLzb3 =  674
   INTEGER, PARAMETER                  :: Spn1FLxb3 =  675
   INTEGER, PARAMETER                  :: Spn1FLyb3 =  676
   INTEGER, PARAMETER                  :: Spn1FLzb3 =  677
   INTEGER, PARAMETER                  :: Spn2FLxb3 =  678
   INTEGER, PARAMETER                  :: Spn2FLyb3 =  679
   INTEGER, PARAMETER                  :: Spn2FLzb3 =  680
   INTEGER, PARAMETER                  :: Spn3FLxb3 =  681
   INTEGER, PARAMETER                  :: Spn3FLyb3 =  682
   INTEGER, PARAMETER                  :: Spn3FLzb3 =  683
   INTEGER, PARAMETER                  :: Spn4FLxb3 =  684
   INTEGER, PARAMETER                  :: Spn4FLyb3 =  685
   INTEGER, PARAMETER                  :: Spn4FLzb3 =  686
   INTEGER, PARAMETER                  :: Spn5FLxb3 =  687
   INTEGER, PARAMETER                  :: Spn5FLyb3 =  688
   INTEGER, PARAMETER                  :: Spn5FLzb3 =  689
   INTEGER, PARAMETER                  :: Spn6FLxb3 =  690
   INTEGER, PARAMETER                  :: Spn6FLyb3 =  691
   INTEGER, PARAMETER                  :: Spn6FLzb3 =  692
   INTEGER, PARAMETER                  :: Spn7FLxb3 =  693
   INTEGER, PARAMETER                  :: Spn7FLyb3 =  694
   INTEGER, PARAMETER                  :: Spn7FLzb3 =  695
   INTEGER, PARAMETER                  :: Spn8FLxb3 =  696
   INTEGER, PARAMETER                  :: Spn8FLyb3 =  697
   INTEGER, PARAMETER                  :: Spn8FLzb3 =  698
   INTEGER, PARAMETER                  :: Spn9FLxb3 =  699
   INTEGER, PARAMETER                  :: Spn9FLyb3 =  700
   INTEGER, PARAMETER                  :: Spn9FLzb3 =  701


  ! Hub and Rotor Loads:

   INTEGER, PARAMETER                  :: LSShftFxa =  702
   INTEGER, PARAMETER                  :: LSShftFya =  703
   INTEGER, PARAMETER                  :: LSShftFza =  704
   INTEGER, PARAMETER                  :: LSShftFys =  705
   INTEGER, PARAMETER                  :: LSShftFzs =  706
   INTEGER, PARAMETER                  :: LSShftMxa =  707
   INTEGER, PARAMETER                  :: LSSTipMya =  708
   INTEGER, PARAMETER                  :: LSSTipMza =  709
   INTEGER, PARAMETER                  :: LSSTipMys =  710
   INTEGER, PARAMETER                  :: LSSTipMzs =  711
   INTEGER, PARAMETER                  :: CThrstAzm =  712
   INTEGER, PARAMETER                  :: CThrstRad =  713
   INTEGER, PARAMETER                  :: RotPwr    =  714
   INTEGER, PARAMETER                  :: RotCq     =  715
   INTEGER, PARAMETER                  :: RotCp     =  716
   INTEGER, PARAMETER                  :: RotCt     =  717


  ! Shaft Strain Gage Loads:

   INTEGER, PARAMETER                  :: LSSGagMya =  718
   INTEGER, PARAMETER                  :: LSSGagMza =  719
   INTEGER, PARAMETER                  :: LSSGagMys =  720
   INTEGER, PARAMETER                  :: LSSGagMzs =  721


  ! Generator and High-Speed Shaft Loads:

   INTEGER, PARAMETER                  :: HSShftTq  =  722
   INTEGER, PARAMETER                  :: HSShftPwr =  723
   INTEGER, PARAMETER                  :: HSShftCq  =  724
   INTEGER, PARAMETER                  :: HSShftCp  =  725
   INTEGER, PARAMETER                  :: GenTq     =  726
   INTEGER, PARAMETER                  :: GenPwr    =  727
   INTEGER, PARAMETER                  :: GenCq     =  728
   INTEGER, PARAMETER                  :: GenCp     =  729
   INTEGER, PARAMETER                  :: HSSBrTq   =  730


  ! Rotor-Furl Bearing Loads:

   INTEGER, PARAMETER                  :: RFrlBrM   =  731


  ! Tail-Furl Bearing Loads:

   INTEGER, PARAMETER                  :: TFrlBrM   =  732


  ! Tail Fin Aerodynamic Loads:

   INTEGER, PARAMETER                  :: TFinAlpha =  733
   INTEGER, PARAMETER                  :: TFinCLift =  734
   INTEGER, PARAMETER                  :: TFinCDrag =  735
   INTEGER, PARAMETER                  :: TFinDnPrs =  736
   INTEGER, PARAMETER                  :: TFinCPFx  =  737
   INTEGER, PARAMETER                  :: TFinCPFy  =  738


  ! Tower-Top / Yaw Bearing Loads:

   INTEGER, PARAMETER                  :: YawBrFxn  =  739
   INTEGER, PARAMETER                  :: YawBrFyn  =  740
   INTEGER, PARAMETER                  :: YawBrFzn  =  741
   INTEGER, PARAMETER                  :: YawBrFxp  =  742
   INTEGER, PARAMETER                  :: YawBrFyp  =  743
   INTEGER, PARAMETER                  :: YawBrMxn  =  744
   INTEGER, PARAMETER                  :: YawBrMyn  =  745
   INTEGER, PARAMETER                  :: YawBrMzn  =  746
   INTEGER, PARAMETER                  :: YawBrMxp  =  747
   INTEGER, PARAMETER                  :: YawBrMyp  =  748


  ! Tower Base Loads:

   INTEGER, PARAMETER                  :: TwrBsFxt  =  749
   INTEGER, PARAMETER                  :: TwrBsFyt  =  750
   INTEGER, PARAMETER                  :: TwrBsFzt  =  751
   INTEGER, PARAMETER                  :: TwrBsMxt  =  752
   INTEGER, PARAMETER                  :: TwrBsMyt  =  753
   INTEGER, PARAMETER                  :: TwrBsMzt  =  754


  ! Local Tower Loads:

   INTEGER, PARAMETER                  :: TwHt1MLxt =  755
   INTEGER, PARAMETER                  :: TwHt1MLyt =  756
   INTEGER, PARAMETER                  :: TwHt1MLzt =  757
   INTEGER, PARAMETER                  :: TwHt2MLxt =  758
   INTEGER, PARAMETER                  :: TwHt2MLyt =  759
   INTEGER, PARAMETER                  :: TwHt2MLzt =  760
   INTEGER, PARAMETER                  :: TwHt3MLxt =  761
   INTEGER, PARAMETER                  :: TwHt3MLyt =  762
   INTEGER, PARAMETER                  :: TwHt3MLzt =  763
   INTEGER, PARAMETER                  :: TwHt4MLxt =  764
   INTEGER, PARAMETER                  :: TwHt4MLyt =  765
   INTEGER, PARAMETER                  :: TwHt4MLzt =  766
   INTEGER, PARAMETER                  :: TwHt5MLxt =  767
   INTEGER, PARAMETER                  :: TwHt5MLyt =  768
   INTEGER, PARAMETER                  :: TwHt5MLzt =  769
   INTEGER, PARAMETER                  :: TwHt6MLxt =  770
   INTEGER, PARAMETER                  :: TwHt6MLyt =  771
   INTEGER, PARAMETER                  :: TwHt6MLzt =  772
   INTEGER, PARAMETER                  :: TwHt7MLxt =  773
   INTEGER, PARAMETER                  :: TwHt7MLyt =  774
   INTEGER, PARAMETER                  :: TwHt7MLzt =  775
   INTEGER, PARAMETER                  :: TwHt8MLxt =  776
   INTEGER, PARAMETER                  :: TwHt8MLyt =  777
   INTEGER, PARAMETER                  :: TwHt8MLzt =  778
   INTEGER, PARAMETER                  :: TwHt9MLxt =  779
   INTEGER, PARAMETER                  :: TwHt9MLyt =  780
   INTEGER, PARAMETER                  :: TwHt9MLzt =  781
   INTEGER, PARAMETER                  :: TwHt1FLxt =  782
   INTEGER, PARAMETER                  :: TwHt1FLyt =  783
   INTEGER, PARAMETER                  :: TwHt1FLzt =  784
   INTEGER, PARAMETER                  :: TwHt2FLxt =  785
   INTEGER, PARAMETER                  :: TwHt2FLyt =  786
   INTEGER, PARAMETER                  :: TwHt2FLzt =  787
   INTEGER, PARAMETER                  :: TwHt3FLxt =  788
   INTEGER, PARAMETER                  :: TwHt3FLyt =  789
   INTEGER, PARAMETER                  :: TwHt3FLzt =  790
   INTEGER, PARAMETER                  :: TwHt4FLxt =  791
   INTEGER, PARAMETER                  :: TwHt4FLyt =  792
   INTEGER, PARAMETER                  :: TwHt4FLzt =  793
   INTEGER, PARAMETER                  :: TwHt5FLxt =  794
   INTEGER, PARAMETER                  :: TwHt5FLyt =  795
   INTEGER, PARAMETER                  :: TwHt5FLzt =  796
   INTEGER, PARAMETER                  :: TwHt6FLxt =  797
   INTEGER, PARAMETER                  :: TwHt6FLyt =  798
   INTEGER, PARAMETER                  :: TwHt6FLzt =  799
   INTEGER, PARAMETER                  :: TwHt7FLxt =  800
   INTEGER, PARAMETER                  :: TwHt7FLyt =  801
   INTEGER, PARAMETER                  :: TwHt7FLzt =  802
   INTEGER, PARAMETER                  :: TwHt8FLxt =  803
   INTEGER, PARAMETER                  :: TwHt8FLyt =  804
   INTEGER, PARAMETER                  :: TwHt8FLzt =  805
   INTEGER, PARAMETER                  :: TwHt9FLxt =  806
   INTEGER, PARAMETER                  :: TwHt9FLyt =  807
   INTEGER, PARAMETER                  :: TwHt9FLzt =  808


  ! Platform Loads:

   INTEGER, PARAMETER                  :: PtfmFxt   =  809
   INTEGER, PARAMETER                  :: PtfmFyt   =  810
   INTEGER, PARAMETER                  :: PtfmFzt   =  811
   INTEGER, PARAMETER                  :: PtfmFxi   =  812
   INTEGER, PARAMETER                  :: PtfmFyi   =  813
   INTEGER, PARAMETER                  :: PtfmFzi   =  814
   INTEGER, PARAMETER                  :: PtfmMxt   =  815
   INTEGER, PARAMETER                  :: PtfmMyt   =  816
   INTEGER, PARAMETER                  :: PtfmMzt   =  817
   INTEGER, PARAMETER                  :: PtfmMxi   =  818
   INTEGER, PARAMETER                  :: PtfmMyi   =  819
   INTEGER, PARAMETER                  :: PtfmMzi   =  820


  ! Mooring Line Loads:

   INTEGER, PARAMETER                  :: Fair1Ten  =  821
   INTEGER, PARAMETER                  :: Fair1Ang  =  822
   INTEGER, PARAMETER                  :: Anch1Ten  =  823
   INTEGER, PARAMETER                  :: Anch1Ang  =  824
   INTEGER, PARAMETER                  :: Fair2Ten  =  825
   INTEGER, PARAMETER                  :: Fair2Ang  =  826
   INTEGER, PARAMETER                  :: Anch2Ten  =  827
   INTEGER, PARAMETER                  :: Anch2Ang  =  828
   INTEGER, PARAMETER                  :: Fair3Ten  =  829
   INTEGER, PARAMETER                  :: Fair3Ang  =  830
   INTEGER, PARAMETER                  :: Anch3Ten  =  831
   INTEGER, PARAMETER                  :: Anch3Ang  =  832
   INTEGER, PARAMETER                  :: Fair4Ten  =  833
   INTEGER, PARAMETER                  :: Fair4Ang  =  834
   INTEGER, PARAMETER                  :: Anch4Ten  =  835
   INTEGER, PARAMETER                  :: Anch4Ang  =  836
   INTEGER, PARAMETER                  :: Fair5Ten  =  837
   INTEGER, PARAMETER                  :: Fair5Ang  =  838
   INTEGER, PARAMETER                  :: Anch5Ten  =  839
   INTEGER, PARAMETER                  :: Anch5Ang  =  840
   INTEGER, PARAMETER                  :: Fair6Ten  =  841
   INTEGER, PARAMETER                  :: Fair6Ang  =  842
   INTEGER, PARAMETER                  :: Anch6Ten  =  843
   INTEGER, PARAMETER                  :: Anch6Ang  =  844
   INTEGER, PARAMETER                  :: Fair7Ten  =  845
   INTEGER, PARAMETER                  :: Fair7Ang  =  846
   INTEGER, PARAMETER                  :: Anch7Ten  =  847
   INTEGER, PARAMETER                  :: Anch7Ang  =  848
   INTEGER, PARAMETER                  :: Fair8Ten  =  849
   INTEGER, PARAMETER                  :: Fair8Ang  =  850
   INTEGER, PARAMETER                  :: Anch8Ten  =  851
   INTEGER, PARAMETER                  :: Anch8Ang  =  852
   INTEGER, PARAMETER                  :: Fair9Ten  =  853
   INTEGER, PARAMETER                  :: Fair9Ang  =  854
   INTEGER, PARAMETER                  :: Anch9Ten  =  855
   INTEGER, PARAMETER                  :: Anch9Ang  =  856


  ! Wave Motions:

   INTEGER, PARAMETER                  :: WaveElev  =  857
   INTEGER, PARAMETER                  :: Wave1Vxi  =  858
   INTEGER, PARAMETER                  :: Wave1Vyi  =  859
   INTEGER, PARAMETER                  :: Wave1Vzi  =  860
   INTEGER, PARAMETER                  :: Wave1Axi  =  861
   INTEGER, PARAMETER                  :: Wave1Ayi  =  862
   INTEGER, PARAMETER                  :: Wave1Azi  =  863
   INTEGER, PARAMETER                  :: Wave2Vxi  =  864
   INTEGER, PARAMETER                  :: Wave2Vyi  =  865
   INTEGER, PARAMETER                  :: Wave2Vzi  =  866
   INTEGER, PARAMETER                  :: Wave2Axi  =  867
   INTEGER, PARAMETER                  :: Wave2Ayi  =  868
   INTEGER, PARAMETER                  :: Wave2Azi  =  869
   INTEGER, PARAMETER                  :: Wave3Vxi  =  870
   INTEGER, PARAMETER                  :: Wave3Vyi  =  871
   INTEGER, PARAMETER                  :: Wave3Vzi  =  872
   INTEGER, PARAMETER                  :: Wave3Axi  =  873
   INTEGER, PARAMETER                  :: Wave3Ayi  =  874
   INTEGER, PARAMETER                  :: Wave3Azi  =  875
   INTEGER, PARAMETER                  :: Wave4Vxi  =  876
   INTEGER, PARAMETER                  :: Wave4Vyi  =  877
   INTEGER, PARAMETER                  :: Wave4Vzi  =  878
   INTEGER, PARAMETER                  :: Wave4Axi  =  879
   INTEGER, PARAMETER                  :: Wave4Ayi  =  880
   INTEGER, PARAMETER                  :: Wave4Azi  =  881
   INTEGER, PARAMETER                  :: Wave5Vxi  =  882
   INTEGER, PARAMETER                  :: Wave5Vyi  =  883
   INTEGER, PARAMETER                  :: Wave5Vzi  =  884
   INTEGER, PARAMETER                  :: Wave5Axi  =  885
   INTEGER, PARAMETER                  :: Wave5Ayi  =  886
   INTEGER, PARAMETER                  :: Wave5Azi  =  887
   INTEGER, PARAMETER                  :: Wave6Vxi  =  888
   INTEGER, PARAMETER                  :: Wave6Vyi  =  889
   INTEGER, PARAMETER                  :: Wave6Vzi  =  890
   INTEGER, PARAMETER                  :: Wave6Axi  =  891
   INTEGER, PARAMETER                  :: Wave6Ayi  =  892
   INTEGER, PARAMETER                  :: Wave6Azi  =  893
   INTEGER, PARAMETER                  :: Wave7Vxi  =  894
   INTEGER, PARAMETER                  :: Wave7Vyi  =  895
   INTEGER, PARAMETER                  :: Wave7Vzi  =  896
   INTEGER, PARAMETER                  :: Wave7Axi  =  897
   INTEGER, PARAMETER                  :: Wave7Ayi  =  898
   INTEGER, PARAMETER                  :: Wave7Azi  =  899
   INTEGER, PARAMETER                  :: Wave8Vxi  =  900
   INTEGER, PARAMETER                  :: Wave8Vyi  =  901
   INTEGER, PARAMETER                  :: Wave8Vzi  =  902
   INTEGER, PARAMETER                  :: Wave8Axi  =  903
   INTEGER, PARAMETER                  :: Wave8Ayi  =  904
   INTEGER, PARAMETER                  :: Wave8Azi  =  905
   INTEGER, PARAMETER                  :: Wave9Vxi  =  906
   INTEGER, PARAMETER                  :: Wave9Vyi  =  907
   INTEGER, PARAMETER                  :: Wave9Vzi  =  908
   INTEGER, PARAMETER                  :: Wave9Axi  =  909
   INTEGER, PARAMETER                  :: Wave9Ayi  =  910
   INTEGER, PARAMETER                  :: Wave9Azi  =  911


  ! Internal Degrees of Freedom:

   INTEGER, PARAMETER                  :: Q_B1E1    =  912
   INTEGER, PARAMETER                  :: Q_B2E1    =  913
   INTEGER, PARAMETER                  :: Q_B3E1    =  914
   INTEGER, PARAMETER                  :: Q_B1F1    =  915
   INTEGER, PARAMETER                  :: Q_B2F1    =  916
   INTEGER, PARAMETER                  :: Q_B3F1    =  917
   INTEGER, PARAMETER                  :: Q_B1F2    =  918
   INTEGER, PARAMETER                  :: Q_B2F2    =  919
   INTEGER, PARAMETER                  :: Q_B3F2    =  920
   INTEGER, PARAMETER                  :: Q_Teet    =  921
   INTEGER, PARAMETER                  :: Q_DrTr    =  922
   INTEGER, PARAMETER                  :: Q_GeAz    =  923
   INTEGER, PARAMETER                  :: Q_RFrl    =  924
   INTEGER, PARAMETER                  :: Q_TFrl    =  925
   INTEGER, PARAMETER                  :: Q_Yaw     =  926
   INTEGER, PARAMETER                  :: Q_TFA1    =  927
   INTEGER, PARAMETER                  :: Q_TSS1    =  928
   INTEGER, PARAMETER                  :: Q_TFA2    =  929
   INTEGER, PARAMETER                  :: Q_TSS2    =  930
   INTEGER, PARAMETER                  :: Q_Sg      =  931
   INTEGER, PARAMETER                  :: Q_Sw      =  932
   INTEGER, PARAMETER                  :: Q_Hv      =  933
   INTEGER, PARAMETER                  :: Q_R       =  934
   INTEGER, PARAMETER                  :: Q_P       =  935
   INTEGER, PARAMETER                  :: Q_Y       =  936
   INTEGER, PARAMETER                  :: QD_B1E1   =  937
   INTEGER, PARAMETER                  :: QD_B2E1   =  938
   INTEGER, PARAMETER                  :: QD_B3E1   =  939
   INTEGER, PARAMETER                  :: QD_B1F1   =  940
   INTEGER, PARAMETER                  :: QD_B2F1   =  941
   INTEGER, PARAMETER                  :: QD_B3F1   =  942
   INTEGER, PARAMETER                  :: QD_B1F2   =  943
   INTEGER, PARAMETER                  :: QD_B2F2   =  944
   INTEGER, PARAMETER                  :: QD_B3F2   =  945
   INTEGER, PARAMETER                  :: QD_Teet   =  946
   INTEGER, PARAMETER                  :: QD_DrTr   =  947
   INTEGER, PARAMETER                  :: QD_GeAz   =  948
   INTEGER, PARAMETER                  :: QD_RFrl   =  949
   INTEGER, PARAMETER                  :: QD_TFrl   =  950
   INTEGER, PARAMETER                  :: QD_Yaw    =  951
   INTEGER, PARAMETER                  :: QD_TFA1   =  952
   INTEGER, PARAMETER                  :: QD_TSS1   =  953
   INTEGER, PARAMETER                  :: QD_TFA2   =  954
   INTEGER, PARAMETER                  :: QD_TSS2   =  955
   INTEGER, PARAMETER                  :: QD_Sg     =  956
   INTEGER, PARAMETER                  :: QD_Sw     =  957
   INTEGER, PARAMETER                  :: QD_Hv     =  958
   INTEGER, PARAMETER                  :: QD_R      =  959
   INTEGER, PARAMETER                  :: QD_P      =  960
   INTEGER, PARAMETER                  :: QD_Y      =  961
   INTEGER, PARAMETER                  :: QD2_B1E1  =  962
   INTEGER, PARAMETER                  :: QD2_B2E1  =  963
   INTEGER, PARAMETER                  :: QD2_B3E1  =  964
   INTEGER, PARAMETER                  :: QD2_B1F1  =  965
   INTEGER, PARAMETER                  :: QD2_B2F1  =  966
   INTEGER, PARAMETER                  :: QD2_B3F1  =  967
   INTEGER, PARAMETER                  :: QD2_B1F2  =  968
   INTEGER, PARAMETER                  :: QD2_B2F2  =  969
   INTEGER, PARAMETER                  :: QD2_B3F2  =  970
   INTEGER, PARAMETER                  :: QD2_Teet  =  971
   INTEGER, PARAMETER                  :: QD2_DrTr  =  972
   INTEGER, PARAMETER                  :: QD2_GeAz  =  973
   INTEGER, PARAMETER                  :: QD2_RFrl  =  974
   INTEGER, PARAMETER                  :: QD2_TFrl  =  975
   INTEGER, PARAMETER                  :: QD2_Yaw   =  976
   INTEGER, PARAMETER                  :: QD2_TFA1  =  977
   INTEGER, PARAMETER                  :: QD2_TSS1  =  978
   INTEGER, PARAMETER                  :: QD2_TFA2  =  979
   INTEGER, PARAMETER                  :: QD2_TSS2  =  980
   INTEGER, PARAMETER                  :: QD2_Sg    =  981
   INTEGER, PARAMETER                  :: QD2_Sw    =  982
   INTEGER, PARAMETER                  :: QD2_Hv    =  983
   INTEGER, PARAMETER                  :: QD2_R     =  984
   INTEGER, PARAMETER                  :: QD2_P     =  985
   INTEGER, PARAMETER                  :: QD2_Y     =  986

   ! Yang: Start of adding the output channels of TMDs
      ! TMDs:

INTEGER(4), PARAMETER                  :: TmdXDxn              =	987
INTEGER(4), PARAMETER                  :: TmdXDyn              =	988
INTEGER(4), PARAMETER                  :: TmdXDzn              =	989
INTEGER(4), PARAMETER                  :: TmdXVxn              =	990
INTEGER(4), PARAMETER                  :: TmdXVyn              =	991
INTEGER(4), PARAMETER                  :: TmdXVzn              =	992
INTEGER(4), PARAMETER                  :: TmdXAxn              =	993
INTEGER(4), PARAMETER                  :: TmdXAyn              =	994
INTEGER(4), PARAMETER                  :: TmdXAzn              =	995
INTEGER(4), PARAMETER                  :: TmdYDxn              =	996
INTEGER(4), PARAMETER                  :: TmdYDyn              =	997
INTEGER(4), PARAMETER                  :: TmdYDzn              =	998
INTEGER(4), PARAMETER                  :: TmdYVxn              =	999
INTEGER(4), PARAMETER                  :: TmdYVyn              =	1000
INTEGER(4), PARAMETER                  :: TmdYVzn              =	1001
INTEGER(4), PARAMETER                  :: TmdYAxn              =	1002
INTEGER(4), PARAMETER                  :: TmdYAyn              =	1003
INTEGER(4), PARAMETER                  :: TmdYAzn              =	1004


!mal End of proposed change.  v6.10a-mal  1-Mar-2009.


!mal Start of proposed change.  v6.40a-mal  10-Oct-2009.
!mal add additional output parameters for the tmd if it is in the platform instead of nacelle

INTEGER(4), PARAMETER                  :: TmdXDxt              =	1005
INTEGER(4), PARAMETER                  :: TmdXDyt              =	1006
INTEGER(4), PARAMETER                  :: TmdXDzt              =	1007
INTEGER(4), PARAMETER                  :: TmdXVxt              =	1008
INTEGER(4), PARAMETER                  :: TmdXVyt              =	1009
INTEGER(4), PARAMETER                  :: TmdXVzt              =	1010
INTEGER(4), PARAMETER                  :: TmdXAxt              =	1011
INTEGER(4), PARAMETER                  :: TmdXAyt              =	1012
INTEGER(4), PARAMETER                  :: TmdXAzt              =	1013
INTEGER(4), PARAMETER                  :: TmdYDxt              =	1014
INTEGER(4), PARAMETER                  :: TmdYDyt              =	1015
INTEGER(4), PARAMETER                  :: TmdYDzt              =	1016
INTEGER(4), PARAMETER                  :: TmdYVxt              =	1017
INTEGER(4), PARAMETER                  :: TmdYVyt              =	1018
INTEGER(4), PARAMETER                  :: TmdYVzt              =	1019
INTEGER(4), PARAMETER                  :: TmdYAxt              =	1020
INTEGER(4), PARAMETER                  :: TmdYAyt              =	1021
INTEGER(4), PARAMETER                  :: TmdYAzt              =	1022

   
     ! The maximum number of output channels which can be output by the code.
   INTEGER, PARAMETER                  :: MaxOutPts =  1022

!End of code generated by Matlab script



   ! Regular variables:

! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE VARIABLE BELOW:
REAL(ReKi)                  ,SAVE :: AllOuts  (0:MaxOutPts)                          ! An array holding the value of all of the calculated (not only selected) output channels.
! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE PREVIOUS VARIABLE:
REAL(DbKi), ALLOCATABLE     ,SAVE :: TimeData (:)                                    ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
REAL(ReKi), ALLOCATABLE     ,SAVE :: AllOutData (:,:)                                ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LinAccES (:,:,:)                                ! Total linear acceleration of a point on a   blade (point S) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: LinAccET (:,:)                                  ! Total linear acceleration of a point on the tower (point T) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: LSNodes  (:,:)                                  ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters).
REAL(ReKi), ALLOCATABLE     ,SAVE :: FrcS0B   (:,:)                                  ! Total force at the blade root (point S(0)) due to the blade.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FTHydro  (:,:)                                  ! Total hydrodynamic force per unit length acting on the tower at point T.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MFHydro  (:,:)                                  ! Total hydrodynamic moment per unit length acting on a tower element (body F) at point T.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MomH0B   (:,:)                                  ! Total moment at the hub (body H) / blade root (point S(0)) due to the blade.
REAL(ReKi)                  ,SAVE :: NcIMUxn                                         ! Downwind distance from the tower-top to the nacelle IMU.
REAL(ReKi)                  ,SAVE :: NcIMUyn                                         ! Lateral  distance from the tower-top to the nacelle IMU.
REAL(ReKi)                  ,SAVE :: NcIMUzn                                         ! Vertical distance from the tower-top to the nacelle IMU.
REAL(ReKi), ALLOCATABLE     ,SAVE :: OutData  (:)                                    ! Array to contain the output data.
REAL(ReKi)                  ,SAVE :: ShftGagL                                        ! Distance from hub or teeter pin to shaft strain gages.
REAL(ReKi)                  ,SAVE :: SttsTime                                        ! Amount of time between screen status messages (sec).
REAL(ReKi)                  ,SAVE :: TStart                                          ! Time to begin tabular output.
REAL(ReKi)                  ,SAVE :: WaveElevxi(1) = (/ 0.0 /)                       ! xi-coordinates for points where the incident wave elevation can be output (meters).
REAL(ReKi)                  ,SAVE :: WaveElevyi(1) = (/ 0.0 /)                       ! yi-coordinates for points where the incident wave elevation can be output (meters).

INTEGER(4)                  ,SAVE :: BldGagNd (9)                                    ! Nodes closest to the blade strain gages.
INTEGER                     ,SAVE :: CurrOutStep                                     ! Time index into the AllOutData arrat
INTEGER(4)                  ,SAVE :: DecFact                                         ! Decimation factor for tabular output.
!JASON: ADD OUTPUTS FOR THE MOORING LINE POSITION AND EFFECTIVE TENSION AT EACH NODE.  USE NAMES SUCH AS: Ln#Nd#Pxi, Ln#Nd#Pyi, Ln#Nd#Pzi, Ln#Nd#Ten WHERE # REPRESENTS THE LINE NUMBER OR NODE NUMBER!!!
INTEGER(4)                  ,SAVE :: LineNodes    = 0                                ! Number of nodes per line where the mooring line position and tension can be output (-).
INTEGER(4)                  ,SAVE :: NBlGages                                        ! Number of blade strain gages.
INTEGER(4)                  ,SAVE :: NTwGages                                        ! Number of tower strain gages.
INTEGER(4)                  ,SAVE :: NumOuts                                         ! Number of parameters in the output list.
INTEGER(IntKi)              ,SAVE :: NOutSteps                                       ! Maximum number of output steps
INTEGER(4)                  ,SAVE :: NWaveElev    = 1                                ! Number of points where the incident wave elevation  can be output (-).
INTEGER(4)                  ,SAVE :: NWaveKin     = 0                                ! Number of points where the incident wave kinematics can be output (-).
INTEGER(4)                  ,SAVE :: TwrGagNd (9)                                    ! Nodes closest to the tower strain gages.
INTEGER(4)                  ,SAVE :: WaveKinNd(9)                                    ! List of tower [not floating] or platform [floating] nodes that have wave kinematics sensors.

INTEGER(B2Ki), PARAMETER                  :: FileFmtID_WithTime    = 1                       ! ID for OutputFileFmtID to specify that the time channel should be included in the output file (use if the output can occur at variable times)
INTEGER(B2Ki), PARAMETER                  :: FileFmtID_WithoutTime = 2                       ! ID for OutputFileFmtID to specify that the time channel does not need to be included in the output file (used only with constant time-step output)
INTEGER(B2Ki)               ,SAVE :: OutputFileFmtID = FileFmtID_WithoutTime         ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)

LOGICAL                     ,SAVE :: WrEcho
LOGICAL                     ,SAVE :: WrBinOutFile  = .true.                          ! Write a binary output file? (.outb)
LOGICAL                     ,SAVE :: WrTxtOutFile  = .true.                          ! Write a text (formatted) output file? (.out)
LOGICAL                     ,SAVE :: TabDelim                                        ! Flag to cause tab-delimited output.

CHARACTER(20)               ,SAVE :: OutFmt                                          ! Output format for tabular data.
CHARACTER(10)               ,SAVE :: OutList  (MaxOutPts)                            ! List of output parameters.
CHARACTER(1024)             ,SAVE :: FileDesc                                        ! Description of run to include in binary output file

TYPE(OutParmType),ALLOCATABLE, SAVE :: OutParam (:)                                    ! Names and units of all output parameters.

   ! Let's make these parameters into arrays so we can loop through them
   
INTEGER, PARAMETER                  :: WaveVxi(9) = (/Wave1Vxi,Wave2Vxi,Wave3Vxi,Wave4Vxi,Wave5Vxi,Wave6Vxi,Wave7Vxi,Wave8Vxi,Wave9Vxi/)
INTEGER, PARAMETER                  :: WaveVyi(9) = (/Wave1Vyi,Wave2Vyi,Wave3Vyi,Wave4Vyi,Wave5Vyi,Wave6Vyi,Wave7Vyi,Wave8Vyi,Wave9Vyi/)
INTEGER, PARAMETER                  :: WaveVzi(9) = (/Wave1Vzi,Wave2Vzi,Wave3Vzi,Wave4Vzi,Wave5Vzi,Wave6Vzi,Wave7Vzi,Wave8Vzi,Wave9Vzi/)
INTEGER, PARAMETER                  :: WaveAxi(9) = (/Wave1Axi,Wave2Axi,Wave3Axi,Wave4Axi,Wave5Axi,Wave6Axi,Wave7Axi,Wave8Axi,Wave9Axi/)
INTEGER, PARAMETER                  :: WaveAyi(9) = (/Wave1Ayi,Wave2Ayi,Wave3Ayi,Wave4Ayi,Wave5Ayi,Wave6Ayi,Wave7Ayi,Wave8Ayi,Wave9Ayi/)
INTEGER, PARAMETER                  :: WaveAzi(9) = (/Wave1Azi,Wave2Azi,Wave3Azi,Wave4Azi,Wave5Azi,Wave6Azi,Wave7Azi,Wave8Azi,Wave9Azi/)
INTEGER, PARAMETER                  :: FairTen(9) = (/Fair1Ten,Fair2Ten,Fair3Ten,Fair4Ten,Fair5Ten,Fair6Ten,Fair7Ten,Fair8Ten,Fair9Ten/)
INTEGER, PARAMETER                  :: FairAng(9) = (/Fair1Ang,Fair2Ang,Fair3Ang,Fair4Ang,Fair5Ang,Fair6Ang,Fair7Ang,Fair8Ang,Fair9Ang/)
INTEGER, PARAMETER                  :: AnchTen(9) = (/Anch1Ten,Anch2Ten,Anch3Ten,Anch4Ten,Anch5Ten,Anch6Ten,Anch7Ten,Anch8Ten,Anch9Ten/)
INTEGER, PARAMETER                  :: AnchAng(9) = (/Anch1Ang,Anch2Ang,Anch3Ang,Anch4Ang,Anch5Ang,Anch6Ang,Anch7Ang,Anch8Ang,Anch9Ang/)

INTEGER, PARAMETER                  :: TipDxc( 3)  = (/TipDxc1,  TipDxc2,  TipDxc3/)
INTEGER, PARAMETER                  :: TipDyc( 3)  = (/TipDyc1,  TipDyc2,  TipDyc3/)
INTEGER, PARAMETER                  :: TipDzc( 3)  = (/TipDzc1,  TipDzc2,  TipDzc3/)
INTEGER, PARAMETER                  :: TipDxb( 3)  = (/TipDxb1,  TipDxb2,  TipDxb3/)
INTEGER, PARAMETER                  :: TipDyb( 3)  = (/TipDyb1,  TipDyb2,  TipDyb3/)
INTEGER, PARAMETER                  :: TipALxb(3)  = (/TipALxb1, TipALxb2, TipALxb3/)
INTEGER, PARAMETER                  :: TipALyb(3)  = (/TipALyb1, TipALyb2, TipALyb3/)
INTEGER, PARAMETER                  :: TipALzb(3)  = (/TipALzb1, TipALzb2, TipALzb3/)
INTEGER, PARAMETER                  :: TipRDxb(3)  = (/TipRDxb1, TipRDxb2, TipRDxb3/)
INTEGER, PARAMETER                  :: TipRDyb(3)  = (/TipRDyb1, TipRDyb2, TipRDyb3/)
INTEGER, PARAMETER                  :: TipRDzc(3)  = (/TipRDzc1, TipRDzc2, TipRDzc3/)
INTEGER, PARAMETER                  :: TipClrnc(3) = (/TipClrnc1,TipClrnc2,TipClrnc3/)
INTEGER, PARAMETER                  :: PtchPMzc(3) = (/PtchPMzc1,PtchPMzc2,PtchPMzc3/)

INTEGER, PARAMETER                  :: RootFxc(3) = (/ RootFxc1,RootFxc2,RootFxc3 /)
INTEGER, PARAMETER                  :: RootFyc(3) = (/ RootFyc1,RootFyc2,RootFyc3 /)
INTEGER, PARAMETER                  :: RootFzc(3) = (/ RootFzc1,RootFzc2,RootFzc3 /)
INTEGER, PARAMETER                  :: RootFxb(3) = (/ RootFxb1,RootFxb2,RootFxb3 /)
INTEGER, PARAMETER                  :: RootFyb(3) = (/ RootFyb1,RootFyb2,RootFyb3 /)
INTEGER, PARAMETER                  :: RootMxc(3) = (/ RootMxc1,RootMxc2,RootMxc3 /)
INTEGER, PARAMETER                  :: RootMyc(3) = (/ RootMyc1,RootMyc2,RootMyc3 /)
INTEGER, PARAMETER                  :: RootMzc(3) = (/ RootMzc1,RootMzc2,RootMzc3 /)
INTEGER, PARAMETER                  :: RootMxb(3) = (/ RootMxb1,RootMxb2,RootMxb3 /)
INTEGER, PARAMETER                  :: RootMyb(3) = (/ RootMyb1,RootMyb2,RootMyb3 /)

INTEGER, PARAMETER                  :: SpnALxb(9, 3) = RESHAPE( (/ &
                                    Spn1ALxb1,Spn2ALxb1,Spn3ALxb1,Spn4ALxb1,Spn5ALxb1,Spn6ALxb1,Spn7ALxb1,Spn8ALxb1,Spn9ALxb1, &
                                    Spn1ALxb2,Spn2ALxb2,Spn3ALxb2,Spn4ALxb2,Spn5ALxb2,Spn6ALxb2,Spn7ALxb2,Spn8ALxb2,Spn9ALxb2, &
                                    Spn1ALxb3,Spn2ALxb3,Spn3ALxb3,Spn4ALxb3,Spn5ALxb3,Spn6ALxb3,Spn7ALxb3,Spn8ALxb3,Spn9ALxb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnALyb(9, 3) = RESHAPE( (/ &
                                    Spn1ALyb1,Spn2ALyb1,Spn3ALyb1,Spn4ALyb1,Spn5ALyb1,Spn6ALyb1,Spn7ALyb1,Spn8ALyb1,Spn9ALyb1, &
                                    Spn1ALyb2,Spn2ALyb2,Spn3ALyb2,Spn4ALyb2,Spn5ALyb2,Spn6ALyb2,Spn7ALyb2,Spn8ALyb2,Spn9ALyb2, &
                                    Spn1ALyb3,Spn2ALyb3,Spn3ALyb3,Spn4ALyb3,Spn5ALyb3,Spn6ALyb3,Spn7ALyb3,Spn8ALyb3,Spn9ALyb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnALzb(9, 3) = RESHAPE( (/ &
                                    Spn1ALzb1,Spn2ALzb1,Spn3ALzb1,Spn4ALzb1,Spn5ALzb1,Spn6ALzb1,Spn7ALzb1,Spn8ALzb1,Spn9ALzb1, &
                                    Spn1ALzb2,Spn2ALzb2,Spn3ALzb2,Spn4ALzb2,Spn5ALzb2,Spn6ALzb2,Spn7ALzb2,Spn8ALzb2,Spn9ALzb2, &
                                    Spn1ALzb3,Spn2ALzb3,Spn3ALzb3,Spn4ALzb3,Spn5ALzb3,Spn6ALzb3,Spn7ALzb3,Spn8ALzb3,Spn9ALzb3  &
                                /), (/9, 3/) )

INTEGER, PARAMETER                  :: SpnFLxb(9,3) = RESHAPE( (/ &
                                    Spn1FLxb1,Spn2FLxb1,Spn3FLxb1,Spn4FLxb1,Spn5FLxb1,Spn6FLxb1,Spn7FLxb1,Spn8FLxb1,Spn9FLxb1, &
                                    Spn1FLxb2,Spn2FLxb2,Spn3FLxb2,Spn4FLxb2,Spn5FLxb2,Spn6FLxb2,Spn7FLxb2,Spn8FLxb2,Spn9FLxb2, &
                                    Spn1FLxb3,Spn2FLxb3,Spn3FLxb3,Spn4FLxb3,Spn5FLxb3,Spn6FLxb3,Spn7FLxb3,Spn8FLxb3,Spn9FLxb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnFLyb(9,3) = RESHAPE( (/ &
                                    Spn1FLyb1,Spn2FLyb1,Spn3FLyb1,Spn4FLyb1,Spn5FLyb1,Spn6FLyb1,Spn7FLyb1,Spn8FLyb1,Spn9FLyb1, &
                                    Spn1FLyb2,Spn2FLyb2,Spn3FLyb2,Spn4FLyb2,Spn5FLyb2,Spn6FLyb2,Spn7FLyb2,Spn8FLyb2,Spn9FLyb2, &
                                    Spn1FLyb3,Spn2FLyb3,Spn3FLyb3,Spn4FLyb3,Spn5FLyb3,Spn6FLyb3,Spn7FLyb3,Spn8FLyb3,Spn9FLyb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnFLzb(9,3) = RESHAPE( (/ &
                                    Spn1FLzb1,Spn2FLzb1,Spn3FLzb1,Spn4FLzb1,Spn5FLzb1,Spn6FLzb1,Spn7FLzb1,Spn8FLzb1,Spn9FLzb1, &
                                    Spn1FLzb2,Spn2FLzb2,Spn3FLzb2,Spn4FLzb2,Spn5FLzb2,Spn6FLzb2,Spn7FLzb2,Spn8FLzb2,Spn9FLzb2, &
                                    Spn1FLzb3,Spn2FLzb3,Spn3FLzb3,Spn4FLzb3,Spn5FLzb3,Spn6FLzb3,Spn7FLzb3,Spn8FLzb3,Spn9FLzb3  &
                                /), (/9, 3/) )
                                
INTEGER, PARAMETER                  :: SpnMLxb(9,3) = RESHAPE( (/ &
                                    Spn1MLxb1,Spn2MLxb1,Spn3MLxb1,Spn4MLxb1,Spn5MLxb1,Spn6MLxb1,Spn7MLxb1,Spn8MLxb1,Spn9MLxb1, &
                                    Spn1MLxb2,Spn2MLxb2,Spn3MLxb2,Spn4MLxb2,Spn5MLxb2,Spn6MLxb2,Spn7MLxb2,Spn8MLxb2,Spn9MLxb2, &
                                    Spn1MLxb3,Spn2MLxb3,Spn3MLxb3,Spn4MLxb3,Spn5MLxb3,Spn6MLxb3,Spn7MLxb3,Spn8MLxb3,Spn9MLxb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnMLyb(9,3) = RESHAPE( (/ &
                                    Spn1MLyb1,Spn2MLyb1,Spn3MLyb1,Spn4MLyb1,Spn5MLyb1,Spn6MLyb1,Spn7MLyb1,Spn8MLyb1,Spn9MLyb1, &
                                    Spn1MLyb2,Spn2MLyb2,Spn3MLyb2,Spn4MLyb2,Spn5MLyb2,Spn6MLyb2,Spn7MLyb2,Spn8MLyb2,Spn9MLyb2, &
                                    Spn1MLyb3,Spn2MLyb3,Spn3MLyb3,Spn4MLyb3,Spn5MLyb3,Spn6MLyb3,Spn7MLyb3,Spn8MLyb3,Spn9MLyb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnMLzb(9,3) = RESHAPE( (/ &
                                    Spn1MLzb1,Spn2MLzb1,Spn3MLzb1,Spn4MLzb1,Spn5MLzb1,Spn6MLzb1,Spn7MLzb1,Spn8MLzb1,Spn9MLzb1, &
                                    Spn1MLzb2,Spn2MLzb2,Spn3MLzb2,Spn4MLzb2,Spn5MLzb2,Spn6MLzb2,Spn7MLzb2,Spn8MLzb2,Spn9MLzb2, &
                                    Spn1MLzb3,Spn2MLzb3,Spn3MLzb3,Spn4MLzb3,Spn5MLzb3,Spn6MLzb3,Spn7MLzb3,Spn8MLzb3,Spn9MLzb3  &
                                /), (/9, 3/) )
                                
INTEGER, PARAMETER                  :: SpnTDxb(9,3) = RESHAPE( (/ &
                                    Spn1TDxb1,Spn2TDxb1,Spn3TDxb1,Spn4TDxb1,Spn5TDxb1,Spn6TDxb1,Spn7TDxb1,Spn8TDxb1,Spn9TDxb1, &
                                    Spn1TDxb2,Spn2TDxb2,Spn3TDxb2,Spn4TDxb2,Spn5TDxb2,Spn6TDxb2,Spn7TDxb2,Spn8TDxb2,Spn9TDxb2, &
                                    Spn1TDxb3,Spn2TDxb3,Spn3TDxb3,Spn4TDxb3,Spn5TDxb3,Spn6TDxb3,Spn7TDxb3,Spn8TDxb3,Spn9TDxb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnTDyb(9,3) = RESHAPE( (/ &
                                    Spn1TDyb1,Spn2TDyb1,Spn3TDyb1,Spn4TDyb1,Spn5TDyb1,Spn6TDyb1,Spn7TDyb1,Spn8TDyb1,Spn9TDyb1, &
                                    Spn1TDyb2,Spn2TDyb2,Spn3TDyb2,Spn4TDyb2,Spn5TDyb2,Spn6TDyb2,Spn7TDyb2,Spn8TDyb2,Spn9TDyb2, &
                                    Spn1TDyb3,Spn2TDyb3,Spn3TDyb3,Spn4TDyb3,Spn5TDyb3,Spn6TDyb3,Spn7TDyb3,Spn8TDyb3,Spn9TDyb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnTDzb(9,3) = RESHAPE( (/ &
                                    Spn1TDzb1,Spn2TDzb1,Spn3TDzb1,Spn4TDzb1,Spn5TDzb1,Spn6TDzb1,Spn7TDzb1,Spn8TDzb1,Spn9TDzb1, &
                                    Spn1TDzb2,Spn2TDzb2,Spn3TDzb2,Spn4TDzb2,Spn5TDzb2,Spn6TDzb2,Spn7TDzb2,Spn8TDzb2,Spn9TDzb2, &
                                    Spn1TDzb3,Spn2TDzb3,Spn3TDzb3,Spn4TDzb3,Spn5TDzb3,Spn6TDzb3,Spn7TDzb3,Spn8TDzb3,Spn9TDzb3  &
                                /), (/9, 3/) )

INTEGER, PARAMETER                  :: SpnRDxb(9,3) = RESHAPE( (/ &
                                    Spn1RDxb1,Spn2RDxb1,Spn3RDxb1,Spn4RDxb1,Spn5RDxb1,Spn6RDxb1,Spn7RDxb1,Spn8RDxb1,Spn9RDxb1, &
                                    Spn1RDxb2,Spn2RDxb2,Spn3RDxb2,Spn4RDxb2,Spn5RDxb2,Spn6RDxb2,Spn7RDxb2,Spn8RDxb2,Spn9RDxb2, &
                                    Spn1RDxb3,Spn2RDxb3,Spn3RDxb3,Spn4RDxb3,Spn5RDxb3,Spn6RDxb3,Spn7RDxb3,Spn8RDxb3,Spn9RDxb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnRDyb(9,3) = RESHAPE( (/ &
                                    Spn1RDyb1,Spn2RDyb1,Spn3RDyb1,Spn4RDyb1,Spn5RDyb1,Spn6RDyb1,Spn7RDyb1,Spn8RDyb1,Spn9RDyb1, &
                                    Spn1RDyb2,Spn2RDyb2,Spn3RDyb2,Spn4RDyb2,Spn5RDyb2,Spn6RDyb2,Spn7RDyb2,Spn8RDyb2,Spn9RDyb2, &
                                    Spn1RDyb3,Spn2RDyb3,Spn3RDyb3,Spn4RDyb3,Spn5RDyb3,Spn6RDyb3,Spn7RDyb3,Spn8RDyb3,Spn9RDyb3  &
                                /), (/9, 3/) )
INTEGER, PARAMETER                  :: SpnRDzb(9,3) = RESHAPE( (/ &
                                    Spn1RDzb1,Spn2RDzb1,Spn3RDzb1,Spn4RDzb1,Spn5RDzb1,Spn6RDzb1,Spn7RDzb1,Spn8RDzb1,Spn9RDzb1, &
                                    Spn1RDzb2,Spn2RDzb2,Spn3RDzb2,Spn4RDzb2,Spn5RDzb2,Spn6RDzb2,Spn7RDzb2,Spn8RDzb2,Spn9RDzb2, &
                                    Spn1RDzb3,Spn2RDzb3,Spn3RDzb3,Spn4RDzb3,Spn5RDzb3,Spn6RDzb3,Spn7RDzb3,Spn8RDzb3,Spn9RDzb3  &
                                /), (/9, 3/) )

                                
INTEGER, PARAMETER                  :: TwHtALxt(9) = (/ &
                                    TwHt1ALxt,TwHt2ALxt,TwHt3ALxt,TwHt4ALxt,TwHt5ALxt,TwHt6ALxt,TwHt7ALxt,TwHt8ALxt,TwHt9ALxt /)
INTEGER, PARAMETER                  :: TwHtALyt(9) = (/ &
                                    TwHt1ALyt,TwHt2ALyt,TwHt3ALyt,TwHt4ALyt,TwHt5ALyt,TwHt6ALyt,TwHt7ALyt,TwHt8ALyt,TwHt9ALyt /)
INTEGER, PARAMETER                  :: TwHtALzt(9) = (/ &
                                    TwHt1ALzt,TwHt2ALzt,TwHt3ALzt,TwHt4ALzt,TwHt5ALzt,TwHt6ALzt,TwHt7ALzt,TwHt8ALzt,TwHt9ALzt /)

INTEGER, PARAMETER                  :: TwHtMLxt(9) = (/ &
                                    TwHt1MLxt,TwHt2MLxt,TwHt3MLxt,TwHt4MLxt,TwHt5MLxt,TwHt6MLxt,TwHt7MLxt,TwHt8MLxt,TwHt9MLxt /)
INTEGER, PARAMETER                  :: TwHtMLyt(9) = (/ &
                                    TwHt1MLyt,TwHt2MLyt,TwHt3MLyt,TwHt4MLyt,TwHt5MLyt,TwHt6MLyt,TwHt7MLyt,TwHt8MLyt,TwHt9MLyt /)
INTEGER, PARAMETER                  :: TwHtMLzt(9) = (/ &
                                    TwHt1MLzt,TwHt2MLzt,TwHt3MLzt,TwHt4MLzt,TwHt5MLzt,TwHt6MLzt,TwHt7MLzt,TwHt8MLzt,TwHt9MLzt /)

INTEGER, PARAMETER                  :: TwHtFLxt(9) = (/ &
                                    TwHt1FLxt,TwHt2FLxt,TwHt3FLxt,TwHt4FLxt,TwHt5FLxt,TwHt6FLxt,TwHt7FLxt,TwHt8FLxt,TwHt9FLxt /)
INTEGER, PARAMETER                  :: TwHtFLyt(9) = (/ &
                                    TwHt1FLyt,TwHt2FLyt,TwHt3FLyt,TwHt4FLyt,TwHt5FLyt,TwHt6FLyt,TwHt7FLyt,TwHt8FLyt,TwHt9FLyt /)
INTEGER, PARAMETER                  :: TwHtFLzt(9) = (/ &
                                    TwHt1FLzt,TwHt2FLzt,TwHt3FLzt,TwHt4FLzt,TwHt5FLzt,TwHt6FLzt,TwHt7FLzt,TwHt8FLzt,TwHt9FLzt /)

INTEGER, PARAMETER                  :: TwHtTDxt(9) = (/ &
                                    TwHt1TDxt,TwHt2TDxt,TwHt3TDxt,TwHt4TDxt,TwHt5TDxt,TwHt6TDxt,TwHt7TDxt,TwHt8TDxt,TwHt9TDxt /)
INTEGER, PARAMETER                  :: TwHtTDyt(9) = (/ &
                                    TwHt1TDyt,TwHt2TDyt,TwHt3TDyt,TwHt4TDyt,TwHt5TDyt,TwHt6TDyt,TwHt7TDyt,TwHt8TDyt,TwHt9TDyt /)
INTEGER, PARAMETER                  :: TwHtTDzt(9) = (/ &
                                    TwHt1TDzt,TwHt2TDzt,TwHt3TDzt,TwHt4TDzt,TwHt5TDzt,TwHt6TDzt,TwHt7TDzt,TwHt8TDzt,TwHt9TDzt /)

INTEGER, PARAMETER                  :: TwHtRDxt(9) = (/ &
                                    TwHt1RDxt,TwHt2RDxt,TwHt3RDxt,TwHt4RDxt,TwHt5RDxt,TwHt6RDxt,TwHt7RDxt,TwHt8RDxt,TwHt9RDxt /)
INTEGER, PARAMETER                  :: TwHtRDyt(9) = (/ &
                                    TwHt1RDyt,TwHt2RDyt,TwHt3RDyt,TwHt4RDyt,TwHt5RDyt,TwHt6RDyt,TwHt7RDyt,TwHt8RDyt,TwHt9RDyt /)
INTEGER, PARAMETER                  :: TwHtRDzt(9) = (/ &
                                    TwHt1RDzt,TwHt2RDzt,TwHt3RDzt,TwHt4RDzt,TwHt5RDzt,TwHt6RDzt,TwHt7RDzt,TwHt8RDzt,TwHt9RDzt /)

INTEGER, PARAMETER                  :: TwHtTPxi(9) = (/ &
                                    TwHt1TPxi,TwHt2TPxi,TwHt3TPxi,TwHt4TPxi,TwHt5TPxi,TwHt6TPxi,TwHt7TPxi,TwHt8TPxi,TwHt9TPxi /)
INTEGER, PARAMETER                  :: TwHtTPyi(9) = (/ &
                                    TwHt1TPyi,TwHt2TPyi,TwHt3TPyi,TwHt4TPyi,TwHt5TPyi,TwHt6TPyi,TwHt7TPyi,TwHt8TPyi,TwHt9TPyi /)
INTEGER, PARAMETER                  :: TwHtTPzi(9) = (/ &
                                    TwHt1TPzi,TwHt2TPzi,TwHt3TPzi,TwHt4TPzi,TwHt5TPzi,TwHt6TPzi,TwHt7TPzi,TwHt8TPzi,TwHt9TPzi /)

INTEGER, PARAMETER                  :: TwHtRPxi(9) = (/ &
                                    TwHt1RPxi,TwHt2RPxi,TwHt3RPxi,TwHt4RPxi,TwHt5RPxi,TwHt6RPxi,TwHt7RPxi,TwHt8RPxi,TwHt9RPxi /)
INTEGER, PARAMETER                  :: TwHtRPyi(9) = (/ &
                                    TwHt1RPyi,TwHt2RPyi,TwHt3RPyi,TwHt4RPyi,TwHt5RPyi,TwHt6RPyi,TwHt7RPyi,TwHt8RPyi,TwHt9RPyi /)
INTEGER, PARAMETER                  :: TwHtRPzi(9) = (/ &
                                    TwHt1RPzi,TwHt2RPzi,TwHt3RPzi,TwHt4RPzi,TwHt5RPzi,TwHt6RPzi,TwHt7RPzi,TwHt8RPzi,TwHt9RPzi /)
   
END MODULE Output
!=======================================================================
MODULE Platform


   ! This MODULE stores input variables for platform loading.


USE                             Precision

REAL(ReKi), ALLOCATABLE     ,SAVE :: LAnchxi  (:)                                    ! xi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LAnchyi  (:)                                    ! yi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LAnchzi  (:)                                    ! zi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LDiam    (:)                                    ! Effective diameter of each mooring line for calculation of the line buoyancy.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LEAStff  (:)                                    ! Extensional stiffness of each mooring line.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LFairxt  (:)                                    ! xt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LFairyt  (:)                                    ! yt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LFairzt  (:)                                    ! zt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LMassDen (:)                                    ! Mass density of each mooring line.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LSeabedCD(:)                                    ! Coefficient of seabed static friction drag of each mooring line (a negative value indicates no seabed).
REAL(ReKi), ALLOCATABLE     ,SAVE :: LTenTol  (:)                                    ! Convergence tolerance within Newton-Raphson iteration of each mooring line specified as a fraction of tension.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LUnstrLen(:)                                    ! Unstretched length of each mooring line.
REAL(ReKi)                  ,SAVE :: MaxLRadAnch  = 0.0                              ! Maximum value of input array LRadAnch.

REAL(ReKi)                  ,SAVE :: PtfmAM (6,6) = 0.0                              ! Platform added mass matrix.
REAL(ReKi)                  ,SAVE :: PtfmCD                                          ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                  ,SAVE :: PtfmDiam                                        ! Effective platform diameter in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                  ,SAVE :: PtfmDraft                                       ! Effective platform draft    in calculation of viscous drag term from Morison's equation.

REAL(ReKi)                  ,SAVE :: PtfmFt   (6) = 0.0                              ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the platform force at the platform reference (point Z) and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the platform moment acting at the platform (body X) / platform reference (point Z) associated with everything but the QD2T()'s.

REAL(ReKi)                  ,SAVE :: PtfmVol0                                        ! Displaced volume of water when the platform is in its undisplaced position.
REAL(ReKi)                  ,SAVE :: RdtnDT                                          ! Time step for wave radiation kernel calculations.
REAL(ReKi)                  ,SAVE :: RdtnTMax     = 0.0                              ! Analysis time for wave radiation kernel calculations.

INTEGER(4)                  ,SAVE :: LineMod                                         ! Mooring line model switch.
INTEGER(4)                  ,SAVE :: NumLines     = 0                                ! Number of mooring lines.

INTEGER(4)                  ,SAVE :: PtfmLdMod    = 0                                ! Platform loading model switch. (Initialized to zero b/c not all models read in PtfmFile)
INTEGER(4)                  ,SAVE :: PtfmNodes                                       ! Number of platform nodes used in calculation of viscous drag term from Morison's equation.

CHARACTER(1024)             ,SAVE :: WAMITFile                                       ! Root name of WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst extension), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1 extension), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3 extension).


END MODULE Platform
!=======================================================================
MODULE RotorFurling


   ! This MODULE stores input variables for rotor-furling.


USE                             Precision


REAL(ReKi)   :: RFrlCDmp  = 0.0                                 ! Rotor-furl rate-independent Coulomb-damping moment. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlDmp   = 0.0                                 ! Rotor-furl damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlDSDmp = 0.0                                 ! Rotor-furl down-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlDSDP  = 0.0                                 ! Rotor-furl down-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlDSSP  = 0.0                                 ! Rotor-furl down-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlDSSpr = 0.0                                 ! Rotor-furl down-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlSpr   = 0.0                                 ! Rotor-furl spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlUSDmp = 0.0                                 ! Rotor-furl up-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlUSDP  = 0.0                                 ! Rotor-furl up-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlUSSP  = 0.0                                 ! Rotor-furl up-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)   :: RFrlUSSpr = 0.0                                 ! Rotor-furl up-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)   :: RFrlMod   = 0                                   ! Rotor-furl spring/damper model switch. (Initialized to zero b/c not all models read in FurlFile)


END MODULE RotorFurling
!=======================================================================
MODULE RtHndSid


   ! This MODULE stores variables used in RtHS.


USE                             Precision


REAL(ReKi)                  ,SAVE :: AngAccEBt(3)                                    ! Portion of the angular acceleration of the base plate                                                (body B) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: AngAccERt(3)                                    ! Portion of the angular acceleration of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: AngAccEXt(3)                                    ! Portion of the angular acceleration of the platform                                                  (body X) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: AngAccEFt(:,:)                                  ! Portion of the angular acceleration of tower element J                                               (body F) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.

REAL(ReKi), ALLOCATABLE     ,SAVE :: AngPosEF (:,:)                                  ! Angular position of the current point on the tower                                (body F) in the inertial frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: AngPosXF (:,:)                                  ! Angular position of the current point on the tower                                (body F) in the platform      (body X          ).
REAL(ReKi), ALLOCATABLE     ,SAVE :: AngPosHM (:,:,:)                                ! Angular position of eleMent J of blade K                                          (body M) in the hub           (body H          ).
REAL(ReKi)                  ,SAVE :: AngPosXB (3)                                    ! Angular position of the base plate                                                (body B) in the platform      (body X          ).
REAL(ReKi)                  ,SAVE :: AngVelEB (3)                                    ! Angular velocity of the base plate                                                (body B) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: AngVelEF  (:,:)                                 ! Angular velocity of the current point on the tower                                (body F) in the inertia frame (body E for earth).
REAL(ReKi)                  ,SAVE :: AngVelER (3)                                    ! Angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth).
REAL(ReKi)                  ,SAVE :: AngVelEX (3)                                    ! Angular velocity of the platform                                                  (body X) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: AugMat   (:,:)                                  ! The augmented matrix used for the solution of the QD2T()s.
REAL(ReKi)                  ,SAVE :: FKAero   (3)                                    ! The tail fin aerodynamic force acting at point K, the center-of-pressure of the tail fin.
REAL(ReKi)                  ,SAVE :: FrcONcRtt(3)                                    ! Portion of the force at yaw bearing         (point O   ) due to the nacelle, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: FrcPRott (3)                                    ! Portion of the force at the teeter pin      (point P   ) due to the rotor associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FrcS0Bt  (:,:)                                  ! Portion of the force at the blade root      (point S(0)) due to the blade associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: FrcT0Trbt(3)                                    ! Portion of the force at tower base          (point T(0)) due to the turbine associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: FrcVGnRtt(3)                                    ! Portion of the force at the rotor-furl axis (point V   ) due to the structure that furls with the rotor, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: FrcWTailt(3)                                    ! Portion of the force at the  tail-furl axis (point W   ) due to the tail associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: FrcZAllt (3)                                    ! Portion of the force at platform reference  (point Z   ) due to everything associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FSAero   (:,:,:)                                ! The aerodynamic force per unit span acting on a blade at point S.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FSTipDrag(:,:)                                  ! The aerodynamic force at a blade tip resulting from tip drag.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FTAero   (:,:)                                  ! The aerodynamic force per unit length acting on the tower at point T.
REAL(ReKi), ALLOCATABLE     ,SAVE :: FTHydrot (:,:)                                  ! Portion of the hydrodynamic force per unit length acting on the tower at point T associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: FZHydrot (3)                                    ! Portion of the platform hydrodynamic force at the platform reference (point Z) associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: GBoxEffFac                                      ! The factor used to apply the gearbox efficiency effects to the equation associated with the generator DOF.
REAL(ReKi)                  ,SAVE :: LinAccEIMUt(3)                                  ! Portion of the linear acceleration of the nacelle IMU      (point IMU) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: LinAccEOt(3)                                    ! Portion of the linear acceleration of the base plate         (point O) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LinAccESt(:,:,:)                                ! Portion of the linear acceleration of a point on a blade     (point S) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LinAccETt(:,:)                                  ! Portion of the linear acceleration of a point on the tower   (point T) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: LinAccEZt (3)                                   ! Portion of the linear acceleration of the platform reference (point Z) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: LinVelEIMU(3)                                   ! Linear velocity of the nacelle IMU  (point IMU) in the inertia frame.
REAL(ReKi)                  ,SAVE :: LinVelEZ  (3)                                   ! Linear velocity of platform reference (point Z) in the inertia frame.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LinVelET  (:,:)                                 ! Linear velocity of current point on the tower         (point T) in the inertia frame.
REAL(ReKi), ALLOCATABLE     ,SAVE :: LinVelESm2 (:)                                  ! The m2-component (closest to tip) of LinVelES.
REAL(ReKi)                  ,SAVE :: MAAero   (3)                                    ! The tail fin aerodynamic moment acting at point K, the center-of-pressure of the tail fin.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MFAero   (:,:)                                  ! The aerodynamic moment per unit length acting on the tower at point T.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MFHydrot (:,:)                                  ! Portion of the hydrodynamic moment per unit length acting on a tower element (body F) at point T associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: MomBNcRtt(3)                                    ! Portion of the moment at the base plate (body B) / yaw bearing                       (point O   ) due to the nacelle, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MomH0Bt  (:,:)                                  ! Portion of the moment at the hub        (body H) / blade root                        (point S(0)) due to the blade associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: MomLPRott(3)                                    ! Portion of the moment at the teeter pin (point P) on the low-speed shaft             (body L    ) due to the rotor associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: MomNGnRtt(3)                                    ! Portion of the moment at the nacelle    (body N) / selected point on rotor-furl axis (point V   ) due the structure that furls with the rotor, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: MomNTailt(3)                                    ! Portion of the moment at the nacelle    (body N) / selected point on  tail-furl axis (point W   ) due the tail associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: MomX0Trbt(3)                                    ! Portion of the moment at the platform   (body X) / tower base                        (point T(0)) due to the turbine associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: MomXAllt (3)                                    ! Portion of the moment at the platform   (body X) / platform reference                (point Z   ) due to everything associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MMAero   (:,:,:)                                ! The aerodynamic moment per unit span acting on a blade at point S.
REAL(ReKi)                  ,SAVE :: MXHydrot (3)                                    ! Portion of the platform hydrodynamic moment acting at the platform (body X) / platform reference (point Z) associated with everything but the  QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEA(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the tail                                                      (body A) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEB(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the base plate                                                (body B) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEF(:,:,:,:)                              ! Partial angular velocity (and its 1st time derivative) of tower element J                                               (body F) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEG(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the generator                                                 (body G) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEH(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the hub                                                       (body H) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEL(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the low-speed shaft                                           (body L) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEM(:,:,:,:,:)                            ! Partial angular velocity (and its 1st time derivative) of eleMent J of blade K                                          (body M) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEN(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the nacelle                                                   (body N) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelER(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PAngVelEX(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the platform                                                (body B) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcONcRt(:,:)                                  ! Partial force at the yaw bearing     (point O   ) due to the nacelle, generator, and rotor.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcPRot (:,:)                                  ! Partial force at the teeter pin      (point P   ) due to the rotor.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcS0B  (:,:,:)                                ! Partial force at the blade root      (point S(0)) due to the blade.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcT0Trb(:,:)                                  ! Partial force at the tower base      (point T(0)) due to the turbine.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcVGnRt(:,:)                                  ! Partial force at the rotor-furl axis (point V   ) due to the structure that furls with the rotor, generator, and rotor.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcWTail(:,:)                                  ! Partial force at the  tail-furl axis (point W   ) due to the tail.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcZAll (:,:)                                  ! Partial force at the platform reference (point Z) due to everything.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFTHydro (:,:,:)                                ! Partial hydrodynamic force per unit length acting on the tower at point T.
REAL(ReKi)                  ,SAVE :: PFZHydro (6,3)                                  ! Partial platform hydrodynamic force at the platform reference (point Z).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEC(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the hub center of mass            (point C) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelED(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the center of mass of the structure that furls with the rotor (not including rotor) (point D) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEI(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the tail boom center of mass                                                        (point I) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEIMU(:,:,:)                              ! Partial linear velocity (and its 1st time derivative) of the nacelle IMU                                                                   (point IMU) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEJ(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the tail fin  center of mass                                                        (point J) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEK(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the tail fin  center of pressure                                                    (point K) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEO(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the base plate                                                                      (point O) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEP(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the teeter pin                                                                      (point P) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEQ(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the apex of rotation                                                                (point Q) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelES(:,:,:,:,:)                            ! Partial linear velocity (and its 1st time derivative) of a point on a blade                                                                  (point S) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelET(:,:,:,:)                              ! Partial linear velocity (and its 1st time derivative) of a point on the tower                                                                (point T) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEU(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the nacelle center of mass                                                          (point U) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEV(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the selected point on the rotor-furl axis                                           (point V) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEW(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the selected point on the  tail-furl axis                                           (point W) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEY(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the platform mass center                                                            (point Y) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelEZ(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the platform reference point                                                        (point Z) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMFHydro (:,:,:)                                ! Partial hydrodynamic moment per unit length acting on a tower element (body F) at point T.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomBNcRt(:,:)                                  ! Partial moment at the base plate (body B) / yaw bearing                       (point O   ) due the nacelle, generator, and rotor.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomH0B  (:,:,:)                                ! Partial moment at the hub        (body H) / blade root                        (point S(0)) due to the blade.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomLPRot(:,:)                                  ! Partial moment at the teeter pin (point P) on the low-speed shaft             (body L    ) due to the rotor.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomNGnRt(:,:)                                  ! Partial moment at the nacelle    (body N) / selected point on rotor-furl axis (point V   ) due the structure that furls with the rotor, generator, and rotor.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomNTail(:,:)                                  ! Partial moment at the nacelle    (body N) / selected point on  tail-furl axis (point W   ) due the tail.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomX0Trb(:,:)                                  ! Partial moment at the platform   (body X) / tower base                        (point T(0)) due to the turbine.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomXAll (:,:)                                  ! Partial moment at the platform   (body X) / platform reference                (point Z   ) due to the everything.
REAL(ReKi)                  ,SAVE :: PMXHydro (6,3)                                  ! Partial platform hydrodynamic moment at the platform (body X) / platform reference (point Z).
REAL(ReKi), ALLOCATABLE     ,SAVE :: QDT      (:)                                    ! Current estimate of QD.
REAL(ReKi), ALLOCATABLE     ,SAVE :: QD2T     (:)                                    ! Solution (acceleration) vector.
REAL(ReKi), ALLOCATABLE     ,SAVE :: QD2TC    (:)                                    ! A copy of the value of QD2T used in SUBROUTINE FixHSSBrTq().
REAL(ReKi), ALLOCATABLE     ,SAVE :: OgnlGeAzRo(:)                                   ! The original elements of AugMat that formed the DOF_GeAz equation before application of known initial conditions.
REAL(ReKi), ALLOCATABLE     ,SAVE :: QT       (:)                                    ! Current estimate of Q for each degree of freedom.
REAL(ReKi)                  ,SAVE :: rO       (3)                                    ! Position vector from inertial frame origin             to tower-top / base plate (point O).
REAL(ReKi), ALLOCATABLE     ,SAVE :: rQS      (:,:,:)                                ! Position vector from the apex of rotation (point Q   ) to a point on a blade (point S).
REAL(ReKi), ALLOCATABLE     ,SAVE :: rS       (:,:,:)                                ! Position vector from inertial frame origin             to a point on a blade (point S).
REAL(ReKi), ALLOCATABLE     ,SAVE :: rS0S     (:,:,:)                                ! Position vector from the blade root       (point S(0)) to a point on a blade (point S).
REAL(ReKi), ALLOCATABLE     ,SAVE :: rT       (:,:)                                  ! Position vector from inertial frame origin to the current node (point T(HNodes(J)).
REAL(ReKi)                  ,SAVE :: rT0O     (3)                                    ! Position vector from the tower base       (point T(0)) to tower-top / base plate (point O).
REAL(ReKi), ALLOCATABLE     ,SAVE :: rT0T     (:,:)                                  ! Position vector from a height of TwrRBHt (base of flexible portion of tower) (point T(0)) to a point on the tower (point T).
REAL(ReKi)                  ,SAVE :: rZ       (3)                                    ! Position vector from inertia frame origin to platform reference (point Z).
REAL(ReKi)                  ,SAVE :: rZO      (3)                                    ! Position vector from platform reference   (point Z   ) to tower-top / base plate (point O).
REAL(ReKi), ALLOCATABLE     ,SAVE :: rZT      (:,:)                                  ! Position vector from platform reference   (point Z   ) to a point on a tower     (point T).
REAL(ReKi), ALLOCATABLE     ,SAVE :: SolnVec  (:)                                    ! Solution vector found by solving the equations of motion.
REAL(ReKi)                  ,SAVE :: TeetAng                                         ! Current teeter angle = QT(DOF_Teet) for 2-blader or 0 for 3-blader (this is used in place of QT(DOF_Teet) throughout RtHS().
REAL(ReKi)                  ,SAVE :: TeetAngVel                                      ! Angular velocity of the teeter motion.

! Yang: Start of adding the kinematic variables and positions vectors of the TMDs 
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelETmdX(:,:,:)                             ! Partial linear velocity (and its 1st time derivative) of the axial TMD                                                        
REAL(ReKi), ALLOCATABLE     ,SAVE :: PLinVelETmdY(:,:,:)                             ! Partial linear velocity (and its 1st time derivative) of the transverse TMD  
REAL(ReKi)                  ,SAVE :: LinVelETmdX  (3)                                ! Linear velocity of TmdX in the inertia frame.
REAL(ReKi)                  ,SAVE :: LinVelETmdY  (3)                                ! Linear velocity of TmdY in the inertia frame.  
REAL(ReKi)                  ,SAVE :: LinAccETmdXt  (3)                               ! Portion of the linear acceleration of the TMDX in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                  ,SAVE :: LinAccETmdYt  (3)                               ! Portion of the linear acceleration of the TMDY in the inertia frame (body E for earth) associated with everything but the QD2T()'s.  
REAL(ReKi)                  ,SAVE :: FrcOTmdXt(3)                                    ! Portion of the force at the nacelle TmdX reference point   (point TmdXRef   ) due to the TmdX associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcOTmdX(:,:)                                  ! Partial force at the nacelle TmdX reference point   (point TmdXRef   ) due to the TmdX.
REAL(ReKi)                  ,SAVE :: MomNOTmdXt(3)                                   ! Portion of the moment at the nacelle (body N) / TmdX reference point       (point TmdXRef  ) due to the TmdX associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomNOTmdX(:,:)                                 ! Partial moment at the nacelle (body N) / TmdX reference point       (point TmdXRef  ) due the TmdX.                                                                                
REAL(ReKi)                  ,SAVE :: FrcOTmdYt(3)                                    ! Portion of the force at TmdY reference point   (point TmdYRef   ) due to the TmdY associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcOTmdY(:,:)                                  ! Partial force at the TmdY reference point   (point TmdYRef   ) due to the TmdY.
REAL(ReKi)                  ,SAVE :: MomNOTmdYt(3)                                   ! Portion of the moment at the nacelle (body N) / TmdY reference point       (point TmdYRef  ) due to the TmdY associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomNOTmdY(:,:)                                 ! Partial moment at the nacelle (body N) / TmdY reference point       (point TmdYRef  ) due the TmdY.

! add kinematic variables in case tmd is located at the platform and not in the nacelle

REAL(ReKi)                  ,SAVE :: FrcZTmdXt(3)                                    ! Portion of the force at the platform TmdX reference point   (point TmdXRef   ) due to the TmdX associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcZTmdX(:,:)                                  ! Partial force at the platform TmdX reference point   (point TmdXRef   ) due to the TmdX.
REAL(ReKi)                  ,SAVE :: MomXZTmdXt(3)                                   ! Portion of the moment at the platform (body Z) / platform TmdX reference point (point TmdXRef  ) due to the TmdX associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomXZTmdX(:,:)                                 ! Partial moment at the platform (body Z) / platform TmdX reference point (point TmdXRef  ) due the TmdX.                                                                                
REAL(ReKi)                  ,SAVE :: FrcZTmdYt(3)                                    ! Portion of the force at the platform TmdY reference point   (point TmdYRef   ) due to the TmdY associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PFrcZTmdY(:,:)                                  ! Partial force at the platform TmdY reference point   (point TmdYRef   ) due to the TmdY.
REAL(ReKi)                  ,SAVE :: MomXZTmdYt(3)                                   ! Portion of the moment at the platform (body Z) / platform TmdY reference point       (point TmdYRef  ) due to the TmdY associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PMomXZTmdY(:,:)                                 ! Partial moment at the platform (body Z) / platform TmdY reference point       (point TmdYRef  ) due the TmdY.

!mal End of proposed change.  v6.40a-mal  10-Oct-2009.

END MODULE RtHndSid
!=======================================================================
MODULE SimCont


   ! This MODULE stores variables for simulation control.


USE                             Precision
!bjj: these variables should be initialized in an intialization subroutine (for Simulink)


REAL(ReKi)                  ,SAVE :: DT                                              ! Integration time step.
REAL(ReKi)                  ,SAVE :: DT24                                            ! DT/24.
REAL(ReKi)                  ,SAVE :: TMax                                            ! Total run time.
REAL(ReKi)                  ,SAVE :: ZTime    = 0.0                                  ! Current simulation time.

REAL(4)                     ,SAVE :: UsrTime1                                        ! User CPU time for simulation initialization.

INTEGER(4)                  ,SAVE :: Step     = 0                                    ! Current simulation time step.


END MODULE SimCont
!=======================================================================
MODULE TailAero


   ! This MODULE stores input variables for tail fin aerodynamics.


USE                             Precision


REAL(ReKi)                  ,SAVE :: SQRTTFinA = 0.0                                 ! = SQRT( TFinArea )
REAL(ReKi)                  ,SAVE :: TFinAOA   = 0.0                                 ! Angle-of-attack between the relative wind velocity and tail fin chordline
REAL(ReKi)                  ,SAVE :: TFinArea  = 0.0                                 ! Tail fin planform area.
REAL(ReKi)                  ,SAVE :: TFinCD    = 0.0                                 ! Tail fin drag            coefficient resulting from current TFinAOA
REAL(ReKi)                  ,SAVE :: TFinCL    = 0.0                                 ! Tail fin lift            coefficient resulting from current TFinAOA
REAL(ReKi)                  ,SAVE :: TFinCM    = 0.0                                 ! Tail fin pitching moment coefficient resulting from current TFinAOA
REAL(ReKi)                  ,SAVE :: TFinKFx   = 0.0                                 ! Aerodynamic force  at the tail fin center-of-pressure (point K) along tail fin chordline pointing toward tail fin trailing edge (N)
REAL(ReKi)                  ,SAVE :: TFinKFy   = 0.0                                 ! Aerodynamic force  at the tail fin center-of-pressure (point K) normal to plane of tail fin pointing towards suction surface    (N)
REAL(ReKi)                  ,SAVE :: TFinKMz   = 0.0                                 ! Aerodynamic moment at the tail fin center-of-pressure (point K) in plane of tail fin normal to chordline and nominally upward   (N-m)
REAL(ReKi)                  ,SAVE :: TFinQ     = 0.0                                 ! Dynamic pressure of the relative wind velocity

INTEGER(4)                  ,SAVE :: TFinMod   = 0                                   ! Tail fin aerodynamics model switch. (Initialized to zero b/c not all models read in FurlFile)
INTEGER(4)                  ,SAVE :: TFinNFoil = 1                                   ! Tail fin airfoil number. (iniated to first airfoil number)

LOGICAL                     ,SAVE :: SubAxInd  = .FALSE.                             ! Subtract average rotor axial induction when computing relative wind-inflow at tail fin?


END MODULE TailAero
!=======================================================================
MODULE TailFurling


   ! This MODULE stores input variables for tail-furling.


USE                             Precision


REAL(ReKi)                  ,SAVE :: TFrlCDmp  = 0.0                                 ! Tail-furl rate-independent Coulomb-damping moment. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlDmp   = 0.0                                 ! Tail-furl damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlDSDmp = 0.0                                 ! Tail-furl down-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlDSDP  = 0.0                                 ! Tail-furl down-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlDSSP  = 0.0                                 ! Tail-furl down-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlDSSpr = 0.0                                 ! Tail-furl down-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlSpr   = 0.0                                 ! Tail-furl spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlUSDmp = 0.0                                 ! Tail-furl up-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlUSDP  = 0.0                                 ! Tail-furl up-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlUSSP  = 0.0                                 ! Tail-furl up-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlUSSpr = 0.0                                 ! Tail-furl up-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)                  ,SAVE :: TFrlMod   = 0                                   ! Tail-furl spring/damper model switch. (Initialized to zero b/c not all models read in FurlFile)


END MODULE TailFurling
!=======================================================================
MODULE TeeterVars


   ! This MODULE stores input variables for rotor teeter.


USE                             Precision


REAL(ReKi)                  ,SAVE :: TeetCDmp = 0.0                                  ! Rotor-teeter rate-independent Coulomb-damping. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                  ,SAVE :: TeetDmp  = 0.0                                  ! Rotor-teeter damping constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                  ,SAVE :: TeetDmpP = 0.0                                  ! Rotor-teeter damper position. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                  ,SAVE :: TeetHSSp = 0.0                                  ! Rotor-teeter hard-stop linear-spring constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                  ,SAVE :: TeetHStP = 0.0                                  ! Rotor-teeter hard-stop position. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                  ,SAVE :: TeetSSSp = 0.0                                  ! Rotor-teeter soft-stop linear-spring constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                  ,SAVE :: TeetSStP = 0.0                                  ! Rotor-teeter soft-stop position. (initiated to zero b/c the 3-blader requires it to be zero)

INTEGER(4)                  ,SAVE :: TeetMod  = 0                                    ! Rotor-teeter spring/damper model switch. (initiated to zero b/c the 3-blader requires it to be zero)


END MODULE TeeterVars
!=======================================================================
MODULE TipBrakes


   ! This MODULE stores input variables for tip brakes.


USE                             Precision


REAL(ReKi)                  ,SAVE :: TBDrCon                                         ! Instantaneous tip-brake drag constant, Cd*Area.
REAL(ReKi)                  ,SAVE :: TBDrConD                                        ! Tip-brake drag constant during fully-deployed operation, Cd*Area.
REAL(ReKi)                  ,SAVE :: TBDrConN                                        ! Tip-brake drag constant during normal operation, Cd*Area.
REAL(ReKi)                  ,SAVE :: TpBrDT                                          ! Time for tip-brake to reach full deployment once released (sec).


END MODULE TipBrakes
!=======================================================================
MODULE Tower


   ! This MODULE stores variables for the tower.


USE                             Precision


REAL(ReKi)                  ,SAVE :: AdjFASt                                          ! Factor to adjust tower fore-aft stiffness.
REAL(ReKi)                  ,SAVE :: AdjSSSt                                          ! Factor to adjust tower side-to-side stiffness.
REAL(ReKi)                  ,SAVE :: AdjTwMa                                          ! Factor to adjust tower mass density.
REAL(ReKi), ALLOCATABLE     ,SAVE :: AxRedTFA  (:,:,:)                                ! The axial-reduction terms for the fore-aft tower mode shapes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: AxRedTSS  (:,:,:)                                ! The axial-reduction terms for the side-to-side tower mode shapes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: CAT       (:)                                   ! Interpolated, normalized hydrodynamic added mass   coefficient in Morison's equation.
REAL(ReKi), ALLOCATABLE     ,SAVE :: CDT       (:)                                   ! Interpolated, normalized hydrodynamic viscous drag coefficient in Morison's equation.
REAL(ReKi), ALLOCATABLE     ,SAVE :: cgOffTFA  (:)                                    ! Interpolated tower fore-aft mass cg offset.
REAL(ReKi), ALLOCATABLE     ,SAVE :: cgOffTSS  (:)                                    ! Interpolated tower side-to-side mass cg offset.
REAL(ReKi)                  ,SAVE :: CTFA      (2,2)                                  ! Generalized damping of tower in fore-aft direction.
REAL(ReKi)                  ,SAVE :: CTSS      (2,2)                                  ! Generalized damping of tower in side-to-side direction.
REAL(ReKi), ALLOCATABLE     ,SAVE :: DHNodes   (:)                                    ! Length of variable-length tower elements
REAL(ReKi), ALLOCATABLE     ,SAVE :: DiamT     (:)                                   ! Interpolated tower diameter in Morison's equation.
REAL(ReKi)                  ,SAVE :: FAStTunr  (2)                                    ! Tower fore-aft modal stiffness tuners.
REAL(ReKi), ALLOCATABLE     ,SAVE :: HNodes    (:)                                    ! Location of variable-spaced tower nodes (relative to the tower rigid base height)
REAL(ReKi), ALLOCATABLE     ,SAVE :: HNodesNorm(:)                                    ! Normalized location of variable-spaced tower nodes (relative to the tower rigid base height) (0 < HNodesNorm(:) < 1)
REAL(ReKi), ALLOCATABLE     ,SAVE :: HtFract   (:)                                    ! Fractional height of the flexible portion of tower for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: InerTFA   (:)                                    ! Interpolated tower fore-aft (about yt-axis) mass inertia per unit length.
REAL(ReKi), ALLOCATABLE     ,SAVE :: InerTSS   (:)                                    ! Interpolated tower side-to-side (about xt-axis) mass inertia per unit length.
REAL(ReKi)                  ,SAVE :: KTFA      (2,2)                                  ! Generalized stiffness of tower in fore-aft direction.
REAL(ReKi)                  ,SAVE :: KTSS      (2,2)                                  ! Generalized stiffness of tower in side-to-side direction.
REAL(ReKi), ALLOCATABLE     ,SAVE :: MassT     (:)                                    ! Interpolated lineal mass density of tower.
REAL(ReKi)                  ,SAVE :: SSStTunr  (2)                                    ! Tower side-to-side modal stiffness tuners.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffTEA  (:)                                    ! Interpolated tower extensional stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffTFA  (:)                                    ! Interpolated fore-aft tower stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffTGJ  (:)                                    ! Interpolated tower torsional stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: StiffTSS  (:)                                    ! Interpolated side-side tower stiffness.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TMassDen  (:)                                    ! Tower mass density for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwEAStif  (:)                                    ! Tower extensional stiffness for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwFAcgOf  (:)                                    ! Tower fore-aft (along the xt-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwFAIner  (:)                                    ! Tower fore-aft (about yt-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwFAStif  (:)                                    ! Tower fore-aft stiffness for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwGJStif  (:)                                    ! Tower torsional stiffness for a given input station.
REAL(ReKi)                  ,SAVE :: TwrAM     (6,6) = 0.0                           ! Added mass matrix of the current tower element per unit length.
REAL(ReKi)                  ,SAVE :: TwrCA    = 0.0                                  ! Normalized hydrodynamic added mass   coefficient in Morison's equation.
REAL(ReKi)                  ,SAVE :: TwrCD    = 0.0                                  ! Normalized hydrodynamic viscous drag coefficient in Morison's equation.
REAL(ReKi)                  ,SAVE :: TwrDiam  = 0.0                                  ! Tower diameter in Morison's equation.
REAL(ReKi)                  ,SAVE :: TwrFADmp  (2)                                    ! Tower fore-aft structural damping ratios.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwrFASF   (:,:,:)                                ! Tower fore-aft shape functions.
! Yang: Add normal tower modal shape functions for earthquake force calculations
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwrFASF_Eq   (:,:,:)                                ! Tower fore-aft shape functions for seismic analysis.
REAL(ReKi)                  ,SAVE :: TwrFlexL                                         ! Height / length of the flexible portion of the tower.
REAL(ReKi)                  ,SAVE :: TwrFt     (6)   = 0.0                            ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force at the current tower element (point T) and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the tower moment acting at the current tower element (body F) / (point T) per unit length associated with everything but the QD2T()'s.

REAL(ReKi)                  ,SAVE :: TwrSSDmp  (2)                                    ! Tower side-to-side structural damping ratios.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwrSSSF   (:,:,:)                                ! Tower side-to-side shape functions.
! Yang: Add normal tower modal shape functions for earthquake force calculations
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwrSSSF_Eq   (:,:,:)                                ! Tower fore-aft shape functions for seismic analysis.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwSScgOf  (:)                                    ! Tower fore-aft (along the yt-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwSSIner  (:)                                    ! Tower side-to-side (about xt-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TwSSStif  (:)                                    ! Tower side-to-side stiffness for a given input station.
! Yang: Add modal participation factors
REAL(ReKi)                  ,SAVE :: ModalPartFct(4)                                 ! Modal participation factors. [FA1,FA2,SS1,SS2]

INTEGER(4)                  ,SAVE :: NTwInpSt                                        ! Number of tower input stations.
INTEGER(4)                  ,SAVE :: TTopNode                                        ! Index of the additional node located at the tower-top = TwrNodes + 1
INTEGER(4)                  ,SAVE :: TwrLdMod = 0                                    ! Tower loading model switch.
INTEGER(4)                  ,SAVE :: TwrNodes                                        ! Number of tower nodes used in the analysis.


END MODULE Tower
!=======================================================================
MODULE TurbConf


   ! This MODULE stores variables for turbine configuration.


USE                             Precision


REAL(ReKi)                  ,SAVE :: AvgNrmTpRd                                      ! Average tip radius normal to the saft.
REAL(ReKi)                  ,SAVE :: AzimB1Up                                        ! Azimuth value to use for I/O when blade 1 points up.
REAL(ReKi)                  ,SAVE :: BoomCMxn  = 0.0                                 ! Downwind distance from tower-top to tail boom CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: BoomCMyn  = 0.0                                 ! Lateral  distance from tower-top to tail boom CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: BoomCMzn  = 0.0                                 ! Vertical distance from tower-top to tail boom CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: CosDel3   = 1.0                                 ! Cosine of the Delta-3 angle for teetering rotors.
REAL(ReKi), ALLOCATABLE     ,SAVE :: CosPreC   (:)                                   ! Cosines of the precone angles.
REAL(ReKi)                  ,SAVE :: CRFrlSkew = 1.0                                 ! Cosine of the rotor-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: CRFrlSkw2 = 1.0                                 ! Cosine-squared of the rotor-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: CRFrlTilt = 1.0                                 ! Cosine of the rotor-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: CRFrlTlt2 = 1.0                                 ! Cosine-squared of the rotor-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: CShftSkew = 1.0                                 ! Cosine of the shaft skew angle.
REAL(ReKi)                  ,SAVE :: CShftTilt                                       ! Cosine of the shaft tilt angle.
REAL(ReKi)                  ,SAVE :: CSRFrlSkw = 0.0                                 ! Cosine*Sine of the rotor-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: CSRFrlTlt = 0.0                                 ! Cosine*Sine of the rotor-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: CSTFrlSkw = 0.0                                 ! Cosine*Sine of the tail-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: CSTFrlTlt = 0.0                                 ! Cosine*Sine of the tail-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: CTFinBank = 1.0                                 ! Cosine of the tail fin planform  bank angle.
REAL(ReKi)                  ,SAVE :: CTFinSkew = 1.0                                 ! Cosine of the tail fin chordline skew angle.
REAL(ReKi)                  ,SAVE :: CTFinTilt = 1.0                                 ! Cosine of the tail fin chordline tilt angle.
REAL(ReKi)                  ,SAVE :: CTFrlSkew = 1.0                                 ! Cosine of the tail-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: CTFrlSkw2 = 1.0                                 ! Cosine-squared of the tail-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: CTFrlTilt = 1.0                                 ! Cosine of the tail-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: CTFrlTlt2 = 1.0                                 ! Cosine-squared of the tail-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: Delta3    = 0.0                                 ! Delta-3 angle for teetering rotors.
REAL(ReKi)                  ,SAVE :: FASTHH                                          ! Hub-height as computed using FAST inputs [= TowerHt + Twr2Shft + OverHang*SIN( ShftTilt ) ].
REAL(ReKi)                  ,SAVE :: HubCM                                           ! Distance from rotor apex to hub mass.
REAL(ReKi)                  ,SAVE :: HubRad                                          ! Preconed hub radius.
REAL(ReKi)                  ,SAVE :: NacCMxn                                         ! Downwind distance from tower-top to nacelle CM.
REAL(ReKi)                  ,SAVE :: NacCMyn                                         ! Lateral  distance from tower-top to nacelle CM.
REAL(ReKi)                  ,SAVE :: NacCMzn                                         ! Vertical distance from tower-top to nacelle CM.
REAL(ReKi)                  ,SAVE :: OverHang                                        ! Distance from yaw axis to rotor apex or teeter pin.
REAL(ReKi), ALLOCATABLE     ,SAVE :: PreCone   (:)                                   ! Rotor precone angle.
REAL(ReKi)                  ,SAVE :: ProjArea                                        ! Swept area of the rotor projected onto the rotor plane (the plane normal to the low-speed shaft).
REAL(ReKi)                  ,SAVE :: PtfmCM    = 0.0                                 ! Downward distance from the ground [onshore] or MSL [offshore] to the platform CM. (Initialized to zero b/c not all models read in PtfmFile)
REAL(ReKi)                  ,SAVE :: PtfmRef   = 0.0                                 ! Downward distance from the ground [onshore] or MSL [offshore] to the platform reference point. (Initialized to zero b/c not all models read in PtfmFile)
REAL(ReKi)                  ,SAVE :: RefTwrHt                                        ! Vertical distance between FAST's undisplaced tower       height (variable TowerHt) and FAST's inertia frame reference point (variable PtfmRef); that is, RefTwrHt = TowerHt + PtfmRef.
REAL(ReKi)                  ,SAVE :: RFrlCMxn  = 0.0                                 ! Downwind distance from tower-top to rotor-furl CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlCMyn  = 0.0                                 ! Lateral  distance from tower-top to rotor-furl CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlCMzn  = 0.0                                 ! Vertical distance from tower-top to rotor-furl CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlPntxn = 0.0                                 ! Downwind distance from tower-top to arbitrary point on rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlPntyn = 0.0                                 ! Lateral  distance from tower-top to arbitrary point on rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlPntzn = 0.0                                 ! Vertical distance from tower-top to arbitrary point on rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlSkew  = 0.0                                 ! Rotor-furl axis skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: RFrlTilt  = 0.0                                 ! Rotor-furl axis tilt angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: rVDxn     = 0.0                                 ! xn-component of position vector rVD.
REAL(ReKi)                  ,SAVE :: rVDyn     = 0.0                                 ! yn-component of position vector rVD.
REAL(ReKi)                  ,SAVE :: rVDzn     = 0.0                                 ! zn-component of position vector rVD.
REAL(ReKi)                  ,SAVE :: rVIMUxn                                         ! xn-component of position vector rVIMU.
REAL(ReKi)                  ,SAVE :: rVIMUyn                                         ! yn-component of position vector rVIMU.
REAL(ReKi)                  ,SAVE :: rVIMUzn                                         ! zn-component of position vector rVIMU.
REAL(ReKi)                  ,SAVE :: rVPxn     = 0.0                                 ! xn-component of position vector rVP.
REAL(ReKi)                  ,SAVE :: rVPyn     = 0.0                                 ! yn-component of position vector rVP.
REAL(ReKi)                  ,SAVE :: rVPzn                                           ! zn-component of position vector rVP (need not be initialized to zero).
REAL(ReKi)                  ,SAVE :: rWIxn     = 0.0                                 ! xn-component of position vector rWI.
REAL(ReKi)                  ,SAVE :: rWIyn     = 0.0                                 ! yn-component of position vector rWI.
REAL(ReKi)                  ,SAVE :: rWIzn     = 0.0                                 ! zn-component of position vector rWI.
REAL(ReKi)                  ,SAVE :: rWJxn     = 0.0                                 ! xn-component of position vector rWJ.
REAL(ReKi)                  ,SAVE :: rWJyn     = 0.0                                 ! yn-component of position vector rWJ.
REAL(ReKi)                  ,SAVE :: rWJzn     = 0.0                                 ! zn-component of position vector rWJ.
REAL(ReKi)                  ,SAVE :: rWKxn     = 0.0                                 ! xn-component of position vector rWK.
REAL(ReKi)                  ,SAVE :: rWKyn     = 0.0                                 ! yn-component of position vector rWK.
REAL(ReKi)                  ,SAVE :: rWKzn     = 0.0                                 ! zn-component of position vector rWK.
REAL(ReKi)                  ,SAVE :: rZT0zt                                          ! zt-component of position vector rZT0.
REAL(ReKi)                  ,SAVE :: rZYzt     = 0.0                                 ! zt-component of position vector rZY.
REAL(ReKi)                  ,SAVE :: ShftSkew  = 0.0                                 ! Rotor shaft skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: ShftTilt                                        ! Rotor shaft tilt angle.
REAL(ReKi)                  ,SAVE :: SinDel3   = 0.0                                 ! Sine of the Delta-3 angle for teetering rotors.
REAL(ReKi), ALLOCATABLE     ,SAVE :: SinPreC   (:)                                   ! Sines of the precone angles.
REAL(ReKi)                  ,SAVE :: SRFrlSkew = 0.0                                 ! Sine of the rotor-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: SRFrlSkw2 = 0.0                                 ! Sine-squared of the rotor-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: SRFrlTilt = 0.0                                 ! Sine of the rotor-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: SRFrlTlt2 = 0.0                                 ! Sine-squared of the rotor-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: SShftSkew = 0.0                                 ! Sine of the shaft skew angle.
REAL(ReKi)                  ,SAVE :: SShftTilt                                       ! Sine of the shaft tilt angle.
REAL(ReKi)                  ,SAVE :: STFinBank = 0.0                                 ! Sine of the tail fin planform  bank angle.
REAL(ReKi)                  ,SAVE :: STFinSkew = 0.0                                 ! Sine of the tail fin chordline skew angle.
REAL(ReKi)                  ,SAVE :: STFinTilt = 0.0                                 ! Sine of the tail fin chordline tilt angle.
REAL(ReKi)                  ,SAVE :: STFrlSkew = 0.0                                 ! Sine of the tail-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: STFrlSkw2 = 0.0                                 ! Sine-squared of the tail-furl axis skew angle.
REAL(ReKi)                  ,SAVE :: STFrlTilt = 0.0                                 ! Sine of the tail-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: STFrlTlt2 = 0.0                                 ! Sine-squared of the tail-furl axis tilt angle.
REAL(ReKi)                  ,SAVE :: TFrlPntxn = 0.0                                 ! Downwind distance from tower-top to arbitrary point on tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlPntyn = 0.0                                 ! Lateral  distance from tower-top to arbitrary point on tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlPntzn = 0.0                                 ! Vertical distance from tower-top to arbitrary point on tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlSkew  = 0.0                                 ! Rotor-furl axis skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFrlTilt  = 0.0                                 ! Rotor-furl axis tilt angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinBank  = 0.0                                 ! Tail fin planform  bank angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinCMxn  = 0.0                                 ! Downwind distance from tower-top to tail fin CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinCMyn  = 0.0                                 ! Lateral  distance from tower-top to tail fin CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinCMzn  = 0.0                                 ! Vertical distance from tower-top to tail fin CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinCPxn  = 0.0                                 ! Downwind distance from tower-top to tail fin CP. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinCPyn  = 0.0                                 ! Lateral  distance from tower-top to tail fin CP. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinCPzn  = 0.0                                 ! Vertical distance from tower-top to tail fin CP. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinSkew  = 0.0                                 ! Tail fin chordline skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TFinTilt  = 0.0                                 ! Tail fin chordline tilt angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                  ,SAVE :: TipRad                                          ! Preconed blade-tip radius.
REAL(ReKi)                  ,SAVE :: TowerHt                                         ! Height of tower above ground level.
REAL(ReKi)                  ,SAVE :: Twr2Shft                                        ! Vertical distance from the tower-top to the rotor shaft.
REAL(ReKi)                  ,SAVE :: TwrDraft  = 0.0                                 ! Downward distance from the ground [onshore] or MSL [offshore] to the tower base platform connection. (Initialized to zero b/c not all models read in PtfmFile)
REAL(ReKi)                  ,SAVE :: TwrRBHt                                         ! Tower rigid base height.
REAL(ReKi)                  ,SAVE :: UndSling                                        ! Undersling length.
REAL(ReKi)                  ,SAVE :: Yaw2Shft  = 0.0                                 ! Lateral distance from the yaw axis to the rotor shaft. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)                  ,SAVE :: NumBl                                           ! Number of blades.
INTEGER(4)                  ,SAVE :: PSpnElN   = 1                                   ! Number of the innermost blade element which is still part of the pitchable portion of the blade for partial-span pitch control.


END MODULE TurbConf
!=======================================================================
MODULE TurbCont


   ! This MODULE stores input variables for turbine control.


USE                             Precision


REAL(ReKi), ALLOCATABLE     ,SAVE :: BlPitch  (:)                                    ! Initial and current blade pitch angles.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BlPitchCom(:)                                   ! Commanded blade pitch angles.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BlPitchF (:)                                    ! Final blade pitch.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BlPitchFrct(:)                                  ! Blade pitch angle fractions used for the override pitch maneuver calculation.
REAL(ReKi), ALLOCATABLE     ,SAVE :: BlPitchI (:)                                    ! Initial blade pitch angles at the start of the override pitch maneuver.
REAL(ReKi)                  ,SAVE :: NacYawF                                         ! Final yaw angle after override yaw maneuver.
REAL(ReKi)                  ,SAVE :: SpdGenOn                                        ! Generator speed to turn on the generator for a startup.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TBDepISp (:)                                    ! Deployment-initiation speed for the tip brakes.
REAL(ReKi)                  ,SAVE :: THSSBrDp                                        ! Time to initiate deployment of the shaft brake.
REAL(ReKi)                  ,SAVE :: THSSBrFl                                        ! Time at which shaft brake is fully deployed.
REAL(ReKi)                  ,SAVE :: TiDynBrk                                        ! Time to initiate deployment of the dynamic generator brake.
REAL(ReKi)                  ,SAVE :: TimGenOf                                        ! Time to turn off generator for braking or modeling a run-away.
REAL(ReKi)                  ,SAVE :: TimGenOn                                        ! Time to turn on generator for startup.
REAL(ReKi)                  ,SAVE :: TPCOn                                           ! Time to enable active pitch control.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TPitManE (:)                                    ! Time to end pitch maneuvers for each blade.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TPitManS (:)                                    ! Time to start pitch maneuvers for each blade.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TTpBrDp  (:)                                    ! Times to initiate deployment of tip brakes.
REAL(ReKi), ALLOCATABLE     ,SAVE :: TTpBrFl  (:)                                    ! Times at which tip brakes are fully deployed.
REAL(ReKi)                  ,SAVE :: TYawManE                                        ! Time to end override yaw maneuver.
REAL(ReKi)                  ,SAVE :: TYawManS                                        ! Time to start override yaw maneuver.
REAL(ReKi)                  ,SAVE :: TYCOn                                           ! Time to enable active yaw control.
REAL(ReKi)                  ,SAVE :: VS_Rgn2K                                        ! Generator torque constant in Region 2 (HSS side), N-m/rpm^2.
REAL(ReKi)                  ,SAVE :: VS_RtGnSp                                       ! Rated generator speed (HSS side), rpm.
REAL(ReKi)                  ,SAVE :: VS_RtTq                                         ! Rated generator torque/constant generator torque in Region 3 (HSS side), N-m.
REAL(ReKi)                  ,SAVE :: VS_Slope                                        ! Torque/speed slope of region 2 1/2 induction generator.
REAL(ReKi)                  ,SAVE :: VS_SlPc                                         ! Rated generator slip percentage in Region 2 1/2, %.
REAL(ReKi)                  ,SAVE :: VS_SySp                                         ! Synchronous speed of region 2 1/2 induction generator.
REAL(ReKi)                  ,SAVE :: VS_TrGnSp                                       ! Transitional generator speed between regions 2 and 2 1/2.
REAL(ReKi)                  ,SAVE :: YawPosCom                                       ! Commanded yaw angle from user-defined routines, rad.
REAL(ReKi)                  ,SAVE :: YawRateCom                                      ! Commanded yaw rate  from user-defined routines, rad/s.

INTEGER(4)                  ,SAVE :: GenModel                                        ! Generator model
INTEGER(4)                  ,SAVE :: HSSBrMode                                       ! HSS brake model.
INTEGER(4)                  ,SAVE :: PCMode                                          ! Pitch control mode
INTEGER(4)                  ,SAVE :: VSContrl                                        ! Variable-speed-generator control switch.
INTEGER(4)                  ,SAVE :: YCMode                                          ! Yaw control mode

LOGICAL,    ALLOCATABLE     ,SAVE :: BegPitMan(:)                                    ! .TRUE. before the override pitch manuever has begun (begin pitch manuever).
LOGICAL                     ,SAVE :: GenTiStp                                        ! Stop generator based upon T: time or F: generator power = 0.
LOGICAL                     ,SAVE :: GenTiStr                                        ! Start generator based upon T: time or F: generator speed.


END MODULE TurbCont
!=======================================================================

! Yang: Add the TMDControl module. 14-Aug-2018

MODULE TMD


   ! This MODULE stores input variables for TMDs.

USE                             NWTC_Library
USE                             Precision

! Variables for active control mode
INTEGER(4)                  ,SAVE :: TmdXCMode                                       ! TmdX control mode
REAL(ReKi)                  ,SAVE :: TTmdXCOn                                        ! time TMDX control starts
INTEGER(4)                  ,SAVE :: TmdYCMode                                       ! TmdY control mode
REAL(ReKi)                  ,SAVE :: TTmdYCOn                                        ! time TMDY control starts
REAL(ReKi)                  ,SAVE :: TmdXSprCom                                      ! Commanded TmdX spring stiffness.
REAL(ReKi)                  ,SAVE :: TmdXDampCom                                     ! Commanded TmdX damping.
REAL(ReKi)                  ,SAVE :: TmdXFextCom                                      ! Commanded TmdX external force from user-defined routines, m.
REAL(ReKi)                  ,SAVE :: TmdYSprCom                                      ! Commanded TmdX spring stiffness.
REAL(ReKi)                  ,SAVE :: TmdYDampCom                                     ! Commanded TmdX damping.
REAL(ReKi)                  ,SAVE :: TmdYFextCom                                      ! Commanded TmdX external force from user-defined routines, m.

! Variables for Tuned mass damper in X drection (fore-aft)
LOGICAL                     ,SAVE :: TmdXDOF   = .FALSE.                             ! Axial TMD DOF. (Initialize to .FALSE. b/c not all models will use a TMD)
INTEGER(4)                  ,SAVE :: TmdXLoc                                         ! TmdX location {1: nacelle, 2: platform} 
REAL(ReKi)                  ,SAVE :: TmdXRefxnt                                      ! Downwind distance from tower-top to TmdX CM.
REAL(ReKi)                  ,SAVE :: TmdXRefynt                                      ! Lateral  distance from tower-top to TmdX CM.
REAL(ReKi)                  ,SAVE :: TmdXRefznt                                      ! Vertical distance from tower-top to TmdX CM.
REAL(ReKi)                  ,SAVE :: TmdXAngle                                       ! TmdX rotation angle about vertical axis. 
REAL(ReKi)                  ,SAVE :: TmdXMass   = 0.0                                ! TmdX mass. (Initialized to zero b/c not all models will use a Tmd)
REAL(ReKi)                  ,SAVE :: TmdXSpr                                         ! Initial and current TmdX spring stiffness.
REAL(ReKi)                  ,SAVE :: TmdXDamp                                        ! Initial and current TmdX damping.
REAL(ReKi)                  ,SAVE :: TmdXFext                                        ! Commanded TmdX external force from user-defined routines, m
REAL(ReKi)                  ,SAVE :: TmdXDWSP                                        ! TmdX stop downwind stop position. 
REAL(ReKi)                  ,SAVE :: TmdXUWSP                                        ! TmdX stop upwind stop position.
REAL(ReKi)                  ,SAVE :: TmdXSDamp                                       ! TmdX stop damping constant. 
REAL(ReKi)                  ,SAVE :: TmdXSSpr                                        ! TmdX stop spring constant.
REAL(ReKi)                  ,SAVE :: TmdXNeut                                        ! Neutral TmdX position.
! Variables for Tuned mass damper in Y drection (fore-aft)
LOGICAL                     ,SAVE :: TmdYDOF   = .FALSE.                             ! Transverse TMD DOF. (Initialize to .FALSE. b/c not all models will use a TMD)
INTEGER(4)                  ,SAVE :: TmdYLoc                                         ! TmdY location {1: nacelle, 2: platform} 
REAL(ReKi)                  ,SAVE :: TmdYRefxnt                                      ! Downwind distance from tower-top to TmdY CM.
REAL(ReKi)                  ,SAVE :: TmdYRefynt                                      ! Lateral  distance from tower-top to TmdY CM.
REAL(ReKi)                  ,SAVE :: TmdYRefznt                                      ! Vertical distance from tower-top to TmdY CM.
REAL(ReKi)                  ,SAVE :: TmdYAngle                                       ! TmdY rotation angle about vertical axis. 
REAL(ReKi)                  ,SAVE :: TmdYMass   = 0.0                                ! TmdY mass. (Initialized to zero b/c not all models will use a Tmd)
REAL(ReKi)                  ,SAVE :: TmdYSpr                                         ! Initial and current TmdY spring stiffness.
REAL(ReKi)                  ,SAVE :: TmdYDamp                                        ! Initial and current TmdY damping.
REAL(ReKi)                  ,SAVE :: TmdYFext                                        ! Commanded TmdY external force from user-defined routines, m.
REAL(ReKi)                  ,SAVE :: TmdYSDamp                                       ! TmdY stop damping constant. 
REAL(ReKi)                  ,SAVE :: TmdYSSpr                                        ! TmdY stop spring constant. 
REAL(ReKi)                  ,SAVE :: TmdYPLSP                                        ! TmdY stop positive-lateral stop position. 
REAL(ReKi)                  ,SAVE :: TmdYNLSP                                        ! TmdY stop negative-lateral stop position. 
REAL(ReKi)                  ,SAVE :: TmdYNeut                                        ! Neutral TmdY position.


CONTAINS
!=======================================================================
! subroutine ReadTMDFile
! Read all variables defined in the TMD control file 
!=======================================================================
SUBROUTINE ReadTMDFile()
USE                                 Precision
USE                                 General,ONLY:UnSc,StrcCtrlFile,Cmpl4SFun
USE                                 NWTC_IO
USE                                 NWTC_Library

INTEGER(4)                 ,SAVE :: i
 ! Open the TMD control input file:
    CALL OpenFInpFile(UnSc, StrcCtrlFile )
    CALL ReadCom ( UnSc, StrcCtrlFile, 'comments1' )
    CALL ReadCom ( UnSc, StrcCtrlFile, 'comments2' )
    CALL ReadCom ( UnSc, StrcCtrlFile, 'comments3' )
    CALL ReadCom ( UnSc, StrcCtrlFile, 'Active control mode' )
    CALL ReadIVar ( UnSc, StrcCtrlFile, TmdXCMode, 'TmdXCMode', 'TmdX control mode' )
    IF ( ( .NOT. Cmpl4SFun ) .AND. ( TmdXCMode == 2 ) )  THEN
       CALL Abort ( ' TmdXCMode can only equal 2 when FAST is interfaced with Simulink.'// &
                '  Set TmdXCMode to 0 or 1 or interface FAST with Simulink.'             )
    ELSEIF ( ( TmdXCMode < 0 ) .OR. ( TmdXCMode > 2 ) )  THEN
           CALL Abort ( ' TmdXCMode must be 0, 1, or 2.' )
    ENDIF
    
    CALL ReadRVar(UnSc, StrcCtrlFile, TTmdXCOn, 'TTmdXCOn', 'Time to enable active TmdX control (s)' )
    
    IF ( Cmpl4SFun .AND. ( TmdXCMode == 2 ) .AND. ( TTmdXCOn /= 0.0 ) )  THEN
       CALL Abort ( ' TmdX control must be enabled at time zero when implemented in Simulink.'//      &
                    '  Set TTmdXCOn to 0.0, set TmdXCMode to 0 or 1, or use the standard version of FAST.'   )
    ELSEIF ( TTmdXCOn < 0.0 )  THEN
       CALL Abort ( ' TTmdXCOn must not be negative.' )
    ENDIF
    
    CALL ReadIVar ( UnSc, StrcCtrlFile, TmdYCMode, 'TmdYCMode', 'TmdY control mode' )
    IF ( ( .NOT. Cmpl4SFun ) .AND. ( TmdYCMode == 2 ) )  THEN
       CALL Abort ( ' TmdYCMode can only equal 2 when FAST is interfaced with Simulink.'// &
                '  Set TmdYCMode to 0 or 1 or interface FAST with Simulink.'             )
    ELSEIF ( ( TmdYCMode < 0 ) .OR. ( TmdYCMode > 2 ) )  THEN
           CALL Abort ( ' TmdYCMode must be 0, 1, or 2.' )
    ENDIF
    CALL ReadRVar(UnSc, StrcCtrlFile, TTmdYCOn, 'TTmdYCOn', 'Time to enable active TmdY control (s)' )
    IF ( Cmpl4SFun .AND. ( TmdYCMode == 2 ) .AND. ( TTmdYCOn /= 0.0 ) )  THEN
       CALL Abort ( ' TmdY control must be enabled at time zero when implemented in Simulink.'//      &
                    '  Set TTmdYCOn to 0.0, set TmdYCMode to 0 or 1, or use the standard version of FAST.'   )
    ELSEIF ( TTmdYCOn < 0.0 )  THEN
       CALL Abort ( ' TTmdYCOn must not be negative.' )
    ENDIF
    
 ! read the configurations of TMD_X   
    CALL ReadCom ( UnSc, StrcCtrlFile, 'Tuned mass damper in X drection (fore-aft)' )
    CALL ReadLVar(UnSc, StrcCtrlFile, TmdXDOF, 'TmdXDOF', 'TmdX DOF (flag)' )
    CALL ReadIVar ( UnSc, StrcCtrlFile, TmdXLoc, 'TmdXLoc', 'Location of the TmdX' )
    IF ( ( TmdXLoc < 1 ) .OR. ( TmdXLoc > 2 ) )  THEN
       CALL Abort ( ' TmdXLoc must be 1 or 2.' )
    ENDIF
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXRefxnt, 'TmdXRefxnt', 'Downwind distance from the tower-top (n) or platform reference (t) to the TmdX' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXRefynt, 'TmdXRefynt', 'Lateral  distance from the tower-top (n) or platform reference (t) to the TmdX' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXRefznt, 'TmdXRefznt', 'Vertical distance from the tower-top (n) or platform reference (t) to the TmdX' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXAngle, 'TmdXAngle', 'Initial or fixed TmdX coordinate system rotation angle about vertical axis (degrees)' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXMass, 'TmdXMass', 'TmdX mass' )
    IF (( TmdXMass == 0.0 ) .AND. (TmdXDOF))  CALL Abort ( ' TmdXMass must not be zero if the TmdXDOF is True.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXSpr, 'TmdXSpr', 'TmdX initial or fixed spring stiffness' )
    IF (( TmdXSpr < 0.0 ) .AND. (TmdXDOF))  CALL ProgAbort ( ' TmdXSpr must not be negative.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXDamp, 'TmdXDamp', 'TmdX initial or fixed damping' )
    IF (( TmdXDamp < 0.0 ) .AND. (TmdXDOF))  CALL ProgAbort ( ' TmdXDamp must not be negative.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXFext, 'TmdXFext', 'TmdX initial or fixed external force' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXDWSP, 'TmdXDWSP', 'TmdX downwind distance to downwind stop position' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXUWSP, 'TmdXUWSP', 'TmdX downwind distance to upwind stop position' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXSSpr, 'TmdXSSpr', 'TmdX initial or fixed stop spring stiffness' )
    IF ( TmdXSSpr < 0.0 )  CALL ProgAbort ( ' TmdXSSpr must not be negative.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdXSDamp, 'TmdXSDamp', 'TmdX initial or fixed stop damping' )
    IF ( TmdXSDamp < 0.0 )  CALL ProgAbort ( ' TmdXSDamp must not be negative.' )
 ! read the configurations of TMD_Y    
    CALL ReadCom (UnSc, StrcCtrlFile, 'Tuned mass damper in Y drection (side-side)' )
    CALL ReadLVar(UnSc, StrcCtrlFile, TmdYDOF, 'TmdYDOF', 'TmdY DOF (flag)' )
    CALL ReadIVar (UnSc, StrcCtrlFile, TmdYLoc, 'TmdYLoc', 'Location of the TmdY' )
    IF ( ( TmdYLoc < 1 ) .OR. ( TmdYLoc > 2 ) )  THEN
       CALL Abort ( ' TmdYLoc must be 1 or 2.' )
    ENDIF
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYRefxnt, 'TmdYRefxnt', 'Downwind distance from the tower-top (n) or platform reference (t) to the TmdY' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYRefynt, 'TmdYRefynt', 'Lateral  distance from the tower-top (n) or platform reference (t) to the TmdY' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYRefznt, 'TmdYRefznt', 'Vertical distance from the tower-top (n) or platform reference (t) to the TmdY' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYAngle, 'TmdYAngle', 'Initial or fixed TmdY coordinate system rotation angle about vertical axis (degrees)' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYMass, 'TmdYMass', 'TmdY mass' )
    IF (( TmdYMass == 0.0 ) .AND. (TmdYDOF))  CALL Abort ( ' TmdYMass must not be zero if the TmdYDOF is True.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYSpr, 'TmdYSpr', 'TmdY initial or fixed spring stiffness' )
    IF (( TmdYSpr < 0.0 ) .AND. (TmdYDOF))  CALL ProgAbort ( ' TmdYSpr must not be negative.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYDamp, 'TmdYDamp', 'TmdY initial or fixed damping' )
    IF (( TmdYDamp < 0.0 ) .AND. (TmdYDOF))  CALL ProgAbort ( ' TmdYDamp must not be negative.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYFext, 'TmdYFext', 'TmdY initial or fixed external force' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYDWSP, 'TmdYDWSP', 'TmdY downwind distance to downwind stop position' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYUWSP, 'TmdYUWSP', 'TmdY downwind distance to upwind stop position' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYSSpr, 'TmdYSSpr', 'TmdY initial or fixed stop spring stiffness' )
    IF ( TmdYSSpr < 0.0 )  CALL ProgAbort ( ' TmdYSSpr must not be negative.' )
    CALL ReadRVar(UnSc, StrcCtrlFile, TmdYSDamp, 'TmdYSDamp', 'TmdY initial or fixed stop damping' )
    IF ( TmdYSDamp < 0.0 )  CALL ProgAbort ( ' TmdYSDamp must not be negative.' )
    
    CLOSE ( UnSc )
END SUBROUTINE ReadTMDFile    

END MODULE TMD

! Yang
MODULE AQWAinFAST


   ! This MODULE stores the variables passing through the DLL of AQWA


USE                             Precision

REAL (Reki)          ,SAVE :: PtfmMotAQWA(6)                    ! Platform motion at the current time step
REAL (Reki)          ,SAVE :: PtfmVelAQWA(6)                    ! Platform velocity at the current time step
REAL (Reki)          ,SAVE :: PtfmAcceAQWA(6)                   ! Platform acceleration at the current time step
REAL (Reki)          ,SAVE :: PtfmMomAQWA(3)                    ! Platform moment due to tower base forces 
REAL (Reki)          ,SAVE :: PtfmMomAQWA_Yaw                   ! Platform moment due to tower base forces 
REAL (Reki)          ,SAVE :: PtfmFrcAQWA(3)                    ! Platform forces due to tower base forces 
LOGICAL              ,SAVE :: FAST2AQWA                         ! Coupling flag between FAST & AQWA 
LOGICAL              ,SAVE :: PrtFASTRslts                      ! Flag of outputting the FAST results at the current time step. .TRUE. for STAGE = 1 and .FALSE. for STAGE = 2.

END MODULE AQWAinFAST