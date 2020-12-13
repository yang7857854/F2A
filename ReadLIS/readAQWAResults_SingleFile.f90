!  readAQWAResults_SingleFile.f90 
!
!  FUNCTIONS:
!  readAQWAResults_SingleFile - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: readAQWAResults_SingleFile
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************


    program readAQWAResults_SingleFile
    USE IFPORT

    implicit none
! Variables
    character (len=1024) :: dirBase,LISfileBase,OutFolder
    real :: TimeDur,TimeStep
    integer :: StartCase,EndCase,NumBodies,NumAddedPts,IDBodyAddedPts
    character(len=1024),allocatable ::NameBodies(:)
    logical,allocatable :: FlagOPMot(:),FlagOPVel(:),FlagOPAcc(:),FlagOPFrc(:)
    integer,allocatable :: IDNodes(:)
    integer :: TdonCncID(2),NumTdon,MoorLineID,NumMoorLines
    logical :: FlagMLType,FlagTdMlTurn,FlagWvElv
    integer :: i,j,k,lenWave
    integer,allocatable :: UnFiles(:)
! Variables in reading
    integer :: icase,LTID
    character(len=10) :: icaseChar
    character(len=1024) :: charTendon,LineStr,charMoor,CharWvEle
    character(len=13) :: tempUnit
    character(len=5) :: charAddedPts
    integer :: res,FEXIST
    character(len=1024) :: Outbase,InputLISFile,charUnit,TempChar
    character (len=41),allocatable :: charNodeOP(:)
    character(len=13) :: CharTension
    integer :: TdonLTID,TdonUTID,status,NumSteps,iStep,iStruct,tempC
    real,allocatable :: Disp(:,:),Vel(:,:),Acce(:,:),Force(:,:),TdTensionLT(:),TdTensionUT(:),AddedPtAcce(:,:),MoorTensionFL(:),MoorTensionAC(:),AddedPtVel(:,:),AddedPtDisp(:,:)
    real :: WaveVel(3),WaveAcce(3)
    real :: SimTime,WvEle
    logical :: wrtFlag
! Time and date information that will be output in the result file.
    character(len=8) :: Currdate
    character(len=10) :: CurrTime
    character(len=50) :: Time_Date


! Open the input file to read the parameters
    open (1, file = 'InputforreadAQWAResults_SingleFile.txt')
     read (1,*) 
     read (1,*) 
     read (1,*)  ! --------- Global parameters -----------------------
     read (1,*)  InputLISFile 		! - Input LIS file
     read (1,*)  OutFolder		! - The target folder to store the result files			
     read (1,*)  TimeDur		! - Duration of each simulation
     read (1,*)  TimeStep		! - Time step of each simulation
     read (1,*)  ! ---------- Details of the LIS file -----------------
     read (1,*)  NumBodies		!- Number of bodies defined in the LIS file.
     allocate (NameBodies(NumBodies))
     allocate (FlagOPMot(NumBodies))
     allocate (FlagOPVel(NumBodies))
     allocate (FlagOPAcc(NumBodies))
     allocate (FlagOPFrc(NumBodies))
     do i=1,NumBodies
        read (1,*)  !----------- Body i ----------------------------------
        read (1,*)  NameBodies(i)		! - Name of the 1st body
        read (1,*)  FlagOPMot(i)		! - Flag of outputting motions at the 1st body's CoG
        read (1,*)  FlagOPVel(i)		! - Flag of outputting velocities at the 1st body's CoG
        read (1,*)  FlagOPAcc(i)		! - Flag of outputting accelerations at the 1st body's CoG
        read (1,*)  FlagOPFrc(i)		! - Flag of outputting forces at the 1st body's CoG
     enddo
     read (1,*) !---------- Added points -----------------------------
     read (1,*) NumAddedPts
     if (NumAddedPts>0) then
        allocate (IDNodes(NumAddedPts))
        do i=1,NumAddedPts
            read (1,*) IDNodes(i)
        enddo
     endif
     read (1,*) ! --------- Tension of cables -----------------------
     read (1,*) FlagMLType                  ! - Tendon type (True: nonlinear catenary, False: linear cable )
     read (1,*) NumMoorLines		        ! - Number of mooring lines
     read (1,*) ! --------- Wave elevation -----------------------
     read (1,*) FlagWvElv		            ! - Output wave elevation
     read (1,*) CharWvEle		            ! - Characteristic string of the wave elevation line.
     close(1)

    lenWave = len(trim(adjustl(CharWvEle)))
! allocate output ID string
    allocate (charNodeOP(NumAddedPts))
    if (NumAddedPts>0) then
       do i=1,NumAddedPts
          write(charNodeOP(i),'(I6)') IDNodes(i)  
          charNodeOP(i) = '                      POSITION NODE'//charNodeOP(i)
       enddo
    endif

! Time and date of generating the current file.
    call date_and_time(DATE=Currdate,Time=CurrTime)
    Time_Date = ' at '//CurrTime(1:2)//':'//CurrTime(3:4)//':'//CurrTime(5:6)
    Time_Date = trim(Time_Date)//' on '//Currdate(7:8)//'-'//Currdate(5:6)//'-'//Currdate(1:4)
! Allocate the file IO matrix (UnFiles)
    allocate (UnFiles(NumBodies*4))
    do i=1,NumBodies*4
       UnFiles (i) = 1000 + i
    enddo
! Creat the folder
    inquire(DIRECTORY=trim(OutFolder),EXIST=FEXIST)
    if(.not.FEXIST) then
       res = MAKEDIRQQ(trim(OutFolder))
    endif
! Read the results in loop
       Outbase      = trim(OutFolder)//'\'//InputLISFile(1:len(trim(InputLISFile))-4)
       ! Creat the files
       do i=1,NumBodies
          ! motion file
          if (FlagOPMot(i)) then
             open (UnFiles((i-1)*4+1),file =trim(Outbase)//'_'//trim(NameBodies(i))//'_Motion.dat' )  
              write(UnFiles((i-1)*4+1),'(A)') 'Motion of the '//trim(NameBodies(i))//' of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
              write(UnFiles((i-1)*4+1),'(A10,6(A13))')'Time','Surge','Sway','Heave','Roll','Pitch','Yaw'
              write(UnFiles((i-1)*4+1),'(A10,6(A13))') '(s)','(m)','(m)','(m)','(deg)','(deg)','(deg)'
          endif
          ! Velocity file
          if (FlagOPVel(i)) then
             open (UnFiles((i-1)*4+2),file =trim(Outbase)//'_'//trim(NameBodies(i))//'_Velocity.dat' ) 
             write(UnFiles((i-1)*4+2),'(A)') 'Velocity of the '//trim(NameBodies(i))//' of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
             write(UnFiles((i-1)*4+2),'(A10,6(A13))')'Time','SurgeVel','SwayVel','HeaveVel','RollVel','PitchVel','YawVel'
             write(UnFiles((i-1)*4+2),'(A10,6(A13))') '(s)','(m/s)','(m/s)','(m/s)','(deg/s)','(deg/s)','(deg/s)' 
          endif
          ! Acceleration file
          if (FlagOPAcc(i)) then
             open (UnFiles((i-1)*4+3),file =trim(Outbase)//'_'//trim(NameBodies(i))//'_Acceleration.dat' ) 
             write(UnFiles((i-1)*4+3),'(A)') 'Acceleration of the '//trim(NameBodies(i))//' of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
             write(UnFiles((i-1)*4+3),'(A10,6(A13))')'Time','SurgeAcce','SwayAcce','HeaveAcce','RollAcce','PitchAcce','YawAcce'
             write(UnFiles((i-1)*4+3),'(A10,6(A13))') '(s)','(m/s^2)','(m/s^2)','(m/s^2)','(deg/s^2)','(deg/s^2)','(deg/s^2)'
          endif
          ! Force file
          if (FlagOPFrc(i)) then
             open (UnFiles((i-1)*4+4),file =trim(Outbase)//'_'//trim(NameBodies(i))//'_Force.dat' )
             write(UnFiles((i-1)*4+4),'(A)') 'Overall forces applied at the '//trim(NameBodies(i))//' of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
             write(UnFiles((i-1)*4+4),'(A10,6(A13))')'Time','Force_x','Force_y','Force_z','Moment_x','Moment_y','Moment_z'
             write(UnFiles((i-1)*4+4),'(A10,6(A13))') '(s)','(N)','(N)','(N)','(N-m)','(N-m)','(N-m)' 
          endif
       enddo
      ! Results of added points 
      if (NumAddedPts>0) then
          do j=1,NumAddedPts
              write(charAddedPts,'(I5)') j
              open (UnFiles(4*NumBodies)+(j-1)*3+1,file = trim(Outbase)//'_'//'AddedPt_'//trim(adjustl(charAddedPts))//'_Acceleration.dat')
               write(UnFiles(4*NumBodies)+(j-1)*3+1,'(A)') 'Acclerations of the added point '//trim(adjustl(charAddedPts))//' of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
               write(UnFiles(4*NumBodies)+(j-1)*3+1,'(A10,3(A13))')'Time','x_Acce','y_Acce','z_Acce'
               write(UnFiles(4*NumBodies)+(j-1)*3+1,'(A10,3(A13))') '(s)','(m/s^2)','(m/s^2)','(m/s^2)'
              open (UnFiles(4*NumBodies)+(j-1)*3+2,file = trim(Outbase)//'_'//'AddedPt_'//trim(adjustl(charAddedPts))//'_Velocity.dat')
               write(UnFiles(4*NumBodies)+(j-1)*3+2,'(A)') 'Velocity of the added point '//trim(adjustl(charAddedPts))//' of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
               write(UnFiles(4*NumBodies)+(j-1)*3+2,'(A10,3(A13))')'Time','x_Vel','y_Vel','z_Vel'
               write(UnFiles(4*NumBodies)+(j-1)*3+2,'(A10,3(A13))') '(s)','(m/s)','(m/s)','(m/s)'
              open (UnFiles(4*NumBodies)+(j-1)*3+3,file = trim(Outbase)//'_'//'AddedPt_'//trim(adjustl(charAddedPts))//'_Displacement.dat')
               write(UnFiles(4*NumBodies)+(j-1)*3+3,'(A)') 'Displacement of the added point '//trim(adjustl(charAddedPts))//' of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
               write(UnFiles(4*NumBodies)+(j-1)*3+3,'(A10,3(A13))')'Time','x_Disp','y_Disp','z_Disp'
               write(UnFiles(4*NumBodies)+(j-1)*3+3,'(A10,3(A13))') '(s)','(m)','(m)','(m)'
          enddo
      endif
      ! Results of the cable tensions
      ! Mooring lines
      charUnit = ''
      charMoor = ''
      do j=1,NumMoorLines
         write(tempUnit,'(I3)') j
         charMoor = trim(charMoor)//'  Mooring'//trim(adjustl(tempUnit))//'_FL'
         write(tempUnit,'(A13)') '(N)'
         charUnit = trim(charUnit)//trim(tempUnit)
      enddo
      do j=1,NumMoorLines
         write(tempUnit,'(I3)') j
         charMoor = trim(charMoor)//'  Mooring'//trim(adjustl(tempUnit))//'_AC'
         write(tempUnit,'(A13)') '(N)'
         charUnit = trim(charUnit)//trim(tempUnit)
      enddo

     ! Mooring lines    
      open (UnFiles(4*NumBodies)+NumAddedPts*3+3, file =trim(Outbase)//'_'//'MooringTension.dat')
       write (UnFiles(4*NumBodies)+NumAddedPts*3+3,'(A)') 'Tension of mooring lines of the '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
       write (UnFiles(4*NumBodies)+NumAddedPts*3+3,'(A10,A)') 'Time',trim(charMoor)  
       write (UnFiles(4*NumBodies)+NumAddedPts*3+3,'(A10,A)') 's',trim(charUnit)
       
    ! Wave elevation
      open (UnFiles(4*NumBodies)+NumAddedPts*3+4, file =trim(Outbase)//'_'//'WaveElevation.dat')
       write (UnFiles(4*NumBodies)+NumAddedPts*3+4,'(A)') 'Wave elevation '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
       write (UnFiles(4*NumBodies)+NumAddedPts*3+4,'(A10,A13)') 'Time','WvEle'
       write (UnFiles(4*NumBodies)+NumAddedPts*3+4,'(A10,A13)')  '(s)','(m)'
      open (UnFiles(4*NumBodies)+NumAddedPts*3+5, file =trim(Outbase)//'_'//'WaveVelocity.dat')
       write (UnFiles(4*NumBodies)+NumAddedPts*3+5,'(A)') 'Wave elevation '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
       write (UnFiles(4*NumBodies)+NumAddedPts*3+5,'(A10,A13)') 'Time','WvEle'
       write (UnFiles(4*NumBodies)+NumAddedPts*3+5,'(A10,A13)')  '(s)','(m)'
      open (UnFiles(4*NumBodies)+NumAddedPts*3+6, file =trim(Outbase)//'_'//'WaveAcceleration.dat')
       write (UnFiles(4*NumBodies)+NumAddedPts*3+6,'(A)') 'Wave elevation '//trim(Outbase)//' case extracted from AQWA results using the code developed by Y.Yang (PDRA in LJMU)'//trim(Time_Date)
       write (UnFiles(4*NumBodies)+NumAddedPts*3+6,'(A10,A13)') 'Time','WvEle'
       write (UnFiles(4*NumBodies)+NumAddedPts*3+6,'(A10,A13)')  '(s)','(m)'
    
    
! Read the results
      NumSteps = ceiling(TimeDur/TimeStep)+1
! Allocate matrix
      allocate (Disp(NumBodies,6))
      allocate (Vel(NumBodies,6))
      allocate (Acce(NumBodies,6))
      allocate (Force(NumBodies,6))
      if (NumAddedPts > 0) then
         allocate (AddedPtAcce(NumAddedPts,3))
         allocate (AddedPtVel(NumAddedPts,3))
         allocate (AddedPtDisp(NumAddedPts,3))
      endif
      allocate (MoorTensionFL(NumMoorLines))
      allocate (MoorTensionAC(NumMoorLines))

      iStep = 1
      SimTime = 0
      wrtFlag = .false.
      write(*,*) 'InputLISFile= ',trim(InputLISFile)
      open(1,file = trim(InputLISFile))
       iStruct = 0
       do    ! time step
          read(1,'(A1024)',IOSTAT = status) LineStr
          if (status/=0) exit
          if (LineStr(1:11) == ' TIME(SECS)') then   ! time step
             iStruct = iStruct + 1
             !wrtFlag = .true.
             if (iStruct>NumBodies) then
                ! write results
          
             ! CoG results
                do i=1,NumBodies
                   ! motion file
                   if (FlagOPMot(i)) then
                       write(UnFiles((i-1)*4+1),'(F10.3,6(ES13.5))') SimTime,(Disp(i,j),j=1,6)
                   endif
                   ! Velocity file
                   if (FlagOPVel(i)) then
                      write(UnFiles((i-1)*4+2),'(F10.3,6(ES13.5))') SimTime,(Vel(i,j),j=1,6)
                   endif
                   ! Acceleration file
                   if (FlagOPAcc(i)) then
                      write(UnFiles((i-1)*4+3),'(F10.3,6(ES13.5))') SimTime,(Acce(i,j),j=1,6)
                   endif
                   ! Force file
                   if (FlagOPFrc(i)) then
                     write(UnFiles((i-1)*4+4),'(F10.3,6(ES13.5))') SimTime,(Force(i,j),j=1,6)
                   endif
                enddo
                 ! Results of added points 
                if (NumAddedPts>0) then
                    do j=1,NumAddedPts
                       write(UnFiles(4*NumBodies)+(j-1)*3+1,'(F10.3,3(ES13.5))')SimTime,(AddedPtAcce(j,k),k=1,3)
                       write(UnFiles(4*NumBodies)+(j-1)*3+2,'(F10.3,3(ES13.5))')SimTime,(AddedPtVel(j,k),k=1,3)
                       write(UnFiles(4*NumBodies)+(j-1)*3+3,'(F10.3,3(ES13.5))')SimTime,(AddedPtDisp(j,k),k=1,3)
                    enddo
                endif
                
                ! Mooring lines
                TempChar = ''
                do i=1,NumMoorLines
                   write (CharTension,'(ES13.5)') MoorTensionFL(i)
                   TempChar = trim(TempChar)//adjustr(CharTension)
                enddo
                do i=1,NumMoorLines
                   write (CharTension,'(ES13.5)') MoorTensionAC(i)
                   TempChar = trim(TempChar)//adjustr(CharTension)
                enddo
                write(UnFiles(4*NumBodies)+NumAddedPts*3+3,'(F10.3,A)') SimTime,trim(TempChar)
                ! Wave elevation
                if (FlagWvElv) then
                   write (UnFiles(4*NumBodies)+NumAddedPts*3+4,'(F10.3, ES13.5)') SimTime,WvEle
                   write (UnFiles(4*NumBodies)+NumAddedPts*3+5,'(F10.3, 3(ES13.5))') SimTime,WaveVel
                   write (UnFiles(4*NumBodies)+NumAddedPts*3+6,'(F10.3, 3(ES13.5))') SimTime,WaveAcce
                endif
                
                
                wrtFlag = .false.
                if (MOD(SimTime,100.0)==0  ) then   ! monitor the progress
                   write (*,'(A,I5,A)') 'Time = ',NINT(SimTime),' s'
                endif
                
                
             
                iStruct = iStruct -NumBodies
                iStep = iStep+1
                SimTime = (iStep -1)*TimeStep
             endif
             !do i=1,7
             !   read(1,'(A1024)')  LineStr ! Skip the 7 lines for the current body
             !enddo
          elseif (LineStr(23:35) == 'POSITION     ') then
             
             read(LineStr(48:1024),*) (Disp(iStruct,j),j=1,6)                     ! Displacement
             read(1,'(A1024)') LineStr
             read(LineStr(48:1024),*) (Vel(iStruct,j),j=1,6)                      ! Velocity
          elseif (LineStr(1:34) == '                      ACCELERATION') then    ! Acceleration
              read(LineStr(48:1024),*) (Acce(iStruct,j),j=1,6)
          elseif (LineStr(1:36) == '                      EXTERNAL FORCE') then     ! Total force
              read(LineStr(48:1024),*) (Force(iStruct,j),j=1,6)
          elseif( LineStr(1:35) == '                      POSITION NODE') then   ! Added points
              if (NumAddedPts > 0) then
                 do i=1,NumAddedPts
                     if (LineStr(1:41) == charNodeOP(i)) then
                         read (LineStr(42:1024),*) (AddedPtDisp(i,j),j=1,3)
                         read (1,'(A1024)') LineStr
                         read (LineStr(42:1024),*) (AddedPtVel(i,j),j=1,3)
                         read (1,'(A1024)') LineStr
                         read (LineStr(42:1024),*) (AddedPtAcce(i,j),j=1,3)
                         
                     endif
                 enddo
              endif
              
          elseif( LineStr(1:lenWave+22) == '                      '//trim(adjustl(CharWvEle))) then   ! Wave elevation
              read (LineStr(75:1024),*) WvEle          
              read(1,'(A1024)',IOSTAT = status) LineStr
              read (LineStr(42:1024),*) (WaveVel(j),j=1,3)
              read(1,'(A1024)',IOSTAT = status) LineStr
              read (LineStr(42:1024),*) (WaveAcce(j),j=1,3)
          elseif ( LineStr(1:32) == '                      FORCE LINE' .and. FlagMLType) then           !Nonlinear catenary modelling for the mooring lines
              read (LineStr(92:1024),*) MoorTensionFL(1)
              do i=2,NumMoorLines
                 read (1,'(A1024)') LineStr
                 read (LineStr(92:1024),*) MoorTensionFL(i) 
              enddo
              
          elseif( .not. FlagMLType ) then  ! Linar cable modelling for the mooring line
              if (LineStr(1:34) == '                      TENSION LINE')  then  !  linear cable modelling             
                 read (LineStr(102:1024),*) MoorTensionFL(1)
                 do i=2,NumMoorLines
                    read (1,'(A1024)') LineStr
                    read (LineStr(102:1024),*) MoorTensionFL(i) 
                 enddo
              endif
          elseif (LineStr(1:34) == '                      TENSION LINE' .and.  FlagMLType  )  then   ! anchor tension
              read (LineStr(92:1024),*) MoorTensionAC(1)
              do i=2,NumMoorLines
                 read(1,'(A1024)') LineStr
                 read (LineStr(92:1024),*) MoorTensionAC(i)
              enddo
              wrtFlag = .true.
          endif ! current step ending the reading
          
          
          if (iStep > NumSteps ) exit       ! exit the loop
        enddo  ! end of reading file loop
       close(1)    ! LIS file
       do i=1,4*NumBodies
          close(UnFiles(i))
       enddo
       if (NumAddedPts>0) then
           do j=1,NumAddedPts
              close(UnFiles(4*NumBodies)+(j-1)*3+1)
              close(UnFiles(4*NumBodies)+(j-1)*3+2)
              close(UnFiles(4*NumBodies)+(j-1)*3+3)
           enddo
       endif
       close(UnFiles(4*NumBodies)+NumAddedPts*3+1)
       close(UnFiles(4*NumBodies)+NumAddedPts*3+2)
       close(UnFiles(4*NumBodies)+NumAddedPts*3+3)
       close(UnFiles(4*NumBodies)+NumAddedPts*3+4)
       close(UnFiles(4*NumBodies)+NumAddedPts*3+5)
       close(UnFiles(4*NumBodies)+NumAddedPts*3+6)
       ! deallocate the matrix
      if (allocated (Disp)) deallocate(Disp)
      if (allocated (Vel)) deallocate(Vel)
      if (allocated (Acce)) deallocate(Acce)
      if (allocated (Force)) deallocate(Force)
      if (allocated (TdTensionLT)) deallocate(TdTensionLT)
      if (allocated (TdTensionUT)) deallocate(TdTensionUT)
      if (allocated (AddedPtAcce)) deallocate(AddedPtAcce)
      if (allocated (AddedPtVel)) deallocate(AddedPtVel)
      if (allocated (AddedPtDisp)) deallocate(AddedPtDisp)
      if (allocated (MoorTensionFL)) deallocate(MoorTensionFL)
      if (allocated (MoorTensionAC)) deallocate(MoorTensionAC)

  end program readAQWAResults_SingleFile

! *************************************************************************
! get
! *************************************************************************



