!--------------------------------------------------------------!
!                                                              !
!        Switch13lade v5.1                                     !
!        Created by: Thomas Halverson                          !
!        Date: 5-14-13                                         !
!                                                              !
!        MAIN SOURCE CODE                                      !
!        1. Contains NO user interface                         !
!        2. Parallel                                           !
!        3. Uses third party marginals from MM.f90             !
!        4. Uses third party PST from PST.f90                  !
!        5. Uses data file for potential energy                !   
!             PE FILE FORMAT                                   !
!             <dimension> <number of terms>                    !
!             <dimension power> | <coeffient>                  !
!                                                              !
!             mPE FILE FORMAT                                  ! 
!             <mX1@ -20*Delta> ... <mXD@ -20*Delta>            ! 
!                  ...                  ...                    !
!             <mX1@ +20*Delta> ... <mXD@ +20*Delta>            !
!                                                              !
!             PST FILE FORMAT                                  !
!             <mtype>   <Energy Cuttoff>   <Basis Size>        !
!             <m1 n1 ... mD nD>                                !
!                                                              !
!        6. Scalapack parallel diagolization                   !
!             -PDSYGVX                                         !
!        7. Allows for both half and whole integer m           !
!--------------------------------------------------------------!

    !+++++++++++++++++++++++++++++++++++++++++
    !++++++++++  Main Program Loop +++++++++++
    !+++++++++++++++++++++++++++++++++++++++++
Program Main
!----Modules----
 Use UserInputs
 Use Scalapack_Params
!---------------

  Include 'mpif.h'

  Double Precision :: t1, t2, tot_time
  Character(len = 5) :: units

  !Total Time for the entire program
  t1 = MPI_Wtime()

! ------BLACS Intialization------
  Call BLACS_PINFO(IAM, NPROCS)

  If (IAM==0) Then
    Print *,'========================================================'
    Print *,'|     Switch13lade v5.1                                |'
    Print *,'|     Created by: Thomas Halverson                     |'
    Print *,'|     Date: 5-14-13                                    |'
    Print *,'|                                                      |'
    Print *,'|  Solving N-D coupled/anharmonic oscillators using a  |'
    Print *,'|  phase space truncated, momentum symmetrized, doubly |'
    Print *,'|  doubly dense gaussian basis                         |'
    Print *,'========================================================'
    Print *
  End If

! -----Gets Grid info from file------
  Call ReadUserInputs()

! -----Grid Intialization-----
  Call SL_INIT(ICTXT, NPROWS, NPCOLS)
  Call BLACS_GRIDINFO(ICTXT, NPROWS, NPCOLS, MYROW, MYCOL)

  Call GetOneDArrays()
  Call CreateMatrix() 
  Call Diagonalize()
  Call WriteVals()

  t2 = MPI_Wtime()

  If(t2-t1 > 18000.0d0) then
    tot_time = (t2-t1)/(60.0d0*60.0d0)
    units = 'hrs.'
  Else if (t2-t1 > 300) then
    tot_time = (t2-t1)/(60.0d0)
    units = 'mins.'
  Else 
    tot_time = (t2-t1)
    units = 'sec.'
  End If

  If(IAM == 0) Then
    Print 9999, tot_time, units
    Print *,'*****************************************************************'
    Print *
    Print *
  End If

  9999 Format ('  Total time taken: ', f10.5, 2x, a)
  
  Call BLACS_GRIDEXIT(ICTXT)
  Call BLACS_EXIT(0)
       
End Program Main
    !+++++++++++++++++++++++++++++++++++++++++
    !++++++  End of Main Program Loop ++++++++
    !+++++++++++++++++++++++++++++++++++++++++


    !************************************************************!
    !********************    SUBROUTINES   **********************!
    !************************************************************!

  Subroutine ReadUserInputs()
    !----Modules----
    Use Constants
    Use UserInputs
    Use Outputs
    Use ParameterArrays
    Use Scalapack_Params
    !---------------
    
    Implicit NONE
    Integer :: i,j
    Character(len = 40) :: input_file
 
    Call GETARG(1, input_file)
    input_file = trim(input_file)
    
    If (IAM == 0) Then
      Print *,'*****************************************************************'
      Print *,'        STEP 1:'
      Print '(x,a,2x,a20)',"Reading data from: ", input_file
      Print *
    End If
 
!   ---- Read in files names (and directories) of user inputs ----  
    Open(UNIT = 15, FILE = input_file, STATUS='OLD',  ACTION ='Read')
    
      Read (15,'(A)') PE_FN
      PE_FN = trim(PE_FN)
      Read (15,'(A)') masses_FN
      masses_FN = trim(masses_FN)
      Read (15,'(A)') SysParams_FN
      SysParams_FN = trim(SysParams_FN)      
      Read (15,'(A)') psp_FN
      psp_FN = trim(psp_FN)   
      Read (15,'(A)') Out_FN
      Out_FN = trim(Out_FN)   

    Close(15)
!   -----------------------------------------


    If(IAM == 0) then
      Print *, ' Use requested files: '
      Print 1002, "PE:                 ", PE_FN
      Print 1002, "Masses:             ", masses_FN
      Print 1002, "Sytem parameters:   ", SysParams_FN
      Print 1002, "Phase space points: ", psp_FN
      Print 1002, "Output:             ", Out_FN
      Print *
      Print *, ' Reading data from suggested files'
    End If  

   
!   ---- Read in user input system paramters from file ----  
    Open(UNIT = 15, FILE = SysParams_FN, STATUS='OLD',  ACTION ='Read')
    
      Read (15,*) NPROWS
      Read (15,*) NPCOLS
      Read (15,*) BLOCK_ROW
      BLOCK_COL = BLOCK_ROW
    
    Close(15)
!   -----------------------------------------



!   -----Read in PE info from file-----
    Open(UNIT = 15, FILE = PE_FN, STATUS='OLD',  ACTION ='Read')
    
      Read(15,*) Dmax, coeffs_len

      Allocate(Coeffs(coeffs_len))
      Allocate(Powers(coeffs_len, Dmax))
     
      Do i = 1, coeffs_len
        Read (15,*) (Powers(i,j), j=1, Dmax), Coeffs(i)
      End Do 
    
    Close(15)
!   -----------------------------------



!   ----Read masses from file-----
    Open(15, file = masses_FN, STATUS='OLD',  ACTION ='Read')
    
      Allocate(Masses(DMax))
     
      Do i = 1, Dmax
        Read(15,*) Masses(i)
      End Do
            
    Close(15)
!   ------------------------------


  
!   ---- Read phase space points from file-----
    Open(15, file = psp_FN, STATUS='OLD',  ACTION ='Read')

      Read(15,*) mtype, Emax, iMax
      Allocate(MandN(iMax,2*Dmax))
    
      Do i = 1, iMax
        Read(15,*) (MandN(i,j), j=1, 2*Dmax)
      End Do

      If(mtype==1) then
        shift=0.50d0
      Else if (mtype==2) then
        shift=1.0d0
      Else
        Print *, "ERROR 1: INVALID SWITCH VALUE"
      End If
!   ------------------------------

    

    If(IAM == 0) Then
      Print *, " Data read: SUCCESS"
      Print *
      Print *, "    _________________________________________________________"
      Print 1000, "Dimension:             ", Dmax
      Print 1000, "Number of parameters:  ", coeffs_len
      Print 1003, "Energy cuttoff:        ", EMax
      Print 1000, "Basis size:            ", iMax
      Print 1000, "Process grid rows:     ", NPROWS
      Print 1000, "Process grid columns:  ", NPCOLS
      If (mtype == 1) then
        Print '(10x,a)',  " Lattice spacing:     half integer"
      Else if (mtype ==2) then
        Print '(10x,a)',  " Lattice spacing:     whole integer"
      Else
        Print *, "ERROR 1: INVALID SWITCH VALUE"
      End If
      Print *
      Print *, "    _________________________________________________________"
      Print *
      Print *,'*****************************************************************'
      Print *
      Print *
      Print *
      
      1000 Format(10x, a, i13)
      1001 Format(7x, a, f13.1)
      1002 Format(7x, a29, a)
      1003 Format(10x, a, f13.5)
    End If 
    
  End Subroutine ReadUserInputs
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine GetOneDArrays()
!   - GetOneDArrays subroutine poplulates the one dimensional arrays
!     from the precreated half integer or whole integer m data files 
!     depending on the user specification. The 1D arrays are used 
!     to create the Hamiltonian and overlap matrices. 
!       
!-----------------------------------------------------------------------------
  Subroutine GetOneDArrays()

   !----Modules-----
    Use Constants
    Use ParameterArrays
    Use UserInputs
    Use Scalapack_Params
   !----------------    

    Implicit NONE
    
    !Declare Variables
    Integer :: u,v,up,vp 	!Indices for the Loop
    Integer :: index
  
    If (mtype == 1) then 
      Allocate(OneDArrays_Array(2*Size,Size,2*Size,Size,6))
      Open(UNIT = 10, FILE = './OneDarrays/SArray_HI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 11, FILE = './OneDarrays/XArray_HI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 12, FILE = './OneDarrays/X2Array_HI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 13, FILE = './OneDarrays/X3Array_HI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 14, FILE = './OneDarrays/X4Array_HI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 15, FILE = './OneDarrays/P2Array_HI.dat',STATUS='OLD',  ACTION ='Read')
    Else if (mtype == 2) then
      Allocate(OneDArrays_Array(2*Size+1,Size,2*Size+1,Size,6))
      Open(UNIT = 10, FILE = './OneDarrays/SArray_WI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 11, FILE = './OneDarrays/XArray_WI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 12, FILE = './OneDarrays/X2Array_WI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 13, FILE = './OneDarrays/X3Array_WI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 14, FILE = './OneDarrays/X4Array_WI.dat',STATUS='OLD',  ACTION ='Read')
      Open(UNIT = 15, FILE = './OneDarrays/P2Array_WI.dat',STATUS='OLD',  ACTION ='Read')
    Else
      Print *, "ERROR 2: MISREAD ARRAYS"
    End IF  
    
    If(Iam == 0) then
      If (mtype == 1) then
        Print *,'*****************************************************************'
        Print *, '        STEP 2:'
        Print *,"   Initializing 1-D arrays for half integer m"
      Else if (mtype == 2) then
        Print *,'*****************************************************************'
        Print *, '        STEP 2:'
        Print *,"   Initializing 1-D arrays for whole integer m"
      End If
    End If
       
    If (mtype==1) then
      Do u = 1,(2*Size)
        Do v = 1, Size
          Do up = 1,(2*Size)
            Do vp = 1, Size 
              Do index = 1, 6

                Read(index+9,*) OneDArrays_Array(u,v,up,vp,index)

              End Do
            End Do
          End Do
        End Do
      End Do

    Else If (mtype==2) then
      Do u = 1,(2*Size) + 1
        Do v = 1, Size
          Do up = 1,(2*Size) + 1
            Do vp = 1, Size 
              Do index = 1, 6

                Read(index+9,*) OneDArrays_Array(u,v,up,vp,index)

              End Do
            End Do
          End Do
        End Do
      End Do
    End If
    
    Close(UNIT = 10)
    Close(UNIT = 11)
    Close(UNIT = 13)
    Close(UNIT = 14)
    Close(UNIT = 15)

    If (IAM == 0) then
      Print*,'   1-D array initialization: SUCCESS'
      Print *,'*****************************************************************'
      Print *
      Print *
      Print *
    End If

End Subroutine GetOneDArrays
!-----------------------------------------------------------------------------





!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine CreateMatrix()
!   - CreateMatrix subroutine builds the Hamiltonian and overlap matrices
!     from the given list of truncated phase space points MandN
!
!-----------------------------------------------------------------------------
  Subroutine CreateMatrix()

   !----Modules-----
    Use Constants
    Use UserInputs
    Use ParameterArrays
    Use Outputs
    Use Scalapack_Params
   !----------------   
  
    Implicit NONE
    
            
    Integer          :: eye, jay, kay   !Indicies for the loops
    Integer          :: index, power
    Double Precision :: KE, PE, dummy
    Double Precision :: t1,t2,tot_time

    Character(len = 5) :: units
    
    !-------------------------BLACS Vars----------------------------
    Integer :: RSRC, CSRC, LLD, INFO
    Integer :: GLOBAL_rindex, GLOBAL_cindex
    Integer :: LOC_ROW, LOC_COL
    Integer :: LOC_roffset, LOC_coffset, LOC_rindex, LOC_cindex
    Integer :: ROW_first, COL_first, ROW_end, COL_end, ROW_start, COL_start
    !---------------------------------------------------------------
  
    !External  and BLACS functions
    Double Precision MPI_Wtime
    Integer numroc, indxg2p
    
    Integer, Dimension(Dmax) :: Yous
    Integer, Dimension(Dmax) :: Vees
    Integer, Dimension(Dmax) :: YouPrimes
    Integer, Dimension(Dmax) :: VeePrimes
    
    Double Precision, Dimension(iMax) :: Work !Needed for PDLAPRT
  
    t1 = MPI_Wtime()
    
    !----------------Intitialize Descriptor-------------------
    !-----------Preparing for block-cylic distribution--------
    !Proc offset
    RSRC = 0
    CSRC = 0
    
    !------------Find local array size----------------
    !-------------------------------------------------
    LOC_COL = numroc(iMax, BLOCK_COL, NPCOLS, CSRC, NPCOLS)
    LOC_COL = max(1,LOC_COL)
    LOC_ROW = numroc(iMax, BLOCK_ROW, NPROWS, RSRC, NPROWS)
    LLD = max(LOC_ROW,1)
    !-------------------------------------------------
  
    Call DESCINIT(DESCH, iMax, iMax, BLOCK_ROW, BLOCK_COL, RSRC, CSRC, ICTXT, LLD, INFO)
    Call DESCINIT(DESCS, iMax, iMax, BLOCK_ROW, BLOCK_COL, RSRC, CSRC, ICTXT, LLD, INFO)    
  
    Allocate(Ham(LOC_ROW, LOC_COL))
    Allocate(Overlap(LOC_ROW, LOC_COL))
  
    
    If (IAM == 0) Then
      Print *,'*****************************************************************'
      Print *, '        STEP 3:'
      Print *, 'Preparing to create and distribute Hamiltonian and overlap matrices'
      Print *
      Print *, ' Distribution type: block-cyclic'
      Print 1000, 'Total number of processors: ', NPROCS
      Print 1001, 'Proccess grid:              ', NPROWS, ' x ', NPCOLS
      Print 1001, 'Global array size:          ', imax, ' x ', imax
      Print 1001, 'Block size:                 ', BLOCK_ROW, ' x ', BLOCK_COL
      Print 1001, 'Local array size:           ', LOC_ROW, ' x ', LOC_COL
      Print *
    End If
  
    1000 Format(2x, a, 6x, i4)
    1001 Format(2x, a, i10, 5x, a, i10)
  
    !---------------Begin Block-Cyclic Distribution---------
    !-------------------------------------------------------
  
    !Establish the starting point on each processor
    If (MYROW >= DESCH(RSRC_)) Then
      ROW_first = (MYROW - DESCH(RSRC_))*DESCH(BLOCK_ROWS_) + 1
    Else
      ROW_first = (MYROW + (NPROWS - DESCH(RSRC_))) * DESCH(BLOCK_ROWS_) + 1
    End If
  
    If (MYCOL >= DESCH(CSRC_)) Then
      COL_first = (MYCOL - DESCH(CSRC_)) * DESCH(BLOCK_COLs_) + 1
    Else
      COL_first = (MYCOL + (NPCOLS - DESCH(CSRC_))) * DESCH(BLOCK_COLS_) + 1
    End If
  
    call blacs_barrier(ICTXT, 'All')
  
    !Move through A by block
    Do COL_start = COL_first, DESCH(COLS_), NPCOLS * DESCH(BLOCK_COLS_)
      Do ROW_start = ROW_first, DESCH(ROWS_), NPROWS * DESCH(BLOCK_ROWS_)
  
        !find the last index in the block
        ROW_end = min( DESCH(ROWS_), ROW_start + DESCH(BLOCK_ROWS_)-1)
        COL_end = min( DESCH(COLS_), COL_start + DESCH(BLOCK_COLS_)-1)
  
        !find the local offset
        call infog2l(ROW_start, COL_start, DESCH, NPROWS, NPCOLS, MYROW, MYCOL, &
          &  LOC_roffset, LOC_coffset, RSRC, CSRC)
  
        !Move through each block
        Do GLOBAL_cindex = COL_start, COL_end
          Do GLOBAL_rindex = ROW_start, ROW_end
  
            LOC_rindex = LOC_roffset + (GLOBAL_rindex - ROW_start)
            LOC_cindex = LOC_coffset + (GLOBAL_cindex - COL_start)
        
            !Fills up the u and v indices arrays
            Do jay = 1, Dmax
              Yous(jay)      = INT(MandN(GLOBAL_cindex,2*jay-1) + Size + shift)
              YouPrimes(jay) = INT(MandN(GLOBAL_rindex,2*jay-1) + Size + shift)    
            
              Vees(jay)      = INT(MandN(GLOBAL_cindex,2*jay) + .5)
              VeePrimes(jay) = INT(MandN(GLOBAL_rindex,2*jay) + .5)
            End Do
  
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Calculates Kenetic Energy  >>>>>>>>>>>>>>>>>>>>>>>>>
            Ham(LOC_rindex,LOC_cindex) = 0
        
            Do jay = 1, Dmax
        
              KE = (hbar**2)/(2.0d0*Masses(jay))*OneDArrays_Array(Yous(jay),Vees(jay), &
                     & YouPrimes(jay),VeePrimes(jay), 6) 
          
              Do kay = 1, jay - 1
                KE = KE * OneDArrays_Array(Yous(kay),Vees(kay),YouPrimes(kay),VeePrimes(kay),1)
              End Do
              Do kay = jay + 1, Dmax            
                KE = KE * OneDArrays_Array(Yous(kay),Vees(kay),YouPrimes(kay),VeePrimes(kay),1)
              End Do
      
              Ham(LOC_rindex,LOC_cindex) =  Ham(LOC_rindex,LOC_cindex) + KE
        
            End Do
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Calculates Potential Energy  >>>>>>>>>>>>>>>>>>>>>>>>>
            PE = 0.0d0
            Do jay = 1, coeffs_len
              dummy = 1.0d0
              Do kay = 1, Dmax
            
                power = Powers(jay,kay) + 1
                dummy = dummy * OneDArrays_Array(Yous(kay),Vees(kay),YouPrimes(kay),VeePrimes(kay), power)
        
              End Do
              PE = PE + Coeffs(jay)*dummy
            End Do
       

!If(GLOBAL_cindex==GLOBAL_rindex .AND. PE<0) then
!  Write(*,'(2x,i3,2x,i3,2x,f7.3)') GLOBAL_rindex,GLOBAL_cindex,PE
!End If

            Ham(LOC_rindex,LOC_cindex) =  Ham(LOC_rindex,LOC_cindex) + PE
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Calculates Overlap Matrix  >>>>>>>>>>>>>>>>>>>>>>>>>
            Overlap(LOC_rindex,LOC_cindex) = 1
            Do jay = 1, Dmax
              Overlap(LOC_rindex,LOC_cindex) = Overlap(LOC_rindex,LOC_cindex) * &
                  &  OneDArrays_Array(Yous(jay),Vees(jay),YouPrimes(jay),VeePrimes(jay),1)
            End Do
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
          End Do
        End Do
      End Do
    End Do
  
    call blacs_barrier(ICTXT, 'All')
  
    t2 = MPI_Wtime()

    If(t2-t1 > 18000.0d0) then
      tot_time = (t2-t1)/(60.0d0*60.0d0)
      units = 'hrs.'
    Else if (t2-t1 > 300) then
      tot_time = (t2-t1)/(60.0d0)
      units = 'mins.'
    Else 
      tot_time = (t2-t1)
      units = 'sec.'
    End If
    
    If (IAM == 0) then
      Print *, 'Distribution: COMPLETE'
      Print *
      Print 9999, tot_time,units
      Print *,'*****************************************************************'  
      Print *
      Print *
      Print *
    End If  

    Deallocate(Coeffs)
    Deallocate(Powers)
    Deallocate(Masses)
    Deallocate(OneDArrays_Array)
    Deallocate(MandN)
  
    9999 Format ('  Time taken to create and distribute: ',f10.5, 2x, a)
  
  End SubRoutine CreateMatrix
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine Diagonalize()
!   - Diagonalize subroutine calculates the eigenvalues of the Hamiltonian
!    using the Scalapack routine PDSYGVX
!
!------------------------------------------------------------------------------
SubRoutine Diagonalize()

   !----Modules-----
    Use UserInputs
    Use Outputs
    Use Scalapack_Params
   !----------------  

	Implicit NONE

	Integer           :: IBTYPE, IA, JA, IB, JB, VL, VU, IL, IU
    Integer           :: M, NZ, IZ, JZ, LWORK, LIWORK, INFO
    Character*1       :: JOBZ, RANGE, UPLO
    Double Precision  :: ABSTOL, ORFAC, t1, t2, tot_time
    Integer           :: NB, NN, NP0, NNP

    Character(len = 5) :: units
    
    Integer,                       Dimension(DLEN)            :: DESCZ
    Integer,                       Dimension(iMax)            :: IFAIL
    Integer,                       Dimension(2*NPROWS*NPCOLS) :: ICLUSTR
	Double Precision, Allocatable, Dimension(:,:)             :: Z
    Double Precision, Allocatable, Dimension(:)               :: WORK
    Integer,          Allocatable, Dimension(:)               :: IWORK
    Double Precision,              Dimension(NPROWS*NPCOLS)   :: GAP
    
	Double Precision PDLAMCH, MPI_Wtime, NUMROC
   
	IBTYPE = 1
    JOBZ = 'N'
    RANGE = 'A'
    UPLO = 'U'
    
	IA = 1
    JA = 1
    IB = 1
    JB = 1
    IZ = 1
    JZ = 1
    VL = 0
    VU = 0
    IL = iMax
    IU = iMax

    NB = DESCH( BLOCK_ROWS_ )
    NN = MAX( iMax, NB, 2 )
    NP0 = NUMROC( NN, NB, 0, 0, NPROWS )
    NNP = MAX( iMax, NPROWS*NPCOLS + 1, 4 )
    
    ABSTOL = -1.0d0
    ORFAC = -1.0d0

    t1 = MPI_Wtime()
 
    If(IAM==0) Then
      Print *,'*****************************************************************'
      Print *, '        STEP: 4'
      Print *, 'Diagonalizing Hamiltonian'
      Print *
      Print *, ' Using direct diagonalization'
      Print *, ' Scalapack routine:   PDSYGVX'
      Print *
    End If

    Allocate(Eigen_Vals(iMax))
    Allocate(Z(DESCH(LLD_),DESCH(LLD_)))
    Call DESCINIT(DESCZ, iMax, iMax, BLOCK_ROW, BLOCK_COL, DESCH(RSRC_), DESCH(CSRC_), &
                   &  ICTXT, DESCH(LLD_), INFO)

    If(iMax < 75000) then
      
      LWORK = -1
      LIWORK = -1

      Allocate(WORK(1))
      Allocate(IWORK(1))

      Call PDSYGVX( IBTYPE, JOBZ, RANGE, UPLO, iMax, Ham, IA, JA, &
                        & DESCH, Overlap, IB, JB, DESCS, VL, VU, IL, IU, &
                        & ABSTOL, M, NZ, Eigen_Vals, ORFAC, Z, IZ, JZ, DESCZ, &
                        & WORK, LWORK, IWORK, LIWORK, IFAIL, ICLUSTR, &
                        & GAP, INFO )
                        
      LWORK = WORK(1)
      LIWORK = IWORK(1)

      Deallocate(WORK)
      Deallocate(IWORK)
      Allocate(WORK(LWORK))
      Allocate(IWORK(LIWORK))


    Else 
    
      LWORK = 10 * iMax + MAX( 5 * NN, NB * ( NP0 + 1 ) )
      LIWORK = 10 * NNP

      Allocate(WORK(LWORK))
      Allocate(IWORK(LIWORK))

    End If
      
    If(IAM==0) Then
      Print *, '       Workspace Sizes:'
      Print *, '------------------------------------'
      Print *, '  Work array length: ',LWORK
      Print *, '  IWork array legnth: ', LIWORK
      Print *, '------------------------------------'
      Print *
    End If
    

    call blacs_barrier(ICTXT, 'All')

    !Print *,'GOTHERE'

    Call PDSYGVX( IBTYPE, JOBZ, RANGE, UPLO, iMax, Ham, IA, JA, &
                        & DESCH, Overlap, IB, JB, DESCS, VL, VU, IL, IU, &
                        & ABSTOL, M, NZ, Eigen_Vals, ORFAC, Z, IZ, JZ, DESCZ, &
                        & WORK, LWORK, IWORK, LIWORK, IFAIL, ICLUSTR, &
                        & GAP, INFO )

    t2 = MPI_Wtime()

    If(t2-t1 > 18000.0d0) then
      tot_time = (t2-t1)/(60.0d0*60.0d0)
      units = 'hrs.'
    Else if (t2-t1 > 300) then
      tot_time = (t2-t1)/(60.0d0)
      units = 'mins.'
    Else 
      tot_time = (t2-t1)
      units = 'sec.'
    End If

    If (IAM == 0) Then
      Print *, 'Diagonalization: COMPLETE'
      Print *
      Print 9999, tot_time, units
      Print *,'*****************************************************************'
      Print *
      Print *
      Print *
    End If
    
    9999 Format ('  Time taken to diagonalize: ',f10.5,2x,a)

End SubRoutine
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine PrintVals(Emax)
!   - PrintVals subroutine allows the user to output a subset of the
!     eigenvalues to the screen. The user inputs a range, <min, max>.
!
!     Emax
!     - Double precision. Input parameter given by the user.  
!
!-----------------------------------------------------------------------------
Subroutine WriteVals()
    
   !----Modules-----
    Use Outputs
    Use UserInputs
    Use Scalapack_Params
   !----------------  

	Implicit NONE
    Integer	:: i
 
    call blacs_barrier(ICTXT, 'All')
   
    If(IAM == 0) Then

      Print *
      Print *,'*****************************************************************'
      Print *,'       SUMMARY'
      Print *
      Print 1000, Dmax
	  Print 1001, EMax
	  Print 1002, iMax

 	  Open(UNIT = 15, FILE = Out_FN, ACTION='WRITE')    
    
	  Do i = 1, iMax
    	  Write (15,*) Eigen_Vals(i)
      End Do

	  Close(UNIT = 15)
    
	  Print 1003, Out_FN
      Print *

      1000 Format (2x,i2, '-D Coupled Anharmonic Oscillator') 
	  1001 Format ('  Energy cut off:              ', 3x, f13.5)
      1002 Format ('  Total number of eigenvalues: ', 3x, i7)
      1003 Format ('  Eigenvalues read to file:    ', a40)

    End If
    
End Subroutine WriteVals
!-----------------------------------------------------------------------------