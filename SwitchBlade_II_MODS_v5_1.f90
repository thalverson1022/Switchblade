!-------------------------------------
! Constants and variables module for
!       SwitchBlade v5.1 (parallel)
!
!     Created by: Tom Halverson
!     Date: 4-14-13
!-------------------------------------
!-------------------------------------------------------------------------------------------
Module Constants
  Implicit NONE
  Double Precision, Parameter :: PI = 3.1415926535897932385d0
  Double Precision, Parameter :: CONV = 219474.6313705d0
  Integer,          Parameter :: SIZE = 20                     !Parameter for the large arrays
  Integer,          Parameter :: Length = 500000               !Maximum length for MandN
  Double Precision, Parameter :: hbar = 1.0d0	
End Module Constants
!-------------------------------------------------------------------------------------------
Module UserInputs
  Implicit NONE  
  Integer          :: DMax         !Problem Dimensionality
  Integer          :: coeffs_len   !Total number of coulping/anharmonic constants  
  Double Precision :: EMax         !Truncation parameter 
  Double Precision :: shift        !m grid shift
  Integer		   :: mtype        !Switch for (1) half or (2) whole int grid

  !Blacs grid info 
  Integer		   :: NPROWS       !Number of processor rows
  Integer		   :: NPCOLS       !Number of processor columns
  Integer		   :: BLOCK_ROW    !Number of rows in each block
  Integer          :: BLOCK_COL    !NUmber of cloumns in each block
  
  Character(len = 128) :: Out_FN       !Output file: eigenvalues
  Character(len = 128) :: PE_FN        !Input file: potential energy
  Character(len = 128) :: masses_FN    !Input file: masses
  Character(len = 128) :: SysParams_FN  !Input file: Scalapack Grid Info
  Character(len = 128) :: psp_FN       !Input file: valid phase space points
  
End Module UserInputs
!-------------------------------------------------------------------------------------------
Module Outputs
  Implicit NONE
  Integer                                       :: iMax         !Numbe of valid phase space points
  Double Precision,	Allocatable, Dimension(:,:) :: MandN   		!Array for Phase Space Tuncation 
  Double Precision, Allocatable, Dimension(:,:)	:: Ham	   		!Hamiltonian Matrix
  Double Precision, Allocatable, Dimension(:,:)	:: Overlap 		!Overlap Matrix
  Double Precision, Allocatable, Dimension(:)	:: Eigen_Vals	!Array for Eigen Value Storage
End Module Outputs
!-------------------------------------------------------------------------------------------
Module ParameterArrays
  Implicit NONE
  Double Precision, Allocatable, Dimension(:)         :: Coeffs
  Integer,          Allocatable, Dimension(:,:)       :: Powers
  Double Precision, Allocatable, Dimension(:,:)       :: Marginal_PE
  Double Precision, Allocatable, Dimension(:)         :: Masses  
  Double Precision, Allocatable, Dimension(:,:,:,:,:) :: OneDArrays_Array
End Module
!-------------------------------------------------------------------------------------------
Module Scalapack_Params
  
  Implicit NONE

  !Blacs grid vars
  Integer            :: ICTXT      !BLACS context variable
  Integer 			 :: IAM, NPROCS
  Integer 			 :: MYROW, MYCOL
  
  Integer, Parameter :: DLEN = 9   !Length of descriptor DESC

  !Variables for each element of DESC_Variable_() array
  !DESC_Var_={DTYPE,ICTXT,ROWS,COLS,BLOCK_ROWS,BLOCK_COLS,RSRC,CSRC,LLD}
  Integer, Parameter :: BLOCK_CYCLIC_2D = 1
  Integer, Parameter :: DTYPE_          = 1 
  Integer, Parameter :: ICTXT_          = 2
  Integer, Parameter :: ROWS_ 			= 3     
  Integer, Parameter :: COLS_ 			= 4
  Integer, Parameter :: BLOCK_ROWS_		= 5
  Integer, Parameter :: BLOCK_COLS_		= 6
  Integer, Parameter :: RSRC_			= 7
  Integer, Parameter :: CSRC_			= 8
  Integer, Parameter :: LLD_			= 9
  
  Integer, Dimension(DLEN) :: DESCH, DESCS

End Module Scalapack_Params
!-------------------------------------------------------------------------------------------
