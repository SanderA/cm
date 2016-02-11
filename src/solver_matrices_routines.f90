!> \file
!> \author Chris Bradley
!> \brief This module handles all solver matrix and rhs routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module handles all solver matrix and rhs routines.
MODULE SOLVER_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters


  !> \addtogroup SOLVER_MATRICES_ROUTINES_SelectMatricesTypes SOLVER_MATRICES_ROUTINES::SelectMatricesTypes
  !> \brief The types of selection available for the solver matrices
  !> \see SOLVER_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_ALL=1 !<Select all the solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
!  redundant when introducing dynamic nonlinear equations
!  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_LINEAR_ONLY=3 !<Select only the linear solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian solver matrix \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual solver vector \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_ONLY=7 !<Select only the RHS solver vector \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_RESIDUAL_ONLY=8 !<Select only the residual and RHS solver vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLVER_MATRICES_ALL,SOLVER_MATRICES_LINEAR_ONLY,SOLVER_MATRICES_NONLINEAR_ONLY, &
    & SOLVER_MATRICES_JACOBIAN_ONLY,SOLVER_MATRICES_RESIDUAL_ONLY,SOLVER_MATRICES_RHS_ONLY, & 
    & SOLVER_MATRICES_RHS_RESIDUAL_ONLY !,SOLVER_MATRICES_DYNAMIC_ONLY

  PUBLIC SOLVER_MATRIX_EQUATIONS_MATRIX_ADD,SOLVER_MATRIX_INTERFACE_MATRIX_ADD,SOLVER_MATRIX_JACOBIAN_MATRIX_ADD
  
  PUBLIC SOLVER_MATRICES_CREATE_FINISH,SOLVER_MATRICES_CREATE_START,SOLVER_MATRICES_DESTROY,SOLVER_MATRICES_LIBRARY_TYPE_SET, &
    & SOLVER_MATRICES_OUTPUT,SOLVER_MATRICES_STORAGE_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating the solver matrices
  SUBROUTINE SOLVER_MATRICES_CREATE_FINISH(SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrixIdx,numberOfNonZeros
    INTEGER(INTG), ALLOCATABLE :: columnIndices(:),rowIndices(:)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    ENTERS("SOLVER_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL FlagError("Solver matrices have already been finished",ERR,ERROR,*998)
      ELSE
        SOLVER_EQUATIONS=>SOLVER_MATRICES%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            !Now create the individual solver matrices
            IF(ASSOCIATED(SOLVER_MAPPING%DOFS_MAPPING)) THEN
              DO matrixIdx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrixIdx)%PTR
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  !!Create the distributed solver matrix
                  CALL DISTRIBUTED_MATRIX_CREATE_START(SOLVER_MAPPING%DOFS_MAPPING,SOLVER_MAPPING%DOFS_MAPPING, &
                    & SOLVER_MATRICES%MATRICES(matrixIdx)%PTR%MATRIX,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_GHOSTING_TYPE_SET(SOLVER_MATRIX%MATRIX,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_LIBRARY_TYPE_SET(SOLVER_MATRIX%MATRIX,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(SOLVER_MATRIX%MATRIX,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(SOLVER_MATRIX%MATRIX,SOLVER_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                  !Calculate and set the matrix structure/sparsity pattern
                  IF(SOLVER_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                    & SOLVER_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                    CALL SOLVER_MATRIX_STRUCTURE_CALCULATE(SOLVER_MATRIX,numberOfNonZeros,rowIndices, &
                      & columnIndices,ERR,ERROR,*999)                  
                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(SOLVER_MATRIX%MATRIX,numberOfNonZeros, &
                      & ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(SOLVER_MATRIX%MATRIX,rowIndices,columnIndices, &
                      & ERR,ERROR,*999)
                    IF(ALLOCATED(rowIndices)) DEALLOCATE(rowIndices)
                    IF(ALLOCATED(columnIndices)) DEALLOCATE(columnIndices)
                  ENDIF
                  CALL DISTRIBUTED_MATRIX_CREATE_FINISH(SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
                  !Allocate the distributed solver vector
                  CALL DISTRIBUTED_VECTOR_CREATE_START(SOLVER_MAPPING%DOFS_MAPPING,SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_GHOSTING_TYPE_SET(SOLVER_MATRIX%SOLVER_VECTOR, &
                    & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRIX%SOLVER_VECTOR,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRIX%SOLVER_VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !matrixIdx
              IF(SOLVER_EQUATIONS%LINEARITY==PROBLEM_SOLVER_NONLINEAR) THEN
                !Allocate the nonlinear matrices and vectors                  
                !Allocate the distributed residual vector
                CALL DISTRIBUTED_VECTOR_CREATE_START(SOLVER_MAPPING%DOFS_MAPPING,SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_GHOSTING_TYPE_SET(SOLVER_MATRICES%RESIDUAL, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRICES%RESIDUAL,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRICES%RESIDUAL,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)                  
              ENDIF
!!TODO: what to do if there is no RHS
              !Allocate the distributed rhs vector
              CALL DISTRIBUTED_VECTOR_CREATE_START(SOLVER_MAPPING%DOFS_MAPPING,SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_GHOSTING_TYPE_SET(SOLVER_MATRICES%RHS_VECTOR, &
                & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRICES%RHS_VECTOR,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRICES%RHS_VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
              !Finish up
              SOLVER_MATRICES%SOLVER_MATRICES_FINISHED=.TRUE.
            ELSE
              CALL FlagError("Row domain mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FlagError("Solver matrices solver equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",ERR,ERROR,*998)
    ENDIF
        
    EXITS("SOLVER_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ALLOCATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ALLOCATED(columnIndices)) DEALLOCATE(columnIndices)
    CALL SOLVER_MATRICES_FINALISE(SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("SOLVER_MATRICES_CREATE_FINISH",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE SOLVER_MATRICES_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating the solver matrices
  SUBROUTINE SOLVER_MATRICES_CREATE_START(SOLVER_EQUATIONS,SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to create the solver matrices for
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<On return, a pointer to the solver matrices. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    ENTERS("SOLVER_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(SOLVER_MATRICES)) THEN
          CALL FlagError("Solver matrices is already associated",ERR,ERROR,*998)
        ELSE
          NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES)
          CALL SOLVER_MATRICES_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
        ENDIF
      ELSE
        CALL FlagError("Solver equations are not finished",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated",ERR,ERROR,*998)
    ENDIF
        
    EXITS("SOLVER_MATRICES_CREATE_START")
    RETURN
999 CALL SOLVER_MATRICES_FINALISE(SOLVER_EQUATIONS%SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("SOLVER_MATRICES_CREATE_START",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_CREATE_START
        
  !
  !================================================================================================================================
  !

  !>Destroy the solver matrices
  SUBROUTINE SOLVER_MATRICES_DESTROY(SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer the solver matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SOLVER_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      CALL SOLVER_MATRICES_FINALISE(SOLVER_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Solver matrices is not associated",ERR,ERROR,*999)
    ENDIF
        
    EXITS("SOLVER_MATRICES_DESTROY")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_DESTROY",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the solver matrices and deallocates all memory
  SUBROUTINE SOLVER_MATRICES_FINALISE(SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx

    ENTERS("SOLVER_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(ALLOCATED(SOLVER_MATRICES%MATRICES)) THEN
        DO matrixIdx=1,SIZE(SOLVER_MATRICES%MATRICES,1)
          CALL SOLVER_MATRIX_FINALISE(SOLVER_MATRICES%MATRICES(matrixIdx)%PTR,ERR,ERROR,*999)
        ENDDO !matrixIdx
        DEALLOCATE(SOLVER_MATRICES%MATRICES)
      ENDIF
      IF(ASSOCIATED(SOLVER_MATRICES%RESIDUAL)) CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLVER_MATRICES%RHS_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(SOLVER_MATRICES)
    ENDIF
        
    EXITS("SOLVER_MATRICES_FINALISE")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Initialises the solver matrices for solver equations
  SUBROUTINE SOLVER_MATRICES_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to initialise the solver matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrixIdx
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("SOLVER_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MATRICES)) THEN
        CALL FlagError("Solver matrices is already associated for this solver equations.",ERR,ERROR,*998)
      ELSE
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          ALLOCATE(SOLVER_EQUATIONS%SOLVER_MATRICES,STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate solver matrices.",ERR,ERROR,*999)
          SOLVER_EQUATIONS%SOLVER_MATRICES%SOLVER_EQUATIONS=>SOLVER_EQUATIONS
          SOLVER_EQUATIONS%SOLVER_MATRICES%SOLVER_MATRICES_FINISHED=.FALSE.
          SOLVER_EQUATIONS%SOLVER_MATRICES%SOLVER_MAPPING=>SOLVER_MAPPING
          SOLVER_EQUATIONS%SOLVER_MATRICES%NUMBER_OF_ROWS=SOLVER_MAPPING%NUMBER_OF_DOFS
          SOLVER_EQUATIONS%SOLVER_MATRICES%NUMBER_OF_GLOBAL_ROWS=SOLVER_MAPPING%NUMBER_OF_GLOBAL_DOFS
          SOLVER_EQUATIONS%SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
          SOLVER_EQUATIONS%SOLVER_MATRICES%NUMBER_OF_MATRICES=SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
          ALLOCATE(SOLVER_EQUATIONS%SOLVER_MATRICES%MATRICES(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate solver matrices matrices.",ERR,ERROR,*999)
          DO matrixIdx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
            NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES%MATRICES(matrixIdx)%PTR)
            CALL SOLVER_MATRIX_INITIALISE(SOLVER_EQUATIONS%SOLVER_MATRICES,matrixIdx,ERR,ERROR,*999)
          ENDDO !matrixIdx
          IF(SOLVER_EQUATIONS%LINEARITY==PROBLEM_SOLVER_NONLINEAR) THEN
            SOLVER_EQUATIONS%SOLVER_MATRICES%UPDATE_RESIDUAL=.TRUE.
          ELSE
            SOLVER_EQUATIONS%SOLVER_MATRICES%UPDATE_RESIDUAL=.FALSE.
          ENDIF
          NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES%RESIDUAL)
          SOLVER_EQUATIONS%SOLVER_MATRICES%UPDATE_RHS_VECTOR=.TRUE.
          NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES%RHS_VECTOR)
        ELSE
          CALL FlagError("Solver equations solver mapping is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver equations is not associated",ERR,ERROR,*998)
    ENDIF
        
    EXITS("SOLVER_MATRICES_INITIALISE")
    RETURN
999 CALL SOLVER_MATRICES_FINALISE(SOLVER_EQUATIONS%SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("SOLVER_MATRICES_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_INITIALISE
  
  !
  !================================================================================================================================
  !
  
  !>Gets the library type for the solver matrices (and vectors)
  SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_GET(SOLVER_MATRICES,LIBRARY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices.
    INTEGER(INTG), INTENT(OUT) :: LIBRARY_TYPE !<On return, the library type of the specified solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("SOLVER_MATRICES_LIBRARY_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        LIBRARY_TYPE=SOLVER_MATRICES%LIBRARY_TYPE
      ELSE
        CALL FlagError("Solver matrices has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MATRICES_LIBRARY_TYPE_GET")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_LIBRARY_TYPE_GET",ERR,ERROR)
    RETURN
  END SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_GET
  
        
  !
  !================================================================================================================================
  !

  !>Sets the library type for the solver matrices (and vectors)
  SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,LIBRARY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices.
    INTEGER(INTG), INTENT(IN) :: LIBRARY_TYPE !<The library type to set \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("SOLVER_MATRICES_LIBRARY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL FlagError("Solver matrices has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(LIBRARY_TYPE)
        CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
          SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE
        CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
          SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
        CASE DEFAULT
          LOCAL_ERROR="The library type of "// TRIM(NUMBER_TO_VSTRING(LIBRARY_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MATRICES_LIBRARY_TYPE_SET")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_LIBRARY_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Outputs the solver matrices
  SUBROUTINE SOLVER_MATRICES_OUTPUT(ID,SELECTION_TYPE,SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The type of matrix selection \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    
    ENTERS("SOLVER_MATRICES_OUTPUT",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"",ERR,ERROR,*999)
        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
!           & SELECTION_TYPE==SOLVER_MATRICES_DYNAMIC_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
          CALL WRITE_STRING(ID,"Solver matrices:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(ID,"Number of matrices = ",SOLVER_MATRICES%NUMBER_OF_MATRICES,ERR,ERROR,*999)
          DO matrixIdx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrixIdx)%PTR
            IF(ASSOCIATED(SOLVER_MATRIX)) THEN
              CALL WRITE_STRING_VALUE(ID,"Solver matrix : ",matrixIdx,ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrixIdx
        ENDIF
        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
          IF(ASSOCIATED(SOLVER_MATRICES%RESIDUAL)) THEN
            CALL WRITE_STRING(ID,"Solver residual vector:",ERR,ERROR,*999)     
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)  
          ENDIF
        ENDIF
        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
!          & SELECTION_TYPE==SOLVER_MATRICES_DYNAMIC_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RHS_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
          IF(ASSOCIATED(SOLVER_MATRICES%RHS_VECTOR)) THEN
            CALL WRITE_STRING(ID,"Solver RHS vector:",ERR,ERROR,*999)     
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Solver matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MATRICES_OUTPUT")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_OUTPUT",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_OUTPUT
  
  !
  !================================================================================================================================
  !
  
  !>Gets the storage type (sparsity) of the solver matrices
  SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_GET(SOLVER_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrixIdx). On return, the storage type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("SOLVER_MATRICES_STORAGE_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        IF(SIZE(STORAGE_TYPE,1)>=SOLVER_MATRICES%NUMBER_OF_MATRICES) THEN
          DO matrixIdx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrixIdx)%PTR
            IF(ASSOCIATED(SOLVER_MATRIX)) THEN
              STORAGE_TYPE(matrixIdx)=SOLVER_MATRIX%STORAGE_TYPE
            ELSE
              CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrixIdx
        ELSE
          LOCAL_ERROR="The size of STORAGE_TYPE is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver matrices have not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MATRICES_STORAGE_TYPE_GET")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_STORAGE_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the solver matrices
  SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrixIdx). The storage type for the matrixIdx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("SOLVER_MATRICES_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL FlagError("Solver matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STORAGE_TYPE,1)==SOLVER_MATRICES%NUMBER_OF_MATRICES) THEN
          DO matrixIdx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrixIdx)%PTR
            IF(ASSOCIATED(SOLVER_MATRIX)) THEN
              SELECT CASE(STORAGE_TYPE(matrixIdx))
              CASE(MATRIX_BLOCK_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
              CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_DIAGONAL_STORAGE_TYPE        
              CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_COLUMN_MAJOR_STORAGE_TYPE
              CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_ROW_MAJOR_STORAGE_TYPE
              CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
              CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
              CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_ROW_COLUMN_STORAGE_TYPE
              CASE DEFAULT
                LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrixIdx),"*",ERR,ERROR))// &
                  & " for the matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrixIdx
        ELSE
          LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MATRICES_STORAGE_TYPE_SET")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_STORAGE_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Adds alpha times the equations matrix into the solver matrix
  SUBROUTINE SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equationsSetIdx,alpha,EQUATIONS_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping that contains the equations matrix to add
    REAL(DP), INTENT(IN) :: alpha !<The multiplicative factor for the equations matrix
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX !<A pointer to the equations matrix to add    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsColumnIdx,equationsColumnNumber,equationsRowNumber,equationsStorageType, &
      & solverColumnNumber,solverRowNumber,rowConstraintRowNumber,rowConstraintColumnNumber,rowConstraintColumnIdx, &
      & columnConstraintRowNumber,columnConstraintColumnNumber,columnConstraintColumnIdx
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:),columnConstraintRowIndices(:),columnConstraintColumnIndices(:), &
      & rowConstraintRowIndices(:),rowConstraintColumnIndices(:)
    REAL(DP) :: VALUE
    REAL(DP), POINTER :: EQUATIONS_MATRIX_DATA(:),rowConstraintMatrixData(:),columnConstraintMatrixData(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX, &
      & rowConstraintMatrix,columnConstraintMatrix
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE), POINTER :: rowSolverVariableMap,columnSolverVariableMap
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("SOLVER_MATRIX_EQUATIONS_MATRIX_ADD",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(SOLVER_MATRIX)) CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_MATRIX)) CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
    IF(ABS(alpha)<ZERO_TOLERANCE) RETURN
    IF(.NOT.(equationsSetIdx>0.AND.equationsSetIdx<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)) THEN
      LOCAL_ERROR="The specified equations set index of "// &
        & TRIM(NUMBER_TO_VSTRING(equationsSetIdx,"*",ERR,ERROR))// &
        & " is invalid. The equations set index needs to be between 1 and "// &
        & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,"*",ERR,ERROR))//"."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    SOLVER_MATRICES=>SOLVER_MATRIX%SOLVER_MATRICES
    IF(.NOT.ASSOCIATED(SOLVER_MATRICES)) CALL FlagError("Solver matrix solver matrices is not associated.",ERR,ERROR,*999)
    IF(.NOT.SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) CALL FlagError("Solver matrices have not been finished.",ERR,ERROR,*999)
    SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(SOLVER_MAPPING)) CALL FlagError("Solver matrices solver mapping is not associated.",ERR,ERROR,*999)
    LINEAR_MATRICES=>EQUATIONS_MATRIX%LINEAR_MATRICES
    DYNAMIC_MATRICES=>EQUATIONS_MATRIX%DYNAMIC_MATRICES
    IF(.NOT.(ASSOCIATED(DYNAMIC_MATRICES).OR.ASSOCIATED(LINEAR_MATRICES))) &
      & CALL FlagError("Equations matrix dynamic or linear matrices is not associated.",ERR,ERROR,*999)
    rowSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%solverVariableMap
    IF(.NOT.ASSOCIATED(rowSolverVariableMap)) CALL FlagError("Row solver variable map is not associated.",ERR,ERROR,*999)
    IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
      EQUATIONS_MATRICES=>DYNAMIC_MATRICES%EQUATIONS_MATRICES
      columnSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)% &
        & equationsToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
        & dynamicMatrixToSolverVariableMap(EQUATIONS_MATRIX%MATRIX_NUMBER)%ptr
    ELSE
      EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
      columnSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)% &
        & equationsToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
        & linearMatrixToSolverVariableMap(EQUATIONS_MATRIX%MATRIX_NUMBER)%ptr
    ENDIF
    IF(.NOT.ASSOCIATED(columnSolverVariableMap)) CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
    rowConstraintMatrix=>rowSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
    NULLIFY(rowConstraintMatrixData)
    NULLIFY(rowConstraintRowIndices)
    NULLIFY(rowConstraintColumnIndices)
    CALL DISTRIBUTED_MATRIX_DATA_GET(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(rowConstraintMatrix,rowConstraintRowIndices, &
      & rowConstraintColumnIndices,ERR,ERROR,*999)
    columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
    NULLIFY(columnConstraintMatrixData)
    NULLIFY(columnConstraintRowIndices)
    NULLIFY(columnConstraintColumnIndices)
    CALL DISTRIBUTED_MATRIX_DATA_GET(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix,columnConstraintRowIndices, &
      & columnConstraintColumnIndices,ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(EQUATIONS_MATRICES)) &
      & CALL FlagError("Dynamic or linear matrices equations matrices is not associated.",ERR,ERROR,*999)
    IF(.NOT.EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) &
      & CALL FlagError("Equations matrices have not been finished.",ERR,ERROR,*999)
    SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
    IF(.NOT.ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) &
      & CALL FlagError("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
    EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
    IF(.NOT.ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) &
      & CALL FlagError("The equations matrix distributed matrix is not associated",ERR,ERROR,*999)
    NULLIFY(EQUATIONS_MATRIX_DATA)
    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX,equationsStorageType,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA,ERR,ERROR,*999)
      
    SELECT CASE(equationsStorageType)
    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
      !Loop over the rows of the equations matrix
      DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
        IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
          solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
          !Loop over the columns of the equations matrix
          DO equationsColumnNumber=1,EQUATIONS_MATRIX%TOTAL_NUMBER_OF_COLUMNS
            IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
              VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsRowNumber+ &
                & (equationsColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            ELSE
              !The column dof is constrained
              columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
              DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsRowNumber+ &
                  & (equationsColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                  & columnConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              END DO !constraintColumnIdx
            END IF
          ENDDO !equationsColumnIdx
        ELSE
          !The row dof is constrained
          rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
          DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
            & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
            rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
            !Loop over the columns of the equations matrix
            DO equationsColumnNumber=1,EQUATIONS_MATRIX%TOTAL_NUMBER_OF_COLUMNS
              IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsRowNumber+ &
                  & (equationsColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                  & rowConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsRowNumber+ &
                    & (equationsColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                    & rowConstraintMatrixData(columnConstraintColumnIdx)* &
                    & columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !columnConstraintColumnIdx
              END IF
            END DO !equationsColumnIdx
          END DO !rowConstraintColumnIdx
        END IF
      END DO !equationsRowNumber
    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
      !Loop over the rows of the equations matrix
      DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
        equationsColumnNumber=equationsRowNumber
        IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
          solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
          IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
            solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
            VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsRowNumber)
            CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
              & ERR,ERROR,*999)
          ELSE
            !The column dof is constrained
            columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
            DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
              & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
              columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
              VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsRowNumber)*columnConstraintMatrixData(columnConstraintColumnIdx)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            END DO !constraintColumnIdx
          END IF
        ELSE
          !The row dof is constrained
          rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
          DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
            & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
            rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
            IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
              VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            ELSE
              !The column dof is constrained
              columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
              DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                  & columnConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              END DO !columnConstraintColumnIdx
            END IF
          END DO !rowConstraintColumnIdx
        END IF
      END DO !equationsRowNumber
    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      NULLIFY(columnIndices)
      NULLIFY(rowIndices)
      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX,rowIndices,columnIndices,ERR,ERROR,*999)
      !Loop over the rows of the equations matrix
      DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
        IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
          solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
          !Loop over the columns of the equations matrix
          DO equationsColumnIdx=rowIndices(equationsRowNumber),rowIndices(equationsRowNumber+1)-1
            equationsColumnNumber=columnIndices(equationsColumnIdx)
            IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
              VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsColumnIdx)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            ELSE
              !The column dof is constrained
              columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
              DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsColumnIdx)*columnConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              END DO !constraintColumnIdx
            END IF
          ENDDO !equationsColumnIdx
        ELSE
          !The row dof is constrained
          rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
          DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
            & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
            rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
            !Loop over the columns of the equations matrix
            DO equationsColumnIdx=rowIndices(equationsRowNumber),rowIndices(equationsRowNumber+1)-1
              equationsColumnNumber=columnIndices(equationsColumnIdx)
              IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha*EQUATIONS_MATRIX_DATA(equationsColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                    & columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !columnConstraintColumnIdx
              END IF
            END DO !equationsColumnIdx
          END DO !rowConstraintColumnIdx
        END IF
      END DO !equationsRowNumber
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="The equations matrix storage type of "// &
        & TRIM(NUMBER_TO_VSTRING(equationsStorageType,"*",ERR,ERROR))//" is invalid."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA,ERR,ERROR,*999)
    
    EXITS("SOLVER_MATRIX_EQUATIONS_MATRIX_ADD")
    RETURN
999 ERRORSEXITS("SOLVER_MATRIX_EQUATIONS_MATRIX_ADD",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRIX_EQUATIONS_MATRIX_ADD

  !
  !================================================================================================================================
  !

  !>Adds alpha times the interface matrix into the solver matrix
  SUBROUTINE SOLVER_MATRIX_INTERFACE_MATRIX_ADD(SOLVER_MATRIX,interfaceConditionIdx,alpha,INTERFACE_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interfaceConditionIdx index in the solver mapping that contains the interface matrix to add
    REAL(DP), INTENT(IN) :: alpha(2) !<The multiplicative factor for the interface matrix
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX !<A pointer to the interface matrix to add    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceColumnIdx,interfaceColumnNumber,interfaceRowNumber,interfaceStorageType, &
      & solverColumnNumber,solverRowNumber,rowConstraintRowNumber,rowConstraintColumnNumber,rowConstraintColumnIdx, &
      & columnConstraintRowNumber,columnConstraintColumnNumber,columnConstraintColumnIdx, &
      & interfaceMatrixNumberOfRows,interfaceMatrixNumberOfColumns
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:),columnConstraintRowIndices(:),columnConstraintColumnIndices(:), &
      & rowConstraintRowIndices(:),rowConstraintColumnIndices(:)
    REAL(DP) :: VALUE
    REAL(DP), POINTER :: INTERFACE_MATRIX_DATA(:),rowConstraintMatrixData(:),columnConstraintMatrixData(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: INTERFACE_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX, &
      & rowConstraintMatrix,columnConstraintMatrix
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE), POINTER :: rowSolverVariableMap,columnSolverVariableMap
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("SOLVER_MATRIX_INTERFACE_MATRIX_ADD",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(SOLVER_MATRIX)) CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(INTERFACE_MATRIX)) CALL FlagError("Interface matrix is not associated.",ERR,ERROR,*999)
    IF(.NOT.(interfaceConditionIdx>0.AND.interfaceConditionIdx<=SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS)) THEN
      LOCAL_ERROR="The specified interface set index of "// &
        & TRIM(NUMBER_TO_VSTRING(interfaceConditionIdx,"*",ERR,ERROR))// &
        & " is invalid. The interface set index needs to be between 1 and "// &
        & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS,"*",ERR,ERROR))//"."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    SOLVER_MATRICES=>SOLVER_MATRIX%SOLVER_MATRICES
    IF(.NOT.ASSOCIATED(SOLVER_MATRICES)) CALL FlagError("Solver matrix solver matrices is not associated.",ERR,ERROR,*999)
    IF(.NOT.SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) CALL FlagError("Solver matrices have not been finished.",ERR,ERROR,*999)
    SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(SOLVER_MAPPING)) CALL FlagError("Solver matrices solver mapping is not associated.",ERR,ERROR,*999)
    INTERFACE_MATRICES=>INTERFACE_MATRIX%INTERFACE_MATRICES
    IF(.NOT.(ASSOCIATED(INTERFACE_MATRICES))) &
      & CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
    IF(.NOT.INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) &
      & CALL FlagError("Interface matrices have not been finished.",ERR,ERROR,*999)
    SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
    IF(.NOT.ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) &
      & CALL FlagError("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
      
    IF(ABS(alpha(1))>ZERO_TOLERANCE) THEN
      rowSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)% &
        & interfaceToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
        & interfaceMatrixToSolverVariableMap(INTERFACE_MATRIX%MATRIX_NUMBER)%ptr
      IF(.NOT.ASSOCIATED(rowSolverVariableMap)) CALL FlagError("Row solver variable map is not associated.",ERR,ERROR,*999)
      !The Lagrange variable
      columnSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%solverVariableMap
      IF(.NOT.ASSOCIATED(columnSolverVariableMap)) CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
      rowConstraintMatrix=>rowSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
      NULLIFY(rowConstraintMatrixData)
      NULLIFY(rowConstraintRowIndices)
      NULLIFY(rowConstraintColumnIndices)
      CALL DISTRIBUTED_MATRIX_DATA_GET(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(rowConstraintMatrix,rowConstraintRowIndices, &
        & rowConstraintColumnIndices,ERR,ERROR,*999)
      columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
      NULLIFY(columnConstraintMatrixData)
      NULLIFY(columnConstraintRowIndices)
      NULLIFY(columnConstraintColumnIndices)
      CALL DISTRIBUTED_MATRIX_DATA_GET(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix,columnConstraintRowIndices, &
      & columnConstraintColumnIndices,ERR,ERROR,*999)
      INTERFACE_DISTRIBUTED_MATRIX=>INTERFACE_MATRIX%MATRIX
      IF(.NOT.ASSOCIATED(INTERFACE_DISTRIBUTED_MATRIX)) &
        & CALL FlagError("The interface matrix distributed matrix is not associated",ERR,ERROR,*999)
      NULLIFY(INTERFACE_MATRIX_DATA)
      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(INTERFACE_DISTRIBUTED_MATRIX,interfaceStorageType,ERR,ERROR,*999)
      CALL DISTRIBUTED_MATRIX_DATA_GET(INTERFACE_DISTRIBUTED_MATRIX,INTERFACE_MATRIX_DATA,ERR,ERROR,*999)
      interfaceMatrixNumberOfRows=INTERFACE_DISTRIBUTED_MATRIX%ROW_DOMAIN_MAPPING%NUMBER_OF_LOCAL
      interfaceMatrixNumberOfColumns=INTERFACE_DISTRIBUTED_MATRIX%COLUMN_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
      SELECT CASE(interfaceStorageType)
      CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
        !Loop over the rows of the interface matrix
        DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
          IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
            !Loop over the columns of the interface matrix
            DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
              IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                  & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                    & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)* &
                    & columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !constraintColumnIdx
              END IF
            ENDDO !interfaceColumnIdx
          ELSE
            !The row dof is constrained
            rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
            DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
              & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
              rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
              solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
              !Loop over the columns of the interface matrix
              DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
                IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                    & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)* &
                    & rowConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                ELSE
                  !The column dof is constrained
                  columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                    & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                    columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                    VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                      & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)* &
                      & rowConstraintMatrixData(columnConstraintColumnIdx)* &
                      & columnConstraintMatrixData(columnConstraintColumnIdx)
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                      & ERR,ERROR,*999)
                  END DO !columnConstraintColumnIdx
                END IF
              END DO !interfaceColumnIdx
            END DO !rowConstraintColumnIdx
          END IF
        END DO !interfaceRowNumber
      CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
        !Loop over the rows of the interface matrix
        DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
          interfaceColumnNumber=interfaceRowNumber
          IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
            IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
              VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceRowNumber)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            ELSE
              !The column dof is constrained
              columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
              DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceRowNumber)*columnConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              END DO !constraintColumnIdx
            END IF
          ELSE
            !The row dof is constrained
            rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
            DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
              & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
              rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
              solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
              IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                    & columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !columnConstraintColumnIdx
              END IF
            END DO !rowConstraintColumnIdx
          END IF
        END DO !interfaceRowNumber
      CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
        NULLIFY(columnIndices)
        NULLIFY(rowIndices)
        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(INTERFACE_DISTRIBUTED_MATRIX,rowIndices,columnIndices,ERR,ERROR,*999)
        !Loop over the rows of the interface matrix
        DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
          IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
            !Loop over the columns of the interface matrix
            DO interfaceColumnIdx=rowIndices(interfaceRowNumber),rowIndices(interfaceRowNumber+1)-1
              interfaceColumnNumber=columnIndices(interfaceColumnIdx)
              IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !constraintColumnIdx
              END IF
            ENDDO !interfaceColumnIdx
          ELSE
            !The row dof is constrained
            rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
            DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
              & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
              rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
              solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
              !Loop over the columns of the interface matrix
              DO interfaceColumnIdx=rowIndices(interfaceRowNumber),rowIndices(interfaceRowNumber+1)-1
                interfaceColumnNumber=columnIndices(interfaceColumnIdx)
                IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                ELSE
                  !The column dof is constrained
                  columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                    & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                    columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                    VALUE=alpha(1)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                      & columnConstraintMatrixData(columnConstraintColumnIdx)
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                      & ERR,ERROR,*999)
                  END DO !columnConstraintColumnIdx
                END IF
              END DO !interfaceColumnIdx
            END DO !rowConstraintColumnIdx
          END IF
        END DO !interfaceRowNumber
      CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The interface matrix storage type of "// &
          & TRIM(NUMBER_TO_VSTRING(interfaceStorageType,"*",ERR,ERROR))//" is invalid."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      CALL DISTRIBUTED_MATRIX_DATA_RESTORE(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
      CALL DISTRIBUTED_MATRIX_DATA_RESTORE(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
      CALL DISTRIBUTED_MATRIX_DATA_RESTORE(INTERFACE_DISTRIBUTED_MATRIX,INTERFACE_MATRIX_DATA,ERR,ERROR,*999)
    END IF
    IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
      IF(ABS(alpha(2))>ZERO_TOLERANCE) THEN
        !The Lagrange variable
        rowSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%solverVariableMap
        IF(.NOT.ASSOCIATED(rowSolverVariableMap)) CALL FlagError("Row solver variable map is not associated.",ERR,ERROR,*999)
        columnSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)% &
          & interfaceToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
          & interfaceMatrixToSolverVariableMap(INTERFACE_MATRIX%MATRIX_NUMBER)%ptr
        IF(.NOT.ASSOCIATED(columnSolverVariableMap)) CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
        rowConstraintMatrix=>rowSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
        NULLIFY(rowConstraintMatrixData)
        NULLIFY(rowConstraintRowIndices)
        NULLIFY(rowConstraintColumnIndices)
        CALL DISTRIBUTED_MATRIX_DATA_GET(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(rowConstraintMatrix,rowConstraintRowIndices, &
          & rowConstraintColumnIndices,ERR,ERROR,*999)
        columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
        NULLIFY(columnConstraintMatrixData)
        NULLIFY(columnConstraintRowIndices)
        NULLIFY(columnConstraintColumnIndices)
        CALL DISTRIBUTED_MATRIX_DATA_GET(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix,columnConstraintRowIndices, &
        & columnConstraintColumnIndices,ERR,ERROR,*999)
        INTERFACE_DISTRIBUTED_MATRIX=>INTERFACE_MATRIX%MATRIX_TRANSPOSE
        IF(.NOT.ASSOCIATED(INTERFACE_DISTRIBUTED_MATRIX)) &
          & CALL FlagError("The interface matrix distributed matrix is not associated",ERR,ERROR,*999)
        NULLIFY(INTERFACE_MATRIX_DATA)
        CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(INTERFACE_DISTRIBUTED_MATRIX,interfaceStorageType,ERR,ERROR,*999)
        CALL DISTRIBUTED_MATRIX_DATA_GET(INTERFACE_DISTRIBUTED_MATRIX,INTERFACE_MATRIX_DATA,ERR,ERROR,*999)
        interfaceMatrixNumberOfRows=INTERFACE_DISTRIBUTED_MATRIX%ROW_DOMAIN_MAPPING%NUMBER_OF_LOCAL
        interfaceMatrixNumberOfColumns=INTERFACE_DISTRIBUTED_MATRIX%COLUMN_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
        SELECT CASE(interfaceStorageType)
        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
          !Loop over the rows of the interface matrix
          DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
            IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
              solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
              !Loop over the columns of the interface matrix
              DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
                IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                    & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                ELSE
                  !The column dof is constrained
                  columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                    & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                    columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                    VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                      & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)* &
                      & columnConstraintMatrixData(columnConstraintColumnIdx)
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                      & ERR,ERROR,*999)
                  END DO !constraintColumnIdx
                END IF
              ENDDO !interfaceColumnIdx
            ELSE
              !The row dof is constrained
              rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
              DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                !Loop over the columns of the interface matrix
                DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
                  IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                      & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)* &
                      & rowConstraintMatrixData(columnConstraintColumnIdx)
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                      & ERR,ERROR,*999)
                  ELSE
                    !The column dof is constrained
                    columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                      & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                      columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                      VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceRowNumber+ &
                        & (interfaceColumnNumber-1)*interfaceMatrixNumberOfRows)* &
                        & rowConstraintMatrixData(columnConstraintColumnIdx)* &
                        & columnConstraintMatrixData(columnConstraintColumnIdx)
                      CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                        & ERR,ERROR,*999)
                    END DO !columnConstraintColumnIdx
                  END IF
                END DO !interfaceColumnIdx
              END DO !rowConstraintColumnIdx
            END IF
          END DO !interfaceRowNumber
        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
          !Loop over the rows of the interface matrix
          DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
            interfaceColumnNumber=interfaceRowNumber
            IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
              solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
              IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceRowNumber)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceRowNumber)*columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !constraintColumnIdx
              END IF
            ELSE
              !The row dof is constrained
              rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
              DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                ELSE
                  !The column dof is constrained
                  columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                    & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                    columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                    VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                      & columnConstraintMatrixData(columnConstraintColumnIdx)
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                      & ERR,ERROR,*999)
                  END DO !columnConstraintColumnIdx
                END IF
              END DO !rowConstraintColumnIdx
            END IF
          END DO !interfaceRowNumber
        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          NULLIFY(columnIndices)
          NULLIFY(rowIndices)
          CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(INTERFACE_DISTRIBUTED_MATRIX,rowIndices,columnIndices,ERR,ERROR,*999)
          !Loop over the rows of the interface matrix
          DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
            IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
              solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
              !Loop over the columns of the interface matrix
              DO interfaceColumnIdx=rowIndices(interfaceRowNumber),rowIndices(interfaceRowNumber+1)-1
                interfaceColumnNumber=columnIndices(interfaceColumnIdx)
                IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                ELSE
                  !The column dof is constrained
                  columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                  DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                    & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                    columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                    VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*columnConstraintMatrixData(columnConstraintColumnIdx)
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                      & ERR,ERROR,*999)
                  END DO !constraintColumnIdx
                END IF
              ENDDO !interfaceColumnIdx
            ELSE
              !The row dof is constrained
              rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
              DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                !Loop over the columns of the interface matrix
                DO interfaceColumnIdx=rowIndices(interfaceRowNumber),rowIndices(interfaceRowNumber+1)-1
                  interfaceColumnNumber=columnIndices(interfaceColumnIdx)
                  IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
                    CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                      & ERR,ERROR,*999)
                  ELSE
                    !The column dof is constrained
                    columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                      & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                      columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                      VALUE=alpha(2)*INTERFACE_MATRIX_DATA(interfaceColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                        & columnConstraintMatrixData(columnConstraintColumnIdx)
                      CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                        & ERR,ERROR,*999)
                    END DO !columnConstraintColumnIdx
                  END IF
                END DO !interfaceColumnIdx
              END DO !rowConstraintColumnIdx
            END IF
          END DO !interfaceRowNumber
        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The interface matrix storage type of "// &
            & TRIM(NUMBER_TO_VSTRING(interfaceStorageType,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        CALL DISTRIBUTED_MATRIX_DATA_RESTORE(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
        CALL DISTRIBUTED_MATRIX_DATA_RESTORE(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
        CALL DISTRIBUTED_MATRIX_DATA_RESTORE(INTERFACE_DISTRIBUTED_MATRIX,INTERFACE_MATRIX_DATA,ERR,ERROR,*999)
      END IF
    END IF
    
    EXITS("SOLVER_MATRIX_INTERFACE_MATRIX_ADD")
    RETURN
999 ERRORSEXITS("SOLVER_MATRIX_INTERFACE_MATRIX_ADD",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRIX_INTERFACE_MATRIX_ADD

  !
  !================================================================================================================================
  !

  !>Adds alpha times the Jacobian matrix into the solver matrix
  SUBROUTINE SOLVER_MATRIX_JACOBIAN_MATRIX_ADD(SOLVER_MATRIX,equationsSetIdx,alpha,JACOBIAN_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping that contains the Jacobian matrix to add
    REAL(DP), INTENT(IN) :: alpha !<The multiplicative factor for the Jacobian matrix
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX !<A pointer to the Jacobian matrix to add    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianColumnIdx,jacobianColumnNumber,jacobianRowNumber,JACOBIAN_STORAGE_TYPE, &
      & solverColumnNumber,solverRowNumber,rowConstraintRowNumber,rowConstraintColumnNumber,rowConstraintColumnIdx, &
      & columnConstraintRowNumber,columnConstraintColumnNumber,columnConstraintColumnIdx
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:),columnConstraintRowIndices(:),columnConstraintColumnIndices(:), &
      & rowConstraintRowIndices(:),rowConstraintColumnIndices(:)
    REAL(DP) :: VALUE
    REAL(DP), POINTER :: JACOBIAN_MATRIX_DATA(:),rowConstraintMatrixData(:),columnConstraintMatrixData(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: JACOBIAN_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX, &
      & rowConstraintMatrix,columnConstraintMatrix
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE), POINTER :: rowSolverVariableMap,columnSolverVariableMap
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("SOLVER_MATRIX_JACOBIAN_MATRIX_ADD",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(SOLVER_MATRIX)) CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(JACOBIAN_MATRIX)) CALL FlagError("Jacobian matrix is not associated.",ERR,ERROR,*999)
    IF(ABS(alpha)<ZERO_TOLERANCE) RETURN
    IF(.NOT.(equationsSetIdx>0.AND.equationsSetIdx<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)) THEN
      LOCAL_ERROR="The specified equations set index of "// &
        & TRIM(NUMBER_TO_VSTRING(equationsSetIdx,"*",ERR,ERROR))// &
        & " is invalid. The equations set index needs to be between 1 and "// &
        & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,"*",ERR,ERROR))//"."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF
    SOLVER_MATRICES=>SOLVER_MATRIX%SOLVER_MATRICES
    IF(.NOT.ASSOCIATED(SOLVER_MATRICES)) CALL FlagError("Solver matrix solver matrices is not associated.",ERR,ERROR,*999)
    IF(.NOT.SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) CALL FlagError("Solver matrices have not been finished.",ERR,ERROR,*999)
    SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(SOLVER_MAPPING)) CALL FlagError("Solver matrices solver mapping is not associated.",ERR,ERROR,*999)
    NONLINEAR_MATRICES=>JACOBIAN_MATRIX%NONLINEAR_MATRICES
    IF(.NOT.(ASSOCIATED(NONLINEAR_MATRICES))) &
      & CALL FlagError("Jacobian matrix nonlinear matrices is not associated.",ERR,ERROR,*999)
    rowSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%solverVariableMap
    IF(.NOT.ASSOCIATED(rowSolverVariableMap)) CALL FlagError("Row solver variable map is not associated.",ERR,ERROR,*999)
    columnSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)% &
      & equationsToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
      & linearMatrixToSolverVariableMap(JACOBIAN_MATRIX%JACOBIAN_NUMBER)%ptr
    IF(.NOT.ASSOCIATED(columnSolverVariableMap)) CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
    rowConstraintMatrix=>rowSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
    NULLIFY(rowConstraintMatrixData)
    NULLIFY(rowConstraintRowIndices)
    NULLIFY(rowConstraintColumnIndices)
    CALL DISTRIBUTED_MATRIX_DATA_GET(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(rowConstraintMatrix,rowConstraintRowIndices, &
      & rowConstraintColumnIndices,ERR,ERROR,*999)
    columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
    NULLIFY(columnConstraintMatrixData)
    NULLIFY(columnConstraintRowIndices)
    NULLIFY(columnConstraintColumnIndices)
    CALL DISTRIBUTED_MATRIX_DATA_GET(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix,columnConstraintRowIndices, &
      & columnConstraintColumnIndices,ERR,ERROR,*999)
    EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
    IF(.NOT.ASSOCIATED(EQUATIONS_MATRICES)) &
      & CALL FlagError("Dynamic or linear matrices equations matrices is not associated.",ERR,ERROR,*999)
    IF(.NOT.EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) &
      & CALL FlagError("Equations matrices have not been finished.",ERR,ERROR,*999)
    SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
    IF(.NOT.ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) &
      & CALL FlagError("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
    JACOBIAN_DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
    IF(.NOT.ASSOCIATED(JACOBIAN_DISTRIBUTED_MATRIX)) &
      & CALL FlagError("The jacobian matrix distributed matrix is not associated",ERR,ERROR,*999)
    NULLIFY(JACOBIAN_MATRIX_DATA)
    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_STORAGE_TYPE,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_DATA_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA,ERR,ERROR,*999)
      
    SELECT CASE(JACOBIAN_STORAGE_TYPE)
    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
      !Loop over the rows of the jacobian matrix
      DO jacobianRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
        IF(rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)>0) THEN
          solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)
          !Loop over the columns of the jacobian matrix
          DO jacobianColumnNumber=1,JACOBIAN_MATRIX%TOTAL_NUMBER_OF_COLUMNS
            IF(columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)>0) THEN
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
              VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianRowNumber+ &
                & (jacobianColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            ELSE
              !The column dof is constrained
              columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
              DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianRowNumber+ &
                  & (jacobianColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                  & columnConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              END DO !constraintColumnIdx
            END IF
          ENDDO !jacobianColumnIdx
        ELSE
          !The row dof is constrained
          rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)
          DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
            & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
            rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
            !Loop over the columns of the jacobian matrix
            DO jacobianColumnNumber=1,JACOBIAN_MATRIX%TOTAL_NUMBER_OF_COLUMNS
              IF(columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
                VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianRowNumber+ &
                  & (jacobianColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                  & rowConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianRowNumber+ &
                    & (jacobianColumnNumber-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                    & rowConstraintMatrixData(columnConstraintColumnIdx)* &
                    & columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !columnConstraintColumnIdx
              END IF
            END DO !jacobianColumnIdx
          END DO !rowConstraintColumnIdx
        END IF
      END DO !jacobianRowNumber
    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
      !Loop over the rows of the jacobian matrix
      DO jacobianRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
        jacobianColumnNumber=jacobianRowNumber
        IF(rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)>0) THEN
          solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)
          IF(columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)>0) THEN
            solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
            VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianRowNumber)
            CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
              & ERR,ERROR,*999)
          ELSE
            !The column dof is constrained
            columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
            DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
              & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
              columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
              VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianRowNumber)*columnConstraintMatrixData(columnConstraintColumnIdx)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            END DO !constraintColumnIdx
          END IF
        ELSE
          !The row dof is constrained
          rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)
          DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
            & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
            rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
            IF(columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)>0) THEN
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
              VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            ELSE
              !The column dof is constrained
              columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
              DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                  & columnConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              END DO !columnConstraintColumnIdx
            END IF
          END DO !rowConstraintColumnIdx
        END IF
      END DO !jacobianRowNumber
    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      NULLIFY(columnIndices)
      NULLIFY(rowIndices)
      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(JACOBIAN_DISTRIBUTED_MATRIX,rowIndices,columnIndices,ERR,ERROR,*999)
      !Loop over the rows of the jacobian matrix
      DO jacobianRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
        IF(rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)>0) THEN
          solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)
          !Loop over the columns of the jacobian matrix
          DO jacobianColumnIdx=rowIndices(jacobianRowNumber),rowIndices(jacobianRowNumber+1)-1
            jacobianColumnNumber=columnIndices(jacobianColumnIdx)
            IF(columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)>0) THEN
              solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
              VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianColumnIdx)
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                & ERR,ERROR,*999)
            ELSE
              !The column dof is constrained
              columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
              DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianColumnIdx)*columnConstraintMatrixData(columnConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              END DO !constraintColumnIdx
            END IF
          ENDDO !jacobianColumnIdx
        ELSE
          !The row dof is constrained
          rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(jacobianRowNumber)
          DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
            & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
            rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
            solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
            !Loop over the columns of the jacobian matrix
            DO jacobianColumnIdx=rowIndices(jacobianRowNumber),rowIndices(jacobianRowNumber+1)-1
              jacobianColumnNumber=columnIndices(jacobianColumnIdx)
              IF(columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)>0) THEN
                solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
                VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                  & ERR,ERROR,*999)
              ELSE
                !The column dof is constrained
                columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(jacobianColumnNumber)
                DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                  & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                  columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                  solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                  VALUE=alpha*JACOBIAN_MATRIX_DATA(jacobianColumnIdx)*rowConstraintMatrixData(rowConstraintColumnIdx)* &
                    & columnConstraintMatrixData(columnConstraintColumnIdx)
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solverRowNumber,solverColumnNumber,VALUE, &
                    & ERR,ERROR,*999)
                END DO !columnConstraintColumnIdx
              END IF
            END DO !jacobianColumnIdx
          END DO !rowConstraintColumnIdx
        END IF
      END DO !jacobianRowNumber
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="The jacobian matrix storage type of "// &
        & TRIM(NUMBER_TO_VSTRING(JACOBIAN_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(rowConstraintMatrix,rowConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(columnConstraintMatrix,columnConstraintMatrixData,ERR,ERROR,*999)
    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA,ERR,ERROR,*999)
    
    EXITS("SOLVER_MATRIX_JACOBIAN_MATRIX_ADD")
    RETURN
999 ERRORSEXITS("SOLVER_MATRIX_JACOBIAN_MATRIX_ADD",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRIX_JACOBIAN_MATRIX_ADD

  !
  !================================================================================================================================
  !

  !>Calculates the structure (sparsity) of the solver matrix from the solution mapping.
  SUBROUTINE SOLVER_MATRIX_STRUCTURE_CALCULATE(SOLVER_MATRIX,numberOfNonZeros,rowIndices,columnIndices,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the solver matrix
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: rowIndices(:) !<On return the row location indices in compressed row format. The calling routine is responsible for deallocation.
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: columnIndices(:) !<On return the column location indices in compressed row format. The calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsColumnIdx,equationsColumnNumber,DUMMY_ERR,equationsMatrixIdx,equationsRowNumber, &
      & equationsSetIdx,equationsStorageType,interfaceColumnIdx,interfaceColumnNumber,interfaceConditionIdx, &
      & interfaceMatrixIdx,interfaceMatrixNumberOfRows,interfaceMatrixNumberOfColumns,interfaceRowNumber,interfaceStorageType, &
      & maxColumnIndices,maxEquationsColumnIndices,maxColumnsPerRow,maxTransposeColumnsPerRow, &
      & numberOfColumns,solverColumnNumber,solverMatrixIdx,solverRowNumber, &
      & rowConstraintRowNumber,rowConstraintColumnNumber,rowConstraintColumnIdx, &
      & columnConstraintRowNumber,columnConstraintColumnNumber,columnConstraintColumnIdx
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    INTEGER(INTG), POINTER :: equationsRowIndices(:),equationsColumnIndices(:),interfaceRowIndices(:),interfaceColumnIndices(:), &
      & columnConstraintRowIndices(:),columnConstraintColumnIndices(:),rowConstraintRowIndices(:),rowConstraintColumnIndices(:)
    REAL(DP) :: SPARSITY
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX, &
      & rowConstraintMatrix,columnConstraintMatrix
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES    
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MAPPING_VARIABLE_TYPE), POINTER :: rowSolverVariableMap,columnSolverVariableMap
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    ENTERS("SOLVER_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(SOLVER_MATRIX)) CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
    SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
    IF(.NOT.ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) &
      & CALL FlagError("Solver matrix distributed matrix is not associated",ERR,ERROR,*999)
    IF(SOLVER_DISTRIBUTED_MATRIX%MATRIX_FINISHED)  &
      & CALL FlagError("The solver distributed matrix has already been finished.",ERR,ERROR,*998)
    SOLVER_MATRICES=>SOLVER_MATRIX%SOLVER_MATRICES
    IF(.NOT.ASSOCIATED(SOLVER_MATRICES)) CALL FlagError("Solver matrix solver matrices is not associated",ERR,ERROR,*999)
    SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(SOLVER_MAPPING)) CALL FlagError("Solver matrices solver mapping is not associated",ERR,ERROR,*999)

    SELECT CASE(SOLVER_MATRIX%STORAGE_TYPE)
    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
      CALL FlagError("Can not calcualte the structure for a block storage matrix.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
      CALL FlagError("Can not calcualte the structure for a diagonal matrix.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
      solverMatrixIdx=SOLVER_MATRIX%MATRIX_NUMBER
      !Find the maximum number of column indices
      maxColumnIndices=0
      DO equationsSetIdx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        maxEquationsColumnIndices=0
        EQUATIONS_MATRICES=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%EQUATIONS%EQUATIONS_MATRICES 
        IF(.NOT.ASSOCIATED(EQUATIONS_MATRICES)) CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          DO equationsMatrixIdx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
            EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(equationsMatrixIdx)%ptr
            IF(.NOT.ASSOCIATED(EQUATIONS_MATRIX)) CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
            DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
            CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,maxColumnsPerRow,ERR,ERROR,*999)
            maxEquationsColumnIndices=MAX(maxEquationsColumnIndices,maxColumnsPerRow)
          END DO !equationsMatrixIdx
        END IF
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          DO equationsMatrixIdx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
            EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(equationsMatrixIdx)%ptr
            IF(.NOT.ASSOCIATED(EQUATIONS_MATRIX)) CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
            DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
            CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,maxColumnsPerRow,ERR,ERROR,*999)
            maxEquationsColumnIndices=maxEquationsColumnIndices+maxColumnsPerRow
          END DO !equationsMatrixIdx
        END IF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          DO equationsMatrixIdx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(equationsMatrixIdx)%ptr
            IF(.NOT.ASSOCIATED(JACOBIAN_MATRIX)) CALL FlagError("Equations jacobian matrix is not associated.",ERR,ERROR,*999)
            DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
            CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,maxColumnsPerRow,ERR,ERROR,*999)
            maxEquationsColumnIndices=maxEquationsColumnIndices+maxColumnsPerRow
          END DO !equationsMatrixIdx
        END IF
        maxColumnIndices=MAX(maxColumnIndices,maxEquationsColumnIndices)
      END DO !equationsSetIdx
      DO interfaceConditionIdx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
        INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%INTERFACE_CONDITION
        SELECT CASE(INTERFACE_CONDITION%METHOD)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          INTERFACE_MATRICES=>INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_MATRICES 
          IF(.NOT.ASSOCIATED(INTERFACE_MATRICES)) CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
          DO interfaceMatrixIdx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
            INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(interfaceMatrixIdx)%ptr
            IF(.NOT.ASSOCIATED(INTERFACE_MATRIX)) CALL FlagError("Interface matrix is not associated.",ERR,ERROR,*999)
            DISTRIBUTED_MATRIX=>INTERFACE_MATRIX%MATRIX
            CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,maxColumnsPerRow,ERR,ERROR,*999)
            maxTransposeColumnsPerRow=0
            IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
              DISTRIBUTED_MATRIX=>INTERFACE_MATRIX%MATRIX_TRANSPOSE
              CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,maxTransposeColumnsPerRow,ERR,ERROR,*999)
            ENDIF
            maxColumnIndices=maxColumnIndices+MAX(maxColumnsPerRow,maxTransposeColumnsPerRow)
          END DO !interfaceMatrixIdx
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The interface condition method of "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      END DO !interfaceConditionIdx
      !Allocate lists
      ALLOCATE(columnIndicesLists(SOLVER_MAPPING%NUMBER_OF_DOFS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate column indices lists.",ERR,ERROR,*999)
      !Allocate row indices
      ALLOCATE(rowIndices(SOLVER_MAPPING%NUMBER_OF_DOFS+1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate row indices.",ERR,ERROR,*999)
      !Set up the column indices lists
      DO solverRowNumber=1,SOLVER_MAPPING%NUMBER_OF_DOFS
        NULLIFY(columnIndicesLists(solverRowNumber)%PTR)
        CALL List_CreateStart(columnIndicesLists(solverRowNumber)%PTR,ERR,ERROR,*999)
        CALL List_DataTypeSet(columnIndicesLists(solverRowNumber)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL List_InitialSizeSet(columnIndicesLists(solverRowNumber)%PTR,maxColumnIndices,ERR,ERROR,*999)
        CALL List_CreateFinish(columnIndicesLists(solverRowNumber)%PTR,ERR,ERROR,*999)
      ENDDO !solverRowNumber
      !Loop over the equations sets
      DO equationsSetIdx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
        EQUATIONS_MATRICES=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%EQUATIONS%EQUATIONS_MATRICES 
        rowSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%solverVariableMap
        IF(.NOT.ASSOCIATED(rowSolverVariableMap)) CALL FlagError("Row solver variable map is not associated.",ERR,ERROR,*999)
        rowConstraintMatrix=>rowSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
        NULLIFY(rowConstraintRowIndices)
        NULLIFY(rowConstraintColumnIndices)
        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(rowConstraintMatrix,rowConstraintRowIndices, &
          & rowConstraintColumnIndices,ERR,ERROR,*999)
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          DO equationsMatrixIdx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
            EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(equationsMatrixIdx)%ptr
            DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,equationsStorageType,ERR,ERROR,*999)
            columnSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)% &
              & equationsToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
              & dynamicMatrixToSolverVariableMap(equationsMatrixIdx)%ptr
            IF(.NOT.ASSOCIATED(columnSolverVariableMap)) &
              & CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
            columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
            NULLIFY(columnConstraintRowIndices)
            NULLIFY(columnConstraintColumnIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix, &
              & columnConstraintRowIndices,columnConstraintColumnIndices,ERR,ERROR,*999)
            SELECT CASE(equationsStorageType)
            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  !Loop over the columns of the equations matrix
                  DO equationsColumnNumber=1,EQUATIONS_MATRIX%TOTAL_NUMBER_OF_COLUMNS
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !equationsColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the equations matrix
                    DO equationsColumnNumber=1,EQUATIONS_MATRIX%TOTAL_NUMBER_OF_COLUMNS
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !equationsColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                equationsColumnNumber=equationsRowNumber
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                    CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                  ELSE
                    !The column dof is constrained
                    columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                    DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                      & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                      columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    END DO !constraintColumnIdx
                  END IF
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !columnConstraintColumnIdx
                    END IF
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              NULLIFY(equationsRowIndices)
              NULLIFY(equationsColumnIndices)
              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX, &
                & equationsRowIndices,equationsColumnIndices,ERR,ERROR,*999)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  !Loop over the columns of the equations matrix
                  DO equationsColumnIdx=equationsRowIndices(equationsRowNumber),equationsRowIndices(equationsRowNumber+1)-1
                    equationsColumnNumber=equationsColumnIndices(equationsColumnIdx)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !equationsColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the equations matrix
                    DO equationsColumnIdx=equationsRowIndices(equationsRowNumber),equationsRowIndices(equationsRowNumber+1)-1
                      equationsColumnNumber=equationsColumnIndices(equationsColumnIdx)
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !equationsColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(equationsStorageType,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          END DO !equationsMatrixIdx
        END IF
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          DO equationsMatrixIdx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
            EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(equationsMatrixIdx)%ptr
            DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,equationsStorageType,ERR,ERROR,*999)
            columnSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)% &
              & equationsToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
              & linearMatrixToSolverVariableMap(equationsMatrixIdx)%ptr
            IF(.NOT.ASSOCIATED(columnSolverVariableMap)) &
              & CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
            columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
            NULLIFY(columnConstraintRowIndices)
            NULLIFY(columnConstraintColumnIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix, &
              & columnConstraintRowIndices,columnConstraintColumnIndices,ERR,ERROR,*999)
            SELECT CASE(equationsStorageType)
            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  !Loop over the columns of the equations matrix
                  DO equationsColumnNumber=1,EQUATIONS_MATRIX%TOTAL_NUMBER_OF_COLUMNS
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !equationsColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the equations matrix
                    DO equationsColumnNumber=1,EQUATIONS_MATRIX%TOTAL_NUMBER_OF_COLUMNS
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !equationsColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                equationsColumnNumber=equationsRowNumber
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                    CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                  ELSE
                    !The column dof is constrained
                    columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                    DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                      & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                      columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    END DO !constraintColumnIdx
                  END IF
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !columnConstraintColumnIdx
                    END IF
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              NULLIFY(equationsRowIndices)
              NULLIFY(equationsColumnIndices)
              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX, &
                & equationsRowIndices,equationsColumnIndices,ERR,ERROR,*999)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  !Loop over the columns of the equations matrix
                  DO equationsColumnIdx=equationsRowIndices(equationsRowNumber),equationsRowIndices(equationsRowNumber+1)-1
                    equationsColumnNumber=equationsColumnIndices(equationsColumnIdx)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !equationsColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the equations matrix
                    DO equationsColumnIdx=equationsRowIndices(equationsRowNumber),equationsRowIndices(equationsRowNumber+1)-1
                      equationsColumnNumber=equationsColumnIndices(equationsColumnIdx)
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !equationsColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(equationsStorageType,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          END DO !equationsMatrixIdx
        END IF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          DO equationsMatrixIdx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(equationsMatrixIdx)%ptr
            DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,equationsStorageType,ERR,ERROR,*999)
            columnSolverVariableMap=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)% &
              & equationsToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
              & linearMatrixToSolverVariableMap(equationsMatrixIdx)%ptr
            IF(.NOT.ASSOCIATED(columnSolverVariableMap)) &
              & CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
            columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
            NULLIFY(columnConstraintRowIndices)
            NULLIFY(columnConstraintColumnIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix, &
              & columnConstraintRowIndices,columnConstraintColumnIndices,ERR,ERROR,*999)
            SELECT CASE(equationsStorageType)
            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  !Loop over the columns of the equations matrix
                  DO equationsColumnNumber=1,JACOBIAN_MATRIX%TOTAL_NUMBER_OF_COLUMNS
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !equationsColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the equations matrix
                    DO equationsColumnNumber=1,JACOBIAN_MATRIX%TOTAL_NUMBER_OF_COLUMNS
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !equationsColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                equationsColumnNumber=equationsRowNumber
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                    CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                  ELSE
                    !The column dof is constrained
                    columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                    DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                      & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                      columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    END DO !constraintColumnIdx
                  END IF
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !columnConstraintColumnIdx
                    END IF
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              NULLIFY(equationsRowIndices)
              NULLIFY(equationsColumnIndices)
              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX, &
                & equationsRowIndices,equationsColumnIndices,ERR,ERROR,*999)
              !Loop over the rows of the equations matrix
              DO equationsRowNumber=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                IF(rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  !Loop over the columns of the equations matrix
                  DO equationsColumnIdx=equationsRowIndices(equationsRowNumber),equationsRowIndices(equationsRowNumber+1)-1
                    equationsColumnNumber=equationsColumnIndices(equationsColumnIdx)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !equationsColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(equationsRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the equations matrix
                    DO equationsColumnIdx=equationsRowIndices(equationsRowNumber),equationsRowIndices(equationsRowNumber+1)-1
                      equationsColumnNumber=equationsColumnIndices(equationsColumnIdx)
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(equationsColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !equationsColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(equationsStorageType,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          END DO !equationsMatrixIdx
        END IF
      END DO !equationsSetIdx
      
      DO interfaceConditionIdx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
        INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%INTERFACE_CONDITION
        SELECT CASE(INTERFACE_CONDITION%METHOD)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          INTERFACE_MATRICES=>INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_MATRICES 
          IF(.NOT.ASSOCIATED(INTERFACE_MATRICES)) CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
          DO interfaceMatrixIdx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
            INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(interfaceMatrixIdx)%ptr
            DISTRIBUTED_MATRIX=>INTERFACE_MATRIX%MATRIX
            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,interfaceStorageType,ERR,ERROR,*999)
            rowSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)% &
              & interfaceToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
              & interfaceMatrixToSolverVariableMap(INTERFACE_MATRIX%MATRIX_NUMBER)%ptr
            IF(.NOT.ASSOCIATED(rowSolverVariableMap)) CALL FlagError("Row solver variable map is not associated.",ERR,ERROR,*999)
            !The Lagrange variable
            columnSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%solverVariableMap
            IF(.NOT.ASSOCIATED(columnSolverVariableMap)) CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
            rowConstraintMatrix=>rowSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
            NULLIFY(rowConstraintRowIndices)
            NULLIFY(rowConstraintColumnIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(rowConstraintMatrix,rowConstraintRowIndices, &
              & rowConstraintColumnIndices,ERR,ERROR,*999)
            columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
            NULLIFY(columnConstraintRowIndices)
            NULLIFY(columnConstraintColumnIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix,columnConstraintRowIndices, &
            & columnConstraintColumnIndices,ERR,ERROR,*999)
            interfaceMatrixNumberOfRows=DISTRIBUTED_MATRIX%ROW_DOMAIN_MAPPING%NUMBER_OF_LOCAL
            interfaceMatrixNumberOfColumns=DISTRIBUTED_MATRIX%COLUMN_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
            SELECT CASE(interfaceStorageType)
            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
              !Loop over the rows of the interface matrix
              DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
                IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  !Loop over the columns of the interface matrix
                  DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !interfaceColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the interface matrix
                    DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !interfaceColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !interfaceRowNumber
            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
              !Loop over the rows of the interface matrix
              DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
                interfaceColumnNumber=interfaceRowNumber
                IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                  ELSE
                    !The column dof is constrained
                    columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                      & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                      columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    END DO !constraintColumnIdx
                  END IF
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !columnConstraintColumnIdx
                    END IF
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !interfaceRowNumber
            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              NULLIFY(interfaceColumnIndices)
              NULLIFY(interfaceRowIndices)
              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX, &
                & interfaceRowIndices,interfaceColumnIndices,ERR,ERROR,*999)
              !Loop over the rows of the interface matrix
              DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
                IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  !Loop over the columns of the interface matrix
                  DO interfaceColumnIdx=interfaceRowIndices(interfaceRowNumber),interfaceRowIndices(interfaceRowNumber+1)-1
                    interfaceColumnNumber=interfaceColumnIndices(interfaceColumnIdx)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !interfaceColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the interface matrix
                    DO interfaceColumnIdx=interfaceRowIndices(interfaceRowNumber),interfaceRowIndices(interfaceRowNumber+1)-1
                      interfaceColumnNumber=interfaceColumnIndices(interfaceColumnIdx)
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !interfaceColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !interfaceRowNumber
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(interfaceStorageType,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
            !The Lagrange variable
            rowSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)%solverVariableMap
            IF(.NOT.ASSOCIATED(rowSolverVariableMap)) &
              & CALL FlagError("Row solver variable map is not associated.",ERR,ERROR,*999)
            columnSolverVariableMap=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interfaceConditionIdx)% &
              & interfaceToSolverVariablesMaps(SOLVER_MATRIX%MATRIX_NUMBER)% &
              & interfaceMatrixToSolverVariableMap(INTERFACE_MATRIX%MATRIX_NUMBER)%ptr
            IF(.NOT.ASSOCIATED(columnSolverVariableMap)) &
              & CALL FlagError("Column solver variable map is not associated.",ERR,ERROR,*999)
            rowConstraintMatrix=>rowSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
            NULLIFY(rowConstraintRowIndices)
            NULLIFY(rowConstraintColumnIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(rowConstraintMatrix,rowConstraintRowIndices, &
              & rowConstraintColumnIndices,ERR,ERROR,*999)
            columnConstraintMatrix=>columnSolverVariableMap%boundaryConditionsVariable%constraints%constraintMatrix
            NULLIFY(columnConstraintRowIndices)
            NULLIFY(columnConstraintColumnIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(columnConstraintMatrix,columnConstraintRowIndices, &
            & columnConstraintColumnIndices,ERR,ERROR,*999)
            DISTRIBUTED_MATRIX=>INTERFACE_MATRIX%MATRIX_TRANSPOSE
            NULLIFY(interfaceColumnIndices)
            NULLIFY(interfaceRowIndices)
            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,interfaceStorageType,ERR,ERROR,*999)
            interfaceMatrixNumberOfRows=DISTRIBUTED_MATRIX%ROW_DOMAIN_MAPPING%NUMBER_OF_LOCAL
            interfaceMatrixNumberOfColumns=DISTRIBUTED_MATRIX%COLUMN_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
            SELECT CASE(interfaceStorageType)
            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
              !Loop over the rows of the interface matrix
              DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
                IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  !Loop over the columns of the interface matrix
                  DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !interfaceColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the interface matrix
                    DO interfaceColumnNumber=1,interfaceMatrixNumberOfColumns
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !interfaceColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !interfaceRowNumber
            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
              !Loop over the rows of the interface matrix
              DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
                interfaceColumnNumber=interfaceRowNumber
                IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                    solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                  ELSE
                    !The column dof is constrained
                    columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                    DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                      & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                      columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    END DO !constraintColumnIdx
                  END IF
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !columnConstraintColumnIdx
                    END IF
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !interfaceRowNumber
            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              NULLIFY(interfaceRowIndices)
              NULLIFY(interfaceColumnIndices)
              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX, &
                & interfaceRowIndices,interfaceColumnIndices,ERR,ERROR,*999)
              !Loop over the rows of the interface matrix
              DO interfaceRowNumber=1,interfaceMatrixNumberOfRows
                IF(rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)>0) THEN
                  solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  !Loop over the columns of the interface matrix
                  DO interfaceColumnIdx=interfaceRowIndices(interfaceRowNumber),interfaceRowIndices(interfaceRowNumber+1)-1
                    interfaceColumnNumber=interfaceColumnIndices(interfaceColumnIdx)
                    IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                      solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                    ELSE
                      !The column dof is constrained
                      columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                      DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                        & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                        columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      END DO !constraintColumnIdx
                    END IF
                  ENDDO !interfaceColumnIdx
                ELSE
                  !The row dof is constrained
                  rowConstraintRowNumber=-rowSolverVariableMap%variableDofToSolverDofMap(interfaceRowNumber)
                  DO rowConstraintColumnIdx=rowConstraintRowIndices(rowConstraintRowNumber), &
                    & rowConstraintRowIndices(rowConstraintRowNumber+1)-1
                    rowConstraintColumnNumber=rowConstraintColumnIndices(rowConstraintColumnIdx)
                    solverRowNumber=rowSolverVariableMap%variableDofToSolverDofMap(rowConstraintColumnNumber)
                    !Loop over the columns of the interface matrix
                    DO interfaceColumnIdx=interfaceRowIndices(interfaceRowNumber),interfaceRowIndices(interfaceRowNumber+1)-1
                      interfaceColumnNumber=interfaceColumnIndices(interfaceColumnIdx)
                      IF(columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)>0) THEN
                        solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                      ELSE
                        !The column dof is constrained
                        columnConstraintRowNumber=-columnSolverVariableMap%variableDofToSolverDofMap(interfaceColumnNumber)
                        DO columnConstraintColumnIdx=columnConstraintRowIndices(columnConstraintRowNumber), &
                          & columnConstraintRowIndices(columnConstraintRowNumber+1)-1
                          columnConstraintColumnNumber=columnConstraintColumnIndices(columnConstraintColumnIdx)
                          solverColumnNumber=columnSolverVariableMap%variableDofToSolverDofMap(columnConstraintColumnNumber)
                          CALL List_ItemAdd(columnIndicesLists(solverRowNumber)%PTR,solverColumnNumber,ERR,ERROR,*999)
                        END DO !columnConstraintColumnIdx
                      END IF
                    END DO !interfaceColumnIdx
                  END DO !rowConstraintColumnIdx
                END IF
              END DO !interfaceRowNumber
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(interfaceStorageType,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          END IF
          END DO !interfaceMatrixIdx
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The interface condition method of "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      END DO !interfaceConditionIdx

      !Loop over the rows to calculate the number of non-zeros and setup the row indicces
      rowIndices(1)=1
      DO solverRowNumber=1,SOLVER_MAPPING%NUMBER_OF_DOFS
        CALL List_RemoveDuplicates(columnIndicesLists(solverRowNumber)%PTR,ERR,ERROR,*999)
        CALL List_NumberOfItemsGet(columnIndicesLists(solverRowNumber)%PTR,numberOfColumns,ERR,ERROR,*999)
        rowIndices(solverRowNumber+1)=rowIndices(solverRowNumber)+numberOfColumns
      END DO !solverRowNumber
      !Allocate and setup the column locations
      numberOfNonZeros=rowIndices(SOLVER_MAPPING%NUMBER_OF_DOFS+1)-1
      ALLOCATE(columnIndices(numberOfNonZeros),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate column indices.",ERR,ERROR,*999)
      DO solverRowNumber=1,SOLVER_MAPPING%NUMBER_OF_DOFS
        CALL List_DetachAndDestroy(columnIndicesLists(solverRowNumber)%PTR,numberOfColumns,columns,ERR,ERROR,*999)
        columnIndices(rowIndices(solverRowNumber):rowIndices(solverRowNumber+1)-1)=columns
        DEALLOCATE(columns)
      END DO !solverRowNumber
    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)                        
    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
      CALL FlagError("Not implemented.",ERR,ERROR,*999)                      
    CASE DEFAULT
      LOCAL_ERROR="The matrix storage type of "// &
        & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix structure:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix number : ",SOLVER_MATRIX%MATRIX_NUMBER, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",SOLVER_MATRICES%NUMBER_OF_ROWS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",SOLVER_MATRIX%NUMBER_OF_GLOBAL_COLUMNS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",numberOfNonZeros,ERR,ERROR,*999)
      IF(SOLVER_MATRICES%NUMBER_OF_ROWS*SOLVER_MATRIX%NUMBER_OF_GLOBAL_COLUMNS/=0) THEN
        SPARSITY=REAL(numberOfNonZeros,DP)/REAL(SOLVER_MATRICES%NUMBER_OF_ROWS*SOLVER_MATRIX%NUMBER_OF_GLOBAL_COLUMNS,DP)*100.0_DP
        CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY,"F6.2", ERR,ERROR,*999)
      ENDIF
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MATRICES%NUMBER_OF_ROWS+1,8,8,rowIndices, &
          & '("  Row indices    :",8(X,I13))','(18X,8(X,I13))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', ERR,ERROR,*999)
      ENDIF
    ENDIF

    EXITS("SOLVER_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ALLOCATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ALLOCATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO solverRowNumber=1,SOLVER_MAPPING%NUMBER_OF_DOFS
        IF(ASSOCIATED(columnIndicesLists(solverRowNumber)%PTR)) &
          & CALL List_Destroy(columnIndicesLists(solverRowNumber)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !solverRowNumber
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("SOLVER_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_STRUCTURE_CALCULATE
        
  !
  !================================================================================================================================
  !

  !>Finalises a solver matrix and deallocates all memory
  SUBROUTINE SOLVER_MATRIX_FINALISE(SOLVER_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("SOLVER_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
      IF(ASSOCIATED(SOLVER_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLVER_MATRIX%SOLVER_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(SOLVER_MATRIX)
    ENDIF
    
    EXITS("SOLVER_MATRIX_FINALISE")
    RETURN
999 ERRORSEXITS("SOLVER_MATRIX_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Forms a solver matrix by initialising the structure of the matrix to zero.
  SUBROUTINE SOLVER_MATRIX_FORM(SOLVER_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("SOLVER_MATRIX_FORM",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
      CALL DISTRIBUTED_MATRIX_FORM(SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MATRIX_FORM")
    RETURN
999 ERRORSEXITS("SOLVER_MATRIX_FORM",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_FORM
        
    !
  !================================================================================================================================
  !

  !>Initialises a solver matrix
  SUBROUTINE SOLVER_MATRIX_INITIALISE(SOLVER_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices to initialise
    INTEGER(INTG), INTENT(IN) :: MATRIX_NUMBER !<The matrix number in the solver matrices to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    ENTERS("SOLVER_MATRIX_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=SOLVER_MATRICES%NUMBER_OF_MATRICES) THEN
        SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          IF(ASSOCIATED(SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
            CALL FlagError("Solver matrix is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate solver matrix.",ERR,ERROR,*999)
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
            SOLVER_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
            SOLVER_MATRIX%SOLVER_MATRICES=>SOLVER_MATRICES
            SOLVER_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
            SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
            SOLVER_MATRIX%NUMBER_OF_COLUMNS=SOLVER_MAPPING%NUMBER_OF_DOFS
            SOLVER_MATRIX%TOTAL_NUMBER_OF_COLUMNS=SOLVER_MAPPING%TOTAL_NUMBER_OF_DOFS
            SOLVER_MATRIX%NUMBER_OF_GLOBAL_COLUMNS=SOLVER_MAPPING%NUMBER_OF_GLOBAL_DOFS
            NULLIFY(SOLVER_MATRIX%SOLVER_VECTOR)
            NULLIFY(SOLVER_MATRIX%MATRIX)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Solver matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("SOLVER_MATRIX_INITIALISE")
    RETURN
999 CALL SOLVER_MATRIX_FINALISE(SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("SOLVER_MATRIX_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_INITIALISE
        
  !
  !================================================================================================================================
  !

END MODULE SOLVER_MATRICES_ROUTINES
