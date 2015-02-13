!> \file
!> \author Chris Bradley
!> \brief This module handles all constraint condition routines.
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
!> Contributor(s): Sander Arens
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

!> This module handles all constraint condition routines.
MODULE CONSTRAINT_CONDITION_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE FIELD_ROUTINES
  USE CONSTRAINT_ROUTINES
  USE CONSTRAINT_CONDITION_CONSTANTS
  USE CONSTRAINT_MATRICES_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE MPI
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC CONSTRAINT_CONDITION_ASSEMBLE
  
  PUBLIC CONSTRAINT_CONDITION_NONLINEAR_RHS_UPDATE
  
  PUBLIC CONSTRAINT_CONDITION_CREATE_START,CONSTRAINT_CONDITION_CREATE_FINISH

  PUBLIC CONSTRAINT_CONDITION_DESTROY

  PUBLIC CONSTRAINT_CONDITIONS_FINALISE,CONSTRAINT_CONDITIONS_INITIALISE

  PUBLIC CONSTRAINT_CONDITION_CONSTRAINT_CREATE_FINISH,CONSTRAINT_CONDITION_CONSTRAINT_CREATE_START

  PUBLIC CONSTRAINT_CONDITION_CONSTRAINT_DESTROY
  
  PUBLIC CONSTRAINT_CONDITION_DEPENDENT_CREATE_START,CONSTRAINT_CONDITION_DEPENDENT_CREATE_FINISH

  PUBLIC CONSTRAINT_CONDITION_DEPENDENT_DESTROY
  
  PUBLIC CONSTRAINT_CONDITION_JACOBIAN_EVALUATE,CONSTRAINT_CONDITION_RESIDUAL_EVALUATE
  
  PUBLIC CONSTRAINT_CONDITION_SOLUTION_METHOD_GET,CONSTRAINT_CONDITION_SOLUTION_METHOD_SET
  
  PUBLIC CONSTRAINT_CONDITION_USER_NUMBER_FIND

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assembles the equations for an constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_CONDITION_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(CONSTRAINT_EQUATIONS%CONSTRAINT_FINISHED) THEN
          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
            SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
            CASE(CONSTRAINT_CONDITION_STATIC)
              SELECT CASE(CONSTRAINT_EQUATIONS%LINEARITY)
              CASE(CONSTRAINT_CONDITION_LINEAR)
                SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
                CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                  CALL CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The constraint condition solution method of "// &
                    & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE(CONSTRAINT_CONDITION_NONLINEAR)
                SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
                CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                  CALL CONSTRAINT_CONDITION_ASSEMBLE_STATIC_NONLINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The constraint condition solution method of "// &
                    & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE(CONSTRAINT_CONDITION_NONLINEAR_BCS)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The equations linearity of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC,CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
              SELECT CASE(CONSTRAINT_EQUATIONS%LINEARITY)
              CASE(CONSTRAINT_CONDITION_LINEAR)
                SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
                CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                  CALL CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The constraint condition solution method of "// &
                    & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE(CONSTRAINT_CONDITION_NONLINEAR)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_NONLINEAR_BCS)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition linearity of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(CONSTRAINT_TIME_STEPPING)
              CALL FLAG_ERROR("Time stepping equations are not assembled.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The constraint condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))// &
              & "is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Constraint equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint condition equations is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_ASSEMBLE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix and rhs for a dynamic linear constraint condition using the finite element method.
  SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    
    CALL ENTERS("CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      LAGRANGE_FIELD=>CONSTRAINT_CONDITION%LAGRANGE%LAGRANGE_FIELD
      IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Initialise the matrices and rhs vector
            CALL CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,CONSTRAINT_MATRICES_LINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>LAGRANGE_FIELD%DECOMPOSITION%DOMAIN(LAGRANGE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for constraint conditionup and initialisation = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for constraint conditionup and initialisation = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal constraint equations assembly = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal constraint equations assembly = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
             ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost constraint equations assembly = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost constraint equations assembly = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for constraint equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for constraint equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_MATRIX_OUTPUT) THEN
              CALL CONSTRAINT_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for constraint equations assembly = ", & 
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for constraint equations assembly = ", & 
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Lagrange field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix and rhs for a linear static constraint condition using the finite element method.
  SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to assemble the constraint equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    
!#ifdef TAUPROF
!    CHARACTER(28) :: CVAR
!    INTEGER :: PHASE(2)= (/ 0, 0 /)
!    SAVE PHASE
!#endif

    CALL ENTERS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      LAGRANGE_FIELD=>CONSTRAINT_CONDITION%LAGRANGE%LAGRANGE_FIELD
      IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Initialise the matrices and rhs vector
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_VALUES_INITIALISE()")
#endif
            CALL CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,CONSTRAINT_MATRICES_LINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_MATRICES_VALUES_INITIALISE()")
#endif
            !Assemble the elements
            !Allocate the element matrices 
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_ELEMENT_INITIALISE()")
#endif
            CALL CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>LAGRANGE_FIELD%DECOMPOSITION%DOMAIN(LAGRANGE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_MATRICES_ELEMENT_INITIALISE()")
#endif
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for constraint conditionup and initialisation = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for constraint conditionup and initialisation = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements

#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("Internal Elements Loop")
#endif
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
!#ifdef TAUPROF
!              WRITE (CVAR,'(a23,i3)') 'Internal Elements Loop ',element_idx
!              CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
!              CALL TAU_PHASE_START(PHASE)
!#endif
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
!#ifdef TAUPROF
!              CALL TAU_PHASE_STOP(PHASE)
!#endif
            ENDDO !element_idx
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("Internal Elements Loop")
#endif

            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal constraint equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal constraint equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
            ENDIF
            !Loop over the boundary and ghost elements
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("Boundary and Ghost Elements Loop")
#endif
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("Boundary and Ghost Elements Loop")
#endif
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost constraint equations assembly = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost constraint equations assembly = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for constraint equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for constraint equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("CONSTRAINT_MATRICES_ELEMENT_FINALISE()")
#endif
            CALL CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_MATRICES_ELEMENT_FINALISE()")
#endif
            !Output constraint equations matrices and vector if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_MATRIX_OUTPUT) THEN
              CALL CONSTRAINT_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for constraint equations assembly = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for constraint equations assembly = ", & 
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint constraint equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint constraint equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Lagrange field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM

  !
  !================================================================================================================================
  !
  
  !>Assembles the constraint equations stiffness matrix, residuals and rhs for a nonlinear static constraint condition using the finite element method.
  SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE_STATIC_NONLINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to assemble the constraint equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    
    CALL ENTERS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_NONLINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      LAGRANGE_FIELD=>CONSTRAINT_CONDITION%LAGRANGE%LAGRANGE_FIELD
      IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
             !Initialise the matrices and rhs vector
            CALL CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES, &
              & CONSTRAINT_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>LAGRANGE_FIELD%DECOMPOSITION%DOMAIN(LAGRANGE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for constraint conditionup and initialisation = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for constraint conditionup and initialisation = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal constraint equations assembly = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal constraint equations assembly = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost constraint equations assembly = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost constraint equations assembly = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for constraint equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for constraint equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            !Output constraint equations matrices and RHS vector if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_MATRIX_OUTPUT) THEN
              CALL CONSTRAINT_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for constraint equations assembly = ", &
                & USER_ELAPSED,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for constraint equations assembly = ", &
                & SYSTEM_ELAPSED,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint constraint equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint constraint equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Lagrange field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_NONLINEAR_FEM")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_NONLINEAR_FEM",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_ASSEMBLE_STATIC_NONLINEAR_FEM")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_ASSEMBLE_STATIC_NONLINEAR_FEM

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equation set on a region. \see OPENCMISS::CMISSConstraintConditionCreateStart
  SUBROUTINE CONSTRAINT_CONDITION_CREATE_FINISH(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_SETUP_TYPE) :: CONSTRAINT_CONDITION_SETUP_INFO
    
    CALL ENTERS("CONSTRAINT_CONDITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      IF(CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Constraint condition has already been finished.",ERR,ERROR,*999)
      ELSE
        !Test various inputs have been set up.
        SELECT CASE(CONSTRAINT_CONDITION%METHOD)
        CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
          CONSTRAINT_DEPENDENT=>CONSTRAINT_CONDITION%DEPENDENT
          IF(ASSOCIATED(CONSTRAINT_DEPENDENT)) THEN
          !Do nothing?
        ELSE
          CALL FLAG_ERROR("Constraint condition dependen is not associated.",ERR,ERROR,*999)
        ENDIF
      CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The constraint condition method of "//TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))// &
          & "is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
        !Finish the constraint condition creation
        CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_CONDITION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_CONDITION_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating an constraint condition defined by USER_NUMBER in the region identified by REGION. \see OPENCMISS::CMISSConstraintConditionCreateStart

  SUBROUTINE CONSTRAINT_CONDITION_CREATE_START(USER_NUMBER,REGION,GEOM_FIBRE_FIELD,CONSTRAINT_CONDITION_CLASS,CONSTRAINT_CONDITION_TYPE_,&
    & CONSTRAINT_CONDITION_SUBTYPE,CONSTRAINT_CONDITION_FIELD_USER_NUMBER,CONSTRAINT_CONDITION_FIELD_FIELD,CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the constraint condition
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the constraint condition on
    TYPE(FIELD_TYPE), POINTER :: GEOM_FIBRE_FIELD !<A pointer to the either the geometry or, if appropriate, the fibre field for the equation set
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_CONDITION_FIELD_USER_NUMBER !<The user number of the constraint condition field
    TYPE(FIELD_TYPE), POINTER :: CONSTRAINT_CONDITION_FIELD_FIELD !<On return, a pointer to the constraint condition field
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<On return, a pointer to the constraint condition
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_CONDITION_CLASS !<The constraint condition class to set
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_CONDITION_TYPE_ !<The constraint condition type to set
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_CONDITION_SUBTYPE !<The constraint condition subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,equations_set_idx
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: NEW_CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_CONDITION_PTR_TYPE), POINTER :: NEW_CONSTRAINT_CONDITIONS(:)
    TYPE(CONSTRAINT_CONDITION_SETUP_TYPE) :: CONSTRAINT_CONDITION_SETUP_INFO
    TYPE(REGION_TYPE), POINTER :: GEOM_FIBRE_FIELD_REGION,CONSTRAINT_CONDITION_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    TYPE(CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_TYPE), POINTER :: CONSTRAINT_CONSTRAINT_CONDITION_FIELD
    TYPE(FIELD_TYPE), POINTER :: FIELD

    NULLIFY(NEW_CONSTRAINT_CONDITION)
    NULLIFY(NEW_CONSTRAINT_CONDITIONS)
    NULLIFY(CONSTRAINT_CONSTRAINT_CONDITION_FIELD)

    CALL ENTERS("CONSTRAINT_CONDITION_CREATE_START",ERR,ERROR,*997)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS)) THEN
        CALL CONSTRAINT_CONDITION_USER_NUMBER_FIND(USER_NUMBER,REGION,NEW_CONSTRAINT_CONDITION,ERR,ERROR,*997)
        IF(ASSOCIATED(NEW_CONSTRAINT_CONDITION)) THEN
          LOCAL_ERROR="Constraint condition user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
        ELSE
          NULLIFY(NEW_CONSTRAINT_CONDITION)
          IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
            IF(GEOM_FIBRE_FIELD%FIELD_FINISHED) THEN
              IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR.GEOM_FIBRE_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
                GEOM_FIBRE_FIELD_REGION=>GEOM_FIBRE_FIELD%REGION
                IF(ASSOCIATED(GEOM_FIBRE_FIELD_REGION)) THEN
                  IF(GEOM_FIBRE_FIELD_REGION%USER_NUMBER==REGION%USER_NUMBER) THEN
                      IF(ASSOCIATED(CONSTRAINT_CONDITION_FIELD_FIELD)) THEN
                        !Check the constraint condition field has been finished
                        IF(CONSTRAINT_CONDITION_FIELD_FIELD%FIELD_FINISHED.eqv..TRUE.) THEN
                          !Check the user numbers match
                          IF(CONSTRAINT_CONDITION_FIELD_USER_NUMBER/=CONSTRAINT_CONDITION_FIELD_FIELD%USER_NUMBER) THEN
                            LOCAL_ERROR="The specified constraint condition field user number of "// &
                              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                              & " does not match the user number of the specified constraint condition field of "// &
                              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION_FIELD_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                          CONSTRAINT_CONDITION_FIELD_REGION=>CONSTRAINT_CONDITION_FIELD_FIELD%REGION
                          IF(ASSOCIATED(CONSTRAINT_CONDITION_FIELD_REGION)) THEN                
                            !Check the field is defined on the same region as the constraint condition
                            IF(CONSTRAINT_CONDITION_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                              LOCAL_ERROR="Invalid region setup. The specified constraint condition field was created on region no. "// &
                                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                                & " and the specified constraint condition has been created on region number "// &
                                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            !Check the specified constraint condition field has the same decomposition as the geometric field
                            IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
                              IF(.NOT.ASSOCIATED(GEOM_FIBRE_FIELD%DECOMPOSITION,CONSTRAINT_CONDITION_FIELD_FIELD%DECOMPOSITION)) THEN
                                CALL FLAG_ERROR("The specified constraint condition field does not have the same decomposition "// &
                                  & "as the geometric field for the specified constraint condition.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("The geom. field is not associated for the specified constraint condition.",ERR,ERROR,*999)
                            ENDIF
                              
                          ELSE
                            CALL FLAG_ERROR("The specified constraint condition field region is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("The specified constraint condition field has not been finished.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        !Check the user number has not already been used for a field in this region.
                        NULLIFY(FIELD)
                        CALL FIELD_USER_NUMBER_FIND(CONSTRAINT_CONDITION_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
                        IF(ASSOCIATED(FIELD)) THEN
                          LOCAL_ERROR="The specified constraint condition field user number of "// &
                            & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                            & "has already been used to create a field on region number "// &
                            & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ENDIF
                      !Initialise the constraint condition materials
!                       CALL CONSTRAINT_CONDITION_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*999)
!                        WRITE(*,'(A)') "constraint condition initialise called"
!                       IF(.NOT.ASSOCIATED(CONSTRAINT_CONDITION_FIELD_FIELD)) THEN
!                         CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_AUTO_CREATED=.TRUE.
!                       ENDIF
!--- tob 1            
                      !Initalise constraint condition
                      CALL CONSTRAINT_CONDITION_INITIALISE(NEW_CONSTRAINT_CONDITION,ERR,ERROR,*999)
                      !Set default constraint condition values
                      NEW_CONSTRAINT_CONDITION%USER_NUMBER=USER_NUMBER
                      NEW_CONSTRAINT_CONDITION%GLOBAL_NUMBER=REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS+1
                      NEW_CONSTRAINT_CONDITION%CONSTRAINT_CONDITIONS=>REGION%CONSTRAINT_CONDITIONS
                      NEW_CONSTRAINT_CONDITION%REGION=>REGION
                      !Set the constraint condition class, type and subtype
                      CALL CONSTRAINT_CONDITION_SPECIFICATION_SET(NEW_CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_CLASS,CONSTRAINT_CONDITION_TYPE_, &
                        & CONSTRAINT_CONDITION_SUBTYPE,ERR,ERROR,*999)
                      NEW_CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FINISHED=.FALSE.
                      !Initialise the setup
                      CALL CONSTRAINT_CONDITION_SETUP_INITIALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
                      CONSTRAINT_CONDITION_SETUP_INFO%SETUP_TYPE=CONSTRAINT_CONDITION_SETUP_INITIAL_TYPE
                      CONSTRAINT_CONDITION_SETUP_INFO%ACTION_TYPE=CONSTRAINT_CONDITION_SETUP_START_ACTION
                      !Here, we get a pointer to the equations_set_field; default is null
                      CONSTRAINT_CONDITION_SETUP_INFO%FIELD_USER_NUMBER=CONSTRAINT_CONDITION_FIELD_USER_NUMBER
                      CONSTRAINT_CONDITION_SETUP_INFO%FIELD=>CONSTRAINT_CONDITION_FIELD_FIELD
                      !Start constraint condition specific setup
                      CALL CONSTRAINT_CONDITION_SETUP(NEW_CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
                      CALL CONSTRAINT_CONDITION_SETUP_FINALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
                      !Set up the constraint condition geometric fields
                      CALL CONSTRAINT_CONDITION_GEOMETRY_INITIALISE(NEW_CONSTRAINT_CONDITION,ERR,ERROR,*999)
                      IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
                        NEW_CONSTRAINT_CONDITION%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD
                        NULLIFY(NEW_CONSTRAINT_CONDITION%GEOMETRY%FIBRE_FIELD)
                      ELSE
                        NEW_CONSTRAINT_CONDITION%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD%GEOMETRIC_FIELD
                        NEW_CONSTRAINT_CONDITION%GEOMETRY%FIBRE_FIELD=>GEOM_FIBRE_FIELD
                      ENDIF
                      CONSTRAINT_CONDITION_SETUP_INFO%SETUP_TYPE=CONSTRAINT_CONDITION_SETUP_GEOMETRY_TYPE
                      CONSTRAINT_CONDITION_SETUP_INFO%ACTION_TYPE=CONSTRAINT_CONDITION_SETUP_START_ACTION
                      CONSTRAINT_CONDITION_SETUP_INFO%FIELD_USER_NUMBER=GEOM_FIBRE_FIELD%USER_NUMBER
                      CONSTRAINT_CONDITION_SETUP_INFO%FIELD=>GEOM_FIBRE_FIELD
                      !Set up constraint condition specific geometry
                      CALL CONSTRAINT_CONDITION_SETUP(NEW_CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
                      !Finalise the setup
                      CALL CONSTRAINT_CONDITION_SETUP_FINALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
                      !Add new constraint condition into list of constraint condition in the region
                      ALLOCATE(NEW_CONSTRAINT_CONDITIONS(REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS+1),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new constraint conditions.",ERR,ERROR,*999)
                      DO equations_set_idx=1,REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS
                        NEW_CONSTRAINT_CONDITIONS(equations_set_idx)%PTR=>REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR
                      ENDDO !equations_set_idx
                      NEW_CONSTRAINT_CONDITIONS(REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS+1)%PTR=>NEW_CONSTRAINT_CONDITION
                      IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)) DEALLOCATE(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)
                      REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS=>NEW_CONSTRAINT_CONDITIONS
                      REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS=REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS+1
                      CONSTRAINT_CONDITION=>NEW_CONSTRAINT_CONDITION
                      CONSTRAINT_CONSTRAINT_CONDITION_FIELD=>CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD
                      !\todo check pointer setup
                      IF(CONSTRAINT_CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_AUTO_CREATED) THEN
                        CONSTRAINT_CONDITION_FIELD_FIELD=>CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_FIELD
                      ELSE
                        CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_FIELD=>CONSTRAINT_CONDITION_FIELD_FIELD
                      ENDIF
                  ELSE
                    LOCAL_ERROR="The geometric field region and the specified region do not match. "// &
                      & "The geometric field was created on region number "// &
                      & TRIM(NUMBER_TO_VSTRING(GEOM_FIBRE_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                      & " and the specified region number is "// &
                      & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The specified geometric fields region is not associated.",ERR,ERROR,*997)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified geometric field is not a geometric or fibre field.",ERR,ERROR,*997)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified geometric field is not finished.",ERR,ERROR,*997)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The specified geometric field is not associated.",ERR,ERROR,*997)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The constraint conditions on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*997)
    ENDIF
    
    CALL EXITS("CONSTRAINT_CONDITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_CONSTRAINT_CONDITION))CALL CONSTRAINT_CONDITION_FINALISE(NEW_CONSTRAINT_CONDITION,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_CONSTRAINT_CONDITIONS)) DEALLOCATE(NEW_CONSTRAINT_CONDITIONS)
997 CALL ERRORS("CONSTRAINT_CONDITION_CREATE_START",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_CREATE_START")
    RETURN 1   
  END SUBROUTINE CONSTRAINT_CONDITION_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroys an constraint condition identified by a user number on the give region and deallocates all memory. \see OPENCMISS::CMISSConstraintConditionDestroy
  SUBROUTINE CONSTRAINT_CONDITION_DESTROY_NUMBER(USER_NUMBER,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the constraint condition to destroy
    TYPE(REGION_TYPE), POINTER :: REGION !<The region of the constraint condition to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,equations_set_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_CONDITION_PTR_TYPE), POINTER :: NEW_CONSTRAINT_CONDITIONS(:)

    NULLIFY(NEW_CONSTRAINT_CONDITIONS)

    CALL ENTERS("CONSTRAINT_CONDITION_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS)) THEN
        
        !Find the constraint condition identified by the user number
        FOUND=.FALSE.
        equations_set_position=0
        DO WHILE(equations_set_position<REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS.AND..NOT.FOUND)
          equations_set_position=equations_set_position+1
          IF(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_position)%PTR%USER_NUMBER==USER_NUMBER)FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          CONSTRAINT_CONDITION=>REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_position)%PTR
          
          !Destroy all the constraint condition components
          CALL CONSTRAINT_CONDITION_FINALISE(CONSTRAINT_CONDITION,ERR,ERROR,*999)
          
          !Remove the constraint condition from the list of constraint condition
          IF(REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS>1) THEN
            ALLOCATE(NEW_CONSTRAINT_CONDITIONS(REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new constraint conditions.",ERR,ERROR,*999)
            DO equations_set_idx=1,REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS
              IF(equations_set_idx<equations_set_position) THEN
                NEW_CONSTRAINT_CONDITIONS(equations_set_idx)%PTR=>REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR
              ELSE IF(equations_set_idx>equations_set_position) THEN
                REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR%GLOBAL_NUMBER=REGION%CONSTRAINT_CONDITIONS% &
                  & CONSTRAINT_CONDITIONS(equations_set_idx)%PTR%GLOBAL_NUMBER-1
                NEW_CONSTRAINT_CONDITIONS(equations_set_idx-1)%PTR=>REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR
              ENDIF
            ENDDO !equations_set_idx
            IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)) DEALLOCATE(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)
            REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS=>NEW_CONSTRAINT_CONDITIONS
            REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS=REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS-1
          ELSE
            DEALLOCATE(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)
            REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Constraint condition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The constraint conditions on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("CONSTRAINT_CONDITION_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_CONSTRAINT_CONDITIONS)) DEALLOCATE(NEW_CONSTRAINT_CONDITIONS)
998 CALL ERRORS("CONSTRAINT_CONDITION_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_DESTROY_NUMBER")
    RETURN 1   
  END SUBROUTINE CONSTRAINT_CONDITION_DESTROY_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Destroys an constraint condition identified by a pointer and deallocates all memory. \see OPENCMISS::CMISSConstraintConditionDestroy
  SUBROUTINE CONSTRAINT_CONDITION_DESTROY(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,equations_set_position
    TYPE(CONSTRAINT_CONDITIONS_TYPE), POINTER :: CONSTRAINT_CONDITIONS
    TYPE(CONSTRAINT_CONDITION_PTR_TYPE), POINTER :: NEW_CONSTRAINT_CONDITIONS(:)

    NULLIFY(NEW_CONSTRAINT_CONDITIONS)

    CALL ENTERS("CONSTRAINT_CONDITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_CONDITIONS=>CONSTRAINT_CONDITION%CONSTRAINT_CONDITIONS
      IF(ASSOCIATED(CONSTRAINT_CONDITIONS)) THEN
        equations_set_position=CONSTRAINT_CONDITION%GLOBAL_NUMBER

        !Destroy all the constraint condition components
        CALL CONSTRAINT_CONDITION_FINALISE(CONSTRAINT_CONDITION,ERR,ERROR,*999)
        
        !Remove the constraint condition from the list of constraint condition
        IF(CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS>1) THEN
          ALLOCATE(NEW_CONSTRAINT_CONDITIONS(CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new constraint conditions.",ERR,ERROR,*999)
          DO equations_set_idx=1,CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS
            IF(equations_set_idx<equations_set_position) THEN
              NEW_CONSTRAINT_CONDITIONS(equations_set_idx)%PTR=>CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR
            ELSE IF(equations_set_idx>equations_set_position) THEN
              CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR%GLOBAL_NUMBER=CONSTRAINT_CONDITIONS% &
                & CONSTRAINT_CONDITIONS(equations_set_idx)%PTR%GLOBAL_NUMBER-1
              NEW_CONSTRAINT_CONDITIONS(equations_set_idx-1)%PTR=>CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR
            ENDIF
          ENDDO !equations_set_idx
          IF(ASSOCIATED(CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)) DEALLOCATE(CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)
          CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS=>NEW_CONSTRAINT_CONDITIONS
          CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS=CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS-1
        ELSE
          DEALLOCATE(CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)
          CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Constraint condition constraint condition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("CONSTRAINT_CONDITION_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_CONSTRAINT_CONDITIONS)) DEALLOCATE(NEW_CONSTRAINT_CONDITIONS)
998 CALL ERRORS("CONSTRAINT_CONDITION_DESTROY",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_DESTROY")
    RETURN 1   
  END SUBROUTINE CONSTRAINT_CONDITION_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalise the constraint condition and deallocate all memory.
  SUBROUTINE CONSTRAINT_CONDITION_FINALISE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CALL CONSTRAINT_CONDITION_GEOMETRY_FINALISE(CONSTRAINT_CONDITION%GEOMETRY,ERR,ERROR,*999)
      CALL CONSTRAINT_CONDITION_DEPENDENT_FINALISE(CONSTRAINT_CONDITION%DEPENDENT,ERR,ERROR,*999)
      CALL CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_FINALISE(CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD,ERR,ERROR,*999)
      IF(ASSOCIATED(CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS))CALL CONSTRAINT_DESTROY(CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS,ERR,ERROR,*999)
      DEALLOCATE(CONSTRAINT_CONDITION)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_FINALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a finite element constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(ELEMENT_VECTOR_TYPE), POINTER :: ELEMENT_VECTOR
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      SELECT CASE(CONSTRAINT_CONDITION%CLASS)
      CASE(CONSTRAINT_CONDITION_ELASTICITY_CLASS)
        CALL ELASTICITY_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_FLUID_MECHANICS_CLASS)
       CALL FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_FITTING_CLASS)
        CALL FITTING_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_BIOELECTRICS_CLASS)
        IF(CONSTRAINT_CONDITION%TYPE == CONSTRAINT_CONDITION_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
          CALL MONODOMAIN_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
        ELSE
          CALL BIOELECTRIC_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
        END IF
      CASE(CONSTRAINT_CONDITION_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_MULTI_PHYSICS_CLASS)
       CALL MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Constraint condition class "//TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_ELEMENT_MATRIX_OUTPUT) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element stiffness matrices:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
            DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
            IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Dynamic matrices:",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",DYNAMIC_MATRICES% &
                & NUMBER_OF_DYNAMIC_MATRICES,ERR,ERROR,*999)
              DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR% &
                  & UPDATE_MATRIX,ERR,ERROR,*999)
                IF(DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                  ELEMENT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                    & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                    & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                    & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                    & '(16X,8(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDDO !matrix_idx
            ENDIF
            LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
            IF(ASSOCIATED(LINEAR_MATRICES)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Linear matrices:",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",LINEAR_MATRICES% &
                & NUMBER_OF_LINEAR_MATRICES,ERR,ERROR,*999)
              DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",LINEAR_MATRICES%MATRICES(matrix_idx)%PTR% &
                  & UPDATE_MATRIX,ERR,ERROR,*999)
                IF(LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                  ELEMENT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                    & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                    & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                    & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                    & '(16X,8(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDDO !matrix_idx
            ENDIF
            RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
            IF(ASSOCIATED(RHS_VECTOR)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element RHS vector :",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",RHS_VECTOR%UPDATE_VECTOR,ERR,ERROR,*999)
              IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                ELEMENT_VECTOR=>RHS_VECTOR%ELEMENT_VECTOR
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                  & '("  Vector(:):",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equation matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF    

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE()")
#endif
       
    CALL EXITS("CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_CONDITION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian for the given element number for a finite element constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
        IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
          NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
          IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              SELECT CASE(NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR%JACOBIAN_CALCULATION_TYPE)
              CASE(CONSTRAINT_JACOBIAN_ANALYTIC_CALCULATED)
                ! None of these routines currently support calculating off diagonal terms for coupled problems,
                ! but when one does we will have to pass through the matrix_idx parameter
                IF(matrix_idx>1) THEN
                  CALL FLAG_ERROR("Analytic off-diagonal Jacobian calculation not implemented.",ERR,ERROR,*999)
                END IF
                SELECT CASE(CONSTRAINT_CONDITION%CLASS)
                CASE(CONSTRAINT_CONDITION_ELASTICITY_CLASS)
                  CALL ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_FLUID_MECHANICS_CLASS)
                  CALL FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_ELECTROMAGNETICS_CLASS)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_CLASSICAL_FIELD_CLASS)
                  CALL CLASSICAL_FIELD_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_BIOELECTRICS_CLASS)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_MODAL_CLASS)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_MULTI_PHYSICS_CLASS)
                  CALL MULTI_PHYSICS_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="Constraint condition class "//TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%CLASS,"*",ERR,ERROR))//" is not valid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE(CONSTRAINT_JACOBIAN_FINITE_DIFFERENCE_CALCULATED)
                CALL Constraint equationsSet_FiniteElementJacobianEvaluateFD(CONSTRAINT_CONDITION,ELEMENT_NUMBER,matrix_idx,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="Jacobian calculation type "//TRIM(NUMBER_TO_VSTRING(NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR% &
                  & JACOBIAN_CALCULATION_TYPE,"*",ERR,ERROR))//" is not valid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            END DO
          ELSE
            CALL FLAG_ERROR("Constraint equations nonlinear matrices is not associated.",ERR,ERROR,*999)
          END IF
        ELSE
          CALL FLAG_ERROR("Constraint equations matrices is not associated.",ERR,ERROR,*999)
        END IF
        IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_ELEMENT_MATRIX_OUTPUT) THEN
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element Jacobian matrix:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element Jacobian:",ERR,ERROR,*999)
          DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Jacobian number = ",matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update Jacobian = ",NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR% &
              & UPDATE_JACOBIAN,ERR,ERROR,*999)
            IF(NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR%UPDATE_JACOBIAN) THEN
              ELEMENT_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR%ELEMENT_JACOBIAN
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
              CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
              CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                & '(16X,8(X,E13.6))',ERR,ERROR,*999)
!!TODO: Write out the element residual???
            END IF
          END DO
        END IF
      ELSE
        CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    END IF

    CALL EXITS("CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1

  END SUBROUTINE CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix entries using finite differencing for a general finite element constraint condition.
  SUBROUTINE Constraint equationsSet_FiniteElementJacobianEvaluateFD(equationsSet,elementNumber,jacobianNumber,err,error,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: equationsSet  !<A pointer to the constraint condition to evaluate the element Jacobian for
    INTEGER(INTG), INTENT(IN) :: elementNumber  !<The element number to calculate the Jacobian for
    INTEGER(INTG), INTENT(IN) :: jacobianNumber  !<The Jacobian number to calculate when there are coupled problems
    INTEGER(INTG), INTENT(OUT) :: err  !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error  !<The error string
    !Local Variables
    TYPE(CONSTRAINT_TYPE), POINTER :: equations
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: elementsTopology
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: parameters
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowVariable,columnVariable
    TYPE(ELEMENT_VECTOR_TYPE) :: elementVector
    INTEGER(INTG) :: componentIdx,localNy,version,derivativeIdx,derivative,nodeIdx,node,column
    INTEGER(INTG) :: componentInterpolationType
    INTEGER(INTG) :: numberOfRows
    REAL(DP) :: delta,origDepVar

    CALL ENTERS("Constraint equationsSet_FiniteElementJacobianEvaluateFD",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(equations)) THEN
        equationsMatrices=>equations%CONSTRAINT_MATRICES
        nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
        nonlinearMapping=>equations%CONSTRAINT_MAPPING%NONLINEAR_MAPPING
        ! The first residual variable is always the row variable, which is the variable the
        ! residual is calculated for
        rowVariable=>nonlinearMapping%RESIDUAL_VARIABLES(1)%PTR
        ! For coupled problems this routine will be called multiple times if multiple Jacobians use finite
        ! differencing, so make sure we only calculate the residual vector once, to save time and because
        ! it would otherwise add together
        IF(nonlinearMatrices%ELEMENT_RESIDUAL_CALCULATED/=elementNumber) THEN
          CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
        END IF
        ! make a temporary copy of the unperturbed residuals
        elementVector=nonlinearMatrices%ELEMENT_RESIDUAL
        IF(jacobianNumber<=nonlinearMatrices%NUMBER_OF_JACOBIANS) THEN
          ! For coupled nonlinear problems there will be multiple Jacobians
          ! For this constraint condition, we calculate the residual for the row variable
          ! while pertubing parameters from the column variable.
          ! For non coupled problems these two variables will be the same
          columnVariable=>nonlinearMapping%RESIDUAL_VARIABLES(jacobianNumber)%PTR
          parameters=>columnVariable%PARAMETER_SETS%PARAMETER_SETS(FIELD_VALUES_SET_TYPE)%PTR%PARAMETERS  ! vector of dependent variables, basically
          numberOfRows=nonlinearMatrices%JACOBIANS(jacobianNumber)%PTR%ELEMENT_JACOBIAN%NUMBER_OF_ROWS
          IF(numberOfRows/=nonlinearMatrices%ELEMENT_RESIDUAL%NUMBER_OF_ROWS) THEN
            CALL FlagError("Element matrix number of rows does not match element residual vector size.",err,error,*999)
          END IF
          ! determine step size
          CALL DistributedVector_L2Norm(parameters,delta,err,error,*999)
          delta=(1.0_DP+delta)*1E-7_DP
          ! the actual finite differencing algorithm is about 4 lines but since the parameters are all
          ! distributed out, have to use proper field accessing routines..
          ! so let's just loop over component, node/el, derivative
          column=0  ! element jacobian matrix column number
          DO componentIdx=1,columnVariable%NUMBER_OF_COMPONENTS
            elementsTopology=>columnVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%ELEMENTS
            componentInterpolationType=columnVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE
            SELECT CASE (componentInterpolationType)
            CASE (FIELD_NODE_BASED_INTERPOLATION)
              basis=>elementsTopology%ELEMENTS(elementNumber)%BASIS
              DO nodeIdx=1,basis%NUMBER_OF_NODES
                node=elementsTopology%ELEMENTS(elementNumber)%ELEMENT_NODES(nodeIdx)
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                  derivative=elementsTopology%ELEMENTS(elementNumber)%ELEMENT_DERIVATIVES(derivativeIdx,nodeIdx)
                  version=elementsTopology%ELEMENTS(elementNumber)%elementVersions(derivativeIdx,nodeIdx)
                  localNy=columnVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                    & DERIVATIVES(derivative)%VERSIONS(version)
                  ! one-sided finite difference
                  CALL DISTRIBUTED_VECTOR_VALUES_GET(parameters,localNy,origDepVar,err,error,*999)
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar+delta,err,error,*999)
                  nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
                  CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar,err,error,*999)
                  column=column+1
                  nonlinearMatrices%JACOBIANS(jacobianNumber)%PTR%ELEMENT_JACOBIAN%MATRIX(1:numberOfRows,column)= &
                      & (nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR(1:numberOfRows)-elementVector%VECTOR(1:numberOfRows))/delta
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            CASE (FIELD_ELEMENT_BASED_INTERPOLATION)
              localNy=columnVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(elementNumber)
              ! one-sided finite difference
              CALL DISTRIBUTED_VECTOR_VALUES_GET(parameters,localNy,origDepVar,err,error,*999)
              CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar+delta,err,error,*999)
              nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
              CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar,err,error,*999)
              column=column+1
              nonlinearMatrices%JACOBIANS(jacobianNumber)%PTR%ELEMENT_JACOBIAN%MATRIX(1:numberOfRows,column)= &
                  & (nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR(1:numberOfRows)-elementVector%VECTOR(1:numberOfRows))/delta
            CASE DEFAULT
              CALL FLAG_ERROR("Unsupported type of interpolation.",err,error,*999)
            END SELECT
          END DO
          ! put the original residual back in
          nonlinearMatrices%ELEMENT_RESIDUAL=elementVector
        ELSE
          CALL FLAG_ERROR("Invalid Jacobian number of "//TRIM(NUMBER_TO_VSTRING(jacobianNumber,"*",err,error))// &
            & ". The number should be <= "//TRIM(NUMBER_TO_VSTRING(nonlinearMatrices%NUMBER_OF_JACOBIANS,"*",err,error))// &
            & ".",err,error,*999)
        END IF
      ELSE
        CALL FLAG_ERROR("Constraint condition equations is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",err,error,*999)
    END IF

    CALL EXITS("Constraint equationsSet_FiniteElementJacobianEvaluateFD")
    RETURN
999 CALL ERRORS("Constraint equationsSet_FiniteElementJacobianEvaluateFD",err,error)
    CALL EXITS("Constraint equationsSet_FiniteElementJacobianEvaluateFD")
    RETURN 1
  END SUBROUTINE Constraint equationsSet_FiniteElementJacobianEvaluateFD

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vector for the given element number for a finite element constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(ELEMENT_VECTOR_TYPE), POINTER :: ELEMENT_VECTOR
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(CONSTRAINT_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      SELECT CASE(CONSTRAINT_CONDITION%CLASS)
      CASE(CONSTRAINT_CONDITION_ELASTICITY_CLASS)
        CALL ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_FLUID_MECHANICS_CLASS)
        CALL FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_BIOELECTRICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_MULTI_PHYSICS_CLASS)
        CALL MULTI_PHYSICS_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Constraint condition class "//TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
        IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
          NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
          IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL_CALCULATED=ELEMENT_NUMBER
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_ELEMENT_MATRIX_OUTPUT) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element residual matrices and vectors:",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
              LINEAR_MATRICES=>CONSTRAINT_MATRICES%LINEAR_MATRICES
              IF(ASSOCIATED(LINEAR_MATRICES)) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Linear matrices:",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",LINEAR_MATRICES% &
                  & NUMBER_OF_LINEAR_MATRICES,ERR,ERROR,*999)
                DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",LINEAR_MATRICES%MATRICES(matrix_idx)%PTR% &
                    & UPDATE_MATRIX,ERR,ERROR,*999)
                  IF(LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                    ELEMENT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                      & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                      & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                      & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                    CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                      & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                      & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                      & '(16X,8(X,E13.6))',ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
              ENDIF
              DYNAMIC_MATRICES=>CONSTRAINT_MATRICES%DYNAMIC_MATRICES
              IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Dynamnic matrices:",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",DYNAMIC_MATRICES% &
                  & NUMBER_OF_DYNAMIC_MATRICES,ERR,ERROR,*999)
                DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR% &
                    & UPDATE_MATRIX,ERR,ERROR,*999)
                  IF(DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                    ELEMENT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                      & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                      & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                      & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                    CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                      & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                      & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                      & '(16X,8(X,E13.6))',ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
              ENDIF
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element residual vector:",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",NONLINEAR_MATRICES%UPDATE_RESIDUAL,ERR,ERROR,*999)
              IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
                ELEMENT_VECTOR=>NONLINEAR_MATRICES%ELEMENT_RESIDUAL
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                  & '("  Vector(:):",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
              RHS_VECTOR=>CONSTRAINT_MATRICES%RHS_VECTOR
              IF(ASSOCIATED(RHS_VECTOR)) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element RHS vector :",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",RHS_VECTOR%UPDATE_VECTOR,ERR,ERROR,*999)
                IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                  ELEMENT_VECTOR=>RHS_VECTOR%ELEMENT_VECTOR
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                    & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                    & '("  Vector(:)    :",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equation nonlinear matrices not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equation matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF    
       
    CALL EXITS("CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Initialises an constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<The pointer to the constraint condition to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONSTRAINT_CONDITION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CALL FLAG_ERROR("Constraint condition is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(CONSTRAINT_CONDITION,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint condition.",ERR,ERROR,*999)
      CONSTRAINT_CONDITION%USER_NUMBER=0
      CONSTRAINT_CONDITION%GLOBAL_NUMBER=0
      CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FINISHED=.FALSE.
      NULLIFY(CONSTRAINT_CONDITION%CONSTRAINT_CONDITIONS)
      NULLIFY(CONSTRAINT_CONDITION%REGION)
      CONSTRAINT_CONDITION%CLASS=CONSTRAINT_CONDITION_NO_CLASS
      CONSTRAINT_CONDITION%TYPE=CONSTRAINT_CONDITION_NO_TYPE
      CONSTRAINT_CONDITION%SUBTYPE=CONSTRAINT_CONDITION_NO_SUBTYPE
      CONSTRAINT_CONDITION%SOLUTION_METHOD=0
      CALL CONSTRAINT_CONDITION_GEOMETRY_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*999)
      CALL CONSTRAINT_CONDITION_DEPENDENT_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*999)
      CALL CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*999)
      NULLIFY(CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS)
      NULLIFY(CONSTRAINT_CONDITION%BOUNDARY_CONDITIONS)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_INITIALISE")
    RETURN
999 CALL CONSTRAINT_CONDITION_FINALISE(CONSTRAINT_CONDITION,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_CONDITION_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise the geometry for an constraint condition
  SUBROUTINE CONSTRAINT_CONDITION_GEOMETRY_FINALISE(CONSTRAINT_CONDITION_GEOMETRY,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_GEOMETRY_TYPE) :: CONSTRAINT_CONDITION_GEOMETRY !<A pointer to the constraint condition geometry to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITION_GEOMETRY_FINALISE",ERR,ERROR,*999)
    
    NULLIFY(CONSTRAINT_CONDITION_GEOMETRY%GEOMETRIC_FIELD)
    NULLIFY(CONSTRAINT_CONDITION_GEOMETRY%FIBRE_FIELD)
       
    CALL EXITS("CONSTRAINT_CONDITION_GEOMETRY_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_GEOMETRY_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_GEOMETRY_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_GEOMETRY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the geometry for an equation set
  SUBROUTINE CONSTRAINT_CONDITION_GEOMETRY_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to initialise the geometry for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_CONDITION_GEOMETRY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_CONDITION%GEOMETRY%CONSTRAINT_CONDITION=>CONSTRAINT_CONDITION
      NULLIFY(CONSTRAINT_CONDITION%GEOMETRY%GEOMETRIC_FIELD)
      NULLIFY(CONSTRAINT_CONDITION%GEOMETRY%FIBRE_FIELD)
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_GEOMETRY_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_GEOMETRY_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_GEOMETRY_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_GEOMETRY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of a dependent variables for an constraint condition. \see OPENCMISS::CMISSConstraintConditionDependentCreateFinish
  SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_CREATE_FINISH(CONSTRAINT_CONDITION,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_SETUP_TYPE) :: CONSTRAINT_CONDITION_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD

    CALL ENTERS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      IF(CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("Constraint condition dependent has already been finished",ERR,ERROR,*999)
      ELSE
        !Initialise the setup
        CALL CONSTRAINT_CONDITION_SETUP_INITIALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        CONSTRAINT_CONDITION_SETUP_INFO%SETUP_TYPE=CONSTRAINT_CONDITION_SETUP_DEPENDENT_TYPE
        CONSTRAINT_CONDITION_SETUP_INFO%ACTION_TYPE=CONSTRAINT_CONDITION_SETUP_FINISH_ACTION
        DEPENDENT_FIELD=>CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          CONSTRAINT_CONDITION_SETUP_INFO%FIELD_USER_NUMBER=DEPENDENT_FIELD%USER_NUMBER
          CONSTRAINT_CONDITION_SETUP_INFO%FIELD=>DEPENDENT_FIELD
          !Finish constraint condition specific setup
          CALL CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Constraint condition dependent dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
        !Finalise the setup
        CALL CONSTRAINT_CONDITION_SETUP_FINALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        !Finish the constraint condition creation
        CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of dependent variables for an constraint condition. \see OPENCMISS::CMISSConstraintConditionDependentCreateStart
  SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_CREATE_START(CONSTRAINT_CONDITION,DEPENDENT_FIELD_USER_NUMBER,DEPENDENT_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to start the creation of a dependent field on
    INTEGER(INTG), INTENT(IN) :: DEPENDENT_FIELD_USER_NUMBER !<The user specified dependent field number
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<If associated on entry, a pointer to the user created dependent field which has the same user number as the specified dependent field user number. If not associated on entry, on exit, a pointer to the created dependent field for the constraint condition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(CONSTRAINT_CONDITION_SETUP_TYPE) :: CONSTRAINT_CONDITION_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,DEPENDENT_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      IF(CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("The constraint condition dependent has been finished.",ERR,ERROR,*999)
      ELSE
        REGION=>CONSTRAINT_CONDITION%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
            !Check the dependent field has been finished
            IF(DEPENDENT_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(DEPENDENT_FIELD_USER_NUMBER/=DEPENDENT_FIELD%USER_NUMBER) THEN
                LOCAL_ERROR="The specified dependent field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not match the user number of the specified dependent field of "// &
                  & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              DEPENDENT_FIELD_REGION=>DEPENDENT_FIELD%REGION
              IF(ASSOCIATED(DEPENDENT_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the constraint condition
                IF(DEPENDENT_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  LOCAL_ERROR="Invalid region setup. The specified dependent field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and the specified constraint condition has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                !Check the specified dependent field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>CONSTRAINT_CONDITION%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,DEPENDENT_FIELD%DECOMPOSITION)) THEN
                    CALL FLAG_ERROR("The specified dependent field does not have the same decomposition as the geometric "// &
                      & "field for the specified constraint condition.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The geometric field is not associated for the specified constraint condition.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified dependent field region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(FIELD)
            CALL FIELD_USER_NUMBER_FIND(DEPENDENT_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            IF(ASSOCIATED(FIELD)) THEN
              LOCAL_ERROR="The specified dependent field user number of "// &
                & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                & " has already been used to create a field on region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
            CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.TRUE.
          ENDIF
          !Initialise the setup
          CALL CONSTRAINT_CONDITION_SETUP_INITIALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
          CONSTRAINT_CONDITION_SETUP_INFO%SETUP_TYPE=CONSTRAINT_CONDITION_SETUP_DEPENDENT_TYPE
          CONSTRAINT_CONDITION_SETUP_INFO%ACTION_TYPE=CONSTRAINT_CONDITION_SETUP_START_ACTION
          CONSTRAINT_CONDITION_SETUP_INFO%FIELD_USER_NUMBER=DEPENDENT_FIELD_USER_NUMBER
          CONSTRAINT_CONDITION_SETUP_INFO%FIELD=>DEPENDENT_FIELD
          !Start the constraint condition specfic solution setup
          CALL CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the setup
          CALL CONSTRAINT_CONDITION_SETUP_FINALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
          !Set pointers
          IF(CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
            DEPENDENT_FIELD=>CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD
          ELSE
            CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD=>DEPENDENT_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equation set region is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint equations_set is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_START")
    RETURN
999 CALL CONSTRAINT_CONDITION_DEPENDENT_FINALISE(CONSTRAINT_CONDITION%DEPENDENT,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_START",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_CREATE_START")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_CREATE_START

  !
  !================================================================================================================================
  !
  
  !>Destroy the dependent variables for an constraint condition. \see OPENCMISS::CMISSConstraintConditionDependentDestroy
  SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_DESTROY(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<The pointer to the constraint condition to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITION_DEPENDENT_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CALL CONSTRAINT_CONDITION_DEPENDENT_FINALISE(CONSTRAINT_CONDITION%DEPENDENT,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_DESTROY")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_DEPENDENT_DESTROY",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_DESTROY")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_DESTROY
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the dependent variables for an equation set and deallocates all memory.
  SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_FINALISE(CONSTRAINT_CONDITION_DEPENDENT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_DEPENDENT_TYPE) :: CONSTRAINT_CONDITION_DEPENDENT !<The pointer to the constraint condition
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITION_DEPENDENT_FINALISE",ERR,ERROR,*999)

    NULLIFY(CONSTRAINT_CONDITION_DEPENDENT%CONSTRAINT_CONDITION)
    CONSTRAINT_CONDITION_DEPENDENT%DEPENDENT_FINISHED=.FALSE.
    CONSTRAINT_CONDITION_DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.FALSE.
    NULLIFY(CONSTRAINT_CONDITION_DEPENDENT%DEPENDENT_FIELD)
    
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_DEPENDENT_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the dependent variables for a constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to initialise the dependent field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_CONDITION_DEPENDENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_CONDITION%DEPENDENT%CONSTRAINT_CONDITION=>CONSTRAINT_CONDITION
      CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FINISHED=.FALSE.
      CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.FALSE.
      NULLIFY(CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD)
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_DEPENDENT_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_DEPENDENT_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_DEPENDENT_INITIALISE

  !
  !================================================================================================================================
  !
  !>Finalises the dependent variables for an equation set and deallocates all memory.
  SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_FINALISE(CONSTRAINT_CONDITION_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_TYPE) :: CONSTRAINT_CONDITION_FIELD !<The pointer to the constraint condition
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_FINALISE",ERR,ERROR,*999)

    NULLIFY(CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION)
    CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_FINISHED=.FALSE.
    CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_AUTO_CREATED=.FALSE.
    NULLIFY(CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_FIELD)
    
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_FINALISE
  
  !
  !================================================================================================================================
  !
  !>Initialises the constraint condition field for a constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_INITIALISE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to initialise the dependent field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION=>CONSTRAINT_CONDITION
      CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_FINISHED=.FALSE.
      CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_AUTO_CREATED=.TRUE.
      NULLIFY(CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FIELD%CONSTRAINT_CONDITION_FIELD_FIELD)
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CONDITION_FIELD_INITIALISE

  !
  !================================================================================================================================
  !



  !>Sets up the specifices for an equation set.
  SUBROUTINE CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to perform the setup on
    TYPE(CONSTRAINT_CONDITION_SETUP_TYPE), INTENT(INOUT) :: CONSTRAINT_CONDITION_SETUP_INFO !<The constraint condition setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_CONDITION_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      SELECT CASE(CONSTRAINT_CONDITION%CLASS)
      CASE(CONSTRAINT_CONDITION_ELASTICITY_CLASS)
        CALL ELASTICITY_CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_FLUID_MECHANICS_CLASS)
        CALL FLUID_MECHANICS_CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_BIOELECTRICS_CLASS)
        IF(CONSTRAINT_CONDITION%TYPE == CONSTRAINT_CONDITION_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
          CALL MONODOMAIN_EQUATION_CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        ELSE
          CALL BIOELECTRIC_CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        END IF
      CASE(CONSTRAINT_CONDITION_FITTING_CLASS)
        CALL FITTING_CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(CONSTRAINT_CONDITION_MULTI_PHYSICS_CLASS)
        CALL MULTI_PHYSICS_CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Constraint condition class "//TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_SETUP")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_SETUP",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_SETUP")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_SETUP

  !
  !================================================================================================================================
  !

 !>Finish the creation of equations for the constraint condition. \see OPENCMISS::CMISSConstraintConditionConstraint equationsCreateFinish
  SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CREATE_FINISH(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to finish the creation of the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_SETUP_TYPE) :: CONSTRAINT_CONDITION_SETUP_INFO
    
    CALL ENTERS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      !Initialise the setup
      CALL CONSTRAINT_CONDITION_SETUP_INITIALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
      CONSTRAINT_CONDITION_SETUP_INFO%SETUP_TYPE=CONSTRAINT_CONDITION_SETUP_CONSTRAINT_TYPE
      CONSTRAINT_CONDITION_SETUP_INFO%ACTION_TYPE=CONSTRAINT_CONDITION_SETUP_FINISH_ACTION
      !Finish the equations specific solution setup.
      CALL CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the setup
      CALL CONSTRAINT_CONDITION_SETUP_FINALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of equations for the equation set. \see CMISSConstraintConditionConstraint equationsCreateStart
  !>Default values set for the EQUATIONS's attributes are:
  !>- OUTPUT_TYPE: 0 (CONSTRAINT_CONDITION_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (CONSTRAINT_CONDITION_SPARSE_MATRICES)
  !>- NONLINEAR_JACOBIAN_TYPE: 0
  !>- INTERPOLATION: null
  !>- LINEAR_DATA: null 
  !>- NONLINEAR_DATA: null
  !>- TIME_DATA: null
  !>- CONSTRAINT_MAPPING:  
  !>- CONSTRAINT_MATRICES:  
  SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CREATE_START(CONSTRAINT_CONDITION,EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to create equations for
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS !<On exit, a pointer to the created equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_SETUP_TYPE) :: CONSTRAINT_CONDITION_SETUP_INFO

    CALL ENTERS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      IF(ASSOCIATED(EQUATIONS)) THEN
        CALL FLAG_ERROR("Constraint equations is already associated.",ERR,ERROR,*999)
      ELSE
        !Initialise the setup
        CALL CONSTRAINT_CONDITION_SETUP_INITIALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        CONSTRAINT_CONDITION_SETUP_INFO%SETUP_TYPE=CONSTRAINT_CONDITION_SETUP_CONSTRAINT_TYPE
        CONSTRAINT_CONDITION_SETUP_INFO%ACTION_TYPE=CONSTRAINT_CONDITION_SETUP_START_ACTION
        !Start the constraint condition specific solution setup
        CALL CONSTRAINT_CONDITION_SETUP(CONSTRAINT_CONDITION,CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        !Finalise the setup
        CALL CONSTRAINT_CONDITION_SETUP_FINALISE(CONSTRAINT_CONDITION_SETUP_INFO,ERR,ERROR,*999)
        !Return the pointer
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_START")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_START",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_CREATE_START")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the equations for an constraint condition. \see OPENCMISS::CMISSConstraintConditionConstraint equationsDestroy
  SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_DESTROY(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to destroy the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITION_CONSTRAINT_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      IF(ASSOCIATED(CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS)) THEN
        CALL CONSTRAINT_FINALISE(CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Constraint condition equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_DESTROY")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_CONSTRAINT_DESTROY",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_CONSTRAINT_DESTROY")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_CONSTRAINT_DESTROY

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(CONSTRAINT_EQUATIONS%CONSTRAINT_FINISHED) THEN
          SELECT CASE(CONSTRAINT_EQUATIONS%LINEARITY)
          CASE(CONSTRAINT_CONDITION_LINEAR)
            SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
            CASE(CONSTRAINT_CONDITION_STATIC)
              SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
              CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                CALL CONSTRAINT_CONDITION_ASSEMBLE_STATIC_LINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition solution method of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC,CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
              SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
              CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                CALL CONSTRAINT_CONDITION_ASSEMBLE_DYNAMIC_LINEAR_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition solution method of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The equations time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(CONSTRAINT_CONDITION_NONLINEAR)
            SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
            CASE(CONSTRAINT_CONDITION_STATIC)
              SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
              CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                CALL CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_STATIC_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition solution method  of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC,CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
! sebk 15/09/09
              SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
              CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                CALL CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_DYNAMIC_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition solution method  of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(CONSTRAINT_TIME_STEPPING)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The constraint condition time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(CONSTRAINT_CONDITION_NONLINEAR_BCS)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations linearity of "// &
              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Constraint equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint condition equations is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for an static constraint condition using the finite element method
  SUBROUTINE CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_STATIC_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
  
    CALL ENTERS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      DEPENDENT_FIELD=>CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
!!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,CONSTRAINT_MATRICES_JACOBIAN_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for constraint conditionup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for constraint conditionup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_MATRIX_OUTPUT) THEN
              CALL CONSTRAINT_MATRICES_JACOBIAN_OUTPUT(GENERAL_OUTPUT_TYPE,CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF            
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_STATIC_FEM")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_STATIC_FEM",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_STATIC_FEM")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_STATIC_FEM

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for an dynamic constraint condition using the finite element method
  SUBROUTINE CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_DYNAMIC_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
  
    CALL ENTERS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_DYNAMIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      DEPENDENT_FIELD=>CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
!!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,CONSTRAINT_MATRICES_JACOBIAN_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for constraint conditionup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for constraint conditionup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_JACOBIAN_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_JACOBIAN_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
             !Output equations matrices and RHS vector if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_MATRIX_OUTPUT) THEN
              CALL CONSTRAINT_MATRICES_JACOBIAN_OUTPUT(GENERAL_OUTPUT_TYPE,CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF            
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_DYNAMIC_FEM")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_DYNAMIC_FEM",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_DYNAMIC_FEM")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_JACOBIAN_EVALUATE_DYNAMIC_FEM

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an constraint condition.
  SUBROUTINE CONSTRAINT_CONDITION_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: residual_variable_idx
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(CONSTRAINT_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: RESIDUAL_PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: RESIDUAL_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(CONSTRAINT_EQUATIONS%CONSTRAINT_FINISHED) THEN
          SELECT CASE(CONSTRAINT_EQUATIONS%LINEARITY)
          CASE(CONSTRAINT_CONDITION_LINEAR)
            CALL FLAG_ERROR("Can not evaluate a residual for linear equations.",ERR,ERROR,*999)
          CASE(CONSTRAINT_CONDITION_NONLINEAR)
            SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
            CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC,CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
              SELECT CASE(CONSTRAINT_CONDITION%SOLUTION_METHOD)
              CASE(CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD)
                CALL CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_DYNAMIC_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition solution method  of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The constraint condition time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(CONSTRAINT_CONDITION_NONLINEAR_BCS)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations linearity of "// &
              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Update the residual parameter set if it exists
          CONSTRAINT_MAPPING=>CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING
          IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
            NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
              DO residual_variable_idx=1,NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
                RESIDUAL_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(residual_variable_idx)%PTR
                IF(ASSOCIATED(RESIDUAL_VARIABLE)) THEN
                  RESIDUAL_PARAMETER_SET=>RESIDUAL_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_RESIDUAL_SET_TYPE)%PTR
                  IF(ASSOCIATED(RESIDUAL_PARAMETER_SET)) THEN
                    !Residual parameter set exists
                    CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
                    IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
                      NONLINEAR_MATRICES=>CONSTRAINT_MATRICES%NONLINEAR_MATRICES
                      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                        !Copy the residual vector to the residuals parameter set.
                        CALL DISTRIBUTED_VECTOR_COPY(NONLINEAR_MATRICES%RESIDUAL,RESIDUAL_PARAMETER_SET%PARAMETERS,1.0_DP, &
                          & ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("Constraint equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Constraint equations equations matrices is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  LOCAL_ERROR="Nonlinear mapping residual variable for residual variable index "// &
                    & TRIM(NUMBER_TO_VSTRING(residual_variable_idx,"*",ERR,ERROR))//" is not associated."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !residual_variable_idx
            ELSE
              CALL FLAG_ERROR("Constraint equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint equations equations mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint condition equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_CONDITION_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an dynamic constraint condition using the finite element method
  SUBROUTINE CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_DYNAMIC_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
 
    CALL ENTERS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_DYNAMIC_FEM",ERR,ERROR,*999)

    NULLIFY(ELEMENTS_MAPPING)
    NULLIFY(EQUATIONS)
    NULLIFY(CONSTRAINT_MATRICES)
    NULLIFY(DEPENDENT_FIELD)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      DEPENDENT_FIELD=>CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,CONSTRAINT_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for constraint conditionup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for constraint conditionup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
             !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_MATRIX_OUTPUT) THEN
              CALL CONSTRAINT_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_DYNAMIC_FEM")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_DYNAMIC_FEM",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_DYNAMIC_FEM")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_DYNAMIC_FEM

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an static constraint condition using the finite element method
  SUBROUTINE CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_STATIC_FEM(CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MATRICES_TYPE), POINTER :: CONSTRAINT_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
 
    CALL ENTERS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      DEPENDENT_FIELD=>CONSTRAINT_CONDITION%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          CONSTRAINT_MATRICES=>CONSTRAINT_EQUATIONS%CONSTRAINT_MATRICES
          IF(ASSOCIATED(CONSTRAINT_MATRICES)) THEN
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL CONSTRAINT_MATRICES_VALUES_INITIALISE(CONSTRAINT_MATRICES,CONSTRAINT_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL CONSTRAINT_MATRICES_ELEMENT_INITIALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for constraint conditionup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for constraint conditionup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL CONSTRAINT_MATRICES_ELEMENT_CALCULATE(CONSTRAINT_MATRICES,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_CONDITION_FINITE_ELEMENT_RESIDUAL_EVALUATE(CONSTRAINT_CONDITION,ne,ERR,ERROR,*999)
              CALL CONSTRAINT_MATRICES_ELEMENT_ADD(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL CONSTRAINT_MATRICES_ELEMENT_FINALISE(CONSTRAINT_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_MATRIX_OUTPUT) THEN
              CALL CONSTRAINT_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,CONSTRAINT_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(CONSTRAINT_EQUATIONS%OUTPUT_TYPE>=CONSTRAINT_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_STATIC_FEM")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_STATIC_FEM",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_STATIC_FEM")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_RESIDUAL_EVALUATE_STATIC_FEM
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for an constraint condition. \see OPENCMISS::CMISSConstraintConditionSolutionMethodSet
  SUBROUTINE CONSTRAINT_CONDITION_SOLUTION_METHOD_SET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The constraint condition solution method to set \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_CONDITION_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      IF(CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Constraint condition has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(CONSTRAINT_CONDITION%CLASS)
        CASE(CONSTRAINT_CONDITION_ELASTICITY_CLASS)
          CALL ELASTICITY_CONSTRAINT_CONDITION_SOLUTION_METHOD_SET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE(CONSTRAINT_CONDITION_FLUID_MECHANICS_CLASS)
          CALL FLUID_MECHANICS_CONSTRAINT_CONDITION_SOLUTION_METHOD_SET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE(CONSTRAINT_CONDITION_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(CONSTRAINT_CONDITION_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_CONSTRAINT_CONDITION_SOLUTION_METHOD_SET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE(CONSTRAINT_CONDITION_BIOELECTRICS_CLASS)
          IF(CONSTRAINT_CONDITION%TYPE == CONSTRAINT_CONDITION_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
            CALL MONODOMAIN_EQUATION_CONSTRAINT_CONDITION_SOLUTION_METHOD_SET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*999)
          ELSE
            CALL BIOELECTRIC_CONSTRAINT_CONDITION_SOLUTION_METHOD_SET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*999)
          END IF
        CASE(CONSTRAINT_CONDITION_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(CONSTRAINT_CONDITION_MULTI_PHYSICS_CLASS)
          CALL MULTI_PHYSICS_CONSTRAINT_CONDITION_SOLUTION_METHOD_SET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Constraint condition class "//TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%CLASS,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_CONDITION_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_SOLUTION_METHOD_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the solution method for an constraint condition. \see OPENCMISS::CMISSConstraintConditionSolutionMethodGet
  SUBROUTINE CONSTRAINT_CONDITION_SOLUTION_METHOD_GET(CONSTRAINT_CONDITION,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<A pointer to the constraint condition to get the solution method for
    INTEGER(INTG), INTENT(OUT) :: SOLUTION_METHOD !<On return, the constraint condition solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITION_SOLUTION_METHOD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
      IF(CONSTRAINT_CONDITION%CONSTRAINT_CONDITION_FINISHED) THEN
        SOLUTION_METHOD=CONSTRAINT_CONDITION%SOLUTION_METHOD
      ELSE
        CALL FLAG_ERROR("Constraint condition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_CONDITION_SOLUTION_METHOD_GET")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_SOLUTION_METHOD_GET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_SOLUTION_METHOD_GET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_SOLUTION_METHOD_GET
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in CONSTRAINT_CONDITION a pointer to the constraint condition identified by USER_NUMBER in the given REGION. If no constraint condition with that USER_NUMBER exists CONSTRAINT_CONDITION is left nullified.
  SUBROUTINE CONSTRAINT_CONDITION_USER_NUMBER_FIND(USER_NUMBER,REGION,CONSTRAINT_CONDITION,ERR,ERROR,*)

    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find the equation set
    TYPE(REGION_TYPE), POINTER :: REGION !<The region to find the constraint condition in
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION !<On return, a pointer to the constraint condition if an constraint condition with the specified user number exists in the given region. If no equation set with the specified number exists a NULL pointer is returned. The pointer must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_CONDITION_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
        CALL FLAG_ERROR("Constraint condition is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(CONSTRAINT_CONDITION)
        IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS)) THEN
          equations_set_idx=1
          DO WHILE(equations_set_idx<=REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS.AND..NOT.ASSOCIATED(CONSTRAINT_CONDITION))
            IF(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
              CONSTRAINT_CONDITION=>REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(equations_set_idx)%PTR
            ELSE
              equations_set_idx=equations_set_idx+1
            ENDIF
          ENDDO
        ELSE
          LOCAL_ERROR="The constraint conditions on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_CONDITION_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITION_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITION_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises all constraint conditions on a region and deallocates all memory.
  SUBROUTINE CONSTRAINT_CONDITIONS_FINALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finalise the problems for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS)) THEN
        DO WHILE(REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS>0)
          CALL CONSTRAINT_CONDITION_DESTROY(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS(1)%PTR,ERR,ERROR,*999)
        ENDDO !problem_idx
        DEALLOCATE(REGION%CONSTRAINT_CONDITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_CONDITIONS_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_CONDITIONS_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITIONS_FINALISE")
    RETURN 1   
  END SUBROUTINE CONSTRAINT_CONDITIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises all constraint conditions on a region.
  SUBROUTINE CONSTRAINT_CONDITIONS_INITIALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the constraint conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_CONDITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS)) THEN
        CALL FLAG_ERROR("Region already has associated constraint conditions",ERR,ERROR,*998)
      ELSE
!!TODO: Inherit any constraint conditions from the parent region???
        ALLOCATE(REGION%CONSTRAINT_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region constraint conditions",ERR,ERROR,*999)
        REGION%CONSTRAINT_CONDITIONS%REGION=>REGION
        REGION%CONSTRAINT_CONDITIONS%NUMBER_OF_CONSTRAINT_CONDITIONS=0
        NULLIFY(REGION%CONSTRAINT_CONDITIONS%CONSTRAINT_CONDITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("CONSTRAINT_CONDITIONS_INITIALISE")
    RETURN
999 IF(ASSOCIATED(REGION%CONSTRAINT_CONDITIONS)) DEALLOCATE(REGION%CONSTRAINT_CONDITIONS)
998 CALL ERRORS("CONSTRAINT_CONDITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CONDITIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CONDITIONS_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE CONSTRAINT_CONDITION_ROUTINES
