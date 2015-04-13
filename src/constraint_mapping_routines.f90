!> \file
!> \author Chris Bradley
!> \brief This module contains all constraint mapping routines.
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

!>This module contains all constraint mapping routines.
MODULE CONSTRAINT_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE CONSTRAINT_CONDITIONS_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables


  INTERFACE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET
    MODULE PROCEDURE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_1
    MODULE PROCEDURE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_2
  END INTERFACE !CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET
  
  INTERFACE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET
    MODULE PROCEDURE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1
    MODULE PROCEDURE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2
  END INTERFACE !CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET
  
  PUBLIC CONSTRAINT_MAPPING_CREATE_FINISH,CONSTRAINT_MAPPING_CREATE_START,CONSTRAINT_MAPPING_DESTROY, &
    & CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET,CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET, &
    & CONSTRAINT_MAPPING_LINEAR_MATRICES_COEFFS_SET,CONSTRAINT_MAPPING_JACOBIAN_COEFF_SET, &
    & CONSTRAINT_MAPPING_RHS_COEFF_SET,CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET
 
  PUBLIC CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET
 
  PUBLIC CONSTRAINT_MAPPING_DYNAMIC_MATRICES_TRANSPOSE_SET,CONSTRAINT_MAPPING_LINEAR_MATRICES_TRANSPOSE_SET, &
    & CONSTRAINT_MAPPING_JACOBIAN_MATRICES_TRANSPOSE_SET


CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates constraint mapping
  SUBROUTINE CONSTRAINT_MAPPING_CALCULATE(CONSTRAINT_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to calculate the mapping for
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: row_idx,column_idx,dof_idx,matrix_idx,number_of_constraint_matrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD,DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,LAGRANGE_VARIABLE
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_DEPENDENT_TYPE), POINTER :: CONSTRAINT_DEPENDENT
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(CONSTRAINT_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MAPPING_CALCULATE",ERR,ERROR,*999)

   IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
     CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
     IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
       CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
       IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
         CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
         SELECT CASE(CONSTRAINT_CONDITION%METHOD)
         CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
           LAGRANGE=>CONSTRAINT_CONDITION%LAGRANGE
           IF(ASSOCIATED(LAGRANGE)) THEN
             CONSTRAINT_DEPENDENT=>CONSTRAINT_CONDITION%DEPENDENT
             IF(ASSOCIATED(CONSTRAINT_DEPENDENT)) THEN
               !Set the Lagrange variable information
               LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
               DEPENDENT_FIELD=>CONSTRAINT_DEPENDENT%DEPENDENT_FIELD
               NULLIFY(LAGRANGE_VARIABLE)
               CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                 & ERR,ERROR,*999)
               CONSTRAINT_MAPPING%LAGRANGE_VARIABLE_TYPE=CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE
               CONSTRAINT_MAPPING%LAGRANGE_VARIABLE=>LAGRANGE_VARIABLE
               !Set the number of columns in the constraint matrices
               CONSTRAINT_MAPPING%NUMBER_OF_COLUMNS=LAGRANGE_VARIABLE%NUMBER_OF_DOFS
               CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS=LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
               CONSTRAINT_MAPPING%NUMBER_OF_GLOBAL_COLUMNS=LAGRANGE_VARIABLE%NUMBER_OF_GLOBAL_DOFS
               !Set the column dofs mapping
               CONSTRAINT_MAPPING%COLUMN_DOFS_MAPPING=>LAGRANGE_VARIABLE%DOMAIN_MAPPING

                !Calculate dynamic mappings
                IF(CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES>0) THEN                  
                  CALL CONSTRAINT_MAPPING_DYNAMIC_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*999)
                  DYNAMIC_MAPPING=>CONSTRAINT_MAPPING%DYNAMIC_MAPPING
                  IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                    DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                    DYNAMIC_MAPPING%STIFFNESS_MATRIX_NUMBER=CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER
                    DYNAMIC_MAPPING%DAMPING_MATRIX_NUMBER=CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER
                    DYNAMIC_MAPPING%MASS_MATRIX_NUMBER=CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER
                    ALLOCATE(DYNAMIC_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Lagrange dof to column map.",ERR,ERROR,*999)
                    !1-1 mapping for now
                    DO dof_idx=1,LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                      column_idx=LAGRANGE_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                      DYNAMIC_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(dof_idx)=column_idx
                    ENDDO
                    ALLOCATE(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP,STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices map.",ERR,ERROR,*999)
                    CALL CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT( &
                      & DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*999)
                    DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_TYPE=CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE
                    CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,CREATE_VALUES_CACHE%DYNAMIC_VARIABLE_TYPE,FIELD_VARIABLE, &
                      & ERR,ERROR,*999)
                    DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE=>FIELD_VARIABLE
                    DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES= &
                      & CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                    IF(CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE/=0) & 
                      & DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=-1
                    !Allocate and initialise the variable to constraint matrices maps
                    FIELD_VARIABLE=>DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      IF(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES==-1) THEN
                        !!TODO: check if this can be removed and just allocate those variables that are actually used
                        ALLOCATE(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP=0
                      ELSE IF(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES>0) THEN
                        ALLOCATE(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS( &
                          & DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) &
                          & CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps constraint matrix numbers.", &
                          & ERR,ERROR,*999)
                        DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS=0
                        IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                          DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS( &
                            & CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)=CREATE_VALUES_CACHE% &
                            & DYNAMIC_STIFFNESS_MATRIX_NUMBER
                        ENDIF
                        IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                          DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS( &
                            & CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)=CREATE_VALUES_CACHE% &
                            & DYNAMIC_DAMPING_MATRIX_NUMBER
                        ENDIF
                        IF(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                          DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS( &
                            & CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER)=CREATE_VALUES_CACHE% &
                            & DYNAMIC_MASS_MATRIX_NUMBER
                        ENDIF
                        ALLOCATE(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                        !1-1 mappings for now.
                          row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                          DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP(dof_idx)=row_idx
                        ENDDO !dof_idx
                      ENDIF
                    ENDIF
                    !Allocate and initialise the constraint matrix to variable maps types
                    ALLOCATE(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(DYNAMIC_MAPPING% &
                      & NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrix rows to variable maps.", &
                      & ERR,ERROR,*999)
                    !Loop over the constraint matrices and calculate the row mappings
                    !The pointers below have been checked for association above.
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      number_of_constraint_matrices=DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      !Number of constraint matrices whose rows/columns are related to Dependent/Lagrange variables and not Lagrange/Lagrange variables (last constraint matrix is Lagrange/Lagrange (Penalty matrix)
                      number_of_constraint_matrices=DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES/2 
                    ENDSELECT
                    DO matrix_idx=1,number_of_constraint_matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP( &
                        & matrix_idx),ERR,ERROR,*999)
                      EQUATIONS_SET=>CONSTRAINT_DEPENDENT%EQUATIONS_SET
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        FIELD_VARIABLE=>DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%EQUATIONS_SET=>EQUATIONS_SET
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE= &
                            & FIELD_VARIABLE%VARIABLE_TYPE
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%MATRIX_COEFFICIENT= &
                            & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(matrix_idx)
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE(matrix_idx)
                           !Set the number of rows
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_DOFS
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row to variable dofs map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)% &
                              & ROW_TO_VARIABLE_DOF_MAP(row_idx)=dof_idx
                          ENDDO !dof_idx
                        ELSE
                          CALL FLAG_ERROR("Dependent variable could not be found.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations set could not be found.",ERR,ERROR,*999)
                      ENDIF
                    ENDDO !matrix_idx
                    !The pointers below have been checked for association above.
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      !Sets up the Lagrange-(Penalty) constraint matrix mapping and calculate the row mappings
                      DO matrix_idx=number_of_constraint_matrices,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                        !Initialise and setup the constraint matrix
                        CALL CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP( &
                          & matrix_idx),ERR,ERROR,*999)
                        CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                          & ERR,ERROR,*999)
                        FIELD_VARIABLE=>LAGRANGE_VARIABLE
                        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
                        IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
                          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%CONSTRAINT_EQUATIONS=> &
                              & CONSTRAINT_EQUATIONS
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE% &
                              & VARIABLE_TYPE
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%MATRIX_COEFFICIENT= &
                              & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(matrix_idx)
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                              & CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE(matrix_idx)
                              !Set the number of rows
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS= &
                              & FIELD_VARIABLE%NUMBER_OF_DOFS
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                              & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                              & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                            !Set the row mapping
                            DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_DOFS_MAPPING=> &
                              & FIELD_VARIABLE%DOMAIN_MAPPING
                            ALLOCATE(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP( &
                              & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                            !1-1 mapping for now
                            DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                              row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                              DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)% &
                                & ROW_TO_VARIABLE_DOF_MAP(row_idx)=dof_idx
                            ENDDO !dof_idx
                          ELSE
                            CALL FLAG_ERROR("Lagrange variable could not be found.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Constraint equations could not be found.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO
                    ENDSELECT
                  ELSE
                    CALL FLAG_ERROR("Dynamic mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
 
               !Calculate linear mappings
               IF(CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES>0) THEN                  
                 CALL CONSTRAINT_MAPPING_LINEAR_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*999)
                 LINEAR_MAPPING=>CONSTRAINT_MAPPING%LINEAR_MAPPING
                 IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                   LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
                   ALLOCATE(LINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Lagrange dof to column map.",ERR,ERROR,*999)
                   !1-1 mapping for now
                   DO dof_idx=1,LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                     column_idx=LAGRANGE_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                     LINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(dof_idx)=column_idx
                   ENDDO
                   ALLOCATE(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP,STAT=ERR)
                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices map.",ERR,ERROR,*999)
                   CALL CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT( &
                     & LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*999)
                    LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_TYPE=CREATE_VALUES_CACHE%LINEAR_VARIABLE_TYPE
                    CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,CREATE_VALUES_CACHE%LINEAR_VARIABLE_TYPE,FIELD_VARIABLE, &
                      & ERR,ERROR,*999)
                    LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE=>FIELD_VARIABLE
                    LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=CREATE_VALUES_CACHE% &
                      & NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
                    IF(CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE/=0) &
                      & LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=-1
                    !Allocate and initialise the variable to constraint matrices maps
                    FIELD_VARIABLE=>LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      IF(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES==-1) THEN
                        !!TODO: check if this can be removed and just allocate those variables that are actually used
                        ALLOCATE(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP=0
                      ELSE IF(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES>0) THEN
                        ALLOCATE(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS( &
                          & LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) &
                          & CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps constraint matrix numbers.", &
                          & ERR,ERROR,*999)
                        LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS=0
 
                        ALLOCATE(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                        !1-1 mappings for now.
                        row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                        LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP(row_idx)=dof_idx
                        ENDDO !dof_idx
                      ENDIF
                    ENDIF
                    !Allocate and initialise the constraint matrix to variable maps types
                    ALLOCATE(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(LINEAR_MAPPING% &
                      & NUMBER_OF_LINEAR_CONSTRAINT_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrix rows to variable maps.", &
                      & ERR,ERROR,*999)
                    !Loop over the constraint matrices and calculate the row mappings
                    !The pointers below have been checked for association above.
                    matrix_idx=1
                    !Initialise and setup the constraint matrix
                    CALL CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP( &
                      & matrix_idx),ERR,ERROR,*999)
                    EQUATIONS_SET=>CONSTRAINT_DEPENDENT%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,CREATE_VALUES_CACHE%LINEAR_VARIABLE_TYPE,FIELD_VARIABLE, &
                        & ERR,ERROR,*999)
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%EQUATIONS_SET=>EQUATIONS_SET
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%MATRIX_COEFFICIENT= &
                          & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(matrix_idx)
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                          & CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE(matrix_idx)
                         !Set the number of rows
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE%NUMBER_OF_DOFS
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                          & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                        !Set the row mapping
                        LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_DOFS_MAPPING=> &
                          & FIELD_VARIABLE%DOMAIN_MAPPING
                        ALLOCATE(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                        !1-1 mapping for now
                        DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP(row_idx)=dof_idx
                        ENDDO !dof_idx
                      ELSE
                        CALL FLAG_ERROR("Dependent variable could not be found.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set could not be found.",ERR,ERROR,*999)
                    ENDIF
                    !The pointers below have been checked for association above.
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      !Sets up the Lagrange-(Penalty) constraint matrix mapping and calculate the row mappings
                      matrix_idx=LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES !last of the constraint matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP( &
                        & matrix_idx),ERR,ERROR,*999)
                      CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                        & ERR,ERROR,*999)
                      FIELD_VARIABLE=>LAGRANGE_VARIABLE
                      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
                      IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%CONSTRAINT_EQUATIONS=>CONSTRAINT_EQUATIONS
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%MATRIX_COEFFICIENT=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(matrix_idx)
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE(matrix_idx)
                            !Set the number of rows
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE%NUMBER_OF_DOFS
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                            LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)% &
                              & ROW_TO_VARIABLE_DOF_MAP(row_idx)=row_idx
                          ENDDO !dof_idx
                        ELSE
                          CALL FLAG_ERROR("Lagrange variable could not be found.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Constraint equations could not be found.",ERR,ERROR,*999)
                      ENDIF
                    ENDSELECT
                 ELSE
                   CALL FLAG_ERROR("Linear mapping is not associated.",ERR,ERROR,*999)
                 ENDIF
               ENDIF
 
                !Calculate nonlinear mappings
                IF(CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES>0) THEN                  
                  CALL CONSTRAINT_MAPPING_NONLINEAR_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*999)
                  NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
                  IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                    NONLINEAR_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES= &
                      & CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES
                    ALLOCATE(NONLINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Lagrange dof to column map.",ERR,ERROR,*999)
                    !1-1 mapping for now
                    DO dof_idx=1,LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                      column_idx=LAGRANGE_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                      NONLINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP(dof_idx)=column_idx
                    ENDDO
                    ALLOCATE(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP,STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint Jacobians map.",ERR,ERROR,*999)
                    CALL CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_INIT( &
                      & NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP,ERR,ERROR,*999)
                    NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_TYPE=CREATE_VALUES_CACHE%JACOBIAN_VARIABLE_TYPE
                    CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,CREATE_VALUES_CACHE%JACOBIAN_VARIABLE_TYPE,FIELD_VARIABLE, &
                      & ERR,ERROR,*999)
                    NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE=>FIELD_VARIABLE
                    NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%NUMBER_OF_CONSTRAINT_JACOBIANS=CREATE_VALUES_CACHE% &
                      & NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES
                    IF(CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE/=0) &
                      & NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%NUMBER_OF_CONSTRAINT_JACOBIANS=-1
                    !Allocate and initialise the variable to constraint matrices maps
                    FIELD_VARIABLE=>NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      IF(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%NUMBER_OF_CONSTRAINT_JACOBIANS==-1) THEN
                        !!TODO: check if this can be removed and just allocate those variables that are actually used
                        ALLOCATE(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_DOF_TO_ROWS_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_DOF_TO_ROWS_MAP=0
                      ELSE IF(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%NUMBER_OF_CONSTRAINT_JACOBIANS>0) THEN
                        ALLOCATE(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%CONSTRAINT_JACOBIAN_NUMBERS( &
                          & NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%NUMBER_OF_CONSTRAINT_JACOBIANS),STAT=ERR)
                        IF(ERR/=0) &
                          & CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps constraint matrix numbers.", &
                          & ERR,ERROR,*999)
                        NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%CONSTRAINT_JACOBIAN_NUMBERS=0
 
                        ALLOCATE(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_DOF_TO_ROWS_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable to constraint matrices maps dof to rows map.", &
                          & ERR,ERROR,*999)
                        DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                        !1-1 mappings for now.
                        row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                        NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_DOF_TO_ROWS_MAP(dof_idx)=row_idx
                        ENDDO !dof_idx
                      ENDIF
                    ENDIF
                    !Allocate and initialise the Jacobian constraint matrix to variable maps types
                    ALLOCATE(NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(NONLINEAR_MAPPING% &
                      & NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint Jacobian rows to variable maps.", &
                      & ERR,ERROR,*999)
                    !Loop over the constraint matrices and calculate the row mappings
                    !The pointers below have been checked for association above.
                    matrix_idx=1
                    !Initialise and setup the constraint matrix
                    CALL CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT(NONLINEAR_MAPPING% &
                      & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx),ERR,ERROR,*999)
                    EQUATIONS_SET=>CONSTRAINT_DEPENDENT%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,CREATE_VALUES_CACHE%JACOBIAN_VARIABLE_TYPE,FIELD_VARIABLE, &
                        & ERR,ERROR,*999)
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%EQUATIONS_SET=>EQUATIONS_SET
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE= &
                          & FIELD_VARIABLE%VARIABLE_TYPE
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%JACOBIAN_COEFFICIENT= &
                          & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%JACOBIAN_COEFFICIENT
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                          & CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE(matrix_idx)
                         !Set the number of rows
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS= &
                          & FIELD_VARIABLE%NUMBER_OF_DOFS
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                          & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                        !Set the row mapping
                        NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%ROW_DOFS_MAPPING=> &
                          & FIELD_VARIABLE%DOMAIN_MAPPING
                        ALLOCATE(NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP( &
                          & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                        !1-1 mapping for now
                        DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP(row_idx)=dof_idx
                        ENDDO !dof_idx
                      ELSE
                        CALL FLAG_ERROR("Dependent variable could not be found.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set could not be found.",ERR,ERROR,*999)
                    ENDIF
                    !The pointers below have been checked for association above.
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      !Sets up the Lagrange-(Penalty) constraint matrix mapping and calculate the row mappings
                      matrix_idx=NONLINEAR_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES !last of the constraint matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT(NONLINEAR_MAPPING% &
                        & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx),ERR,ERROR,*999)
                      NULLIFY(LAGRANGE_VARIABLE)
                      CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                        & ERR,ERROR,*999)
                      NULLIFY(CONSTRAINT_EQUATIONS)
                      NULLIFY(FIELD_VARIABLE)
                      FIELD_VARIABLE=>LAGRANGE_VARIABLE
                      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
                      IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%CONSTRAINT_EQUATIONS=> &
                            & CONSTRAINT_EQUATIONS
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE% &
                            & VARIABLE_TYPE
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%JACOBIAN_COEFFICIENT= &
                            & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%JACOBIAN_COEFFICIENT
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE(matrix_idx)
                            !Set the number of rows
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE% &
                            & NUMBER_OF_DOFS
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            row_idx=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                            NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)% &
                              & ROW_TO_VARIABLE_DOF_MAP(row_idx)=dof_idx
                          ENDDO !dof_idx
                        ELSE
                          CALL FLAG_ERROR("Lagrange variable could not be found.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Constraint equations could not be found.",ERR,ERROR,*999)
                      ENDIF
                    ENDSELECT
                  ELSE
                    CALL FLAG_ERROR("Nonlinear mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
 
               !Calculate RHS mappings
               IF(CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE/=0) THEN
                 CALL CONSTRAINT_MAPPING_RHS_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*999)
                 RHS_MAPPING=>CONSTRAINT_MAPPING%RHS_MAPPING
                 IF(ASSOCIATED(RHS_MAPPING)) THEN
                   RHS_MAPPING%RHS_VARIABLE_TYPE=CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE
                   LAGRANGE_VARIABLE=>LAGRANGE_FIELD%VARIABLE_TYPE_MAP(CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE)%PTR
                   RHS_MAPPING%RHS_VARIABLE=>LAGRANGE_VARIABLE
                   RHS_MAPPING%RHS_VARIABLE_MAPPING=>LAGRANGE_VARIABLE%DOMAIN_MAPPING
                   RHS_MAPPING%RHS_COEFFICIENT=CREATE_VALUES_CACHE%RHS_COEFFICIENT
                   !Allocate and set up the row mappings
                   ALLOCATE(RHS_MAPPING%RHS_DOF_TO_CONSTRAINT_ROW_MAP(LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate rhs dof to constraint row map.",ERR,ERROR,*999)
                   ALLOCATE(RHS_MAPPING%CONSTRAINT_ROW_TO_RHS_DOF_MAP(CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS),STAT=ERR)
                   IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint row to dof map.",ERR,ERROR,*999)
                   DO dof_idx=1,LAGRANGE_VARIABLE%TOTAL_NUMBER_OF_DOFS
                     !1-1 mapping for now  
                     row_idx=LAGRANGE_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(dof_idx)
                     RHS_MAPPING%RHS_DOF_TO_CONSTRAINT_ROW_MAP(dof_idx)=row_idx
                     RHS_MAPPING%CONSTRAINT_ROW_TO_RHS_DOF_MAP(row_idx)=dof_idx
                   ENDDO !dof_idx
                 ELSE
                   CALL FLAG_ERROR("RHS mapping is not associated.",ERR,ERROR,*999)
                 ENDIF
               ENDIF
             ELSE
               CALL FLAG_ERROR("Constraint condition dependent is not associated.",ERR,ERROR,*999)
             ENDIF
           ELSE
             CALL FLAG_ERROR("Constraint condition Lagrange is not associated.",ERR,ERROR,*999)
           ENDIF
         CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE DEFAULT
           LOCAL_ERROR="The constraint condition method of "// &
             & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
           CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
         END SELECT
       ELSE
         CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",ERR,ERROR,*999)
       ENDIF
     ELSE
       CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
     ENDIF
   ELSE
     CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*999)
   ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Constraint mappings:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",CONSTRAINT_MAPPING%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of columns = ",CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global columns = ",CONSTRAINT_MAPPING%NUMBER_OF_GLOBAL_COLUMNS, &
        & ERR,ERROR,*999)
      DYNAMIC_MAPPING=>CONSTRAINT_MAPPING%DYNAMIC_MAPPING
      IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Dynamic mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dynamic equations matrices = ",DYNAMIC_MAPPING% &
          & NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic stiffness matrix number = ",DYNAMIC_MAPPING% &
          & STIFFNESS_MATRIX_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic damping matrix number = ",DYNAMIC_MAPPING% &
          & DAMPING_MATRIX_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic mass matrix number = ",DYNAMIC_MAPPING% &
          & MASS_MATRIX_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dynamic variable type = ",DYNAMIC_MAPPING% &
          & VARIABLE_TYPE,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ", DYNAMIC_MAPPING%VARIABLE_TYPE,ERR,ERROR,*999)
        IF(ASSOCIATED(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",DYNAMIC_MAPPING% &
            & VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",DYNAMIC_MAPPING% &
            & VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES,ERR,ERROR,*999)
          IF(DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & NUMBER_OF_CONSTRAINT_MATRICES,4,4,DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & CONSTRAINT_MATRIX_NUMBERS,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',ERR,ERROR,*999) 
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Variable DOF to row maps :",ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & VARIABLE_DOF_TO_ROWS_MAP, &
              & '("        Row numbers :",5(X,I13))','(24X,5(X,I13))',ERR,ERROR,*999) 
          ENDIF
        ENDIF
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",ERR,ERROR,*999)
        DO matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrix_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",DYNAMIC_MAPPING% &
            & CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows = ",DYNAMIC_MAPPING% &
            & CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",DYNAMIC_MAPPING% &
            & CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%MATRIX_COEFFICIENT,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)% &
            & NUMBER_OF_ROWS,5,5,DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP, &
            & '("        Variable row to DOF maps :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999) 
        ENDDO !matrix_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Column mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS,5,5, &
          & DYNAMIC_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP,'("    DOF to column map :",5(X,I13))','(21X,5(X,I13))', &
          & ERR,ERROR,*999) 
      ENDIF
      LINEAR_MAPPING=>CONSTRAINT_MAPPING%LINEAR_MAPPING
      IF(ASSOCIATED(LINEAR_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Linear mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of linear constraint matrices = ",LINEAR_MAPPING% &
          & NUMBER_OF_LINEAR_CONSTRAINT_MATRICES,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ",LINEAR_MAPPING%VARIABLE_TYPE,ERR,ERROR,*999)
        IF(ASSOCIATED(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",LINEAR_MAPPING% &
            & VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",LINEAR_MAPPING% &
            & VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES,ERR,ERROR,*999)
          IF(LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & NUMBER_OF_CONSTRAINT_MATRICES,4,4,LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & CONSTRAINT_MATRIX_NUMBERS,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',ERR,ERROR,*999) 
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Variable DOF to row maps :",ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP% &
              & VARIABLE_DOF_TO_ROWS_MAP, &
              & '("        Row numbers :",5(X,I13))','(24X,5(X,I13))',ERR,ERROR,*999) 
          ENDIF
        ENDIF
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",ERR,ERROR,*999)
        DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrix_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",LINEAR_MAPPING% &
            & CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows = ",LINEAR_MAPPING% &
            & CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",LINEAR_MAPPING% &
            & CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%MATRIX_COEFFICIENT,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)% &
            & NUMBER_OF_ROWS,5,5,LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP, &
            & '("        Variable row to DOF maps :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999) 
        ENDDO !matrix_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Column mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS,5,5, &
          & LINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP,'("    DOF to column map :",5(X,I13))','(21X,5(X,I13))', &
          & ERR,ERROR,*999) 
      ENDIF
      NONLINEAR_MAPPING=>CONSTRAINT_MAPPING%NONLINEAR_MAPPING
      IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Nonlinear mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of nonlinear constraint matrices = ",NONLINEAR_MAPPING% &
          & NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Variable to matrices mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type : ",NONLINEAR_MAPPING%VARIABLE_TYPE,ERR,ERROR,*999)
        IF(ASSOCIATED(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Total number of DOFs = ",NONLINEAR_MAPPING% &
            & VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE%TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of equations matrices = ",NONLINEAR_MAPPING% &
            & VAR_TO_CONSTRAINT_JACOBIAN_MAP%NUMBER_OF_CONSTRAINT_JACOBIANS,ERR,ERROR,*999)
          IF(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP%NUMBER_OF_CONSTRAINT_JACOBIANS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP% &
              & NUMBER_OF_CONSTRAINT_JACOBIANS,4,4,NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP% &
              & CONSTRAINT_JACOBIAN_NUMBERS,'("      Matrix numbers :",4(X,I12))','(22X,4(X,I12))',ERR,ERROR,*999) 
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Variable DOF to row maps :",ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP% &
              & VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5,NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP% &
              & VARIABLE_DOF_TO_ROWS_MAP, &
              & '("        Row numbers :",5(X,I13))','(24X,5(X,I13))',ERR,ERROR,*999) 
          ENDIF
        ENDIF
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Matrix to variable mappings:",ERR,ERROR,*999)
        DO matrix_idx=1,NONLINEAR_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Matrix number : ",matrix_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Variable type = ",NONLINEAR_MAPPING% &
            & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%VARIABLE_TYPE,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows = ",NONLINEAR_MAPPING% &
            & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%NUMBER_OF_ROWS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Matrix coefficient = ",NONLINEAR_MAPPING% &
            & CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%JACOBIAN_COEFFICIENT,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)% &
            & NUMBER_OF_ROWS,5,5,NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx)%ROW_TO_VARIABLE_DOF_MAP, &
            & '("        Variable row to DOF maps :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999) 
        ENDDO !matrix_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Column mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS,5,5, &
          & NONLINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP,'("    DOF to column map :",5(X,I13))','(21X,5(X,I13))', &
          & ERR,ERROR,*999) 
      ENDIF
      RHS_MAPPING=>CONSTRAINT_MAPPING%RHS_MAPPING
      IF(ASSOCIATED(RHS_MAPPING)) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  RHS mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    RHS variable type = ", &
          & RHS_MAPPING%RHS_VARIABLE_TYPE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of RHS DOFs = ",RHS_MAPPING%RHS_VARIABLE% &
          & TOTAL_NUMBER_OF_DOFS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    RHS coefficient = ",RHS_MAPPING%RHS_COEFFICIENT,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Row mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,RHS_MAPPING%RHS_VARIABLE%TOTAL_NUMBER_OF_DOFS,5,5, &
          & RHS_MAPPING%RHS_DOF_TO_CONSTRAINT_ROW_MAP,'("    DOF to row mappings :",5(X,I13))','(25X,5(X,I13))',ERR,ERROR,*999) 
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS,5,5, &
          & RHS_MAPPING%CONSTRAINT_ROW_TO_RHS_DOF_MAP,'("    Row to DOF mappings :",5(X,I13))','(25X,5(X,I13))',ERR,ERROR,*999) 
      ENDIF

    ENDIF

    CALL EXITS("CONSTRAINT_MAPPING_CALCULATE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_CALCULATE")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of constraint mapping
  SUBROUTINE CONSTRAINT_MAPPING_CREATE_FINISH(CONSTRAINT_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to finish the creation of.
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_MAPPING_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        !Calculate the equations mapping and clean up
        CALL CONSTRAINT_MAPPING_CALCULATE(CONSTRAINT_MAPPING,ERR,ERROR,*999)
        CALL CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
        CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED=.TRUE.
      ENDIF
   ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !
  !>Starts the process of creating the constraint mapping for constraint equations.
  SUBROUTINE CONSTRAINT_MAPPING_CREATE_START(CONSTRAINT_EQUATIONS,CONSTRAINT_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS !<A pointer to the constraint equations to create the mapping for.
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<On exit, a pointer to the created constraint mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINT_MAPPING_CREATE_START",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
      IF(CONSTRAINT_EQUATIONS%CONSTRAINT_EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
          CALL FLAG_ERROR("Constraint mapping is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(CONSTRAINT_MAPPING)
          CALL CONSTRAINT_MAPPING_INITIALISE(CONSTRAINT_EQUATIONS,ERR,ERROR,*999)
          CONSTRAINT_MAPPING=>CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint equations have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CREATE_START",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_START")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_CREATE_START
 
  !
  !================================================================================================================================
  !
 
  !>Finalises a constraint mapping create values cache and deallocates all memory
  SUBROUTINE CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE)) DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE)) DEALLOCATE(CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE)) DEALLOCATE(CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE
 
 !
 !================================================================================================================================
 !
  !>Initialises a constraint mapping create values cache
  SUBROUTINE CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to create the create values cache for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_type_idx,variable_type_idx2
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD,CONSTRAINT_DEPENDENT_FIELD
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_DEPENDENT_TYPE), POINTER :: CONSTRAINT_DEPENDENT
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(ASSOCIATED(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Constraint mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
          CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
          IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
            !Allocate and initialise the create values cache
            ALLOCATE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint mapping create values cache.",ERR,ERROR,*999)
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=0
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=0
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=0
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=0
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%JACOBIAN_COEFFICIENT=1.0_DP
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=0
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=0
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=1.0_DP
            !First calculate how many constraint matrices we have and set the variable types
            SELECT CASE(CONSTRAINT_CONDITION%METHOD)
            CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
              LAGRANGE=>CONSTRAINT_CONDITION%LAGRANGE
              IF(ASSOCIATED(LAGRANGE)) THEN
                LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                  CONSTRAINT_DEPENDENT=>CONSTRAINT_CONDITION%DEPENDENT
                  IF(ASSOCIATED(CONSTRAINT_DEPENDENT)) THEN
                    CONSTRAINT_DEPENDENT_FIELD=>CONSTRAINT_DEPENDENT%DEPENDENT_FIELD
                    IF(ASSOCIATED(CONSTRAINT_DEPENDENT_FIELD)) THEN
                      !The pointers below have been checked for association above.
                      !Default the Lagrange variable to the first Lagrange variable
                      CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=0
                      DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                        IF(ASSOCIATED(LAGRANGE_FIELD%VARIABLE_TYPE_MAP(variable_type_idx)%PTR)) THEN
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=variable_type_idx
                          EXIT
                        ENDIF
                      ENDDO !variable_type_idx
                      IF(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE==0) &
                        & CALL FLAG_ERROR("Could not find a Lagrange variable type in the Lagrange field.",ERR,ERROR,*999)
                      !Default the RHS Lagrange variable to the second Lagrange variable
                      DO variable_type_idx2=variable_type_idx+1,FIELD_NUMBER_OF_VARIABLE_TYPES
                        IF(ASSOCIATED(LAGRANGE_FIELD%VARIABLE_TYPE_MAP(variable_type_idx2)%PTR)) THEN
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=variable_type_idx2
                          EXIT
                        ENDIF
                      ENDDO !variable_type_idx2
                      IF(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE==0) &
                        & CALL FLAG_ERROR("Could not find a RHS Lagrange variable type in the Lagrange field.",ERR,ERROR,*999)
                      ALLOCATE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(CONSTRAINT_MAPPING% &
                        & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix coefficients.",ERR,ERROR,*999)
                      !Set the default constraint mapping in the create values cache
                      SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
                      CASE(CONSTRAINT_CONDITION_STATIC)
                        SELECT CASE(CONSTRAINT_EQUATIONS%LINEARITY)
                        CASE(CONSTRAINT_CONDITION_LINEAR,CONSTRAINT_CONDITION_NONLINEAR_BCS)
                          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=1
                          CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables plus a single
                            !Lagrange variable
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=2
                          END SELECT
                        CASE(CONSTRAINT_CONDITION_NONLINEAR)
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=0
                          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                            !Default the number of Jacobia constraint matrices to the number of added dependent variables
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES=1
                          CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                            !Default the number of Jacobian constraint matrices to the number of added dependent variables plus a
                            !single Lagrange variable
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES=2
                          END SELECT
                        CASE DEFAULT
                          LOCAL_ERROR="The constraint linearity type of "// &
                            & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC,CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
                        IF(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE==CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC) THEN
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=2
                          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=2
                          CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables plus a single Lagrange
                            !variable
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=4
                          END SELECT
                        ELSE
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=1
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=2
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=3
                          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=3
                          CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables plus a single Lagrange variable
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=6
                          END SELECT
                        ENDIF
                        CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=0
                      CASE DEFAULT
                        LOCAL_ERROR="The constraint time dependence type of "// &
                          & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                      !Allocate the dynamic matrix coefficients and set their values
                      IF(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES>0) THEN
                        ALLOCATE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CONSTRAINT_MAPPING% &
                          & CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL &
                          & FLAG_ERROR("Could not allocate constraint mapping create values cache dynamic matrix coefficients.", &
                          & ERR,ERROR,*999)
                        CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS=1.0_DP !Equations matrices are added by default
                        ALLOCATE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE(CONSTRAINT_MAPPING% &
                          & CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache dynamic has transpose.",ERR,ERROR,*999)
                        !Default the dynamic constraint matrices to all have a transpose
                        CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE=.TRUE.
                        !The pointers below have been checked for association above.
                        SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                        CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                          !Default the constraint matrix (Penalty) to have no transpose
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE(CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES)=.FALSE.
                        END SELECT
                      ENDIF
                      !Allocate the Jacobian matrix variable types and Jacobian matrix coefficients and set their values
                      IF(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES>0) THEN
                        ALLOCATE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE(CONSTRAINT_MAPPING% &
                          & CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache Jacobian has transpose.",ERR,ERROR,*999)
                        !Default the Jacobia constraint matrices to all have a transpose
                        CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE=.TRUE.
                        !The pointers below have been checked for association above.
                        SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                        CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                          !Default the constraint matrix (Penalty) to have no transpose
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE(CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES)=.FALSE.
                        END SELECT
                      ENDIF
                      !Allocate the linear matrix variable types and linear matrix coefficients and set their values
                      IF(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES>0) THEN
                        ALLOCATE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(CONSTRAINT_MAPPING% &
                          & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL &
                          & FLAG_ERROR("Could not allocate constraint mapping create values cache linear matrix coefficients.", &
                          & ERR,ERROR,*999)
                        !Default the linear constraint matrices coefficients to add.
                        CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS=1.0_DP
                        ALLOCATE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE(CONSTRAINT_MAPPING% &
                          & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache linear has transpose.",ERR,ERROR,*999)
                        !Default the linear constraint matrices to all have a transpose
                        CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE=.TRUE.
                        !The pointers below have been checked for association above.
                        SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                        CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                          !Default the constraint matrix (Penalty) to have no transpose
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE(CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES)=.FALSE.
                        END SELECT
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Constraint condition dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Constraint condition dependent is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Constraint condition Lagrange field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Constraint condition Lagrange is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The constraint equations method of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping constraint equations is not associated.",ERR,ERROR,*998)
        ENDIF
        ENDIF
      ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
    ENDIF
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_INITIALISE
 
  !
  !================================================================================================================================
  !
 
  !>Destroys a constraint mapping.
  SUBROUTINE CONSTRAINT_MAPPING_DESTROY(CONSTRAINT_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer the constraint mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_DESTROY",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      CALL CONSTRAINT_MAPPING_FINALISE(CONSTRAINT_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("CONSTRAINT_MAPPING_DESTROY")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_DESTROY",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_DESTROY
 
 
 !
 !================================================================================================================================
 !
 
  !>Finalise a constraint matrix to variable maps and deallocate all memory.
  SUBROUTINE CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE(CONSTRAINT_MATRIX_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRIX_TO_VAR_MAP_TYPE) :: CONSTRAINT_MATRIX_TO_VAR_MAP !<The constraint matrix to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR,*999)
 
    IF(ALLOCATED(CONSTRAINT_MATRIX_TO_VAR_MAP%ROW_TO_VARIABLE_DOF_MAP)) &
      & DEALLOCATE(CONSTRAINT_MATRIX_TO_VAR_MAP%ROW_TO_VARIABLE_DOF_MAP)
    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE
 
  !
  !================================================================================================================================
  !
 
  !>Initialise an constraint matrix to variable maps.
  SUBROUTINE CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT(CONSTRAINT_MATRIX_TO_VAR_MAP,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MATRIX_TO_VAR_MAP_TYPE) :: CONSTRAINT_MATRIX_TO_VAR_MAP !<The constraint matrix to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT",ERR,ERROR,*999)
 
    CONSTRAINT_MATRIX_TO_VAR_MAP%MATRIX_NUMBER=0
    CONSTRAINT_MATRIX_TO_VAR_MAP%VARIABLE_TYPE=0
    NULLIFY(CONSTRAINT_MATRIX_TO_VAR_MAP%VARIABLE)
    CONSTRAINT_MATRIX_TO_VAR_MAP%NUMBER_OF_ROWS=0
    CONSTRAINT_MATRIX_TO_VAR_MAP%MATRIX_COEFFICIENT=1.0_DP !Matrices in constraint condition are added by default
    NULLIFY(CONSTRAINT_MATRIX_TO_VAR_MAP%ROW_DOFS_MAPPING)
    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_INIT
 
  !
  !================================================================================================================================
  !
 
  !>Finalises the constraint mapping dynamic  mapping and deallocates all memory
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE(DYNAMIC_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING !<A pointer to the dynamic mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE",ERR,ERROR,*999)
 
    IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
      CALL CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE(DYNAMIC_MAPPING% &
        & VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*999)
      IF(ALLOCATED(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP)) THEN
        DO matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
          CALL CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE( &
            & DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP)
      ENDIF
      IF(ALLOCATED(DYNAMIC_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)) &
        & DEALLOCATE(DYNAMIC_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)
      DEALLOCATE(DYNAMIC_MAPPING)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE
 
  !
  !================================================================================================================================
  !
 
  !>Initialises the constraint mapping dynamic mapping
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to initialise the dynamic mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(ASSOCIATED(CONSTRAINT_MAPPING%DYNAMIC_MAPPING)) THEN
        CALL FLAG_ERROR("Constraint mapping dynamic mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONSTRAINT_MAPPING%DYNAMIC_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint mapping dynamic mapping.",ERR,ERROR,*999)
        CONSTRAINT_MAPPING%DYNAMIC_MAPPING%CONSTRAINT_MAPPING=>CONSTRAINT_MAPPING
        CONSTRAINT_MAPPING%DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=0
        CONSTRAINT_MAPPING%DYNAMIC_MAPPING%STIFFNESS_MATRIX_NUMBER=0
        CONSTRAINT_MAPPING%DYNAMIC_MAPPING%DAMPING_MATRIX_NUMBER=0
        CONSTRAINT_MAPPING%DYNAMIC_MAPPING%MASS_MATRIX_NUMBER=0
        CONSTRAINT_MAPPING%DYNAMIC_MAPPING%VARIABLE_TYPE=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE(CONSTRAINT_MAPPING%DYNAMIC_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MAPPING_INITIALISE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MAPPING_INITIALISE
 
  !
  !================================================================================================================================
  !
 
  !>Sets/changes the matrices involved in dynamic constraint mapping
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL(CONSTRAINT_MAPPING,MASS_MATRIX,DAMPING_MATRIX,STIFFNESS_MATRIX, &
    & ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set the atrices for
    LOGICAL, INTENT(IN) :: MASS_MATRIX !<Is .TRUE. if the mass matrix is in the constraint mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: DAMPING_MATRIX !<Is .TRUE. if the damping matrix is in the constraint mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: STIFFNESS_MATRIX !<Is .TRUE. if the stiffness matrix is in the constraint mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_DYNAMIC_DAMPING_MATRIX_NUMBER,NEW_DYNAMIC_MASS_MATRIX_NUMBER,NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER, &
      & NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
    REAL(DP), ALLOCATABLE :: OLD_DYNAMIC_MATRIX_COEFFICIENTS(:)
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
      IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
        CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
        IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
          CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
          IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
            SELECT CASE(CONSTRAINT_EQUATIONS%LINEARITY)
            CASE(CONSTRAINT_CONDITION_LINEAR,CONSTRAINT_CONDITION_NONLINEAR_BCS)
              NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=0
              NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=0
              NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=0
              NEW_DYNAMIC_MASS_MATRIX_NUMBER=0
              IF(STIFFNESS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES+1
                NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
              ENDIF
              IF(DAMPING_MATRIX) THEN
                NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES+1
                NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
              ENDIF
              IF(MASS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES+1
                NEW_DYNAMIC_MASS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
              ENDIF
              IF(NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES>0) THEN
                ALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old dynamic matrix coefficients.",ERR,ERROR,*999)
                OLD_DYNAMIC_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES)= &
                  & CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES)
                DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)
                ALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic matrix coefficients.",ERR,ERROR,*999)
                IF(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=NEW_DYNAMIC_DAMPING_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=NEW_DYNAMIC_MASS_MATRIX_NUMBER
                IF(ALLOCATED(OLD_DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS)
              ELSE
                CALL FLAG_ERROR("Invalid dynamic matrices set up. There are no dynamic constraint matrices.",ERR,ERROR,*999)
              ENDIF
            CASE(CONSTRAINT_CONDITION_NONLINEAR)
              NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=0
              NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=0
              NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=0
              NEW_DYNAMIC_MASS_MATRIX_NUMBER=0
              IF(STIFFNESS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES+1
                NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
              ENDIF
              IF(DAMPING_MATRIX) THEN
                NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES+1
                NEW_DYNAMIC_DAMPING_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
              ENDIF
              IF(MASS_MATRIX) THEN
                NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES+1
                NEW_DYNAMIC_MASS_MATRIX_NUMBER=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
              ENDIF
              IF(NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES>0) THEN
                ALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old dynamic matrix coefficients.",ERR,ERROR,*999)
                OLD_DYNAMIC_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES)= &
                  & CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES)
                DEALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS)
                ALLOCATE(CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dynamic matrix coefficients.",ERR,ERROR,*999)
                IF(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                IF(NEW_DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                  IF(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER==0) THEN
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)=1.0_DP
                  ELSE
                    CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(NEW_DYNAMIC_MASS_MATRIX_NUMBER)= &
                      & OLD_DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER)
                  ENDIF
                ENDIF
                CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER=NEW_DYNAMIC_STIFFNESS_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=NEW_DYNAMIC_DAMPING_MATRIX_NUMBER
                CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER=NEW_DYNAMIC_MASS_MATRIX_NUMBER
                IF(ALLOCATED(OLD_DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS)
              ELSE
                CALL FLAG_ERROR("Invalid dynamic matrices set up. There are no dynamic constraint matrices.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The equations linearity type of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations constraint condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint mapping equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL_ORDER")
    RETURN
999 IF(ALLOCATED(OLD_DYNAMIC_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_DYNAMIC_MATRIX_COEFFICIENTS)    
    CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL
 
  !
  !================================================================================================================================
  !
 
  !>Sets/changes the matrices involved in a first order dynamic constraint mapping
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_1(CONSTRAINT_MAPPING,DAMPING_MATRIX,STIFFNESS_MATRIX,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set the first order matrices for
    LOGICAL, INTENT(IN) :: DAMPING_MATRIX !<Is .TRUE. if the damping matrix is in the constraint mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: STIFFNESS_MATRIX !<Is .TRUE. if the stiffness matrix is in the constraint mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_1",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
          SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
          CASE(CONSTRAINT_CONDITION_STATIC)
            CALL FLAG_ERROR("Can not set dynamic matrices for static equations.",ERR,ERROR,*999)
          CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC)
            IF(.NOT.DAMPING_MATRIX) CALL FLAG_WARNING("No damping matrix for first order dynamic equations.",ERR,ERROR,*999)
            CALL CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL(CONSTRAINT_MAPPING,.FALSE.,DAMPING_MATRIX,STIFFNESS_MATRIX, &
              ERR,ERROR,*999)
          CASE(CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
            CALL FLAG_ERROR("Need to specify three matrices to set for second order dynamic equations.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations time dependence type of "// &
              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Constraint mapping equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_1")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_1",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_1")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_1
 
  !
  !================================================================================================================================
  !
 
  !>Sets/changes the matrices involved in a second order dynamic constraint mapping
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_2(CONSTRAINT_MAPPING,MASS_MATRIX,DAMPING_MATRIX,STIFFNESS_MATRIX, &
    & ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set the first order matrices for
    LOGICAL, INTENT(IN) :: MASS_MATRIX !<Is .TRUE. if the mass matrix is in the constraint mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: DAMPING_MATRIX !<Is .TRUE. if the damping matrix is in the constraint mapping, .FALSE. if not
    LOGICAL, INTENT(IN) :: STIFFNESS_MATRIX !<Is .TRUE. if the stiffness matrix is in the constraint mapping, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_2",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
        IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
          SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
          CASE(CONSTRAINT_CONDITION_STATIC)
            CALL FLAG_ERROR("Can not set dynamic matrices for static equations.",ERR,ERROR,*999)
          CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC)
            IF(MASS_MATRIX) THEN
              CALL FLAG_ERROR("The mass matrix cannot be present for first order dynamic equations.",ERR,ERROR,*999)
            ELSE
              IF(.NOT.DAMPING_MATRIX) CALL FLAG_WARNING("No damping matrix for a first order dynamic system.",ERR,ERROR,*999)
              CALL CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL(CONSTRAINT_MAPPING,.FALSE.,DAMPING_MATRIX,STIFFNESS_MATRIX, &
                ERR,ERROR,*999)
            ENDIF
          CASE(CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
            IF(.NOT.MASS_MATRIX) CALL FLAG_WARNING("No mass matrix for a second order dynamic system.",ERR,ERROR,*999)
            CALL CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_ALL(CONSTRAINT_MAPPING,MASS_MATRIX,DAMPING_MATRIX, &
              & STIFFNESS_MATRIX,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations time dependence type of "// &
              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Constraint mapping equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_2")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_2",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_2")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_SET_2
 
  !
  !================================================================================================================================
  !
 
  !>Sets/changes the matrix coefficients in a first order dynamic constraint mapping
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1(CONSTRAINT_MAPPING,DAMPING_MATRIX_COEFFICIENT, &
    & STIFFNESS_MATRIX_COEFFICIENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set 
    REAL(DP), INTENT(IN) :: DAMPING_MATRIX_COEFFICIENT !<The damping matrix coefficient
    REAL(DP), INTENT(IN) :: STIFFNESS_MATRIX_COEFFICIENT !<The stiffness matrix coefficient
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
          IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
            SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
            CASE(CONSTRAINT_CONDITION_STATIC)
              CALL FLAG_ERROR("Can not set dynamic matrix coefficients for static equations.",ERR,ERROR,*999)
            CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC)
              IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                  & STIFFNESS_MATRIX_COEFFICIENT
              ENDIF
              IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                  & DAMPING_MATRIX_COEFFICIENT
              ENDIF
            CASE(CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
              CALL FLAG_ERROR("Need to specify three matrix coefficients for second order dynamic equations.", &
                & ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Constraint mapping equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_1
 
  !
  !================================================================================================================================
  !
 
  !>Sets/changes the matrix coefficients in a second order dynamic constraint mapping
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2(CONSTRAINT_MAPPING,MASS_MATRIX_COEFFICIENT, &
    & DAMPING_MATRIX_COEFFICIENT,STIFFNESS_MATRIX_COEFFICIENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set 
    REAL(DP), INTENT(IN) :: MASS_MATRIX_COEFFICIENT !<The mass matrix coefficient
    REAL(DP), INTENT(IN) :: DAMPING_MATRIX_COEFFICIENT !<The damping matrix coefficient
    REAL(DP), INTENT(IN) :: STIFFNESS_MATRIX_COEFFICIENT !<The stiffness matrix coefficient
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
          IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
            SELECT CASE(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE)
            CASE(CONSTRAINT_CONDITION_STATIC)
              CALL FLAG_ERROR("Can not set dynamic matrices for static equations.",ERR,ERROR,*999)
            CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC)
              CALL FLAG_ERROR("Need to specify two matrix coefficients for second order dynamic equations.",ERR,ERROR,*999)
            CASE(CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
              IF(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_STIFFNESS_MATRIX_NUMBER)= &
                  & STIFFNESS_MATRIX_COEFFICIENT
              ENDIF
              IF(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER)= &
                  & DAMPING_MATRIX_COEFFICIENT
              ENDIF
              IF(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER/=0) THEN
                CREATE_VALUES_CACHE%DYNAMIC_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%DYNAMIC_MASS_MATRIX_NUMBER)= &
                  & MASS_MATRIX_COEFFICIENT
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The equations time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Constraint mapping equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_COEFFS_SET_2
 
  !
  !================================================================================================================================
  !
 
 
  !>Finalises a variable to constraint Jacobian map and deallocates all memory.
  SUBROUTINE CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_FINALISE(CONSTRAINT_JACOBIAN_TO_VAR_MAP,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_JACOBIAN_TO_VAR_MAP_TYPE) :: CONSTRAINT_JACOBIAN_TO_VAR_MAP !<The constraint Jacobian to variable map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(CONSTRAINT_JACOBIAN_TO_VAR_MAP%ROW_TO_VARIABLE_DOF_MAP)) &
      & DEALLOCATE(CONSTRAINT_JACOBIAN_TO_VAR_MAP%ROW_TO_VARIABLE_DOF_MAP)
    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_FINALISE
 
  !
  !================================================================================================================================
  !
 
  !>Initialises a variable to equations Jacobian map.
  SUBROUTINE CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT(CONSTRAINT_JACOBIAN_TO_VAR_MAP,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_JACOBIAN_TO_VAR_MAP_TYPE) :: CONSTRAINT_JACOBIAN_TO_VAR_MAP !<The constraint Jacobian to variable map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT",ERR,ERROR,*999)
    
    CONSTRAINT_JACOBIAN_TO_VAR_MAP%VARIABLE_TYPE=0
    NULLIFY(CONSTRAINT_JACOBIAN_TO_VAR_MAP%VARIABLE)
    NULLIFY(CONSTRAINT_JACOBIAN_TO_VAR_MAP%CONSTRAINT_JACOBIAN)
    CONSTRAINT_JACOBIAN_TO_VAR_MAP%NUMBER_OF_ROWS=0
    CONSTRAINT_JACOBIAN_TO_VAR_MAP%JACOBIAN_COEFFICIENT=0
    NULLIFY(CONSTRAINT_JACOBIAN_TO_VAR_MAP%ROW_DOFS_MAPPING)    
    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_INIT
 
  !
  !================================================================================================================================
  !
 
  !>Finalises the constraint mapping and deallocates all memory.
  SUBROUTINE CONSTRAINT_MAPPING_FINALISE(CONSTRAINT_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_FINALISE",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      NULLIFY(CONSTRAINT_MAPPING%COLUMN_DOFS_MAPPING)
      CALL CONSTRAINT_MAPPING_DYNAMIC_MAPPING_FINALISE(CONSTRAINT_MAPPING%DYNAMIC_MAPPING,ERR,ERROR,*999)
      CALL CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE(CONSTRAINT_MAPPING%LINEAR_MAPPING,ERR,ERROR,*999)
      CALL CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE(CONSTRAINT_MAPPING%NONLINEAR_MAPPING,ERR,ERROR,*999)
      CALL CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE(CONSTRAINT_MAPPING%RHS_MAPPING,ERR,ERROR,*999)
      CALL CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_FINALISE(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(CONSTRAINT_MAPPING)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_FINALISE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_FINALISE
  !
  !================================================================================================================================
  !

 !>Initialises the constraint mapping and deallocates all memory.
 SUBROUTINE CONSTRAINT_MAPPING_INITIALISE(CONSTRAINT_EQUATIONS,ERR,ERROR,*)

   !Argument variables
   TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS !<A pointer to the constraint equations to initialise the constraint mapping for
   INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
   TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   !Local Variables
   INTEGER(INTG) :: DUMMY_ERR
   TYPE(VARYING_STRING) :: DUMMY_ERROR

   CALL ENTERS("CONSTRAINT_MAPPING_INITIALISE",ERR,ERROR,*998)

   IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
     IF(ASSOCIATED(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING)) THEN
       CALL FLAG_ERROR("Constraint mapping is already associated.",ERR,ERROR,*998)
     ELSE
       ALLOCATE(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING,STAT=ERR)
       IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint equations constraint mapping.",ERR,ERROR,*999)
       CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS=>CONSTRAINT_EQUATIONS
       CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED=.FALSE.
       CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%LAGRANGE_VARIABLE_TYPE=0
       CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%NUMBER_OF_COLUMNS=0
       CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS=0
       CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%NUMBER_OF_GLOBAL_COLUMNS=0
       NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%LAGRANGE_VARIABLE)
       NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%COLUMN_DOFS_MAPPING)
       NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%DYNAMIC_MAPPING)
       NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%LINEAR_MAPPING)
       NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%NONLINEAR_MAPPING)
       NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%RHS_MAPPING)
       NULLIFY(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%CREATE_VALUES_CACHE)
       CALL CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_INITIALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING,ERR,ERROR,*999)
     ENDIF
   ELSE
     CALL FLAG_ERROR("Constraint equations is not associated.",ERR,ERROR,*998)
   ENDIF
      
   CALL EXITS("CONSTRAINT_MAPPING_INITIALISE")
   RETURN
999 CALL CONSTRAINT_MAPPING_FINALISE(CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MAPPING_INITIALISE",ERR,ERROR)
   CALL EXITS("CONSTRAINT_MAPPING_INITIALISE")
   RETURN 1
   
 END SUBROUTINE CONSTRAINT_MAPPING_INITIALISE

 !
 !================================================================================================================================
 !

 !>Finalises the constraint mapping linear  mapping and deallocates all memory
 SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE(LINEAR_MAPPING,ERR,ERROR,*)

   !Argument variables
   TYPE(CONSTRAINT_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING !<A pointer to the linear mapping to finalise
   INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
   TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   !Local Variables
   INTEGER(INTG) :: matrix_idx

   CALL ENTERS("CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE",ERR,ERROR,*999)

   IF(ASSOCIATED(LINEAR_MAPPING)) THEN
     CALL CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE(LINEAR_MAPPING% &
       & VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*999)
     IF(ALLOCATED(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP)) THEN
       DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
         CALL CONSTRAINT_MAPPING_CONSTR_MATRIX_TO_VAR_MAP_FINALISE(LINEAR_MAPPING% &
           & CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP(matrix_idx),ERR,ERROR,*999)
       ENDDO !matrix_idx
       DEALLOCATE(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAP)
     ENDIF
     IF(ALLOCATED(LINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)) &
       & DEALLOCATE(LINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)
     DEALLOCATE(LINEAR_MAPPING)
   ENDIF
      
   CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE")
   RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE",ERR,ERROR)
   CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE")
   RETURN 1
   
 END SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE

 !
 !================================================================================================================================
 !
 
  !>Initialises the constraint mapping linear mapping
  SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to initialise the linear mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_LINEAR_MAPPING_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(ASSOCIATED(CONSTRAINT_MAPPING%LINEAR_MAPPING)) THEN
        CALL FLAG_ERROR("Constraint mapping linear mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONSTRAINT_MAPPING%LINEAR_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint mapping linear mapping.",ERR,ERROR,*999)
        CONSTRAINT_MAPPING%LINEAR_MAPPING%CONSTRAINT_MAPPING=>CONSTRAINT_MAPPING       
        CONSTRAINT_MAPPING%LINEAR_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MAPPING_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MAPPING_LINEAR_MAPPING_FINALISE(CONSTRAINT_MAPPING%LINEAR_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MAPPING_LINEAR_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MAPPING_INITIALISE
 
  !
  !================================================================================================================================
  !
 
  !>Sets the coefficients for the linear constraint matrices in an equation set. 
  SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MATRICES_COEFFS_SET(CONSTRAINT_MAPPING,LINEAR_MATRIX_COEFFICIENTS,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    REAL(DP), INTENT(IN) :: LINEAR_MATRIX_COEFFICIENTS(:) !<The linear matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_LINEAR_MATRICES_COEFFS_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping is finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN          
          IF(SIZE(LINEAR_MATRIX_COEFFICIENTS,1)==CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES) THEN
            CREATE_VALUES_CACHE%LINEAR_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES)= &
              & LINEAR_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES)
          ELSE
            LOCAL_ERROR="Invalid size of linear matrix coefficeints. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(LINEAR_MATRIX_COEFFICIENTS,1),"*",ERR,ERROR))// &
              & ") must match the number of linear constraint matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE% &
              & NUMBER_OF_LINEAR_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MATRICES_COEFFS_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_LINEAR_MATRICES_COEFFS_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MATRICES_COEFFS_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MATRICES_COEFFS_SET
 
  !
  !================================================================================================================================
  !
 
  !>Finalises the constraint mapping nonlinear mapping and deallocates all memory
  SUBROUTINE CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE(NONLINEAR_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING !<A pointer to the nonlinear mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) matrix_idx
 
    CALL ENTERS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE",ERR,ERROR,*999)
 
    IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
      CALL CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_FINALISE(NONLINEAR_MAPPING%VAR_TO_CONSTRAINT_JACOBIAN_MAP,ERR,ERROR,*999)
      IF(ALLOCATED(NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP)) THEN
        DO matrix_idx=1,NONLINEAR_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES
          CALL CONSTRAINT_MAPPING_CONSTR_JACOBIAN_TO_VAR_MAP_FINALISE( &
            & NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP(matrix_idx),ERR,ERROR,*999)
        ENDDO
        DEALLOCATE(NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAP)
      ENDIF
      IF(ALLOCATED(NONLINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)) &
        & DEALLOCATE(NONLINEAR_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)
      DEALLOCATE(NONLINEAR_MAPPING)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE
 
  !
  !================================================================================================================================
  !
 
  !>Initialises the constraint mapping nonlinear mapping
  SUBROUTINE CONSTRAINT_MAPPING_NONLINEAR_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to initialise the nonlinear mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(ASSOCIATED(CONSTRAINT_MAPPING%NONLINEAR_MAPPING)) THEN
        CALL FLAG_ERROR("Constraint mapping nonlinear mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONSTRAINT_MAPPING%NONLINEAR_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint mapping nonlinear mapping.",ERR,ERROR,*999)
        CONSTRAINT_MAPPING%NONLINEAR_MAPPING%CONSTRAINT_MAPPING=>CONSTRAINT_MAPPING
        CONSTRAINT_MAPPING%NONLINEAR_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES=0
        CONSTRAINT_MAPPING%NONLINEAR_MAPPING%JACOBIAN_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
 
    CALL EXITS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MAPPING_NONLINEAR_MAPPING_FINALISE(CONSTRAINT_MAPPING%NONLINEAR_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_NONLINEAR_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_NONLINEAR_MAPPING_INITIALISE
 
  !
  !================================================================================================================================
  !
 
  !>Sets the coefficient applied to the constraint condition residual vector.
  SUBROUTINE CONSTRAINT_MAPPING_JACOBIAN_COEFF_SET(CONSTRAINT_MAPPING,JACOBIAN_COEFFICIENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set
    REAL(DP), INTENT(IN) :: JACOBIAN_COEFFICIENT!<The coefficient applied to the constraint condition residual vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_JACOBIAN_COEFF_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping have been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE)) THEN
          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%JACOBIAN_COEFFICIENT=JACOBIAN_COEFFICIENT
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_JACOBIAN_COEFF_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_JACOBIAN_COEFF_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_JACOBIAN_COEFF_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_JACOBIAN_COEFF_SET
 
  !
  !================================================================================================================================
  !
 
  !>Finalises a variable to constraint Jacobian map and deallocates all memory.
  SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_FINALISE(VAR_TO_CONSTRAINT_JACOBIAN_MAP,ERR,ERROR,*)
 
    !Argument variables
    TYPE(VAR_TO_CONSTRAINT_JACOBIAN_MAP_TYPE) :: VAR_TO_CONSTRAINT_JACOBIAN_MAP !<The variable to constraint Jacobian map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_DOF_TO_ROWS_MAP)) &
      & DEALLOCATE(VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_DOF_TO_ROWS_MAP)
    IF(ALLOCATED(VAR_TO_CONSTRAINT_JACOBIAN_MAP%CONSTRAINT_JACOBIAN_NUMBERS)) &
      & DEALLOCATE(VAR_TO_CONSTRAINT_JACOBIAN_MAP%CONSTRAINT_JACOBIAN_NUMBERS)
    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_FINALISE
 
  !
  !================================================================================================================================
  !
 
  !>Initialises a variable to equations Jacobian map
  SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_INIT(VAR_TO_CONSTRAINT_JACOBIAN_MAP,ERR,ERROR,*)
 
    !Argument variables
    TYPE(VAR_TO_CONSTRAINT_JACOBIAN_MAP_TYPE) :: VAR_TO_CONSTRAINT_JACOBIAN_MAP !<The variable to constraint Jacobian map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_INIT",ERR,ERROR,*999)
    
    VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE_TYPE=0
    NULLIFY(VAR_TO_CONSTRAINT_JACOBIAN_MAP%VARIABLE)
    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_INIT")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_INIT",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_INIT")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_JACOBIAN_MAP_INIT
 
  !
  !================================================================================================================================
  !
 
  !>Finalises a variable to constraint matrices map and deallocates all memory.
  SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE(VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*)
 
    !Argument variables
    TYPE(VAR_TO_CONSTRAINT_MATRICES_MAP_TYPE) :: VAR_TO_CONSTRAINT_MATRICES_MAP !<The variable to constraint matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS)) &
      & DEALLOCATE(VAR_TO_CONSTRAINT_MATRICES_MAP%CONSTRAINT_MATRIX_NUMBERS)
    IF(ALLOCATED(VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP)) &
      & DEALLOCATE(VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP)
    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_FINALISE
 
  !
  !================================================================================================================================
  !
 
  !>Initialise an constraint mapping constraint matrix map.
  SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT(VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*)
 
    !Argument variables
    TYPE(VAR_TO_CONSTRAINT_MATRICES_MAP_TYPE) :: VAR_TO_CONSTRAINT_MATRICES_MAP !<The variable to constraint matrices map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT",ERR,ERROR,*999)
 
    VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_TYPE=0
    NULLIFY(VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE)
    VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=0
    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT",ERR,ERROR)    
    CALL EXITS("CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT")
    RETURN 1
   
  END SUBROUTINE CONSTRAINT_MAPPING_VAR_TO_CONSTR_MATRICES_MAP_INIT
 
  !
  !================================================================================================================================
  !
 
  !>Sets the Lagrange variable type for the constraint mapping. 
  SUBROUTINE CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET(CONSTRAINT_MAPPING,LAGRANGE_VARIABLE_TYPE,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    INTEGER(INTG), INTENT(IN) :: LAGRANGE_VARIABLE_TYPE !<The Lagrange variable type to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: LAGRANGE_VARIABLE
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
          IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
            CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
            IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
              SELECT CASE(CONSTRAINT_CONDITION%METHOD)
              CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
                LAGRANGE=>CONSTRAINT_CONDITION%LAGRANGE
                IF(ASSOCIATED(LAGRANGE)) THEN
                  IF(LAGRANGE%LAGRANGE_FINISHED) THEN
                    LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                    NULLIFY(LAGRANGE_VARIABLE)
                    CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE,ERR,ERROR,*999)
                    CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE=LAGRANGE_VARIABLE_TYPE
                  ELSE
                    CALL FLAG_ERROR("Constraint condition Lagrange field has not been finished.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Constraint condition Lagrange is not associated.",ERR,ERROR,*999)
                ENDIF
              CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mapping constraint equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET
 
  !
  !================================================================================================================================
  !
 
  !>Sets the transpose flag for the constraint matrices. 
  SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MATRICES_TRANSPOSE_SET(CONSTRAINT_MAPPING,LINEAR_MATRIX_TRANSPOSE,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    LOGICAL, INTENT(IN) :: LINEAR_MATRIX_TRANSPOSE(:) !<LINEAR_MATRIX_TRANSPOSE(matrix_idx). The constraint linear matrix transpose flag for the matrix_idx'th constraint matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_LINEAR_MATRICES_TRANSPOSE_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
           CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
          IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
            CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
            IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
              SELECT CASE(CONSTRAINT_CONDITION%METHOD)
              CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
                !Check that the number of supplied coefficients matches the number of constraint matrices
                IF(SIZE(LINEAR_MATRIX_TRANSPOSE,1)==CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES) THEN
                  CREATE_VALUES_CACHE%LINEAR_HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES)= &
                    LINEAR_MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES)
                ELSE
                  LOCAL_ERROR="Invalid size of linear matrix tranpose. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(LINEAR_MATRIX_TRANSPOSE,1),"*",ERR,ERROR))// &
                    & ") must match the number of linear constraint matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mapping constraint equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MATRICES_TRANSPOSE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_LINEAR_MATRICES_TRANSPOSE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_LINEAR_MATRICES_TRANSPOSE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_LINEAR_MATRICES_TRANSPOSE_SET
 
  !
  !================================================================================================================================
  !
 
  !>Sets the transpose flag for the constraint matrices. 
  SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_TRANSPOSE_SET(CONSTRAINT_MAPPING,DYNAMIC_MATRIX_TRANSPOSE,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    LOGICAL, INTENT(IN) :: DYNAMIC_MATRIX_TRANSPOSE(:) !<DYNAMIC_MATRIX_TRANSPOSE(matrix_idx). The constraint dynamic matrix transpose flag for the matrix_idx'th constraint matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_TRANSPOSE_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
           CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
          IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
            CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
            IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
              SELECT CASE(CONSTRAINT_CONDITION%METHOD)
              CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
                !Check that the number of supplied coefficients matches the number of constraint matrices
                IF(SIZE(DYNAMIC_MATRIX_TRANSPOSE,1)==CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES) THEN
                  CREATE_VALUES_CACHE%DYNAMIC_HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES)= &
                    DYNAMIC_MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES)
                ELSE
                  LOCAL_ERROR="Invalid size of dynamic matrix tranpose. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(DYNAMIC_MATRIX_TRANSPOSE,1),"*",ERR,ERROR))// &
                    & ") must match the number of dynamic constraint matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mapping constraint equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_TRANSPOSE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_TRANSPOSE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_DYNAMIC_MATRICES_TRANSPOSE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_DYNAMIC_MATRICES_TRANSPOSE_SET
 
  !
  !================================================================================================================================
  !
 
  !>Sets the transpose flag for the constraint matrices. 
  SUBROUTINE CONSTRAINT_MAPPING_JACOBIAN_MATRICES_TRANSPOSE_SET(CONSTRAINT_MAPPING,JACOBIAN_MATRIX_TRANSPOSE,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    LOGICAL, INTENT(IN) :: JACOBIAN_MATRIX_TRANSPOSE(:) !<JACOBIAN_MATRIX_TRANSPOSE(matrix_idx). The constraint jacobian matrix transpose flag for the matrix_idx'th constraint matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_JACOBIAN_MATRICES_TRANSPOSE_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
           CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
          IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
            CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
            IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
              SELECT CASE(CONSTRAINT_CONDITION%METHOD)
              CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
                !Check that the number of supplied coefficients matches the number of constraint matrices
                IF(SIZE(JACOBIAN_MATRIX_TRANSPOSE,1)==CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES) THEN
                  CREATE_VALUES_CACHE%JACOBIAN_HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES)= &
                    JACOBIAN_MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES)
                ELSE
                  LOCAL_ERROR="Invalid size of jacobian matrix tranpose. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(JACOBIAN_MATRIX_TRANSPOSE,1),"*",ERR,ERROR))// &
                    & ") must match the number of jacobian constraint matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The constraint condition method of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mapping constraint equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_JACOBIAN_MATRICES_TRANSPOSE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_JACOBIAN_MATRICES_TRANSPOSE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_JACOBIAN_MATRICES_TRANSPOSE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_JACOBIAN_MATRICES_TRANSPOSE_SET
 
  !
  !================================================================================================================================
  !
 
  !>Sets the coefficient applied to the constraint RHS vector.
  SUBROUTINE CONSTRAINT_MAPPING_RHS_COEFF_SET(CONSTRAINT_MAPPING,RHS_COEFFICIENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set the RHS coefficent for
    REAL(DP), INTENT(IN) :: RHS_COEFFICIENT !<The coefficient applied to the constraint RHS vector.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_RHS_COEFF_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE)) THEN
          IF(CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE/=0) THEN
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RHS_COEFFICIENT=RHS_COEFFICIENT
          ELSE
            CALL FLAG_ERROR("The constraint mapping RHS Lagrange variable type has not been set.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_RHS_COEFF_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_RHS_COEFF_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_RHS_COEFF_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_RHS_COEFF_SET
  
  !
  !================================================================================================================================
  !
 
  !>Finalises the constraint mapping RHS mapping and deallocates all memory
  SUBROUTINE CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE(RHS_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING !<A pointer to the RHS mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE",ERR,ERROR,*999)
 
    IF(ASSOCIATED(RHS_MAPPING)) THEN
      IF(ALLOCATED(RHS_MAPPING%RHS_DOF_TO_CONSTRAINT_ROW_MAP)) DEALLOCATE(RHS_MAPPING%RHS_DOF_TO_CONSTRAINT_ROW_MAP)
      IF(ALLOCATED(RHS_MAPPING%CONSTRAINT_ROW_TO_RHS_DOF_MAP)) DEALLOCATE(RHS_MAPPING%CONSTRAINT_ROW_TO_RHS_DOF_MAP)
      DEALLOCATE(RHS_MAPPING)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE
 
  !
  !================================================================================================================================
  !
 
  !>Initialises the constraint mapping RHS mapping
  SUBROUTINE CONSTRAINT_MAPPING_RHS_MAPPING_INITIALISE(CONSTRAINT_MAPPING,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to initialise the RHS mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(ASSOCIATED(CONSTRAINT_MAPPING%RHS_MAPPING)) THEN
        CALL FLAG_ERROR("Constraint mapping RHS mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONSTRAINT_MAPPING%RHS_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint mapping RHS mapping.",ERR,ERROR,*999)
        CONSTRAINT_MAPPING%RHS_MAPPING%CONSTRAINT_MAPPING=>CONSTRAINT_MAPPING        
        CONSTRAINT_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE=0
        NULLIFY(CONSTRAINT_MAPPING%RHS_MAPPING%RHS_VARIABLE)
        NULLIFY(CONSTRAINT_MAPPING%RHS_MAPPING%RHS_VARIABLE_MAPPING)
        CONSTRAINT_MAPPING%RHS_MAPPING%RHS_COEFFICIENT=1.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN
999 CALL CONSTRAINT_MAPPING_RHS_MAPPING_FINALISE(CONSTRAINT_MAPPING%RHS_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONSTRAINT_MAPPING_RHS_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_RHS_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_RHS_MAPPING_INITIALISE
 
  !
  !================================================================================================================================
  !
 
  !>Sets the mapping between a Lagrange field variable and the constraint rhs vector.
  SUBROUTINE CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET(CONSTRAINT_MAPPING,RHS_VARIABLE_TYPE,ERR,ERROR,*)
 
    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set the RHS variable type for.
    INTEGER(INTG), INTENT(IN) :: RHS_VARIABLE_TYPE !<The variable type associated with the constraint rhs vector. If the constraint condition does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_LAGRANGE_TYPE), POINTER :: CONSTRAINT_LAGRANGE
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR,*999)
 
    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          IF(RHS_VARIABLE_TYPE==0) THEN
            CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=0
          ELSE
            CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
            IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
              CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
              IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
                SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
                  CONSTRAINT_LAGRANGE=>CONSTRAINT_CONDITION%LAGRANGE
                  IF(ASSOCIATED(CONSTRAINT_LAGRANGE)) THEN
                    LAGRANGE_FIELD=>CONSTRAINT_LAGRANGE%LAGRANGE_FIELD
                    IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                      !Check the RHS variable type is not being by the constraint matrices
                      IF(CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE==RHS_VARIABLE_TYPE) THEN
                        LOCAL_ERROR="The specified RHS variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " is the same as the Lagrange variable type for the constraint matrices."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Check the RHS variable number is defined on the Lagrange field
                      IF(RHS_VARIABLE_TYPE>=1.AND.RHS_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                        IF(ASSOCIATED(LAGRANGE_FIELD%VARIABLE_TYPE_MAP(RHS_VARIABLE_TYPE)%PTR)) THEN
                          CREATE_VALUES_CACHE%RHS_LAGRANGE_VARIABLE_TYPE=RHS_VARIABLE_TYPE
                        ELSE
                          LOCAL_ERROR="The specified RHS variable type of "// &
                            & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " is not defined on the Lagrange field."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="The specified RHS variable type of "// &
                          & TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " is invalid. The number must either be zero or >= 1 and <= "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Lagrange field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Constraint Lagrange is not associated.",ERR,ERROR,*999)
                  ENDIF
                CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The constraint condition method of "// &
                    & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_CONDITION%METHOD,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Constraint equations constraint condition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint mapping constraint equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mapping create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

END MODULE CONSTRAINT_MAPPING_ROUTINES
