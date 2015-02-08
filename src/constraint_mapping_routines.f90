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

  PUBLIC CONSTRAINT_MAPPING_CREATE_FINISH,CONSTRAINT_MAPPING_CREATE_START

  PUBLIC CONSTRAINT_MAPPING_DESTROY

  PUBLIC CONSTRAINT_MAPPING_LAGRANGE_VARIABLE_TYPE_SET

  PUBLIC CONSTRAINT_MAPPING_MATRICES_COEFFS_SET

  PUBLIC CONSTRAINT_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET,CONSTRAINT_MAPPING_MATRICES_ROW_MESH_INDICES_SET

  PUBLIC CONSTRAINT_MAPPING_MATRICES_NUMBER_SET

  PUBLIC CONSTRAINT_MAPPING_MATRICES_TRANSPOSE_SET

  PUBLIC CONSTRAINT_MAPPING_RHS_COEFF_SET

  PUBLIC CONSTRAINT_MAPPING_RHS_VARIABLE_TYPE_SET

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
    INTEGER(INTG) :: column_idx,dof_idx,matrix_idx,mesh_idx,variable_idx,number_of_constraint_matrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,LAGRANGE_VARIABLE
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(DEPENDENT_TYPE), POINTER :: DEPENDENT
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
              DEPENDENT=>CONSTRAINT_CONDITION%DEPENDENT
              IF(ASSOCIATED(DEPENDENT)) THEN
                !Set the Lagrange variable information
                LAGRANGE_FIELD=>LAGRANGE%LAGRANGE_FIELD
                DEPENDENT_FIELD=>DEPENDENT%DEPENDENT_FIELD
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
                    CALL CONSTRAINT_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE( &
                      & DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*999)
                    DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_TYPE=CREATE_VALUES_CACHE%VARIABLE_TYPE
                    DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE=>CREATE_VALUES_CACHE%VARIABLE
                    DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                    IF(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=-1
                    !Allocate and initialise the variable to constraint matrices maps
                    FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE
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
                        row_idx=dof_idx
                        DYNAMIC_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP(dof_idx)=row_idx
                        ENDDO !dof_idx
                      ENDIF
                    ENDIF
                    !Allocate and initialise the constraint matrix to variable maps types
                    ALLOCATE(DYNAMIC_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(DYNAMIC_MAPPING% &
                      & NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrix rows to variable maps.", &
                      & ERR,ERROR,*999)
                    !Loop over the constraint matrices and calculate the row mappings
                    !The pointers below have been checked for association above.
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      number_of_constraint_matrices=CONSTRAINT_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      !Number of constraint matrices whose rows/columns are related to Dependent/Lagrange variables and not Lagrange/Lagrange variables (last constraint matrix is Lagrange/Lagrange (Penalty matrix)
                      number_of_constraint_matrices=CONSTRAINT_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES/2 
                    ENDSELECT
                    DO matrix_idx=1,number_of_constraint_matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_MATRIX_ROW_TO_VAR_MAP_INITIALISE(CONSTRAINT_MAPPING,matrix_idx,ERR,ERROR,*999)
                      EQUATIONS_SET=>DEPENDENT%EQUATIONS_SET
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        FIELD_VARIABLE=>DEPENDENT%FIELD_VARIABLE
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%EQUATIONS_SET=>EQUATIONS_SET
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE= &
                            & FIELD_VARIABLE%VARIABLE_TYPE
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT= &
                            & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                           !Set the number of rows
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)% &
                              & VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
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
                      DO matrix_idx=number_of_constraint_matrices,CONSTRAINT_MAPPING%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES
                        !Initialise and setup the constraint matrix
                        CALL CONSTRAINT_MAPPING_MATRIX_ROW_TO_VAR_MAP_INITIALISE(CONSTRAINT_MAPPING,matrix_idx,ERR,ERROR,*999)
                        CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                          & ERR,ERROR,*999)
                        FIELD_VARIABLE=>LAGRANGE_VARIABLE
                        CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
                        IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
                          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%CONSTRAINT_EQUATIONS=> &
                              & CONSTRAINT_EQUATIONS
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE% &
                              & VARIABLE_TYPE
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT= &
                              & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                              & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                              !Set the number of rows
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS= &
                              & FIELD_VARIABLE%NUMBER_OF_DOFS
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                              & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                              & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                            !Set the row mapping
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                              & FIELD_VARIABLE%DOMAIN_MAPPING
                            ALLOCATE(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                              & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                            !1-1 mapping for now
                            DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                              CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)% &
                                & VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
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
                    CALL CONSTRAINT_MAPPING_VAR_TO_EQUATS_MATRICES_MAP_INITIALISE( &
                      & LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP,ERR,ERROR,*999)
                    LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_TYPE=CREATE_VALUES_CACHE%VARIABLE_TYPE
                    LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE=>CREATE_VALUES_CACHE%VARIABLE
                    LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
                    IF(CREATE_VALUES_CACHE%RHS_VARIABLE_TYPE/=0) LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%NUMBER_OF_CONSTRAINT_MATRICES=-1
                    !Allocate and initialise the variable to constraint matrices maps
                    FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE
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
                        row_idx=dof_idx
                        LINEAR_MAPPING%VAR_TO_CONSTRAINT_MATRICES_MAP%VARIABLE_DOF_TO_ROWS_MAP(dof_idx)=row_idx
                        ENDDO !dof_idx
                      ENDIF
                    ENDIF
                    !Allocate and initialise the constraint matrix to variable maps types
                    ALLOCATE(LINEAR_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(LINEAR_MAPPING% &
                      & NUMBER_OF_LINEAR_CONSTRAINT_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint matrix rows to variable maps.", &
                      & ERR,ERROR,*999)
                    !Loop over the constraint matrices and calculate the row mappings
                    !The pointers below have been checked for association above.
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      number_of_constraint_matrices=CONSTRAINT_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      !Number of constraint matrices whose rows/columns are related to Dependent/Lagrange variables and not Lagrange/Lagrange variables (last constraint matrix is Lagrange/Lagrange (Penalty matrix)
                      number_of_constraint_matrices=CONSTRAINT_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES-1 
                    ENDSELECT
                    DO matrix_idx=1,number_of_constraint_matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_MATRIX_ROW_TO_VAR_MAP_INITIALISE(CONSTRAINT_MAPPING,matrix_idx,ERR,ERROR,*999)
                      EQUATIONS_SET=>DEPENDENT%EQUATIONS_SET
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        FIELD_VARIABLE=>DEPENDENT%FIELD_VARIABLE
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%EQUATIONS_SET=>EQUATIONS_SET
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                           !Set the number of rows
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE%NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
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
                      matrix_idx = CONSTRAINT_MAPPING%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES !last of the constraint matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_MATRIX_ROW_TO_VAR_MAP_INITIALISE(CONSTRAINT_MAPPING,matrix_idx,ERR,ERROR,*999)
                      CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                        & ERR,ERROR,*999)
                      FIELD_VARIABLE=>LAGRANGE_VARIABLE
                      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
                      IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%CONSTRAINT_EQUATIONS=>CONSTRAINT_EQUATIONS
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                            !Set the number of rows
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE%NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)% &
                              & VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
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
                    !Allocate and initialise the Jacobian constraint matrix to variable maps types
                    ALLOCATE(NONLINEAR_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(LINEAR_MAPPING% &
                      & NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint Jacobian rows to variable maps.", &
                      & ERR,ERROR,*999)
                    !Loop over the constraint matrices and calculate the row mappings
                    !The pointers below have been checked for association above.
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      number_of_constraint_matrices=CONSTRAINT_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      !Number of constraint matrices whose rows/columns are related to Dependent/Lagrange variables and not Lagrange/Lagrange variables (last constraint matrix is Lagrange/Lagrange (Penalty matrix)
                      number_of_constraint_matrices=CONSTRAINT_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES-1 
                    ENDSELECT
                    DO matrix_idx=1,number_of_constraint_matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_JACOBIAN_ROWS_TO_VAR_MAP_INITIALISE(CONSTRAINT_MAPPING,matrix_idx,ERR,ERROR,*999)
                      EQUATIONS_SET=>DEPENDENT%EQUATIONS_SET
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        FIELD_VARIABLE=>DEPENDENT%FIELD_VARIABLE
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%EQUATIONS_SET=>EQUATIONS_SET
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE= &
                            & FIELD_VARIABLE%VARIABLE_TYPE
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%RESIDUAL_COEFFICIENT= &
                            & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENT
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                           !Set the number of rows
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
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
                      matrix_idx = CONSTRAINT_MAPPING%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES !last of the constraint matrices
                      !Initialise and setup the constraint matrix
                      CALL CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE(CONSTRAINT_MAPPING,matrix_idx,ERR,ERROR,*999)
                      NULLIFY(LAGRANGE_VARIABLE)
                      CALL FIELD_VARIABLE_GET(LAGRANGE_FIELD,CREATE_VALUES_CACHE%LAGRANGE_VARIABLE_TYPE,LAGRANGE_VARIABLE, &
                        & ERR,ERROR,*999)
                      NULLIFY(CONSTRAINT_EQUATIONS)
                      NULLIFY(FIELD_VARIABLE)
                      FIELD_VARIABLE=>LAGRANGE_VARIABLE
                      CONSTRAINT_EQUATIONS=>CONSTRAINT_CONDITION%CONSTRAINT_EQUATIONS
                      IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%CONSTRAINT_EQUATIONS=> &
                            & CONSTRAINT_EQUATIONS
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=FIELD_VARIABLE% &
                            & VARIABLE_TYPE
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE=>FIELD_VARIABLE
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%RESIDUAL_COEFFICIENT= &
                            & CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENT
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=CONSTRAINT_MAPPING% &
                            & CREATE_VALUES_CACHE%HAS_TRANSPOSE(matrix_idx)
                            !Set the number of rows
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS=FIELD_VARIABLE% &
                            & NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS= &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS= &
                            & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS
                          !Set the row mapping
                          CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING=> &
                            & FIELD_VARIABLE%DOMAIN_MAPPING
                          ALLOCATE(CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP( &
                            & FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dof to row map.",ERR,ERROR,*999)
                          !1-1 mapping for now
                          DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                            CONSTRAINT_MAPPING%CONSTRAINT_JACOBIAN_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_DOF_TO_ROW_MAP(dof_idx)=dof_idx
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
                      column_idx=dof_idx
                      RHS_MAPPING%RHS_DOF_TO_CONSTRAINT_ROW_MAP(dof_idx)=column_idx
                    ENDDO !dof_idx
                    DO column_idx=1,CONSTRAINT_MAPPING%TOTAL_NUMBER_OF_COLUMNS
                      !1-1 mapping for now
                      dof_idx=column_idx
                      RHS_MAPPING%CONSTRAINT_ROW_TO_RHS_DOF_MAP(column_idx)=dof_idx
                    ENDDO !column_idx
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
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%HAS_TRANSPOSE)) DEALLOCATE(CREATE_VALUES_CACHE%HAS_TRANSPOSE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES))  &
        & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COL_FIELD_VARIABLE_INDICES)) &
        & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COL_FIELD_VARIABLE_INDICES)
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
    INTEGER(INTG) :: DUMMY_ERR,variable_idx,variable_type_idx,variable_type_idx2
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
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
            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%RESIDUAL_COEFFICIENT=1.0_DP
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
                    CONSTRAINT_DEPENDENT_FIELD=>CONSTRAINT_DEPENDENT=>DEPENDENT_FIELD
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
                        & CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix coefficients.",ERR,ERROR,*999)
                      !Set the default constraint mapping in the create values cache
                      SELECT CASE(CONSTRAINT_EQATIONS%TIME_DEPENDENCE)
                      CASE(CONSTRAINT_CONDITION_STATIC)
                        SELECT CASE(CONSTRAINT_EQATIONS%LINEARITY)
                        CASE(CONSTRAINT_CONDITION_LINEAR,CONSTRAINT_CONDITION_NONLINEAR_BCS)
                          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=1
                          CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables plus a single Lagrange variable
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=2
                          END SELECT
                        CASE(CONSTRAINT_CONDITION_NONLINEAR)
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_LINEAR_CONSTRAINT_MATRICES=0
                          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                            !Default the number of Jacobia constraint matrices to the number of added dependent variables
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES=1
                          CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                            !Default the number of Jacobia constraint matrices to the number of added dependent variables plus a single Lagrange variable
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_JACOBIAN_CONSTRAINT_MATRICES=2
                          END SELECT
                        CASE DEFAULT
                          LOCAL_ERROR="The constraint linearity type of "// &
                            & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      CASE(CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC,CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC)
                        IF(CONSTRAINT_EQATIONS%TIME_DEPENDENCE==CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC) THEN
                          CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%DYNAMIC_DAMPING_MATRIX_NUMBER=2
                          SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                          CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables
                            CONSTRAINT_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_DYNAMIC_CONSTRAINT_MATRICES=2
                          CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                            !Default the number of constraint matrices to the number of added dependent variables plus a single Lagrange variable
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
                          & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_EQATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
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

  !>Finalises the constraint mapping and deallocates all memory.
  SUBROUTINE CONSTRAINT_MAPPING_FINALISE(CONSTRAINT_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("CONSTRAINT_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(ALLOCATED(CONSTRAINT_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)) DEALLOCATE(CONSTRAINT_MAPPING%LAGRANGE_DOF_TO_COLUMN_MAP)
      IF(ALLOCATED(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS)) THEN
        DO matrix_idx=1,SIZE(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS,1)
          CALL CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_FINALISE(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx), &
            & ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS)
      ENDIF
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
        CONSTRAINT_EQUATIONS%CONSTRAINT_MAPPING%NUMBER_OF_CONSTRAINT_MATRICES=0
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

  !>Finalises a constraint matrix to variable map and deallocates all memory.
  SUBROUTINE CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_FINALISE(CONSTRAINT_MATRIX_TO_VAR_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MATRIX_TO_VAR_MAP_TYPE) :: CONSTRAINT_MATRIX_TO_VAR_MAP !<The constraint matrix to var map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    CALL ENTERS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(CONSTRAINT_MATRIX_TO_VAR_MAP%VARIABLE_DOF_TO_ROW_MAP)) &
      & DEALLOCATE(CONSTRAINT_MATRIX_TO_VAR_MAP%VARIABLE_DOF_TO_ROW_MAP)
    
    CALL EXITS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a constraint matrix to variable map.
  SUBROUTINE CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE(CONSTRAINT_MAPPING,matrix_idx,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to initialise the matrix to variable map for a given matrix index.
    INTEGER(INTG), INTENT(IN) :: matrix_idx !<The matrix index to intialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
  
    CALL ENTERS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(ALLOCATED(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS)) THEN
        IF(matrix_idx>0.AND.matrix_idx<=CONSTRAINT_MAPPING%NUMBER_OF_CONSTRAINT_MATRICES) THEN
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_NUMBER=matrix_idx
          NULLIFY(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%CONSTRAINT_MATRIX)
          NULLIFY(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%EQUATIONS_SET)
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE_TYPE=0
          NULLIFY(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE)
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MESH_INDEX=0
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%MATRIX_COEFFICIENT=0.0_DP
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%HAS_TRANSPOSE=.FALSE.
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_ROWS=0
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%TOTAL_NUMBER_OF_ROWS=0
          CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%NUMBER_OF_GLOBAL_ROWS=0
          NULLIFY(CONSTRAINT_MAPPING%CONSTRAINT_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING)          
        ELSE
          LOCAL_ERROR="The specified matrix index of "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is invalid. The index must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(CONSTRAINT_MAPPING%NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint mapping matrix rows to var maps is not allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_MATRIX_TO_VAR_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the constraint matrices. 
  SUBROUTINE CONSTRAINT_MAPPING_MATRICES_COEFFS_SET(CONSTRAINT_MAPPING,MATRIX_COEFFICIENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    REAL(DP), INTENT(IN) :: MATRIX_COEFFICIENTS(:) !<The constraint matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MAPPING_MATRICES_COEFFS_SET",ERR,ERROR,*999)

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
                IF(SIZE(MATRIX_COEFFICIENTS,1)==CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES) THEN
                  CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)= &
                    & MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)
                ELSE
                  LOCAL_ERROR="Invalid size of matrix coefficeints. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_COEFFICIENTS,1),"*",ERR,ERROR))// &
                    & ") must match the number of constraint matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
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
    
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_COEFFS_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_MATRICES_COEFFS_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_COEFFS_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_MATRICES_COEFFS_SET

  !
  !================================================================================================================================
  !

  !>Sets the column mesh indices for the constraint matrices. 
  SUBROUTINE CONSTRAINT_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET(CONSTRAINT_MAPPING,COLUMN_MESH_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    INTEGER(INTG), INTENT(IN) :: COLUMN_MESH_INDICES(:) !<The constraint matrix column mesh indices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET",ERR,ERROR,*999)

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
              CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                CALL FLAG_ERROR("Can not set the column mesh indices when using the Lagrange multipliers "// &
                  "constraint condition method.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
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
    
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_MATRICES_COLUMN_MESH_INDICES_SET

  !
  !================================================================================================================================
  !

  !>Sets the number of constraint matrices for a constraint mapping.
  SUBROUTINE CONSTRAINT_MAPPING_MATRICES_NUMBER_SET(CONSTRAINT_MAPPING,NUMBER_OF_CONSTRAINT_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping to set the number of linear matrices for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_CONSTRAINT_MATRICES !<The number of constraint matrices to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx,matrix_idx2,variable_idx,number_of_dependent_variables
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(:)
    LOGICAL :: FOUND
    LOGICAL, ALLOCATABLE :: OLD_MATRIX_TRANSPOSE(:)
    REAL(DP), ALLOCATABLE :: OLD_MATRIX_COEFFICIENTS(:)
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_DEPENDENT_TYPE), POINTER :: CONSTRAINT_DEPENDENT
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MAPPING)) THEN
      IF(CONSTRAINT_MAPPING%CONSTRAINT_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mapping has already been finished.",ERR,ERROR,*999)
      ELSE
        CREATE_VALUES_CACHE=>CONSTRAINT_MAPPING%CREATE_VALUES_CACHE
        IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
          CONSTRAINT_EQUATIONS=>CONSTRAINT_MAPPING%CONSTRAINT_EQUATIONS
          IF(ASSOCIATED(CONSTRAINT_EQUATIONS)) THEN
            CONSTRAINT_CONDITION=>CONSTRAINT_EQUATIONS%CONSTRAINT_CONDITION
            IF(ASSOCIATED(CONSTRAINT_CONDITION)) THEN
              SELECT CASE(CONSTRAINT_CONDITION%METHOD)
              CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,CONSTRAINT_CONDITION_PENALTY_METHOD)
                !Check the number of constraint matrices
                IF(NUMBER_OF_CONSTRAINT_MATRICES>0) THEN
                  CONSTRAINT_DEPENDENT=>CONSTRAINT_CONDITION%DEPENDENT
                  IF(ASSOCIATED(CONSTRAINT_DEPENDENT)) THEN
                    SELECT CASE(CONSTRAINT_CONDITION%METHOD)
                    CASE(CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                      number_of_dependent_variables=CONSTRAINT_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    CASE(CONSTRAINT_CONDITION_PENALTY_METHOD)
                      number_of_dependent_variables=CONSTRAINT_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1
                    END SELECT
                    IF(NUMBER_OF_CONSTRAINT_MATRICES<=number_of_dependent_variables) THEN
                      !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
                      IF(NUMBER_OF_CONSTRAINT_MATRICES/=CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES) THEN
                        ALLOCATE(OLD_MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix coefficients.",ERR,ERROR,*999)
                        ALLOCATE(OLD_MATRIX_TRANSPOSE(CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix transpose.",ERR,ERROR,*999)
                        ALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix row field indexes.",ERR,ERROR,*999)
                        OLD_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)= &
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)
                        OLD_MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)= &
                          & CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)
                        OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)= &
                          & CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE% &
                          & NUMBER_OF_CONSTRAINT_MATRICES)
                        IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)) &
                          & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
                        IF(ALLOCATED(CREATE_VALUES_CACHE%HAS_TRANSPOSE)) &
                          & DEALLOCATE(CREATE_VALUES_CACHE%HAS_TRANSPOSE)
                        IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES)) &
                          & DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES)
                        ALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix coefficients.",ERR,ERROR,*999)
                        ALLOCATE(CREATE_VALUES_CACHE%HAS_TRANSPOSE(NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix tranpose.",ERR,ERROR,*999)
                        ALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(NUMBER_OF_CONSTRAINT_MATRICES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix row field variable indexes.",ERR,ERROR,*999)
                        IF(NUMBER_OF_CONSTRAINT_MATRICES>CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES) THEN
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)= &
                            & OLD_MATRIX_COEFFICIENTS(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES+1: &
                            & NUMBER_OF_CONSTRAINT_MATRICES)=1.0_DP
                          CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)= &
                            & OLD_MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)
                          CREATE_VALUES_CACHE%HAS_TRANSPOSE(CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES+1: &
                            & NUMBER_OF_CONSTRAINT_MATRICES)=.TRUE.
                          CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE% &
                            & NUMBER_OF_CONSTRAINT_MATRICES)=OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:CREATE_VALUES_CACHE% &
                            & NUMBER_OF_CONSTRAINT_MATRICES)
                          !Loop through in mesh index order and set the default matrix to variable map to be in mesh index order
                          DO matrix_idx=CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES+1,NUMBER_OF_CONSTRAINT_MATRICES
                            CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)=0
                            DO variable_idx=1,CONSTRAINT_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                              FOUND=.FALSE.
                              DO matrix_idx2=1,CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES
                                IF(CONSTRAINT_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)==CREATE_VALUES_CACHE% &
                                  MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx2)) THEN
                                  FOUND=.TRUE.
                                  EXIT
                                ENDIF
                              ENDDO !matrix_idx2
                              IF(.NOT.FOUND) THEN
                                CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)=CONSTRAINT_DEPENDENT% &
                                  & VARIABLE_MESH_INDICES(variable_idx)
                              ENDIF
                            ENDDO !variable_idx2
                            IF(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(matrix_idx)==0) THEN
                              LOCAL_ERROR="Could not map a constraint mesh index for constraint matrix "// &
                                & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//"."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDDO !matrix_idx
                        ELSE
                          CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:NUMBER_OF_CONSTRAINT_MATRICES)= &
                            & OLD_MATRIX_COEFFICIENTS(1:NUMBER_OF_CONSTRAINT_MATRICES)
                          CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:NUMBER_OF_CONSTRAINT_MATRICES)= &
                            & OLD_MATRIX_TRANSPOSE(1:NUMBER_OF_CONSTRAINT_MATRICES)
                          CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_VARIABLE_INDICES(1:NUMBER_OF_CONSTRAINT_MATRICES)= &
                            & OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES(1:NUMBER_OF_CONSTRAINT_MATRICES)
                        ENDIF
                        IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
                        IF(ALLOCATED(OLD_MATRIX_TRANSPOSE)) DEALLOCATE(OLD_MATRIX_TRANSPOSE)
                        IF(ALLOCATED(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)) DEALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The specified number of constraint matrices of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))// &
                        & " is invalid. The number must be <= the number of added dependent variables of "// &
                        & TRIM(NUMBER_TO_VSTRING(number_of_dependent_variables,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Constraint condition dependent is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified number of constraint matrices of "// &
                    & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))// &
                    & " is invalid. The number must be > 0."
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
      CALL FLAG_ERROR("Constraint mapping is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
    IF(ALLOCATED(OLD_MATRIX_TRANSPOSE)) DEALLOCATE(OLD_MATRIX_TRANSPOSE)
    IF(ALLOCATED(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)) DEALLOCATE(OLD_MATRIX_ROW_FIELD_VARIABLE_INDICES)
    CALL ERRORS("CONSTRAINT_MAPPING_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_NUMBER_SET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MAPPING_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Sets the transpose flag for the constraint matrices. 
  SUBROUTINE CONSTRAINT_MAPPING_MATRICES_TRANSPOSE_SET(CONSTRAINT_MAPPING,MATRIX_TRANSPOSE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MAPPING_TYPE), POINTER :: CONSTRAINT_MAPPING !<A pointer to the constraint mapping.
    LOGICAL, INTENT(IN) :: MATRIX_TRANSPOSE(:) !<MATRIX_TRANSPOSE(matrix_idx). The constraint matrix transpose flag for the matrix_idx'th constraint matrix.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: CONSTRAINT_CONDITION
    TYPE(CONSTRAINT_EQUATIONS_TYPE), POINTER :: CONSTRAINT_EQUATIONS
    TYPE(CONSTRAINT_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MAPPING_MATRICES_TRANSPOSE_SET",ERR,ERROR,*999)

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
                IF(SIZE(MATRIX_TRANSPOSE,1)==CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES) THEN
                  CREATE_VALUES_CACHE%HAS_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)= &
                    MATRIX_TRANSPOSE(1:CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES)
                ELSE
                  LOCAL_ERROR="Invalid size of matrix tranpose. The size of the supplied array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_TRANSPOSE,1),"*",ERR,ERROR))// &
                    & ") must match the number of constraint matrices ("// &
                    & TRIM(NUMBER_TO_VSTRING(CREATE_VALUES_CACHE%NUMBER_OF_CONSTRAINT_MATRICES,"*",ERR,ERROR))//")."
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
    
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_TRANSPOSE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MAPPING_MATRICES_TRANSPOSE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MAPPING_MATRICES_TRANSPOSE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MAPPING_MATRICES_TRANSPOSE_SET

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
