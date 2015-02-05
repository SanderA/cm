!> \file
!> \author Chris Bradley
!> \brief This module contains all constraint routines.
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
!> Contributor(s): David Nordsletten, Thiranja Prasad Babarenda Gamage
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delte
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module contains all constraint routines.
MODULE CONSTRAINT_ROUTINES

  USE BASE_ROUTINES
  USE DATA_POINT_ROUTINES
  USE DATA_PROJECTION_ROUTINES
  USE FIELD_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE CONSTRAINT_CONDITIONS_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Constraints

  CONSTRAINT CONSTRAINT_LABEL_GET
    MODULE PROCEDURE CONSTRAINT_LABEL_GET_C
    MODULE PROCEDURE CONSTRAINT_LABEL_GET_VS
  END CONSTRAINT !CONSTRAINT_LABEL_GET
  
  CONSTRAINT CONSTRAINT_LABEL_SET
    MODULE PROCEDURE CONSTRAINT_LABEL_SET_C
    MODULE PROCEDURE CONSTRAINT_LABEL_SET_VS
  END CONSTRAINT !CONSTRAINT_LABEL_SET
  
  PUBLIC CONSTRAINT_MESH_ADD

  PUBLIC CONSTRAINT_CREATE_START, CONSTRAINT_CREATE_FINISH
  
  PUBLIC CONSTRAINT_COORDINATE_SYSTEM_SET,CONSTRAINT_COORDINATE_SYSTEM_GET
  
  PUBLIC CONSTRAINT_DATA_POINTS_GET

  PUBLIC CONSTRAINT_DESTROY, CONSTRAINT_MESH_CONNECTIVITY_DESTROY, ConstraintPointsConnectivity_Destroy

  PUBLIC CONSTRAINT_LABEL_GET,CONSTRAINT_LABEL_SET

  PUBLIC CONSTRAINT_NODES_GET

  PUBLIC CONSTRAINT_MESH_CONNECTIVITY_CREATE_START, CONSTRAINT_MESH_CONNECTIVITY_CREATE_FINISH

  PUBLIC CONSTRAINT_USER_NUMBER_FIND

  PUBLIC CONSTRAINTS_FINALISE,CONSTRAINTS_INITIALISE

  PUBLIC CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_XI_SET, CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET
  
  PUBLIC CONSTRAINT_MESH_CONNECTIVITY_NODE_NUMBER_SET

  PUBLIC CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET
  
  PUBLIC ConstraintPointsConnectivity_CreateStart,ConstraintPointsConnectivity_CreateFinish
  
  PUBLIC ConstraintPointsConnectivity_DataReprojection
  
  PUBLIC ConstraintPointsConnectivity_ElementNumberGet,ConstraintPointsConnectivity_ElementNumberSet
  
  PUBLIC ConstraintPointsConnectivity_PointXiGet,ConstraintPointsConnectivity_PointXiSet
  
  PUBLIC ConstraintPointsConnectivity_UpdateFromProjection
  
CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE CONSTRAINT_MESH_ADD(CONSTRAINT,MESH,MESH_INDEX,ERR,ERROR,*)   

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to add a mesh to
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to add to the constraint
    INTEGER(INTG), INTENT(OUT) :: MESH_INDEX !<On return, the index of the added mesh in the list of meshes in the constraint
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx
    LOGICAL :: MESH_ALREADY_COUPLED
    TYPE(MESH_TYPE), POINTER :: COUPLED_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_COUPLED_MESHES(:)
    TYPE(REGION_TYPE), POINTER :: COUPLED_MESH_REGION,MESH_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_COUPLED_MESHES)
    
    CALL ENTERS("CONSTRAINT_MESH_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN
        CALL FLAG_ERROR("Constraint has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(MESH)) THEN
          IF(MESH%MESH_FINISHED) THEN
            MESH_REGION=>MESH%REGION
            IF(ASSOCIATED(MESH_REGION)) THEN
              ALLOCATE(NEW_COUPLED_MESHES(CONSTRAINT%NUMBER_OF_COUPLED_MESHES+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new coupled meshes.",ERR,ERROR,*999)
              !Check that the mesh is not already in the list of meshes for the constraint.
              IF(CONSTRAINT%NUMBER_OF_COUPLED_MESHES>0) THEN
                IF(ASSOCIATED(CONSTRAINT%COUPLED_MESHES)) THEN
                  MESH_ALREADY_COUPLED=.FALSE.
                  DO mesh_idx=1,CONSTRAINT%NUMBER_OF_COUPLED_MESHES
                    COUPLED_MESH=>CONSTRAINT%COUPLED_MESHES(mesh_idx)%PTR
                    IF(ASSOCIATED(COUPLED_MESH)) THEN
                      COUPLED_MESH_REGION=>COUPLED_MESH%REGION
                      IF(ASSOCIATED(COUPLED_MESH_REGION)) THEN
                        IF(MESH_REGION%USER_NUMBER==COUPLED_MESH_REGION%USER_NUMBER) THEN
                          IF(MESH%USER_NUMBER==COUPLED_MESH%USER_NUMBER) THEN
                            MESH_ALREADY_COUPLED=.TRUE.
                            EXIT
                          ENDIF
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Coupled constraint mesh region for mesh index "// &
                          & TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))//" is not associated."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Coupled constraint mesh for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))// &
                        & " is not associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    NEW_COUPLED_MESHES(mesh_idx)%PTR=>CONSTRAINT%COUPLED_MESHES(mesh_idx)%PTR
                  ENDDO !mesh_idx
                  IF(MESH_ALREADY_COUPLED) THEN
                    LOCAL_ERROR="The supplied mesh has already been added to the list of coupled meshes at mesh index "// &
                      & TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  DEALLOCATE(CONSTRAINT%COUPLED_MESHES)
                ELSE
                  CALL FLAG_ERROR("Constraint coupled meshes is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              !Add the mesh to the list of coupled meshes
              NEW_COUPLED_MESHES(CONSTRAINT%NUMBER_OF_COUPLED_MESHES+1)%PTR=>MESH
              CONSTRAINT%COUPLED_MESHES=>NEW_COUPLED_MESHES
              !Increment the number of coupled meshes and return the index
              CONSTRAINT%NUMBER_OF_COUPLED_MESHES=CONSTRAINT%NUMBER_OF_COUPLED_MESHES+1
              MESH_INDEX=CONSTRAINT%NUMBER_OF_COUPLED_MESHES
            ELSE
              CALL FLAG_ERROR("Mesh region is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MESH_ADD")
    RETURN
999 IF(ASSOCIATED(NEW_COUPLED_MESHES)) DEALLOCATE(NEW_COUPLED_MESHES)
    CALL ERRORS("CONSTRAINT_MESH_ADD",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_ADD")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MESH_ADD

  !
  !================================================================================================================================
  !

  !>Finishes the creation of an constraint. \see OPENCMISS::CMISSConstraintCreateFinish
  SUBROUTINE CONSTRAINT_CREATE_FINISH(CONSTRAINT,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("CONSTRAINT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN
        CALL FLAG_ERROR("Constraint has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONSTRAINT%NUMBER_OF_COUPLED_MESHES<2) THEN
          LOCAL_ERROR="Invalid mesh coupling. Only "//TRIM(NUMBER_TO_VSTRING(CONSTRAINT%NUMBER_OF_COUPLED_MESHES,"*",ERR,ERROR))// &
            & " have been coupled. The number of coupled meshes must be >= 2."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
        CONSTRAINT%CONSTRAINT_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Constraint :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  User number = ",CONSTRAINT%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global number = ",CONSTRAINT%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",CONSTRAINT%LABEL,ERR,ERROR,*999)
      IF(ASSOCIATED(CONSTRAINT%CONSTRAINTS)) THEN
        IF(ASSOCIATED(CONSTRAINT%CONSTRAINTS%PARENT_REGION)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",CONSTRAINT%CONSTRAINTS% &
            & PARENT_REGION%USER_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",CONSTRAINT%CONSTRAINTS% &
            & PARENT_REGION%LABEL,ERR,ERROR,*999)        
        ELSE
          CALL FLAG_ERROR("Constraints parent region is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint constraints is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    CALL EXITS("CONSTRAINT_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CONSTRAINT_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of an constraint on a parent region. \see OPENCMISS::CMISSConstraintCreateStart
  SUBROUTINE CONSTRAINT_CREATE_START(USER_NUMBER,PARENT_REGION,CONSTRAINT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the constraint to create
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer to the parent region to create the constraint on.
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<On exit, a pointer to the created constraint. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: constraint_idx
    TYPE(CONSTRAINT_TYPE), POINTER :: NEW_CONSTRAINT
    TYPE(CONSTRAINT_PTR_TYPE), POINTER :: NEW_CONSTRAINTS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_CONSTRAINT)
    NULLIFY(NEW_CONSTRAINTS)
    
    CALL ENTERS("CONSTRAINT_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(PARENT_REGION)) THEN
      IF(ASSOCIATED(CONSTRAINT)) THEN
        CALL FLAG_ERROR("Constraint is already associated.",ERR,ERROR,*998)
      ELSE
        NULLIFY(CONSTRAINT)
        CALL CONSTRAINT_USER_NUMBER_FIND(USER_NUMBER,PARENT_REGION,CONSTRAINT,ERR,ERROR,*998)
        IF(ASSOCIATED(CONSTRAINT)) THEN
          LOCAL_ERROR="Constraint number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ELSE        
          NULLIFY(CONSTRAINT)
          !Allocate and set default constraint properties.
          CALL CONSTRAINT_INITIALISE(NEW_CONSTRAINT,ERR,ERROR,*999)
          NEW_CONSTRAINT%USER_NUMBER=USER_NUMBER
          NEW_CONSTRAINT%GLOBAL_NUMBER=PARENT_REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS+1
          LOCAL_STRING="Constraint_"//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
          NEW_CONSTRAINT%LABEL=CHAR(LOCAL_STRING)
          IF(ERR/=0) GOTO 999
          NEW_CONSTRAINT%CONSTRAINTS=>PARENT_REGION%CONSTRAINTS
          NEW_CONSTRAINT%PARENT_REGION=>PARENT_REGION
          !Add new initerface into list of constraints in the parent region
          ALLOCATE(NEW_CONSTRAINTS(PARENT_REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS+1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new constraints.",ERR,ERROR,*999)
          DO constraint_idx=1,PARENT_REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS
            NEW_CONSTRAINTS(constraint_idx)%PTR=>PARENT_REGION%CONSTRAINTS%CONSTRAINTS(constraint_idx)%PTR
          ENDDO !constraint_idx
          NEW_CONSTRAINTS(PARENT_REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS+1)%PTR=>NEW_CONSTRAINT
          IF(ASSOCIATED(PARENT_REGION%CONSTRAINTS%CONSTRAINTS)) DEALLOCATE(PARENT_REGION%CONSTRAINTS%CONSTRAINTS)
          PARENT_REGION%CONSTRAINTS%CONSTRAINTS=>NEW_CONSTRAINTS
          PARENT_REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS=PARENT_REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS+1
          CONSTRAINT=>NEW_CONSTRAINT
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Parent region is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("CONSTRAINT_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_CONSTRAINTS)) DEALLOCATE(NEW_CONSTRAINTS)
    CALL CONSTRAINT_FINALISE(CONSTRAINT,ERR,ERROR,*998)
998 CALL ERRORS("CONSTRAINT_CREATE_START",ERR,ERROR)
    CALL EXITS("CONSTRAINT_CREATE_START")
    RETURN 1
  END SUBROUTINE CONSTRAINT_CREATE_START

  !
  !================================================================================================================================
  !
  
  !>Returns the coordinate system of an constraint. \see OPENCMISS::CMISSConstraint_CoordinateSystemGet
  SUBROUTINE CONSTRAINT_COORDINATE_SYSTEM_GET(CONSTRAINT,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<On exit, the coordinate system for the specified constraint. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINT_COORDINATE_SYSTEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          CALL FLAG_ERROR("Coordinate system is already associated.",ERR,ERROR,*999)
        ELSE
          COORDINATE_SYSTEM=>CONSTRAINT%COORDINATE_SYSTEM
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_COORDINATE_SYSTEM_GET")
    RETURN
999 CALL ERRORS("CONSTRAINT_COORDINATE_SYSTEM_GET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_COORDINATE_SYSTEM_GET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_COORDINATE_SYSTEM_GET
  
  
  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of an constraint.  \see OPENCMISS::CMISSConstraint_CoordinateSystemSet
  SUBROUTINE CONSTRAINT_COORDINATE_SYSTEM_SET(CONSTRAINT,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to set the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_COORDINATE_SYSTEM_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN
        CALL FLAG_ERROR("Constraint has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
            CONSTRAINT%COORDINATE_SYSTEM=>COORDINATE_SYSTEM
          ELSE
            CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_COORDINATE_SYSTEM_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_COORDINATE_SYSTEM_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_COORDINATE_SYSTEM_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINT_COORDINATE_SYSTEM_SET
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the data points for a region. \see OPENCMISS::CMISSRegionDataPointsGet
  SUBROUTINE CONSTRAINT_DATA_POINTS_GET(CONSTRAINT,DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the region to get the data points for
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<On exit, a pointer to the data points for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("REGION_DATA_POINTS_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN 
        IF(ASSOCIATED(DATA_POINTS)) THEN
          CALL FLAG_ERROR("Data points is already associated.",ERR,ERROR,*998)
        ELSE
          DATA_POINTS=>CONSTRAINT%DATA_POINTS
          IF(.NOT.ASSOCIATED(DATA_POINTS)) CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONSTRAINT_DATA_POINTS_GET")
    RETURN
999 NULLIFY(DATA_POINTS)
998 CALL ERRORS("CONSTRAINT_DATA_POINTS_GET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_DATA_POINTS_GET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_DATA_POINTS_GET    

  !
  !================================================================================================================================
  !

  !>Destroys an constraint. \see OPENCMISS::CMISSConstraintDestroy
  SUBROUTINE CONSTRAINT_DESTROY(CONSTRAINT,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: constraint_idx,constraint_position
    TYPE(CONSTRAINT_PTR_TYPE), POINTER :: NEW_CONSTRAINTS(:)
    TYPE(CONSTRAINTS_TYPE), POINTER :: CONSTRAINTS
     
    NULLIFY(NEW_CONSTRAINTS)

    CALL ENTERS("CONSTRAINT_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      CONSTRAINTS=>CONSTRAINT%CONSTRAINTS
      IF(ASSOCIATED(CONSTRAINTS)) THEN
        constraint_position=CONSTRAINT%GLOBAL_NUMBER

        !Destroy all the constraint condition components
        CALL CONSTRAINT_FINALISE(CONSTRAINT,ERR,ERROR,*999)
        
        !Remove the constraint from the list of constraints
        IF(CONSTRAINTS%NUMBER_OF_CONSTRAINTS>1) THEN
          ALLOCATE(NEW_CONSTRAINTS(CONSTRAINTS%NUMBER_OF_CONSTRAINTS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new constraint conditions.",ERR,ERROR,*999)
          DO constraint_idx=1,CONSTRAINTS%NUMBER_OF_CONSTRAINTS
            IF(constraint_idx<constraint_position) THEN
              NEW_CONSTRAINTS(constraint_idx)%PTR=>CONSTRAINTS%CONSTRAINTS(constraint_idx)%PTR
            ELSE IF(constraint_idx>constraint_position) THEN
              CONSTRAINTS%CONSTRAINTS(constraint_idx)%PTR%GLOBAL_NUMBER=CONSTRAINTS%CONSTRAINTS(constraint_idx)%PTR%GLOBAL_NUMBER-1
              NEW_CONSTRAINTS(constraint_idx-1)%PTR=>CONSTRAINTS%CONSTRAINTS(constraint_idx)%PTR
            ENDIF
          ENDDO !constraint_idx
          IF(ASSOCIATED(CONSTRAINTS%CONSTRAINTS)) DEALLOCATE(CONSTRAINTS%CONSTRAINTS)
          CONSTRAINTS%CONSTRAINTS=>NEW_CONSTRAINTS
          CONSTRAINTS%NUMBER_OF_CONSTRAINTS=CONSTRAINTS%NUMBER_OF_CONSTRAINTS-1
        ELSE
          DEALLOCATE(CONSTRAINTS%CONSTRAINTS)
          CONSTRAINTS%NUMBER_OF_CONSTRAINTS=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Constraint constraints is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("CONSTRAINT_DESTROY")
    RETURN
999 CALL ERRORS("CONSTRAINT_DESTROY",ERR,ERROR)
    CALL EXITS("CONSTRAINT_DESTROY")
    RETURN 1
  END SUBROUTINE CONSTRAINT_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises an constraint and deallocates all memory.
  SUBROUTINE CONSTRAINT_FINALISE(CONSTRAINT,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("CONSTRAINT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(ASSOCIATED(CONSTRAINT%COUPLED_MESHES)) DEALLOCATE(CONSTRAINT%COUPLED_MESHES)
      IF(ASSOCIATED(CONSTRAINT%MESH_CONNECTIVITY)) &
         & CALL CONSTRAINT_MESH_CONNECTIVITY_FINALISE(CONSTRAINT%MESH_CONNECTIVITY,ERR,ERROR,*999)
      IF(ASSOCIATED(CONSTRAINT%pointsConnectivity)) &
        & CALL ConstraintPointsConnectivity_Finalise(CONSTRAINT%pointsConnectivity,ERR,ERROR,*999)
      IF(ASSOCIATED(CONSTRAINT%NODES)) CALL NODES_DESTROY(CONSTRAINT%NODES,ERR,ERROR,*999)
      CALL MESHES_FINALISE(CONSTRAINT%MESHES,ERR,ERROR,*999)
      CALL FIELDS_FINALISE(CONSTRAINT%FIELDS,ERR,ERROR,*999)
      CALL CONSTRAINT_CONDITIONS_FINALISE(CONSTRAINT%CONSTRAINT_CONDITIONS,ERR,ERROR,*999)
      DEALLOCATE(CONSTRAINT)
    ENDIF
    
    CALL EXITS("CONSTRAINT_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an constraint.
  SUBROUTINE CONSTRAINT_INITIALISE(CONSTRAINT,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("CONSTRAINT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      CALL FLAG_ERROR("Constraint is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(CONSTRAINT,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint.",ERR,ERROR,*999)
      CONSTRAINT%USER_NUMBER=0
      CONSTRAINT%GLOBAL_NUMBER=0
      CONSTRAINT%CONSTRAINT_FINISHED=.FALSE.
      CONSTRAINT%LABEL=""
      NULLIFY(CONSTRAINT%CONSTRAINTS)
      NULLIFY(CONSTRAINT%PARENT_REGION)
      CONSTRAINT%NUMBER_OF_COUPLED_MESHES=0
      NULLIFY(CONSTRAINT%COUPLED_MESHES)
      NULLIFY(CONSTRAINT%MESH_CONNECTIVITY)
      NULLIFY(CONSTRAINT%pointsConnectivity)
      NULLIFY(CONSTRAINT%NODES)
      NULLIFY(CONSTRAINT%MESHES)
      NULLIFY(CONSTRAINT%GENERATED_MESHES)
      NULLIFY(CONSTRAINT%FIELDS)
      NULLIFY(CONSTRAINT%CONSTRAINT_CONDITIONS)
      NULLIFY(CONSTRAINT%COORDINATE_SYSTEM)
      NULLIFY(CONSTRAINT%DATA_POINTS)
      CALL MESHES_INITIALISE(CONSTRAINT,ERR,ERROR,*999)
      CALL GENERATED_MESHES_INITIALISE(CONSTRAINT,ERR,ERROR,*999)
      CALL FIELDS_INITIALISE(CONSTRAINT,ERR,ERROR,*999)
      CALL CONSTRAINT_CONDITIONS_INITIALISE(CONSTRAINT,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the label of an constraint for a character label. \see OPENCMISS::CMISSConstraintLabelGet
  SUBROUTINE CONSTRAINT_LABEL_GET_C(CONSTRAINT,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On return the constraint label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: C_LENGTH,VS_LENGTH

    CALL ENTERS("CONSTRAINT_LABEL_GET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      C_LENGTH=LEN(LABEL)
      VS_LENGTH=LEN_TRIM(CONSTRAINT%LABEL)
      IF(C_LENGTH>VS_LENGTH) THEN
        LABEL=CHAR(CONSTRAINT%LABEL,VS_LENGTH)
      ELSE
        LABEL=CHAR(CONSTRAINT%LABEL,C_LENGTH)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_LABEL_GET_C")
    RETURN
999 CALL ERRORS("CONSTRAINT_LABEL_GET_C",ERR,ERROR)
    CALL EXITS("CONSTRAINT_LABEL_GET_C")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_LABEL_GET_C

   !
  !================================================================================================================================
  !

  !>Returns the label of an constraint for a varying string label. \see OPENCMISS::CMISSConstraintLabelGet
  SUBROUTINE CONSTRAINT_LABEL_GET_VS(CONSTRAINT,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On return the constraint label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_LABEL_GET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      !\todo The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
      LABEL=VAR_STR(CHAR(CONSTRAINT%LABEL))
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_LABEL_GET_VS")
    RETURN
999 CALL ERRORS("CONSTRAINT_LABEL_GET_VS",ERR,ERROR)
    CALL EXITS("CONSTRAINT_LABEL_GET_VS")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Sets the label of an constraint for a character label. \see OPENCMISS::CMISSConstraintLabelSet
  SUBROUTINE CONSTRAINT_LABEL_SET_C(CONSTRAINT,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_LABEL_SET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN
        CALL FLAG_ERROR("Constraint has been finished.",ERR,ERROR,*999)
      ELSE
        CONSTRAINT%LABEL=LABEL
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_LABEL_SET_C")
    RETURN
999 CALL ERRORS("CONSTRAINT_LABEL_SET_C",ERR,ERROR)
    CALL EXITS("CONSTRAINT_LABEL_SET_C")
    RETURN 1
  END SUBROUTINE CONSTRAINT_LABEL_SET_C

  !
  !================================================================================================================================
  !

  !>Sets the label of an constraint for a varying string label. \see OPENCMISS::CMISSConstraintLabelSet
  SUBROUTINE CONSTRAINT_LABEL_SET_VS(CONSTRAINT,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_LABEL_SET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN
        CALL FLAG_ERROR("Constraint has been finished.",ERR,ERROR,*999)
      ELSE
        CONSTRAINT%LABEL=LABEL
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_LABEL_SET_VS")
    RETURN
999 CALL ERRORS("CONSTRAINT_LABEL_SET_VS",ERR,ERROR)
    CALL EXITS("CONSTRAINT_LABEL_SET_VS")
    RETURN 1
  END SUBROUTINE CONSTRAINT_LABEL_SET_VS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the nodes for a constraint. \see OPENCMISS::CMISSConstraintNodesGet
  SUBROUTINE CONSTRAINT_NODES_GET(CONSTRAINT,NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to get the nodes for
    TYPE(NODES_TYPE), POINTER :: NODES !<On exit, a pointer to the nodes for the constraint. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONSTRAINT_NODES_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN 
        IF(ASSOCIATED(NODES)) THEN
          CALL FLAG_ERROR("Nodes is already associated.",ERR,ERROR,*998)
        ELSE
          NODES=>CONSTRAINT%NODES
          IF(.NOT.ASSOCIATED(NODES)) CALL FLAG_ERROR("Nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONSTRAINT_NODES_GET")
    RETURN
999 NULLIFY(NODES)
998 CALL ERRORS("CONSTRAINT_NODES_GET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_NODES_GET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_NODES_GET

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity for an constraint.
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_CREATE_FINISH(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: CONSTRAINT_MESH_CONNECTIVITY !<A pointer to the constraint meshes connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ConstraintElementIdx,CoupledMeshIdx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MESH_CONNECTIVITY)) THEN
      IF(CONSTRAINT_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Constraint meshes connectivity has already been finished.",ERR,ERROR,*999)
      ELSE
        !Check if connectivity from each constraint element to an appropriate element in the coupled meshes has been setup
        DO ConstraintElementIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS
          DO CoupledMeshIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
            IF (CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)% &
              & COUPLED_MESH_ELEMENT_NUMBER==0) THEN
              LOCAL_ERROR="The connectivity from constraint element " &
                //TRIM(NUMBER_TO_VSTRING(ConstraintElementIdx,"*",ERR,ERROR))//" to an element in coupled mesh " & 
                //TRIM(NUMBER_TO_VSTRING(CoupledMeshIdx,"*",ERR,ERROR))//" has not been defined."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !CoupledMeshIdx
        ENDDO !ConstraintElementIdx
        !Calculate line or face numbers for coupled mesh elements that are connected to the constraint mesh
        CALL CONSTRAINT_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*999)
        CONSTRAINT_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint meshes connectivity is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Initialises a meshes connectivity for an constraint.
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_CREATE_START(CONSTRAINT,MESH,CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to create the meshes connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: CONSTRAINT_MESH_CONNECTIVITY !<On return, a pointer to the created meshes connectivity
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(CONSTRAINT%CONSTRAINT_FINISHED) THEN
        IF(ASSOCIATED(CONSTRAINT%MESH_CONNECTIVITY)) THEN
          CALL FLAG_ERROR("The constraint already has a meshes connectivity associated.",ERR,ERROR,*999)
        ELSE
          !Initialise the meshes connectivity
          CALL CONSTRAINT_MESH_CONNECTIVITY_INITIALISE(CONSTRAINT,MESH,ERR,ERROR,*999)
          !Return the pointer
          CONSTRAINT_MESH_CONNECTIVITY=>CONSTRAINT%MESH_CONNECTIVITY
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_START")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_START",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_CREATE_START")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity and deallocates all memory
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_DESTROY(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: CONSTRAINT_MESH_CONNECTIVITY !<A pointer to the constraint meshes connectivity to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MESH_CONNECTIVITY)) THEN
      CALL CONSTRAINT_MESH_CONNECTIVITY_FINALISE(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Constraint meshes connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_DESTROY")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_DESTROY",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_DESTROY")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_DESTROY

  !
  !================================================================================================================================
  !

  !>Sets the constraint mesh connectivity basis
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET(CONSTRAINT_MESH_CONNECTIVITY,BASIS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: CONSTRAINT_MESH_CONNECTIVITY !<A pointer to constraint mesh connectivity to set the element number of elements for.
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    INTEGER(INTG) :: ConstraintElementIdx,CoupledMeshIdx,NumberOfConstraintElementNodes,NumberOfCoupledMeshXiDirections

    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MESH_CONNECTIVITY)) THEN
      IF(CONSTRAINT_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mesh connectivity already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(CONSTRAINT_MESH_CONNECTIVITY%BASIS)) THEN
          CALL FLAG_ERROR("Mesh connectivity basis already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(BASIS)) THEN
            CONSTRAINT_MESH_CONNECTIVITY%BASIS=>BASIS
            !Now that the mesh connectivity basis is set the number of constraint element nodes can be determined and now MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)%XI can be allocated
            !\todo NumberOfCoupledMeshXiDirections currently set to the number of constraint mesh xi directions + 1. Restructure ELEMENT_CONNECTIVITY type see below
            NumberOfCoupledMeshXiDirections=CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT_MESH%NUMBER_OF_DIMENSIONS+1
            NumberOfConstraintElementNodes=CONSTRAINT_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES
            DO ConstraintElementIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS
              DO CoupledMeshIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
                ELEMENT_CONNECTIVITY=>CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)
                !\todo Update mesh component index to look at the number of mesh components in each element. 
                !\todo Currently this defaults to the first mesh component ie %XI(NumberOfConstraintMeshXi,1,NumberOfConstraintElementNodes)). 
                !\todo The constraint mesh types will also need to be restructured.
                !eg. CONSTRAINT%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)%MESH_COMPONENT(MeshComponentIdx)%XI(NumberOfCoupledMeshXiDirections,NumberOfConstraintElementNodes) and adding appropriate initialize/finialize routines
                ALLOCATE(ELEMENT_CONNECTIVITY%XI(NumberOfCoupledMeshXiDirections,1,NumberOfConstraintElementNodes),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint element connectivity.",ERR,ERROR,*999)
                ELEMENT_CONNECTIVITY%XI=0.0_DP
              ENDDO !CoupledMeshIdx
            ENDDO !ConstraintElementIdx
          ELSE
            CALL FLAG_ERROR("Basis to set mesh connectivity not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mesh connectivity is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET

  !
  !================================================================================================================================
  !
    
  !>Sets the mapping from an xi position of a coupled mesh element to a node of an constraint mesh element
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_XI_SET(CONSTRAINT_MESH_CONNECTIVITY,CONSTRAINT_MESH_ELEMENT_NUMBER, &
    & COUPLED_MESH_INDEX,COUPLED_MESH_ELEMENT_NUMBER,CONSTRAINT_MESH_LOCAL_NODE_NUMBER,CONSTRAINT_MESH_COMPONENT_NUMBER,XI, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: CONSTRAINT_MESH_CONNECTIVITY !<A pointer to the constraint mesh connectivity for the constraint mesh
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_MESH_ELEMENT_NUMBER !<The constraint mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_INDEX !<The index of the coupled mesh at the constraint to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_ELEMENT_NUMBER !<The coupled mesh element to define the element xi connectivity from
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_MESH_LOCAL_NODE_NUMBER !<The constraint mesh node to assign the coupled mesh element xi to
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_MESH_COMPONENT_NUMBER !<The constraint mesh node's component to assign the coupled mesh element xi to
    REAL(DP), INTENT(IN) :: XI(:) !<XI(xi_idx). The xi value for the xi_idx'th xi direction in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: coupled_mesh_number_of_xi
    TYPE(CONSTRAINT_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_XI_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MESH_CONNECTIVITY)) THEN
      IF(CONSTRAINT_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mesh connectivity already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ALLOCATED(CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          IF((CONSTRAINT_MESH_ELEMENT_NUMBER>0).AND. &
            & (CONSTRAINT_MESH_ELEMENT_NUMBER<=CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS)) THEN
            IF((COUPLED_MESH_INDEX>0).AND.(COUPLED_MESH_INDEX<=CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES)) THEN
              IF((COUPLED_MESH_ELEMENT_NUMBER>0).AND.(COUPLED_MESH_ELEMENT_NUMBER<= &
                & CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT%COUPLED_MESHES(COUPLED_MESH_INDEX)%PTR%NUMBER_OF_ELEMENTS))THEN
                IF((CONSTRAINT_MESH_COMPONENT_NUMBER>0).AND. &
                  & (CONSTRAINT_MESH_COMPONENT_NUMBER<=CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT_MESH%NUMBER_OF_COMPONENTS)) THEN
                  IF((CONSTRAINT_MESH_LOCAL_NODE_NUMBER>0).AND.(CONSTRAINT_MESH_LOCAL_NODE_NUMBER<= &
                    & CONSTRAINT_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES))THEN
                    ELEMENT_CONNECTIVITY=>CONSTRAINT_MESH_CONNECTIVITY% &
                      & ELEMENT_CONNECTIVITY(CONSTRAINT_MESH_ELEMENT_NUMBER,COUPLED_MESH_INDEX)
                    IF(ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER==COUPLED_MESH_ELEMENT_NUMBER)THEN
                      coupled_mesh_number_of_xi = CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT%COUPLED_MESHES(COUPLED_MESH_INDEX)%PTR% & 
                        & TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(COUPLED_MESH_ELEMENT_NUMBER)%BASIS%NUMBER_OF_XI
                      IF(SIZE(XI)==coupled_mesh_number_of_xi) THEN
                        !\todo the ELEMENT_CONNECTIVITY%XI array needs to be restructured to efficiently utilize memory when coupling bodies with 2xi directions to bodies with 3xi directions using an constraint.
                        ELEMENT_CONNECTIVITY%XI(1:coupled_mesh_number_of_xi,CONSTRAINT_MESH_COMPONENT_NUMBER, &
                          & CONSTRAINT_MESH_LOCAL_NODE_NUMBER)=XI(1:coupled_mesh_number_of_xi)
                      ELSE
                        CALL FLAG_ERROR("The size of the xi array provided does not match the coupled mesh element's' number"// &
                          & " of xi.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Coupled mesh element number doesn't match that set to the constraint.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Constraint local node number is out of range.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Constraint component number is out of range.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Coupled mesh element number out of range.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint coupled mesh index number out of range.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint mesh element number out of range.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint elements connectivity array not allocated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mesh connectivity is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_XI_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_XI_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_XI_SET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_XI_SET

  !
  !================================================================================================================================
  !
  
  !>Sets the connectivity between an element in a coupled mesh to an element in the constraint mesh
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET(CONSTRAINT_MESH_CONNECTIVITY,CONSTRAINT_MESH_ELEMENT_NUMBER, &
      & COUPLED_MESH_INDEX,COUPLED_MESH_ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: CONSTRAINT_MESH_CONNECTIVITY !<A pointer to the constraint mesh connectivity for the constraint mesh
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_MESH_ELEMENT_NUMBER !<The constraint mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_INDEX !<The index of the coupled mesh at the constraint to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_ELEMENT_NUMBER !<The coupled mesh element to be connected to the constraint
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT_MESH_CONNECTIVITY)) THEN
      IF(CONSTRAINT_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FLAG_ERROR("Constraint mesh connectivity has already been finished.",ERR,ERROR,*999)
      ELSE
        IF((CONSTRAINT_MESH_ELEMENT_NUMBER>0).AND.(CONSTRAINT_MESH_ELEMENT_NUMBER<= &
          & CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS)) THEN
          IF((COUPLED_MESH_INDEX>0).AND.(COUPLED_MESH_INDEX<=CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES)) THEN
            IF (ALLOCATED(CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
              ELEMENT_CONNECTIVITY=>CONSTRAINT_MESH_CONNECTIVITY% &
                & ELEMENT_CONNECTIVITY(CONSTRAINT_MESH_ELEMENT_NUMBER,COUPLED_MESH_INDEX)
              ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER=COUPLED_MESH_ELEMENT_NUMBER
            ELSE
              CALL FLAG_ERROR("Constraint elements connectivity array not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint coupled mesh index number out of range.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mesh element number out of range.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint mesh connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_NUMBER_SET

  !
  !================================================================================================================================
  !
  
  !>Sets the connectivity between an element in a coupled mesh to an element in the constraint mesh
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_NODE_NUMBER_SET(NODES,CONSTRAINT_MESH_NODE_NUMBERS, &
      & FIRST_COUPLED_MESH_INDEX,FIRST_COUPLED_MESH_NODE_NUMBERS,SECOND_COUPLED_MESH_INDEX,SECOND_COUPLED_MESH_NODE_NUMBERS, &
      & ERR,ERROR,*)        
        
    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the constraint mesh connectivity for the constraint mesh
    INTEGER(INTG), INTENT(IN) :: CONSTRAINT_MESH_NODE_NUMBERS(:) !<The constraint mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: FIRST_COUPLED_MESH_INDEX,SECOND_COUPLED_MESH_INDEX !<The index of the coupled mesh at the constraint to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: FIRST_COUPLED_MESH_NODE_NUMBERS(:),SECOND_COUPLED_MESH_NODE_NUMBERS(:) !<The coupled mesh element to be connected to the constraint
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIndex
    
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_NODE_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%CONSTRAINT%MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        PRINT *, 'CHECK how to circumvent! constraint_routines.f90:1133'
        CALL FLAG_ERROR("Constraint mesh connectivity has already been finished.",ERR,ERROR,*999)
      ELSE
        !Default to two coupled meshes
        ALLOCATE(NODES%COUPLED_NODES(2,SIZE(CONSTRAINT_MESH_NODE_NUMBERS(:))))
        DO nodeIndex=1,SIZE(CONSTRAINT_MESH_NODE_NUMBERS(:))
          NODES%COUPLED_NODES(FIRST_COUPLED_MESH_INDEX,nodeIndex)=FIRST_COUPLED_MESH_NODE_NUMBERS(nodeIndex)
          NODES%COUPLED_NODES(SECOND_COUPLED_MESH_INDEX,nodeIndex)=SECOND_COUPLED_MESH_NODE_NUMBERS(nodeIndex)
        ENDDO
      ENDIF
    ELSE
      CALL FLAG_ERROR("Nodes are not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_NODE_NUMBER_SET")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_NODE_NUMBER_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_NODE_NUMBER_SET")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_NODE_NUMBER_SET

  !
  !================================================================================================================================
  !
  
  !>Calculate line or face numbers for coupled mesh elements that are connected to the constraint mesh
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE), POINTER :: CONSTRAINT_MESH_CONNECTIVITY !<A pointer to constraint mesh connectivity to calculate line or face numbers for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: XiDifference
    INTEGER(INTG) :: CoupledMeshIdx,ConstraintElementIdx,CoupledMeshXiIdx,NumberOfConstraintElementNodes,NumberOfConstraintMeshXi
    TYPE(CONSTRAINT_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE",ERR,ERROR,*999)
    
    DO CoupledMeshIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
      DO ConstraintElementIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS
        ELEMENT_CONNECTIVITY=>CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)
        NumberOfConstraintElementNodes=CONSTRAINT_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES
        NumberOfConstraintMeshXi=CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT_MESH%NUMBER_OF_DIMENSIONS
        SELECTCASE(NumberOfConstraintMeshXi)
        CASE(1) !Lines
          DO CoupledMeshXiIdx=1,CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT%COUPLED_MESHES(CoupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
            ! Calculate difference between first node and last node of an element
            XiDifference=ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,1)- &
              & ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfConstraintElementNodes)
            IF(ABS(XiDifference)<ZERO_TOLERANCE) THEN
              IF(ABS(ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfConstraintElementNodes))<ZERO_TOLERANCE) THEN
                ELEMENT_CONNECTIVITY%CONNECTED_LINE=3-(CoupledMeshXiIdx-1)*2
              ELSE
                ELEMENT_CONNECTIVITY%CONNECTED_LINE=4-(CoupledMeshXiIdx-1)*2
              ENDIF
            ENDIF
          ENDDO
        CASE(2) !Faces
          DO CoupledMeshXiIdx=1,CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT%COUPLED_MESHES(CoupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
            XiDifference=ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,1)- &
              & ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfConstraintElementNodes)
            IF(ABS(XiDifference)<ZERO_TOLERANCE) THEN
              IF(ABS(ELEMENT_CONNECTIVITY%XI(CoupledMeshXiIdx,1,NumberOfConstraintElementNodes))<ZERO_TOLERANCE) THEN
                ELEMENT_CONNECTIVITY%CONNECTED_FACE=(CoupledMeshXiIdx-1)*2+2
              ELSE
                ELEMENT_CONNECTIVITY%CONNECTED_FACE=(CoupledMeshXiIdx-1)*2+1
              ENDIF
            ENDIF
          ENDDO
        CASE DEFAULT 
          LOCAL_ERROR="Number of constraint mesh dimension of "//TRIM(NUMBER_TO_VSTRING(NumberOfConstraintMeshXi,"*",ERR,ERROR))// &
            & "is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDSELECT
      ENDDO
    ENDDO
    
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE")
    RETURN
    
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE")
    RETURN 1
  
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_CONNECTED_LINES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_FINALISE(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE) :: CONSTRAINT_MESH_CONNECTIVITY !<The constraint mesh connectivity to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_FINALISE",ERR,ERROR,*999)

    CALL CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_FINALISE(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*999)
    NULLIFY(CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT)
    NULLIFY(CONSTRAINT_MESH_CONNECTIVITY%CONSTRAINT_MESH)
    NULLIFY(CONSTRAINT_MESH_CONNECTIVITY%BASIS)
    CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS=0
    CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES=0
       
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_FINALISE")
    RETURN 1
    
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the constraint mesh connectivity.
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_INITIALISE(CONSTRAINT,MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to initialise the mesh connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(ASSOCIATED(CONSTRAINT%MESH_CONNECTIVITY)) THEN
        CALL FLAG_ERROR("Constraint mesh connectivity is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(CONSTRAINT%MESH_CONNECTIVITY,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint mesh connectivity.",ERR,ERROR,*999)
        CONSTRAINT%MESH_CONNECTIVITY%CONSTRAINT=>CONSTRAINT
        CONSTRAINT%MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED=.FALSE.
        CONSTRAINT%MESH_CONNECTIVITY%CONSTRAINT_MESH=>MESH
        NULLIFY(CONSTRAINT%MESH_CONNECTIVITY%BASIS)
        CALL CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_INITIALISE(CONSTRAINT,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculate the coupled mesh elements that are connected to each constraint element
  SUBROUTINE ConstraintPointsConnectivity_CoupledElementsCalculate(constraintPointsConnectivity,coupledMeshIdx,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity !<A pointer to the constraint points connectivity to calculate coupled elements for
    INTEGER(INTG), INTENT(IN) :: coupledMeshIdx !<The coupled mesh index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LIST_TYPE), POINTER :: elementNumbersList
    INTEGER(INTG) :: elementIdx,dataPointIdx,globalDataPointNumber,globalElementNumber,numberOfElementDataPoints, &
      & numberOfCoupledElements,coupledElementIdx
    INTEGER(INTG), ALLOCATABLE :: elementNumbers(:)
  
    CALL ENTERS("ConstraintPointsConnectivity_CoupledElementsCalculate",err,error,*999)

    IF(ASSOCIATED(constraintPointsConnectivity)) THEN
      IF(ALLOCATED(constraintPointsConnectivity%coupledElements)) THEN
        constraintPointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)=0; !Initialise the number of coupled mesh elements 
        DO elementIdx=1,SIZE(constraintPointsConnectivity%coupledElements,1) !Number of constraint elements
          numberOfElementDataPoints=constraintPointsConnectivity%constraintMesh%TOPOLOGY(1)%PTR%dataPoints% &
            & elementDataPoint(elementIdx)%numberOfProjectedData !Get the number of data points in constraint mesh element
          !Set up list
          NULLIFY(elementNumbersList)
          CALL LIST_CREATE_START(elementNumbersList,ERR,ERROR,*999)
          CALL LIST_DATA_TYPE_SET(elementNumbersList,LIST_INTG_TYPE,ERR,ERROR,*999)
          CALL LIST_INITIAL_SIZE_SET(elementNumbersList,numberOfElementDataPoints,ERR,ERROR,*999)
          CALL LIST_CREATE_FINISH(elementNumbersList,ERR,ERROR,*999)
          DO dataPointIdx=1,numberOfElementDataPoints
            globalDataPointNumber=constraintPointsConnectivity%constraintMesh%TOPOLOGY(1)%PTR%dataPoints% &
              & elementDataPoint(elementIdx)%dataIndices(dataPointIdx)%globalNumber
            globalElementNumber=constraintPointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
              & coupledMeshElementNumber
            CALL LIST_ITEM_ADD(elementNumbersList,globalElementNumber,ERR,ERROR,*999)
          ENDDO !dataPointIdx
          CALL LIST_REMOVE_DUPLICATES(elementNumbersList,ERR,ERROR,*999)
          CALL LIST_DETACH_AND_DESTROY(elementNumbersList,numberOfCoupledElements,elementNumbers, &
            & ERR,ERROR,*999)
          IF(ALLOCATED(constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers)) & 
            & DEALLOCATE(constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers) !for updating coupledElements after projection
          ALLOCATE(constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)% &
            & elementNumbers(numberOfCoupledElements),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate coupled mesh element numbers.",ERR,ERROR,*999)
          DO coupledElementIdx=1,numberOfCoupledElements
            constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers(coupledElementIdx)= &
              & elementNumbers(coupledElementIdx)
          ENDDO
          constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%numberOfCoupledElements=numberOfCoupledElements
          IF(ALLOCATED(elementNumbers)) DEALLOCATE(elementNumbers)
          IF(constraintPointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)<numberOfCoupledElements) THEN
            constraintPointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)=numberOfCoupledElements
          ENDIF
        ENDDO !elementIdx
      ELSE
        CALL FLAG_ERROR("Constraint points connectivity coupled elements is not allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_CoupledElementsCalculate")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_CoupledElementsCalculate",err,error)
    CALL EXITS("ConstraintPointsConnectivity_CoupledElementsCalculate")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_CoupledElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalise the points connectivity coupled mesh elements 
  SUBROUTINE ConstraintPointsConnectivity_CoupledElementsFinalise(constraintPointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity !<A pointer to the constraint points connectivity whose coupled mesh elements is to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,coupledMeshIdx
  
    CALL ENTERS("ConstraintPointsConnectivity_CoupledElementsFinalise",err,error,*999)

    IF(ASSOCIATED(constraintPointsConnectivity)) THEN
      IF(ALLOCATED(constraintPointsConnectivity%coupledElements)) THEN
        DO coupledMeshIdx=1,SIZE(constraintPointsConnectivity%coupledElements,2)
          DO elementIdx=1,SIZE(constraintPointsConnectivity%coupledElements,1)
            constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%numberOfCoupledElements=0
            IF(ALLOCATED(constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers)) &
              & DEALLOCATE(constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers)
          ENDDO
        ENDDO
        DEALLOCATE(constraintPointsConnectivity%coupledElements)
      END IF
      IF(ALLOCATED(constraintPointsConnectivity%maxNumberOfCoupledElements)) &
        & DEALLOCATE(constraintPointsConnectivity%maxNumberOfCoupledElements)
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_CoupledElementsFinalise")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_CoupledElementsFinalise",err,error)
    CALL EXITS("ConstraintPointsConnectivity_CoupledElementsFinalise")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_CoupledElementsFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialise the coupled mesh elements for points connectivity
  SUBROUTINE ConstraintPointsConnectivity_CoupledElementsInitialise(constraintPointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity !<A pointer to the constraint points connectivity to initialise coupled elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfConstraintElements,numberOfCoupledMeshes,coupledMeshIdx,elementIdx
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
  
    CALL ENTERS("ConstraintPointsConnectivity_CoupledElementsInitialise",err,error,*999)

     IF(ASSOCIATED(constraintPointsConnectivity)) THEN
       IF(ALLOCATED(constraintPointsConnectivity%coupledElements)) THEN
         CALL FLAG_ERROR("Constraint points connectivity coupled elements is already allocated.",err,error,*999)
       ELSE
         IF(ASSOCIATED(constraintPointsConnectivity%constraint)) THEN
           IF(ASSOCIATED(constraintPointsConnectivity%constraintMesh)) THEN
             numberOfConstraintElements=constraintPointsConnectivity%constraintMesh%NUMBER_OF_ELEMENTS
             numberOfCoupledMeshes=constraintPointsConnectivity%constraint%NUMBER_OF_COUPLED_MESHES
             ALLOCATE(constraintPointsConnectivity%coupledElements(numberOfConstraintElements,numberOfCoupledMeshes),STAT=ERR)
             IF(ERR/=0) CALL FLAG_ERROR("Could not allocate points connectivity coupled element.",ERR,ERROR,*999)
             DO coupledMeshIdx=1,numberOfCoupledMeshes
               DO elementIdx=1,numberOfConstraintElements
                 constraintPointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%numberOfCoupledElements=0
               ENDDO !elementIdx
             ENDDO !coupledMeshIdx
           ELSE
             CALL FLAG_ERROR("Constraint mesh is not associated.",err,error,*999)
           ENDIF
         ELSE
            CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
         ENDIF
       ENDIF
     ELSE
       CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
     ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_CoupledElementsInitialise")
    RETURN
999 CALL ConstraintPointsConnectivity_CoupledElementsFinalise(constraintPointsConnectivity,dummyErr,dummyError,*998) 
998 CALL ERRORS("ConstraintPointsConnectivity_CoupledElementsInitialise",err,error)
    CALL EXITS("ConstraintPointsConnectivity_CoupledElementsInitialise")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_CoupledElementsInitialise
  
  !
  !================================================================================================================================
  !

  !>Finish create constraint points connectivity
  SUBROUTINE ConstraintPointsConnectivity_CreateFinish(constraintPointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity !<A pointer to the constraint points connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONSTRAINT_TYPE), POINTER :: constraint
    INTEGER(INTG) :: coupledMeshIdx
  
    CALL ENTERS("ConstraintPointsConnectivity_CreateFinish",err,error,*999)

    IF(ASSOCIATED(constraintPointsConnectivity)) THEN
      IF(constraintPointsConnectivity%pointsConnectivityFinished) THEN
        CALL FLAG_ERROR("Constraint points connectivity has already been finished.",err,error,*999)
      ELSE
        CALL ConstraintPointsConnectivity_ReducedXiCalculate(constraintPointsConnectivity,err,error,*999) 
        constraint=>constraintPointsConnectivity%constraint
        IF(ASSOCIATED(constraint)) THEN
          CALL ConstraintPointsConnectivity_CoupledElementsInitialise(constraintPointsConnectivity,err,error,*999) 
          DO coupledMeshIdx=1,constraintPointsConnectivity%constraint%NUMBER_OF_COUPLED_MESHES
            CALL ConstraintPointsConnectivity_CoupledElementsCalculate(constraintPointsConnectivity,coupledMeshIdx,err,error,*999) 
          ENDDO  
        ELSE
          CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
        ENDIF
        constraintPointsConnectivity%pointsConnectivityFinished=.TRUE.   
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_CreateFinish")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_CreateFinish",err,error)
    CALL EXITS("ConstraintPointsConnectivity_CreateFinish")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start create constraint points connectivity
  SUBROUTINE ConstraintPointsConnectivity_CreateStart(constraint,constraintMesh,constraintPointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: constraint !<A pointer to the constraint to create the points connectivity for
    TYPE(MESH_TYPE), POINTER :: constraintMesh !<A pointer to the constraint mesh for which the points connectivity is created
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity !<On return, a pointer to the created points connectivity
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("ConstraintPointsConnectivity_CreateStart",err,error,*999)

    IF(ASSOCIATED(constraint)) THEN
      IF(constraint%CONSTRAINT_FINISHED) THEN
        IF(ASSOCIATED(constraint%pointsConnectivity)) THEN
          CALL FLAG_ERROR("The constraint already has a points connectivity associated.",err,error,*999)
        ELSE
          IF(ASSOCIATED(constraintMesh)) THEN
            CALL ConstraintPointsConnectivity_Initialise(constraint,constraintMesh,err,error,*999)
            !Return the pointer
            constraintPointsConnectivity=>constraint%pointsConnectivity
          ELSE
            CALL FLAG_ERROR("Constraint mesh is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_CreateStart")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_CreateStart",err,error)
    CALL EXITS("ConstraintPointsConnectivity_CreateStart")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_CreateStart
  
  !
  !================================================================================================================================
  !

  !>Reproject data points for points connectivity
  SUBROUTINE ConstraintPointsConnectivity_DataReprojection(constraint,constraintCondition,err,error,*) 

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: constraint !<A pointer to the constraint where data reprojection is performed
    TYPE(CONSTRAINT_CONDITION_TYPE), POINTER :: constraintCondition !<A pointer to the constraint where data reprojection is performed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity
    TYPE(FIELD_TYPE), POINTER :: dependentFieldFixed,dependentFieldProjection
    TYPE(DATA_POINTS_TYPE), POINTER :: dataPoints
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    INTEGER(INTG) :: fixedBodyIdx,projectionBodyIdx,dataPointIdx
    INTEGER(INTG) :: elementNumber,numberOfGeometricComponents
    INTEGER(INTG) :: coupledMeshFaceLineNumber,component
  
    CALL ENTERS("ConstraintPointsConnectivity_DataReprojection",err,error,*999)
    
    NULLIFY(interpolatedPoints)
    NULLIFY(interpolationParameters)
    fixedBodyIdx=2 !\todo: need to generalise
    projectionBodyIdx=1

    IF(ASSOCIATED(constraint)) THEN
      IF(ASSOCIATED(constraintCondition)) THEN
        constraintPointsConnectivity=>constraint%pointsConnectivity
        IF(ASSOCIATED(constraintPointsConnectivity)) THEN
          dataPoints=>constraint%DATA_POINTS
          IF(ASSOCIATED(dataPoints)) THEN

            !Evaluate data points positions
            dependentFieldFixed=>constraintCondition%DEPENDENT%FIELD_VARIABLES(fixedBodyIdx)%PTR%FIELD
            IF(ASSOCIATED(dependentFieldFixed)) THEN
              numberOfGeometricComponents=dependentFieldFixed%GEOMETRIC_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS
              CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(dependentFieldFixed,interpolationParameters,err,error,*999, &
                & FIELD_GEOMETRIC_COMPONENTS_TYPE)
              CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
                & FIELD_GEOMETRIC_COMPONENTS_TYPE)
              interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
              DO dataPointIdx=1,dataPoints%NUMBER_OF_DATA_POINTS
                elementNumber=constraintPointsConnectivity%pointsConnectivity(dataPointIdx,fixedBodyIdx)% &
                  & coupledMeshElementNumber
                coupledMeshFaceLineNumber=dependentFieldFixed%DECOMPOSITION%TOPOLOGY%ELEMENTS% &
                  & ELEMENTS(elementNumber)% &
                  & ELEMENT_FACES(constraintPointsConnectivity%pointsConnectivity(dataPointIdx,fixedBodyIdx)% &
                  & elementLineFaceNumber)
                CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,coupledMeshFaceLineNumber, &
                  & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,constraintPointsConnectivity%pointsConnectivity(dataPointIdx, &
                  & fixedBodyIdx)%reducedXi(:),interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
                DO component=1,numberOfGeometricComponents
                  dataPoints%DATA_POINTS(dataPointIdx)%position(component) = interpolatedPoint%VALUES(component,NO_PART_DERIV)
                ENDDO !component
              ENDDO !dataPointIdx
            ELSE
              CALL FLAG_ERROR("Fixed dependent field is not associated.",err,error,*999)
            ENDIF
            
            !Data reprojection and update points connectivity information with the projection results
            dataProjection=>dataPoints%DATA_PROJECTIONS(projectionBodyIdx+1)%PTR 
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ProjectedBodyDataProjectionLabel",ERR,ERROR,*999)
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,dataProjection%label,ERR,ERROR,*999)
            IF(ASSOCIATED(dataProjection)) THEN
              dependentFieldProjection=>constraintCondition%DEPENDENT%FIELD_VARIABLES(projectionBodyIdx)%PTR%FIELD
              IF(ASSOCIATED(dependentFieldProjection)) THEN
                !Projection the data points (with know spatial positions) on the projection dependent field 
                CALL DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE(dataProjection,dependentFieldProjection,err,error,*999)
                CALL ConstraintPointsConnectivity_UpdateFromProjection(ConstraintPointsConnectivity,dataProjection, &
                  & projectionBodyIdx,err,error,*999) 
              ELSE
                CALL FLAG_ERROR("Projection dependent field is not associated.",err,error,*999)
              ENDIF       
            ELSE
              CALL FLAG_ERROR("Constraint data projection is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint data points is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint condition is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_DataReprojection")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_DataReprojection",err,error)
    CALL EXITS("ConstraintPointsConnectivity_DataReprojection")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_DataReprojection
  
  !
  !================================================================================================================================
  !

  !>Destroy constraint points connectivity
  SUBROUTINE ConstraintPointsConnectivity_Destroy(constraintPointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity !<A pointer to constraint points connectivity to be destroyed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("ConstraintPointsConnectivity_Destroy",err,error,*999)

    IF(ASSOCIATED(constraintPointsConnectivity)) THEN
      CALL ConstraintPointsConnectivity_Finalise(constraintPointsConnectivity,err,error,*999) 
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_Destroy")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_Destroy",err,error)
    CALL EXITS("ConstraintPointsConnectivity_Destroy")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_Destroy
  
  !
  !================================================================================================================================
  !
  
  !>Gets the number of coupled mesh elements which are linked to a specific constraint element.
  SUBROUTINE ConstraintPointsConnectivity_ElementNumberGet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
      & meshComponentNumber,coupledMeshUserElementNumber,err,error,*)
      
    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to constraint points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the constraint to set element number for
    INTEGER(INTG), INTENT(OUT) :: coupledMeshUserElementNumber !<The coupled mesh element user number
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the constraint mesh that points connectivity is associated to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: constraintMesh
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists

    CALL ENTERS("ConstraintPointsConnectivity_ElementNumberGet",err,error,*999)
    
    IF(ASSOCIATED(pointsConnectivity)) THEN 
      CALL DATA_POINT_CHECK_EXISTS(pointsConnectivity%constraint%DATA_POINTS,dataPointUserNumber,dataPointExists, &
        & dataPointGlobalNumber,err,error,*999)
      IF(dataPointExists) THEN
        constraintMesh=>pointsConnectivity%constraintMesh
        coupledMeshUserElementNumber=pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)% &
          & coupledMeshElementNumber
      ELSE
        CALL FLAG_ERROR("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
          & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_ElementNumberGet")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_ElementNumberGet",err,error)
    CALL EXITS("ConstraintPointsConnectivity_ElementNumberGet")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_ElementNumberGet
  
  !
  !================================================================================================================================
  !
  
  !>Sets the number of coupled mesh elements which are linked to a specific constraint element.
  SUBROUTINE ConstraintPointsConnectivity_ElementNumberSet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
      & coupledMeshUserElementNumber,meshComponentNumber,err,error,*)
      
    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to constraint points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the constraint to set element number for
    INTEGER(INTG), INTENT(IN) :: coupledMeshUserElementNumber !<The coupled mesh element user number
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the constraint mesh that points connectivity is associated to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,elementMeshNumber
    LOGICAL :: dataPointExists,elementExists

    CALL ENTERS("ConstraintPointsConnectivity_ElementNumberSet",err,error,*999)
    
    IF(ASSOCIATED(pointsConnectivity)) THEN
      IF(pointsConnectivity%pointsConnectivityFinished) THEN
        CALL FLAG_ERROR("Constraint points connectivity has already been finished.",err,error,*999)
      ELSE
        CALL DATA_POINT_CHECK_EXISTS(pointsConnectivity%constraint%DATA_POINTS,dataPointUserNumber,dataPointExists, &
          & dataPointGlobalNumber,err,error,*999)
        IF(dataPointExists) THEN
          IF ((coupledMeshIndexNumber<=pointsConnectivity%constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS).OR. &
              & (coupledMeshIndexNumber>0)) THEN
            IF (ALLOCATED(pointsConnectivity%pointsConnectivity)) THEN
              CALL MeshTopologyElementCheckExists(pointsConnectivity%CONSTRAINT%COUPLED_MESHES(coupledMeshIndexNumber)%PTR, &
                & meshComponentNumber,coupledMeshUserElementNumber,elementExists,elementMeshNumber,err,error,*999) !Make sure user element exists       
              IF(elementExists) THEN
                pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)%coupledMeshElementNumber= &
                  & elementMeshNumber
              ELSE
                CALL FLAG_ERROR("Element with user number ("//TRIM(NUMBER_TO_VSTRING &
                  & (coupledMeshUserElementNumber,"*",err,error))//") does not exist.",err,error,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Constraint points connectivity array not allocated.",err,error,*999)
            END IF
          ELSE
            CALL FLAG_ERROR("Constraint coupled mesh index number out of range.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
              & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
        ENDIF
      ENDIF  
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_ElementNumberSet")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_ElementNumberSet",err,error)
    CALL EXITS("ConstraintPointsConnectivity_ElementNumberSet")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_ElementNumberSet
  
  !
  !================================================================================================================================
  !

  !>Finalise constraint points connectivity
  SUBROUTINE ConstraintPointsConnectivity_Finalise(constraintPointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: constraintPointsConnectivity !<A pointer to constraint points connectivity to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledMeshIdx,dataPointIdx
    
    CALL ENTERS("ConstraintPointsConnectivity_Finalise",err,error,*999)
    
    IF(ASSOCIATED(constraintPointsConnectivity)) THEN
      IF(ALLOCATED(constraintPointsConnectivity%PointsConnectivity)) THEN
        DO coupledMeshIdx=1,size(constraintPointsConnectivity%pointsConnectivity,2)
          DO dataPointIdx=1,size(constraintPointsConnectivity%pointsConnectivity,1) !Deallocate memory for each data point
            CALL ConstraintPointsConnectivity_PointFinalise(constraintPointsConnectivity%pointsConnectivity(dataPointIdx, &
              & coupledMeshIdx),err,error,*999) 
          ENDDO
        ENDDO
        DEALLOCATE(constraintPointsConnectivity%pointsConnectivity)  
      ENDIF
      CALL ConstraintPointsConnectivity_CoupledElementsFinalise(constraintPointsConnectivity,err,error,*999) 
      NULLIFY(constraintPointsConnectivity%constraint)
      NULLIFY(constraintPointsConnectivity%constraintMesh) 
      IF(ALLOCATED(constraintPointsConnectivity%coupledElements)) DEALLOCATE(constraintPointsConnectivity%coupledElements)
      IF(ALLOCATED(constraintPointsConnectivity%maxNumberOfCoupledElements))  &
        & DEALLOCATE(constraintPointsConnectivity%maxNumberOfCoupledElements)
      DEALLOCATE(constraintPointsConnectivity)
    ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_Finalise")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_Finalise",err,error)
    CALL EXITS("ConstraintPointsConnectivity_Finalise")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_Finalise
  
  !
  !================================================================================================================================
  !

  !>Finalise constraint point connectivity
  SUBROUTINE ConstraintPointsConnectivity_PointFinalise(constraintPointConnectivity,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointConnectivityType) :: constraintPointConnectivity !< An constraint point connectivity to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    CALL ENTERS("ConstraintPointsConnectivity_PointFinalise",err,error,*999)
    
    constraintPointConnectivity%coupledMeshElementNumber=0
    constraintPointConnectivity%elementLineFaceNumber=0
    IF(ALLOCATED(constraintPointConnectivity%xi)) DEALLOCATE(constraintPointConnectivity%xi)
    IF(ALLOCATED(constraintPointConnectivity%reducedXi)) DEALLOCATE(constraintPointConnectivity%reducedXi)
    
    CALL EXITS("ConstraintPointsConnectivity_PointFinalise")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_PointFinalise",err,error)
    CALL EXITS("ConstraintPointsConnectivity_PointFinalise")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_PointFinalise
  
  !
  !================================================================================================================================
  !

  !>Calculate full xi locations in points connectivity from reduced xi
  SUBROUTINE ConstraintPointsConnectivity_FullXiCalculate(ConstraintPointsConnectivity,coupledMeshIdx,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: ConstraintPointsConnectivity !<A pointer to the constraint points connectivity to finish creating
    INTEGER(INTG), INTENT(IN) :: coupledMeshIdx !<Coupled mesh index to calculate the full xi location for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,constraintMeshDimensions,coupledMeshDimensions
    TYPE(CONSTRAINT_TYPE), POINTER :: constraint
    TYPE(ConstraintPointConnectivityType), POINTER :: pointConnectivity
    TYPE(MESH_TYPE), POINTER :: constraintMesh,coupledMesh
    TYPE(VARYING_STRING) :: localError
  
    CALL ENTERS("ConstraintPointsConnectivity_FullXiCalculate",err,error,*999)

     IF(ASSOCIATED(ConstraintPointsConnectivity)) THEN
       constraintMesh=>ConstraintPointsConnectivity%constraintMesh
       IF(ASSOCIATED(constraintMesh)) THEN
         constraintMeshDimensions=ConstraintPointsConnectivity%constraintMesh%NUMBER_OF_DIMENSIONS
         constraint=>ConstraintPointsConnectivity%constraint
         IF(ASSOCIATED(constraint)) THEN
           coupledMesh=>constraint%COUPLED_MESHES(coupledMeshIdx)%PTR
           IF(ASSOCIATED(coupledMesh)) THEN
             coupledMeshDimensions=coupledMesh%NUMBER_OF_DIMENSIONS
             IF(constraintMeshDimensions==coupledMeshDimensions) THEN !e.g. If 1D-2D, 2D-3D coupling, constraint dimension is 1D and 2D respectively for 1st body
               DO dataPointIdx=1,constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS
                 ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%xi(:)= &
                   & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%reducedXi(:)
               ENDDO !dataPointIdx
             ELSE
               !Update full xi location from reduced xi and element face/line number
               SELECT CASE(coupledMeshDimensions)
               CASE(2)
                 DO dataPointIdx=1,constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS
                   pointConnectivity=>ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)
                   SELECT CASE(pointConnectivity%elementLineFaceNumber)
                   CASE(1)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=0.0_DP
                   CASE(2)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=1.0_DP
                   CASE(3)
                     pointConnectivity%xi(1)=0.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                   CASE(4)
                     pointConnectivity%xi(1)=1.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                   CASE DEFAULT
                     localError="Invalid local line number "// &
                       & TRIM(NUMBER_TO_VSTRING(pointConnectivity%elementLineFaceNumber,"*",err,error))//" ."
                     CALL FLAG_ERROR(localError,err,error,*999)
                   END SELECT
                 ENDDO !dataPointIdx
               CASE(3)
                 DO dataPointIdx=1,constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS
                   pointConnectivity=>ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)
                   SELECT CASE(pointConnectivity%elementLineFaceNumber)
                   CASE(1)
                     pointConnectivity%xi(1)=1.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(2)
                     pointConnectivity%xi(1)=0.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(3)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=1.0_DP
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(4)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=0.0_DP
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(5)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(2)
                     pointConnectivity%xi(3)=1.0_DP
                   CASE(6)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(2)
                     pointConnectivity%xi(3)=0.0_DP
                   CASE DEFAULT
                     localError="Invalid local face number "// &
                       & TRIM(NUMBER_TO_VSTRING(pointConnectivity%elementLineFaceNumber,"*",err,error))//" ."
                     CALL FLAG_ERROR(localError,err,error,*999)
                   END SELECT
                 ENDDO !dataPointIdx
               CASE DEFAULT
                 localError="Invalid coupled mesh dimension "// &
                   & TRIM(NUMBER_TO_VSTRING(pointConnectivity%elementLineFaceNumber,"*",err,error))//" ."
                 CALL FLAG_ERROR(localError,err,error,*999)
               END SELECT
             ENDIF
           ELSE
             CALL FLAG_ERROR("Coupled mesh is not associated.",err,error,*999)
           ENDIF
         ELSE
           CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
         ENDIF
       ELSE
         CALL FLAG_ERROR("Constraint mesh is not associated.",err,error,*999)
       ENDIF
     ELSE
       CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
     ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_FullXiCalculate")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_FullXiCalculate",err,error)
    CALL EXITS("ConstraintPointsConnectivity_FullXiCalculate")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_FullXiCalculate
  
    !
  !================================================================================================================================
  !

  !>Initialises the constraint mesh connectivity.
  SUBROUTINE ConstraintPointsConnectivity_Initialise(constraint,constraintMesh,err,error,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: constraint !<A pointer to the constraint to initialise points connectivity for
    TYPE(MESH_TYPE), POINTER :: constraintMesh !<A pointer to the constraint mesh to initialise points connectivity with
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,meshIdx,coupledMeshDimension,constraintMeshDimension
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
     
    CALL ENTERS("ConstraintPointsConnectivity_Initialise",err,error,*999)

    IF(ASSOCIATED(constraint)) THEN
      IF(ASSOCIATED(constraint%pointsConnectivity)) THEN
        CALL FLAG_ERROR("Constraint has already got points connectivity associated.",err,error,*998)
      ELSE
        !Initialise the poins connectivity
        ALLOCATE(constraint%pointsConnectivity,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint points connectivity.",err,error,*998)
        constraint%pointsConnectivity%constraint=>constraint
        constraint%pointsConnectivity%pointsConnectivityFinished=.FALSE.
        constraint%pointsConnectivity%constraintMesh=>constraintMesh
        constraintMeshDimension=constraintMesh%NUMBER_OF_DIMENSIONS
        IF(constraint%NUMBER_OF_COUPLED_MESHES>0) THEN 
          IF(ASSOCIATED(constraint%DATA_POINTS)) THEN
            IF(constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS>0) THEN
              ALLOCATE(constraint%pointsConnectivity%pointsConnectivity(constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS, &
                & constraint%NUMBER_OF_COUPLED_MESHES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint point connectivity.",err,error,*999)
              DO dataPointIdx=1,constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS
                DO meshIdx=1,constraint%NUMBER_OF_COUPLED_MESHES
                  coupledMeshDimension=constraint%COUPLED_MESHES(meshIdx)%PTR%NUMBER_OF_DIMENSIONS
                  CALL ConstraintPointsConnectivity_PointInitialise(constraint%pointsConnectivity%pointsConnectivity(dataPointIdx, &
                    & meshIdx),coupledMeshDimension,constraintMeshDimension,err,error,*999)
                ENDDO!meshIdx
              ENDDO!dataPointIdx
            ELSE
              CALL FLAG_ERROR("Number of constraint data points must be > 0.",err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint data points are not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Number of coupled meshes in the constraint must be > 0.",err,error,*999)
        ENDIF 
        ALLOCATE(constraint%pointsConnectivity%maxNumberOfCoupledElements(constraint%NUMBER_OF_COUPLED_MESHES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint max number of coupled mesh elements.",err,error,*999)
        constraint%pointsConnectivity%maxNumberOfCoupledElements=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("ConstraintPointsConnectivity_Initialise")
    RETURN
999 CALL ConstraintPointsConnectivity_Finalise(constraint%pointsConnectivity,dummyErr,dummyError,*998) 
998 CALL ERRORS("ConstraintPointsConnectivity_Initialise",err,error)
    CALL EXITS("ConstraintPointsConnectivity_Initialise")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_Initialise
  
  !
  !================================================================================================================================
  !

  !>Initialises the constraint mesh connectivity.
  SUBROUTINE ConstraintPointsConnectivity_PointInitialise(constraintPointConnectivity,coupledMeshDimension,constraintMeshDimension, &
      & err,error,*)

    !Argument variables
    TYPE(ConstraintPointConnectivityType) :: constraintPointConnectivity !<An constraint point connectivity to be initliased for
    INTEGER(INTG), INTENT(IN) :: coupledMeshDimension,constraintMeshDimension
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
     
    CALL ENTERS("ConstraintPointsConnectivity_PointInitialise",err,error,*999)

    constraintPointConnectivity%coupledMeshElementNumber=0
    constraintPointConnectivity%elementLineFaceNumber=0
    !Allocate memory for coupled mesh full and reduced xi location
    ALLOCATE(constraintPointConnectivity%xi(coupledMeshDimension),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint point connectivity full xi.",err,error,*999)
    constraintPointConnectivity%xi=0.0_DP
    ALLOCATE(constraintPointConnectivity%reducedXi(constraintMeshDimension),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint point connectivity reduced xi.",err,error,*999)
    constraintPointConnectivity%reducedXi=0.0_DP
       
    CALL EXITS("ConstraintPointsConnectivity_PointInitialise")
    RETURN
999 CALL ConstraintPointsConnectivity_PointFinalise(constraintPointConnectivity,dummyErr,dummyError,*998)
998 CALL ERRORS("ConstraintPointsConnectivity_PointInitialise",err,error)
    CALL EXITS("ConstraintPointsConnectivity_PointInitialise")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_PointInitialise
  
  !
  !================================================================================================================================
  !
    
  !>Gets the xi coordinate mapping between the data points in constraint and xi coordinates in a coupled region mesh
  SUBROUTINE ConstraintPointsConnectivity_PointXiGet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
    & xi,err,error,*)

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to constraint points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the constraint to set the number of elements for.
    REAL(DP), INTENT(OUT) :: xi(:) !<xi(xi_idx). The full xi location in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists
    
    CALL ENTERS("ConstraintPointsConnectivity_PointXiGet",err,error,*999)
    
    ! Preliminary error checks to verify user input information
    IF(ASSOCIATED(pointsConnectivity)) THEN
      IF (ALLOCATED(pointsConnectivity%pointsConnectivity)) THEN
        CALL DATA_POINT_CHECK_EXISTS(pointsConnectivity%CONSTRAINT%DATA_POINTS,dataPointUserNumber,dataPointExists, &
          & dataPointGlobalNumber,err,error,*999)
        IF(dataPointExists) THEN
          IF(SIZE(xi)>=SIZE(pointsConnectivity%pointsConnectivity(dataPointUserNumber,coupledMeshIndexNumber)%xi,1)) THEN
            xi=pointsConnectivity%pointsConnectivity(dataPointUserNumber,coupledMeshIndexNumber)%xi(:)
          ELSE
            CALL FLAG_ERROR("Input xi array size is smaller than points connectivity xi array.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
              & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint points connectivity array not allocated.",err,error,*999)
      ENDIF 
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF   

    CALL EXITS("ConstraintPointsConnectivity_PointXiGet")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_PointXiGet",err,error)
    CALL EXITS("ConstraintPointsConnectivity_PointXiGet")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_PointXiGet
  
  !
  !================================================================================================================================
  !
    
  !>Sets the xi coordinate mapping between the data points in constraint and xi coordinates in a coupled region mesh
  SUBROUTINE ConstraintPointsConnectivity_PointXiSet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
    & xi,err,error,*)

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to constraint points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the constraint to set the number of elements for.
    REAL(DP), INTENT(IN) :: xi(:) !<xi(xi_idx). The full xi location in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists
    
    CALL ENTERS("ConstraintPointsConnectivity_PointXiSet",err,error,*999)
    
    ! Preliminary error checks to verify user input information
    IF(ASSOCIATED(pointsConnectivity)) THEN
      IF(pointsConnectivity%pointsConnectivityFinished) THEN
        CALL FLAG_ERROR("Constraint mesh connectivity already been finished.",err,error,*999)
      ELSE
        IF (ALLOCATED(pointsConnectivity%pointsConnectivity)) THEN
          CALL DATA_POINT_CHECK_EXISTS(pointsConnectivity%CONSTRAINT%DATA_POINTS,dataPointUserNumber,dataPointExists, &
            & dataPointGlobalNumber,err,error,*999)
          IF(dataPointExists) THEN
            IF ((coupledMeshIndexNumber>pointsConnectivity%constraint%NUMBER_OF_COUPLED_MESHES).OR.(coupledMeshIndexNumber<0)) THEN
              CALL FLAG_ERROR("Constraint coupled mesh index number out of range.",err,error,*999)
            ELSE
              IF(SIZE(pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)%xi,1)== &
                  & SIZE(xi,1)) THEN
                pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)% &
                  & xi(:)=xi(:) 
              ELSE
                CALL FLAG_ERROR("Input xi dimension does not match full coupled mesh xi dimension.",err,error,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
                & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint points connectivity array not allocated.",err,error,*999)
        ENDIF 
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF   

    CALL EXITS("ConstraintPointsConnectivity_PointXiSet")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_PointXiSet",err,error)
    CALL EXITS("ConstraintPointsConnectivity_PointXiSet")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_PointXiSet
  
  !
  !================================================================================================================================
  !
  
  !>Update points connectivity with projection results
  SUBROUTINE ConstraintPointsConnectivity_UpdateFromProjection(ConstraintPointsConnectivity,dataProjection, &
      & coupledMeshIndex,err,error,*) 
  
    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: ConstraintPointsConnectivity !<A pointer to the constraint points connectivity to finish creating
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection !<The data projection that points connectivity update with
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndex !<The mesh index of the the points connectivity to be updated
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
    TYPE(DATA_PROJECTION_RESULT_TYPE), POINTER :: dataProjectionResult
    
    CALL ENTERS("ConstraintPointsConnectivity_UpdateFromProjection",err,error,*999)
    
    IF(ASSOCIATED(ConstraintPointsConnectivity)) THEN
      IF(ASSOCIATED(dataProjection)) THEN
        IF(dataProjection%DATA_PROJECTION_FINISHED) THEN
          WRITE(*,*) "ConstraintPointsConnectivity_UpdateFromProjection"
          DO dataPointIdx=1,SIZE(dataProjection%DATA_PROJECTION_RESULTS,1) !Update reduced xi location, projection element number and element face/line number with projection result
            dataProjectionResult=>dataProjection%DATA_PROJECTION_RESULTS(dataPointIdx)
            ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%reducedXi(:)= &
              & dataProjectionResult%XI
            ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%coupledMeshElementNumber= &
              & dataProjectionResult%ELEMENT_NUMBER
            IF(dataProjectionResult%ELEMENT_LINE_NUMBER/=0) THEN
              ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%elementLineFaceNumber= &
                & dataProjectionResult%ELEMENT_LINE_NUMBER
            ELSEIF(dataProjectionResult%ELEMENT_FACE_NUMBER/=0) THEN
              ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%elementLineFaceNumber= &
                & dataProjectionResult%ELEMENT_FACE_NUMBER
            ENDIF
          ENDDO
          CALL ConstraintPointsConnectivity_FullXiCalculate(ConstraintPointsConnectivity,coupledMeshIndex, &
            & err,error,*999) 
          !Update points connectivity coupledElement information
          CALL ConstraintPointsConnectivity_CoupledElementsCalculate(ConstraintPointsConnectivity,coupledMeshIndex,err,error,*999) 
        ELSE
          CALL FLAG_ERROR("Data projection is not finished.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
    ENDIF
     
    CALL EXITS("ConstraintPointsConnectivity_UpdateFromProjection")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_UpdateFromProjection",err,error)
    CALL EXITS("ConstraintPointsConnectivity_UpdateFromProjection")
    RETURN 1
    
  END SUBROUTINE ConstraintPointsConnectivity_UpdateFromProjection
  
  !
  !================================================================================================================================
  !

  !>Calculate reduced xi locations in points connectivity from full xi
  SUBROUTINE ConstraintPointsConnectivity_ReducedXiCalculate(ConstraintPointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(ConstraintPointsConnectivityType), POINTER :: ConstraintPointsConnectivity !<A pointer to the constraint points connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshIdx,dataPointIdx,xiIdx,constraintMeshDimensions,coupledMeshDimensions
    TYPE(CONSTRAINT_TYPE), POINTER :: constraint
  
    CALL ENTERS("ConstraintPointsConnectivity_ReducedXiCalculate",err,error,*999)

     IF(ASSOCIATED(ConstraintPointsConnectivity)) THEN
       IF(ConstraintPointsConnectivity%pointsConnectivityFinished) THEN
         CALL FLAG_ERROR("Constraint points connectivity has already been finished.",err,error,*999)
       ELSE
         constraint=>ConstraintPointsConnectivity%constraint
         constraintMeshDimensions=ConstraintPointsConnectivity%constraintMesh%NUMBER_OF_DIMENSIONS
         IF(ASSOCIATED(constraint)) THEN
           DO meshIdx=1,constraint%NUMBER_OF_COUPLED_MESHES
             coupledMeshDimensions=constraint%COUPLED_MESHES(meshIdx)%PTR%NUMBER_OF_DIMENSIONS
             IF(constraintMeshDimensions==coupledMeshDimensions) THEN !e.g. If 1D-2D, 2D-3D coupling, constraint dimension is 1D and 2D respectively for 1st body
               DO dataPointIdx=1,constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS
                 ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(:)= &
                   & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(:)
                 ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=1 !The only local face/line for the body with lower dimension 
               ENDDO !dataPointIdx
             ELSE
               SELECT CASE(coupledMeshDimensions)
               CASE(2)
                 DO dataPointIdx=1,constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS
                   DO xiIdx=1,coupledMeshDimensions
                     IF(ABS(ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)) &
                       & < ZERO_TOLERANCE) THEN
                       ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=4-(xiIdx-1)*2 !Calculate line number
                       EXIT
                     ELSEIF(ABS(ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)-1.0_DP) &
                       & < ZERO_TOLERANCE) THEN
                       ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=3-(xiIdx-1)*2 !Calculate line number
                       EXIT
                     ENDIF
                   ENDDO !xiIdx
                   ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi= &
                     & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)% &
                     & xi(OTHER_XI_DIRECTIONS2(xiIdx))  !Populate reducedXi 
                 ENDDO !dataPointIdx
               CASE(3)
                 DO dataPointIdx=1,constraint%DATA_POINTS%NUMBER_OF_DATA_POINTS
                   DO xiIdx=1,coupledMeshDimensions
                     IF(ABS(ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)) &
                       & < ZERO_TOLERANCE) THEN
                       ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=(xiIdx-1)*2+2 !Calculate face number
                       EXIT
                     ELSE IF(ABS(ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)-1.0_DP) &
                       & > ZERO_TOLERANCE) THEN
                       ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=(xiIdx-1)*2+1 !Calculate face number
                       EXIT
                     ENDIF
                   ENDDO !xiIdx
                   SELECT CASE(xiIdx) !Populate reducedXi 
                   CASE(1)
                     ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                      & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(2)
                     ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                      & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(3)
                   CASE(2)
                     ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                      & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(1)
                     ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                      & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(3)
                   CASE(3)
                     ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                      & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(1)
                     ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                      & ConstraintPointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(2)
                   END SELECT
                 ENDDO !dataPointIdx
               CASE DEFAULT
                 ! Do nothing
               END SELECT
             ENDIF
           ENDDO !meshIdx
         ELSE
           CALL FLAG_ERROR("Constraint is not associated.",err,error,*999)
         ENDIF
       ENDIF
     ELSE
       CALL FLAG_ERROR("Constraint points connectivity is not associated.",err,error,*999)
     ENDIF
    
    CALL EXITS("ConstraintPointsConnectivity_ReducedXiCalculate")
    RETURN
999 CALL ERRORS("ConstraintPointsConnectivity_ReducedXiCalculate",err,error)
    CALL EXITS("ConstraintPointsConnectivity_ReducedXiCalculate")
    RETURN 1
  END SUBROUTINE ConstraintPointsConnectivity_ReducedXiCalculate

  !
  !================================================================================================================================
  !

  !>Finds and returns in CONSTRAINT a pointer to the constraint identified by USER_NUMBER in the given PARENT_REGION. If no constraint with that USER_NUMBER exists CONSTRAINT is left nullified.
  SUBROUTINE CONSTRAINT_USER_NUMBER_FIND(USER_NUMBER,PARENT_REGION,CONSTRAINT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<The parent region to find the constraint in    
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<On return a pointer to the constraint with the given user number. If no constraint with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: constraint_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINT_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(PARENT_REGION)) THEN
      IF(ASSOCIATED(CONSTRAINT)) THEN
        CALL FLAG_ERROR("Constraint is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(CONSTRAINT)
        IF(ASSOCIATED(PARENT_REGION%CONSTRAINTS)) THEN
          constraint_idx=1
          DO WHILE(constraint_idx<=PARENT_REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS.AND..NOT.ASSOCIATED(CONSTRAINT))
            IF(PARENT_REGION%CONSTRAINTS%CONSTRAINTS(constraint_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
              CONSTRAINT=>PARENT_REGION%CONSTRAINTS%CONSTRAINTS(constraint_idx)%PTR
            ELSE
              constraint_idx=constraint_idx+1
            ENDIF
          ENDDO
        ELSE
          LOCAL_ERROR="The constraints on parent region number "// &
            & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//" are not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Parent region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("CONSTRAINT_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("CONSTRAINT_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE CONSTRAINT_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises constraints and deallocates all memory.
  SUBROUTINE CONSTRAINTS_FINALISE(CONSTRAINTS,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINTS_TYPE), POINTER :: CONSTRAINTS !<A pointer to the constraints to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT
     
    CALL ENTERS("CONSTRAINTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINTS)) THEN
      DO WHILE(CONSTRAINTS%NUMBER_OF_CONSTRAINTS>0)
        CONSTRAINT=>CONSTRAINTS%CONSTRAINTS(1)%PTR
        CALL CONSTRAINT_DESTROY(CONSTRAINT,ERR,ERROR,*999)
      ENDDO
      IF(ASSOCIATED(CONSTRAINTS%CONSTRAINTS)) DEALLOCATE(CONSTRAINTS%CONSTRAINTS)
      DEALLOCATE(CONSTRAINTS)
    ENDIF
    
    CALL EXITS("CONSTRAINTS_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINTS_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises constraints for a region.
  SUBROUTINE CONSTRAINTS_INITIALISE(REGION,ERR,ERROR,*) 

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the constraints for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    CALL ENTERS("CONSTRAINTS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%CONSTRAINTS)) THEN
        LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has constraints associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        ALLOCATE(REGION%CONSTRAINTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region constraints.",ERR,ERROR,*999)
        REGION%CONSTRAINTS%PARENT_REGION=>REGION
        REGION%CONSTRAINTS%NUMBER_OF_CONSTRAINTS=0
        NULLIFY(REGION%CONSTRAINTS%CONSTRAINTS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINTS_INITIALISE")
    RETURN
999 CALL CONSTRAINTS_FINALISE(REGION%CONSTRAINTS,ERR,ERROR,*998)
998 CALL ERRORS("CONSTRAINTS_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the constraint element connectivity.
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_INITIALISE(CONSTRAINT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONSTRAINT_TYPE), POINTER :: CONSTRAINT !<A pointer to the constraint to initialise the mesh connectivity for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ConstraintElementIdx,CoupledMeshIdx
    TYPE(CONSTRAINT_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
     
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONSTRAINT)) THEN
      IF(ASSOCIATED(CONSTRAINT%MESH_CONNECTIVITY)) THEN
        IF(ALLOCATED(CONSTRAINT%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          CALL FLAG_ERROR("Constraint mesh element connectivity is already allocated.",ERR,ERROR,*999)
        ELSE
          IF(CONSTRAINT%NUMBER_OF_COUPLED_MESHES>0) THEN
            IF(CONSTRAINT%MESH_CONNECTIVITY%CONSTRAINT_MESH%NUMBER_OF_ELEMENTS>0) THEN
              CONSTRAINT%MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS=CONSTRAINT%MESHES%MESHES(1)%PTR%NUMBER_OF_ELEMENTS
              CONSTRAINT%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES=CONSTRAINT%NUMBER_OF_COUPLED_MESHES
              ALLOCATE(CONSTRAINT%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(CONSTRAINT%MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS, &
                & CONSTRAINT%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate constraint element connectivity.",ERR,ERROR,*999)
              DO ConstraintElementIdx=1,CONSTRAINT%MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS
                DO CoupledMeshIdx=1,CONSTRAINT%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
                  ELEMENT_CONNECTIVITY=>CONSTRAINT%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)
                  ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER=0
                  ELEMENT_CONNECTIVITY%CONNECTED_FACE=0
                  ELEMENT_CONNECTIVITY%CONNECTED_LINE=0
                  ! Note that ELEMENT_CONNECTIVITY%XI(NumberOfCoupledMeshXiDirections,1,NumberOfConstraintElementNodes) is allocated after the basis for the mesh connectivity has been specified in CONSTRAINT_MESH_CONNECTIVITY_BASIS_SET where the number of NumberOfConstraintElementNodes can be determined.
                  !\todo see corresponding todo in regarding updating the structure of ELEMENT_CONNECTIVITY%XI
                ENDDO !CoupledMeshIdx
              ENDDO !ConstraintElementIdx
            ELSE
              CALL FLAG_ERROR("Constraint coupled meshes are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Constraint coupled meshes are not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Constraint mesh connectivity is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Constraint is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises an constraint element connectivity and deallocates all memory.
  SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_FINALISE(CONSTRAINT_MESH_CONNECTIVITY,ERR,ERROR,*) 

    !Argument variables
    TYPE(CONSTRAINT_MESH_CONNECTIVITY_TYPE) :: CONSTRAINT_MESH_CONNECTIVITY !<The constraint element connectivity to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ConstraintElementIdx,CoupledMeshIdx
     
    CALL ENTERS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_FINALISE",ERR,ERROR,*999)

    DO ConstraintElementIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_CONSTRAINT_ELEMENTS  
      DO CoupledMeshIdx=1,CONSTRAINT_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
        IF(ALLOCATED(CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)%COUPLED_MESH_ELEMENT_NUMBER=0
          IF(ALLOCATED(CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)%XI)) THEN
            DEALLOCATE(CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ConstraintElementIdx,CoupledMeshIdx)%XI)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Constraint mesh connectivity element connectivity is being deallocated before allocation.", &
            & ERR,ERROR,*999)
        ENDIF
      ENDDO !ConstraintElementIdx
    ENDDO !CoupledMeshIdx

    DEALLOCATE(CONSTRAINT_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)
    
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_FINALISE")
    RETURN 1
  END SUBROUTINE CONSTRAINT_MESH_CONNECTIVITY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

END MODULE CONSTRAINT_ROUTINES
