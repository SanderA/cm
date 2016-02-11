!> \file
!> \author Chris Bradley
!> \brief This module handles all mesh (node and element) routines.
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

!> This module handles all mesh (node and element) routines.
MODULE MESH_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS_MPI
  USE CMISS_PARMETIS
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE DOMAIN_MAPPINGS
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE LISTS
  USE MPI
  USE NODE_ROUTINES
  USE SORTING
  USE STRINGS
  USE TREES
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup MESH_ROUTINES_DecompositionTypes MESH_ROUTINES::DecompositionTypes
  !> \brief The Decomposition types parameters
  !> \see MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_ALL_TYPE=1 !<The decomposition contains all elements. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_CALCULATED_TYPE=2 !<The element decomposition is calculated by graph partitioning. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_USER_DEFINED_TYPE=3 !<The user will set the element decomposition. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Starts the process of creating a mesh
  INTERFACE MESH_CREATE_START
    MODULE PROCEDURE MESH_CREATE_START_INTERFACE
    MODULE PROCEDURE MESH_CREATE_START_REGION
  END INTERFACE !MESH_CREATE_START

  !>Initialises the meshes for a region or interface.
  INTERFACE MESHES_INITIALISE
    MODULE PROCEDURE MESHES_INITIALISE_INTERFACE
    MODULE PROCEDURE MESHES_INITIALISE_REGION
  END INTERFACE !MESHES_INITIALISE

  INTERFACE MESH_USER_NUMBER_FIND
    MODULE PROCEDURE MESH_USER_NUMBER_FIND_INTERFACE
    MODULE PROCEDURE MESH_USER_NUMBER_FIND_REGION
  END INTERFACE !MESH_USER_NUMBER_FIND

  INTERFACE MeshTopologyElementCheckExists
    MODULE PROCEDURE MeshTopologyElementCheckExistsMesh
    MODULE PROCEDURE MeshTopologyElementCheckExistsMeshElements
  END INTERFACE MeshTopologyElementCheckExists
  
  PUBLIC DECOMPOSITION_ALL_TYPE,DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE

  PUBLIC DECOMPOSITIONS_INITIALISE,DECOMPOSITIONS_FINALISE

  PUBLIC DECOMPOSITION_CREATE_START,DECOMPOSITION_CREATE_FINISH

  PUBLIC DECOMPOSITION_DESTROY

  PUBLIC DECOMPOSITION_MESH_COMPONENT_NUMBER_GET,DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  PUBLIC DECOMPOSITION_NUMBER_OF_DOMAINS_GET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET
  
  PUBLIC DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS,DecompositionTopology_DataPointCheckExists
  
  PUBLIC DecompositionTopology_DataProjectionCalculate

  PUBLIC DecompositionTopology_ElementDataPointLocalNumberGet

  PUBLIC DecompositionTopology_ElementDataPointUserNumberGet

  PUBLIC DecompositionTopology_NumberOfElementDataPointsGet
  
  PUBLIC DECOMPOSITION_TYPE_GET,DECOMPOSITION_TYPE_SET
  
  PUBLIC DECOMPOSITION_USER_NUMBER_FIND,DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION
  
  PUBLIC DECOMPOSITION_CALCULATE_LINES_SET,DECOMPOSITION_CALCULATE_FACES_SET

  PUBLIC DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS

  PUBLIC DomainTopology_ElementBasisGet
  
  PUBLIC MeshTopologyElementCheckExists
  
  PUBLIC MESH_CREATE_START,MESH_CREATE_FINISH

  PUBLIC MESH_DESTROY
  
  PUBLIC MESH_NUMBER_OF_COMPONENTS_GET,MESH_NUMBER_OF_COMPONENTS_SET

  PUBLIC MESH_NUMBER_OF_ELEMENTS_GET,MESH_NUMBER_OF_ELEMENTS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_CREATE_START,MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH

  PUBLIC MESH_TOPOLOGY_ELEMENTS_DESTROY

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET,MeshElements_ElementNodeVersionSet

  PUBLIC MESH_TOPOLOGY_ELEMENTS_GET

  PUBLIC MeshElements_ElementUserNumberGet,MeshElements_ElementUserNumberSet
  
  PUBLIC MeshTopologyElementsUserNumbersAllSet

  PUBLIC MeshTopologyDataPointsCalculateProjection

  PUBLIC MESH_USER_NUMBER_FIND,MESH_USER_NUMBER_TO_MESH

  PUBLIC MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  PUBLIC MESH_EMBEDDING_CREATE,MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  PUBLIC MESH_EMBEDDING_SET_GAUSS_POINT_DATA

  PUBLIC MESHES_INITIALISE,MESHES_FINALISE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalises a decomposition adjacent element information and deallocates all memory
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
    IF(ALLOCATED(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises a decomposition adjacent element information.
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition on a given mesh. \see OPENCMISS::CMISSDecompositionCreateFinish
  SUBROUTINE DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(MESH_TYPE), POINTER :: MESH

    ENTERS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(DECOMPOSITION)) CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
      !Calculate which elements belong to which domain
      CALL Decomposition_ElementDomainsCalculate(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the topology information for this decomposition
      CALL DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the domain for this computational node
      CALL DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Calculate the domain mappings and topology from the new element domains.
      CALL Domain_Calculate(decomposition,err,error,*999) 
      !Calculate the decomposition topology
      CALL DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      DECOMPOSITION%DECOMPOSITION_FINISHED=.TRUE.

    IF(DIAGNOSTICS1) THEN
      MESH=>DECOMPOSITION%MESH
      IF(.NOT.ASSOCIATED(MESH)) CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Mesh = ",MESH%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of decompositions = ", &
        & MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS,ERR,ERROR,*999)
      DO decomposition_no=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Decomposition number = ",decomposition_no,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
          & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%GLOBAL_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
          & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%USER_NUMBER,ERR,ERROR,*999)
      ENDDO !decomposition_no
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a given mesh. \see OPENCMISS::CMISSDecompositionCreateStart
  SUBROUTINE DECOMPOSITION_CREATE_START(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition 
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to decompose
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On return a pointer to the created decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: NEW_DECOMPOSITION
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITION)
    NULLIFY(NEW_DECOMPOSITIONS)

    ENTERS("DECOMPOSITION_CREATE_START",ERR,ERROR,*999)

    NULLIFY(DECOMPOSITION)
    
    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
          IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
            CALL DECOMPOSITION_USER_NUMBER_FIND(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
            IF(ASSOCIATED(DECOMPOSITION)) THEN
              LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
                & " has already been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              !\todo Split this into an initialise and create start.
              ALLOCATE(NEW_DECOMPOSITION,STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decomposition.",ERR,ERROR,*999)
              !Set default decomposition properties
              NEW_DECOMPOSITION%GLOBAL_NUMBER=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1
              NEW_DECOMPOSITION%USER_NUMBER=USER_NUMBER
              NEW_DECOMPOSITION%DECOMPOSITION_FINISHED=.FALSE.
              NEW_DECOMPOSITION%CALCULATE_LINES=.TRUE. !Default
              NEW_DECOMPOSITION%CALCULATE_FACES=.FALSE. !Default
              NEW_DECOMPOSITION%DECOMPOSITIONS=>MESH%DECOMPOSITIONS
              NEW_DECOMPOSITION%MESH=>MESH
              !By default, the process of decompostion was done on the first mesh components. But the decomposition is the same for all mesh components, since the decomposition is element-based.
              NEW_DECOMPOSITION%MESH_COMPONENT_NUMBER=1
              !Default decomposition is all the mesh with one domain.
              NEW_DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
              NEW_DECOMPOSITION%NUMBER_OF_DOMAINS=1          
              NULLIFY(NEW_DECOMPOSITION%elementDomains)
              CALL Decomposition_ElementDomainsInitialise(NEW_DECOMPOSITION%elementDomains,err,error,*999)
              !Nullify the domain
              NULLIFY(NEW_DECOMPOSITION%DOMAIN)
              !Nullify the topology
              NULLIFY(NEW_DECOMPOSITION%TOPOLOGY)
              !\todo change this to use move alloc.
              !Add new decomposition into list of decompositions on the mesh
              ALLOCATE(NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
              DO decomposition_no=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
                NEW_DECOMPOSITIONS(decomposition_no)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR
              ENDDO !decomposition_no
              NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1)%PTR=>NEW_DECOMPOSITION
              IF(ASSOCIATED(MESH%DECOMPOSITIONS%DECOMPOSITIONS)) DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
              MESH%DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
              MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1        
              DECOMPOSITION=>NEW_DECOMPOSITION
            ENDIF
          ELSE
            LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
              & " are not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITION)) THEN
      CALL Decomposition_ElementDomainsFinalise(NEW_DECOMPOSITION%elementDomains,err,error,*999) 
      DEALLOCATE(NEW_DECOMPOSITION)
    ENDIF
    IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    NULLIFY(DECOMPOSITION)
    ERRORSEXITS("DECOMPOSITION_CREATE_START",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a user number and deallocates all memory. \see OPENCMISS::CMISSDecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY_NUMBER(USER_NUMBER,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to destroy.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh containing the decomposition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITIONS)

    ENTERS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN

        !Find the decomposition identified by the user number
        FOUND=.FALSE.
        decomposition_position=0
        DO WHILE(decomposition_position<MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.FOUND)
          decomposition_position=decomposition_position+1
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          DECOMPOSITION=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR

          !Destroy all the decomposition components          
          CALL Decomposition_ElementDomainsFinalise(DECOMPOSITION%elementDomains,err,error,*999)
          CALL DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
          CALL DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
          
          DEALLOCATE(DECOMPOSITION)

          !Remove the decomposition from the list of decompositions
          IF(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
            ALLOCATE(NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
            DO decomposition_idx=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
              IF(decomposition_idx<decomposition_position) THEN
                NEW_DECOMPOSITIONS(decomposition_idx)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ELSE IF(decomposition_idx>decomposition_position) THEN
                MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                  & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
                NEW_DECOMPOSITIONS(decomposition_idx-1)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ENDIF
            ENDDO !decomposition_idx
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
          ELSE
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    ERRORSEXITS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a pointer and deallocates all memory. \see OPENCMISS::CMISSDecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: DECOMPOSITIONS
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITIONS)

    ENTERS("DECOMPOSITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      DECOMPOSITIONS=>DECOMPOSITION%DECOMPOSITIONS
      IF(ASSOCIATED(DECOMPOSITIONS)) THEN
        decomposition_position=DECOMPOSITION%GLOBAL_NUMBER

        !Destroy all the decomposition components          
        CALL Decomposition_ElementDomainsFinalise(DECOMPOSITION%elementDomains,err,error,*999) 
        CALL DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
        CALL DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
        
        DEALLOCATE(DECOMPOSITION)
        
        !Remove the decomposition from the list of decompositions
        IF(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
          ALLOCATE(NEW_DECOMPOSITIONS(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
          DO decomposition_idx=1,DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
            IF(decomposition_idx<decomposition_position) THEN
              NEW_DECOMPOSITIONS(decomposition_idx)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
            ELSE IF(decomposition_idx>decomposition_position) THEN
              DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                & DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
              NEW_DECOMPOSITIONS(decomposition_idx-1)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
            ENDIF
          ENDDO !decomposition_idx
          DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
          DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
          DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
        ELSE
          DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
          DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        ENDIF
      ELSE
        CALL FlagError("Decomposition decompositions is not associated.",ERR,ERROR,*999)
       ENDIF
    ELSE
      CALL FlagError("Decompositions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    ERRORSEXITS("DECOMPOSITION_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the (local) partition for a decomposition of a mesh.
  SUBROUTINE Decomposition_ElementDomainsFinalise(elementDomains,err,error,*)

    !Argument variables
    TYPE(DecompositionElementDomainsType), POINTER :: elementDomains !<A pointer to the decomposition element domains to finalise for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_ElementDomainsFinalise",err,error,*999)

    IF(ASSOCIATED(elementDomains)) THEN 
      IF(ALLOCATED(elementDomains%offsets)) DEALLOCATE(elementDomains%offsets)
      IF(ALLOCATED(elementDomains%domains)) DEALLOCATE(elementDomains%domains)
      DEALLOCATE(elementDomains)
    END IF

    EXITS("Decomposition_ElementDomainsFinalise")
    RETURN
999 ERRORSEXITS("Decomposition_ElementDomainsFinalise",err,error)
    RETURN 1
  END SUBROUTINE Decomposition_ElementDomainsFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the (local) partition for a decomposition of a mesh.
  SUBROUTINE Decomposition_ElementDomainsInitialise(elementDomains,err,error,*)

    !Argument variables
    TYPE(DecompositionElementDomainsType), POINTER :: elementDomains !<A pointer to the decomposition element domains to initialise for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_ElementDomainsInitialise",err,error,*999)

    IF(ASSOCIATED(elementDomains)) CALL FlagError("Decomposition element domains is already associated.",err,error,*999)
    ALLOCATE(elementDomains,stat=err) 
    IF(err/=0) CALL FlagError("Coult not allocate decomposition element domains.",err,error,*999)

    EXITS("Decomposition_ElementDomainsInitialise")
    RETURN
999 ERRORSEXITS("Decomposition_ElementDomainsInitialise",err,error)
    RETURN 1
  END SUBROUTINE Decomposition_ElementDomainsInitialise
  
  !
  !================================================================================================================================
  !

  !>Calculates the (local) partition for a decomposition of a mesh.
  SUBROUTINE Decomposition_ElementDomainsCalculate(decomposition,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<A pointer to the decomposition to calculate the partition for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,myDomain,numberOfDomains,domainIdx,meshComponentNumber,numberOfCutEdges,userNodeNumber, &
      & elementNodeIdx,numberOfLocalElements,localElementNumber,globalElementNumber,numberOfAdjacentDomains, &
      & adjacentElementIdx,adjacentElementNumber,numberOfBoundaryElements,boundaryElementIdx,start,finish,mid,elementOwner, &
      & receiveBufferIdx,otherDomain,adjacentDomainIdx,adjacentDomainNumber,insertStatus,numberOfElementDomains, &
      & elementDomainIdx,dummyErr,MPI_IERROR
    INTEGER(INTG) :: elementWeight(1),adjacentWeight(1),weightFlag,numberFlag,numberOfConstraints,numberOfCommonNodes, &
      & parmetisOptions(3)
    INTEGER(INTG), ALLOCATABLE :: elementNodeOffsets(:),elementNodeNumbers(:),domainOffsets(:),partition(:), &
      & boundaryElements(:),adjacentDomainNumbers(:),sendCounts(:),receiveCounts(:),requests(:),statuses(:,:)
    INTEGER(INTG), POINTER :: adjacentElementOffsets(:),adjacentElements(:)
    REAL(DP) :: ubvec(1)
    REAL(DP), ALLOCATABLE :: tpwgts(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DecompositionElementDomainsType), POINTER :: elementDomains
    TYPE(INTEGER_INTG_ALLOC_TYPE), ALLOCATABLE :: sendBuffers(:),receiveBuffers(:)
    TYPE(LIST_TYPE), POINTER :: boundaryElementsList,elementDomainsList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: adjacentDomainsLists(:),newAdjacentDomainsLists(:),sendBufferLists(:)
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(TREE_TYPE), POINTER :: ghostOwnerTree
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Decomposition_ElementDomainsCalculate",err,error,*999)

    IF(ASSOCIATED(decomposition)) THEN
      mesh=>decomposition%mesh
      IF(ASSOCIATED(mesh)) THEN
        IF(ASSOCIATED(mesh%topology)) THEN
          SELECT CASE(decomposition%DECOMPOSITION_TYPE)
          CASE(DECOMPOSITION_ALL_TYPE)
            !Do nothing. Decomposition checked below.
          CASE(DECOMPOSITION_CALCULATED_TYPE)
            meshComponentNumber=decomposition%MESH_COMPONENT_NUMBER
            myDomain=COMPUTATIONAL_NODE_NUMBER_GET(err,error)+1
            IF(err/=0) GOTO 999
            numberOfDomains=decomposition%NUMBER_OF_DOMAINS
            elementDomains=>decomposition%elementDomains
            IF(.NOT.ASSOCIATED(elementDomains)) CALL FlagError("Decomposition element domains is not associated.",err,error,*999)
            !Calculate the general decomposition
            IF(numberOfDomains==1) THEN
              ALLOCATE(elementDomains%offsets(mesh%NUMBER_OF_ELEMENTS+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate element domain offsets.",err,error,*999)
              elementDomains%offsets=[(elementIdx,elementIdx=1,mesh%NUMBER_OF_ELEMENTS+1)]
              ALLOCATE(elementDomains%domains(mesh%NUMBER_OF_ELEMENTS),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate element domains.",err,error,*999)
              elementDomains%domains=1
            ELSE
              meshElements=>mesh%topology(meshComponentNumber)%ptr%elements

              ALLOCATE(domainOffsets(numberOfDomains+1),stat=err)
              IF(ERR/=0) CALL FlagError("Could not allocate domain elements offset.",err,error,*999)
              !Uniform distribution of elements over the domains.
              !Here NUMBER_OF_ELEMENTS is the global number of elements.
              domainOffsets=[((mesh%NUMBER_OF_ELEMENTS*(domainIdx-1))/numberOfDomains+1,domainIdx=1,numberOfDomains+1)]
              numberOfLocalElements=domainOffsets(myDomain+1)-domainOffsets(myDomain)

              ALLOCATE(elementNodeOffsets(numberOfLocalElements+1),stat=err)
              IF(ERR/=0) CALL FlagError("Could not allocate element nodes offset.",err,error,*999)
              elementNodeOffsets(1)=1
              DO elementIdx=1,numberOfLocalElements
                globalElementNumber=domainOffsets(myDomain)+elementIdx-1
                basis=>meshElements%elements(globalElementNumber)%basis
                elementNodeOffsets(elementIdx+1)=elementNodeOffsets(elementIdx)+basis%NUMBER_OF_NODES
              ENDDO !elementIdx

              ALLOCATE(elementNodeNumbers(elementNodeOffsets(numberOfLocalElements+1)-1),stat=err)
              IF(ERR/=0) CALL FlagError("Could not allocate element node numbers.",err,error,*999)

              !Question: can parmetis work with user numbers if the node numbering is non-contiguous, e.g. with 4 nodes the
              !numbering is 1,999,2000,9999? If so, there should be either a constraint on the node numbering, or we should
              !calculate the dual graph ourselves if possible.
              
              DO elementIdx=1,numberOfLocalElements
                globalElementNumber=domainOffsets(myDomain)+elementIdx-1
                basis=>meshElements%elements(globalElementNumber)%basis
                DO elementNodeIdx=1,basis%NUMBER_OF_NODES
                  userNodeNumber=meshElements%elements(globalElementNumber)%USER_ELEMENT_NODES(elementNodeIdx)
                  elementNodeNumbers(elementNodeOffsets(elementIdx)+elementNodeIdx-1)=userNodeNumber
                END DO !elementNodeIdx
              END DO !elementIdx
              
              !Set up ParMETIS variables
              weightFlag=0 !No weights
              elementWeight(1)=1 !Isn't used due to weight flag
              adjacentWeight(1)=1 !Isn't used due to weight flag
              numberFlag=1 !Use Fortran numbering
              numberOfConstraints=1
              !Has to be 1, or we can't calculate all the ghost elements.
              numberOfCommonNodes=1
              
              ALLOCATE(tpwgts(numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate tpwgts.",err,error,*999)
              tpwgts=1.0_DP/numberOfDomains
              ubvec=1.05_DP
              parmetisOptions(1)=1 !If zero, defaults are used, otherwise next two values are used
              parmetisOptions(2)=7 !Level of information to output
              parmetisOptions(3)=CMISS_RANDOM_SEEDS(1) !Seed for random number generator

              !Calculate the dual mesh. This gives us information about which elements are adjacent (share nodes with) with which
              !other elements. From this we can also calculate which elements are ghost in the initial (uniform) partitioning.
              NULLIFY(adjacentElementOffsets)
              NULLIFY(adjacentElements)

              CALL ParMETIS_Mesh2Dual(domainOffsets,elementNodeOffsets,elementNodeNumbers,numberFlag,numberOfCommonNodes, &
                & adjacentElementOffsets,adjacentElements,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,err,error,*999)

              IF(ALLOCATED(elementNodeOffsets)) DEALLOCATE(elementNodeOffsets)
              IF(ALLOCATED(elementNodeNumbers)) DEALLOCATE(elementNodeNumbers)

              !Call ParMETIS to calculate the partitioning of the dual mesh.
              ALLOCATE(partition(numberOfLocalElements),stat=err)
              IF(ERR/=0) CALL FlagError("Could not allocate element domain.",err,error,*999)
              CALL ParMETIS_PartKway(domainOffsets,adjacentElementOffsets,adjacentElements,elementWeight,adjacentWeight,weightFlag, &
                & numberFlag,numberOfConstraints,numberOfDomains,tpwgts,ubvec,parmetisOptions,numberOfCutEdges,partition, &
                & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,err,error,*999)

              IF(ALLOCATED(tpwgts)) DEALLOCATE(tpwgts)
              
              !Find the boundary elements of the initial partitioning.
              NULLIFY(boundaryElementsList)
              CALL List_CreateStart(boundaryElementsList,err,error,*999)
              CALL List_DataTypeSet(boundaryElementsList,LIST_INTG_TYPE,err,error,*999)
              CALL List_InitialSizeSet(boundaryElementsList,numberOfLocalElements,err,error,*999)
              CALL List_CreateFinish(boundaryElementsList,err,error,*999)

              DO elementIdx=1,numberOfLocalElements
                DO adjacentElementIdx=adjacentElementOffsets(elementIdx),adjacentElementOffsets(elementIdx+1)-1
                  adjacentElementNumber=adjacentElements(adjacentElementIdx)
                  IF(adjacentElementNumber<domainOffsets(myDomain).OR.domainOffsets(myDomain+1)<=adjacentElementNumber) THEN
                    !The adjacent element is a ghost element, so we have found a boundary element.
                    CALL List_ItemAdd(boundaryElementsList,elementIdx,err,error,*999)
                    EXIT
                  END IF
                END DO !adjacentElementIdx
              END DO !elementIdx
              CALL List_RemoveDuplicates(boundaryElementsList,err,error,*999)
              CALL List_DetachAndDestroy(boundaryElementsList,numberOfBoundaryElements,boundaryElements,err,error,*999)

              ALLOCATE(adjacentDomainsLists(numberOfBoundaryElements),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate adjacent domains lists.",err,error,*999)

              !Calculate the adjacent domains of the boundary elements in the initial (uniform) partitioning of the mesh.
              DO boundaryElementIdx=1,numberOfBoundaryElements
                localElementNumber=boundaryElements(boundaryElementIdx)
                NULLIFY(adjacentDomainsLists(boundaryElementIdx)%ptr)
                CALL List_CreateStart(adjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                CALL List_DataTypeSet(adjacentDomainsLists(boundaryElementIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                CALL List_InitialSizeSet(adjacentDomainsLists(boundaryElementIdx)%ptr,9,err,error,*999)
                CALL List_CreateFinish(adjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                DO adjacentElementIdx=adjacentElementOffsets(localElementNumber),adjacentElementOffsets(localElementNumber+1)-1
                  adjacentElementNumber=adjacentElements(adjacentElementIdx)
                  IF(adjacentElementNumber<domainOffsets(myDomain).OR.domainOffsets(myDomain+1)<=adjacentElementNumber) THEN
                    !The adjacent element is a ghost element. Determine which process it owns by searching for the last domain index
                    !for which the domain offset is not greater than the adjacent element number (just in case for when there would
                    !be zero elements in a domain). The searching is done with the binary search algorithm.
                    start=LBOUND(domainOffsets,1)
                    finish=UBOUND(domainOffsets,1)
                    DO WHILE(finish-start>0)
                      mid=start+(finish-start+1)/2 !Round up
                      IF(domainOffsets(mid)>adjacentElementNumber) THEN
                        finish=mid-1
                      ELSE
                        start=mid
                      END IF
                    END DO
                    elementOwner=start
                    CALL List_ItemAdd(adjacentDomainsLists(boundaryElementIdx)%ptr,elementOwner,err,error,*999)
                  END IF
                END DO !adjacentElementIdx
              END DO !elementIdx

              !Send the global numbers and the new domain owners of the boundary elements to the adjacent domains.
              ALLOCATE(sendBuffers(numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate send buffers.",err,error,*999)
              ALLOCATE(receiveBuffers(numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate receive buffers.",err,error,*999)

              ALLOCATE(receiveCounts(numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate receive counts.",err,error,*999)
              receiveCounts=0
              ALLOCATE(sendCounts(numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate send counts.",err,error,*999)
              sendCounts=0
              
              ALLOCATE(sendBufferLists(numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate send buffer lists.",err,error,*999)
              DO domainIdx=1,numberOfDomains
                NULLIFY(sendBufferLists(domainIdx)%ptr)
                CALL List_CreateStart(sendBufferLists(domainIdx)%ptr,err,error,*999)
                CALL List_DataTypeSet(sendBufferLists(domainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                CALL List_InitialSizeSet(sendBufferLists(domainIdx)%ptr,2*numberOfBoundaryElements,err,error,*999)
                CALL List_CreateFinish(sendBufferLists(domainIdx)%ptr,err,error,*999)
              END DO !domainIdx

              DO boundaryElementIdx=1,numberOfBoundaryElements
                localElementNumber=boundaryElements(boundaryElementIdx)
                globalElementNumber=domainOffsets(myDomain)+localElementNumber-1
                CALL List_RemoveDuplicates(adjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                CALL List_NumberOfItemsGet(adjacentDomainsLists(boundaryElementIdx)%ptr,numberOfAdjacentDomains,err,error,*999)
                CALL List_DetachAndDestroy(adjacentDomainsLists(boundaryElementIdx)%ptr, &
                  & numberOfAdjacentDomains,adjacentDomainNumbers,err,error,*999)
                DO adjacentDomainIdx=1,numberOfAdjacentDomains
                  adjacentDomainNumber=adjacentDomainNumbers(adjacentDomainIdx)
                  CALL List_ItemAdd(sendBufferLists(adjacentDomainNumber)%ptr,globalElementNumber,err,error,*999)
                  CALL List_ItemAdd(sendBufferLists(adjacentDomainNumber)%ptr,partition(localElementNumber),err,error,*999)
                END DO
                IF(ALLOCATED(adjacentDomainNumbers)) DEALLOCATE(adjacentDomainNumbers)
              END DO
              IF(ALLOCATED(adjacentDomainsLists)) DEALLOCATE(adjacentDomainsLists)
              IF(ALLOCATED(boundaryElements)) DEALLOCATE(boundaryElements)

              DO domainIdx=1,numberOfDomains
                CALL List_DetachAndDestroy(sendBufferLists(domainIdx)%ptr,sendCounts(domainIdx), &
                  & sendBuffers(domainIdx)%array,err,error,*999)
              END DO
              IF(ALLOCATED(sendBufferLists)) DEALLOCATE(sendBufferLists)

              !Collect the receive counts from all domains.
              CALL MPI_ALLTOALL(sendCounts,1,MPI_INTEGER,receiveCounts,1,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_ALLTOALL",MPI_IERROR,err,error,*999)

              !Allocate all the receive buffers.
              DO domainIdx=1,numberOfDomains
                ALLOCATE(receiveBuffers(domainIdx)%array(receiveCounts(domainIdx)),stat=err)
                IF(err/=0) CALL FlagError("Could not allocate receive buffer.",err,error,*999)
                receiveBuffers(domainIdx)%array=0
              END DO

              ALLOCATE(requests(2*numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate mpi receive requests.",err,error,*999)
              requests=MPI_REQUEST_NULL

              !Get all the receive buffers from the other domains.
              DO domainIdx=1,numberOfDomains
                IF(receiveCounts(domainIdx)>0) THEN
                  CALL MPI_IRECV(receiveBuffers(domainIdx)%array,receiveCounts(domainIdx),MPI_INTEGER,domainIdx-1, &
                    & MPI_GHOST_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,requests(domainIdx),MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,err,error,*999)
                END IF
              END DO
              DO domainIdx=1,numberOfDomains
                IF(sendCounts(domainIdx)>0) THEN
                  CALL MPI_ISEND(sendBuffers(domainIdx)%array,sendCounts(domainIdx),MPI_INTEGER,domainIdx-1, &
                    & MPI_GHOST_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,requests(numberOfDomains+domainIdx),MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,err,error,*999)
                END IF
              END DO

              !Create a tree with as key global element number and value the new owner.
              NULLIFY(ghostOwnerTree)
              CALL Tree_CreateStart(ghostOwnerTree,err,error,*999)
              CALL Tree_InsertTypeSet(ghostOwnerTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
              CALL Tree_CreateFinish(ghostOwnerTree,err,error,*999)

              ALLOCATE(statuses(MPI_STATUS_SIZE,2*numberOfDomains),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate mpi statuses.",err,error,*999)

              CALL MPI_WAITALL(2*numberOfDomains,requests,statuses,MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,err,error,*999)
              
              IF(ALLOCATED(statuses)) DEALLOCATE(statuses)
              IF(ALLOCATED(requests)) DEALLOCATE(requests)

              !Unpack the received data.
              DO domainIdx=1,numberOfDomains
                DO receiveBufferIdx=1,receiveCounts(domainIdx),2
                  globalElementNumber=receiveBuffers(domainIdx)%array(receiveBufferIdx)
                  elementOwner=receiveBuffers(domainIdx)%array(receiveBufferIdx+1)
                  CALL Tree_ItemInsert(ghostOwnerTree,globalElementNumber,elementOwner,insertStatus,err,error,*999)
                END DO !receiveBufferIdx
              END DO !domainIdx

              IF(ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
              IF(ALLOCATED(receiveCounts)) DEALLOCATE(receiveCounts)

              DO domainIdx=1,numberOfDomains
                IF(ALLOCATED(receiveBuffers(domainIdx)%array)) DEALLOCATE(receiveBuffers(domainIdx)%array)
                IF(ALLOCATED(sendBuffers(domainIdx)%array)) DEALLOCATE(sendBuffers(domainIdx)%array)
              END DO
              IF(ALLOCATED(receiveBuffers)) DEALLOCATE(receiveBuffers)
              IF(ALLOCATED(sendBuffers)) DEALLOCATE(sendBuffers)

              !Find out the boundary elements in the new partitioning
              NULLIFY(boundaryElementsList)
              CALL List_CreateStart(boundaryElementsList,err,error,*999)
              CALL List_DataTypeSet(boundaryElementsList,LIST_INTG_TYPE,err,error,*999)
              CALL List_InitialSizeSet(boundaryElementsList,numberOfLocalElements,err,error,*999)
              CALL List_CreateFinish(boundaryElementsList,err,error,*999)

              DO elementIdx=1,numberOfLocalElements
                DO adjacentElementIdx=adjacentElementOffsets(elementIdx),adjacentElementOffsets(elementIdx+1)-1
                  adjacentElementNumber=adjacentElements(adjacentElementIdx)
                  IF(adjacentElementNumber<domainOffsets(myDomain).OR.domainOffsets(myDomain+1)-1<adjacentElementNumber) THEN
                    !The adjacent element is a ghost element in the initial partitioning.
                    !Find ghost element owner from tree
                    NULLIFY(treeNode)
                    CALL Tree_Search(ghostOwnerTree,adjacentElementNumber,treeNode,err,error,*999)
                    CALL Tree_NodeValueGet(ghostOwnerTree,treeNode,otherDomain,err,error,*999)
                  ELSE
                    !The adjacent element is local in the initial partitioning.
                    otherDomain=partition(adjacentElementNumber-domainOffsets(myDomain)+1)
                  END IF
                  IF(partition(elementIdx)/=otherDomain) THEN
                    !This local element is a boundary element for the new domain.
                    CALL List_ItemAdd(boundaryElementsList,elementIdx,err,error,*999)
                    EXIT
                  END IF
                END DO !adjacentElementIdx
              END DO !elementIdx
              CALL List_RemoveDuplicates(boundaryElementsList,err,error,*999)
              CALL List_DetachAndDestroy(boundaryElementsList,numberOfBoundaryElements,boundaryElements,err,error,*999)

              ALLOCATE(newAdjacentDomainsLists(numberOfBoundaryElements),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate new adjacent domains lists.",err,error,*999)

              DO boundaryElementIdx=1,numberOfBoundaryElements
                localElementNumber=boundaryElements(boundaryElementIdx)
                NULLIFY(newAdjacentDomainsLists(boundaryElementIdx)%ptr)
                CALL List_CreateStart(newAdjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                CALL List_DataTypeSet(newAdjacentDomainsLists(boundaryElementIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                CALL List_InitialSizeSet(newAdjacentDomainsLists(boundaryElementIdx)%ptr,9,err,error,*999)
                CALL List_CreateFinish(newAdjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                DO adjacentElementIdx=adjacentElementOffsets(localElementNumber),adjacentElementOffsets(localElementNumber+1)-1
                  adjacentElementNumber=adjacentElements(adjacentElementIdx)
                  IF(adjacentElementNumber<domainOffsets(myDomain).OR.domainOffsets(myDomain+1)<=adjacentElementNumber) THEN
                    !The adjacent element is a ghost element in the initial partitioning.
                    !Find ghost element owner from tree
                    NULLIFY(treeNode)
                    CALL Tree_Search(ghostOwnerTree,adjacentElementNumber,treeNode,err,error,*999)
                    CALL Tree_NodeValueGet(ghostOwnerTree,treeNode,otherDomain,err,error,*999)
                  ELSE
                    !The adjacent element is local in the initial partitioning.
                    otherDomain=partition(adjacentElementNumber-domainOffsets(myDomain)+1)
                  END IF
                  IF(partition(localElementNumber)/=otherDomain) THEN
                    !This local element is a boundary element for the new domain.
                    CALL List_ItemAdd(newAdjacentDomainsLists(boundaryElementIdx)%ptr,otherDomain,err,error,*999)
                  END IF
                END DO !adjacentElementIdx
              END DO !localElementNumber
              IF(ASSOCIATED(adjacentElementOffsets)) DEALLOCATE(adjacentElementOffsets)
              IF(ASSOCIATED(adjacentElements)) DEALLOCATE(adjacentElements)
              CALL Tree_Destroy(ghostOwnerTree,err,error,*999)

              NULLIFY(elementDomainsList)
              CALL List_CreateStart(elementDomainsList,err,error,*999)
              CALL List_DataTypeSet(elementDomainsList,LIST_INTG_TYPE,err,error,*999)
              CALL List_InitialSizeSet(elementDomainsList,numberOfLocalElements+numberOfBoundaryElements*9,err,error,*999)
              CALL List_CreateFinish(elementDomainsList,err,error,*999)

              ALLOCATE(elementDomains%offsets(numberOfLocalElements+1),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate element domain offsets.",err,error,*999)

              elementDomains%offsets(1)=1
              DO elementIdx=1,numberOfLocalElements
                elementDomains%offsets(elementIdx+1)=elementDomains%offsets(elementIdx)+1
                CALL List_ItemAdd(elementDomainsList,partition(elementIdx),err,error,*999)
                CALL List_Search(boundaryElements,elementIdx,boundaryElementIdx,err,error,*999)
                IF(boundaryElementIdx>0) THEN
                  CALL List_RemoveDuplicates(newAdjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                  CALL List_NumberOfItemsGet(newAdjacentDomainsLists(boundaryElementIdx)%ptr,numberOfAdjacentDomains,err,error,*999)
                  elementDomains%offsets(elementIdx+1)=elementDomains%offsets(elementIdx+1)+numberOfAdjacentDomains
                  CALL List_AppendList(elementDomainsList,newAdjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                  CALL List_Destroy(newAdjacentDomainsLists(boundaryElementIdx)%ptr,err,error,*999)
                END IF
              END DO !elementIdx
              IF(ALLOCATED(partition)) DEALLOCATE(partition)
              IF(ALLOCATED(boundaryElements)) DEALLOCATE(boundaryElements)
              IF(ALLOCATED(newAdjacentDomainsLists)) DEALLOCATE(newAdjacentDomainsLists)
              CALL List_DetachAndDestroy(elementDomainsList,numberOfElementDomains,elementDomains%domains,err,error,*999)
            END IF
          CASE(DECOMPOSITION_USER_DEFINED_TYPE)
            !The user should give for each initial local element a domain number and based on that ghost elements should be
            !calculated and eventually all the domain numbers.
            CALL FlagError("Not implemented.",err,error,*999)            
          CASE DEFAULT
            CALL FlagError("Invalid domain decomposition type.",err,error,*999)            
          END SELECT
        ELSE
          CALL FlagError("Decomposition mesh topology is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Decomposition is not associated.",err,error,*999)
    END IF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition for mesh number ",decomposition%mesh%USER_NUMBER, &
        & err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains = ", decomposition%NUMBER_OF_DOMAINS, &
        & err,error,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Partition:",err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposition type = ", decomposition%DECOMPOSITION_TYPE, &
        & err,error,*999)
      IF(decomposition%DECOMPOSITION_TYPE==DECOMPOSITION_CALCULATED_TYPE) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of edges cut = ",numberOfCutEdges,err,error,*999)
      END IF
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local elements = ",numberOfLocalElements, &
        & err,error,*999)
      DO elementIdx=1,numberOfLocalElements
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global element number = ",domainOffsets(myDomain)+elementIdx-1, &
          & err,error,*999)
        elementDomains=>decomposition%elementDomains
        IF(.NOT.ASSOCIATED(elementDomains)) CALL FlagError("Decomposition element domains is not associated.",err,error,*999)
        DO elementDomainIdx=elementDomains%offsets(elementIdx),elementDomains%offsets(elementIdx+1)-1        
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Domain = ",elementDomains%domains(elementDomainIdx), &
            & err,error,*999)
        END DO !elementDomainIdx
      END DO !elementIdx
    END IF
    IF(ALLOCATED(domainOffsets)) DEALLOCATE(domainOffsets)
    
    EXITS("Decomposition_ElementDomainsCalculate")
    RETURN
999 IF(ALLOCATED(domainOffsets)) DEALLOCATE(domainOffsets)
    IF(ALLOCATED(elementNodeOffsets)) DEALLOCATE(elementNodeOffsets)
    IF(ALLOCATED(elementNodeNumbers)) DEALLOCATE(elementNodeNumbers)
    IF(ASSOCIATED(adjacentElementOffsets)) DEALLOCATE(adjacentElementOffsets)
    IF(ASSOCIATED(adjacentElements)) DEALLOCATE(adjacentElements)
    IF(ALLOCATED(boundaryElements)) DEALLOCATE(boundaryElements)
    IF(ALLOCATED(adjacentDomainNumbers)) DEALLOCATE(adjacentDomainNumbers)
    IF(ALLOCATED(partition)) DEALLOCATE(partition)
    IF(ALLOCATED(tpwgts)) DEALLOCATE(tpwgts)
    IF(ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
    IF(ALLOCATED(receiveCounts)) DEALLOCATE(receiveCounts)
    IF(ALLOCATED(requests)) DEALLOCATE(requests)
    IF(ALLOCATED(statuses)) DEALLOCATE(statuses)
    IF(ALLOCATED(receiveBuffers)) THEN
      DO domainIdx=1,SIZE(receiveBuffers)
        IF(ALLOCATED(receiveBuffers(domainIdx)%array)) DEALLOCATE(receiveBuffers(domainIdx)%array)
      END DO
      DEALLOCATE(receiveBuffers)
    END IF
    IF(ALLOCATED(sendBuffers)) THEN
      DO domainIdx=1,SIZE(sendBuffers)
        IF(ALLOCATED(sendBuffers(domainIdx)%array)) DEALLOCATE(sendBuffers(domainIdx)%array)
      END DO
      DEALLOCATE(sendBuffers)
    END IF
    IF(ASSOCIATED(boundaryElementsList)) CALL List_Destroy(boundaryElementsList,dummyErr,dummyError,*998)
998 IF(ASSOCIATED(elementDomainsList)) CALL List_Destroy(elementDomainsList,dummyErr,dummyError,*997)
997 IF(ASSOCIATED(ghostOwnerTree)) CALL Tree_Destroy(ghostOwnerTree,dummyErr,dummyError,*996)
996 IF(ALLOCATED(adjacentDomainsLists)) THEN
      DO boundaryElementIdx=1,SIZE(adjacentDomainsLists)
        CALL List_Destroy(adjacentDomainsLists(boundaryElementIdx)%ptr,dummyErr,dummyError,*995)
      END DO
      DEALLOCATE(adjacentDomainsLists)
    END IF
995 IF(ALLOCATED(newAdjacentDomainsLists)) THEN
      DO boundaryElementIdx=1,SIZE(newAdjacentDomainsLists)
        CALL List_Destroy(newAdjacentDomainsLists(boundaryElementIdx)%ptr,dummyErr,dummyError,*994)
      END DO
      DEALLOCATE(newAdjacentDomainsLists)
    END IF
994 IF(ALLOCATED(sendBufferLists)) THEN
      DO domainIdx=1,SIZE(sendBufferLists)
        CALL List_Destroy(sendBufferLists(domainIdx)%ptr,dummyErr,dummyError,*993)
      END DO
      DEALLOCATE(sendBufferLists)
    END IF
993 ERRORSEXITS("Decomposition_ElementDomainsCalculate",err,error)
    RETURN 1
  END SUBROUTINE Decomposition_ElementDomainsCalculate
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::CMISSDecompositionMeshComponentGet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the mesh component for
    INTEGER(INTG), INTENT(OUT) :: MESH_COMPONENT_NUMBER !<On return, the mesh component number to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET
  
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::CMISSDecompositionMeshComponentSet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS) THEN
            DECOMPOSITION%MESH_COMPONENT_NUMBER=MESH_COMPONENT_NUMBER
          ELSE
            LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
              & "is invalid. The component number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the number of domains for a decomposition. \see OPENCMISS::CMISSDecompositionNumberOfDomainsGet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the number of domains for.
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_DOMAINS !<On return, the number of domains to get.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        NUMBER_OF_DOMAINS=DECOMPOSITION%NUMBER_OF_DOMAINS
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition. \see OPENCMISS::CMISSDecompositionNumberOfDomainsSet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the number of domains for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !<The number of domains to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          IF(NUMBER_OF_DOMAINS==1) THEN
            DECOMPOSITION%NUMBER_OF_DOMAINS=1
          ELSE
            CALL FlagError("Can only have one domain for all decomposition type.",ERR,ERROR,*999)
          ENDIF
        CASE(DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE)
          IF(NUMBER_OF_DOMAINS>=1) THEN
            !wolfye???<=?
            IF(NUMBER_OF_DOMAINS<=DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS) THEN
              DECOMPOSITION%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS             
            ELSE
              LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINS,"*",ERR,ERROR))// &
                & ") must be <= the number of global elements ("// &
                & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//") in the mesh."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
             CALL FlagError("Number of domains must be >= 1.",ERR,ERROR,*999)
           ENDIF
         CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%DECOMPOSITION_TYPE,"*",ERR,ERROR))// &
            & " is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  !
  !================================================================================================================================
  !

  !>Calculates the topology for a decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
      !Calculate the elements topology
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      !Calculate the line topology
      IF(DECOMPOSITION%CALCULATE_LINES)THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      !Calculate the face topology
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      meshComponentNumber=DECOMPOSITION%MESH_COMPONENT_NUMBER
      IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsCalculate(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF   
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DecompositionTopology_DataPointsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: localElement,globalElement,dataPointIdx,localData,meshComponentNumber
    INTEGER(INTG) :: INSERT_STATUS,NUMBER_OF_COMPUTATIONAL_NODES
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(MeshDataPointsType), POINTER :: meshData

    ENTERS("DecompositionTopology_DataPointsCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      decompositionData=>TOPOLOGY%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        decomposition=>decompositionData%DECOMPOSITION
        IF(ASSOCIATED(decomposition)) THEN
         decompositionElements=>TOPOLOGY%ELEMENTS
         IF(ASSOCIATED(decompositionElements)) THEN
           elementsMapping=>decomposition%DOMAIN(decomposition%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
           IF(ASSOCIATED(elementsMapping)) THEN
              meshComponentNumber=decomposition%MESH_COMPONENT_NUMBER
              meshData=>decomposition%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints
              IF(ASSOCIATED(meshData)) THEN
                NUMBER_OF_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
                IF(ERR/=0) GOTO 999
                ALLOCATE(decompositionData%elementDataPoint(decompositionElements%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate decomposition element data points.",ERR,ERROR,*999)
                CALL Tree_CreateStart(decompositionData%dataPointsTree,ERR,ERROR,*999)
                CALL Tree_InsertTypeSet(decompositionData%dataPointsTree,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                CALL Tree_CreateFinish(decompositionData%dataPointsTree,ERR,ERROR,*999)
                decompositionData%numberOfGlobalDataPoints=meshData%totalNumberOfProjectedData 
                localData=0
                DO localElement=1,decompositionElements%TOTAL_NUMBER_OF_ELEMENTS
                  globalElement=decompositionElements%ELEMENTS(localElement)%GLOBAL_NUMBER
                  decompositionData%elementDataPoint(localElement)%numberOfProjectedData= &
                    & meshData%elementDataPoint(globalElement)%numberOfProjectedData
                  decompositionData%elementDataPoint(localElement)%globalElementNumber=globalElement
                  IF(localElement<elementsMapping%NUMBER_OF_LOCAL) THEN
                    decompositionData%numberOfDataPoints=decompositionData%numberOfDataPoints+ &
                      & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  ENDIF               
                  decompositionData%totalNumberOfDataPoints=decompositionData%totalNumberOfDataPoints+ &
                    & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  ALLOCATE(decompositionData%elementDataPoint(localElement)%dataIndices(decompositionData% &
                    & elementDataPoint(localElement)%numberOfProjectedData),STAT=ERR)
                  DO dataPointIdx=1,decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%userNumber
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%globalNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%globalNumber
                    localData=localData+1
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%localNumber=localData
                    CALL Tree_ItemInsert(decompositionData%dataPointsTree,decompositionData% &
                      & elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber,localData, &
                      & INSERT_STATUS,ERR,ERROR,*999)
                  ENDDO !dataPointIdx
                ENDDO !localElement   
              ELSE
                CALL FlagError("Mesh data points topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Element mapping  is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition elements topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition data points topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsCalculate")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsCalculate",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataPointsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology for a data projection (for data projections on fields).
  SUBROUTINE DecompositionTopology_DataProjectionCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DecompositionTopology_DataProjectionCalculate",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      CALL DecompositionTopology_DataPointsInitialise(decompositionTopology,err,error,*999)
      CALL DecompositionTopology_DataPointsCalculate(decompositionTopology,err,error,*999)
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN
999 ERRORS("DecompositionTopology_DataProjectionCalculate",err,error)
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataProjectionCalculate

  !
  !================================================================================================================================
  !

  !>Gets the local data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet(decompositionTopology,elementNumber,dataPointIndex, &
       & dataPointLocalNumber,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointLocalNumber !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        numberOfDataPoints = decompositionData%elementDataPoint(elementNumber)%numberOfProjectedData
        IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
          dataPointLocalNumber = decompositionData%elementDataPoint(elementNumber)%dataIndices(dataPointIndex)%localNumber
        ELSE
          localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
           & " out of range for element "//TRIM(NUMBER_TO_VSTRING(elementNumber,"*",ERR,ERROR))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the user (global) data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet(decompositionTopology,userElementNumber,dataPointIndex, &
       & dataPointUserNumber,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointUserNumber !<The data point user number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints,decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointUserNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints = decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
            IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
              dataPointUserNumber = decompositionData%elementDataPoint(decompositionLocalElementNumber)% &
                & dataIndices(dataPointIndex)%userNumber
            ELSE
              localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
               & " out of range for element "//TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",ERR,ERROR))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointUserNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of data points projected on an element
  SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet(decompositionTopology,userElementNumber, &
       & numberOfDataPoints,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(OUT) :: numberOfDataPoints !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_NumberOfElementDataPointsGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints=decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN
999 ERRORS("DecompositionTopology_NumberOfElementDataPointsGet",err,error)
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DecompositionTopology_DataPointCheckExists(decompositionTopology,userDataPointNumber,userDataPointExists, &
        & decompositionLocalDataPointNumber,ghostDataPoint,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to check the data point exists on
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: userDataPointExists !<On exit, is .TRUE. if the data point user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: decompositionLocalDataPointNumber !<On exit, if the data point exists the local number corresponding to the user data point number. If the data point does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostDataPoint !<On exit, is .TRUE. if the local data point (if it exists) is a ghost data point, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("DecompositionTopology_DataPointCheckExists",ERR,error,*999)

    userDataPointExists=.FALSE.
    decompositionLocalDataPointNumber=0
    ghostDataPoint=.FALSE.
    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        NULLIFY(treeNode)
        CALL Tree_Search(decompositionData%dataPointsTree,userDataPointNumber,treeNode,err,error,*999)
        IF(ASSOCIATED(treeNode)) THEN
          CALL Tree_NodeValueGet(decompositionData%dataPointsTree,treeNode,decompositionLocalDataPointNumber,err,error,*999)
          userDataPointExists=.TRUE.
          ghostDataPoint=decompositionLocalDataPointNumber>decompositionData%numberOfDataPoints
        ENDIF
      ELSE
        CALL FlagError("Decomposition data point topology is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN
999 ERRORS("DecompositionTopology_DataPointCheckExists",err,error)
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointCheckExists
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(DECOMPOSITION_TOPOLOGY,USER_ELEMENT_NUMBER,ELEMENT_EXISTS, &
    & DECOMPOSITION_LOCAL_ELEMENT_NUMBER,GHOST_ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: DECOMPOSITION_TOPOLOGY !<A pointer to the decomposition topology to check the element exists on
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: ELEMENT_EXISTS !<On exit, is .TRUE. if the element user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DECOMPOSITION_LOCAL_ELEMENT_NUMBER !<On exit, if the element exists the local number corresponding to the user element number. If the element does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_ELEMENT !<On exit, is .TRUE. if the local element (if it exists) is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR,*999)

    ELEMENT_EXISTS=.FALSE.
    DECOMPOSITION_LOCAL_ELEMENT_NUMBER=0
    GHOST_ELEMENT=.FALSE.
    IF(ASSOCIATED(DECOMPOSITION_TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>DECOMPOSITION_TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        NULLIFY(TREE_NODE)
        CALL Tree_Search(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL Tree_NodeValueGet(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,TREE_NODE,DECOMPOSITION_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
          ELEMENT_EXISTS=.TRUE.
          GHOST_ELEMENT=DECOMPOSITION_LOCAL_ELEMENT_NUMBER>DECOMPOSITION_ELEMENTS%NUMBER_OF_ELEMENTS
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS

  !
  !================================================================================================================================
  !

  !>Get the basis for an element in the domain identified by its local number
  SUBROUTINE DomainTopology_ElementBasisGet(domainTopology,userElementNumber, &
      & basis,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the element basis for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element user number to get the basis for
    TYPE(BASIS_TYPE), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    LOGICAL :: userElementExists,ghostElement
    INTEGER(INTG) :: localElementNumber

    ENTERS("DomainTopology_ElementBasisGet",err,error,*999)

    NULLIFY(basis)

    IF(ASSOCIATED(domainTopology)) THEN
      domainElements=>domainTopology%elements
      IF(ASSOCIATED(domainElements)) THEN
        decompositionTopology=>domainTopology%domain%decomposition%topology
        IF(ASSOCIATED(decompositionTopology)) THEN
          CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
            & userElementExists,localElementNumber,ghostElement,err,error,*999)
          IF(.NOT.userElementExists) THEN
            CALL FlagError("The specified user element number of "// &
              & TRIM(NumberToVstring(userElementNumber,"*",err,error))// &
              & " does not exist in the domain decomposition.",err,error,*999)
          END IF
          basis=>domainElements%elements(localElementNumber)%basis
        ELSE
          CALL FlagError("Decomposition topology is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Domain topology elements is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Domain topology is not associated.",err,error,*999)
    END IF

    EXITS("DomainTopology_ElementBasisGet")
    RETURN
    999 ERRORSEXITS("DomainTopology_ElementBasisGet",err,error)
    RETURN 1

    END SUBROUTINE DomainTopology_ElementBasisGet

  !
  !================================================================================================================================
  !

  !>Finalises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic !\todo add comment

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) THEN
      DO nic=LBOUND(ELEMENT%ADJACENT_ELEMENTS,1),UBOUND(ELEMENT%ADJACENT_ELEMENTS,1)
        CALL DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(ELEMENT%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    ENDIF
    IF(ALLOCATED(ELEMENT%ELEMENT_LINES)) DEALLOCATE(ELEMENT%ELEMENT_LINES)
    IF(ALLOCATED(ELEMENT%ELEMENT_FACES)) DEALLOCATE(ELEMENT%ELEMENT_FACES)
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%LOCAL_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    ELEMENT%BOUNDARY_ELEMENT=.FALSE.
  
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers adjacent to an element in a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the adjacent elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),NODE_POSITION_INDEX(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: candidate_idx,face_node_idx,node_idx,surrounding_el_idx,candidate_el,idx
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4),numberSurroundingElements
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:), surroundingElements(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST, surroundingElementsList
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    ENTERS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                  !Loop over the elements in the decomposition
                  DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                    !Create a list for every xi direction (plus and minus)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
                      CALL List_CreateStart(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                      CALL List_DataTypeSet(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL List_InitialSizeSet(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
                      CALL List_CreateFinish(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                    ENDDO !nic
                    NUMBER_OF_NODES_XIC=1
                    NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)= &
                      & BASIS%NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)
                    !Place the current element in the surrounding list
                    CALL List_ItemAdd(ADJACENT_ELEMENTS_LIST(0)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER, &
                      & ERR,ERROR,*999)
                    SELECT CASE(BASIS%TYPE)
                    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
!!TODO: Calculate this and set it as part of the basis type
                      !Determine the collapsed "faces" if any
                      NODE_POSITION_INDEX=1
                      !Loop over the face normals of the element
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the face xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Reset the node_position_index in this xi direction
                        NODE_POSITION_INDEX(ni)=1
                        !Loop over the two faces with this normal
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni
                          FACE_COLLAPSED(xi_direction)=.FALSE.
                          DO j=1,2
                            xi_dir_check=FACE_XI(j)
                            IF(xi_dir_check<=BASIS%NUMBER_OF_XI) THEN
                              xi_dir_search=FACE_XI(3-j)
                              NODE_POSITION_INDEX(xi_dir_search)=1
                              XI_COLLAPSED=.TRUE.
                              DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                                !Get the first local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=1
                                nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                !Get the second local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=2
                                nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                IF(nn1/=0.AND.nn2/=0) THEN
                                  IF(DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn1)/= &
                                    & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn2)) XI_COLLAPSED=.FALSE.
                                ENDIF
                                NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                              ENDDO !xi_dir_search
                              IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                            ENDIF
                          ENDDO !j
                          NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                        ENDDO !direction_index
                      ENDDO !ni
                      !Loop over the xi directions and calculate the surrounding elements
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Loop over the two faces
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni                  
                          !Find nodes in the element on the appropriate face/line/point
                          NULLIFY(NODE_MATCH_LIST)
                          CALL List_CreateStart(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL List_DataTypeSet(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL List_InitialSizeSet(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                          CALL List_CreateFinish(NODE_MATCH_LIST,ERR,ERROR,*999)
                          IF(direction_index==-1) THEN
                            NODE_POSITION_INDEX(ni)=1
                          ELSE
                            NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                          ENDIF
                          !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is
                          !also collapsed. This may indicate that we have a funny element in non-rc coordinates that goes around the
                          !central axis back to itself
                          IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                            !Do nothing - the match lists are already empty
                          ELSE
                            !Find the nodes to match and add them to the node match list
                            SELECT CASE(BASIS%NUMBER_OF_XI)
                            CASE(1)
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),1,1,1)
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL List_ItemAdd(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            CASE(2)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),1,1)
                                IF(nn/=0) THEN
                                  np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                  CALL List_ItemAdd(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                ENDIF
                              ENDDO !nn1
                            CASE(3)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2)),NUMBER_OF_NODES_XIC(FACE_XI(2))-1
                                  NODE_POSITION_INDEX(FACE_XI(2))=nn2
                                  nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                    & NODE_POSITION_INDEX(3),1)
                                  IF(nn/=0) THEN
                                    np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                    CALL List_ItemAdd(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !nn2
                              ENDDO !nn1
                            CASE DEFAULT
                              LOCAL_ERROR="The number of xi directions in the basis of "// &
                                & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//" is invalid."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ENDIF
                          CALL List_RemoveDuplicates(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL List_DetachAndDestroy(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=0
                          IF(NUMBER_NODE_MATCHES>0) THEN
                            !NODE_MATCHES now contain the list of corner nodes in the current face with normal_xi of ni.
                            !Look at the surrounding elements of each of these nodes, if there is a repeated element that
                            !is not the current element ne, it's an adjacent element.
                            candidate_idx=0
                            NULLIFY(surroundingElementsList)
                            CALL List_CreateStart(surroundingElementsList,ERR,ERROR,*999)
                            CALL List_DataTypeSet(surroundingElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                            CALL List_InitialSizeSet(surroundingElementsList,2,ERR,ERROR,*999)
                            CALL List_CreateFinish(surroundingElementsList,ERR,ERROR,*999)
                            DO face_node_idx=1,NUMBER_NODE_MATCHES
                              !Dump all the surrounding elements into an array, see if any are repeated
                              node_idx=NODE_MATCHES(face_node_idx)
                              DO surrounding_el_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                                candidate_el=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(surrounding_el_idx)
                                IF(candidate_el/=ne) THEN
                                  candidate_idx=candidate_idx+1
                                  CALL List_ItemAdd(surroundingElementsList,candidate_el,ERR,ERROR,*999)
                                ENDIF
                              ENDDO
                            ENDDO !face_node_idx
                            CALL List_DetachAndDestroy(surroundingElementsList,numberSurroundingElements,surroundingElements, &
                              & ERR,ERROR,*999)
                            DO idx=1,candidate_idx
                              ne1=surroundingElements(idx)
                              IF(COUNT(surroundingElements(1:numberSurroundingElements)==ne1)>=BASIS%NUMBER_OF_XI) THEN
                                !Found it, just exit
                                CALL List_ItemAdd(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                                NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                                EXIT
                              ENDIF
                            ENDDO
                          ENDIF
                          IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                          IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
                        ENDDO !direction_index
                      ENDDO !ni
                    CASE(BASIS_SIMPLEX_TYPE)
                      !Loop over the xi coordinates and calculate the surrounding elements
                      DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                        !Find the other coordinates of the face/line/point
                        FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                        FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                        FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                        !Find nodes in the element on the appropriate face/line/point
                        NULLIFY(NODE_MATCH_LIST)
                        CALL List_CreateStart(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL List_DataTypeSet(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                        CALL List_InitialSizeSet(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                        CALL List_CreateFinish(NODE_MATCH_LIST,ERR,ERROR,*999)
                        NODE_POSITION_INDEX(nic)=1 !Furtherest away from node with the nic'th coordinate
                        !Find the nodes to match and add them to the node match list
                        DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                          NODE_POSITION_INDEX(FACE_XIC(1))=nn1
                          DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                            NODE_POSITION_INDEX(FACE_XIC(2))=nn2
                            DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                              NODE_POSITION_INDEX(FACE_XIC(3))=nn3
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                & NODE_POSITION_INDEX(3),NODE_POSITION_INDEX(4))
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL List_ItemAdd(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !nn3
                          ENDDO !nn2
                        ENDDO !nn1
                        CALL List_RemoveDuplicates(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL List_DetachAndDestroy(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                        IF(NUMBER_NODE_MATCHES>0) THEN
                          !Find list of elements surrounding those nodes
                          DO node_idx=1,NUMBER_NODE_MATCHES
                            np1=NODE_MATCHES(node_idx)
                            DO nep1=1,DOMAIN_NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                              ne1=DOMAIN_NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                              IF(ne1/=ne) THEN !Don't want the current element
                                ! grab the nodes list for current and this surrouding elements
                                ! current face : NODE_MATCHES
                                ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES 
                                ! if all of current face belongs to the candidate element, we will have found the neighbour
                                CALL List_SubsetOf(NODE_MATCHES(1:NUMBER_NODE_MATCHES),DOMAIN_ELEMENTS%ELEMENTS(ne1)% &
                                  & ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                                IF(SUBSET) THEN
                                  CALL List_ItemAdd(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                            ENDDO !nep1
                          ENDDO !node_idx
                        ENDIF
                        IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                      ENDDO !nic
                    CASE(BASIS_SERENDIPITY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_AUXILLIARY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_B_SPLINE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)                     
                    END SELECT
                    !Set the surrounding elements for this element
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
                      BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements.",ERR,ERROR,*999)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      CALL DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic), &
                        & ERR,ERROR,*999)
                      CALL List_DetachAndDestroy(ADJACENT_ELEMENTS_LIST(nic)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS,ERR,ERROR,*999)
                      ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS( &
                        DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate element adjacent elements.",ERR,ERROR,*999)
                      DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(1: &
                        & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)= &
                        ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)
                      IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
                    ENDDO !nic
                  ENDDO !ne           
                ELSE
                  CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not allocated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ",DECOMPOSITION_ELEMENTS% &
        & TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number : ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ERR,ERROR,*999)
          IF(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
              & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,8,8,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)% &
              & ADJACENT_ELEMENTS,'("        Adjacent elements :",8(X,I6))','(30x,8(X,I6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF

    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL List_Destroy(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(surroundingElementsList)) CALL List_Destroy(surroundingElementsList,DUMMY_ERR,DUMMY_ERROR,*997)
997 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL List_Destroy(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*996)
    ENDDO !ni
996 ERRORS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR)
    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the boundary nodes and elements for a decomposition topology. 
  SUBROUTINE DecompositionTopology_BoundaryCalculate(topology,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: topology !<A pointer to the decomposition topology to calculate the boundary for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LOGICAL_ALLOC_TYPE), ALLOCATABLE :: sendBuffers(:),receiveBuffers(:)
    INTEGER(INTG) :: elementIdx,elementNodeIdx,matchIndex,nodeIdx,xiCoordIdx,adjacentDomainIdx, &
      & ghostSendIdx,ghostReceiveIdx,adjacentDomainNumber,meshComponentIdx,MPI_IERROR
    INTEGER(INTG), ALLOCATABLE :: requests(:),statuses(:,:)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition 
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DecompositionTopology_BoundaryCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      decomposition=>topology%decomposition
      IF(ASSOCIATED(decomposition)) THEN
        decompositionElements=>topology%elements
        IF(ASSOCIATED(decompositionElements)) THEN
          domainElements=>decomposition%domain(decomposition%MESH_COMPONENT_NUMBER)%ptr%topology%elements
          IF(.NOT.ASSOCIATED(domainElements)) &
            & CALL FlagError("Domain elements is not associated.",err,error,*999)
          !Only loop over local elements, not ghost elements!
          DO elementIdx=1,decompositionElements%NUMBER_OF_ELEMENTS
            basis=>domainElements%elements(elementIdx)%basis
            SELECT CASE(basis%type)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              DO xiCoordIdx=-basis%NUMBER_OF_XI_COORDINATES,basis%NUMBER_OF_XI_COORDINATES
                IF(xiCoordIdx/=0) THEN
                  IF(decompositionElements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)% &
                    & NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                    decompositionElements%elements(elementIdx)%BOUNDARY_ELEMENT=.TRUE.
                    EXIT
                  ENDIF
                ENDIF
              ENDDO !xiCoordIdx            
            CASE(BASIS_SIMPLEX_TYPE)
              DO xiCoordIdx=1,basis%NUMBER_OF_XI_COORDINATES
                IF(decompositionElements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)% &
                  & NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                  decompositionElements%elements(elementIdx)%BOUNDARY_ELEMENT=.TRUE.
                  EXIT
                ENDIF
              ENDDO !xiCoordIdx
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !elementIdx
           
          elementsMapping=>decomposition%domain(decomposition%MESH_COMPONENT_NUMBER)%ptr%mappings%elements
          IF(.NOT.ASSOCIATED(elementsMapping)) CALL FlagError("Elements mapping is not associated.",err,error,*999)

          !Inform the ghost elements if they are boundary.                
          ALLOCATE(sendBuffers(elementsMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate send buffers.",err,error,*999)

          !Pack the send buffers
          DO adjacentDomainIdx=1,elementsMapping%NUMBER_OF_ADJACENT_DOMAINS
            ALLOCATE(sendBuffers(adjacentDomainIdx)%array( &
              & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS),stat=err)
            IF(err/=0) CALL FlagError("Could not allocate send buffer.",err,error,*999)
            DO ghostSendIdx=1,elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS
              elementIdx=elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES(ghostSendIdx)
              sendBuffers(adjacentDomainIdx)%array(ghostSendIdx)=decompositionElements%elements(elementIdx)%BOUNDARY_ELEMENT
            END DO !ghostSendIdx
          END DO !adjacentDomainIdx

          ALLOCATE(statuses(MPI_STATUS_SIZE,2*elementsMapping%NUMBER_OF_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate mpi statuses.",err,error,*999)
          ALLOCATE(requests(2*elementsMapping%NUMBER_OF_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate mpi receive requests.",err,error,*999)
          requests=MPI_REQUEST_NULL

          ALLOCATE(receiveBuffers(elementsMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate receive buffers.",err,error,*999)
          
          !Get all the receive buffers from the adjacent domains.
          DO adjacentDomainIdx=1,elementsMapping%NUMBER_OF_ADJACENT_DOMAINS
            ALLOCATE(receiveBuffers(adjacentDomainIdx)%array( &
              & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS),stat=err)
            IF(err/=0) CALL FlagError("Could not allocate receive buffer.",err,error,*999)
            adjacentDomainNumber=elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
            CALL MPI_IRECV(receiveBuffers(adjacentDomainIdx)%array, &
              & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS,MPI_LOGICAL, &
              & adjacentDomainNumber-1,MPI_BOUNDARY_ELEMENT_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
              & requests(adjacentDomainNumber),MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,err,error,*999)
          END DO
          DO adjacentDomainIdx=1,elementsMapping%NUMBER_OF_ADJACENT_DOMAINS
            adjacentDomainNumber=elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
            CALL MPI_ISEND(sendBuffers(adjacentDomainIdx)%array, &
              & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS,MPI_LOGICAL, &
              & adjacentDomainNumber-1,MPI_BOUNDARY_ELEMENT_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
              & requests(elementsMapping%NUMBER_OF_DOMAINS+adjacentDomainNumber),MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,err,error,*999)
          END DO

          CALL MPI_WAITALL(2*elementsMapping%NUMBER_OF_DOMAINS,requests,statuses,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,err,error,*999)

          IF(ALLOCATED(statuses)) DEALLOCATE(statuses)
          IF(ALLOCATED(requests)) DEALLOCATE(requests)

          !Unpack the receive buffers
          DO adjacentDomainIdx=1,elementsMapping%NUMBER_OF_ADJACENT_DOMAINS
            DO ghostReceiveIdx=1,elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS
              elementIdx=elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES(ghostReceiveIdx)
              decompositionElements%elements(elementIdx)%BOUNDARY_ELEMENT=receiveBuffers(adjacentDomainIdx)%array(ghostReceiveIdx)
            END DO !ghostReceiveIdx
          END DO !adjacentDomainIdx

          IF(ALLOCATED(sendBuffers)) THEN
            DO adjacentDomainIdx=1,elementsMapping%NUMBER_OF_ADJACENT_DOMAINS
              IF(ALLOCATED(sendBuffers(adjacentDomainIdx)%array)) DEALLOCATE(sendBuffers(adjacentDomainIdx)%array)
            END DO
            DEALLOCATE(sendBuffers)
          END IF
          IF(ALLOCATED(receiveBuffers)) THEN
            DO adjacentDomainIdx=1,elementsMapping%NUMBER_OF_ADJACENT_DOMAINS
              IF(ALLOCATED(receiveBuffers(adjacentDomainIdx)%array)) DEALLOCATE(receiveBuffers(adjacentDomainIdx)%array)
            END DO
            DEALLOCATE(receiveBuffers)
          END IF
          
          !Now calculate the boundary nodes for all mesh components
          DO meshComponentIdx=1,SIZE(decomposition%domain)
            domainNodes=>decomposition%domain(meshComponentIdx)%ptr%topology%nodes
            domainElements=>decomposition%domain(meshComponentIdx)%ptr%topology%elements
            DO elementIdx=1,decompositionElements%TOTAL_NUMBER_OF_ELEMENTS
              IF(decompositionElements%elements(elementIdx)%BOUNDARY_ELEMENT) THEN
                basis=>domainElements%elements(elementIdx)%basis
                SELECT CASE(basis%type)
                CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                  DO xiCoordIdx=-basis%NUMBER_OF_XI_COORDINATES,basis%NUMBER_OF_XI_COORDINATES
                    IF(decompositionElements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)% &
                      & NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                      IF(xiCoordIdx/=0) THEN
                        IF(xiCoordIdx<0) THEN
                          matchIndex=1
                        ELSE
                          matchIndex=BASIS%NUMBER_OF_NODES_XIC(xiCoordIdx)
                        ENDIF                    
                        DO elementNodeIdx=1,BASIS%NUMBER_OF_NODES
                          IF(basis%NODE_POSITION_INDEX(elementNodeIdx,ABS(xiCoordIdx))==matchIndex) THEN
                            nodeIdx=domainElements%elements(elementIdx)%ELEMENT_NODES(elementNodeIdx)
                            domainNodes%nodes(nodeIdx)%BOUNDARY_NODE=.TRUE.
                          ENDIF
                        ENDDO !nn
                      ENDIF
                    END IF
                  ENDDO !xiCoordIdx            
                CASE(BASIS_SIMPLEX_TYPE)
                  DO xiCoordIdx=1,basis%NUMBER_OF_XI_COORDINATES
                    IF(decompositionElements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)% &
                      & NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                      DO elementNodeIdx=1,basis%NUMBER_OF_NODES
                        IF(basis%NODE_POSITION_INDEX(elementNodeIdx,xiCoordIdx)==1) THEN
                          nodeIdx=domainElements%elements(elementIdx)%ELEMENT_NODES(elementNodeIdx)
                          domainNodes%nodes(nodeIdx)%BOUNDARY_NODE=.TRUE.
                        ENDIF
                      ENDDO !elementNodeIdx                  
                    ENDIF
                  ENDDO !xiCoordIdx
                CASE(BASIS_SERENDIPITY_TYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(BASIS_AUXILLIARY_TYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(BASIS_B_SPLINE_TP_TYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              END IF
            END DO !elementIdx
          END DO !meshComponentIdx
        ELSE
          CALL FlagError("Decomposition topology elements is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary elements:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ", &
        & decompositionElements%TOTAL_NUMBER_OF_ELEMENTS,err,error,*999)
      DO elementIdx=1,decompositionElements%TOTAL_NUMBER_OF_ELEMENTS
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Element : ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary element = ",decompositionElements%elements(elementIdx)% &
          & BOUNDARY_ELEMENT,err,error,*999)        
      ENDDO !elementIdx
    ENDIF
 
    EXITS("DecompositionTopology_BoundaryCalculate")
    RETURN
999 IF(ALLOCATED(statuses)) DEALLOCATE(statuses)
    IF(ALLOCATED(requests)) DEALLOCATE(requests)
    IF(ALLOCATED(sendBuffers)) THEN
      DO adjacentDomainIdx=1,SIZE(sendBuffers)
        IF(ALLOCATED(sendBuffers(adjacentDomainIdx)%array)) DEALLOCATE(sendBuffers(adjacentDomainIdx)%array)
      END DO
      DEALLOCATE(sendBuffers)
    END IF
    IF(ALLOCATED(receiveBuffers)) THEN
      DO adjacentDomainIdx=1,SIZE(receiveBuffers)
        IF(ALLOCATED(receiveBuffers(adjacentDomainIdx)%array)) DEALLOCATE(receiveBuffers(adjacentDomainIdx)%array)
      END DO
      DEALLOCATE(receiveBuffers)
    END IF
    ERRORSEXITS("DecompositionTopology_BoundaryCalculate",err,error)

    RETURN 1
  END SUBROUTINE DecompositionTopology_BoundaryCalculate

  !
  !===============================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: global_element,local_element
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        IF(ASSOCIATED(DECOMPOSITION)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
              IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                  DOMAIN_ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS_MAPPING)) THEN
                    MESH=>DECOMPOSITION%MESH
                    IF(ASSOCIATED(MESH)) THEN
                      MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
                      IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
                        MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(MESH_ELEMENTS)) THEN
                          !Allocate the element topology arrays
                          ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate decomposition elements elements.",ERR,ERROR,*999)
                          DECOMPOSITION_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS
                          DO local_element=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                            CALL DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(local_element), &
                              & ERR,ERROR,*999)
                            global_element=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_element)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%USER_NUMBER=MESH_ELEMENTS%ELEMENTS(global_element)% &
                              & USER_NUMBER
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%LOCAL_NUMBER=local_element
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%GLOBAL_NUMBER=global_element
                          ENDDO !local_element
                          !Calculate the elements surrounding the elements in the decomposition topology
                          CALL DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*999)
                          CALL DecompositionTopology_BoundaryCalculate(TOPOLOGY,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given decomposition topology. \todo Pass in the decomposition elements pointer.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)) CALL Tree_Destroy(TOPOLOGY%ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Decomposition already has topology elements associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements.",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given decomposition. \todo Pass in a pointer to the decomposition topology
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      IF(DECOMPOSITION%CALCULATE_LINES) THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      DEALLOCATE(DECOMPOSITION%TOPOLOGY)
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
        CALL FlagError("Decomposition already has topology associated.",ERR,ERROR,*999)
      ELSE
        !Allocate decomposition topology
        ALLOCATE(DECOMPOSITION%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Decomposition topology could not be allocated.",ERR,ERROR,*999)
        DECOMPOSITION%TOPOLOGY%DECOMPOSITION=>DECOMPOSITION
        NULLIFY(DECOMPOSITION%TOPOLOGY%ELEMENTS)
        NULLIFY(DECOMPOSITION%TOPOLOGY%LINES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%FACES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%dataPoints)
        !Initialise the topology components
        CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        IF(DECOMPOSITION%CALCULATE_LINES) THEN !Default is currently true
          CALL DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF
        IF(DECOMPOSITION%CALCULATE_FACES) THEN !Default is currently false
          CALL DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF      
        meshComponentNumber=DECOMPOSITION%MESH_COMPONENT_NUMBER
        IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsInitialise(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF   
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(LINE%SURROUNDING_ELEMENTS)) DEALLOCATE(LINE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(LINE%ELEMENT_LINES)) DEALLOCATE(LINE%ELEMENT_LINES)
    LINE%ADJACENT_LINES=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a decomposition topology line.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    LINE%ADJACENT_LINES=0
    LINE%BOUNDARY_LINE=.FALSE.

    EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the lines in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,element_idx,surrounding_element_idx,basis_local_line_idx, &
      & surrounding_element_basis_local_line_idx,element_local_node_idx,basis_local_line_node_idx,derivative_idx,version_idx, &
      & local_line_idx,surrounding_element_local_line_idx,node_idx,local_node_idx,elem_idx,line_end_node_idx,basis_node_idx, &
      & NODES_IN_LINE(4),NUMBER_OF_LINES,MAX_NUMBER_OF_LINES,NEW_MAX_NUMBER_OF_LINES,LINE_NUMBER,COUNT,APPROX_DIMENSION
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_LINES(:)
    INTEGER(INTG), POINTER :: TEMP_LINES(:,:),NEW_TEMP_LINES(:,:)
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_LINE_TYPE), POINTER :: DECOMPOSITION_LINE,DECOMPOSITION_LINE2
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: DECOMPOSITION_LINES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE,DOMAIN_LINE2
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_LINES)
    NULLIFY(NEW_TEMP_LINES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_LINES=>TOPOLOGY%LINES
      IF(ASSOCIATED(DECOMPOSITION_LINES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish line
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Guestimate the number of lines
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      MAX_NUMBER_OF_LINES=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    CASE(2)
                      APPROX_DIMENSION=SQRT(REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP))
                      !This should give the maximum and will over estimate the number of lines for a "square mesh" by approx 33%
                      MAX_NUMBER_OF_LINES=NINT(3.0_DP*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE(3)
                      !This should give the maximum and will over estimate the number of lines for a "cube mesh" by approx 73%
                      APPROX_DIMENSION=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS**(1.0_DP/3.0_DP)
                      MAX_NUMBER_OF_LINES=NINT(11.0_DP*APPROX_DIMENSION*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain.",ERR,ERROR,*999)
                    END SELECT
                    DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES
                    IF(ASSOCIATED(DOMAIN_LINES)) THEN
                      ALLOCATE(TEMP_LINES(4,MAX_NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate temporary lines array.",ERR,ERROR,*999)
                      ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                      NODES_NUMBER_OF_LINES=0
                      NUMBER_OF_LINES=0
                      TEMP_LINES=0
                      !Loop over the elements in the topology
                      DO element_idx=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_LINES(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate element element lines.",ERR,ERROR,*999)
                        !Loop over the local lines of the element
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          !Calculate the topology node numbers that make up the line
                          NODES_IN_LINE=0
                          DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                            NODES_IN_LINE(basis_local_line_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                              & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                          ENDDO !basis_local_line_node_idx
                          !Try and find a previously created line that matches in the adjacent elements
                          FOUND=.FALSE.
                          node_idx=NODES_IN_LINE(1)
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            IF(surrounding_element_idx/=element_idx) THEN
                              IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_LINES)) THEN
                                BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                DO surrounding_element_basis_local_line_idx=1,BASIS2%NUMBER_OF_LOCAL_LINES
                                  local_line_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)% &
                                    & ELEMENT_LINES(surrounding_element_basis_local_line_idx)
                                  IF(ALL(NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx))== &
                                    & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx),local_line_idx))) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                  ENDIF
                                ENDDO !surrounding_element_basis_local_line_idx
                                IF(FOUND) EXIT
                              ENDIF
                            ENDIF
                          ENDDO !elem_idx
                          IF(FOUND) THEN
                            !Line has already been created
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=local_line_idx
                          ELSE
                            !Line has not been created
                            IF(NUMBER_OF_LINES==MAX_NUMBER_OF_LINES) THEN
                              !We are at maximum. Reallocate the LINES array to be 20% bigger and try again.
                              NEW_MAX_NUMBER_OF_LINES=NINT(1.20_DP*MAX_NUMBER_OF_LINES,INTG)
                              ALLOCATE(NEW_TEMP_LINES(4,NEW_MAX_NUMBER_OF_LINES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate new number of lines.",ERR,ERROR,*999)
                              NEW_TEMP_LINES(:,1:NUMBER_OF_LINES)=TEMP_LINES(:,1:NUMBER_OF_LINES)
                              NEW_TEMP_LINES(:,NUMBER_OF_LINES+1:NEW_MAX_NUMBER_OF_LINES)=0
                              DEALLOCATE(TEMP_LINES)
                              TEMP_LINES=>NEW_TEMP_LINES
                              NULLIFY(NEW_TEMP_LINES)
                              MAX_NUMBER_OF_LINES=NEW_MAX_NUMBER_OF_LINES
                            ENDIF
                            NUMBER_OF_LINES=NUMBER_OF_LINES+1
                            TEMP_LINES(:,NUMBER_OF_LINES)=NODES_IN_LINE
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=NUMBER_OF_LINES
                            DO basis_local_line_node_idx=1,SIZE(NODES_IN_LINE,1)
                              IF(NODES_IN_LINE(basis_local_line_node_idx)/=0) &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))= &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))+1
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      !Allocate the line arrays and set them from the LINES and NODE_LINES arrays
                      DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                        ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate node lines array.",ERR,ERROR,*999)
                        DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                      ENDDO !node_idx
                      DEALLOCATE(NODES_NUMBER_OF_LINES)
                      ALLOCATE(DECOMPOSITION_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology lines.",ERR,ERROR,*999)
                      DECOMPOSITION_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      ALLOCATE(DOMAIN_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate domain topology lines.",ERR,ERROR,*999)
                      DOMAIN_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                        CALL DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(DECOMPOSITION_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        DO basis_local_line_node_idx=1,SIZE(TEMP_LINES,1)
                          IF(TEMP_LINES(basis_local_line_node_idx,local_line_idx)/=0) THEN
                            node_idx=TEMP_LINES(basis_local_line_node_idx,local_line_idx)
                            DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES+1
                            DOMAIN_NODES%NODES(node_idx)%NODE_LINES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES)= &
                              & local_line_idx
                          ENDIF
                        ENDDO !basis_local_line_node_idx  
                      ENDDO !local_line_idx
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DOMAIN_LINE=>DOMAIN_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          IF(.NOT.ASSOCIATED(DOMAIN_LINE%BASIS)) THEN
                            DECOMPOSITION_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%ELEMENT_NUMBER=element_idx !Needs checking
                            DECOMPOSITION_LINE%XI_DIRECTION=BASIS%LOCAL_LINE_XI_DIRECTION(basis_local_line_idx)
                            DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                            ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line nodes in line.",ERR,ERROR,*999)
                            ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                              & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line derivatives in line.",ERR,ERROR,*999)
                            DOMAIN_LINE%NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx))= &
                              & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx),LINE_NUMBER)
                            DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                              !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                              version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                & basis_local_line_node_idx,basis_local_line_idx))
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                              IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                  & BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx), &
                                  & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                  & basis_local_line_node_idx,basis_local_line_idx),BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                  & basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
                              ENDIF
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      DEALLOCATE(TEMP_LINES)
                      !Calculate adjacent lines and the surrounding elements for each line
                      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                        BASIS=>DOMAIN_LINE%BASIS
                        IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                          DECOMPOSITION_LINE%BOUNDARY_LINE=.TRUE.
                          DOMAIN_LINE%BOUNDARY_LINE=.TRUE.
                        ENDIF
                        !Allocate the elements surrounding the line
                        ALLOCATE(DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line surrounding elements.",ERR,ERROR,*999)
                        ALLOCATE(DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line element lines.",ERR,ERROR,*999)
                        DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
                        DECOMPOSITION_LINE%ADJACENT_LINES=0
                        !Loop over the nodes at each end of the line
                        DO line_end_node_idx=0,1
                          FOUND=.FALSE.
                          node_idx=DOMAIN_LINE%NODES_IN_LINE(line_end_node_idx*(BASIS%NUMBER_OF_NODES-1)+1)
                          !Loop over the elements surrounding the node.
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                            DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                            !Loop over the local lines of the element
                            DO basis_local_line_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_LINES
                              surrounding_element_local_line_idx=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                              IF(surrounding_element_local_line_idx/=local_line_idx) THEN
                                DECOMPOSITION_LINE2=>DECOMPOSITION_LINES%LINES(surrounding_element_local_line_idx)
                                DOMAIN_LINE2=>DOMAIN_LINES%LINES(surrounding_element_local_line_idx)
                                IF(DECOMPOSITION_LINE2%XI_DIRECTION==DECOMPOSITION_LINE%XI_DIRECTION) THEN
                                  !Lines run in the same direction.
                                  BASIS2=>DOMAIN_LINE2%BASIS
                                  IF(line_end_node_idx==0) THEN
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES)
                                  ELSE
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(1)
                                  ENDIF
                                  IF(local_node_idx==node_idx) THEN
                                    !The node at the 'other' end of this line matches the node at the current end of the line.
                                    !Check it is not a coexistant line running the other way
                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
                                      COUNT=0
                                      DO basis_node_idx=1,BASIS%NUMBER_OF_NODES
                                        IF(DOMAIN_LINE2%NODES_IN_LINE(basis_node_idx)== &
                                          & DOMAIN_LINE%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES-basis_node_idx+1)) &
                                          & COUNT=COUNT+1
                                      ENDDO !basis_node_idx
                                      IF(COUNT<BASIS%NUMBER_OF_NODES) THEN
                                        FOUND=.TRUE.
                                        EXIT
                                      ENDIF
                                    ELSE
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDIF
                            ENDDO !basis_local_line_idx
                            IF(FOUND) EXIT
                          ENDDO !element_idx
                          IF(FOUND) DECOMPOSITION_LINE%ADJACENT_LINES(line_end_node_idx)=surrounding_element_local_line_idx
                        ENDDO !line_end_node_idx
                      ENDDO !local_line_idx
                      !Set the surrounding elements
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=element_idx
                          DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=basis_local_line_idx
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                    ELSE
                      CALL FlagError("Domain topology lines is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain lines
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES                      
                          IF(ASSOCIATED(DOMAIN_LINES)) THEN
                            ALLOCATE(DOMAIN_LINES%LINES(DECOMPOSITION_LINES%NUMBER_OF_LINES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain lines lines.",ERR,ERROR,*999)
                            DOMAIN_LINES%NUMBER_OF_LINES=DECOMPOSITION_LINES%NUMBER_OF_LINES
                            ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                            NODES_NUMBER_OF_LINES=0
                            !Loop over the lines in the topology
                            DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                              DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                element_idx=DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(1)
                                basis_local_line_idx=DECOMPOSITION_LINE%ELEMENT_LINES(1)
                                CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                                DOMAIN_LINE%NUMBER=local_line_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_LINE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                                DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                                ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in line.",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in line.",ERR,ERROR,*999)
                                DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx, &
                                    & basis_local_line_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                    & basis_local_line_node_idx,basis_local_line_idx))
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                                  IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx),element_local_node_idx)
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx),BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx))
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
                                  ENDIF
                                  NODES_NUMBER_OF_LINES(node_idx)=NODES_NUMBER_OF_LINES(node_idx)+1
                                ENDDO !basis_local_line_node_idx
                              ELSE
                                CALL FlagError("Line is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !local_line_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node lines.",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_LINES)
                            DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              BASIS=>DOMAIN_LINE%BASIS
                              DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_LINES=DOMAIN_NODE%NUMBER_OF_NODE_LINES+1
                                DOMAIN_NODE%NODE_LINES(DOMAIN_NODE%NUMBER_OF_NODE_LINES)=local_line_idx
                              ENDDO !basis_local_line_node_idx
                            ENDDO !local_line_idx
                          ELSE
                            CALL FlagError("Domain lines is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology lines is not associated.",ERR,ERROR,*999)

      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology lines:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of lines = ",DECOMPOSITION_LINES%NUMBER_OF_LINES,ERR,ERROR,*999)
      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Line number = ",DECOMPOSITION_LINE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",DECOMPOSITION_LINE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%ELEMENT_LINES,'("      Element lines        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_LINE%ADJACENT_LINES, &
          & '("      Adjacent lines       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary line = ",DECOMPOSITION_LINE%BOUNDARY_LINE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_LINE=>DOMAIN%TOPOLOGY%LINES%LINES(local_line_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_LINE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_LINE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in lines = ",DOMAIN_LINE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_LINE%NODES_IN_LINE, &
            & '("        Nodes in line        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_line_node_idx=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_line_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(1,:,basis_local_line_node_idx),'("            Derivatives in line  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(2,:,basis_local_line_node_idx), &
              & '("            Derivatives Versions in line  :",4(X,I8))','(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_line_node_idx
        ENDDO !component_idx
      ENDDO !local_line_idx
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_LINES)) DEALLOCATE(TEMP_LINES)
    IF(ASSOCIATED(NEW_TEMP_LINES)) DEALLOCATE(NEW_TEMP_LINES)
    IF(ALLOCATED(NODES_NUMBER_OF_LINES)) DEALLOCATE(NODES_NUMBER_OF_LINES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given decomposition topology. \todo Pass in the topology lines
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DECOMPOSITION_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines.",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_DataPointsInitialise(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DecompositionTopology_DataPointsInitialise",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Decomposition already has topology data points associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points.",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%numberOfDataPoints=0
        TOPOLOGY%dataPoints%totalNumberOfDataPoints=0
        TOPOLOGY%dataPoints%numberOfGlobalDataPoints=0
        NULLIFY(TOPOLOGY%dataPoints%dataPointsTree)
        TOPOLOGY%dataPoints%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION      
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsInitialise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsInitialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataPointsInitialise
  
  !
  !================================================================================================================================
  !
 
  !>Finalises a face in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(FACE%SURROUNDING_ELEMENTS)) DEALLOCATE(FACE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(FACE%ELEMENT_FACES)) DEALLOCATE(FACE%ELEMENT_FACES)
!    FACE%ADJACENT_FACES=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a decomposition topology face.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!    FACE%ADJACENT_FACES=0
    FACE%BOUNDARY_FACE=.FALSE.
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,ne,surrounding_element_idx,basis_local_face_idx,surrounding_element_basis_local_face_idx, &
      & element_local_node_idx,basis_local_face_node_idx,basis_local_face_derivative_idx,derivative_idx,version_idx,face_idx, &
      & node_idx,elem_idx,NODES_IN_FACE(16),NUMBER_OF_FACES,MAX_NUMBER_OF_FACES,NEW_MAX_NUMBER_OF_FACES,FACE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_FACES(:)
    INTEGER(INTG), POINTER :: TEMP_FACES(:,:),NEW_TEMP_FACES(:,:)
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMPOSITION_FACE!,DECOMPOSITION_FACE2
    TYPE(DECOMPOSITION_FACES_TYPE), POINTER :: DECOMPOSITION_FACES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE!,DOMAIN_FACE2
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_FACES)
    NULLIFY(NEW_TEMP_FACES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_FACES=>TOPOLOGY%FACES
      IF(ASSOCIATED(DECOMPOSITION_FACES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish face
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Estimate the number of faces
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      ! Faces not calculated in 1D 
                    CASE(2)
                      ! Faces not calculated in 2D
                    CASE(3)
                      !This should give the maximum and will over estimate the number of faces for a "cube mesh" by approx 33%
                      MAX_NUMBER_OF_FACES= &
                        & NINT((DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS*5.0_DP+1.0_DP)*(4.0_DP/3.0_DP),INTG)

                      DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                      IF(ASSOCIATED(DOMAIN_FACES)) THEN
                        ALLOCATE(TEMP_FACES(16,MAX_NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate temporary faces array",ERR,ERROR,*999)
                        ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                        NODES_NUMBER_OF_FACES=0
                        NUMBER_OF_FACES=0
                        TEMP_FACES=0
                        !Loop over the elements in the topology and fill temp_faces with node numbers for each element
                        DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_FACES(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate element faces of element",ERR,ERROR,*999)
                          !Loop over the local faces of the element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            !Calculate the topology node numbers that make up the face
                            NODES_IN_FACE=0
                            !Check whether face has already been read out
                            DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                              !Read out node numbers of local face from ELEMENT_NODES
                              NODES_IN_FACE(basis_local_face_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                                & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                            ENDDO !basis_local_face_node_idx
                            !Try and find a previously created face that matches in the adjacent elements
                            FOUND=.FALSE.
                            node_idx=NODES_IN_FACE(1)
                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                              surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                              IF(surrounding_element_idx/=ne) THEN
                                IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES)) THEN
                                  BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                  DO surrounding_element_basis_local_face_idx=1,BASIS2%NUMBER_OF_LOCAL_FACES
                                    face_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES( &
                                      & surrounding_element_basis_local_face_idx)
                                    IF(ALL(NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx))== &
                                      & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx),face_idx))) THEN
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !surrounding_element_basis_local_face_idx
                                  IF(FOUND) EXIT
                                ENDIF
                              ENDIF
                            ENDDO !elem_idx
                            IF(FOUND) THEN
                              !Face has already been created
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=face_idx
                            ELSE
                              !Face has not been created
                              IF(NUMBER_OF_FACES==MAX_NUMBER_OF_FACES) THEN
                                !We are at maximum. Reallocate the FACES array to be 20% bigger and try again.
                                NEW_MAX_NUMBER_OF_FACES=NINT(1.20_DP*MAX_NUMBER_OF_FACES,INTG)
                                !\todo: Change 16 to a variable and above for NODES_IN_FACE
                                ALLOCATE(NEW_TEMP_FACES(16,NEW_MAX_NUMBER_OF_FACES),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate new number of faces",ERR,ERROR,*999)
                                NEW_TEMP_FACES(:,1:NUMBER_OF_FACES)=TEMP_FACES(:,1:NUMBER_OF_FACES)
                                NEW_TEMP_FACES(:,NUMBER_OF_FACES+1:NEW_MAX_NUMBER_OF_FACES)=0
                                DEALLOCATE(TEMP_FACES)
                                TEMP_FACES=>NEW_TEMP_FACES
                                NULLIFY(NEW_TEMP_FACES)
                                MAX_NUMBER_OF_FACES=NEW_MAX_NUMBER_OF_FACES
                              ENDIF
                              NUMBER_OF_FACES=NUMBER_OF_FACES+1
                              TEMP_FACES(:,NUMBER_OF_FACES)=NODES_IN_FACE(:)
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=NUMBER_OF_FACES
                              DO basis_local_face_node_idx=1,SIZE(NODES_IN_FACE,1)
                                IF(NODES_IN_FACE(basis_local_face_node_idx)/=0) &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))= &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))+1
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        !Allocate the face arrays and set them from the FACES and NODE_FACES arrays
                        DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                          ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate node faces array",ERR,ERROR,*999)
                          DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                        ENDDO !node_idx
                        DEALLOCATE(NODES_NUMBER_OF_FACES)
                        ALLOCATE(DECOMPOSITION_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology faces",ERR,ERROR,*999)
                        DECOMPOSITION_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        ALLOCATE(DOMAIN_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate domain topology faces",ERR,ERROR,*999)
                        DOMAIN_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                          CALL DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(DECOMPOSITION_FACES%FACES(face_idx),ERR,ERROR,*999)
                          CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                          DO basis_local_face_node_idx=1,SIZE(TEMP_FACES,1)
                            IF(TEMP_FACES(basis_local_face_node_idx,face_idx)/=0) THEN
                              node_idx=TEMP_FACES(basis_local_face_node_idx,face_idx)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES+1
                              DOMAIN_NODES%NODES(node_idx)%NODE_FACES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES)=face_idx
                            ENDIF
                          ENDDO !basis_local_face_node_idx  
                        ENDDO !face_idx

                        !Set nodes in face and derivatives of nodes in face for domain faces
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          !Loop over local faces of element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DOMAIN_FACE=>DOMAIN_FACES%FACES(FACE_NUMBER)
                            DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                            IF(.NOT.ASSOCIATED(DOMAIN_FACE%BASIS)) THEN
                              DECOMPOSITION_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%ELEMENT_NUMBER=ne !! Needs checking
!                              DECOMPOSITION_FACE%ELEMENT_NUMBER=DECOMPOSITION_ELEMENT%NUMBER
!                              DOMAIN_FACE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                              DECOMPOSITION_FACE%XI_DIRECTION=BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx)
                              DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(ABS(DECOMPOSITION_FACE%XI_DIRECTION))%PTR
                              ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                & STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face nodes in face",ERR,ERROR,*999)
                              ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face derivatives in face",ERR,ERROR,*999)
                              DOMAIN_FACE%DERIVATIVES_IN_FACE=0
                              !Set nodes in face based upon face number
                              DOMAIN_FACE%NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx))= &
                                & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx),FACE_NUMBER)
                              !Set derivatives of nodes in domain face from derivatives of nodes in element
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx, &
                                  & basis_local_face_idx)
                                !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_FACE( &
                                  & basis_local_face_node_idx,basis_local_face_idx))
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                  DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(1,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(2,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=version_idx
                                  ENDDO !basis_local_face_derivative_idx
                                ENDIF
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        DEALLOCATE(TEMP_FACES)
                        !\todo Note: Adjacency will be left out of faces calculation for the time being
                        !Calculate adjacent faces and the surrounding elements for each face
                        DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                          DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                          DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                          BASIS=>DOMAIN_FACE%BASIS
                          IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                            DECOMPOSITION_FACE%BOUNDARY_FACE=.TRUE.
                            DOMAIN_FACE%BOUNDARY_FACE=.TRUE.
                          ENDIF
                          !Allocate the elements surrounding the face
                          ALLOCATE(DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face surrounding elements",ERR,ERROR,*999)

                          ALLOCATE(DECOMPOSITION_FACE%ELEMENT_FACES(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face element faces",ERR,ERROR,*999)
!                          DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!                          DECOMPOSITION_FACE%ADJACENT_FACES=0

                           !Loop over the nodes at each end of the face
!                          DO node_idx1=0,1
!                           DO node_idx2=0,1
!                            FOUND=.FALSE.
!                            node_idx=DOMAIN_FACE%NODES_IN_FACE((node_idx2*BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION*(BASIS%NUMBER_OF_FACES-1))&
!                                                                             &+(node_idx1*(BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
                             !Loop over the elements surrounding the node.
!                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
!                              ne=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
!                              DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
!                              DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                               !Loop over the local faces of the element
!                              DO basis_local_face_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_FACES
!                                nf2=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
!                                IF(nf2/=face_idx) THEN
!                                  DECOMPOSITION_FACE2=>DECOMPOSITION_FACES%FACES(nf2)
!                                  DOMAIN_FACE2=>DOMAIN_FACES%FACES(nf2)
                                   !Check whether XI of face have same direction
!                                  IF ((OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),2,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),2,1)).OR.&
!                                     &(OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),3,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),3,1))) THEN
                                     !Loop over nodes in face of surrounding element
!                                    BASIS2=>DOMAIN_FACE2%BASIS
!                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
!                                      NODE_COUNT=0
!                                      DO node_idx3=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                        DO node_idx4=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                          np2=DOMAIN_FACE2%NODES_IN_FACE((node_idx4*(BASIS2%NUMBER_OF_FACES-1))&
!                                                                      &+(node_idx3*(BASIS2%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
!                                          IF(np2==node_idx) NODE_COUNT=NODE_COUNT+1
!                                        ENDDO !node_idx4
!                                      ENDDO !node_idx3
!                                      IF(NODE_COUNT<BASIS%NUMBER_OF_NODES) THEN
!                                        FOUND=.TRUE.
!                                        EXIT
!                                      ENDIF
!                                    ENDIF
!                                  ENDIF
!                                ENDIF
!                              ENDDO !basis_local_face_idx
!                                IF(FOUND) EXIT
!                            ENDDO !elem_idx
!                            IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx2)=nf2
!                           ENDDO !node_idx2
!                           IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx1)=nf2
!                          ENDDO !node_idx1
                        ENDDO !face_idx

                        !Set the surrounding elements
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DO face_idx=1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS
                              DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(face_idx)=ne
                              DECOMPOSITION_FACE%ELEMENT_FACES(face_idx)=basis_local_face_idx
                            ENDDO
                          ENDDO !basis_local_face_idx
                        ENDDO !ne
                      ELSE
                        CALL FlagError("Domain topology faces is not associated",ERR,ERROR,*999)
                      ENDIF
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain",ERR,ERROR,*999)
                    END SELECT
                 ELSE
                    CALL FlagError("Domain topology elements is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain faces
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                          IF(ASSOCIATED(DOMAIN_FACES)) THEN
                            ALLOCATE(DOMAIN_FACES%FACES(DECOMPOSITION_FACES%NUMBER_OF_FACES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain faces faces",ERR,ERROR,*999)
                            DOMAIN_FACES%NUMBER_OF_FACES=DECOMPOSITION_FACES%NUMBER_OF_FACES
                            ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                            NODES_NUMBER_OF_FACES=0
                            !Loop over the faces in the topology
                            DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                              DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                ne=DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(1)
                                basis_local_face_idx=DECOMPOSITION_FACE%ELEMENT_FACES(1)
                                CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                                DOMAIN_FACE%NUMBER=face_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(ABS(DECOMPOSITION_FACE%XI_DIRECTION))%PTR
                                ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in face",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in face",ERR,ERROR,*999)
                                !Set derivatives of nodes in domain face from derivatives of nodes in element
                                DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx, &
                                    & basis_local_face_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_FACE( &
                                    & basis_local_face_node_idx,basis_local_face_idx))
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                  IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                      derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%DERIVATIVES_IN_FACE(1,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=derivative_idx
                                      version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%DERIVATIVES_IN_FACE(2,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=version_idx
                                    ENDDO !basis_local_face_derivative_idx
                                  ENDIF
                                  NODES_NUMBER_OF_FACES(node_idx)=NODES_NUMBER_OF_FACES(node_idx)+1
                                ENDDO !basis_local_face_node_idx
                              ELSE
                                CALL FlagError("Face is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !face_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node faces",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_FACES)
                            DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              BASIS=>DOMAIN_FACE%BASIS
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_FACES=DOMAIN_NODE%NUMBER_OF_NODE_FACES+1
                                !Set the face numbers a node is on
                                DOMAIN_NODE%NODE_FACES(DOMAIN_NODE%NUMBER_OF_NODE_FACES)=face_idx
                              ENDDO !basis_local_face_node_idx
                            ENDDO !face_idx
                          ELSE
                            CALL FlagError("Domain faces is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FlagError("Topology decomposition is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology faces is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology faces:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of faces = ",DECOMPOSITION_FACES%NUMBER_OF_FACES,ERR,ERROR,*999)
      DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
        DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
        DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Face number = ",DECOMPOSITION_FACE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction (Normal to Face) = &
                                                                         &",DECOMPOSITION_FACE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%ELEMENT_FACES,'("      Element faces        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
!        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_FACE%ADJACENT_FACES, &
!          & '("      Adjacent faces       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary face = ",DECOMPOSITION_FACE%BOUNDARY_FACE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_FACE=>DOMAIN%TOPOLOGY%FACES%FACES(face_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_FACE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_FACE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in faces = ",DOMAIN_FACE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_FACE%NODES_IN_FACE, &
            & '("        Nodes in face        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_face_node_idx=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_face_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(1,:,basis_local_face_node_idx),'("            Derivatives in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(2,:,basis_local_face_node_idx),'("            Derivatives Versions in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_face_node_idx
        ENDDO !component_idx
      ENDDO !face_idx
    ENDIF

    EXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_FACES)) DEALLOCATE(TEMP_FACES)
    IF(ASSOCIATED(NEW_TEMP_FACES)) DEALLOCATE(NEW_TEMP_FACES)
    IF(ALLOCATED(NODES_NUMBER_OF_FACES)) DEALLOCATE(NODES_NUMBER_OF_FACES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DECOMPOSITION_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the decomposition type for a decomposition. \see OPENCMISS::CMISSDecompositionTypeGet
  SUBROUTINE DECOMPOSITION_TYPE_GET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the type for
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the decomposition type for the specified decomposition \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        TYPE=DECOMPOSITION%DECOMPOSITION_TYPE
      ELSE
        CALL FlagError("Decomposition has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_GET",ERR,ERROR)
    RETURN
  END SUBROUTINE DECOMPOSITION_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the decomposition type for a decomposition.  \see OPENCMISS::CMISSDecompositionTypeSet
  SUBROUTINE DECOMPOSITION_TYPE_SET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The decomposition type to set \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          !heye: three types for decomposition--decompostion_all_type means no decomposition
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
        CASE(DECOMPOSITION_CALCULATED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_CALCULATED_TYPE
        CASE(DECOMPOSITION_USER_DEFINED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_USER_DEFINED_TYPE
        CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether lines should be calculated in the the decomposition. \see OPENCMISS::CMISSDecompositionCalculateLinesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET(DECOMPOSITION,CALCULATE_LINES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_LINES_FLAG !<The boolean flag to determine whether the lines should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_LINES=CALCULATE_LINES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_LINES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether faces should be calculated in the the decomposition. \see OPENCMISS::CMISSDecompositionCalculateFacesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET(DECOMPOSITION,CALCULATE_FACES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_FACES_FLAG !<The boolean flag to determine whether the faces should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_FACES=CALCULATE_FACES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_FACES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in DECOMPOSITION a pointer to the decomposition identified by USER_NUMBER in the given MESH. If no decomposition with that USER_NUMBER exists DECOMPOSITION is left nullified.
  SUBROUTINE DECOMPOSITION_USER_NUMBER_FIND(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to find
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh containing the decomposition to find
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On return a pointer to the decomposition with the specified user number. If no decomposition with that user number exists then DECOMPOSITION is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_USER_NUMBER_FIND",ERR,ERROR,*999)

    NULLIFY(DECOMPOSITION)
    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        decomposition_idx=1
        DO WHILE(decomposition_idx<=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.ASSOCIATED(DECOMPOSITION))
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            DECOMPOSITION=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
          ELSE
            decomposition_idx=decomposition_idx+1
          ENDIF
        ENDDO
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_USER_NUMBER_FIND")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_USER_NUMBER_FIND",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises the domain decompositions for a given mesh. \todo Pass in a pointer to the decompositions.
  SUBROUTINE DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise the decomposition for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        DO WHILE(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>0)
          CALL DECOMPOSITION_DESTROY(MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR,ERR,ERROR,*999)
        ENDDO !no_decomposition
       DEALLOCATE(MESH%DECOMPOSITIONS)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain decompositions for a given mesh.
  SUBROUTINE DECOMPOSITIONS_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise the decompositions for

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        CALL FlagError("Mesh already has decompositions associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(MESH%DECOMPOSITIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Mesh decompositions could not be allocated.",ERR,ERROR,*999)
        MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        NULLIFY(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
        MESH%DECOMPOSITIONS%MESH=>MESH
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the domain for a given decomposition and deallocates all memory. \todo Pass in a pointer to the domain
  SUBROUTINE DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the domain for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    
    ENTERS("DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          DO component_idx=1,DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS
            CALL DOMAIN_MAPPINGS_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)        
            CALL DOMAIN_TOPOLOGY_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR)
          ENDDO !component_idx
          DEALLOCATE(DECOMPOSITION%DOMAIN)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain for a given decomposition.
  SUBROUTINE DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the domain for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx

    ENTERS("DOMAIN_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          CALL FlagError("Decomposition already has a domain associated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Decomposition domain could not be allocated.",ERR,ERROR,*999)
          DO component_idx=1,DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS !Mesh component
            ALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Decomposition domain component could not be allocated.",ERR,ERROR,*999)
            DECOMPOSITION%DOMAIN(component_idx)%PTR%DECOMPOSITION=>DECOMPOSITION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH=>DECOMPOSITION%MESH
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH_COMPONENT_NUMBER=component_idx
            DECOMPOSITION%DOMAIN(component_idx)%PTR%REGION=>DECOMPOSITION%MESH%REGION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_DIMENSIONS=DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%MAPPINGS)
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%TOPOLOGY)
            CALL DOMAIN_MAPPINGS_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            CALL DOMAIN_TOPOLOGY_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
          ENDDO !component_idx
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the mappings in the given domain. \todo pass in the domain mappings
  SUBROUTINE DOMAIN_MAPPINGS_FINALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%MAPPINGS)
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the element mapping in the given domain mapping. \todo pass in the domain mappings elements
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the elements for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%ELEMENTS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL FlagError("Domain elements mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings elements.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%ELEMENTS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the mappings for a domain decomposition. \todo finalise on error.
  SUBROUTINE DOMAIN_MAPPINGS_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to initialise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        CALL FlagError("Domain already has mappings associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN%MAPPINGS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings.",ERR,ERROR,*999)
        DOMAIN%MAPPINGS%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%MAPPINGS%ELEMENTS)
        NULLIFY(DOMAIN%MAPPINGS%NODES)
        !Calculate the node and element mappings
        CALL DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the node mapping in the given domain mappings. \todo pass in the nodes mapping
  SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mapping to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%NODES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the node mapping in the given domain mapping. \todo finalise on error
  SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL FlagError("Domain nodes mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings nodes.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%NODES,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the local domain topology and domain mappings.
  SUBROUTINE Domain_Calculate(decomposition,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<A pointer to the decomposition to calculate the domain for. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: found
    INTEGER(INTG) :: derivativeIdx,nodeIdx,numberOfDomains,myDomain,domainOffset,numberOfLocalElements,globalElementNumber, &
      & numberOfAdjacentDomains,adjacentDomainIdx,numberOfLocal,numberOfGhost,numberOfGhostElements, &
      & totalNumberOfLocalElements,domainIdx,domainNumber,localElementCount,ghostElementCount,receiveBufferIdx, &
      & elementNumberOfAdjacentDomains,elementDomainIdx,elementOwner,elementIdx,meshComponentIdx,elementNodeIdx, &
      & localNodeCount,localNodeNumber,globalNodeNumber,userNodeNumber,insertStatus,nodeDerivativeIdx,dummyErr,MPI_IERROR
    INTEGER(INTG), ALLOCATABLE :: sendCounts(:),receiveCounts(:),requests(:),statuses(:,:),adjacentDomainNumbers(:), &
      & globalElementNumbers(:),userNodeNumbers(:),elementNodeOffsets(:),elementNodes(:),userElementNumbers(:)
    TYPE(INTEGER_INTG_ALLOC_TYPE), ALLOCATABLE :: sendBuffers(:),receiveBuffers(:),elementDomainsArrays(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DecompositionElementDomainsType), POINTER :: elementDomains
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements 
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(LIST_TYPE), POINTER :: adjacentDomainsList,userNodeNumbersList,ghostUserNodeNumbersList,elementNodesList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: sendBufferLists(:),elementGhostSendLists(:),elementGhostReceiveLists(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping,nodesMapping
    TYPE(MeshElementsType), POINTER ::  meshElements
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(MESH_TYPE), POINTER ::  mesh
    TYPE(TREE_TYPE), POINTER :: elementsTree,nodesTree
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Domain_Calculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    mesh=>decomposition%mesh
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Decomposition mesh is not associated.",err,error,*999)
    myDomain=COMPUTATIONAL_NODE_NUMBER_GET(err,error)+1
    IF(err/=0) GOTO 999
    numberOfDomains=decomposition%NUMBER_OF_DOMAINS
    meshElements=>mesh%topology(decomposition%MESH_COMPONENT_NUMBER)%ptr%elements
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
    elementDomains=>decomposition%elementDomains

    !Redistribute the mesh by sending the relevant data to their new domains.

    !The send buffer for a domain should contain 1)number of local elements, 2)number of ghost elements and further 
    !for each element contain a)the global element number, b)the number of adjacent domains, c)the adjacent domains.

    ALLOCATE(sendBufferLists(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate send buffer lists.",err,error,*999)
    
    DO domainIdx=1,numberOfDomains
      NULLIFY(sendBufferLists(domainIdx)%ptr)
      CALL List_CreateStart(sendBufferLists(domainIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(sendBufferLists(domainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(sendBufferLists(domainIdx)%ptr,2+2*COUNT(elementDomains%domains==domainIdx),err,error,*999)
      CALL List_MutableSet(sendBufferLists(domainIdx)%ptr,.TRUE.,err,error,*999)
      CALL List_CreateFinish(sendBufferLists(domainIdx)%ptr,err,error,*999)
      !The number of local elements on the new domain.
      CALL List_ItemAdd(sendBufferLists(domainIdx)%ptr,0,err,error,*999)
      !The number of ghost elements on the new domain.
      CALL List_ItemAdd(sendBufferLists(domainIdx)%ptr,0,err,error,*999)
    END DO !domainIdx
      
    !See Decomposition_ElementDomainsCalculate for how the elements were initially distributed.
    domainOffset=(mesh%NUMBER_OF_ELEMENTS*(myDomain-1))/numberOfDomains+1
    numberOfLocalElements=(mesh%NUMBER_OF_ELEMENTS*myDomain)/numberOfDomains- &
      & (mesh%NUMBER_OF_ELEMENTS*(myDomain-1))/numberOfDomains
    
    DO elementIdx=1,numberOfLocalElements
      !Note: the element owner is also counted as adjacent domain.
      numberOfAdjacentDomains=elementDomains%offsets(elementIdx+1)-elementDomains%offsets(elementIdx)
      !The element is a local element on the new domain.
      domainNumber=elementDomains%domains(elementDomains%offsets(elementIdx))
      globalElementNumber=domainOffset+elementIdx-1
      CALL List_ItemAdd(sendBufferLists(domainNumber)%ptr,globalElementNumber,err,error,*999)
      CALL List_ItemAdd(sendBufferLists(domainNumber)%ptr,numberOfAdjacentDomains,err,error,*999)
      DO elementDomainIdx=elementDomains%offsets(elementIdx),elementDomains%offsets(elementIdx+1)-1
        CALL List_ItemAdd(sendBufferLists(domainNumber)%ptr,elementDomains%domains(elementDomainIdx),err,error,*999)
      END DO !elementDomainIdx
      CALL List_ItemGet(sendBufferLists(domainNumber)%ptr,1,numberOfLocal,err,error,*999)
      CALL List_ItemSet(sendBufferLists(domainNumber)%ptr,1,numberOfLocal+1,err,error,*999)
      !The element is a ghost element on the new domain.
      DO adjacentDomainIdx=2,numberOfAdjacentDomains
        domainNumber=elementDomains%domains(elementDomains%offsets(elementIdx)+adjacentDomainIdx-1)
        CALL List_ItemAdd(sendBufferLists(domainNumber)%ptr,globalElementNumber,err,error,*999)
        CALL List_ItemAdd(sendBufferLists(domainNumber)%ptr,numberOfAdjacentDomains,err,error,*999)
        DO elementDomainIdx=elementDomains%offsets(elementIdx),elementDomains%offsets(elementIdx+1)-1
          CALL List_ItemAdd(sendBufferLists(domainNumber)%ptr,elementDomains%domains(elementDomainIdx),err,error,*999)
        END DO !elementDomainIdx
        CALL List_ItemGet(sendBufferLists(domainNumber)%ptr,2,numberOfGhost,err,error,*999)
        CALL List_ItemSet(sendBufferLists(domainNumber)%ptr,2,numberOfGhost+1,err,error,*999)
      END DO !adjacentDomainIdx
    END DO !elementIdx

    !Get rid of the old element domains, so that we can calculate the new ones
    CALL Decomposition_ElementDomainsFinalise(elementDomains,err,error,*999)

    ALLOCATE(sendCounts(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate send buffers.",err,error,*999)
    ALLOCATE(sendBuffers(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate send buffers.",err,error,*999)

    ALLOCATE(receiveCounts(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate receive counts.",err,error,*999)
    ALLOCATE(receiveBuffers(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate receive buffers.",err,error,*999)
    
    DO domainIdx=1,numberOfDomains
      CALL List_DetachAndDestroy(sendBufferLists(domainIdx)%ptr,sendCounts(domainIdx), &
        & sendBuffers(domainIdx)%array,err,error,*999)
    END DO
    IF(ALLOCATED(sendBufferLists)) DEALLOCATE(sendBufferLists)

    !Collect the receive counts from all domains.
    CALL MPI_ALLTOALL(sendCounts,1,MPI_INTEGER,receiveCounts,1,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLTOALL",MPI_IERROR,err,error,*999)

    !Allocate all the receive buffers.
    DO domainIdx=1,numberOfDomains
      ALLOCATE(receiveBuffers(domainIdx)%array(receiveCounts(domainIdx)),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate receive buffer.",err,error,*999)
    END DO

    !We already know the receive buffer for this domain.
    receiveBuffers(myDomain)%array=sendBuffers(myDomain)%array
    ALLOCATE(statuses(MPI_STATUS_SIZE,2*numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate mpi statuses.",err,error,*999)
    ALLOCATE(requests(2*numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate mpi receive requests.",err,error,*999)
    requests=MPI_REQUEST_NULL

    !Get all the receive buffers from the other domains.
    DO domainIdx=1,numberOfDomains
      IF(domainIdx/=myDomain.AND.receiveCounts(domainIdx)>0) THEN
        CALL MPI_IRECV(receiveBuffers(domainIdx)%array,receiveCounts(domainIdx),MPI_INTEGER,domainIdx-1, &
          & MPI_ELEMENT_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,requests(domainIdx),MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,err,error,*999)
      END IF
    END DO
    DO domainIdx=1,numberOfDomains
      IF(domainIdx/=myDomain.AND.sendCounts(domainIdx)>0) THEN
        CALL MPI_ISEND(sendBuffers(domainIdx)%array,sendCounts(domainIdx),MPI_INTEGER,domainIdx-1, &
          & MPI_ELEMENT_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,requests(numberOfDomains+domainIdx),MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,err,error,*999)
      END IF
    END DO

    CALL MPI_WAITALL(2*numberOfDomains,requests,statuses,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,err,error,*999)

    IF(DIAGNOSTICS5) THEN
      DO domainIdx=1,numberOfDomains
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI IRECV call posted:",err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive count = ",receiveCounts(domainIdx),err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive datatype = ",MPI_INTEGER,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive source = ",domainIdx-1,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",MPI_ELEMENT_EXCHANGE_TAG,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive comm = ",COMPUTATIONAL_ENVIRONMENT%MPI_COMM,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive request = ",requests(domainIdx),err,error,*999)                
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,receiveCounts(domainIdx),8,8,receiveBuffers(domainIdx)%array, &
          & '("      Receive data       :",8(X,I10))','(39X,8(X,I10))',err,error,*999)      
      ENDDO !domainIdx
      DO domainIdx=1,numberOfDomains
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI ISEND call posted:",err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send count = ",sendCounts(domainIdx),err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send datatype = ",MPI_INTEGER,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send dest = ",domainIdx-1,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",MPI_ELEMENT_EXCHANGE_TAG,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send comm = ",COMPUTATIONAL_ENVIRONMENT%MPI_COMM,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send request = ",requests(numberOfDomains+domainIdx),err,error,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,sendCounts(domainIdx),8,8,sendBuffers(domainIdx)%array, &
          & '("      Send data       :",8(X,I10))','(39X,8(X,I10))',err,error,*999)      
      ENDDO !adjacentDomainIdx
    END IF

    IF(ALLOCATED(requests)) DEALLOCATE(requests)
    IF(ALLOCATED(statuses)) DEALLOCATE(statuses)

    !Count the number of local and ghost elements in the new distributed mesh.
    numberOfLocalElements=0
    numberOfGhostElements=0
    DO domainIdx=1,numberOfDomains
      IF(receiveCounts(domainIdx)>0) THEN
        numberOfLocalElements=numberOfLocalElements+receiveBuffers(domainIdx)%array(1)
        numberOfGhostElements=numberOfGhostElements+receiveBuffers(domainIdx)%array(2)
      END IF
    END DO
    totalNumberOfLocalElements=numberOfLocalElements+numberOfGhostElements

    NULLIFY(elementsTree) 
    CALL Tree_CreateStart(elementsTree,err,error,*999)
    CALL Tree_InsertTypeSet(elementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(elementsTree,err,error,*999)
    
    ALLOCATE(globalElementNumbers(totalNumberOfLocalElements),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate global element numbers.",err,error,*999)
    ALLOCATE(userElementNumbers(totalNumberOfLocalElements),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate user element numbers.",err,error,*999)
    
    ALLOCATE(elementDomainsArrays(totalNumberOfLocalElements),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate element domains arrays.",err,error,*999)
    
    !Unpack the received data
    localElementCount=0
    ghostElementCount=0
    DO domainIdx=1,numberOfDomains
      receiveBufferIdx=2
      DO WHILE(receiveBufferIdx<receiveCounts(domainIdx))
        receiveBufferIdx=receiveBufferIdx+1
        globalElementNumber=receiveBuffers(domainIdx)%array(receiveBufferIdx)
        receiveBufferIdx=receiveBufferIdx+1
        elementNumberOfAdjacentDomains=receiveBuffers(domainIdx)%array(receiveBufferIdx)
        receiveBufferIdx=receiveBufferIdx+1
        elementOwner=receiveBuffers(domainIdx)%array(receiveBufferIdx)
        IF(elementOwner==myDomain) THEN
          !We have a local element
          localElementCount=localElementCount+1
          elementIdx=localElementCount
        ELSE
          !We have a ghost element
          ghostElementCount=ghostElementCount+1
          elementIdx=ghostElementCount+numberOfLocalElements
        END IF
        globalElementNumbers(elementIdx)=globalElementNumber
        userElementNumbers(elementIdx)=meshElements%elements(globalElementNumber)%USER_NUMBER
        CALL Tree_ItemInsert(elementsTree,userElementNumbers(elementIdx),elementIdx,insertStatus,err,error,*999)
        ALLOCATE(elementDomainsArrays(elementIdx)%array(elementNumberOfAdjacentDomains),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate element domains array.",err,error,*999)
        !Store the adjacent domains for later use
        elementDomainsArrays(elementIdx)%array= &
          & receiveBuffers(domainIdx)%array(receiveBufferIdx:receiveBufferIdx+elementNumberOfAdjacentDomains-1)
        receiveBufferIdx=receiveBufferIdx+elementNumberOfAdjacentDomains-1
      END DO
    END DO !domainIdx

    !Deallocate the receive and send buffers.
    IF(ALLOCATED(receiveCounts)) DEALLOCATE(receiveCounts)
    IF(ALLOCATED(receiveBuffers)) THEN
      DO domainIdx=1,numberOfDomains
        IF(ALLOCATED(receiveBuffers(domainIdx)%array)) DEALLOCATE(receiveBuffers(domainIdx)%array)
      END DO
      DEALLOCATE(receiveBuffers)
    END IF
    IF(ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
    IF(ALLOCATED(sendBuffers)) THEN
      DO domainIdx=1,numberOfDomains
        IF(ALLOCATED(sendBuffers(domainIdx)%array)) DEALLOCATE(sendBuffers(domainIdx)%array)
      END DO
      DEALLOCATE(sendBuffers)
    END IF

    decomposition%topology%elements%ELEMENTS_TREE=>elementsTree
   
    !Calculate the new element domains.
    NULLIFY(decomposition%elementDomains)
    CALL Decomposition_ElementDomainsInitialise(decomposition%elementDomains,err,error,*999)
    elementDomains=>decomposition%elementDomains
    ALLOCATE(elementDomains%offsets(totalNumberOfLocalElements+1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate element domain offsets.",err,error,*999)
    elementDomains%offsets(1)=1
    DO elementIdx=1,totalNumberOfLocalElements
      elementDomains%offsets(elementIdx+1)=elementDomains%offsets(elementIdx)+SIZE(elementDomainsArrays(elementIdx)%array)
    END DO !elementIdx
    ALLOCATE(elementDomains%domains(elementDomains%offsets(totalNumberOfLocalElements+1)-1),stat=err)
    DO elementIdx=1,totalNumberOfLocalElements
      elementDomains%domains(elementDomains%offsets(elementIdx):elementDomains%offsets(elementIdx+1)-1)= &
        & elementDomainsArrays(elementIdx)%array
      DEALLOCATE(elementDomainsArrays(elementIdx)%array)
    END DO !elementIdx
    IF(err/=0) CALL FlagError("Could not allocate element domains.",err,error,*999)
    DEALLOCATE(elementDomainsArrays)

    !Determine the adjacent domains for the elements in the new distributed mesh.
    NULLIFY(adjacentDomainsList)
    CALL List_CreateStart(adjacentDomainsList,err,error,*999)
    CALL List_DataTypeSet(adjacentDomainsList,LIST_INTG_TYPE,err,error,*999)
    CALL List_InitialSizeSet(adjacentDomainsList,numberOfDomains,err,error,*999)
    CALL List_CreateFinish(adjacentDomainsList,err,error,*999)
    !Add all adjacent domains of the boundary elements. 
    DO elementIdx=1,numberOfLocalElements
      elementOwner=elementDomains%domains(elementDomains%offsets(elementIdx))
      DO elementDomainIdx=elementDomains%offsets(elementIdx)+1,elementDomains%offsets(elementIdx+1)-1
        domainNumber=elementDomains%domains(elementDomainIdx)  
        CALL List_ItemAdd(adjacentDomainsList,domainNumber,err,error,*999)
      END DO !elementDomainIdx
    END DO !domainIdx
    CALL List_RemoveDuplicates(adjacentDomainsList,err,error,*999)
    CALL List_DetachAndDestroy(adjacentDomainsList,numberOfAdjacentDomains,adjacentDomainNumbers,err,error,*999)

    !Calculate the element ghost send and receive indices.
    ALLOCATE(elementGhostSendLists(numberOfAdjacentDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost send lists.",err,error,*999)
    ALLOCATE(elementGhostReceiveLists(numberOfAdjacentDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost send lists.",err,error,*999)

    DO adjacentDomainIdx=1,numberOfAdjacentDomains
      NULLIFY(elementGhostSendLists(adjacentDomainIdx)%ptr)
      CALL List_CreateStart(elementGhostSendLists(adjacentDomainIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(elementGhostSendLists(adjacentDomainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(elementGhostSendLists(adjacentDomainIdx)%ptr,numberOfGhostElements,err,error,*999)
      CALL List_CreateFinish(elementGhostSendLists(adjacentDomainIdx)%ptr,err,error,*999)
      NULLIFY(elementGhostReceiveLists(adjacentDomainIdx)%ptr)
      CALL List_CreateStart(elementGhostReceiveLists(adjacentDomainIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(elementGhostReceiveLists(adjacentDomainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(elementGhostReceiveLists(adjacentDomainIdx)%ptr,numberOfGhostElements,err,error,*999)
      CALL List_CreateFinish(elementGhostReceiveLists(adjacentDomainIdx)%ptr,err,error,*999)
    END DO !adjacentDomainIdx
   
    !Loop over the local elements
    DO elementIdx=1,numberOfLocalElements
      DO elementDomainIdx=elementDomains%offsets(elementIdx)+1,elementDomains%offsets(elementIdx+1)-1
        domainNumber=elementDomains%domains(elementDomainIdx)
        CALL List_Search(adjacentDomainNumbers,domainNumber,adjacentDomainIdx,err,error,*999)
        IF(adjacentDomainIdx<=0) CALL FlagError("Adjacent domain index not found.",err,error,*999)
        CALL List_ItemAdd(elementGhostSendLists(adjacentDomainIdx)%ptr,elementIdx,err,error,*999)
      END DO !elementDomainIdx
    END DO !elementIdx
    !Loop over the ghost elements
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      elementOwner=elementDomains%domains(elementDomains%offsets(elementIdx))
      CALL List_Search(adjacentDomainNumbers,elementOwner,adjacentDomainIdx,err,error,*999)
      IF(adjacentDomainIdx<=0) CALL FlagError("Adjacent domain index not found.",err,error,*999)
      CALL List_ItemAdd(elementGhostReceiveLists(adjacentDomainIdx)%ptr,elementIdx,err,error,*999)
    END DO !elementIdx

    DO adjacentDomainIdx=1,numberOfAdjacentDomains
      CALL List_RemoveDuplicates(elementGhostSendLists(adjacentDomainIdx)%ptr,err,error,*999)
      CALL List_RemoveDuplicates(elementGhostReceiveLists(adjacentDomainIdx)%ptr,err,error,*999)
    END DO !adjacentDomainIdx
   
    !Set up the element domain mappings 
    DO meshComponentIdx=1,mesh%NUMBER_OF_COMPONENTS          
      domain=>decomposition%domain(meshComponentIdx)%ptr
      meshElements=>mesh%topology(meshComponentIdx)%ptr%elements
      !Calculate elements mapping
      elementsMapping=>domain%mappings%elements
      ALLOCATE(elementsMapping%ADJACENT_DOMAINS(numberOfAdjacentDomains),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate adjacent domains.",err,error,*999)
      DO adjacentDomainIdx=1,numberOfAdjacentDomains
        CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx),err,error,*999)
        elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER=adjacentDomainNumbers(adjacentDomainIdx)
        CALL List_Detach(elementGhostSendLists(adjacentDomainIdx)%ptr, &
          & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS, &
          & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES,err,error,*999)
        CALL List_Detach(elementGhostReceiveLists(adjacentDomainIdx)%ptr, &
          & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS, &
          & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES,err,error,*999)
      END DO !adjacentDomainIdx
      elementsMapping%NUMBER_OF_ADJACENT_DOMAINS=numberOfAdjacentDomains
      elementsMapping%NUMBER_OF_LOCAL=numberOfLocalElements
      elementsMapping%TOTAL_NUMBER_OF_LOCAL=totalNumberOfLocalElements
      elementsMapping%NUMBER_OF_GLOBAL=mesh%NUMBER_OF_ELEMENTS

      !When mesh reading/creation is fully parallel, use the DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE routine
      ALLOCATE(elementsMapping%LOCAL_TO_GLOBAL_MAP(totalNumberOfLocalElements),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate element local to global map.",err,error,*999)
      elementsMapping%LOCAL_TO_GLOBAL_MAP=globalElementNumbers

      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Domain mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains  = ",elementsMapping%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global = ",elementsMapping%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of local = ",elementsMapping%NUMBER_OF_LOCAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of local = ",elementsMapping%TOTAL_NUMBER_OF_LOCAL, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Domain numbers:",ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Domain list:",ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map:",ERR,ERROR,*999)
        DO elementIdx=1,elementsMapping%TOTAL_NUMBER_OF_LOCAL
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local index : ",elementIdx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global index = ",elementsMapping% &
            & LOCAL_TO_GLOBAL_MAP(elementIdx),ERR,ERROR,*999)
        ENDDO !elementIdx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
          & elementsMapping%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
        IF(elementsMapping%NUMBER_OF_ADJACENT_DOMAINS>0) THEN
          DO adjacentDomainIdx=1,elementsMapping%NUMBER_OF_ADJACENT_DOMAINS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",adjacentDomainIdx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
              & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
              & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)% &
              & NUMBER_OF_SEND_GHOSTS,8,8,elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES, &
              & '("      Local send ghost indices       :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)      
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of receive ghosts = ", &
              & elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)% &
              & NUMBER_OF_RECEIVE_GHOSTS,8,8,elementsMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES, &
              & '("      Local receive ghost indices    :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)              
          ENDDO !adjacentDomainIdx
        ENDIF
      ENDIF
    END DO !meshComponentIdx
    IF(ALLOCATED(adjacentDomainNumbers)) DEALLOCATE(adjacentDomainNumbers)
    IF(ALLOCATED(elementGhostSendLists)) THEN
      DO adjacentDomainIdx=1,numberOfAdjacentDomains
        CALL List_Destroy(elementGhostSendLists(adjacentDomainIdx)%ptr,err,error,*999)
      END DO !adjacentDomainIdx
      DEALLOCATE(elementGhostSendLists)
    END IF
    IF(ALLOCATED(elementGhostReceiveLists)) THEN
      DO adjacentDomainIdx=1,numberOfAdjacentDomains
        CALL List_Destroy(elementGhostReceiveLists(adjacentDomainIdx)%ptr,err,error,*999)
      END DO !adjacentDomainIdx
      DEALLOCATE(elementGhostReceiveLists)
    END IF

    DO meshComponentIdx=1,mesh%NUMBER_OF_COMPONENTS          
      !The tree/mapping for mapping the user nodes to the local nodes 
      NULLIFY(nodesTree) 
      CALL Tree_CreateStart(nodesTree,err,error,*999)
      CALL Tree_InsertTypeSet(nodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
      CALL Tree_CreateFinish(nodesTree,err,error,*999)

      !Calculate the element node offsets
      ALLOCATE(elementNodeOffsets(totalNumberOfLocalElements+1),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate element node offsets.",err,error,*999)
      elementNodeOffsets(1)=1
      DO elementIdx=1,totalNumberOfLocalElements
        globalElementNumber=globalElementNumbers(elementIdx)
        basis=>meshElements%elements(globalElementNumber)%basis
        elementNodeOffsets(elementIdx+1)=elementNodeOffsets(elementIdx)+basis%NUMBER_OF_NODES
      END DO !elementIdx
      
      !Calculate the local element nodes
      localNodeCount=1
      DO elementIdx=1,totalNumberOfLocalElements
        globalElementNumber=globalElementNumbers(elementIdx)
        basis=>meshElements%elements(globalElementNumber)%basis
        DO elementNodeIdx=1,basis%NUMBER_OF_NODES
          userNodeNumber=meshElements%elements(globalElementNumber)%USER_ELEMENT_NODES(elementNodeIdx)
          CALL Tree_ItemInsert(nodesTree,userNodeNumber,localNodeCount,insertStatus,err,error,*999)
          IF(insertStatus==TREE_NODE_INSERT_SUCESSFUL) localNodeCount=localNodeCount+1
        END DO !elementNodeIdx
      END DO !elementIdx
      ALLOCATE(elementNodes(elementNodeOffsets(totalNumberOfLocalElements+1)-1),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
      !Assign a local number to the element nodes.
      DO elementIdx=1,totalNumberOfLocalElements
        globalElementNumber=globalElementNumbers(elementIdx)
        basis=>meshElements%elements(globalElementNumber)%basis
        DO elementNodeIdx=1,basis%NUMBER_OF_NODES
          userNodeNumber=meshElements%elements(globalElementNumber)%USER_ELEMENT_NODES(elementNodeIdx)
          NULLIFY(treeNode)
          CALL Tree_Search(nodesTree,userNodeNumber,treeNode,err,error,*999)
          IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node was not found.",err,error,*999)
          CALL Tree_NodeValueGet(nodesTree,treeNode,localNodeNumber,err,error,*999) 
          elementNodes(elementNodeOffsets(elementIdx)+elementNodeIdx-1)=localNodeNumber
        END DO !elementNodeIdx
      END DO !elementIdx
      
      !Calculate the nodes mapping
      nodesMapping=>domain%mappings%nodes
      CALL DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE(nodesMapping,elementNodes,elementNodeOffsets, &
        & elementDomains%domains,elementDomains%offsets,err,error,*999)

      ALLOCATE(userNodeNumbers(nodesMapping%TOTAL_NUMBER_OF_LOCAL),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate element nodes.",err,error,*999)
      !The element nodes may have been changed, so the nodes tree needs to be adjusted as well.
      DO elementIdx=1,totalNumberOfLocalElements
        globalElementNumber=globalElementNumbers(elementIdx)
        basis=>meshElements%elements(globalElementNumber)%basis
        DO elementNodeIdx=1,basis%NUMBER_OF_NODES
          userNodeNumber=meshElements%elements(globalElementNumber)%USER_ELEMENT_NODES(elementNodeIdx)
          NULLIFY(treeNode)
          CALL Tree_Search(nodesTree,userNodeNumber,treeNode,err,error,*999)
          IF(.NOT.ASSOCIATED(treeNode)) CALL FlagError("Tree node was not found.",err,error,*999)
          localNodeNumber=elementNodes(elementNodeOffsets(elementIdx)+elementNodeIdx-1)
          CALL Tree_NodeValueSet(nodesTree,treeNode,localNodeNumber,err,error,*999)
          userNodeNumbers(localNodeNumber)=userNodeNumber
        END DO !elementNodeIdx
      END DO !elementIdx

      !Calculate the domain nodes
      nodes=>mesh%region%nodes
      domainNodes=>domain%topology%nodes
      ALLOCATE(domainNodes%NODES(nodesMapping%TOTAL_NUMBER_OF_LOCAL),stat=err)
      IF(ERR/=0) CALL FlagError("Could not allocate domain nodes nodes.",err,error,*999)
      domainNodes%NUMBER_OF_NODES=nodesMapping%NUMBER_OF_LOCAL
      domainNodes%TOTAL_NUMBER_OF_NODES=nodesMapping%TOTAL_NUMBER_OF_LOCAL
      domainNodes%NUMBER_OF_GLOBAL_NODES=nodesMapping%NUMBER_OF_GLOBAL
      domainNodes%NODES_TREE=>nodesTree

      DO nodeIdx=1,domainNodes%TOTAL_NUMBER_OF_NODES
        CALL DOMAIN_TOPOLOGY_NODE_INITIALISE(domainNodes%nodes(nodeIdx),err,error,*999)
        domainNodes%nodes(nodeIdx)%LOCAL_NUMBER=nodeIdx
        domainNodes%nodes(nodeIdx)%MESH_NUMBER=nodesMapping%LOCAL_TO_GLOBAL_MAP(nodeIdx)
        userNodeNumber=userNodeNumbers(nodeIdx)
        domainNodes%nodes(nodeIdx)%USER_NUMBER=userNodeNumber
        NULLIFY(treeNode)
        CALL Tree_Search(nodes%NODES_TREE,userNodeNumber,treeNode,err,error,*999)
        CALL Tree_NodeValueGet(nodes%NODES_TREE,treeNode,globalNodeNumber,err,error,*999)
        domainNodes%nodes(nodeIdx)%GLOBAL_NUMBER=globalNodeNumber
      END DO !nodeIdx
      IF(ALLOCATED(userNodeNumbers)) DEALLOCATE(userNodeNumbers)

      !Calculate the domain elements.
      domainElements=>domain%topology%elements          
      ALLOCATE(domainElements%elements(totalNumberOfLocalElements),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate domain elements elements.",err,error,*999)
      domainElements%NUMBER_OF_ELEMENTS=elementsMapping%NUMBER_OF_LOCAL
      domainElements%TOTAL_NUMBER_OF_ELEMENTS=elementsMapping%TOTAL_NUMBER_OF_LOCAL
      domainElements%NUMBER_OF_GLOBAL_ELEMENTS=elementsMapping%NUMBER_OF_GLOBAL

      DO elementIdx=1,domainElements%TOTAL_NUMBER_OF_ELEMENTS
        CALL DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(domainElements%elements(elementIdx),err,error,*999)
        globalElementNumber=elementsMapping%LOCAL_TO_GLOBAL_MAP(elementIdx)
        basis=>meshElements%elements(globalElementNumber)%basis
        domainElements%elements(elementIdx)%NUMBER=elementIdx
        domainElements%elements(elementIdx)%basis=>basis
        ALLOCATE(domainElements%elements(elementIdx)%ELEMENT_NODES(basis%NUMBER_OF_NODES),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate domain elements element nodes.",err,error,*999)
        domainElements%elements(elementIdx)%ELEMENT_NODES= &
          & elementNodes(elementNodeOffsets(elementIdx):elementNodeOffsets(elementIdx+1)-1)
        ALLOCATE(domainElements%elements(elementIdx)%ELEMENT_DERIVATIVES(basis%MAXIMUM_NUMBER_OF_DERIVATIVES, &
          & basis%NUMBER_OF_NODES),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate domain elements element derivatives.",err,error,*999)
        ALLOCATE(domainElements%elements(elementIdx)%elementVersions(basis%MAXIMUM_NUMBER_OF_DERIVATIVES, &
          & basis%NUMBER_OF_NODES),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate domain elements element versions.",err,error,*999)
        domainElements%elements(elementIdx)%elementVersions= &
          & meshElements%elements(globalElementNumber)%USER_ELEMENT_NODE_VERSIONS
      END DO !elementIdx

      IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
      IF(ALLOCATED(elementNodeOffsets)) DEALLOCATE(elementNodeOffsets)

      CALL DomainTopology_Calculate(domain%topology,err,error,*999)

      !Calculate the element derivatives
      DO elementIdx=1,domainElements%TOTAL_NUMBER_OF_ELEMENTS
        basis=>domainElements%elements(elementIdx)%basis
        DO elementNodeIdx=1,basis%NUMBER_OF_NODES
          localNodeNumber=domainElements%elements(elementIdx)%ELEMENT_NODES(elementNodeIdx)
          DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(elementNodeIdx)
            !Find equivalent node derivative by matching partial derivative index
            !/todo Take a look at this - is it needed any more?
            found=.FALSE.
            DO nodeDerivativeIdx=1,domainNodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES
              IF(domainNodes%nodes(localNodeNumber)%DERIVATIVES(nodeDerivativeIdx)%PARTIAL_DERIVATIVE_INDEX== &
                & basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,elementNodeIdx)) THEN
                found=.TRUE.
                EXIT
              ENDIF
            ENDDO !nodeDerivativeIdx
            IF(found) THEN
              domainElements%elements(elementIdx)%ELEMENT_DERIVATIVES(derivativeIdx,elementNodeIdx)=nodeDerivativeIdx
            ELSE
              CALL FlagError("Could not find equivalent node derivative",err,error,*999)
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !elementNodeIdx
      END DO !elementIdx

      IF(ALLOCATED(elementNodes)) DEALLOCATE(elementNodes)
      IF(ALLOCATED(elementNodeOffsets)) DEALLOCATE(elementNodeOffsets)

      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Domain topology elements:",err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain elements = ", &
          & domainElements%TOTAL_NUMBER_OF_ELEMENTS,err,error,*999)
        DO elementIdx=1,domainElements%TOTAL_NUMBER_OF_ELEMENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",domainElements%ELEMENTS(elementIdx)%NUMBER, &
            & err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis user number = ", &
            & domainElements%ELEMENTS(elementIdx)%BASIS%USER_NUMBER,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local nodes = ", &
            & domainElements%ELEMENTS(elementIdx)%BASIS%NUMBER_OF_NODES,err,error,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainElements%ELEMENTS(elementIdx)%BASIS%NUMBER_OF_NODES,8,8, &
            & domainElements%ELEMENTS(elementIdx)%ELEMENT_NODES,'("    Element nodes(nn) :",8(X,I9))','(23X,8(X,I9))', &
            & ERR,error,*999)
          DO elementNodeIdx=1,domainElements%ELEMENTS(elementIdx)%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Element node number : ",elementNodeIdx,err,error,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainElements%ELEMENTS(elementIdx)%BASIS% &
              & NUMBER_OF_DERIVATIVES(elementNodeIdx),8,8, &
              & domainElements%ELEMENTS(elementIdx)%ELEMENT_DERIVATIVES(:,elementNodeIdx), &
              & '("        Element derivatives :",8(X,I2))','(29X,8(X,I2))',err,error,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainElements%ELEMENTS(elementIdx)%BASIS% &
              & NUMBER_OF_DERIVATIVES(elementNodeIdx),8,8, &
              & domainElements%ELEMENTS(elementIdx)%elementVersions(:,elementNodeIdx), &
              & '("        Element versions    :",8(X,I2))','(29X,8(X,I2))',err,error,*999)
          ENDDO !elementNodeIdx
        ENDDO !elementIdx
      ENDIF
    END DO !meshComponentIdx
    IF(ALLOCATED(adjacentDomainNumbers)) DEALLOCATE(adjacentDomainNumbers)
    IF(ALLOCATED(elementGhostSendLists)) THEN
      DO adjacentDomainIdx=1,numberOfAdjacentDomains
        CALL List_Destroy(elementGhostSendLists(adjacentDomainIdx)%ptr,err,error,*999)
      END DO !adjacentDomainIdx
      DEALLOCATE(elementGhostSendLists)
    END IF
    IF(ALLOCATED(elementGhostReceiveLists)) THEN
      DO adjacentDomainIdx=1,numberOfAdjacentDomains
        CALL List_Destroy(elementGhostReceiveLists(adjacentDomainIdx)%ptr,err,error,*999)
      END DO !adjacentDomainIdx
      DEALLOCATE(elementGhostReceiveLists)
    END IF
    IF(ALLOCATED(globalElementNumbers)) DEALLOCATE(globalElementNumbers)

    EXITS("Domain_Calculate")
    RETURN
999 IF(ALLOCATED(receiveBuffers)) THEN
      DO domainIdx=1,SIZE(receiveBuffers)
        IF(ALLOCATED(receiveBuffers(domainIdx)%array)) DEALLOCATE(receiveBuffers(domainIdx)%array)
      END DO
      DEALLOCATE(receiveBuffers)
    END IF
    IF(ALLOCATED(sendBuffers)) THEN
      DO domainIdx=1,SIZE(sendBuffers)
        IF(ALLOCATED(sendBuffers(domainIdx)%array)) DEALLOCATE(sendBuffers(domainIdx)%array)
      END DO
      DEALLOCATE(sendBuffers)
    END IF
    IF(ALLOCATED(elementDomainsArrays)) THEN
      DO elementIdx=1,SIZE(elementDomainsArrays)
        IF(ALLOCATED(elementDomainsArrays(elementIdx)%array)) DEALLOCATE(elementDomainsArrays(elementIdx)%array)
      END DO
      DEALLOCATE(elementDomainsArrays)
    END IF
    IF(ALLOCATED(sendBufferLists)) THEN
      DO domainIdx=1,SIZE(sendBufferLists)
        CALL List_Destroy(sendBufferLists(domainIdx)%ptr,dummyErr,dummyError,*998)
      END DO
      DEALLOCATE(sendBufferLists)
    END IF
998 IF(ALLOCATED(elementGhostSendLists)) THEN
      DO adjacentDomainIdx=1,SIZE(elementGhostSendLists)
        CALL List_Destroy(elementGhostSendLists(adjacentDomainIdx)%ptr,dummyErr,dummyError,*997)
      END DO
      DEALLOCATE(elementGhostSendLists)
    END IF
997 IF(ALLOCATED(elementGhostReceiveLists)) THEN
      DO adjacentDomainIdx=1,SIZE(elementGhostReceiveLists)
        CALL List_Destroy(elementGhostReceiveLists(adjacentDomainIdx)%ptr,dummyErr,dummyError,*996)
      END DO
      DEALLOCATE(elementGhostReceiveLists)
    END IF
996 CALL List_Destroy(adjacentDomainsList,dummyErr,dummyError,*995)
995 CALL List_Destroy(userNodeNumbersList,dummyErr,dummyError,*994)
994 CALL List_Destroy(ghostUserNodeNumbersList,dummyErr,dummyError,*993)
993 CALL List_Destroy(elementNodesList,dummyErr,dummyError,*992)
992 ERRORSEXITS("Domain_Calculate",err,error)
    RETURN 1   
  END SUBROUTINE Domain_Calculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the domain topology.
  SUBROUTINE DomainTopology_Calculate(topology,err,error,*)    

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology !<A pointer to the domain topology to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainTopology_Calculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      !Calculate the elements surrounding the nodes in a domain
      CALL DomainTopology_SurroundingElementsCalculate(topology,err,error,*999)
      !Calculate the derivatives and versions at each node in a domain
      CALL DomainTopology_NodesDerivativesVersionsCalculate(topology,err,error,*999)
    ELSE
      CALL FlagError("Topology is not associated",err,error,*999)
    ENDIF
    
    EXITS("DomainTopology_Calculate")
    RETURN
999 ERRORSEXITS("DomainTopology_Calculate",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_Calculate
  
  !
  !===============================================================================================================================
  !

  !>Calculates the the derivatives and versions at each node in a topology.
  SUBROUTINE DomainTopology_NodesDerivativesVersionsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology !<A pointer to the domain topology to calculate the derivates at each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,element,elementIdx,globalDerivative,elementNodeIdx,maxNumberOfDerivatives,nodeIdx, &
      & numberOfDerivatives,numberOfVersions,versionIdx,dofIdx,adjacentDomainIdx,ghostSendIdx,localNodeNumber, &
      & partialDerivativeIndex,versionNumber,adjacentDomainNumber,receiveBufferIdx,ghostReceiveIdx,MPI_IERROR
    INTEGER(INTG), ALLOCATABLE :: derivatives(:),versions(:),sendCounts(:),receiveCounts(:),requests(:),statuses(:,:)
    TYPE(INTEGER_INTG_ALLOC_TYPE), ALLOCATABLE :: sendBuffers(:),receiveBuffers(:)
    LOGICAL :: found
    TYPE(LIST_TYPE), POINTER :: nodeDerivativeList,ghostSendList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: nodeVersionList(:)
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: elements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: nodes
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: nodesMapping
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DomainTopology_NodesDerivativesVersionsCalculate",err,ERROR,*999)

    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes
        IF(ASSOCIATED(nodes)) THEN
          !Loop over the local domain nodes (not the ghosts!)
          DO nodeIdx=1,nodes%NUMBER_OF_NODES
            NULLIFY(nodeDerivativeList)
            CALL List_CreateStart(nodeDerivativeList,err,error,*999)
            CALL List_DataTypeSet(nodeDerivativeList,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(nodeDerivativeList,8,err,error,*999)
            CALL List_CreateFinish(nodeDerivativeList,err,error,*999)
            maxNumberOfDerivatives=-1
            DO elementIdx=1,nodes%nodes(nodeIdx)%NUMBER_OF_SURROUNDING_ELEMENTS
              element=nodes%nodes(nodeIdx)%SURROUNDING_ELEMENTS(elementIdx)
              basis=>elements%elements(element)%basis
              !Find the local node corresponding to this node
              found=.FALSE.
              DO elementNodeIdx=1,basis%NUMBER_OF_NODES
                IF(elements%elements(element)%ELEMENT_NODES(elementNodeIdx)==nodeIdx) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !nn
              IF(found) THEN
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(elementNodeIdx)
                  CALL List_ItemAdd(nodeDerivativeList,basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,elementNodeIdx),err,error,*999)
                ENDDO !derivativeIdx
                maxNumberOfDerivatives=MAX(basis%NUMBER_OF_DERIVATIVES(elementNodeIdx),maxNumberOfDerivatives)
              ELSE
                CALL FlagError("Could not find local node.",err,error,*999)
              ENDIF
            ENDDO !elem_idx
            CALL List_RemoveDuplicates(nodeDerivativeList,err,error,*999)
            CALL List_DetachAndDestroy(nodeDerivativeList,numberOfDerivatives,derivatives,err,error,*999)
            IF(numberOfDerivatives==maxNumberOfDerivatives) THEN
              !Set up the node derivatives.
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(maxNumberOfDerivatives),STAT=err)
              nodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES=maxNumberOfDerivatives
              DO derivativeIdx=1,numberOfDerivatives
                CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx),err,error,*999)
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%PARTIAL_DERIVATIVE_INDEX=derivatives(derivativeIdx)
                globalDerivative=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(derivatives(derivativeIdx))
                IF(globalDerivative/=0) THEN
                  nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX=globalDerivative
                ELSE
                  localError="The partial derivative index of "//TRIM(NumberToVstring(derivatives(derivativeIdx),"*", &
                    & err,error))//" for derivative number "//TRIM(NumberToVstring(derivativeIdx,"*",err,error))// &
                    & " does not have a corresponding global derivative."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
              DEALLOCATE(derivatives)

              !Calculate the number of versions at each node.
              ALLOCATE(nodeVersionList(numberOfDerivatives),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate node version list.",err,error,*999)
              DO derivativeIdx=1,nodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                NULLIFY(nodeVersionList(derivativeIdx)%ptr)
                CALL List_CreateStart(nodeVersionList(derivativeIdx)%ptr,err,error,*999)
                CALL List_DataTypeSet(nodeVersionList(derivativeIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
                CALL List_InitialSizeSet(nodeVersionList(derivativeIdx)%ptr,8,err,error,*999)
                CALL List_CreateFinish(nodeVersionList(derivativeIdx)%ptr,err,error,*999)
              ENDDO !derivativeIdx
              DO elementIdx=1,nodes%nodes(nodeIdx)%NUMBER_OF_SURROUNDING_ELEMENTS
                element=nodes%nodes(nodeIdx)%SURROUNDING_ELEMENTS(elementIdx)
                basis=>elements%elements(element)%basis
                !Find the local node corresponding to this node
                found=.FALSE.
                DO elementNodeIdx=1,basis%NUMBER_OF_NODES
                  IF(elements%elements(element)%ELEMENT_NODES(elementNodeIdx)==nodeIdx) THEN
                    found=.TRUE.
                    EXIT
                  ENDIF
                ENDDO !nn
                IF(found) THEN
                  DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(elementNodeIdx)
                    CALL List_ItemAdd(nodeVersionList(derivativeIdx)%ptr,elements%elements(element)%elementVersions( &
                      & derivativeIdx,elementNodeIdx),err,error,*999)
                  ENDDO !derivativeIdx
                ELSE
                  CALL FlagError("Could not find local node.",err,error,*999)
                ENDIF
              ENDDO !elementIdx
              DO derivativeIdx=1,nodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                CALL List_RemoveDuplicates(nodeVersionList(derivativeIdx)%ptr,err,error,*999)
                CALL List_DetachAndDestroy(nodeVersionList(derivativeIdx)%ptr,numberOfVersions,versions,err,error,*999)
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions=numberOfVersions
                ALLOCATE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(numberOfVersions),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate node global derivative index.",err,error,*999)
                DO versionIdx=1,numberOfVersions 
                  nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx)=versions(versionIdx)
                ENDDO !versionIdx
                IF(ALLOCATED(versions)) DEALLOCATE(versions)
              ENDDO!derivativeIdx
              IF(ALLOCATED(nodeVersionList)) DEALLOCATE(nodeVersionList)
            ELSE
              localError="Invalid domain configuration. User node "// &
                & TRIM(NumberToVstring(nodes%nodes(nodeIdx)%USER_NUMBER,"*",err,error))// &
                & " has inconsistent derivative directions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !nodeIdx

          !Send the partial derivative indices and number of versions the to ghost nodes, as some ghost nodes may not have all the information locally.          
          nodesMapping=>topology%domain%mappings%nodes
          IF(.NOT.ASSOCIATED(nodesMapping)) CALL FlagError("Nodes mapping is not associated.",err,error,*999)
          ALLOCATE(sendCounts(nodesMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate send counts.",err,error,*999)
          ALLOCATE(sendBuffers(nodesMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate send buffers.",err,error,*999)

          !Pack the send buffers
          DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
            NULLIFY(ghostSendList)
            CALL List_CreateStart(ghostSendList,err,error,*999)
            CALL List_DataTypeSet(ghostSendList,LIST_INTG_TYPE,err,error,*999)
            CALL List_InitialSizeSet(ghostSendList,8,err,error,*999)
            CALL List_CreateFinish(ghostSendList,err,error,*999)
            DO ghostSendIdx=1,nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS
              localNodeNumber=nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES(ghostSendIdx)              
              !Add the number of derivatives to the send buffer
              numberOfDerivatives=nodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES
              CALL List_ItemAdd(ghostSendList,numberOfDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfDerivatives
                !Add the partial derivative indices to the send buffer
                partialDerivativeIndex=nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%PARTIAL_DERIVATIVE_INDEX
                CALL List_ItemAdd(ghostSendList,partialDerivativeIndex,err,error,*999)
                !Add the number of versions to the send buffer
                numberOfVersions=nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%numberOfVersions
                CALL List_ItemAdd(ghostSendList,partialDerivativeIndex,err,error,*999)
                DO versionIdx=1,numberOfVersions
                  !Add the version number to the send buffer
                  versionNumber=nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx)
                  CALL List_ItemAdd(ghostSendList,versionNumber,err,error,*999)
                END DO !versionIdx
              END DO !derivativeIdx
            END DO !ghostSendIdx
            CALL List_DetachAndDestroy(ghostSendList,sendCounts(adjacentDomainIdx),sendBuffers(adjacentDomainIdx)%array, &
              & err,error,*999)
          END DO !adjacentDomainIdx

          ALLOCATE(receiveCounts(nodesMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate receive counts.",err,error,*999)
          ALLOCATE(receiveBuffers(nodesMapping%NUMBER_OF_ADJACENT_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate receive buffers.",err,error,*999)
          
          ALLOCATE(statuses(MPI_STATUS_SIZE,2*nodesMapping%NUMBER_OF_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate mpi statuses.",err,error,*999)
          ALLOCATE(requests(2*nodesMapping%NUMBER_OF_DOMAINS),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate mpi receive requests.",err,error,*999)
          requests=MPI_REQUEST_NULL

          !Get the receive counts from the adjacent domains.
          DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
            adjacentDomainNumber=nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
            CALL MPI_IRECV(receiveCounts(adjacentDomainIdx),1,MPI_INTEGER,adjacentDomainNumber-1, &
              & MPI_NODE_DERIVATIVE_COUNT_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,requests(adjacentDomainNumber),MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,err,error,*999)
          END DO
          DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
            adjacentDomainNumber=nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
            CALL MPI_ISEND(sendCounts(adjacentDomainIdx),1,MPI_INTEGER,adjacentDomainNumber-1, &
              & MPI_NODE_DERIVATIVE_COUNT_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
              & requests(adjacentDomainNumber+nodesMapping%NUMBER_OF_DOMAINS),MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,err,error,*999)
          END DO

          CALL MPI_WAITALL(2*nodesMapping%NUMBER_OF_DOMAINS,requests,statuses,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,err,error,*999)
          
          !Allocate all the receive buffers
          DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
            ALLOCATE(receiveBuffers(adjacentDomainIdx)%array(receiveCounts(adjacentDomainIdx)),stat=err)
            IF(err/=0) CALL FlagError("Could not allocate receive buffer.",err,error,*999)
          END DO !adjacentDomainIdx

          requests=MPI_REQUEST_NULL
          !Now get all the receive buffers from the adjacent domains.
          DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
            adjacentDomainNumber=nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
            CALL MPI_IRECV(receiveBuffers(adjacentDomainIdx)%array,receiveCounts(adjacentDomainIdx),MPI_INTEGER, &
              & adjacentDomainNumber-1,MPI_NODE_DERIVATIVE_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
              & requests(adjacentDomainNumber),MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,err,error,*999)
          END DO
          DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
            adjacentDomainNumber=nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER
            CALL MPI_ISEND(sendBuffers(adjacentDomainIdx)%array,sendCounts(adjacentDomainIdx),MPI_INTEGER, &
              & adjacentDomainNumber-1,MPI_NODE_DERIVATIVE_EXCHANGE_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM, &
              & requests(nodesMapping%NUMBER_OF_DOMAINS+adjacentDomainNumber),MPI_IERROR)
            CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,err,error,*999)
          END DO

          CALL MPI_WAITALL(2*nodesMapping%NUMBER_OF_DOMAINS,requests,statuses,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,err,error,*999)

          IF(ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
          IF(ALLOCATED(receiveCounts)) DEALLOCATE(receiveCounts)
          IF(ALLOCATED(statuses)) DEALLOCATE(statuses)
          IF(ALLOCATED(requests)) DEALLOCATE(requests)

          !Unpack the receive buffers
          DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
            receiveBufferIdx=0
            DO ghostReceiveIdx=1,nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS
              localNodeNumber=nodesMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES(ghostReceiveIdx)              
              !Retrieve the number of derivatives of the received node.
              receiveBufferIdx=receiveBufferIdx+1
              numberOfDerivatives=receiveBuffers(adjacentDomainIdx)%array(receiveBufferIdx)
              ALLOCATE(nodes%nodes(localNodeNumber)%derivatives(numberOfDerivatives),STAT=err)
              nodes%nodes(localNodeNumber)%NUMBER_OF_DERIVATIVES=numberOfDerivatives
              DO derivativeIdx=1,numberOfDerivatives
                CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(nodes%nodes(localNodeNumber)%derivatives(derivativeIdx),err,error,*999)
                !Retreive the partial derivative index of the derivativeIdx'th partial derivative of the received node.
                receiveBufferIdx=receiveBufferIdx+1
                partialDerivativeIndex=receiveBuffers(adjacentDomainIdx)%array(receiveBufferIdx)
                nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%PARTIAL_DERIVATIVE_INDEX=partialDerivativeIndex
                globalDerivative=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(partialDerivativeIndex)
                IF(globalDerivative/=0) THEN
                  nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX=globalDerivative
                ELSE
                  localError="The partial derivative index of "//TRIM(NumberToVstring(partialDerivativeIndex,"*", &
                    & err,error))//" for derivative number "//TRIM(NumberToVstring(derivativeIdx,"*",err,error))// &
                    & " does not have a corresponding global derivative."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Retreive the number of versions
                receiveBufferIdx=receiveBufferIdx+1
                numberOfVersions=receiveBuffers(adjacentDomainIdx)%array(receiveBufferIdx)
                nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%numberOfVersions=numberOfVersions
                ALLOCATE(nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%userVersionNumbers(numberOfVersions),stat=err)
                IF(err/=0) CALL FlagError("Could not allocate user version numbers.",err,error,*999)
                DO versionIdx=1,nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%numberOfVersions 
                  !Retreive the version numbers 
                  receiveBufferIdx=receiveBufferIdx+1
                  versionNumber=receiveBuffers(adjacentDomainIdx)%array(receiveBufferIdx)
                  nodes%nodes(localNodeNumber)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx)=versionNumber
                ENDDO !versionIdx
              ENDDO !derivativeIdx
            END DO !ghostReceiveIdx
          END DO !adjacentDomainIdx

          IF(ALLOCATED(sendBuffers)) THEN
            DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
              IF(ALLOCATED(sendBuffers(adjacentDomainIdx)%array)) DEALLOCATE(sendBuffers(adjacentDomainIdx)%array)
            END DO
            DEALLOCATE(sendBuffers)
          END IF
          IF(ALLOCATED(receiveBuffers)) THEN
            DO adjacentDomainIdx=1,nodesMapping%NUMBER_OF_ADJACENT_DOMAINS
              IF(ALLOCATED(receiveBuffers(adjacentDomainIdx)%array)) DEALLOCATE(receiveBuffers(adjacentDomainIdx)%array)
            END DO
            DEALLOCATE(receiveBuffers)
          END IF

          !Calculate the dof indices
          dofIdx=0
          DO nodeIdx=1,nodes%TOTAL_NUMBER_OF_NODES
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
              numberOfVersions=nodes%nodes(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
              ALLOCATE(nodes%nodes(nodeIdx)%DERIVATIVES(derivativeIdx)%DOF_INDEX(numberOfVersions),stat=err)
              IF(err/=0) CALL FlagError("Could not allocate node dervative versions dof index.",err,error,*999)
              DO versionIdx=1,numberOfVersions
                dofIdx=dofIdx+1
                nodes%nodes(nodeIdx)%DERIVATIVES(derivativeIdx)%DOF_INDEX(versionIdx)=dofIdx
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
          IF(DIAGNOSTICS1) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Domain topology nodes:",err,error,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain nodes = ",nodes%TOTAL_NUMBER_OF_NODES, &
              & ERR,error,*999)
            DO nodeIdx=1,nodes%TOTAL_NUMBER_OF_NODES
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",nodes%NODES(nodeIdx)%LOCAL_NUMBER, &
                & ERR,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",nodes%NODES(nodeIdx)%MESH_NUMBER, &
                & ERR,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",nodes%NODES(nodeIdx)%GLOBAL_NUMBER, &
                & ERR,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",nodes%NODES(nodeIdx)%USER_NUMBER, &
                & ERR,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
                & nodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES,err,error,*999)
              DO derivativeIdx=1,nodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node local derivative number = ",derivativeIdx,err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Global derivative index = ", &
                  & nodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX,err,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Partial derivative index = ", &
                  & nodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%PARTIAL_DERIVATIVE_INDEX,err,error,*999)
                CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
                  & nodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions,4,4, &
                  & nodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%DOF_INDEX, &
                  & '("        Degree-of-freedom index(version_idx)  :",4(X,I9))','(36X,4(X,I9))',err,error,*999)
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDIF
        ELSE
          CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DomainTopology_NodesDerivativesVersionsCalculate")
    RETURN
999 IF(ALLOCATED(derivatives)) DEALLOCATE(derivatives)
    IF(ALLOCATED(versions)) DEALLOCATE(versions)
    IF(ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
    IF(ALLOCATED(receiveCounts)) DEALLOCATE(receiveCounts)
    IF(ALLOCATED(requests)) DEALLOCATE(requests)
    IF(ALLOCATED(statuses)) DEALLOCATE(statuses)
    IF(ALLOCATED(receiveBuffers)) THEN
      DO adjacentDomainIdx=1,SIZE(receiveBuffers)
        IF(ALLOCATED(receiveBuffers(adjacentDomainIdx)%array)) DEALLOCATE(receiveBuffers(adjacentDomainIdx)%array)
      END DO
      DEALLOCATE(receiveBuffers)
    END IF
    IF(ALLOCATED(sendBuffers)) THEN
      DO adjacentDomainIdx=1,SIZE(sendBuffers)
        IF(ALLOCATED(sendBuffers(adjacentDomainIdx)%array)) DEALLOCATE(sendBuffers(adjacentDomainIdx)%array)
      END DO
      DEALLOCATE(sendBuffers)
    END IF
    IF(ALLOCATED(nodeVersionList)) THEN
      DO derivativeIdx=1,SIZE(nodeVersionList)
        IF(ASSOCIATED(nodeVersionList(derivativeIdx)%ptr)) CALL List_Destroy(nodeVersionList(derivativeIdx)%ptr,err,error,*998)
      END DO !derivativeIdx
      DEALLOCATE(nodeVersionList)
    END IF
998 IF(ASSOCIATED(nodeDerivativeList)) CALL List_Destroy(nodeDerivativeList,err,error,*997)
997 ERRORSEXITS("DomainTopology_NodesDerivativesVersionsCalculate",err,error)
    RETURN 1

  END SUBROUTINE DomainTopology_NodesDerivativesVersionsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a domain.
  SUBROUTINE DomainTopology_SurroundingElementsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology !<A pointer to the domain topology to calculate the elements surrounding each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: element,elementIdx,insertPosition,elementNodeIdx,node,nodeIdx,surroundingElementNumber,surroundingElementIdx
    INTEGER(INTG), ALLOCATABLE :: newSurroundingElements(:)
    LOGICAL :: foundElement
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: elements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: nodes

    ENTERS("DomainTopology_SurroundingElementsCalculate",err,error,*999)
    
    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements      
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes       
        IF(ASSOCIATED(nodes)) THEN
          IF(ALLOCATED(nodes%nodes)) THEN
            DO elementIdx=1,elements%TOTAL_NUMBER_OF_ELEMENTS
              basis=>elements%elements(elementIdx)%basis
              DO elementNodeIdx=1,basis%NUMBER_OF_NODES
                node=elements%elements(elementIdx)%ELEMENT_NODES(elementNodeIdx)
                foundElement=.FALSE.
                element=1
                insertPosition=1
                DO WHILE(element<=nodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS.AND..NOT.foundElement)
                  surroundingElementNumber=nodes%nodes(node)%SURROUNDING_ELEMENTS(element)
                  IF(surroundingElementNumber==elementIdx) THEN
                    foundElement=.TRUE.
                  ENDIF
                  element=element+1
                  IF(elementIdx>=surroundingElementNumber) THEN
                    insertPosition=element
                  ENDIF
                ENDDO
                IF(.NOT.foundElement) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(newSurroundingElements(nodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS+1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate new surrounding elements.",err,error,*999)
                  IF(nodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                    newSurroundingElements(1:insertPosition-1)=nodes%nodes(node)%SURROUNDING_ELEMENTS(1:insertPosition-1)
                    newSurroundingElements(insertPosition+1:nodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS+1)= &
                      & nodes%nodes(node)%SURROUNDING_ELEMENTS(insertPosition:nodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS)
                  END IF
                  newSurroundingElements(insertPosition)=elementIdx
                  CALL MOVE_ALLOC(newSurroundingElements,nodes%nodes(node)%SURROUNDING_ELEMENTS)
                  nodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS=nodes%nodes(node)%NUMBER_OF_SURROUNDING_ELEMENTS+1
                ENDIF
              ENDDO !elementNodeIdx
            ENDDO !elementIdx
            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Nodes surrounding elements :",err,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain nodes = ",nodes%TOTAL_NUMBER_OF_NODES, &
                & ERR,error,*999)
              DO nodeIdx=1,nodes%TOTAL_NUMBER_OF_NODES
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",nodes%NODES(nodeIdx)%LOCAL_NUMBER, &
                  & ERR,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",nodes%NODES(nodeIdx)%MESH_NUMBER, &
                  & ERR,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",nodes%NODES(nodeIdx)%GLOBAL_NUMBER, &
                  & ERR,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",nodes%NODES(nodeIdx)%USER_NUMBER, &
                  & ERR,error,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
                  & nodes%NODES(nodeIdx)%NUMBER_OF_SURROUNDING_ELEMENTS,err,error,*999)
                DO surroundingElementIdx=1,nodes%NODES(nodeIdx)%NUMBER_OF_SURROUNDING_ELEMENTS
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Surrounding element idx = ",surroundingElementIdx,err,error,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ", &
                    & nodes%NODES(nodeIdx)%SURROUNDING_ELEMENTS(surroundingElementIdx),err,error,*999)
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDIF
          ELSE
            CALL FlagError("Domain topology nodes nodes have not been allocated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain topology nodes are not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology not associated.",err,error,*999)
    ENDIF

    EXITS("DomainTopology_SurroundingElementsCalculate")
    RETURN
999 IF(ALLOCATED(newSurroundingElements)) DEALLOCATE(newSurroundingElements)
    ERRORSEXITS("DomainTopology_SurroundingElementsCalculate",err,error)
    RETURN 1   
  END SUBROUTINE DomainTopology_SurroundingElementsCalculate
  
  !
  !===============================================================================================================================
  !

  !>Finalises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ELEMENT_NODES)) DEALLOCATE(ELEMENT%ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%ELEMENT_DERIVATIVES)) DEALLOCATE(ELEMENT%ELEMENT_DERIVATIVES)
    IF(ALLOCATED(ELEMENT%elementVersions)) DEALLOCATE(ELEMENT%elementVersions)
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%NUMBER=0
    NULLIFY(ELEMENT%BASIS)
  
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_INITALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_INITALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given domain topology. \todo pass in the domain elements
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DOMAIN_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ALLOCATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Decomposition already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DOMAIN=>TOPOLOGY%DOMAIN
        TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=0
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE
  
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given domain. \todo pass in domain topology
  SUBROUTINE DOMAIN_TOPOLOGY_FINALISE(DOMAIN,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_TOPOLOGY_NODES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_LINES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_FACES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%TOPOLOGY)
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given domain. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !A pointer to the domain to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        CALL FlagError("Domain already has topology associated",ERR,ERROR,*999)
      ELSE
        !Allocate domain topology
        ALLOCATE(DOMAIN%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Domain topology could not be allocated",ERR,ERROR,*999)
        DOMAIN%TOPOLOGY%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%TOPOLOGY%ELEMENTS)
        NULLIFY(DOMAIN%TOPOLOGY%NODES)
        NULLIFY(DOMAIN%TOPOLOGY%LINES)
        NULLIFY(DOMAIN%TOPOLOGY%FACES)
        !Initialise the topology components
        CALL DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_NODES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_LINES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_FACES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    IF(ALLOCATED(LINE%NODES_IN_LINE)) DEALLOCATE(LINE%NODES_IN_LINE)
    IF(ALLOCATED(LINE%DERIVATIVES_IN_LINE)) DEALLOCATE(LINE%DERIVATIVES_IN_LINE)
 
    EXITS("DOMAIN_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a domain topology line.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    LINE%BOUNDARY_LINE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given domain topology. \todo pass in domain lines
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DOMAIN_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !
  !>Finalises a face in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    IF(ALLOCATED(FACE%NODES_IN_FACE)) DEALLOCATE(FACE%NODES_IN_FACE)
    IF(ALLOCATED(FACE%DERIVATIVES_IN_FACE)) DEALLOCATE(FACE%DERIVATIVES_IN_FACE)
 
    EXITS("DOMAIN_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a domain topology face.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    FACE%BOUNDARY_FACE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given domain topology. \todo pass in domain faces
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DOMAIN_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a domain nodes topolgoy. 
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(DOMAIN_TOPOLOGY,USER_NODE_NUMBER,NODE_EXISTS,DOMAIN_LOCAL_NODE_NUMBER, &
    & GHOST_NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY !<A pointer to the domain topology to check the node exists on
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: NODE_EXISTS !<On exit, is .TRUE. if the node user number exists in the domain nodes topolgoy (even if it is a ghost node), .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_LOCAL_NODE_NUMBER !<On exit, if the node exists the local number corresponding to the user node number. If the node does not exist then global number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_NODE !<On exit, is .TRUE. if the local node (if it exists) is a ghost node, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    ENTERS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    DOMAIN_LOCAL_NODE_NUMBER=0
    GHOST_NODE=.FALSE.
    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
      IF(ASSOCIATED(DOMAIN_NODES)) THEN
        NULLIFY(TREE_NODE)
        CALL Tree_Search(DOMAIN_NODES%NODES_TREE,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL Tree_NodeValueGet(DOMAIN_NODES%NODES_TREE,TREE_NODE,DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
          NODE_EXISTS=.TRUE.
          GHOST_NODE=DOMAIN_LOCAL_NODE_NUMBER>DOMAIN_NODES%NUMBER_OF_NODES
        ENDIF
      ELSE
        CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node derivative and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE_DERIVATIVE%userVersionNumbers)) DEALLOCATE(NODE_DERIVATIVE%userVersionNumbers)
    IF(ALLOCATED(NODE_DERIVATIVE%DOF_INDEX)) DEALLOCATE(NODE_DERIVATIVE%DOF_INDEX)

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR,*999)

    NODE_DERIVATIVE%numberOfVersions=0
    NODE_DERIVATIVE%GLOBAL_DERIVATIVE_INDEX=0
    NODE_DERIVATIVE%PARTIAL_DERIVATIVE_INDEX=0

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx

    ENTERS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE%DERIVATIVES)) THEN
      DO derivative_idx=1,NODE%NUMBER_OF_DERIVATIVES
        CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
      ENDDO !derivative_idx
      DEALLOCATE(NODE%DERIVATIVES)
    ENDIF
    IF(ALLOCATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(NODE%NODE_LINES)) DEALLOCATE(NODE%NODE_LINES)
 
    EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%LOCAL_NUMBER=0
    NODE%MESH_NUMBER=0
    NODE%GLOBAL_NUMBER=0
    NODE%USER_NUMBER=0
    NODE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    NODE%NUMBER_OF_NODE_LINES=0
    NODE%BOUNDARY_NODE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the nodees in the given domain topology. \todo pass in domain nodes
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np

    ENTERS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
          CALL DOMAIN_TOPOLOGY_NODE_FINALISE(TOPOLOGY%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        IF(ALLOCATED(TOPOLOGY%NODES%NODES)) DEALLOCATE(TOPOLOGY%NODES%NODES)
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES_TREE)) CALL Tree_Destroy(TOPOLOGY%NODES%NODES_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%NODES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        CALL FlagError("Decomposition already has topology nodes associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology nodes",ERR,ERROR,*999)
        TOPOLOGY%NODES%NUMBER_OF_NODES=0
        TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES=0
        TOPOLOGY%NODES%NUMBER_OF_GLOBAL_NODES=0
        TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=0
        TOPOLOGY%NODES%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%NODES%NODES_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a mesh. \see OPENCMISS::CMISSMeshCreateFinish
  SUBROUTINE MESH_CREATE_FINISH(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    LOGICAL :: FINISHED
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
        !Check that the mesh component elements have been finished
        FINISHED=.TRUE.
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
            IF(.NOT.MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS_FINISHED) THEN
              LOCAL_ERROR="The elements for mesh component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " have not been finished"
              FINISHED=.FALSE.
              EXIT
            ENDIF
          ELSE
            LOCAL_ERROR="The elements for mesh topology component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " are not associated"
            FINISHED=.FALSE.
            EXIT
          ENDIF
        ENDDO !component_idx
        IF(.NOT.FINISHED) CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        MESH%MESH_FINISHED=.TRUE.
      ELSE
        CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh user number       = ",MESH%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",MESH%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ",MESH%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_CREATE_FINISH",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh 
  SUBROUTINE MESH_CREATE_START_GENERIC(MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<The pointer to the meshes
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,mesh_idx
    TYPE(MESH_TYPE), POINTER :: NEW_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(NEW_MESH)
    NULLIFY(NEW_MESHES)

    ENTERS("MESH_CREATE_START_GENERIC",ERR,ERROR,*997)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*997)
      ELSE
        CALL MESH_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Set default mesh values
        NEW_MESH%USER_NUMBER=USER_NUMBER
        NEW_MESH%GLOBAL_NUMBER=MESHES%NUMBER_OF_MESHES+1
        NEW_MESH%MESHES=>MESHES
        NEW_MESH%NUMBER_OF_DIMENSIONS=NUMBER_OF_DIMENSIONS
        NEW_MESH%NUMBER_OF_COMPONENTS=1
        NEW_MESH%SURROUNDING_ELEMENTS_CALCULATE=.true. !default true
        !Initialise mesh topology and decompositions
        CALL MeshTopologyInitialise(NEW_MESH,ERR,ERROR,*999)
        CALL DECOMPOSITIONS_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Add new mesh into list of meshes 
        ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES+1),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
        DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
          NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
        ENDDO !mesh_idx
        NEW_MESHES(MESHES%NUMBER_OF_MESHES+1)%PTR=>NEW_MESH
        IF(ASSOCIATED(MESHES%MESHES)) DEALLOCATE(MESHES%MESHES)
        MESHES%MESHES=>NEW_MESHES
        MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES+1
        MESH=>NEW_MESH
      ENDIF
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*997)
    ENDIF
      
    EXITS("MESH_CREATE_START_GENERIC")
    RETURN
999 CALL MESH_FINALISE(NEW_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    NULLIFY(MESH)    
997 ERRORSEXITS("MESH_CREATE_START_GENERIC",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_GENERIC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in an interface. \see OPENCMISS::CMISSMeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_INTERFACE(USER_NUMBER,INTERFACE,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(INTERFACE%MESHES)) THEN
          CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on interface number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(INTERFACE%INTERFACES)) THEN
              PARENT_REGION=>INTERFACE%INTERFACES%PARENT_REGION
              IF(ASSOCIATED(PARENT_REGION)) THEN
                IF(ASSOCIATED(PARENT_REGION%COORDINATE_SYSTEM)) THEN                  
                  IF(NUMBER_OF_DIMENSIONS>0) THEN
                    IF(NUMBER_OF_DIMENSIONS<=PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                      CALL MESH_CREATE_START_GENERIC(INTERFACE%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                      MESH%INTERFACE=>INTERFACE
                    ELSE
                      LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & ") must be <= number of parent region dimensions ("// &
                        & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Parent region coordinate system is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interfaces parent region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface interfaces is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_START_INTERFACE")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_INTERFACE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in the region identified by REGION. \see OPENCMISS::CMISSMeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_REGION(USER_NUMBER,REGION,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(REGION%MESHES)) THEN
          CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
              IF(NUMBER_OF_DIMENSIONS>0) THEN
                IF(NUMBER_OF_DIMENSIONS<=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                  CALL MESH_CREATE_START_GENERIC(REGION%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                  MESH%REGION=>REGION
                ELSE
                  LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                    & ") must be <= number of region dimensions ("// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The coordinate system on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
                & " are not associated."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_CREATE_START_REGION")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_REGION",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_REGION

  !
  !================================================================================================================================
  !

  !>Destroys the mesh identified by a user number on the given region and deallocates all memory. \see OPENCMISS::CMISSMeshDestroy
  SUBROUTINE MESH_DESTROY_NUMBER(USER_NUMBER,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to destroy
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN

!!TODO: have a mesh_destory_ptr and mesh_destroy_number
        
        !Find the problem identified by the user number
        FOUND=.FALSE.
        mesh_position=0
        DO WHILE(mesh_position<REGION%MESHES%NUMBER_OF_MESHES.AND..NOT.FOUND)
          mesh_position=mesh_position+1
          IF(REGION%MESHES%MESHES(mesh_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          MESH=>REGION%MESHES%MESHES(mesh_position)%PTR

          CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

          !Remove the mesh from the list of meshes
          IF(REGION%MESHES%NUMBER_OF_MESHES>1) THEN
            ALLOCATE(NEW_MESHES(REGION%MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
            DO mesh_idx=1,REGION%MESHES%NUMBER_OF_MESHES
              IF(mesh_idx<mesh_position) THEN
                NEW_MESHES(mesh_idx)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ELSE IF(mesh_idx>mesh_position) THEN
                REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
                NEW_MESHES(mesh_idx-1)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ENDIF
            ENDDO !mesh_idx
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%MESHES=>NEW_MESHES
            REGION%MESHES%NUMBER_OF_MESHES=REGION%MESHES%NUMBER_OF_MESHES-1
          ELSE
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%NUMBER_OF_MESHES=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys the mesh and deallocates all memory. \see OPENCMISS::CMISSMeshDestroy
  SUBROUTINE MESH_DESTROY(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    TYPE(MESHES_TYPE), POINTER :: MESHES
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      MESHES=>MESH%MESHES
      IF(ASSOCIATED(MESHES)) THEN
        mesh_position=MESH%GLOBAL_NUMBER
          
        CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

        !Remove the mesh from the list of meshes
        IF(MESHES%NUMBER_OF_MESHES>1) THEN
          ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new meshes.",ERR,ERROR,*999)
          DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
            IF(mesh_idx<mesh_position) THEN
              NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ELSE IF(mesh_idx>mesh_position) THEN
              MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
              NEW_MESHES(mesh_idx-1)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ENDIF
          ENDDO !mesh_idx
          DEALLOCATE(MESHES%MESHES)
          MESHES%MESHES=>NEW_MESHES
          MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES-1
        ELSE
          DEALLOCATE(MESHES%MESHES)
          MESHES%NUMBER_OF_MESHES=0
        ENDIF
      ELSE
        CALL FlagError("The mesh meshes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a mesh and deallocates all memory.
  SUBROUTINE MESH_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL MeshTopologyFinalise(MESH,ERR,ERROR,*999)
      CALL DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*999)
!      IF(ASSOCIATED(MESH%INTF)) CALL INTERFACE_MESH_FINALISE(MESH,ERR,ERROR,*999)  ! <<??>>
      DEALLOCATE(MESH)
    ENDIF
 
    EXITS("MESH_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a region nodes pointer corresponding to the mesh global nodes accounting for interfaces.
  SUBROUTINE MeshGlobalNodesGet(mesh,nodes,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the global nodes for
    TYPE(NODES_TYPE), POINTER :: nodes !<On return, the nodes pointer corresponding to the global nodes for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshGlobalNodesGet",err,error,*999)
    
    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(nodes)) THEN
        CALL FlagError("Nodes is already associated.",err,error,*999)
      ELSE
        NULLIFY(nodes)
        region=>mesh%region
        IF(ASSOCIATED(region)) THEN
          nodes=>region%nodes
        ELSE
          INTERFACE=>mesh%INTERFACE
          IF(ASSOCIATED(interface)) THEN
            nodes=>interface%nodes
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have an associated region or interface."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
        IF(.NOT.ASSOCIATED(nodes)) THEN
          IF(ASSOCIATED(region)) THEN
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have any nodes associated with the mesh region."
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have any nodes associated with the mesh interface."
          ENDIF
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
        
    EXITS("MeshGlobalNodesGet")
    RETURN
999 ERRORSEXITS("MeshGlobalNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshGlobalNodesGet
  
  !
  !================================================================================================================================
  !

  !>Initialises a mesh.
  SUBROUTINE MESH_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(MESH,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate new mesh.",ERR,ERROR,*999)
      MESH%USER_NUMBER=0
      MESH%GLOBAL_NUMBER=0
      MESH%MESH_FINISHED=.FALSE.
      NULLIFY(MESH%MESHES)
      NULLIFY(MESH%REGION)
      NULLIFY(MESH%INTERFACE)
      NULLIFY(MESH%GENERATED_MESH)
      MESH%NUMBER_OF_DIMENSIONS=0
      MESH%NUMBER_OF_COMPONENTS=0
      MESH%MESH_EMBEDDED=.FALSE.
      NULLIFY(MESH%EMBEDDING_MESH)
      MESH%NUMBER_OF_EMBEDDED_MESHES=0
      NULLIFY(MESH%EMBEDDED_MESHES)
      MESH%NUMBER_OF_ELEMENTS=0
      NULLIFY(MESH%TOPOLOGY)
      NULLIFY(MESH%DECOMPOSITIONS)
    ENDIF
    
    EXITS("MESH_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the number of mesh components for a mesh identified by a pointer. \see OPENCMISS::CMISSMeshNumberOfComponentsGet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of components for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COMPONENTS !<On return, the number of components in the specified mesh.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_COMPONENTS=MESH%NUMBER_OF_COMPONENTS
      ELSE
        CALL FlagError("Mesh has not finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of mesh components for a mesh. \see OPENCMISS::CMISSMeshNumberOfComponentsSet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of components for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MeshComponentTopologyPtrType), POINTER :: NEW_TOPOLOGY(:)

    NULLIFY(NEW_TOPOLOGY)
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_COMPONENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FlagError("Mesh has been finished",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_COMPONENTS/=MESH%NUMBER_OF_COMPONENTS) THEN
            ALLOCATE(NEW_TOPOLOGY(NUMBER_OF_COMPONENTS),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new topology",ERR,ERROR,*999)
            IF(NUMBER_OF_COMPONENTS<MESH%NUMBER_OF_COMPONENTS) THEN
              DO component_idx=1,NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
            ELSE !NUMBER_OF_COMPONENTS>MESH%NUMBER_OF_COMPONENTS
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
!!TODO \todo sort out mesh_topology initialise/finalise so that they allocate and deal with this below then call that routine
              DO component_idx=MESH%NUMBER_OF_COMPONENTS+1,NUMBER_OF_COMPONENTS
                ALLOCATE(NEW_TOPOLOGY(component_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate new topology component",ERR,ERROR,*999)
                NEW_TOPOLOGY(component_idx)%PTR%mesh=>mesh
                NEW_TOPOLOGY(component_idx)%PTR%meshComponentNumber=component_idx
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%elements)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%dataPoints)
                !Initialise the topology components
                CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
              ENDDO !component_idx
            ENDIF
            IF(ASSOCIATED(MESH%TOPOLOGY)) DEALLOCATE(MESH%TOPOLOGY)
            MESH%TOPOLOGY=>NEW_TOPOLOGY
            MESH%NUMBER_OF_COMPONENTS=NUMBER_OF_COMPONENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of mesh components ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
          & ") is illegal. You must have >0 mesh components"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_SET")
    RETURN
!!TODO: tidy up memory deallocation on error
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET

  !
  !================================================================================================================================
  !
  
  !>Gets the number of elements for a mesh identified by a pointer. \see OPENCMISS::CMISSMeshNumberOfElementsGet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_ELEMENTS !<On return, the number of elements in the specified mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS
      ELSE
        CALL FlagError("Mesh has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of elements for a mesh. \see OPENCMISS::CMISSMeshNumberOfElementsSet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of elements for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS !<The number of elements to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_ELEMENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_ELEMENTS/=MESH%NUMBER_OF_ELEMENTS) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
                  IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
                    IF(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS>0) THEN
!!TODO: Reallocate the elements and copy information. 
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  CALL FlagError("Mesh topology component pointer is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
            ENDIF
            MESH%NUMBER_OF_ELEMENTS=NUMBER_OF_ELEMENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of elements ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ELEMENTS,"*",ERR,ERROR))// &
          & ") is invalid. You must have > 0 elements."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET

  !
  !================================================================================================================================
  !

  !>Returns the region for a mesh accounting for regions and interfaces
  SUBROUTINE MeshRegionGet(mesh,region,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the region for
    TYPE(REGION_TYPE), POINTER :: region !<On return, the meshes region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: parentRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshRegionGet",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(region)) THEN
        CALL FlagError("Region is already associated.",err,error,*999)
      ELSE
        NULLIFY(region)
        NULLIFY(interface)
        region=>mesh%region
        IF(.NOT.ASSOCIATED(region)) THEN
          interface=>mesh%interface
          IF(ASSOCIATED(interface)) THEN
            parentRegion=>interface%PARENT_REGION
            IF(ASSOCIATED(parentRegion)) THEN
              region=>parentRegion
            ELSE
              localError="The parent region not associated for mesh number "// &
                & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//" of interface number "// &
                & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The region or interface is not associated for mesh number "// &
              & TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF

    EXITS("MeshRegionGet")
    RETURN
999 ERRORSEXITS("MeshRegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshRegionGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the surrounding elements calculate flag. \see OPENCMISS::CMISSMeshSurroundingElementsCalculateSet
  SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET(MESH,SURROUNDING_ELEMENTS_CALCULATE_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the surrounding elements calculate flag for
    LOGICAL, INTENT(IN) :: SURROUNDING_ELEMENTS_CALCULATE_FLAG !<The surrounding elements calculate flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
      ELSE
        MESH%SURROUNDING_ELEMENTS_CALCULATE=SURROUNDING_ELEMENTS_CALCULATE_FLAG
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET")
    RETURN
999 ERRORSEXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET
  
  !
  !===============================================================================================================================
  !

  !>Finishes the process of creating elements for a specified mesh component in a mesh topology. \see OPENCMISS::CMISSMeshElementsCreateFinish
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology

    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Mesh elements have already been finished.",ERR,ERROR,*999)
      ELSE        
        ELEMENTS%ELEMENTS_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      meshComponentTopology=>elements%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(MESH)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of global elements = ",MESH%NUMBER_OF_ELEMENTS, &
            & ERR,ERROR,*999)
          DO ne=1,MESH%NUMBER_OF_ELEMENTS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element = ",ne,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER, &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ",ELEMENTS%ELEMENTS(ne)%USER_NUMBER, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(ELEMENTS%ELEMENTS(ne)%BASIS)) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis number         = ",ELEMENTS%ELEMENTS(ne)%BASIS% &
                & USER_NUMBER,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES)) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)% BASIS%NUMBER_OF_NODES,8,8, &
                & ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES,'("    User element nodes   =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("User element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES)) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
                & ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES,'("    Global element nodes =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("Global element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !ne
        ELSE
          CALL FlagError("Mesh component topology mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh elements mesh component topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating elements in the mesh component identified by MESH and component_idx. The elements will be created with a default basis of BASIS. ELEMENTS is the returned pointer to the MESH_ELEMENTS data structure. \see OPENCMISS::CMISSMeshElementsCreateStart
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,MESH_COMPONENT_NUMBER,BASIS,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to start creating the elements on
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the default basis to use
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the created mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,INSERT_STATUS,ne
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(MESH)) THEN     
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FlagError("Elements is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
              IF(ASSOCIATED(ELEMENTS%ELEMENTS)) THEN
                CALL FlagError("Mesh topology already has elements associated",ERR,ERROR,*998)
              ELSE
                IF(ASSOCIATED(BASIS)) THEN
                  MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%meshComponentNumber=MESH_COMPONENT_NUMBER
                  ALLOCATE(ELEMENTS%ELEMENTS(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate individual elements",ERR,ERROR,*999)
                  ELEMENTS%NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS !Psuedo inheritance of the number of elements
                  CALL Tree_CreateStart(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  CALL Tree_InsertTypeSet(ELEMENTS%ELEMENTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                  CALL Tree_CreateFinish(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  ELEMENTS%ELEMENTS_FINISHED=.FALSE.
                  !Set up the default values and allocate element structures
                  DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                    CALL MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER=ne
                    ELEMENTS%ELEMENTS(ne)%USER_NUMBER=ne
                    CALL Tree_ItemInsert(ELEMENTS%ELEMENTS_TREE,ne,ne,INSERT_STATUS,ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%BASIS=>BASIS
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate user element nodes",ERR,ERROR,*999)
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES=1
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES=1
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                      & BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes versions",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS = 1
                  ENDDO !ne
                ELSE
                  CALL FlagError("Basis is not associated",ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Mesh topology elements is not associated",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Mesh topology is not associated",ERR,ERROR,*998)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START")
    RETURN
999 CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,DUMMY_ERR,DUMMY_ERROR,*998)
998 NULLIFY(ELEMENTS)
    ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys the elements in a mesh topology. \todo as this is a user routine it should take a mesh pointer like create start and finish? Split this into destroy and finalise?
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to destroy 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODE_VERSIONS)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODE_VERSIONS)
    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODES)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%GLOBAL_ELEMENT_NODES)) DEALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES)
    
    EXITS("MESH_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh elements for a given mesh component.
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_GET(MESH,MESH_COMPONENT_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the elements for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to get the elements for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("MESH_TOPOLOGY_ELEMENTS_GET",ERR,ERROR,*998)
    
    IF(ASSOCIATED(MESH)) THEN
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FlagError("Elements is already associated.",ERR,ERROR,*998)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
            ELSE
              CALL FlagError("Mesh topology elements is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_GET")
    RETURN
999 NULLIFY(ELEMENTS)
998 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_GET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    NULLIFY(ELEMENT%BASIS)
    
    EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: Take user number
  
  !>Gets the basis for a mesh element identified by a given global number. \todo should take user number \see OPENCMISS::CMISSMeshElementsBasisGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to get the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to get the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<On return, a pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
          BASIS=>ELEMENT%BASIS
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the basis for a mesh element identified by a given global number. \todo should take user number
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: NEW_USER_ELEMENT_NODES(:),NEW_GLOBAL_ELEMENT_NODES(:),NEW_USER_ELEMENT_NODE_VERSIONS(:,:)
    INTEGER(INTG) :: OVERLAPPING_NUMBER_NODES,OVERLAPPING_NUMBER_DERIVATIVES
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(ASSOCIATED(BASIS)) THEN
            ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
            IF(ELEMENT%BASIS%NUMBER_OF_NODES/=BASIS%NUMBER_OF_NODES.OR. &
                & ELEMENT%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES/=BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES) THEN
              !Allocate new user and global element nodes
              ALLOCATE(NEW_USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                & BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element node versions",ERR,ERROR,*999)

              OVERLAPPING_NUMBER_NODES=MIN(BASIS%NUMBER_OF_NODES,ELEMENT%BASIS%NUMBER_OF_NODES)
              OVERLAPPING_NUMBER_DERIVATIVES=MIN(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,ELEMENT%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES)

              !Set default values
              NEW_USER_ELEMENT_NODE_VERSIONS=1
              NEW_USER_ELEMENT_NODES(OVERLAPPING_NUMBER_NODES+1:)=0
              NEW_GLOBAL_ELEMENT_NODES(OVERLAPPING_NUMBER_NODES+1:)=0
              !Copy previous values
              NEW_USER_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%USER_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)
              NEW_GLOBAL_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%GLOBAL_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)
              NEW_USER_ELEMENT_NODE_VERSIONS(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)= &
                & ELEMENT%USER_ELEMENT_NODE_VERSIONS(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)

              !Replace arrays with new ones
              CALL MOVE_ALLOC(NEW_USER_ELEMENT_NODE_VERSIONS,ELEMENT%USER_ELEMENT_NODE_VERSIONS)
              CALL MOVE_ALLOC(NEW_USER_ELEMENT_NODES,ELEMENT%USER_ELEMENT_NODES)
              CALL MOVE_ALLOC(NEW_GLOBAL_ELEMENT_NODES,ELEMENT%GLOBAL_ELEMENT_NODES)
            ENDIF            
            ELEMENT%BASIS=>BASIS
          ELSE
            CALL FlagError("Basis is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET

  !
  !================================================================================================================================
  !

  !>Gets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(OUT) :: USER_ELEMENT_NODES(:) !<On return, USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)>=SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1)) THEN
            USER_ELEMENT_NODES=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES
          ELSE
            LOCAL_ERROR="The size of USER_ELEMENT_NODES is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(USER_ELEMENT_NODES,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1),"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesSet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODES(:) !<USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nn,NUMBER_OF_BAD_NODES,GLOBAL_NODE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:),BAD_NODES(:)
    LOGICAL :: ELEMENT_NODES_OK,NODE_EXISTS
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)==ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>MESH%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>MESH%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%PARENT_REGION
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(ASSOCIATED(NODES)) THEN
                  ELEMENT_NODES_OK=.TRUE.
                  ALLOCATE(GLOBAL_ELEMENT_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate global element nodes.",ERR,ERROR,*999)
                  ALLOCATE(BAD_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate bad nodes.",ERR,ERROR,*999)
                  NUMBER_OF_BAD_NODES=0
                  DO nn=1,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES
                    CALL NODE_CHECK_EXISTS(NODES,USER_ELEMENT_NODES(nn),NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                    IF(NODE_EXISTS) THEN
                      GLOBAL_ELEMENT_NODES(nn)=GLOBAL_NODE_NUMBER
                    ELSE
                      NUMBER_OF_BAD_NODES=NUMBER_OF_BAD_NODES+1
                      BAD_NODES(NUMBER_OF_BAD_NODES)=USER_ELEMENT_NODES(nn)
                      ELEMENT_NODES_OK=.FALSE.
                    ENDIF
                  ENDDO !nn
                  IF(ELEMENT_NODES_OK) THEN
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES=USER_ELEMENT_NODES
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES=GLOBAL_ELEMENT_NODES
                  ELSE
                    IF(NUMBER_OF_BAD_NODES==1) THEN
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in region "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
                          & " of parent region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))
                      DO nn=2,NUMBER_OF_BAD_NODES-1
                        LOCAL_ERROR=LOCAL_ERROR//","//TRIM(NUMBER_TO_VSTRING(BAD_NODES(nn),"*",ERR,ERROR))
                      ENDDO !nn
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" of parent region number "// &
                          &  TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ENDIF
                    ENDIF
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  IF(ASSOCIATED(REGION)) THEN                   
                    CALL FlagError("The elements mesh region does not have any associated nodes.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("The elements mesh interface does not have any associated nodes.",ERR,ERROR,*999) 
                  ENDIF
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of element nodes does not match number of basis nodes for this element.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets an element node's version for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesSet
  SUBROUTINE MeshElements_ElementNodeVersionSet(GLOBAL_NUMBER,ELEMENTS,VERSION_NUMBER,DERIVATIVE_NUMBER, &
      & USER_ELEMENT_NODE_INDEX,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The version number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODE_INDEX !< The node index of the specified element node to set a version for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementNodeVersionSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(USER_ELEMENT_NODE_INDEX>=1.AND.USER_ELEMENT_NODE_INDEX<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>mesh%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>mesh%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%PARENT_REGION
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(DERIVATIVE_NUMBER>=1.AND.DERIVATIVE_NUMBER<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS% &
                  & NUMBER_OF_DERIVATIVES(USER_ELEMENT_NODE_INDEX)) THEN !Check if the specified derivative exists
                  IF(VERSION_NUMBER>=1) THEN !Check if the specified version is greater than 1
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODE_VERSIONS(DERIVATIVE_NUMBER,USER_ELEMENT_NODE_INDEX) & 
                      & = VERSION_NUMBER
                    !\todo : There is redunancy in USER_ELEMENT_NODE_VERSIONS since it was allocated in MESH_TOPOLOGY_ELEMENTS_CREATE_START based on MAXIMUM_NUMBER_OF_DERIVATIVES for that elements basis:ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
                  ELSE
                    LOCAL_ERROR="The specified node version number of "//TRIM(NUMBER_TO_VSTRING(VERSION_NUMBER,"*", & 
                      & ERR,ERROR))//" is invalid. The element node index should be greater than 1."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified node derivative number of "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*", & 
                    & ERR,ERROR))//" is invalid. The element node derivative index should be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_DERIVATIVES( &
                    & USER_ELEMENT_NODE_INDEX),"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",err,error,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified element node index of "//TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NODE_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The element node index should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementNodeVersionSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementNodeVersionSet",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementNodeVersionSet

  !
  !================================================================================================================================
  !

  !>Finalises the elements data structures for a mesh topology and deallocates any memory. \todo pass in elements
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
        CALL MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
      ENDDO !ne
      DEALLOCATE(ELEMENTS%ELEMENTS)
      IF(ASSOCIATED(ELEMENTS%ELEMENTS_TREE)) CALL Tree_Destroy(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
      DEALLOCATE(ELEMENTS)
    ENDIF
 
    EXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Mesh already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%meshComponentTopology=>TOPOLOGY
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Mesh already has topology data points associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%totalNumberOfProjectedData=0
        TOPOLOGY%dataPoints%meshComponentTopology=>topology
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: ditto.
  
  !>Gets the user number for a global element identified by a given global number. \todo Check that the user number doesn't already exist. \see OPENCMISS::CMISSMeshElementsUserNumberGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the elements to get the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET

  !
  !================================================================================================================================
  !

  !>Returns the user number for a global element identified by a given global number. \see OPENCMISS::CMISSMeshElementsUserNumberGet
  SUBROUTINE MeshElements_ElementUserNumberGet(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberGet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN          
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Elements have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberGet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberGet",ERR,ERROR)
    RETURN 1

   
  END SUBROUTINE MeshElements_ElementUserNumberGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a global element identified by a given global number. \see OPENCMISS::CMISSMeshElementsUserNumberSet
  SUBROUTINE MeshElements_ElementUserNumberSet(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to set.
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the element to set
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,INSERT_STATUS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          NULLIFY(TREE_NODE)
          CALL Tree_Search(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
          IF(ASSOCIATED(TREE_NODE)) THEN
            CALL Tree_NodeValueGet(ELEMENTS%ELEMENTS_TREE,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
            LOCAL_ERROR="Element user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " is already used by global element number "// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            CALL Tree_ItemDelete(ELEMENTS%ELEMENTS_TREE,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
            CALL Tree_ItemInsert(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
            ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberSet",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MeshElements_ElementUserNumberSet
  
  !
  !================================================================================================================================
  !

  !>Changes/sets the user numbers for all elements.
  SUBROUTINE MeshTopologyElementsUserNumbersAllSet(elements,userNumbers,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: elements !<A pointer to the elements to set all the user numbers for 
    INTEGER(INTG), INTENT(IN) :: userNumbers(:) !<The user numbers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,insertStatus
    TYPE(TREE_TYPE), POINTER :: newElementsTree
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newElementsTree)

    ENTERS("MeshTopologyElementsUserNumbersAllSet",err,error,*999)

    IF(ASSOCIATED(elements)) THEN
      IF(elements%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",err,error,*999)
      ELSE
        IF(elements%NUMBER_OF_ELEMENTS==SIZE(userNumbers,1)) THEN
          !Check the users numbers to ensure that there are no duplicates          
          CALL Tree_CreateStart(newElementsTree,err,error,*999)
          CALL Tree_InsertTypeSet(newElementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
          CALL Tree_CreateFinish(newElementsTree,err,error,*999)
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            CALL Tree_ItemInsert(newElementsTree,userNumbers(elementIdx),elementIdx,insertStatus,err,error,*999)
            IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) THEN
              localError="The specified user number of "//TRIM(NumberToVstring(userNumbers(elementIdx),"*",err,error))// &
                & " for global element number "//TRIM(NUMBER_TO_VSTRING(elementIdx,"*",err,error))// &
                & " is a duplicate. The user element numbers must be unique."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !elementIdx
          CALL Tree_Destroy(elements%ELEMENTS_TREE,err,error,*999)
          elements%ELEMENTS_TREE=>newElementsTree
          NULLIFY(newElementsTree)
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            elements%ELEMENTS(elementIdx)%GLOBAL_NUMBER=elementIdx
            elements%ELEMENTS(elementIdx)%USER_NUMBER=userNumbers(elementIdx)
          ENDDO !elementIdx
        ELSE
          localError="The number of specified element user numbers ("// &
            TRIM(NumberToVstring(SIZE(userNumbers,1),"*",err,error))// &
            ") does not match number of elements ("// &
            TRIM(NumberToVstring(elements%NUMBER_OF_ELEMENTS,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyElementsUserNumbersAllSet")
    RETURN
999 IF(ASSOCIATED(newElementsTree)) CALL Tree_Destroy(newElementsTree,err,error,*998)
998 ERRORSEXITS("MeshTopologyElementsUserNumbersAllSet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopologyElementsUserNumbersAllSet
  
  !
  !================================================================================================================================
  !

  !>Calculates the data points in the given mesh topology.
  SUBROUTINE MeshTopologyDataPointsCalculateProjection(mesh,dataProjection,err,error,*)
  
    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh topology to calcualte the data projection for
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection !<A pointer to the data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DATA_POINTS_TYPE), POINTER :: dataPoints !<A pointer to the data points
    TYPE(MeshDataPointsType), POINTER :: dataPointsTopology
    TYPE(DATA_PROJECTION_RESULT_TYPE), POINTER :: dataProjectionResult
    TYPE(MeshElementsType), POINTER :: elements
    INTEGER(INTG) :: dataPointIdx,elementIdx,countIdx,projectionNumber,globalCountIdx,elementNumber

    ENTERS("MeshTopologyDataPointsCalculateProjection",ERR,ERROR,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(dataProjection%DATA_PROJECTION_FINISHED) THEN 
        dataPoints=>dataProjection%DATA_POINTS
        !Default the first mesh component topology to contain data points ! \TODO: need to be changed once the data points topology is moved under meshTopologyType.
        dataPointsTopology=>mesh%TOPOLOGY(1)%PTR%dataPoints
        !Extract the global number of the data projection 
        projectionNumber=dataProjection%GLOBAL_NUMBER
        !Hard code the first mesh component since element topology is the same for all mesh components
        !\TODO: need to be changed once the elements topology is moved under meshTopologyType.
        elements=>mesh%TOPOLOGY(1)%PTR%ELEMENTS
        ALLOCATE(dataPointsTopology%elementDataPoint(elements%NUMBER_OF_ELEMENTS),STAT=ERR)     
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology element.",ERR,ERROR,*999)
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          dataPointsTopology%elementDataPoint(elementIdx)%elementNumber=elements%ELEMENTS(elementIdx)%GLOBAL_NUMBER
          dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData=0
        ENDDO        
        !Calculate number of projected data points on an element
        DO dataPointIdx=1,dataPoints%NUMBER_OF_DATA_POINTS
          dataProjectionResult=>dataProjection%DATA_PROJECTION_RESULTS(dataPointIdx)
          elementNumber=dataProjectionResult%ELEMENT_NUMBER
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
              dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData= &
                & dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData+1;
            ENDIF
          ENDDO !elementIdx
        ENDDO       
        !Allocate memory to store data indices and initialise them to be zero   
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          ALLOCATE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(dataPointsTopology% &
            & elementDataPoint(elementIdx)%numberOfProjectedData),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate data points topology element data points.",ERR,ERROR,*999)
          DO countIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=0
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=0
          ENDDO
        ENDDO     
        !Record the indices of the data that projected on the elements 
        globalCountIdx=0
        dataPointsTopology%totalNumberOfProjectedData=0
        DO dataPointIdx=1,dataPoints%NUMBER_OF_DATA_POINTS 
          dataProjectionResult=>dataProjection%DATA_PROJECTION_RESULTS(dataPointIdx)
          elementNumber=dataProjectionResult%ELEMENT_NUMBER
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            countIdx=1         
            IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
              globalCountIdx=globalCountIdx+1
              !Find the next data point index in this element
              DO WHILE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber/=0)
                countIdx=countIdx+1
              ENDDO
              dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=dataPointIdx
              dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=dataPointIdx!globalCountIdx (used this if only projected data are taken into account)
              dataPointsTopology%totalNumberOfProjectedData=dataPointsTopology%totalNumberOfProjectedData+1
            ENDIF             
          ENDDO !elementIdx
        ENDDO !dataPointIdx
        !Allocate memory to store total data indices in ascending order and element map
        ALLOCATE(dataPointsTopology%dataPoints(dataPointsTopology%totalNumberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology data points.",ERR,ERROR,*999)
        !The global number for the data points will be looping through elements.
        countIdx=1  
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          DO dataPointIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%dataPoints(countIdx)%userNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%userNumber
             dataPointsTopology%dataPoints(countIdx)%globalNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%globalNumber
             dataPointsTopology%dataPoints(countIdx)%elementNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
               & elementNumber
             countIdx=countIdx+1
          ENDDO !dataPointIdx
        ENDDO !elementIdx                      
      ELSE
        CALL FlagError("Data projection is not finished.",err,error,*999)
      ENDIF     
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyDataPointsCalculateProjection")
    RETURN
999 ERRORSEXITS("MeshTopologyDataPointsCalculateProjection",err,error)
    RETURN 1
  END SUBROUTINE MeshTopologyDataPointsCalculateProjection

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. \todo pass in the mesh topology
  SUBROUTINE MeshTopologyFinalise(mesh,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx

    ENTERS("MeshTopologyFinalise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
        CALL MeshTopologyComponentFinalise(mesh%topology(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      DEALLOCATE(mesh%topology)
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
 
    EXITS("MeshTopologyFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopologyFinalise

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. 
  SUBROUTINE MeshTopologyComponentFinalise(meshComponent,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: meshComponent !<A pointer to the mesh component to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopologyComponentFinalise",err,error,*999)

    IF(ASSOCIATED(meshComponent)) THEN
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(meshComponent%elements,err,error,*999)
      DEALLOCATE(meshComponent)
    ENDIF
 
    EXITS("MeshTopologyComponentFinalise")
    RETURN
999 ERRORSEXITS("MeshTopologyComponentFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopologyComponentFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given mesh. \todo finalise on error
  SUBROUTINE MeshTopologyInitialise(mesh,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to initialise the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    
    ENTERS("MeshTopologyInitialise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(mesh%topology)) THEN
        CALL FlagError("Mesh already has topology associated.",err,error,*999)
      ELSE
        !Allocate mesh topology
        ALLOCATE(mesh%topology(mesh%NUMBER_OF_COMPONENTS),STAT=err)
        IF(err/=0) CALL FlagError("Mesh topology could not be allocated.",err,error,*999)
        DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
          ALLOCATE(mesh%topology(componentIdx)%ptr,STAT=err)
          IF(err/=0) CALL FlagError("Mesh topology component could not be allocated.",err,error,*999)
          mesh%topology(componentIdx)%ptr%mesh=>mesh
          NULLIFY(mesh%topology(componentIdx)%ptr%elements)
          NULLIFY(mesh%topology(componentIdx)%ptr%dataPoints)
          !Initialise the topology components
          CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
        ENDDO !componentIdx
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyInitialise")
    RETURN
999 ERRORSEXITS("MeshTopologyInitialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopologyInitialise
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh component. 
  SUBROUTINE MeshTopologyElementCheckExistsMesh(mesh,meshComponentNumber,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to check the element exists on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component to check the element exits on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshElementsType), POINTER :: elements

    NULLIFY(elements)
    
    ENTERS("MeshTopologyElementCheckExistsMesh",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%MESH_FINISHED) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_GET(mesh,meshComponentNumber,elements,err,error,*999)
        CALL MeshTopologyElementCheckExistsMeshElements(elements,userElementNumber,elementExists,globalElementNumber,err,error,*999)
      ELSE
        CALL FlagError("Mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyElementCheckExistsMesh")
    RETURN
999 ERRORSEXITS("MeshTopologyElementCheckExistsMesh",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyElementCheckExistsMesh
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh elements. 
  SUBROUTINE MeshTopologyElementCheckExistsMeshElements(meshElements,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh elements, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("MeshTopologyElementCheckExistsMesh",err,error,*999)

    elementExists=.FALSE.
    globalElementNumber=0
    IF(ASSOCIATED(meshElements)) THEN
      NULLIFY(treeNode)
      CALL Tree_Search(meshElements%ELEMENTS_TREE,userElementNumber,treeNode,err,error,*999)
      IF(ASSOCIATED(treeNode)) THEN
        CALL Tree_NodeValueGet(meshElements%ELEMENTS_TREE,treeNode,globalElementNumber,err,error,*999)
        elementExists=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopologyElementCheckExistsMeshElements")
    RETURN
999 ERRORSEXITS("MeshTopologyElementCheckExistsMeshElements",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopologyElementCheckExistsMeshElements
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given list of MESHES. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,MESHES,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(MESHES_TYPE), POINTER :: MESHES !<The list of meshes containing the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx

    ENTERS("MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        mesh_idx=1
        DO WHILE(mesh_idx<=MESHES%NUMBER_OF_MESHES.AND..NOT.ASSOCIATED(MESH))
          IF(MESHES%MESHES(mesh_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            MESH=>MESHES%MESHES(mesh_idx)%PTR
          ELSE
            mesh_idx=mesh_idx+1
          ENDIF
        ENDDO
      ENDIF
    ELSE
      CALL FlagError("Meshes is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_USER_NUMBER_FIND_GENERIC")
    RETURN
999 ERRORSEXITS("MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_USER_NUMBER_FIND_GENERIC

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given INTERFACE. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_INTERFACE(USER_NUMBER,INTERFACE,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    ENTERS("MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_USER_NUMBER_FIND_INTERFACE")
    RETURN
999 ERRORSEXITS("MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MESH_USER_NUMBER_FIND_INTERFACE

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given REGION. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_REGION(USER_NUMBER,REGION,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: REGION !<The region containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("MESH_USER_NUMBER_FIND_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_USER_NUMBER_FIND_REGION")
    RETURN
999 ERRORSEXITS("MESH_USER_NUMBER_FIND_REGION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_USER_NUMBER_FIND_REGION

  !
  !================================================================================================================================
  !

  !>Finalises the meshes and deallocates all memory
  SUBROUTINE MESHES_FINALISE(MESHES,ERR,ERROR,*)

   !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to finalise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH
 
    ENTERS("MESHES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      DO WHILE(MESHES%NUMBER_OF_MESHES>0)
        MESH=>MESHES%MESHES(1)%PTR
        CALL MESH_DESTROY(MESH,ERR,ERROR,*999)
      ENDDO !mesh_idx
      DEALLOCATE(MESHES)
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("MESHES_FINALISE")
    RETURN
999 ERRORSEXITS("MESHES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESHES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the generic meshes.
  SUBROUTINE MESHES_INITIALISE_GENERIC(MESHES,ERR,ERROR,*)

    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("MESHES_INITIALISE_GENERIC",ERR,ERROR,*998)

    IF(ASSOCIATED(MESHES)) THEN
      CALL FlagError("Meshes is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(MESHES,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Meshes could not be allocated",ERR,ERROR,*999)
      NULLIFY(MESHES%REGION)
      NULLIFY(MESHES%INTERFACE)
      MESHES%NUMBER_OF_MESHES=0
      NULLIFY(MESHES%MESHES)
    ENDIF
    
    EXITS("MESHES_INITIALISE_GENERIC")
    RETURN
999 CALL MESHES_FINALISE(MESHES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("MESHES_INITIALISE_GENERIC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_GENERIC

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given interface.
  SUBROUTINE MESHES_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESHES)) THEN
        LOCAL_ERROR="Interface number "//TRIM(NumberToVString(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(INTERFACE%MESHES,ERR,ERROR,*999)
        INTERFACE%MESHES%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_INTERFACE")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_INTERFACE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_INTERFACE

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given region.
  SUBROUTINE MESHES_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NumberToVString(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(REGION%MESHES,ERR,ERROR,*999)
        REGION%MESHES%REGION=>REGION
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_REGION")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_REGION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_REGION

  !
  !================================================================================================================================
  !

  !>Initialises the embedded meshes.
  SUBROUTINE EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*)

    !Argument variables
    !TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to initialise
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EMBEDDED_MESH_INITIALISE",ERR,ERROR,*998)
    
    ALLOCATE(MESH_EMBEDDING,STAT=ERR)
    NULLIFY(MESH_EMBEDDING%PARENT_MESH)
    NULLIFY(MESH_EMBEDDING%CHILD_MESH)
    
    EXITS("EMBEDDED_MESH_INITIALISE")
    RETURN
!999 CALL EMBEDDED_MESH_FINALISE(MESH_EMBEDDING,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EMBEDDED_MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EMBEDDED_MESH_INITIALISE

  !
  !================================================================================================================================
  !

  !>Creates an embedding of one mesh in another
  SUBROUTINE MESH_EMBEDDING_CREATE(MESH_EMBEDDING, PARENT_MESH, CHILD_MESH,ERR,ERROR,*)
!    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: PARENT_MESH !<The parent mesh
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: CHILD_MESH  !<The child mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: NGP = 0, ne

    ENTERS("MESH_EMBEDDING_CREATE",ERR,ERROR,*999) 
    
    WRITE(*,*) 'parent mesh', PARENT_MESH%NUMBER_OF_ELEMENTS
    WRITE(*,*) 'child mesh', child_MESH%NUMBER_OF_ELEMENTS
    CALL EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*999)

    DO ne=1,PARENT_MESH%NUMBER_OF_ELEMENTS
      NGP = MAX(NGP,PARENT_MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%&
        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS)
    ENDDO !ne

    MESH_EMBEDDING%PARENT_MESH => PARENT_MESH
    MESH_EMBEDDING%CHILD_MESH  => CHILD_MESH
    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate child node positions.",ERR,ERROR,*999)
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(NGP,PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate gauss point positions.",ERR,ERROR,*999)
    
    EXITS("MESH_EMBEDDING_CREATE")
    RETURN 

999 ERRORSEXITS("MESH_EMBEDDING_CREATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_CREATE

  !
  !================================================================================================================================
  !

  !>Sets the positions of nodes in the child mesh for one element in the parent mesh
  SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION(MESH_EMBEDDING, ELEMENT_NUMBER, NODE_NUMBERS, XI_COORDS,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER  !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBERS(:) !<NODE_NUMBERS(node_idx) Node numbers in child mesh for the node_idx'th embedded node in the ELEMENT_NUMBER'th element of the parent mesh
    REAL(DP), INTENT(IN)      :: XI_COORDS(:,:)  !<XI_COORDS(:,node_idx) Xi coordinates of the node_idx'th embedded node in the ELEMENT_NUMBER'th

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR,*999)

    IF(ELEMENT_NUMBER<1 .OR. ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Element number out of range",ERR,ERROR,*999)
    ENDIF

    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NUMBER_OF_NODES = SIZE(NODE_NUMBERS)

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(SIZE(NODE_NUMBERS)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(1:SIZE(NODE_NUMBERS)) = NODE_NUMBERS(1:SIZE(NODE_NUMBERS))

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(SIZE(XI_COORDS,1),SIZE(XI_COORDS,2)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2)) = &
      & XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2))
    
    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  !
  !================================================================================================================================
  !

  !>Sets the positions of a Gauss point of the parent mesh in terms of element/xi coordinate in the child mesh
  SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA(MESH_EMBEDDING, PARENT_ELEMENT_NUMBER, GAUSSPT_NUMBER,&
    & PARENT_XI_COORD, CHILD_ELEMENT_NUMBER, CHILD_XI_COORD,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING   !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: PARENT_ELEMENT_NUMBER           !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: GAUSSPT_NUMBER                  !<Gauss point number in this element
    REAL(DP), INTENT(IN) :: PARENT_XI_COORD(:)              !<Xi coordinate in parent element

    INTEGER(INTG), INTENT(IN) :: CHILD_ELEMENT_NUMBER !<Element number in the child mesh
    REAL(DP), INTENT(IN) :: CHILD_XI_COORD(:)    !<Xi coordinate in child element

    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string

    ENTERS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR,*999)

    IF(PARENT_ELEMENT_NUMBER<1 .OR. PARENT_ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Parent element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(CHILD_ELEMENT_NUMBER<1 .OR. CHILD_ELEMENT_NUMBER > MESH_EMBEDDING%CHILD_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Child element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(GAUSSPT_NUMBER<1 .OR. GAUSSPT_NUMBER > SIZE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION,1)) THEN
      CALL FlagError("Gauss point number out of range",ERR,ERROR,*999)
    ENDIF

    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %PARENT_XI_COORD(SIZE(PARENT_XI_COORD)))
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %CHILD_XI_COORD(SIZE(CHILD_XI_COORD)))


    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD)) = &
      & PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD)) = &
      & CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%ELEMENT_NUMBER = CHILD_ELEMENT_NUMBER

    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA


  !
  !================================================================================================================================
  !

  !> Find the mesh with the given user number, or throw an error if it does not exist.
  SUBROUTINE MESH_USER_NUMBER_TO_MESH( USER_NUMBER, REGION, MESH, ERR, ERROR, * )
    !Arguments
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: REGION !<The region containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the mesh with the specified user number.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_USER_NUMBER_TO_MESH", ERR, ERROR, *999 )

    NULLIFY( MESH )
    CALL MESH_USER_NUMBER_FIND( USER_NUMBER, REGION, MESH, ERR, ERROR, *999 )
    IF( .NOT.ASSOCIATED( MESH ) ) THEN
      LOCAL_ERROR = "A mesh with an user number of "//TRIM(NumberToVString( USER_NUMBER, "*", ERR, ERROR ))// &
        & " does not exist on region number "//TRIM(NumberToVString( REGION%USER_NUMBER, "*", ERR, ERROR ))//"."
      CALL FlagError( LOCAL_ERROR, ERR, ERROR, *999 )
    ENDIF

    EXITS( "MESH_USER_NUMBER_TO_MESH" )
    RETURN
999 ERRORSEXITS( "MESH_USER_NUMBER_TO_MESH", ERR, ERROR )
    RETURN 1

  END SUBROUTINE MESH_USER_NUMBER_TO_MESH

  !
  !================================================================================================================================
  !
!!\todo THIS SHOULD REALLY BE MESH_USER_NUMBER_TO_DECOMPOSITION
  
  !> Find the decomposition with the given user number, or throw an error if it does not exist.
  SUBROUTINE DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION( USER_NUMBER, MESH, DECOMPOSITION, ERR, ERROR, * )
    !Arguments
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to find
    TYPE(MESH_TYPE), POINTER :: MESH !<The mesh containing the decomposition
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On exit, a pointer to the decomposition with the specified user number.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION", ERR, ERROR, *999 )

    NULLIFY( DECOMPOSITION )
    CALL DECOMPOSITION_USER_NUMBER_FIND( USER_NUMBER, MESH, DECOMPOSITION, ERR, ERROR, *999 )
    IF( .NOT.ASSOCIATED( DECOMPOSITION ) ) THEN
      LOCAL_ERROR = "A decomposition with an user number of "//TRIM(NumberToVString( USER_NUMBER, "*", ERR, ERROR ))// &
        & " does not exist on mesh number "//TRIM(NumberToVString( MESH%USER_NUMBER, "*", ERR, ERROR ))//"."
      CALL FlagError( LOCAL_ERROR, ERR, ERROR, *999 )
    ENDIF

    EXITS( "DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION" )
    RETURN
999 ERRORS( "DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION", ERR, ERROR )
    EXITS( "DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION")
    RETURN 1

  END SUBROUTINE DECOMPOSITION_USER_NUMBER_TO_DECOMPOSITION

  !
  !================================================================================================================================
  !

END MODULE MESH_ROUTINES
