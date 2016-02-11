!> \file
!> \author Chris Bradley
!> \brief This module contains all computational environment variables.
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

!> This module contains all computational environment variables.

MODULE COMP_ENVIRONMENT

  USE BASE_ROUTINES
  USE CMISS_MPI
  USE CmissPetsc
  USE CONSTANTS
  USE KINDS
  USE MPI
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
 
#include "macros.h"

  IMPLICIT NONE

  !Module parameters

  !Module types
  
  !>!>pointer type to COMPUTATIONAL_WORK_GROUP_TYPE
  TYPE, PUBLIC :: COMPUTATIONAL_WORK_GROUP_PTR_TYPE
     TYPE(COMPUTATIONAL_WORK_GROUP_TYPE), POINTER :: PTR
  END TYPE COMPUTATIONAL_WORK_GROUP_PTR_TYPE
  
  !>Contains information on logical working groups (added by Robert on 01/04/2010)
  TYPE, PUBLIC :: COMPUTATIONAL_WORK_GROUP_TYPE
    INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES !<size of the total compurational nodes belonging to this group
    INTEGER(INTG) :: NUMBER_SUB_WORK_GROUPS !<size of sub working grous
    TYPE(COMPUTATIONAL_WORK_GROUP_TYPE), POINTER:: PARENT !<Parent of this working groups
    TYPE(COMPUTATIONAL_WORK_GROUP_PTR_TYPE), ALLOCATABLE:: SUB_WORK_GROUPS(:) !<non-leaf node: The sub working groups
    TYPE(COMPUTATIONAL_ENVIRONMENT_TYPE), POINTER :: COMP_ENV !<pointer to the actual working environment
    LOGICAL :: COMP_ENV_FINISHED !<Is .TURE. if the actual working environment has been generated, .FALSE. if not
  END TYPE COMPUTATIONAL_WORK_GROUP_TYPE

  PRIVATE

  !>Contains information on a cache hierarchy
  TYPE CACHE_TYPE
    INTEGER(INTG) :: NUMBER_LEVELS !<The number of levels in the cache hierarchy
    INTEGER(INTG),ALLOCATABLE :: SIZE(:) !<SIZE(level_idx). The size of the level_idx'th cache level.
  END TYPE CACHE_TYPE

  !>Contains information on a computational node containing a number of processors
  TYPE COMPUTATIONAL_NODE_TYPE
    INTEGER(INTG) :: NUMBER_PROCESSORS !<The number of processors for this computational node
    INTEGER(INTG) :: RANK !<The MPI rank of this computational node
   !TYPE(CACHE_TYPE) :: CACHE 
    INTEGER(INTG) :: NODE_NAME_LENGTH !<The length of the name of the computational node
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: NODE_NAME !<The name of the computational node
  END TYPE COMPUTATIONAL_NODE_TYPE

  PUBLIC COMPUTATIONAL_NODE_TYPE

  !>Contains information on the MPI type to transfer information about a computational node
  TYPE MPI_COMPUTATIONAL_NODE_TYPE
    INTEGER(INTG) :: MPI_TYPE !<The MPI data type
    INTEGER(INTG) :: NUM_BLOCKS !<The number of blocks in the MPI data type. This will be equal to 4.
    INTEGER(INTG) :: BLOCK_LENGTHS(4) !<The length of each block.
    INTEGER(INTG) :: TYPES(4) !<The data types of each block.
    INTEGER(MPI_ADDRESS_KIND) :: DISPLACEMENTS(4) !<The address displacements to each block.
  END TYPE MPI_COMPUTATIONAL_NODE_TYPE

  !>Contains information on the computational environment the program is running in.
  TYPE COMPUTATIONAL_ENVIRONMENT_TYPE
    INTEGER(INTG) :: MPI_COMM !<The MPI communicator for cmiss
    INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES !<The number of computational nodes
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER !<The index of the running process
    TYPE(COMPUTATIONAL_NODE_TYPE), ALLOCATABLE :: COMPUTATIONAL_NODES(:) !<COMPUTATIONAL_NODES(node_idx). Contains information on the node_idx'th computational node. 
  END TYPE COMPUTATIONAL_ENVIRONMENT_TYPE

  !Module variables

  TYPE(COMPUTATIONAL_ENVIRONMENT_TYPE), TARGET :: COMPUTATIONAL_ENVIRONMENT !<The computational environment the program is running in.
  TYPE(MPI_COMPUTATIONAL_NODE_TYPE) :: MPI_COMPUTATIONAL_NODE_TYPE_DATA !<The MPI data on the computational nodes.

  !Interfaces
  ! Access specifiers for subroutines and interfaces(if any)
  PUBLIC COMPUTATIONAL_ENVIRONMENT_TYPE
  PUBLIC COMPUTATIONAL_ENVIRONMENT
  PUBLIC COMPUTATIONAL_ENVIRONMENT_INITIALISE,COMPUTATIONAL_ENVIRONMENT_FINALISE,COMPUTATIONAL_NODES_NUMBER_GET, &
    & COMPUTATIONAL_NODE_NUMBER_GET
  PUBLIC COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD, COMPUTATIONAL_WORK_GROUP_CREATE_START, COMPUTATIONAL_WORK_GROUP_CREATE_FINISH

CONTAINS
  !
  !================================================================================================================================
  !

  !>Add the work sub-group to the parent group based on the computational requirements (called by user)
  SUBROUTINE COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD(PARENT_WORK_GROUP, NUMBER_COMPUTATIONAL_NODES, &
   & ADDED_WORK_GROUP,ERR,ERROR,*)

    !Argument Variables
    TYPE(COMPUTATIONAL_WORK_GROUP_TYPE),POINTER, INTENT(INOUT) :: PARENT_WORK_GROUP
    TYPE(COMPUTATIONAL_WORK_GROUP_TYPE),POINTER, INTENT(INOUT) :: ADDED_WORK_GROUP
    INTEGER(INTG),INTENT(IN) :: NUMBER_COMPUTATIONAL_NODES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(COMPUTATIONAL_WORK_GROUP_PTR_TYPE) NEW_WORK_GROUP
    TYPE(COMPUTATIONAL_WORK_GROUP_TYPE),POINTER ::  TMP_PARENT_WORK_GROUP
    TYPE(COMPUTATIONAL_WORK_GROUP_PTR_TYPE), ALLOCATABLE :: SUB_WORK_GROUPS(:)
    INTEGER(INTG):: I

    ALLOCATE(NEW_WORK_GROUP%PTR)
    NEW_WORK_GROUP%PTR%NUMBER_COMPUTATIONAL_NODES = NUMBER_COMPUTATIONAL_NODES
    NEW_WORK_GROUP%PTR%NUMBER_SUB_WORK_GROUPS = 0

    IF(ASSOCIATED(PARENT_WORK_GROUP)) THEN 
      ALLOCATE(SUB_WORK_GROUPS(PARENT_WORK_GROUP%NUMBER_SUB_WORK_GROUPS+1))
      DO I=1,PARENT_WORK_GROUP%NUMBER_SUB_WORK_GROUPS
        SUB_WORK_GROUPS(I)%PTR=>PARENT_WORK_GROUP%SUB_WORK_GROUPS(I)%PTR
      ENDDO
      !SUB_WORK_GROUPS(1:PARENT_WORK_GROUP%NUMBER_SUB_WORK_GROUPS)=>PARENT_WORK_GROUP%SUB_WORK_GROUPS(:)
      
      IF(ALLOCATED(PARENT_WORK_GROUP%SUB_WORK_GROUPS)) THEN 
        DEALLOCATE(PARENT_WORK_GROUP%SUB_WORK_GROUPS)
      ENDIF
      SUB_WORK_GROUPS(1+PARENT_WORK_GROUP%NUMBER_SUB_WORK_GROUPS)%PTR=>NEW_WORK_GROUP%PTR
      ALLOCATE(PARENT_WORK_GROUP%SUB_WORK_GROUPS(SIZE(SUB_WORK_GROUPS,1)))
      DO I=1,SIZE(SUB_WORK_GROUPS,1)
        PARENT_WORK_GROUP%SUB_WORK_GROUPS(I)%PTR => SUB_WORK_GROUPS(I)%PTR
      ENDDO
      !PARENT_WORK_GROUP%SUB_WORK_GROUPS(:) => SUB_WORK_GROUPS(:)
      
      DEALLOCATE(SUB_WORK_GROUPS)
      PARENT_WORK_GROUP%NUMBER_SUB_WORK_GROUPS = 1+PARENT_WORK_GROUP%NUMBER_SUB_WORK_GROUPS
      NEW_WORK_GROUP%PTR%PARENT => PARENT_WORK_GROUP
      TMP_PARENT_WORK_GROUP => PARENT_WORK_GROUP 
      DO WHILE(ASSOCIATED(TMP_PARENT_WORK_GROUP)) !Update the computational number of its ancestors
        TMP_PARENT_WORK_GROUP%NUMBER_COMPUTATIONAL_NODES = TMP_PARENT_WORK_GROUP%NUMBER_COMPUTATIONAL_NODES &
          & + NEW_WORK_GROUP%PTR%NUMBER_COMPUTATIONAL_NODES
        TMP_PARENT_WORK_GROUP => TMP_PARENT_WORK_GROUP%PARENT
      ENDDO
    ELSE !Top level group
      CALL FlagError('PARENT_WORK_GROUP is not associated, call COMPUTATIONAL_WORK_GROUP_CREATE_START first',&
      & ERR,ERROR,*999)
    ENDIF
    ADDED_WORK_GROUP => NEW_WORK_GROUP%PTR

    EXITS("COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD

  !
  !================================================================================================================================
  !

  !>Create the highest level work group (Default: GROUP_WORLD)
  SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_START(WORLD_WORK_GROUP,NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*)

    !Argument Variables
    TYPE(COMPUTATIONAL_WORK_GROUP_TYPE),POINTER, INTENT(INOUT) :: WORLD_WORK_GROUP
    INTEGER(INTG),INTENT(IN) :: NUMBER_COMPUTATIONAL_NODES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(COMPUTATIONAL_WORK_GROUP_PTR_TYPE) NEW_WORK_GROUP

    IF(ASSOCIATED(WORLD_WORK_GROUP)) THEN 
      CALL FlagError('WORLD_WORK_GROUP is already associated', ERR, ERROR, *999)
    ELSE
      ALLOCATE(NEW_WORK_GROUP%PTR)
      NEW_WORK_GROUP%PTR%NUMBER_COMPUTATIONAL_NODES = NUMBER_COMPUTATIONAL_NODES
      NEW_WORK_GROUP%PTR%NUMBER_SUB_WORK_GROUPS = 0   
      NULLIFY(NEW_WORK_GROUP%PTR%PARENT) !It is the highest level work group already 
      NULLIFY(NEW_WORK_GROUP%PTR%COMP_ENV) !Generate this later in COMPUTATIONAL_WORK_GROUP_CREATE_FINISH
      WORLD_WORK_GROUP=>NEW_WORK_GROUP%PTR
    ENDIF
    
    EXITS("COMPUTATIONAL_WORK_GROUP_CREATE_START")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_WORK_GROUP_CREATE_START",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_START

  !
  !================================================================================================================================
  !

  !>Generate computational environment for current level work group tree and all it's subgroups recursively 
  RECURSIVE SUBROUTINE Computational_WorkGroupGenerateCompEnviron(WORK_GROUP,AVAILABLE_RANK_LIST,ERR,ERROR,*)

    !Argument Variables
    TYPE(COMPUTATIONAL_WORK_GROUP_TYPE),POINTER, INTENT(INOUT) :: WORK_GROUP
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: AVAILABLE_RANK_LIST(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: I,MPI_IERROR,RANK,ORIGINAL_GROUP,NEW_GROUP
    INTEGER(INTG), ALLOCATABLE :: NEW_AVAILABLE_RANK_LIST(:)
    
    ENTERS("Computational_WorkGroupGenerateCompEnviron",ERR,ERROR,*999)
    
    ALLOCATE(WORK_GROUP%COMP_ENV)

    !Set size of computational nodes in this communicator
    WORK_GROUP%COMP_ENV%NUMBER_COMPUTATIONAL_NODES = WORK_GROUP%NUMBER_COMPUTATIONAL_NODES
    
    !Determine my processes rank
    CALL MPI_COMM_RANK(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,RANK,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)
    WORK_GROUP%COMP_ENV%MY_COMPUTATIONAL_NODE_NUMBER=RANK
    
    !Fill in the information for every computational node in this group
    ALLOCATE(WORK_GROUP%COMP_ENV%COMPUTATIONAL_NODES(WORK_GROUP%COMP_ENV%NUMBER_COMPUTATIONAL_NODES))
    I=SIZE(AVAILABLE_RANK_LIST,1)
    IF(SIZE(AVAILABLE_RANK_LIST,1)-WORK_GROUP%COMP_ENV%NUMBER_COMPUTATIONAL_NODES < 0) THEN
      CALL FlagError("NOT ENOUGH RANKS", ERR, ERROR, *999)
      GOTO 999
    ENDIF
    DO I=1,WORK_GROUP%COMP_ENV%NUMBER_COMPUTATIONAL_NODES, 1
      WORK_GROUP%COMP_ENV%COMPUTATIONAL_NODES(I) = COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(AVAILABLE_RANK_LIST(I))
    ENDDO

    !Create a communicator
    !CALL MPI_COMM_DUP(MPI_COMM_WORLD,WORK_GROUP%COMP_ENV%MPI_COMM,MPI_IERROR)
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD,ORIGINAL_GROUP,MPI_IERROR);
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)    
    CALL MPI_GROUP_INCL(ORIGINAL_GROUP,I-1,AVAILABLE_RANK_LIST(1:I-1),NEW_GROUP,MPI_IERROR)  !Choose the first I-1 ranks
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)    
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD,NEW_GROUP,WORK_GROUP%COMP_ENV%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)    
    CALL MPI_GROUP_FREE(ORIGINAL_GROUP,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)    
    CALL MPI_GROUP_FREE(NEW_GROUP,MPI_IERROR) 
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)    

    !Shrink the AVAILABLE_RANK_LIST
    ALLOCATE(NEW_AVAILABLE_RANK_LIST(SIZE(AVAILABLE_RANK_LIST,1)-WORK_GROUP%COMP_ENV%NUMBER_COMPUTATIONAL_NODES))  
    NEW_AVAILABLE_RANK_LIST(1:SIZE(NEW_AVAILABLE_RANK_LIST)) = AVAILABLE_RANK_LIST(I:SIZE(AVAILABLE_RANK_LIST,1))
    DEALLOCATE(AVAILABLE_RANK_LIST)
    ALLOCATE(AVAILABLE_RANK_LIST(SIZE(NEW_AVAILABLE_RANK_LIST,1)))
    AVAILABLE_RANK_LIST(:) = NEW_AVAILABLE_RANK_LIST(:)

    WORK_GROUP%COMP_ENV_FINISHED = .TRUE.

    !Recursively do this to all its subgroups
    DO I=1,WORK_GROUP%NUMBER_SUB_WORK_GROUPS,1
      CALL Computational_WorkGroupGenerateCompEnviron(WORK_GROUP%SUB_WORK_GROUPS(I)%PTR,&
        & AVAILABLE_RANK_LIST,ERR,ERROR,*999)      
    ENDDO

    EXITS("Computational_WorkGroupGenerateCompEnviron")
    RETURN
999 ERRORSEXITS("Computational_WorkGroupGenerateCompEnviron",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Computational_WorkGroupGenerateCompEnviron

  !
  !================================================================================================================================
  !

  !>Generate the hierarchy computational environment based on work group tree
  SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_FINISH(WORLD_WORK_GROUP,ERR,ERROR,*)

    !Argument Variables
    TYPE(COMPUTATIONAL_WORK_GROUP_TYPE),POINTER,INTENT(INOUT) :: WORLD_WORK_GROUP
    INTEGER(INTG),INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG),ALLOCATABLE:: AVAILABLE_RANK_LIST(:)
    INTEGER(INTG) :: I

    ENTERS("COMPUTATIONAL_WORK_GROUP_CREATE_FINISH",ERR,ERROR,*999)

    !set the computational environment of the world work group to be the global COMPUTATIONAL_ENVIRONMENT (the default communicator in OpenCMISS)
    WORLD_WORK_GROUP%COMP_ENV => COMPUTATIONAL_ENVIRONMENT 
    WORLD_WORK_GROUP%COMP_ENV_FINISHED = .TRUE.

    !generate the communicators for subgroups if any
    ALLOCATE(AVAILABLE_RANK_LIST(WORLD_WORK_GROUP%COMP_ENV%NUMBER_COMPUTATIONAL_NODES))
    DO I=0,SIZE(AVAILABLE_RANK_LIST,1)-1
      AVAILABLE_RANK_LIST(I+1) = I
    END DO
    DO I=1,WORLD_WORK_GROUP%NUMBER_SUB_WORK_GROUPS,1
      CALL Computational_WorkGroupGenerateCompEnviron(WORLD_WORK_GROUP%SUB_WORK_GROUPS(I)%PTR,AVAILABLE_RANK_LIST,ERR,ERROR,*999)
    END DO

    EXITS("COMPUTATIONAL_WORK_GROUP_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_WORK_GROUP_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finalises the computational node data structures and deallocates all memory.
  SUBROUTINE COMPUTATIONAL_NODE_FINALISE(COMPUTATIONAL_NODE,ERR,ERROR,*)
  
    !Argument Variables
    TYPE(COMPUTATIONAL_NODE_TYPE),INTENT(INOUT) :: COMPUTATIONAL_NODE !<The computational node to finalise
    INTEGER(INTG),INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("COMPUTATIONAL_NODE_FINALISE",ERR,ERROR,*999)

    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=0
    COMPUTATIONAL_NODE%RANK=-1
    COMPUTATIONAL_NODE%NODE_NAME_LENGTH=0
    COMPUTATIONAL_NODE%NODE_NAME=""    

    EXITS("COMPUTATIONAL_NODE_FINALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the computational node data structures.
  SUBROUTINE COMPUTATIONAL_NODE_INITIALISE(COMPUTATIONAL_NODE,RANK,ERR,ERROR,*)
  
    !Argument Variables
    TYPE(COMPUTATIONAL_NODE_TYPE), INTENT(OUT) :: COMPUTATIONAL_NODE !<The computational node to initialise
    INTEGER(INTG), INTENT(IN) :: RANK !<The MPI rank of the computational node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MPI_IERROR

    ENTERS("COMPUTATIONAL_NODE_INITIALISE",ERR,ERROR,*999)

!    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=COMP_DETECT_NUMBER_PROCESSORS(ERR)
!    IF(ERR/=0) GOTO 999
    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=1
    COMPUTATIONAL_NODE%RANK=RANK
    CALL MPI_GET_PROCESSOR_NAME(COMPUTATIONAL_NODE%NODE_NAME,COMPUTATIONAL_NODE%NODE_NAME_LENGTH,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_PROCESSOR_NAME",MPI_IERROR,ERR,ERROR,*999)
    
    EXITS("COMPUTATIONAL_NODE_INITIALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the data structure containing the MPI type information for the COMPUTATIONAL_NODE_TYPE.
  SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(ERR,ERROR,*)
  
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,MPI_IERROR

    ENTERS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE",ERR,ERROR,*999)

    DO i=1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES(i)=0
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS(i)=0
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(i)=0
    ENDDO !i
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS=0

    IF(MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE/=MPI_DATATYPE_NULL) THEN
      CALL MPI_TYPE_FREE(MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_TYPE_FREE",MPI_IERROR,ERR,ERROR,*999)
    ENDIF

    EXITS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the data structure containing the MPI type information for the COMPUTATIONAL_NODE_TYPE.
  SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE(COMPUTATIONAL_NODE,ERR,ERROR,*)
  
    !Argument Variables
    TYPE(COMPUTATIONAL_NODE_TYPE), INTENT(IN) :: COMPUTATIONAL_NODE !<The computational node containing the MPI type to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: I,MPI_IERROR

    ENTERS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE",ERR,ERROR,*999)

    MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE=MPI_DATATYPE_NULL
    
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS=4
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES=(/MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_CHARACTER/)
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS=(/1,1,1,MPI_MAX_PROCESSOR_NAME/)
	
	
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NUMBER_PROCESSORS,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(1),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%RANK,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(2),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NODE_NAME_LENGTH,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(3),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)
    !CPB 19/02/07 AIX compiler complains about the type of the first parameter i.e., the previous 3 have been integers
    !and this one is not so cast the type.
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NODE_NAME,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(4),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)

    DO i=4,1,-1
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(I)=MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(I)- &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(1)
    ENDDO !i

    CALL MPI_TYPE_CREATE_STRUCT(MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS, &
      & MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS,MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES, &
      & MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_TYPE_CREATE_STRUCT",MPI_IERROR,ERR,ERROR,*999)

    CALL MPI_TYPE_COMMIT(MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",MPI_IERROR,ERR,ERROR,*999)
    
    IF(DIAGNOSTICS3) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI Computational Node Type Data:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  MPI type = ",MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number blocks  = ",MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,4,4, &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES,'("  Block types =",4(X,I15))','(15X,4(X,I15))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,8,8, &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS,'("  Block lengths =",8(X,I5))','(17X,8(X,I5))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,8,8, &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS,'("  Displacements =",8(X,I5))','(17X,8(X,I5))',ERR,ERROR,*999)
    ENDIF

    EXITS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE")
    RETURN
999 CALL COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(ERR,ERROR,*998)
998 ERRORSEXITS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the computational environment data structures and deallocates all memory.
  SUBROUTINE COMPUTATIONAL_ENVIRONMENT_FINALISE(ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: COMPUTATIONAL_NODE,MPI_IERROR

    ENTERS("COMPUTATIONAL_ENVIRONMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)) THEN
       DO COMPUTATIONAL_NODE=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
          CALL COMPUTATIONAL_NODE_FINALISE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(COMPUTATIONAL_NODE),ERR,ERROR,*999)
       ENDDO
       DEALLOCATE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)
    ENDIF
    COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES=0

    CALL COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(ERR,ERROR,*999)

    CALL MPI_COMM_FREE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_FREE",MPI_IERROR,ERR,ERROR,*999)

    !Finalise PetSc
    !Call this after MPI_COMM_FREE as PETSc routines are called when some
    !MPI comm attributes are freed.
    !CALL Petsc_LogView(PETSC_COMM_WORLD,"OpenCMISSTest.petsc",ERR,ERROR,*999)
    CALL Petsc_Finalise(ERR,ERROR,*999)

    CALL MPI_FINALIZE(MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_FINALIZE",MPI_IERROR,ERR,ERROR,*999)

    EXITS("COMPUTATIONAL_ENVIRONMENT_FINALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_ENVIRONMENT_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_ENVIRONMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the computational environment data structures.
  SUBROUTINE COMPUTATIONAL_ENVIRONMENT_INITIALISE(ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,DUMMY_ERR,MPI_IERROR,RANK
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("COMPUTATIONAL_ENVIRONMENT_INITIALISE",ERR,ERROR,*999)

    !Initialise the MPI environment
    CALL MPI_INIT(MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_INIT",MPI_IERROR,ERR,ERROR,*999)

    !Create a (private) communicator for cmiss. For now just duplicate MPI_COMM_WORLD
    CALL MPI_COMM_DUP(MPI_COMM_WORLD,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_DUP",MPI_IERROR,ERR,ERROR,*999)
    
    !Determine the number of ranks/computational nodes we have in our computational environment
    CALL MPI_COMM_SIZE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_SIZE",MPI_IERROR,ERR,ERROR,*999)

    !Allocate the computational node data structures
    ALLOCATE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=ERR)
    IF(ERR /=0) CALL FlagError("Could not allocate computational nodes",ERR,ERROR,*999)

    !Determine my processes rank
    CALL MPI_COMM_RANK(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,RANK,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)
    COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER=RANK
    
#ifdef TAUPROF
    CALL TAU_PROFILE_SET_NODE(rank)
#endif

    !Create the MPI type information for the COMPUTATIONAL_NODE_TYPE
    CALL COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(RANK),ERR,ERROR,*999)
    !Fill in all the computational node data structures for this rank at the root position (will be changed later with an
    !allgather call)
    CALL COMPUTATIONAL_NODE_INITIALISE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0),RANK,ERR,ERROR,*999)

!     !Now transfer all the computational node information to the other computational nodes so that each rank has all the
!     !information.
! !!    CALL MPI_ALLGATHER(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0),1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, &
! !!      & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0),1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, &
! !!      & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
!     CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, &
!       & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0),1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, &
!       & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
!     CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Calling MPI_ERROR_CHECK...",ERR,ERROR,*999)
!     CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)

    !Initialise node numbers in base routines.
    CALL COMPUTATIONAL_NODE_NUMBERS_SET(COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER,COMPUTATIONAL_ENVIRONMENT% &
      & NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)
    
    !Initialise PETSc
    CALL Petsc_Initialise(PETSC_NULL_CHARACTER,ERR,ERROR,*999)
    
    IF(DIAGNOSTICS1) THEN
      !Just let the master node write out this information
      IF(RANK==0) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Computational environment:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of computational nodes = ", &
          & COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  My computational node number = ", &
          & COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER,ERR,ERROR,*999)
        IF(DIAGNOSTICS2) THEN
          DO i=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Computational Node:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of Processors = ", &
              & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(i)%NUMBER_PROCESSORS,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    MPI rank = ", &
             & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(i)%RANK,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node Name = ", &
              & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(i)%NODE_NAME,ERR,ERROR,*999)
          ENDDO !i
        ENDIF
      ENDIF
    ENDIF

    EXITS("COMPUTATIONAL_ENVIRONMENT_INITIALISE")
    RETURN
999 CALL COMPUTATIONAL_ENVIRONMENT_FINALISE(DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("COMPUTATIONAL_ENVIRONMENT_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_ENVIRONMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the number/rank of the computational nodes.
  FUNCTION COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: COMPUTATIONAL_NODE_NUMBER_GET !<On exit the computational node number/rank of the current process.
    !Local Variables

    ENTERS("COMPUTATIONAL_NODE_NUMBER_GET",ERR,ERROR,*999)

    IF(ALLOCATED(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)) THEN
      COMPUTATIONAL_NODE_NUMBER_GET=COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER
    ELSE
      CALL FlagError("Computational environment not initialised",ERR,ERROR,*999)
    ENDIF
    
    EXITS("COMPUTATIONAL_NODE_NUMBER_GET")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_NUMBER_GET",ERR,ERROR)
    RETURN 
  END FUNCTION COMPUTATIONAL_NODE_NUMBER_GET

  !
  !================================================================================================================================
  !
  
  !>Returns the number of computational nodes.
  FUNCTION COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
     
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: COMPUTATIONAL_NODES_NUMBER_GET !<On exit, the number of computational nodes for the program.
    !Local Variables

    ENTERS("COMPUTATIONAL_NODES_NUMBER_GET",ERR,ERROR,*999)

    IF(ALLOCATED(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)) THEN
      COMPUTATIONAL_NODES_NUMBER_GET=COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
    ELSE
      CALL FlagError("Computational environment not initialised",ERR,ERROR,*999)
    ENDIF
    
    EXITS("COMPUTATIONAL_NODES_NUMBER_GET")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODES_NUMBER_GET",ERR,ERROR)
    RETURN 
  END FUNCTION COMPUTATIONAL_NODES_NUMBER_GET

  !
  !================================================================================================================================
  !

END MODULE COMP_ENVIRONMENT
