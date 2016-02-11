!> \file
!> \author Chris Bradley
!> \brief This module handles all domain mappings routines.
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

!> This module handles all domain mappings routines.
MODULE DOMAIN_MAPPINGS

  USE BASE_ROUTINES
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MPI
  USE STRINGS
  USE TYPES

#include "macros.h"  


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  INTERFACE DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE
    !This one will be deleted eventually
    MODULE PROCEDURE DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE
    MODULE PROCEDURE DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE_ELEMENT
  END INTERFACE !DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE

  PUBLIC DOMAIN_MAPPINGS_MAPPING_FINALISE,DOMAIN_MAPPINGS_MAPPING_INITIALISE,DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE, &
    & DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE

CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Finalises the adjacent domain and deallocates all memory for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(ADJACENT_DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE) :: ADJACENT_DOMAIN !<The adjacent domain to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ADJACENT_DOMAIN%LOCAL_GHOST_SEND_INDICES)) DEALLOCATE(ADJACENT_DOMAIN%LOCAL_GHOST_SEND_INDICES)
    IF(ALLOCATED(ADJACENT_DOMAIN%LOCAL_GHOST_RECEIVE_INDICES)) DEALLOCATE(ADJACENT_DOMAIN%LOCAL_GHOST_RECEIVE_INDICES)
    ADJACENT_DOMAIN%NUMBER_OF_SEND_GHOSTS=0
    ADJACENT_DOMAIN%NUMBER_OF_RECEIVE_GHOSTS=0
    ADJACENT_DOMAIN%DOMAIN_NUMBER=0
    
    EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialise the adjacent domain for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(ADJACENT_DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE) :: ADJACENT_DOMAIN !<The adjacent domain to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",ERR,ERROR,*999)

    ADJACENT_DOMAIN%NUMBER_OF_SEND_GHOSTS=0
    ADJACENT_DOMAIN%NUMBER_OF_RECEIVE_GHOSTS=0
    ADJACENT_DOMAIN%DOMAIN_NUMBER=0
    
    EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the domain mappings local to global map for element-based dofs, i.e. nodes, Gauss points, etc.
  SUBROUTINE DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE_ELEMENT(domainMapping,elementDofs,elementDofOffsets, &
      & elementDomains,elementDomainOffsets,err,error,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: domainMapping !<The domain mapping to calculate the local mappings
    INTEGER(INTG), INTENT(INOUT) :: elementDofs(:) !<elementDofs(elementDofOffsets(elementIdx)+elementDofIdx-1). On entry, the map from the elementDofIdx'th element dof index in the elementIdx'th local element to a unique local dof number. On exit this map is renumbered such that the locally owned dofs are smaller than or equal to the number of local dofs.
    INTEGER(INTG), INTENT(IN) :: elementDofOffsets(:) !<elementDofOffsets(elementIdx). The offset of the elementIdx'th local element into the elementDofs array.
    INTEGER(INTG), INTENT(IN) :: elementDomains(:) !<elementDomains(elementDomainOffsets(elementIdx)+elementDomainIdx-1). The elementDomainIdx'th domain that shares the elementIdx'th local element. The first elementDomainIdx for each elementIdx is the owner.
    INTEGER(INTG), INTENT(IN) :: elementDomainOffsets(:) !<elementDomainOffsets(elementIdx). The offset of the elementIdx'th local element into the elementDomains array.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: localFlag=0,ghostFlag=-1
    INTEGER(INTG) :: myDomain,numberOfDomains,numberOfLocalElements,totalNumberOfLocalElements,elementIdx,elementOwner, &
      & elementDomainIdx,elementDofIdx,numberOfLocal,totalNumberOfLocal,numberOfGlobal,numberOfGhost,localCount, &
      & domainIdx,domainNumber,numberOfElementDofs,numberOfAdjacentDomains,adjacentDomainIdx,sendBufferIdx,receiveBufferIdx, &
      & localNumber,globalNumber,ghostIdx,dummyErr,MPI_IERROR
    INTEGER(INTG), ALLOCATABLE :: localMap(:),globalMap(:),requests(:),statuses(:,:),numberOfLocalPerDomain(:),domainOffsets(:), &
      & sendCounts(:),receiveCounts(:),adjacentDomainNumbers(:)
    TYPE(INTEGER_INTG_ALLOC_TYPE), ALLOCATABLE :: sendBuffers(:),receiveBuffers(:)
    TYPE(LIST_TYPE), POINTER :: ghostSendList,ghostReceiveList
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE_ELEMENT",err,error,*999)

    IF(.NOT.ASSOCIATED(domainMapping)) CALL FlagError("Domain mapping is not associated.",err,error,*999)
    
    myDomain=COMPUTATIONAL_NODE_NUMBER_GET(err,error)+1
    IF(err/=0) CALL FlagError("Could not get computation node number.",err,error,*999)
    numberOfDomains=domainMapping%NUMBER_OF_DOMAINS

    !The element dofs are contiguously numbered from 1 to the total number of local dofs.
    totalNumberOfLocal=MAXVAL(elementDofs)
    !Store the map from the dofs as given in the element dofs array to the new local numbering, such that all the local dofs come before the ghost dofs.
    ALLOCATE(localMap(totalNumberOfLocal),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate local map.",err,error,*999)
    !Initialise with the ghost flag. These represent dofs that are only on ghost elements (and not shared with local elements).
    localMap=ghostFlag
    !Store the map from the dofs as given in the element dofs array to the global numbering
    ALLOCATE(globalMap(totalNumberOfLocal),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate global map.",err,error,*999)
    globalMap=ghostFlag

    totalNumberOfLocalElements=SIZE(elementDomainOffsets)-1
    !Count the number of elements that has myDomain as element owner
    numberOfLocalElements=0
    DO elementIdx=1,totalNumberOfLocalElements
      elementOwner=elementDomains(elementDomainOffsets(elementIdx))
      IF(elementOwner==myDomain) THEN
        numberOfLocalElements=numberOfLocalElements+1
      ELSE
        EXIT
      END IF
    END DO !elementIdx

    !Flag the dofs on local elements with the local flag.
    DO elementIdx=1,numberOfLocalElements
      DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
        localMap(elementDofs(elementDofIdx))=localFlag
      END DO !elementDofIdx
    END DO !elementIdx
    !Flag the dofs that are shared with ghost elements whose element owner is lower than my domain number to the ghost flag again.
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      elementOwner=elementDomains(elementDomainOffsets(elementIdx))
      IF(elementOwner<myDomain) THEN 
        DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
          localMap(elementDofs(elementDofIdx))=ghostFlag
        END DO !elementDofIdx
      END IF
    END DO !elementIdx
    !Number the local dofs contiguously
    !First number all the local dofs
    localCount=0
    DO elementIdx=1,numberOfLocalElements
      DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
        IF(localMap(elementDofs(elementDofIdx))==localFlag) THEN
          localCount=localCount+1
          localMap(elementDofs(elementDofIdx))=localCount
        END IF
      END DO !elementDofIdx
    END DO !elementIdx

    numberOfLocal=numberOfLocal
    numberOfGhost=totalNumberOfLocal-numberOfLocal

    !Now number all the ghost dofs
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
        IF(localMap(elementDofs(elementDofIdx))==ghostFlag) THEN
          localCount=localCount+1
          localMap(elementDofs(elementDofIdx))=localCount
        END IF
      END DO !elementDofIdx
    END DO !elementIdx

    !Gather number of local for all domains.
    ALLOCATE(numberOfLocalPerDomain(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate number of local.",err,error,*999)
    CALL MPI_ALLGATHER(numberOfLocal,1,MPI_INTEGER,numberOfLocalPerDomain,1,MPI_INTEGER, &
      & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,err,error,*999)
    !Need this to calculate the ghost send and receive indices.
    ALLOCATE(domainOffsets(domainMapping%NUMBER_OF_DOMAINS+1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate domain offsets.",err,error,*999)
    domainOffsets(1)=1
    DO domainIdx=1,numberOfDomains
      domainOffsets(domainIdx+1)=domainOffsets(domainIdx)+numberOfLocalPerDomain(domainIdx)
    END DO !domainIdx
    !\todo: correct this for all constant dofs

    IF(ALLOCATED(numberOfLocalPerDomain)) DEALLOCATE(numberOfLocalPerDomain)
    !Calculate all the local (excluding ghost) numbers based on my domain offset.
    WHERE(localMap<=numberOfLocal) globalMap=localMap+domainOffsets(myDomain)-1
    numberOfGlobal=domainOffsets(numberOfDomains+1)-1

    !Calculate the send counts
    ALLOCATE(sendCounts(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate send counts.",err,error,*999)
    sendCounts=0
    !Count all the dofs from the boundary elements from this domain.
    DO elementIdx=1,numberOfLocalElements
      numberOfElementDofs=elementDofOffsets(elementIdx+1)-elementDofOffsets(elementIdx)
      DO elementDomainIdx=elementDomainOffsets(elementIdx)+1,elementDomainOffsets(elementIdx+1)-1
        domainNumber=elementDomains(elementDomainIdx)
        sendCounts(domainNumber)=sendCounts(domainNumber)+numberOfElementDofs
      END DO !elementDomainIdx
    END DO !elementIdx
    !Count all the dofs in the ghost elements whose element owners are greater than this domain.
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      elementOwner=elementDomains(elementDomainOffsets(elementIdx))
      numberOfElementDofs=elementDofOffsets(elementIdx+1)-elementDofOffsets(elementIdx)
      DO elementDomainIdx=elementDomainOffsets(elementIdx)+1,elementDomainOffsets(elementIdx+1)-1
        domainNumber=elementDomains(elementDomainIdx)
        IF(myDomain<elementOwner.AND.domainNumber/=myDomain) &
          & sendCounts(domainNumber)=sendCounts(domainNumber)+numberOfElementDofs
      END DO !elementDomainIdx
    END DO !elementIdx

    !Calculate the receive counts
    ALLOCATE(receiveCounts(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate receive counts.",err,error,*999)
    receiveCounts=0
    !Count all the dofs coming from the element owners of the ghost elements. 
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      numberOfElementDofs=elementDofOffsets(elementIdx+1)-elementDofOffsets(elementIdx)
      domainNumber=elementDomains(elementDomainOffsets(elementIdx))
      receiveCounts(domainNumber)=receiveCounts(domainNumber)+numberOfElementDofs
    END DO !elementIdx
    !Count all the dofs coming from domain numbers that are smaller than the element owner of the ghost elements. 
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      elementOwner=elementDomains(elementDomainOffsets(elementIdx))
      numberOfElementDofs=elementDofOffsets(elementIdx+1)-elementDofOffsets(elementIdx)
      DO elementDomainIdx=elementDomainOffsets(elementIdx)+1,elementDomainOffsets(elementIdx+1)-1
        domainNumber=elementDomains(elementDomainIdx)
        IF(domainNumber<elementOwner.AND.domainNumber/=myDomain) &
          & receiveCounts(domainNumber)=receiveCounts(domainNumber)+numberOfElementDofs
      END DO !elementDomainIdx
    END DO !elementIdx

    !Allocate and fill the send buffers 
    ALLOCATE(sendBuffers(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate send buffers.",err,error,*999)
    DO domainIdx=1,numberOfDomains
      ALLOCATE(sendBuffers(domainIdx)%array(sendCounts(domainIdx)),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate send buffer.",err,error,*999)
    END DO !domainIdx

    sendCounts=0
    DO elementIdx=1,numberOfLocalElements
      DO elementDomainIdx=elementDomainOffsets(elementIdx)+1,elementDomainOffsets(elementIdx+1)-1
        domainNumber=elementDomains(elementDomainIdx)
        DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
          sendCounts(domainNumber)=sendCounts(domainNumber)+1
          sendBuffers(domainNumber)%array(sendCounts(domainNumber))=globalMap(elementDofs(elementDofIdx))
        END DO !elementDofIdx
      END DO !elementDomainIdx
    END DO !elementIdx
    !Send all the dofs in the ghost elements whose element owners are greater than this domain.
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      elementOwner=elementDomains(elementDomainOffsets(elementIdx))
      DO elementDomainIdx=elementDomainOffsets(elementIdx)+1,elementDomainOffsets(elementIdx+1)-1
        domainNumber=elementDomains(elementDomainIdx)
        IF(myDomain<elementOwner.AND.domainNumber/=myDomain) THEN
          DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
            sendCounts(domainNumber)=sendCounts(domainNumber)+1
            sendBuffers(domainNumber)%array(sendCounts(domainNumber))=globalMap(elementDofs(elementDofIdx))
          END DO !elementDofIdx
        END IF
      END DO !elementDomainIdx
    END DO !elementIdx

    !Allocate all the receive buffers.
    ALLOCATE(receiveBuffers(numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate receive buffers.",err,error,*999)
    DO domainIdx=1,numberOfDomains
      ALLOCATE(receiveBuffers(domainIdx)%array(receiveCounts(domainIdx)),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate receive buffer.",err,error,*999)
      receiveBuffers(domainIdx)%array=0 
    END DO

    ALLOCATE(statuses(MPI_STATUS_SIZE,2*numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate mpi statuses.",err,error,*999)
    ALLOCATE(requests(2*numberOfDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate mpi receive requests.",err,error,*999)
    requests=MPI_REQUEST_NULL

    !Send and receive the buffers.
    DO domainIdx=1,numberOfDomains
      IF(receiveCounts(domainIdx)>0) THEN
        CALL MPI_IRECV(receiveBuffers(domainIdx)%array,receiveCounts(domainIdx),MPI_INTEGER,domainIdx-1, &
          & MPI_LOCAL_TO_GLOBAL_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,requests(domainIdx),MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,err,error,*999)
      END IF
    END DO
    DO domainIdx=1,numberOfDomains
      IF(sendCounts(domainIdx)>0) THEN
        CALL MPI_ISEND(sendBuffers(domainIdx)%array,sendCounts(domainIdx),MPI_INTEGER,domainIdx-1, &
          & MPI_LOCAL_TO_GLOBAL_TAG,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,requests(numberOfDomains+domainIdx),MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,err,error,*999)
      END IF
    END DO

    !While sending and receiving we can set up some parts of the adjacent domains.
    numberOfAdjacentDomains=COUNT((sendCounts+receiveCounts)>0)
    ALLOCATE(adjacentDomainNumbers(numberOfAdjacentDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate adjacent domain numbers.",err,error,*999)
    numberOfAdjacentDomains=0
    DO domainIdx=1,numberOfDomains
      IF((sendCounts(domainIdx)+receiveCounts(domainIdx))>0) THEN
        numberOfAdjacentDomains=numberOfAdjacentDomains+1
        adjacentDomainNumbers(numberOfAdjacentDomains)=domainIdx
      END IF
    END DO !domainIdx
    
    domainMapping%NUMBER_OF_LOCAL=numberOfLocal
    domainMapping%TOTAL_NUMBER_OF_LOCAL=totalNumberOfLocal
    domainMapping%NUMBER_OF_GLOBAL=numberOfGlobal

    domainMapping%NUMBER_OF_ADJACENT_DOMAINS=numberOfAdjacentDomains
    ALLOCATE(domainMapping%ADJACENT_DOMAINS(numberOfAdjacentDomains),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate domain mapping adjacent domains.",err,error,*999)
    DO adjacentDomainIdx=1,numberOfAdjacentDomains
      domainNumber=adjacentDomainNumbers(adjacentDomainIdx)
      CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx),err,error,*999)
      domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER=domainNumber
      NULLIFY(ghostSendList)
      CALL List_CreateStart(ghostSendList,err,error,*999)
      CALL List_DataTypeSet(ghostSendList,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(ghostSendList,numberOfGhost,err,error,*999)
      CALL List_CreateFinish(ghostSendList,err,error,*999)
      DO sendBufferIdx=1,sendCounts(domainNumber)
        globalNumber=sendBuffers(domainNumber)%array(sendBufferIdx)
        IF(globalNumber>0) THEN
          localNumber=globalNumber-domainOffsets(myDomain)+1
          CALL List_ItemAdd(ghostSendList,localNumber,err,error,*999)
        END IF
      END DO !sendBufferIdx
      CALL List_RemoveDuplicates(ghostSendList,err,error,*999)
      CALL List_DetachAndDestroy(ghostSendList,domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS, &
        & domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES,err,error,*999)
    END DO !adjacentDomainIdx

    CALL MPI_WAITALL(2*numberOfDomains,requests,statuses,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,err,error,*999)

    IF(DIAGNOSTICS5) THEN
      DO domainIdx=1,numberOfDomains
        IF(receiveCounts(domainIdx)>0) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI IRECV call posted:",err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive count = ",receiveCounts(domainIdx),err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive datatype = ",MPI_INTEGER,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive source = ",domainIdx-1,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",MPI_LOCAL_TO_GLOBAL_TAG,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive comm = ",COMPUTATIONAL_ENVIRONMENT%MPI_COMM,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive request = ",requests(domainIdx),err,error,*999)                
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,receiveCounts(domainIdx),8,8,receiveBuffers(domainIdx)%array, &
            & '("      Receive data       :",8(X,I10))','(39X,8(X,I10))',err,error,*999)      
        END IF
      ENDDO !domainIdx
      DO domainIdx=1,numberOfDomains
        IF(sendCounts(domainIdx)>0) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI ISEND call posted:",err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send count = ",sendCounts(domainIdx),err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send datatype = ",MPI_INTEGER,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send dest = ",domainIdx-1,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",MPI_LOCAL_TO_GLOBAL_TAG,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send comm = ",COMPUTATIONAL_ENVIRONMENT%MPI_COMM,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send request = ",requests(numberOfDomains+domainIdx),err,error,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,sendCounts(domainIdx),8,8,sendBuffers(domainIdx)%array, &
            & '("      Send data       :",8(X,I10))','(39X,8(X,I10))',err,error,*999)      
        END IF
      ENDDO !domainIdx
    END IF

    IF(ALLOCATED(requests)) DEALLOCATE(requests)
    IF(ALLOCATED(statuses)) DEALLOCATE(statuses)

    !Complete the global map

    receiveCounts=0
    !Get all the dofs coming from the element owners of the ghost elements. 
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      domainNumber=elementDomains(elementDomainOffsets(elementIdx))
      DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
        receiveCounts(domainNumber)=receiveCounts(domainNumber)+1
        IF(receiveBuffers(domainNumber)%array(receiveCounts(domainNumber))>0) &
          & globalMap(elementDofs(elementDofIdx))=receiveBuffers(domainNumber)%array(receiveCounts(domainNumber))
      END DO !elementDofIdx
    END DO !elementIdx
    !Get all the dofs coming from domain numbers that are smaller than the element owner of the ghost elements. 
    DO elementIdx=numberOfLocalElements+1,totalNumberOfLocalElements
      elementOwner=elementDomains(elementDomainOffsets(elementIdx))
      DO elementDomainIdx=elementDomainOffsets(elementIdx)+1,elementDomainOffsets(elementIdx+1)-1
        domainNumber=elementDomains(elementDomainIdx)
        IF(domainNumber<elementOwner.AND.domainNumber/=myDomain) THEN
          DO elementDofIdx=elementDofOffsets(elementIdx),elementDofOffsets(elementIdx+1)-1
            receiveCounts(domainNumber)=receiveCounts(domainNumber)+1
            IF(receiveBuffers(domainNumber)%array(receiveCounts(domainNumber))>0) &
              & globalMap(elementDofs(elementDofIdx))=receiveBuffers(domainNumber)%array(receiveCounts(domainNumber))
          END DO !elementDofIdx
        END IF
      END DO !elementDomainIdx
    END DO !elementIdx

    ALLOCATE(domainMapping%LOCAL_TO_GLOBAL_MAP(totalNumberOfLocal),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate domain local to global map.",err,error,*999)
    domainMapping%LOCAL_TO_GLOBAL_MAP(localMap)=globalMap
    IF(ALLOCATED(globalMap)) DEALLOCATE(globalMap)

    !Calculate the ghost receive indices
    DO adjacentDomainIdx=1,numberOfAdjacentDomains
      domainNumber=adjacentDomainNumbers(adjacentDomainIdx)
      NULLIFY(ghostReceiveList)
      CALL List_CreateStart(ghostReceiveList,err,error,*999)
      CALL List_DataTypeSet(ghostReceiveList,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(ghostReceiveList,numberOfGhost,err,error,*999)
      CALL List_CreateFinish(ghostReceiveList,err,error,*999)
      DO receiveBufferIdx=1,receiveCounts(domainNumber)
        globalNumber=receiveBuffers(domainNumber)%array(receiveBufferIdx)
        IF(globalNumber>0) THEN
          !Search for the global number in the local to global map (only the part with the ghost numbers is needed).
          CALL List_Search(domainMapping%LOCAL_TO_GLOBAL_MAP(numberOfLocal+1:),globalNumber,ghostIdx,err,error,*999)
          CALL List_ItemAdd(ghostReceiveList,numberOfLocal+ghostIdx,err,error,*999)
        END IF
      END DO !receiveBufferIdx
      CALL List_RemoveDuplicates(ghostReceiveList,err,error,*999)
      CALL List_DetachAndDestroy(ghostReceiveList,domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS, &
        & domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES,err,error,*999)
    END DO !adjacentDomainIdx

    IF(ALLOCATED(receiveCounts)) DEALLOCATE(receiveCounts)
    IF(ALLOCATED(sendCounts)) DEALLOCATE(sendCounts)
    !Deallocate send and receive buffers.
    IF(ALLOCATED(sendBuffers)) THEN
      DO domainIdx=1,SIZE(sendBuffers)
        IF(ALLOCATED(sendBuffers(domainIdx)%array)) DEALLOCATE(sendBuffers(domainIdx)%array)
      END DO ! domainIdx
      DEALLOCATE(sendBuffers)
    ENDIF
    IF(ALLOCATED(receiveBuffers)) THEN
      DO domainIdx=1,SIZE(receiveBuffers)
        IF(ALLOCATED(receiveBuffers(domainIdx)%array)) DEALLOCATE(receiveBuffers(domainIdx)%array)
      END DO ! domainIdx
      DEALLOCATE(receiveBuffers)
    ENDIF

    !Renumber the element dofs such that all local dofs come before ghost dofs.
    elementDofs=localMap(elementDofs)

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Domain mappings:",err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains  = ",domainMapping%NUMBER_OF_DOMAINS,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global = ",domainMapping%NUMBER_OF_GLOBAL,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of local = ",domainMapping%NUMBER_OF_LOCAL,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of local = ",domainMapping%TOTAL_NUMBER_OF_LOCAL, &
        & err,error,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map:",err,error,*999)
      DO localNumber=1,domainMapping%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local index : ",localNumber,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global index = ",domainMapping%LOCAL_TO_GLOBAL_MAP(localNumber), &
          & err,error,*999)
      ENDDO !localNumber
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains:",err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & domainMapping%NUMBER_OF_ADJACENT_DOMAINS,err,error,*999)
      DO adjacentDomainIdx=1,domainMapping%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",adjacentDomainIdx,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS,err,error,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)% &
          & NUMBER_OF_SEND_GHOSTS,8,8,domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indices       :",8(X,I10))','(39X,8(X,I10))',err,error,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of receive ghosts = ", &
          & domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS,err,error,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,8,8,domainMapping%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indices    :",8(X,I10))','(39X,8(X,I10))',err,error,*999)              
      ENDDO !adjacentDomainIdx
    ENDIF

    EXITS("DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE")
    RETURN
999 IF(ALLOCATED(sendBuffers)) THEN
      DO domainIdx=1,SIZE(sendBuffers)
        IF(ALLOCATED(sendBuffers(domainIdx)%array)) DEALLOCATE(sendBuffers(domainIdx)%array)
      END DO ! domainIdx
      DEALLOCATE(sendBuffers)
    ENDIF
    IF(ALLOCATED(receiveBuffers)) THEN
      DO domainIdx=1,SIZE(receiveBuffers)
        IF(ALLOCATED(receiveBuffers(domainIdx)%array)) DEALLOCATE(receiveBuffers(domainIdx)%array)
      END DO ! domainIdx
      DEALLOCATE(receiveBuffers)
    ENDIF
    CALL List_Destroy(ghostSendList,dummyErr,dummyError,*998)
998 CALL List_Destroy(ghostReceiveList,dummyErr,dummyError,*997)
997 ERRORSEXITS("DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE_ELEMENT",err,error)
    RETURN 1
    
  END SUBROUTINE DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE_ELEMENT
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the domain mappings local to global map.
  SUBROUTINE DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE(DOMAIN_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The domain mapping to calculate the local mappings
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,adjacentDomainIdx,localNumber,NUMBER_OF_GHOST_RECEIVE,NUMBER_OF_GHOST_SEND,domainOffset,myDomain,MPI_IERROR
    INTEGER(INTG), ALLOCATABLE :: REQUESTS(:),STATUSES(:,:),DOMAIN_NUMBER_OF_LOCAL(:)
    TYPE(INTEGER_INTG_ALLOC_TYPE), ALLOCATABLE :: SEND_BUFFERS(:),RECEIVE_BUFFERS(:)

    ENTERS("DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      myDomain=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)+1
      IF(ERR/=0) GOTO 999
    
      !Check that the number of local and total number of local is set.
      IF(DOMAIN_MAPPING%NUMBER_OF_LOCAL<=0.OR.DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL<=0) THEN
        CALL FlagError("Domain mapping number of local and total number of local is not set.",ERR,ERROR,*999)
      ELSE IF(DOMAIN_MAPPING%NUMBER_OF_LOCAL>DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
        CALL FlagError("Domain mapping number of local is greater than total number of local.",ERR,ERROR,*999)
      END IF

      !Count the local (excluding ghosts) numbers. There should be no number owned by more than one domain (this should be handled
      !during the mapping calculation of the nodes and dofs).
      ALLOCATE(DOMAIN_NUMBER_OF_LOCAL(DOMAIN_MAPPING%NUMBER_OF_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate number of local.",ERR,ERROR,*999)
      DOMAIN_NUMBER_OF_LOCAL=0

      !Gather number of local for all domains.
      CALL MPI_ALLGATHER(DOMAIN_MAPPING%NUMBER_OF_LOCAL,1,MPI_INTEGER,DOMAIN_NUMBER_OF_LOCAL,1,MPI_INTEGER, &
        & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
      DOMAIN_MAPPING%NUMBER_OF_GLOBAL=SUM(DOMAIN_NUMBER_OF_LOCAL)

      ALLOCATE(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate domain local to global map.",ERR,ERROR,*999)
      DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP=0

      !Calculate the local to global map.
      !First calculate all the local (excluding ghost) numbers and then communicate them to their adjacent domains.
      domainOffset=SUM(DOMAIN_NUMBER_OF_LOCAL(:myDomain-1))
      DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(:DOMAIN_MAPPING%NUMBER_OF_LOCAL)=[(i+domainOffset,i=1,DOMAIN_MAPPING%NUMBER_OF_LOCAL)]

      IF(ALLOCATED(DOMAIN_NUMBER_OF_LOCAL)) DEALLOCATE(DOMAIN_NUMBER_OF_LOCAL)

      !Send and receive the global numbers for the ghosts
      ALLOCATE(SEND_BUFFERS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate number of local.",ERR,ERROR,*999)
      ALLOCATE(RECEIVE_BUFFERS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate number of local.",ERR,ERROR,*999)
      ALLOCATE(REQUESTS(2*DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate mpi receive requests.",ERR,ERROR,*999)
      REQUESTS=MPI_REQUEST_NULL
      ALLOCATE(STATUSES(MPI_STATUS_SIZE,2*DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate mpi statuses.",ERR,ERROR,*999)

      !Check that the adjacent domains are allocated
      IF(.NOT.ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS)) &
        & CALL FlagError("Domain mapping adjacent domains is not allocated.",ERR,ERROR,*999)

      !Post all the receive calls first and then the send calls.
      DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        NUMBER_OF_GHOST_RECEIVE=DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS
        ALLOCATE(RECEIVE_BUFFERS(adjacentDomainIdx)%array(NUMBER_OF_GHOST_RECEIVE),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate receive buffer pointer.",ERR,ERROR,*999)
        !Check that the ghost receive indices array is allocated.
        IF(.NOT.ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES)) &
          & CALL FlagError("Domain mapping adjacent domains local ghost receive indices is not allocated.",ERR,ERROR,*999)
        RECEIVE_BUFFERS(adjacentDomainIdx)%array=0
        CALL MPI_IRECV(RECEIVE_BUFFERS(adjacentDomainIdx)%array,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS, &
          & MPI_INTEGER,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER-1,MPI_LOCAL_TO_GLOBAL_TAG, &
          & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,REQUESTS(adjacentDomainIdx),MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,ERR,ERROR,*999)
      ENDDO !adjacentDomainIdx
      !Post all the send calls.
      DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        NUMBER_OF_GHOST_SEND=DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS
        ALLOCATE(SEND_BUFFERS(adjacentDomainIdx)%array(NUMBER_OF_GHOST_SEND),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate send buffer pointer.",ERR,ERROR,*999)
        !Check that the ghost send indices array is allocated.
        IF(.NOT.ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES)) &
          & CALL FlagError("Domain mapping adjacent domains local ghost send indices is not allocated.",ERR,ERROR,*999)
        SEND_BUFFERS(adjacentDomainIdx)%array=DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
          & LOCAL_GHOST_SEND_INDICES)
        CALL MPI_ISEND(SEND_BUFFERS(adjacentDomainIdx)%array,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS, &
          & MPI_INTEGER,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER-1,MPI_LOCAL_TO_GLOBAL_TAG, &
          & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,REQUESTS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS+adjacentDomainIdx),MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,ERR,ERROR,*999)
      ENDDO !adjacentDomainIdx

      CALL MPI_WAITALL(2*DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,REQUESTS,STATUSES,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)

      IF(DIAGNOSTICS5) THEN
        DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI IRECV call posted:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive count = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive datatype = ",MPI_INTEGER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive source = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER-1,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",MPI_LOCAL_TO_GLOBAL_TAG,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive comm = ", &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive request = ", &
            & REQUESTS(adjacentDomainIdx),ERR,ERROR,*999)                
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
            & NUMBER_OF_RECEIVE_GHOSTS,8,8,RECEIVE_BUFFERS(adjacentDomainIdx)%array, &
            & '("      Receive data       :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)      
        ENDDO !adjacentDomainIdx
        DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI ISEND call posted:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send count = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send datatype = ",MPI_INTEGER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send dest = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER-1,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Receive tag = ",MPI_LOCAL_TO_GLOBAL_TAG,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send comm = ",COMPUTATIONAL_ENVIRONMENT%MPI_COMM,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Send request = ", &
            & REQUESTS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS+adjacentDomainIdx),ERR,ERROR,*999)                
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
            & NUMBER_OF_SEND_GHOSTS,8,8,SEND_BUFFERS(adjacentDomainIdx)%array, &
            & '("      Send data       :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)      
        ENDDO !adjacentDomainIdx
      END IF
      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Domain mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains  = ",DOMAIN_MAPPING%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global = ",DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of local = ",DOMAIN_MAPPING%NUMBER_OF_LOCAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of local = ",DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
          & DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
        DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",adjacentDomainIdx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
            & NUMBER_OF_SEND_GHOSTS,8,8,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES, &
            & '("      Local send ghost indices       :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)      
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of receive ghosts = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
            & NUMBER_OF_RECEIVE_GHOSTS,8,8,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES, &
            & '("      Local receive ghost indices    :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)              
        ENDDO !adjacentDomainIdx
      ENDIF

      !Copy the receive buffers back to the ghost positions in the local to global map.
      DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES)= &
          RECEIVE_BUFFERS(adjacentDomainIdx)%array
      ENDDO !adjacentDomainIdx


      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Domain mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains  = ",DOMAIN_MAPPING%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global = ",DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of local = ",DOMAIN_MAPPING%NUMBER_OF_LOCAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of local = ",DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map:",ERR,ERROR,*999)
        DO localNumber=1,DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local index : ",localNumber,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global index = ",DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localNumber), &
            & ERR,ERROR,*999)
        ENDDO !localNumber
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
          & DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
        DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",adjacentDomainIdx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%DOMAIN_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
            & NUMBER_OF_SEND_GHOSTS,8,8,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES, &
            & '("      Local send ghost indices       :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)      
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of receive ghosts = ", &
            & DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)% &
            & NUMBER_OF_RECEIVE_GHOSTS,8,8,DOMAIN_MAPPING%ADJACENT_DOMAINS(adjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES, &
            & '("      Local receive ghost indices    :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)              
        ENDDO !adjacentDomainIdx
      ENDIF

      !Deallocate send and receive buffers.
      DO adjacentDomainIdx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        IF(ALLOCATED(SEND_BUFFERS(adjacentDomainIdx)%array)) DEALLOCATE(SEND_BUFFERS(adjacentDomainIdx)%array)
        IF(ALLOCATED(RECEIVE_BUFFERS(adjacentDomainIdx)%array)) DEALLOCATE(RECEIVE_BUFFERS(adjacentDomainIdx)%array)
      END DO !adjacentDomainIdx

      IF(ALLOCATED(SEND_BUFFERS)) DEALLOCATE(SEND_BUFFERS)
      IF(ALLOCATED(RECEIVE_BUFFERS)) DEALLOCATE(RECEIVE_BUFFERS)
      IF(ALLOCATED(REQUESTS)) DEALLOCATE(REQUESTS)
      IF(ALLOCATED(STATUSES)) DEALLOCATE(STATUSES)
      
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAIN_NUMBER_OF_LOCAL)) DEALLOCATE(DOMAIN_NUMBER_OF_LOCAL)
    IF(ALLOCATED(REQUESTS)) DEALLOCATE(REQUESTS)
    IF(ALLOCATED(STATUSES)) DEALLOCATE(STATUSES)
    IF(ALLOCATED(SEND_BUFFERS)) THEN
      DO adjacentDomainIdx=1,SIZE(SEND_BUFFERS)
        IF(ALLOCATED(SEND_BUFFERS(adjacentDomainIdx)%array)) DEALLOCATE(SEND_BUFFERS(adjacentDomainIdx)%array)
      END DO ! adjacentDomainIdx
      DEALLOCATE(SEND_BUFFERS)
    ENDIF
    IF(ALLOCATED(RECEIVE_BUFFERS)) THEN
      DO adjacentDomainIdx=1,SIZE(RECEIVE_BUFFERS)
        IF(ALLOCATED(RECEIVE_BUFFERS(adjacentDomainIdx)%array)) DEALLOCATE(RECEIVE_BUFFERS(adjacentDomainIdx)%array)
      END DO ! adjacentDomainIdx
      DEALLOCATE(RECEIVE_BUFFERS)
    ENDIF
    ERRORSEXITS("DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DOMAIN_MAPPINGS_LOCAL_TO_GLOBAL_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the mapping for a domain mappings mapping and deallocates all memory.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: idx

    ENTERS("DOMAIN_MAPPINGS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(ALLOCATED(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP)) DEALLOCATE(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP)
      IF(ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS)) THEN
        DO idx=1,SIZE(DOMAIN_MAPPING%ADJACENT_DOMAINS,1)
          CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(DOMAIN_MAPPING%ADJACENT_DOMAINS(idx),ERR,ERROR,*999)
        ENDDO !idx
        DEALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS)
      ENDIF
      DEALLOCATE(DOMAIN_MAPPING)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the mapping for a domain mappings mapping.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPING,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to initialise the mappings for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !<The number of domains 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_MAPPING_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(NUMBER_OF_DOMAINS>0) THEN
        DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL=0
        DOMAIN_MAPPING%NUMBER_OF_LOCAL=0
        DOMAIN_MAPPING%NUMBER_OF_GLOBAL=0
        DOMAIN_MAPPING%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
        DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS=0
      ELSE
        LOCAL_ERROR="The specified number of domains of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINS,"*",ERR,ERROR))// &
          & " is invalid. The number of domains must be > 0."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_MAPPING_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_MAPPING_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE DOMAIN_MAPPINGS
