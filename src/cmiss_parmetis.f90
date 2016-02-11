!> \file
!> \author Chris Bradley
!> \brief This module is a CMISS buffer module to the ParMETIS library.
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

!> This module is a CMISS buffer module to the ParMETIS library.
MODULE CMISS_PARMETIS
  
  USE BASE_ROUTINES
  USE KINDS
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING

#define idx_t C_INT
#define real_t C_DOUBLE
#include "macros.h"

  IMPLICIT NONE
 
  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  INTERFACE

    FUNCTION ParMETIS_V3_Mesh2Dual(elmdist,eptr,eind,numflag,ncommonnodes,xadj,adjncy,comm) BIND(C) 
      USE ISO_C_BINDING, ONLY: C_INT,C_PTR,idx_t
      INTEGER(C_INT) :: ParMETIS_V3_Mesh2Dual
      INTEGER(idx_t) :: elmdist(*)
      INTEGER(idx_t) :: eptr(*)
      INTEGER(idx_t) :: eind(*)
      INTEGER(idx_t) :: numflag
      INTEGER(idx_t) :: ncommonnodes
      TYPE(C_PTR) :: xadj
      TYPE(C_PTR) :: adjncy
      INTEGER(C_INT) :: comm
    END FUNCTION ParMETIS_V3_Mesh2Dual

    FUNCTION ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,vwgt,adjwgt,wgtflag,numflag,ncon,nparts,tpwgts,ubvec,options,edgecut, &
        & part,comm) BIND(C)
      USE ISO_C_BINDING, ONLY: C_INT,idx_t,real_t
      INTEGER(C_INT) :: ParMETIS_V3_PartKway
      INTEGER(idx_t) :: vtxdist(*)
      INTEGER(idx_t) :: xadj(*)
      INTEGER(idx_t) :: adjncy(*)
      INTEGER(idx_t) :: vwgt(*)
      INTEGER(idx_t) :: adjwgt(*)
      INTEGER(idx_t) :: wgtflag
      INTEGER(idx_t) :: numflag
      INTEGER(idx_t) :: ncon
      INTEGER(idx_t) :: nparts
      REAL(real_t) :: tpwgts(*)
      REAL(real_t) :: ubvec(*)
      INTEGER(idx_t) :: options(*)
      INTEGER(idx_t) :: edgecut
      INTEGER(idx_t) :: part(*)
      INTEGER(C_INT) :: comm
    END FUNCTION ParMETIS_V3_PartKway

    FUNCTION ParMETIS_V3_PartMeshKway(elmdist,eptr,eind,elmwgt,wgtflag,numflag,ncon,ncommonnodes,nparts,tpwgts, &
      & ubvec,options,edgecut,part,comm) BIND(C)
      USE ISO_C_BINDING, ONLY: C_INT,idx_t,real_t
      INTEGER(C_INT) :: ParMETIS_V3_PartMeshKway
      INTEGER(idx_t) :: elmdist(*)
      INTEGER(idx_t) :: eptr(*)
      INTEGER(idx_t) :: eind(*)
      INTEGER(idx_t) :: elmwgt(*)
      INTEGER(idx_t) :: wgtflag
      INTEGER(idx_t) :: numflag
      INTEGER(idx_t) :: ncon
      INTEGER(idx_t) :: ncommonnodes
      INTEGER(idx_t) :: nparts
      REAL(real_t) :: tpwgts(*)
      REAL(real_t) :: ubvec(*)
      INTEGER(idx_t) :: options(*)
      INTEGER(idx_t) :: edgecut
      INTEGER(idx_t) :: part(*)
      INTEGER(C_INT) :: comm
    END FUNCTION ParMETIS_V3_PartMeshKway
    
  END INTERFACE

  PUBLIC ParMETIS_Mesh2Dual,PARMETIS_PARTKWAY,PARMETIS_PARTMESHKWAY
  
CONTAINS

  !>Buffer routine to the ParMetis ParMETIS_V3_Mesh2Dual routine.
  SUBROUTINE ParMETIS_Mesh2Dual(elmdist,eptr,eind,numflag,ncommonnodes,xadj,adjncy,comm,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: elmdist(:)
    INTEGER(INTG), INTENT(IN) :: eptr(:)
    INTEGER(INTG), INTENT(IN) :: eind(:)
    INTEGER(INTG), INTENT(IN) :: numflag
    INTEGER(INTG), INTENT(IN) :: ncommonnodes
    INTEGER(INTG), POINTER :: xadj(:)
    INTEGER(INTG), POINTER :: adjncy(:)
    INTEGER(INTG), INTENT(IN) :: comm      
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG), POINTER :: xadjFPtr(:),adjncyFPtr(:)
    TYPE(C_PTR) :: xadjCPtr,adjncyCPtr

    ENTERS("ParMETIS_Mesh2Dual",err,error,*999)

    IF(ASSOCIATED(xadj)) CALL FlagError("The xadj pointer is already associated.",err,error,*999)
    IF(ASSOCIATED(adjncy)) CALL FlagError("The adjncy pointer is already associated.",err,error,*999)

    IF(ParMETIS_V3_Mesh2Dual(elmdist,eptr,eind,numflag,ncommonnodes,xadjCPtr,adjncyCPtr,comm)/=1) THEN
      CALL FlagError("ParMetis error in ParMETIS_V3_Mesh2Dual",err,error,*999)
    ENDIF

    CALL C_F_POINTER(xadjCPtr,xadjFptr,[SIZE(eptr)])
    CALL C_F_POINTER(adjncyCPtr,adjncyFptr,[xadjFptr(SIZE(xadjFptr))-1])
    
    SELECT CASE(numflag)
    CASE(0)
      !C-style numbering
      xadj(0:SIZE(xadjFptr)-1)=>xadjFPtr
      adjncy(0:SIZE(adjncyFptr)-1)=>adjncyFPtr
    CASE(1)
      !F-style numbering
      xadj=>xadjFPtr
      adjncy=>adjncyFPtr
    CASE DEFAULT
      CALL FlagError("Invalid number flag. Must be either 0 or 1.",err,error,*999)
    END SELECT
    
    EXITS("ParMETIS_Mesh2Dual")
    RETURN
999 ERRORSEXITS("ParMETIS_Mesh2Dual",err,error)
    RETURN 1
  END SUBROUTINE ParMETIS_Mesh2Dual

  !
  !================================================================================================================================
  !

  !>Buffer routine to the ParMetis ParMETIS_V3_PartKway routine.
  SUBROUTINE PARMETIS_PARTKWAY(VERTEX_DISTANCE,XADJ,ADJNCY,VERTEX_WEIGHT,ADJ_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
    & NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR,ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: VERTEX_DISTANCE(:)
    INTEGER(INTG), INTENT(IN) :: XADJ(:)
    INTEGER(INTG), INTENT(IN) :: ADJNCY(:)
    INTEGER(INTG), INTENT(IN) :: VERTEX_WEIGHT(:)
    INTEGER(INTG), INTENT(IN) :: ADJ_WEIGHT(:)
    INTEGER(INTG), INTENT(IN) :: WEIGHT_FLAG
    INTEGER(INTG), INTENT(IN) :: NUM_FLAG
    INTEGER(INTG), INTENT(IN) :: NCON
    INTEGER(INTG), INTENT(IN) :: NUMBER_PARTS
    REAL(DP), INTENT(IN) :: TP_WEIGHTS(:)
    REAL(DP), INTENT(IN) :: UB_VEC(:)
    INTEGER(INTG), INTENT(IN) :: OPTIONS(:)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_EDGES_CUT
    INTEGER(INTG), INTENT(OUT) :: PARTITION(:)
    INTEGER(INTG), INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    ENTERS("PARMETIS_PARTKWAY",ERR,ERROR,*999)

    IF(ParMETIS_V3_PartKway(VERTEX_DISTANCE,XADJ,ADJNCY,VERTEX_WEIGHT,ADJ_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
      & NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR)/=1) THEN
      CALL FlagError("ParMetis error in ParMETIS_V3_PartKway",ERR,ERROR,*999)
    ENDIF
    
    EXITS("PARMETIS_PARTKWAY")
    RETURN
999 ERRORSEXITS("PARMETIS_PARTKWAY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE PARMETIS_PARTKWAY

  !
  !================================================================================================================================
  !

  !>Buffer routine to the ParMetis ParMETIS_V3_PartMeshKway routine.
  SUBROUTINE PARMETIS_PARTMESHKWAY(ELEMENT_DISTANCE,ELEMENT_PTR,ELEMENT_INDEX,ELEMENT_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
    & NUMBER_COMMON_NODES,NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR,ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: ELEMENT_DISTANCE(:)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_PTR(:)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_INDEX(:)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_WEIGHT(:)
    INTEGER(INTG), INTENT(IN) :: WEIGHT_FLAG
    INTEGER(INTG), INTENT(IN) :: NUM_FLAG
    INTEGER(INTG), INTENT(IN) :: NCON
    INTEGER(INTG), INTENT(IN) :: NUMBER_COMMON_NODES
    INTEGER(INTG), INTENT(IN) :: NUMBER_PARTS
    REAL(DP), INTENT(IN) :: TP_WEIGHTS(:)
    REAL(DP), INTENT(IN) :: UB_VEC(:)
    INTEGER(INTG), INTENT(IN) :: OPTIONS(:)
    INTEGER(INTG), INTENT(OUT) :: NUMBER_EDGES_CUT
    INTEGER(INTG), INTENT(OUT) :: PARTITION(:)
    INTEGER(INTG), INTENT(IN) :: COMMUNICATOR
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    ENTERS("PARMETIS_PARTMESHKWAY",ERR,ERROR,*999)
    
    IF(ParMETIS_V3_PartMeshKway(ELEMENT_DISTANCE,ELEMENT_PTR,ELEMENT_INDEX,ELEMENT_WEIGHT,WEIGHT_FLAG,NUM_FLAG,NCON, &
      & NUMBER_COMMON_NODES,NUMBER_PARTS,TP_WEIGHTS,UB_VEC,OPTIONS,NUMBER_EDGES_CUT,PARTITION,COMMUNICATOR)/=1) THEN
      CALL FlagError("ParMetis error in ParMETIS_V3_PartMeshKway",ERR,ERROR,*999)
    ENDIF
    
    EXITS("PARMETIS_PARTMESHKWAY")
    RETURN
999 ERRORSEXITS("PARMETIS_PARTMESHKWAY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE PARMETIS_PARTMESHKWAY

  !
  !================================================================================================================================
  !
    
END MODULE CMISS_PARMETIS
