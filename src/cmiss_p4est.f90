!> \file
!> \author Sander Arens
!> \brief This module is a CMISS buffer module to the p4est library.
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

!> This module is a CMISS buffer module to the p4est library.
MODULE CMISS_P4EST
  
  USE BASE_ROUTINES
  USE KINDS
  USE ISO_VARYING_STRING
  USE ISO_C_BINDING

  IMPLICIT NONE
 
  PRIVATE

  !Module parameters

  !Module types
  TYPE P4estType
    TYPE(C_PTR) :: p4est
  END TYPE P4estType

  TYPE P4estConnectivityType
    TYPE(C_PTR) :: p4estConnectivity
  END TYPE P4estConnectivityType

  !Module variables

  !Interfaces
  INTERFACE
    FUNCTION p4est_connectivity_new(num_vertices,num_trees,num_corners,num_ctt) BIND(C,NAME='p4est_connectivity_new')
      USE ISO_C_BINDING, ONLY:p4est_topidx_t=>C_INT32_T,C_PTR
      TYPE(C_PTR) :: p4est_connectivity_new
      INTEGER(p4est_topidx_t) :: num_vertices
      INTEGER(p4est_topidx_t) :: num_trees
      INTEGER(p4est_topidx_t) :: num_corners
      INTEGER(p4est_topidx_t) :: num_ctt
    END FUNCTION
  END INTERFACE

!  PUBLIC 
  
CONTAINS

  !
  !================================================================================================================================
  !


  !
  !================================================================================================================================
  !
    
END MODULE CMISS_P4EST
