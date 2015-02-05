!> \file
!> \author Chris Bradley
!> \brief This module handles all constants shared across constraint condition routines.
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

!> This module defines all constants shared across constraint condition routines.
MODULE CONSTRAINT_CONDITIONS_CONSTANTS

  USE KINDS

  IMPLICIT NONE

  !Module parameters

  !> \addtogroup CONSTRAINT_CONDITIONS_CONSTANTS_Methods CONSTRAINT_CONDITIONS_CONSTANTS::Methods
  !> \brief Constraint condition methods.
  !> \see CONSTRAINT_CONDITIONS_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_LAGRANGE_MULTIPLIERS_METHOD=1 !<Lagrange multipliers constraint condition method. \see CONSTRAINT_CONDITIONS_CONSTANTS_Methods,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_AUGMENTED_LAGRANGE_METHOD=2 !<Augmented Lagrange multiplers constraint condition method. \see CONSTRAINT_CONDITIONS_CONSTANTS_Methods,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_PENALTY_METHOD=3 !<Penalty constraint condition method. \see CONSTRAINT_CONDITIONS_CONSTANTS_Methods,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_POINT_TO_POINT_METHOD=4 !<Point to point constraint condition method. \see CONSTRAINT_CONDITIONS_CONSTANTS_Methods,CONSTRAINT_CONDITIONS_CONSTANTS
  !>@}

  !> \addtogroup CONSTRAINT_CONDITIONS_Operators CONSTRAINT_CONDITIONS_CONSTANTS::Operators
  !> \brief Constraint condition operators.
  !> \see CONSTRAINT_CONDITIONS_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_FE_INCOMPRESSIBILITY_OPERATOR=1 !<Finite elasticity incompressibility operator, i.e., lambda.(J-1). \see CONSTRAINT_CONDITIONS_CONSTANTS_Operators,CONSTRAINT_CONDITIONS_CONSTANTS 
  !>@}

  !> \addtogroup CONSTRAINT_CONDITIONS_CONSTANTS_LinearityTypes CONSTRAINT_CONDITIONS_CONSTANTS::LinearityTypes
  !> \brief The constraint condition linearity types
  !> \see CONSTRAINT_CONDITIONS_CONSTANTS,OPENCMISS_EquationsLinearityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_CONSTRAINT_CONDITION_LINEARITY_TYPES=3 !<The number of constraint conditions linearity types defined. \see CONSTRAINT_CONDITIONS_CONSTANTS_LinearityTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_LINEAR=1 !<The constraint conditions are linear. \see CONSTRAINT_CONDITIONS_CONSTANTS_LinearityTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_NONLINEAR=2 !<The constraint conditions are non-linear. \see CONSTRAINT_CONDITIONS_CONSTANTS_LinearityTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_NONLINEAR_BCS=3 !<The constraint conditions have non-linear boundary conditions. \see CONSTRAINT_CONDITIONS_CONSTANTS_LinearityTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  !>@}

  !> \addtogroup CONSTRAINT_CONDITIONS_CONSTANTS_TimeDepedenceTypes CONSTRAINT_CONDITIONS_CONSTANTS::TimeDepedenceTypes
  !> \brief The constraint condition time dependence type parameters
  !> \see CONSTRAINT_CONDITIONS_CONSTANTS,OPENCMISS_EquationsTimeDependenceTypes
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_CONSTRAINT_CONDITION_TIME_TYPES=5 !<The number of constraint conditions time dependence types defined. \see CONSTRAINT_CONDITIONS_CONSTANTS_TimeDependenceTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_STATIC=1 !<The constraint conditions are static and have no time dependence. \see CONSTRAINT_CONDITIONS_CONSTANTS_TimeDependenceTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_QUASISTATIC=2 !<The constraint conditions are quasi-static. \see CONSTRAINT_CONDITIONS_CONSTANTS_TimeDependenceTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_FIRST_ORDER_DYNAMIC=3 !<The constraint conditions are first order dynamic. \see CONSTRAINT_CONDITIONS_CONSTANTS_TimeDependenceTypes,CONSTRAINT_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_SECOND_ORDER_DYNAMIC=4 !<The constraint conditions are a second order dynamic. \see CONSTRAINT_CONDITIONS_CONSTANTS_TimeDependenceTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_TIME_STEPPING=5 !<The constraint conditions are for time stepping. \see CONSTRAINT_CONDITIONS_CONSTANTS_TimeDependenceTypes,EQUATIONS_ROUTINES
  !>@}
  
  !> \addtogroup CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods CONSTRAINT_CONDITION_CONSTANTS::SolutionMethods
  !> \brief The solution method parameters
  !> \see CONSTRAINT_CONDITION_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_CONSTRAINT_CONDITION_SOLUTION_METHODS=7 !<The number of solution methods defined. \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_FEM_SOLUTION_METHOD=1 !<Finite Element Method solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_BEM_SOLUTION_METHOD=2 !<Boundary Element Method solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_FD_SOLUTION_METHOD=3 !<Finite Difference solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_FV_SOLUTION_METHOD=4 !<Finite Volume solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_GFEM_SOLUTION_METHOD=5 !<Grid-based Finite Element Method solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_GFD_SOLUTION_METHOD=6 !<Grid-based Finite Difference Method solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  INTEGER(INTG), PARAMETER :: CONSTRAINT_CONDITION_GFV_SOLUTION_METHOD=7 !<Grid-based Finite Volume solution method \see CONSTRAINT_CONDITION_CONSTANTS_SolutionMethods,CONSTRAINT_CONDITION_CONSTANTS
  !>@}

END MODULE CONSTRAINT_CONDITIONS_CONSTANTS
