!> \file
!> \authors Sander Arens 
!> \brief This module handles all constraints routines.
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

!> This module handles all constraints class routines.
MODULE CONSTRAINTS_ROUTINES

  USE BASE_ROUTINES
  USE INCOMPRESSIBILITY_EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  
  PUBLIC CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_GET,CONSTRAINTS_PROBLEM_CLASS_TYPE_GET

  PUBLIC CONSTRAINTS_FINITE_ELEMENT_JACOBIAN_EVALUATE,CONSTRAINTS_FINITE_ELEMENT_RESIDUAL_EVALUATE
  
  PUBLIC CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_SET,CONSTRAINTS_FINITE_ELEMENT_CALCULATE, &
    & CONSTRAINTS_EQUATIONS_SET_SETUP,CONSTRAINTS_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & CONSTRAINTS_PROBLEM_CLASS_TYPE_SET,CONSTRAINTS_PROBLEM_SETUP, &
    & CONSTRAINTS_POST_SOLVE,CONSTRAINTS_PRE_SOLVE,CONSTRAINTS_CONTROL_LOOP_PRE_LOOP, &
    & CONSTRAINTS_CONTROL_LOOP_POST_LOOP
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a constraints equation set class.
  SUBROUTINE CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_GET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_TYPE !<On return, the equation type
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SUBTYPE !<On return, the equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_CONSTRAINTS_CLASS) THEN
        EQUATIONS_TYPE=EQUATIONS_SET%TYPE
        EQUATIONS_SUBTYPE=EQUATIONS_SET%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Equations set is not the constraints type",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_GET")
    RETURN
999 CALL ERRORS("CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_GET",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_GET")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a constraints equation set class.
  SUBROUTINE CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_TYPE,EQUATIONS_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_TYPE !<The equation type
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SUBTYPE !<The equation subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_EQUATIONS_SET_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_TYPE)
      CASE(EQUATIONS_SET_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SUBTYPE,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_EQUATIONS_SET_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a constraints class finite element equation set.
  SUBROUTINE CONSTRAINTS_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_INCOMPRESSIBILITY_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("CONSTRAINTS_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a constraints class finite element equation set.
  SUBROUTINE CONSTRAINTS_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("CONSTRAINTS_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vectors for the given element number for a constraints class finite element equation set.
  SUBROUTINE CONSTRAINTS_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("CONSTRAINTS_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a constraints equations set class.
  SUBROUTINE CONSTRAINTS_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("CONSTRAINTS_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_EQUATIONS_SET_SETUP
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a constraints equation set class.
  SUBROUTINE CONSTRAINTS_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_EQUATIONS_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("CONSTRAINTS_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Gets the problem type and subtype for a constraints problem class.
  SUBROUTINE CONSTRAINTS_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<On return, the problem type
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<On return, the proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("CONSTRAINTS_PROBLEM_CLASS_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%CLASS==PROBLEM_CONSTRAINTS_CLASS) THEN
        PROBLEM_EQUATION_TYPE=PROBLEM%TYPE
        PROBLEM_SUBTYPE=PROBLEM%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Problem is not constraints class",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_PROBLEM_CLASS_TYPE_GET")
    RETURN
999 CALL ERRORS("CONSTRAINTS_PROBLEM_CLASS_TYPE_GET",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_PROBLEM_CLASS_TYPE_GET")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_PROBLEM_CLASS_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a constraints problem class.
  SUBROUTINE CONSTRAINTS_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE !<The problem type
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The proboem subtype
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_PROBLEM_CLASS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_EQUATION_TYPE)
      CASE(PROBLEM_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem equation type "//TRIM(NUMBER_TO_VSTRING(PROBLEM_EQUATION_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_PROBLEM_CLASS_TYPE_SET")
    RETURN
999 CALL ERRORS("CONSTRAINTS_PROBLEM_CLASS_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_PROBLEM_CLASS_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_PROBLEM_CLASS_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a constraints problem class.
  SUBROUTINE CONSTRAINTS_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%TYPE)
      CASE(PROBLEM_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("CONSTRAINTS_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a constraints problem class.
  SUBROUTINE CONSTRAINTS_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%PROBLEM%TYPE)
      CASE(PROBLEM_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_POST_SOLVE")
    RETURN
999 CALL ERRORS("CONSTRAINTS_POST_SOLVE",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_POST_SOLVE")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a constraints problem class.
  SUBROUTINE CONSTRAINTS_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONSTRAINTS_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%PROBLEM%TYPE)
      CASE(PROBLEM_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONSTRAINTS_PRE_SOLVE")
    RETURN
999 CALL ERRORS("CONSTRAINTS_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE CONSTRAINTS_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINTS_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%PROBLEM%TYPE)
      CASE(PROBLEM_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSIBILITY_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINTS_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("CONSTRAINTS_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE CONSTRAINTS_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("CONSTRAINTS_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%PROBLEM%TYPE)
      CASE(PROBLEM_INCOMPRESSIBILITY_TYPE)
        CALL INCOMPRESSBILITY_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a constraints problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CONSTRAINTS_CONTROL_LOOP_POST_LOOP")
    RETURN
999 CALL ERRORS("CONSTRAINTS_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    CALL EXITS("CONSTRAINTS_CONTROL_LOOP_POST_LOOP")
    RETURN 1
  END SUBROUTINE CONSTRAINTS_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

END MODULE CONSTRAINTS_ROUTINES

