PROGRAM FORTRAN_TEST

!To solve the problem of heat transfer in a plate using RELAXATION TECHNIQUE to solve the elliptic equation

100 REAL, DIMENSION(2000, 2000)::T, T_OLD    !Solution matrix


REAL DX, DY, X, Y, Z, ERR, LX, LY, T1, T2, T3, T4, T0 ! Variables to solve the problem


INTEGER I, J, K, M, N, COUNT, MAX_COUNT, CHOICE, DUMMY !Iteration variable

CALL INITIALIZATION(LX, LY, M, N, DX, DY, X, Y,Z) !SUBROUTINE TO INTIALIZE THE THE GRID VARIABLE

CALL BOUNDARY_CONDITIONS(T1, T2, T3, T4, T0) !SUBROUTINE FOR BOUNDARY CONDITIONS

CALL BC_INIT(I, J, M, N, T1, T2, T3, T4, T0, T) !SUBROUTINE FOR INITIALIZING BOUNDARY CONDITIONS IN THE SOLUTION MATRIX

COUNT = 0                     !INITIALIZING VARIABLES FOR
MAX_COUNT = (M-2)*(N-2)		  !ENSURING CONVERGENCE 		


CALL SOLVER(I, J, K, M, N, ERR, COUNT, MAX_COUNT, T, X, Y, Z, DUMMY, T_OLD) !SUBROUTINE FOR RUNNING THE SOLVER UPTO CONVERGENCE


PRINT*, "DO YOU WANT TO RUN THE PROBLEM AGAIN(1/0)"
READ*, CHOICE

!ASKING THE USER IF THEY WANT TO RUN THE PROBLEM AGAIN

IF (CHOICE == 1) THEN
  GOTO 100
END IF   
END PROGRAM


!SUBROUTINE DEFINITION TO INTIALIZE THE THE GRID VARIABLES


SUBROUTINE INITIALIZATION(LX, LY, M, N, DX, DY , X, Y, Z)

INTEGER M, N
REAL LX, LY, DX, DY

PRINT*,"ENTER THE LENGTH AND BREADTH OF THE PLATE, RESPECTIVELY:" !Asking the user to enter plate dimensions
READ*, LX, LY

!Now we ask the user to enter the number of divisions the want on both length and breadth

PRINT*,"ENTER THE NUMBER OF DIVISIONS YOU WANT ON THE LENGTH SIDE OF THE PLATE (MAX 2000):" 
READ*, M

PRINT*,"ENTER THE THE NUMBER OF DIVISIONS YOU WANT ON THE BREADTH SIDE OF THE PLATE (MAX 2000):"
READ*, N


DX = LX/M !Calculating the value of smallest gird size in x direction

DY = LY/N !Calculating the value of smallest grid size in y direction


!We calculate the values required for calculating the iteration formula

X = DX**2 

Y = DY**2

Z = (X*Y)/(2*(X + Y))

RETURN

END SUBROUTINE


!SUBROUTINE DEFINITION FOR BOUNDARY CONDITIONS FROM USER 


SUBROUTINE BOUNDARY_CONDITIONS(T1, T2, T3, T4, T0)

REAL T1, T2, T3, T4, T0

!Asking the users to enter the boundary conditions

PRINT*,"ENTER THE TEMPARATURE AT THE LEFT BOUNDARY:"
READ*, T1

PRINT*,"ENTER THE TEMPARATURE AT THE RIGHT BOUNDARY:"
READ*, T2

PRINT*,"ENTER THE TEMPARATURE AT THE BOTTOM BOUNDARY:"
READ*, T3

PRINT*,"ENTER THE TEMPARATURE AT THE TOP BOUNDARY:"
READ*, T4


T0 = (T1 + T2 + T3 + T4)*0.25 !Calculating the value of initial guess using Boundary Conditions

RETURN 

END SUBROUTINE


!SUBROUTINE DEFINITION TO INTIALIZE THE BOUNDARY CONDTIONS IN THE SOLUTION MATRIX

SUBROUTINE BC_INIT(I, J, M, N, T1, T2, T3, T4, T0, T)

REAL, DIMENSION(2000, 2000)::T
REAL T1, T2, T3, T4, T0
INTEGER I, J, M, N

!Now we energize the solution matrix(i.e. the grid points) with initial value
!and boundary conditions and write this data in a file

OPEN( unit = 100, file = 'LAPLACE_INITIAL.dat', status = 'unknown') !DEFINITION FOR OPENING THE FILE NAMED LAPLACE_INITIAL

 DO I = 1, M
  DO J = 1, N
    IF(I == 1) THEN
      T(1, J) = T3     
    ELSE IF(I == M) THEN
      T(M, J) = T4
    ELSE IF(J == 1) THEN   
      T(I, 1) = T1      
    ELSE IF(J == N) THEN
      T(I, N) = T2
    ELSE
      T(I, J) = T0
    END IF  
  END DO
 END DO

!Now we write the above matrix in the file that we opened above
    
DO I = 1, M
  DO J = 1, N
   WRITE(100, *)I,J,T(I , J)
  END DO
   WRITE(100, *)
END DO   
CLOSE(100)

RETURN

END SUBROUTINE


!SUBROUTINE DEFINITION TO RUN THE SOLVER UNTIL CONVERGENCE


SUBROUTINE SOLVER(I, J, K, M, N, ERR, COUNT, MAX_COUNT, T, X, Y, Z, DUMMY, T_OLD)

INTEGER I, J, K, M, N, COUNT, MAX_COUNT, DUMMY
REAL X, Y, Z, ERR
REAL, DIMENSION(2000, 2000)::T, T_OLD

!Now we solve the problem by running the  solver 

K = 0

OPEN( unit = 222, file = 'LAPLACE_SOLUTION.dat', status = 'unknown') !DEFINITION FOR OPENING THE FILE NAMED LAPLACE_SOLUTION

DO WHILE(K < 200000)
  DO I = 2, M-1
    DO J = 2, N-1
       T_OLD(I, J) = T(I, J)
       T(I, J) = Z*((T(I-1, J) + T(I+1, J))/X + (T(I, J-1) + T(I, J+1))/Y)
       
	   !The variable ERR calculates the difference between T(i,j) at present iterations
       !and previous iteration, if the value of ERR is less than 10^-3 then 1 is added
       !to COUNT, this means that the present node has converged to an acceptable value 
       
       ERR = T(I, J) - T_OLD(I, J)
      
       IF (ABS(ERR) < 1E-6) THEN
         COUNT = COUNT + 1
       END IF  
    END DO
  END DO
  
  !If the value of COUNT is equal to the total number of nodes in the
  !calculation domain that means the solution has converged for the 
  !whole domain and now the solver can stop, if not the count is equal to
  !zero, and the solver moves to next iteration
  DUMMY = COUNT
  IF (COUNT == MAX_COUNT ) THEN
    
    GOTO 100
  ELSE
    COUNT = 0  
  END IF
  K = K + 1     
END DO


!Now we print at which iteration the solution has stopped and write the solution
!in a file

100 PRINT*, "THE SOLUTION HAS CONVERGED AT ALL THE NODES AFTER", K," ITERATION"

DO I = 1, M
 DO J = 1, N
   WRITE(222, *)I,J,T(I , J)
 END DO
 WRITE(222, *)
END DO   

PRINT*," "

CLOSE(222)

RETURN

END SUBROUTINE