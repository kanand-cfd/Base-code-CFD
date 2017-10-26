PROGRAM ELLIPTIC

!To solve the problem of heat transfer in a plate using RELAXATION TECHNIQUE to solve the elliptic equation

REAL, DIMENSION(500, 200, 200)::T    !Solution matrix

REAL DX, DY, X, Y, Z, ERR, LX, LY, T1, T2, T3, T4, T0 ! Variables to solve the problem

INTEGER I, J, K, M, N, A, COUNT, MAX_COUNT !Iteration variable

PRINT*,"ENTER THE LENGTH AND BREADTH OF THE PLATE, RESPECTIVELY:" !Asking the user to enter plate dimensions
READ*, LX, LY

!Now we ask the user to enter the number of divisions the want on both length and breadth

PRINT*,"ENTER THE NUMBER OF DIVISIONS YOU WANT ON THE LENGTH SIDE OF THE PLATE:" 
READ*, M

PRINT*,"ENTER THE THE NUMBER OF DIVISIONS YOU WANT ON THE BREADTH SIDE OF THE PLATE:"
READ*, N

!Asking the users to enter the boundary conditions

PRINT*,"ENTER THE TEMPARATURE AT THE LEFT BOUNDARY:"
READ*, T1

PRINT*,"ENTER THE TEMPARATURE AT THE RIGHT BOUNDARY:"
READ*, T2

PRINT*,"ENTER THE TEMPARATURE AT THE BOTTOM BOUNDARY:"
READ*, T3

PRINT*,"ENTER THE TEMPARATURE AT THE TOP BOUNDARY:"
READ*, T4

!We ask the user to enter the maximum number of iterations the solver should run

PRINT*,"ENTER THE MAXIMUM NUMBER OF ITERATION THE SOLVER SHOULD RUN:"
READ*, A

DX = LX/M !Calculating the value of smallest gird size in x direction

DY = LY/N !Calculating the value of smallest grid size in y direction

T0 = (T1 + T2 + T3 + T4)*0.25 !Calculating the value of initial guess using Boundary Conditions

!We calculate the values required for calculating the iteration formula

X = DX**2 

Y = DY**2

Z = (X*Y)/(2*(X + Y))

!Now we initialize the variables required for deciding the convergence of the problem

COUNT = 0

MAX_COUNT = (M-2)*(N-2)

!Now we energize the solution matrix(i.e. the grid points) with initial value
!and boundary conditions and write this data in a file

OPEN( unit = 100, file = 'ELLIPTIC_INITIAL.dat', status = 'unknown')

DO K = 1, A
 DO I = 1, M
  DO J = 1, N
    IF(I == 1) THEN
      T(K, 1, J) = T3     
    ELSE IF(I == M) THEN
      T(K, M, J) = T4
    ELSE IF(J == 1) THEN   
      T(K, I, 1) = T1      
    ELSE IF(J == N) THEN
      T(K, I, N) = T2
    ELSE
      T(K, I, J) = T0
    END IF  
  END DO
 END DO
END DO

    
DO I = 1, M
  DO J = 1, N
   WRITE(100, *)I,J,T(1, I , J)
  END DO
   WRITE(100, *)
END DO   
CLOSE(100)

!Now we solve the problem by running the  solver 

OPEN( unit = 222, file = 'ELLIPTIC_SOLUTION.dat', status = 'unknown')

DO K = 2, A
  DO I = 2, M-1
    DO J = 2, N-1
       T(K, I, J) = Z*((T(K, I-1, J) + T(K-1, I+1, J))/X + (T(K, I, J-1) + T(K-1, I, J+1))/Y)
       
	   !The variable ERR calculates the difference between T(i,j) at present iterations
       !and previous iteration, if the value of ERR is less than 10^-3 then 1 is added
       !to COUNT, this means that the present node has converged to an acceptable value 
       
       ERR = T(K, I, J) - T(K-1, I, J)
      
       IF (ABS(ERR) < 1E-3) THEN
         COUNT = COUNT + 1
       END IF  
    END DO
  END DO
  
  !If the value of COUNT is equal to the total number of nodes in the
  !calculation domain that means the solution has converged for the 
  !whole domain and now the solver can stop, if not the count is equal to
  !zero, and the solver moves to next iteration
  
  IF (COUNT == MAX_COUNT ) THEN
    GOTO 100
  ELSE
    COUNT = 0  
  END IF     
END DO
GOTO 200

!Now we print at which iteration the solution has stopped and write the solution
!in a file

100 PRINT*, "THE SOLUTION HAS CONVERGED AT ALL THE NODES AFTER", K,"ITERATION"

DO I = 1, M
 DO J = 1, N
   WRITE(222, *)I,J,T(K, I , J)
 END DO
 WRITE(222, *)
END DO   
GOTO 300

200 PRINT*, "THE SOLUTION HAS CONVERGED AFTER", K-1,"ITERATION"

DO I = 1, M
 DO J = 1, N
   WRITE(222, *)I,J,T(K-1, I , J)
 END DO
 WRITE(222, *)
END DO   

300 PRINT*," "

CLOSE(222)


END PROGRAM