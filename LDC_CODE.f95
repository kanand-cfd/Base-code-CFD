PROGRAM LID_DRIVEN_CAVITY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!																																					   !	 
! THE FOLLOWING CODE AIMS TO SOLVE FOR THE VELOCITY AND PRESSURE FIELD IN A LID DRIVEN CAVITY FLOW. IT HAS BEEN WRITTEN FOR A 21X21 GRID FOR A CREEPING! 
! FLOW CONDITION i.e. REYNOLDS NUMBER IS 1. THE TIME STEP SIZE HAS BEEN TAKEN AS 0.001. THE SOLVER RUNS FOR 1000 TIME STEPS AND SOLVES THE X-MOMENTUM  ! 
! AND Y-MOMENTUM EQUATION ALONG WITH THE PRESSURE POISSON EQUATION TO SOLVE FOR PRESSURE. THE CODE CONTAINS TWO SUBROUTINES, ONE TO SOLVE THE PRESSURE !
! POISSON EQUATION AND THE OTHER TO CALCULATE THE VALUE OF A SOURCE TERM OCCURING THE PRESSURE POISSON EQUATION. BOTH THE SUBROUTINES HAVE BEEN CALLED !
! INSIDE THE DO-WHILE LOOP OF THE MAIN PROGRAM.																										   !	
!																																					   !	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




DOUBLE PRECISION, DIMENSION(100, 100):: P, UN, VN, ST, VELOCITY    ! INITIALIZATION OF SOLUTION MATRICES

DOUBLE PRECISION DX, DY, DT, NU, A , B, C, D, X, Y, T              ! INTIALIZATION OF GRID VARIABLES, TIME STEP SIZE,  
INTEGER RHO, NX, NY, NIT, I, J		 							   ! ITERATING VARIABLES, DENSITY AND VISCOSITY	


NX = 21															   ! SPATIAL DISCRETIZATION WITH 21 NODES 	
NY = 21															   ! ON BOTH X AND Y AXIS	
X = 1
Y= 1
DX = X/(NX-1)                                                      ! UNIFORM GRID SPACINGS OF 0.05
DY = Y/(NY-1)                                                      ! FOR BOTH X AND Y
DT = 0.001
T = 0
RHO = 1         												   ! DENSITY = 1			
NU = 1                                                             ! KINEMATIC VISCOSITY = 1
NIT = 1                                                            ! FOR A TIME STEP OF 0.001 SECONDS AND A TOTAL 
																   ! OF 1000 TIMESTEPS THE CALCULATION RUNS FOR 1 SECOND IN REAL TIME


A = DT/DX 														   ! VARIABLES NEEDED FOR CALCULATION IN 
B = DT/DY														   ! THE MOMENTUM EQUATION EXPRESSION 						
C = DT/(DX**2)
D = DT/(DY**2)


! BOTH THE LOOPS IN THE FOLLOWING LINES ARE FOR INITIALIZING 
! BOUNDARY CONDITIONS AND INITIAL VALUES IN THE SOLUTION MATRICES

DO I = 1, NX													   		
  DO J= 1, NY													   		
    IF (J == NY) THEN											  	
      UN(I, NY) = 1													
      VN(I, NY) = 0
      P(I, NY) = 0
      ST(I, NY) = 0 
     END IF 
  END DO
END DO

DO I = 1, NX
  DO J = 1, NY-1
    UN(I, J) = 0 
      VN(I, J) = 0
      P(I, J) = 0
      ST(I, J) = 0
  END DO
END DO
      
!BEGINNING OF THE SOLVER LOOP 

DO WHILE (T < NIT)
  DO I = 2, NX-1
    DO J = 2, NY-1

      CALL SOURCE_TERM(RHO, DX, DY, DT, UN, VN, I, J, ST)   ! CALLING THE SOURCE TERM SUBROUTINE TO BE USED IN THE PRESSURE POISSON EQUATION
      CALL PRESS_POISSON(DX, DY, ST, I, J, P, NX, NY)         ! CALLING THE PRESSURE POISSON SUBROUTINE TO CALCULATE PRESSURE
	  
      ! FOLLOWING IS THE EXPRESSION TO CALCULATE U_VELOCITY EXPILCITLY USING THE DISCRETIZED X-MOMENTUM EQUATION   

      UN(I,J) = UN(I,J) - (UN(I,J)*A*(UN(I,J) - UN(I-1,J))) - (VN(I,J)*B*(UN(I,J) - UN(I,J-1))) - &             
      & ((A/(2*RHO))*(P(I+1,J) - P(I-1,J))) + (NU*(UN(I+1,J) - (2*UN(I,J)) + &
      &  UN(I-1,J))*C) + (NU*(UN(I,J+1) - (2*UN(I,J)) + UN(I,J-1))*D)

      ! FOLLOWING IS THE EXPRESSION TO CALCULATE V_VELOCITY EXPILCITLY USING THE DISCRETIZED Y-MOMENTUM EQUATION
      
      VN(I,J) = VN(I,J) - (UN(I,J)*A*(VN(I,J) - VN(I-1,J))) - &
      & (VN(I,J)*B*(VN(I,J) - VN(I,J-1))) - ((B/(2*RHO))*(P(I+1,J) - P(I-1,J))) + (NU*(VN(I+1,J) - &
      & (2*VN(I,J)) + VN(I-1,J))*C) + (NU*(VN(I,J+1) - (2*VN(I,J)) + VN(I,J-1))*D)
	  
	  
     END DO
    END DO
  T = T + DT      ! INCREMENTING THE TOTAL TIME BY TIME STEP VALUE 

END DO

! THE RESULTS ARE TO BE WRITTEN IN A FILE 

OPEN(unit = 100, file = 'U_VELOCITY.TXT', status = 'unknown')           !  OPENING A FILE TO STORE X VELOCITY COMPONENTS 
OPEN(unit = 200, file = 'V_VELOCITY.TXT', status = 'unknown')           !  OPENING A FILE TO STORE Y VELOCITY COMPONENTS
OPEN(unit = 300, file = 'PRESSURE.TXT', status = 'unknown')             !  OPENING A FILE TO STORE PRESSURE
OPEN(unit = 400, file = 'VELOCITY.TXT', status = 'unknown')             !  OPENING A FILE TO STORE MAGNITUDE OF VELOCITY 

! LOOP TO WRITE THE VALUES IN THE FILES 
    
DO I = 1, NX
  DO J = 1, NY
   VELOCITY(I,J) = ((UN(I,J))**2 + (VN(I,J))**2)**0.5 
   WRITE(100, *)I, J, UN(I , J)
   WRITE(200, *)I, J, VN(I , J)
   WRITE(300, *)I, J, P(I , J)
   WRITE(400, *)I, J, VELOCITY(I,J)
  END DO
   WRITE(100, *)
   WRITE(200, *)
   WRITE(300, *)
   WRITE(400, *)
END DO

!CLOSING ALL THE OPENED FILES
   
CLOSE(100)
CLOSE(200)
CLOSE(300)
CLOSE(400)


END PROGRAM    ! END OF THE MAIN PROGRAM


! SUBROUTINE TO CALCULATE THE SOURCE TERM REQUIRED IN PRESSURE POISSON EQUATION

SUBROUTINE SOURCE_TERM(RHO, DX, DY, DT, UN, VN, I, J, ST)

DOUBLE PRECISION, DIMENSION(100, 100):: ST, UN, VN
DOUBLE PRECISION DX, DY, DT, CF, C
INTEGER RHO, I, J

CF = (RHO*(DX**2)*(DY**2))/(2*((DX**2) + (DY**2)))
C = 1/DT

ST(I,J) = CF*((C*(((UN(I+1,J) - UN(I-1,J))/(2*DX)) + ((VN(I,J+1) - VN(I,J-1))/(2*DY)))) - &
& (((UN(I+1,J) - UN(I-1, J))/(2*DX))**2) - (((VN(I,J+1) - VN(I,J-1))/(2*DY))**2) - &
& ((UN(I,J+1) - UN(I, J-1))/(2*DY))*((VN(I+1,J) - VN(I-1, J))/(2*DX)))

RETURN

END SUBROUTINE SOURCE_TERM ! END OF THE SOURCE TERM SUBROUTINE

! SUBROUTINE TO CALCULATE THE PRESSURE, EXPLICITLY USING PRESSURE POISSON EQUATION

SUBROUTINE PRESS_POISSON(DX, DY, ST, I, J, P, NX, NY)

DOUBLE PRECISION, DIMENSION(100, 100):: P, ST
DOUBLE PRECISION DX, DY
INTEGER I, J

P(I, J) = ((P(I+1, J) + P(I-1, J))*(DY**2) + (P(I, J+1) + P(I, J-1))*(DX**2))/(2*((DX**2) + (DY**2))) - ST(I, J)

 P(1, J) = P(2, J)
 P(NX, J) = P(NX-1, J)
 P(I, 1) = P(I, 2)
 P(I, NY) = P(I, NY-1)

RETURN

END SUBROUTINE PRESS_POISSON    ! END OF THE PRESSURE POISSON EQUATION SUBROUTINE