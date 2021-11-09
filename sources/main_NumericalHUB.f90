!***************************************************
!* Book:  How to learn Applied maths
!***************************************************    
program main_NumericalHUB

       
       use API_Example_Systems_of_Equations
       use API_Example_Lagrange_Interpolation
       use API_Example_Cauchy_Problem
       use API_Example_Finite_Differences
       use API_Example_Boundary_Value_Problem
       use API_Example_Initial_Boundary_Value_Problem
       use API_Example_IBVP_and_BVP
       
       use my_examples
       use MUSE_2020
       use API_Example_Fourier_series
       use API_Example_IBVP_Fourier
      
       use API_examples_dependencies
          
       implicit none 
       integer :: option = 1  
     
     
       
       !call myexampleA
       !call myexampleB
       !call myexampleC 
       !call myexampleD
 
           
do while (option>0) 
    
     write(*,*) "Welcome to NumericalHUB" 
     
     write(*,*) " select an option " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Systems of equations  "
     write(*,*) " 2. Lagrange interpolation  " 
     write(*,*) " 3. Finite difference   "
     write(*,*) " 4. ODE Cauchy problems   "
     write(*,*) " 5. Boundary value problems  "
     write(*,*) " 6. Initial-boundary value problems  "
     write(*,*) " 7. Mixed problems: IBVP+BVP  "
     write(*,*) " 8. Advanced methods ODE methods "
     write(*,*) " 9. Advanced methods and problems "
     
     read(*,*) option 
     
     select case(option)
     case(1) 
         call Systems_of_Equations_examples
         
     case(2) 
         call Lagrange_Interpolation_examples 
         
     case(3) 
         call Finite_difference_examples
       
     case(4)
         call Cauchy_problem_examples
      
     case(5) 
         call BVP_examples
      
     case(6) 
         call IBVP_examples 
      
     case(7) 
         call Nonlinear_Plate_Vibrations
       
     case(8) 
         call Advanced_Cauchy_problem_examples
         
     case(9) 
         call Advanced_problems 
         
         case default
              
     end select 
     
end do
  

contains 
    
subroutine Advanced_problems 

integer :: option = 1  
           
do while (option>0) 
    
     write(*,*) "Welcome to Advanced methods" 
     
     write(*,*) " select an option " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Orbits and numerical methods (MUSE 2020)  "
     write(*,*) " 2. SVD applications (not yet implemented)  "
     write(*,*) " 3. Fourier problems  "
     write(*,*) " 4. Chebyshev problems (not yet implemented)  "
     
     read(*,*) option 
     
     select case(option)
     case(1) 
           call Orbits_and_Numerical_Methods
         
     case(2) 
         
         
     case(3) 
           call Fourier_interpolation_examples
           call Fourier_IBVP_examples
         
     case(4)
        
         case default
              
     end select 
     
end do
    
end subroutine 
    
 
end program  

