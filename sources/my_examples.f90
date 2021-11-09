
module my_examples

        use dislin 
        use Linear_systems
        use Cauchy_Problem
        use Temporal_Schemes
        use plots 
        implicit none 

       real, parameter :: PI = 4 * atan(1d0) 
       
    contains  
  


!**************************************************************************************************
! Plot a simple graph 
!*****************************************************************************************
subroutine myexampleA
 
 integer, parameter :: N=200
 real :: x(0:N), y(0:N)
 integer :: i 
 real :: a = 0, b = 2 * PI 
 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 

   y = sin(x ) 
 
   call scrmod("reverse")
   call qplot(x, y, N+1) 
 
end subroutine


!**************************************************************************************************
! plot a simple trayectory 
!*****************************************************************************************
subroutine myexampleB

        
   
    integer, parameter :: N = 200  !Time steps
    real :: Time(0:N), U(0:N,1)
    real :: t0 = 0, tf = 8
    integer :: i 

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
    U(0,1) =  1
   
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F,& 
                          Solution = U, Scheme = Euler )
    
    call scrmod("reverse")
    call qplot(Time, U, N+1) 
  
contains

function F( U, t ) 
    real :: U(:), t 
    real :: F(size(U)) 
    
    F(1) = - U(1)
  
end function 

end subroutine
 
!**************************************************************************************************
! Plot a simple graph and generates latex file in a specific location 
!*****************************************************************************************
subroutine myexampleC 
 
    integer, parameter :: N=200, Np = 3 
    real :: x(0:N), y(0:N, Np), a = 0, b = 2 * PI 
    integer :: i 
    character(len=100) :: path(4) =                    & 
    ["./results/myexampleCa", "./results/myexampleCb", & 
     "./results/myexampleCc", "./results/myexampleCd"  ]
    
    x  = [ (a + (b-a)*i/N, i=0, N) ] 
    y(:, 1)  = sin(x); y(:, 2)  = cos(x); y(:, 3)  = sin(2*x)

   call plot_parametrics( x, y, ["$\sin x$", "$\cos x$", "$\sin 2x$"], & 
                         "$x$", "$y$", "(a)", path(1) ) 
   call plot_parametrics( y(:,1), y(:,:), ["O1", "O2", "O3"],          & 
                         "$y_2$", "$y_1$", "(b)", path(2) ) 
   call plot_parametrics( y(:,1), y(:,2:2), ["O2"],  "$y_2$", "$y_1$", & 
                          "(c)", path(3) ) 
   call plot_parametrics( y(:,1), y(:,3:3), ["O3"], "$y_2$", "$y_1$",  & 
                          "(d)", path(4) ) 
end subroutine    
    
!**************************************************************************************************
! Plot a contour graph and generates latex file in a specific location 
!*****************************************************************************************
subroutine myexampleD 
 
    integer, parameter :: N=20, Nl = 29 
    real :: x(0:N), y(0:N), z(0:N, 0:N)
    real :: levels(0:Nl), a = 0, b = 2 * PI  
    integer :: i 
    character(len=100) :: path(2) =  ["./results/myexampleDa", & 
                                      "./results/myexampleDb"  ] 
    x  = [ (a + (b-a)*i/N, i=0, N) ] 
    y  = [ (a + (b-a)*i/N, i=0, N) ]
    a = -1; b = 1 
    levels  = [ (a + (b-a)*i/Nl, i=0, Nl) ]
    z = Tensor_product( sin(x), sin(y) ) 
   
    call plot_contour(x, y, z, "x", "y", levels, "(a)",path(1),"color") 
    call plot_contour(x, y, z, "x", "y", levels, "(b)",path(2),"isolines") 
end subroutine       




end module  

