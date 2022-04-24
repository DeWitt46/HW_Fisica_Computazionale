PROGRAM CLT
IMPLICIT NONE
!-----------------------------------------------------
!Code to numerically prove central limit theory 
!starting from uniform and exponential distribution 
!-----------------------------------------------------

INTEGER, PARAMETER  :: nkind = 4                                               !4 for single precision, 8 for double

REAL(kind=nkind), PARAMETER     :: pi = 3.14159265359
INTEGER, PARAMETER  :: n_num = 8                                               !number of different N

INTEGER(kind=nkind)                             :: i, j, k                         !Counters
INTEGER(kind=nkind), DIMENSION(:), ALLOCATABLE  :: N, nrun                         !Number of steps and points
INTEGER(kind=nkind)                             :: nbin                            !Number of bins
REAL(kind=nkind), DIMENSION(:,:), ALLOCATABLE   :: P_N, H_u, H_exp, H_L            !Histograms from gaussian, uniform, 
REAL(kind=nkind), DIMENSION(:,:), ALLOCATABLE   :: L1, L2                          !exponential distribution

REAL(kind=nkind), DIMENSION(:), ALLOCATABLE     :: x, y, z                         !Average variables
REAL(kind=nkind)                                :: dr                              !Bin length
REAL(kind=nkind)                                :: mu_n, sigma_n                   !Numerical gaussian parameters for x_N

REAL(kind=nkind)                                :: zero = 0., uno = 1.             !Normal gaussian parameters with
                                                                                   !the right kind

REAL(kind=nkind), PARAMETER                     :: mu_a_u    = float(0)            !Analytical gaussian parameters for x_N 
REAL(kind=nkind), PARAMETER                     :: s_u       = sqrt(float(1)/3)    !UNIFORM DISTRIBUTION IN [-1;1]
REAL(kind=nkind)                                :: sigma_a_u = sqrt(float(1)/3)    !UNIFORM DISTRIBUTION IN [-1;1]

REAL(kind=nkind), PARAMETER                     :: mu_a_e    = float(1)            !Analytical gaussian parameters for x_N
REAL(kind=nkind), PARAMETER                     :: s_e       = sqrt(float(1))      !EXPONENTIAL DISTRIBUTION
REAL(kind=nkind)                                :: sigma_a_e = sqrt(float(1))      !EXPONENTIAL DISTRIBUTION

REAL(kind=nkind), PARAMETER                     :: mu_a_L    = float(0)            !Analytical gaussian parameters for x_N
REAL(kind=nkind), PARAMETER                     :: s_L       = sqrt(2*(2 - atan(2.))/pi)    !LORENTZ DISTRIBUTION      
REAL(kind=nkind)                                :: sigma_a_L = 0.                  !LORENTZ DISTRIBUTION      

REAL(kind=nkind), DIMENSION(:), ALLOCATABLE     :: rnd                             !Array random numbers




INTEGER                            :: sizer, values(1:8)                       !Seed variables
INTEGER, DIMENSION(:), ALLOCATABLE :: seed, box                                !Seed variables






!--------------------------------------------
!Seed generation and saving
!--------------------------------------------

CALL random_seed(sizer)
ALLOCATE(seed(sizer))
ALLOCATE(box(sizer))       
CALL date_and_time(values=values)
seed(:) = values(8)

PRINT*,'Here the seed has ',sizer,'components'  
CALL random_seed(put=seed)
CALL random_seed(get=box)
OPEN (unit=1,file="seed.dat",status="replace",action="write")
DO i = 1 , sizer
   WRITE(unit=1,fmt=*)box(i)
END DO
CLOSE(1)
PRINT*,'seed data stored in "seed.dat"'





!-----------------------------------------------
!Read and set initial value
!-----------------------------------------------
PRINT*,'Numerical proof for central limit theorem'
PRINT*,""
PRINT*, "Enter number of steps for average>"
!READ*, N
PRINT*, "Enter number of x_N needed>"
!READ*, nrun
PRINT*, "Enter number of bin needed>"
!READ*, nbin

nbin = 350

ALLOCATE(N(n_num))
   N(1) = 500
   N(2) = 850
   N(3) = 1000
   N(4) = 1200
   N(5) = 1350
   N(6) = 1500
   N(7) = 1700
   N(8) = 2000
ALLOCATE(nrun(1))
nrun = 1000

ALLOCATE(x(nrun(1)), y(nrun(1)), z(nrun(1)))

OPEN(7,file="es_2.dat",STATUS="REPLACE", ACTION="WRITE")
OPEN(1,file="hist.dat",STATUS="REPLACE", ACTION="WRITE")

ALLOCATE(P_N(n_num, -500000:500000), H_u(n_num, -nbin/2:nbin/2),&
H_exp(n_num, -2*nbin:2*nbin), H_L(n_num, -500000:500000), &
L1(n_num, -500000:500000), L2(n_num, -500000:500000))

dr      = float(2)/nbin                   !Interval length uniform distribution is 2
H_u     = 0
H_exp   = 0
H_L     = 0
L1      = 0
L2      = 0


WRITE(7,*)"#CENTRAL LIMIT RESULTS FOR UNIFORM DISTRIBUTION"           
WRITE(7,*)"#N, mu_a, mu_n, sigma_a, sigma_n, <z^4>, 3<z^2>^2"           
DO k = 1, n_num
    WRITE(1,*)"#INDEX ", k-1
    ALLOCATE(rnd(N(k)))
    x         = 0                                     !Set the accumulator
    sigma_a_u = s_u/sqrt(float(N(k)))          !Sigma is sigma/N
    
    DO i = 1, nrun(1)
        CALL random_number(rnd)
        rnd = rnd*2 - 1.                 !Shift to uniform distribution in [-1;1] 
        x(i)  = SUM(rnd/N(k))            !Don't work with Sum/N
                                         !bc SUM can be too big for big N
    END DO
    
    mu_n    = SUM(x/nrun(1))
    sigma_n = sqrt(SUM(x**2/nrun(1)) - mu_n**2)
    z       = (x - mu_a_u)/sigma_a_u
    DEALLOCATE(rnd)
    WRITE(7,*)N(k), mu_a_u, mu_n, sigma_a_u, sigma_n, SUM(z**4/nrun(1)), 3*(SUM(z**2/nrun(1))**2)       !Calculate z moments
    
    !HERE WE FILL HISTOGRAMS
    
    DO i = 1, nrun(1)
       j         = int(x(i)/dr)          !See in which bin falls x(i) and fill Histograms
       H_u(k, j) = H_u(k, j) + 1          
    END DO
    DO j = -nbin/2, nbin/2
       CALL gauss(mu_a_u, sigma_a_u**2, (j*dr), P_N(k, j))    !Standard Gaussian
       WRITE(1,*)(j*dr), H_u(k, j)/(dr*nrun(1)), P_N(k, j)         
    END DO
    WRITE(1,*)""
    WRITE(1,*)""
    
END DO
WRITE(7,*)""
WRITE(7,*)""



x       = 0                !Set value to zero, not needed but just in case
z       = 0
mu_n    = 0 
sigma_n = 0                !Set value to zero, not needed but just in case

WRITE(7,*)"#CENTRAL LIMIT RESULTS FOR EXPONENTIAL DISTRIBUTION"           
WRITE(7,*)"#N, mu_a, mu_n, sigma_a, sigma_n, <z^4>, 3<z^2>^2"           
DO k = 1, n_num
    WRITE(1,*)"#INDEX ", n_num + k-1
    ALLOCATE(rnd(N(k)))
    x         = 0                                     !Set the accumulator
    sigma_a_e = s_e/sqrt(float(N(k)))          !Sigma is sigma/N
    
    DO i = 1, nrun(1)
        DO j = 1, N(k)
            CALL expdev(rnd(j))
        END DO    
        x(i)  = SUM(rnd/N(k))            !Don't work with Sum/N
                                         !bc SUM can be too big for big N 
    END DO
    
    mu_n    = SUM(x/nrun(1))
    sigma_n = sqrt(SUM(x**2/nrun(1)) - mu_n**2)
    z       = (x - mu_a_e)/sigma_a_e 
    DEALLOCATE(rnd)
    WRITE(7,*)N(k), mu_a_e, mu_n, sigma_a_e, sigma_n, SUM(z**4/nrun(1)), 3*(SUM(z**2/nrun(1))**2)       !Calculate z moments


    !HERE WE FILL HISTOGRAMS

    DO i = 1, nrun(1)
       j      = int(x(i)/dr)               !See in which bin falls x(i) and fill Histograms
       H_exp(k, j) = H_exp(k, j) + 1          
    END DO

    DO j = 0, nbin
       CALL gauss(mu_a_e, sigma_a_e**2, (j*dr), P_N(k, j))    !Standard Gaussian
       WRITE(1,*)(j*dr), H_exp(k, j)/(dr*nrun(1)), P_N(k, j)         
    END DO
    WRITE(1,*)""
    WRITE(1,*)""

END DO
WRITE(7,*)""
WRITE(7,*)""


x       = 0                !Set value to zero, not needed but just in case
z       = 0
mu_n    = 0 
sigma_n = 0                !Set value to zero, not needed but just in case

WRITE(7,*)"#CENTRAL LIMIT RESULTS FOR LORENTZ DISTRIBUTION"           
WRITE(7,*)"#N, mu_a, mu_n, sigma_a, sigma_n, <z^4>, 3<z^2>^2"           
DO k = 1, n_num
    WRITE(1,*)"#INDEX ", 2*n_num + k -1
    ALLOCATE(rnd(N(k)))
    x         = 0                                     !Set the accumulator
    sigma_a_L = s_L/sqrt(float(N(k)))          !Sigma is sigma/N
    
    DO i = 1, nrun(1)
        DO j = 1, N(k)
            CALL cauchy(rnd(j))
        END DO    
        x(i)  = SUM(rnd/N(k))            !Don't work with Sum/N
                                         !bc SUM can be too big for big N 
    END DO
    
    mu_n    = SUM(x/nrun(1))
    sigma_n = sqrt(SUM(x**2/nrun(1)) - mu_n**2)
    z       = (x - mu_a_e)/sigma_a_e 
    WRITE(7,*)N(k), mu_a_L, mu_n, sigma_a_L, sigma_n, SUM(z**4/nrun(1)), 3*(SUM(z**2/nrun(1))**2)       !Calculate z moments


    !HERE WE FILL HISTOGRAMS

    DO i = 1, nrun(1)
       j         = int(x(i)/dr)               !See in which bin falls x(i) and fill Histograms
       H_L(k, j) = H_L(k, j) + 1 
    END DO

    DO i = 1, N(k)
       j       = int(rnd(i)/dr)
       L1(k, j) = L1(k, j) + 1
    END DO
    
    DO j = -5*nbin, 5*nbin
       CALL gauss(mu_a_L, sigma_a_L**2, (j*dr), P_N(k, j))    !Standard Gaussian
       CALL lorentz((j*dr), L2(k, j))                     !LORENTZ
       WRITE(1,*)(j*dr), H_L(k, j)/(dr*nrun(1)), P_N(k, j), L1(k, j)/(dr*N(k)), L2(k, j)         
    END DO
    WRITE(1,*)""
    WRITE(1,*)""
    
    DEALLOCATE(rnd)
END DO
WRITE(7,*)""
WRITE(7,*)""
PRINT*,'Results stored in "es_2.dat"'
PRINT*,'Histogram data stored in "hist.dat"'
CLOSE(1)
CLOSE(7)    


STOP

DEALLOCATE(P_N, H_u, H_exp, nrun)









CONTAINS



!-----------------------------------------------------
!Subroutine to have point x from e^-x distribution
!-----------------------------------------------------


SUBROUTINE expdev(x)
REAL(kind=nkind), INTENT(OUT) :: x
REAL(kind=nkind)              :: r
    DO
        CALL random_number(r)
        IF(r > 0) exit
    END DO
    x = -log(r)
END SUBROUTINE expdev





!-----------------------------------------------------
!Subroutine to have point y from Gauss distribution
!-----------------------------------------------------


SUBROUTINE gauss(mu, var, x, y)
!INTEGER, INTENT(IN) :: x
REAL(kind=nkind), INTENT(IN)    :: mu, var, x
REAL(kind=nkind), INTENT(OUT)   :: y
REAL(kind=nkind), PARAMETER     :: pi = 3.14159265359


y = ((1/(2.*pi*var))**0.5)*exp(-(mu - x)**2/(2*var))

END SUBROUTINE gauss


!-----------------------------------------------------
!Subroutine to have point y from Lorentz distribution
!-----------------------------------------------------


SUBROUTINE lorentz(x, y)
!INTEGER, INTENT(IN) :: x
REAL(kind=nkind), INTENT(IN)    :: x
REAL(kind=nkind), INTENT(OUT)   :: y
REAL(kind=nkind), PARAMETER     :: pi = 3.14159265359


y = float(1)/(pi*(1 + x**2))

END SUBROUTINE lorentz



!------------------------------------------------------------------------
!!Generate a random deviate from the standard Cauchy distribution
!------------------------------------------------------------------------

SUBROUTINE cauchy(fn_val)
REAL(kind=nkind), INTENT(OUT)     :: fn_val
REAL(kind=nkind)                  :: v(2)
!Variables needed for test
REAL(kind=nkind)                  :: half, one, two, vsmalL

half   = 0.5
one    = 1.0
two    = 2.0
vsmall = TINY(1.0)


DO
    DO
      CALL random_number(v)
      v = two*(v - half)
      IF (abs(v(2)) < vsmall) CYCLE               ! Test for zero
      IF (v(1)**2 + v(2)**2 < one) EXIT
    END DO
    fn_val = v(1) / v(2)
    IF(fn_val < 1.85) then
        IF(fn_val > -1.85) EXIT
    END IF
END DO

END SUBROUTINE cauchy


END PROGRAM CLT
