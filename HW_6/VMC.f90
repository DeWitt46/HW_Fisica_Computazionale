PROGRAM VMC
!-----------------------------------------------------------------------------
!Code to compare GS of harmonic and anharmonic oscillator esteemed through
!VMC with exponential and parabolic trial function
!-----------------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER, PARAMETER        :: dp=selected_real_kind(13)

   INTEGER, PARAMETER        :: N = 100000, b_step=350, a_step=20
   REAL(kind=dp), PARAMETER  :: b_i = 0.56, b_f = 0.6, a_i=1.05, a_f=3.05
   INTEGER                   :: j
   REAL(kind=dp)             :: a, beta, y, sigma, x0 = 0._dp, delta=6._dp, g
   LOGICAL, SAVE             :: ISO = .false. !Flag to look at harmonic or anharmonic oscillator
   

IF(ISO) THEN
   OPEN (unit=1,file="Metropolis_exp.dat",status="replace",action="write")
   WRITE(unit=1,fmt=*)"#N, delta, beta, etot, sigma, acc"     
   OPEN (unit=7,file="Metropolis_par.dat",status="replace",action="write")
   WRITE(7,fmt=*)"#N, delta, a, etot, sigma, acc"
   CLOSE(1)
   CLOSE(7)

   DO j = 0, b_step
      beta  = b_i + j*(b_f-b_i)/b_step
      sigma = sqrt(1._dp/(4*beta))
      CALL metropolis_exp(N, sigma, x0, delta)
      PRINT*,j, "iteration done"
   END DO

   DO j = 0, a_step
      a = a_i + j*(a_f-a_i)/a_step
      CALL metropolis_par(a, N, x0, delta)
      PRINT*,j, "iteration done"
   END DO
   
   !For Potential plots
   OPEN(unit=5,file="pot.dat",status="replace",action="write")
   WRITE(unit=5,fmt=*)"#x, V(x), psi_exp(x), psi_par(x)"
   a = 2.0533333460489906_dp
   sigma = sqrt(1._dp/4._dp*0.49999999999999989_dp)
   DO j = 0, 300
      y = -a + (j*2._dp*a)/300._dp
      CALL gauss(0._dp, sigma**2, y, g)
      WRITE(unit=5,fmt=*)y, y**2/4._dp, g, sqrt(15._dp/(16._dp*a**5))*(a**2 - y**2)
   END DO 
   CLOSE(5)   
ELSE
   OPEN (unit=1,file="Metropolis_aniso_exp.dat",status="replace",action="write")
   WRITE(unit=1,fmt=*)"#N, delta, beta, etot, sigma, acc"     
   OPEN (unit=7,file="Metropolis_aniso_par.dat",status="replace",action="write")
   WRITE(7,fmt=*)"#N, delta, a, etot, sigma, acc"
   CLOSE(1)
   CLOSE(7)

   DO j = 0, b_step
      beta  = b_i + j*(b_f-b_i)/b_step
      sigma = sqrt(1._dp/(4*beta))
      CALL metropolis_aniso_exp(N, sigma, x0, delta)
      PRINT*,j, "iteration done"
   END DO

   DO j = 0, a_step
      a = a_i + j*(a_f-a_i)/a_step
      CALL metropolis_aniso_par(a, N, x0, delta)
      PRINT*,j, "iteration done"
   END DO
   
   !For Potential plots
   OPEN(unit=3,file="Aniso_pot.dat",status="replace",action="write")
   WRITE(unit=3,fmt=*)"#x, V(x), psi_exp(x), psi_par(x)"
   a = 1.8166666666666667_dp
   sigma = sqrt(1._dp/4._dp*0.62426667372385658_dp)
   DO j = 0, 300
      y = -a + (j*2._dp*a)/300._dp
      CALL gauss(0._dp, sigma**2, y, g)
      WRITE(unit=3,fmt=*)y, y**2/4._dp + y**4/8._dp, g, sqrt(15._dp/(16._dp*a**5))*(a**2 - y**2)
   END DO
   CLOSE(3)
END IF   
   
   
   
CONTAINS

!-------------------------------------------------------------------
!Subroutine to have point y from Gauss distribution re-normalized
!-------------------------------------------------------------------
SUBROUTINE gauss(mu, var, x, y)
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN)    :: mu, var, x
   REAL(kind=dp), INTENT(OUT)   :: y
   REAL(kind=dp), PARAMETER     :: pi = 3.14159265359

   y = ((1._dp/(2._dp*pi*var))**0.5)*exp(-(mu - x)**2/(2*var))

END SUBROUTINE gauss






!---------------------------------------------------------------------
! METROPOLIS sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2),
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = exp(-x^2/(4\sigma^2))
! \sigma=1 => psi^2(x) = costant * standard gaussian
!  P(x) =  exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!----------------------------------------------------------------------
SUBROUTINE metropolis_exp(N, sigma, x0, delta)
   IMPLICIT NONE
   INTEGER, PARAMETER        :: dp=selected_real_kind(13)
   INTEGER                   :: i
   INTEGER, INTENT(IN)       :: N
   REAL(kind=dp), INTENT(IN) :: sigma, x0, delta
   REAL(kind=dp)             :: etot,etot2,ekin,ekinL,epot,epotL,rnd
   REAL(kind=dp)             :: x,x1,x2,xp,expx,expxp,p,acc
   CHARACTER(len=13), SAVE   :: format1 = "(a7,2x,2f9.5)"


   INTEGER                             :: sizer, values(1:8)                !Seed variables
   INTEGER, DIMENSION(:), ALLOCATABLE  :: seed                              !Seed variables


!--------------------------------------------
!Seed generation
!--------------------------------------------

   CALL random_seed(sizer)
   ALLOCATE(seed(sizer))
   CALL date_and_time(values=values)
   seed(:) = values(8)
   CALL random_seed(put=seed)


   acc   = 0.0_dp
   x     = x0
   x1    = 0.0_dp
   x2    = 0.0_dp
   ekin  = 0.0_dp
   epot  = 0.0_dp
   etot  = 0.0_dp
   etot2 = 0.0_dp   
   
   DO i = 1, N
      ekinL = - 0.5_dp * ((x/(2*sigma**2))**2 - 1/(2*sigma**2))
      epotL = 0.5_dp * x**2
      ekin  = ekin + ekinL
      epot  = epot + epotL
      etot  = etot + ekinL + epotL
      etot2 = etot2 + (ekinL + epotL)**2
      
      x1 = x1 + x
      x2 = x2 + x**2
      !-------------------------------
      expx = - x**2 /(2*sigma**2)    !
      CALL random_number(rnd)        !
      xp = x + delta * (rnd-0.5_dp)  !
      expxp = - xp**2 /(2*sigma**2)  !   metropolis
      p = exp (expxp-expx)           !   algorithm
      CALL random_number(rnd)        !
      IF (p > rnd) THEN              !
         x = xp                      !
      !-------------------------------
         acc = acc + 1.0_dp                
      END IF
   END DO
  
!   WRITE(unit=*,fmt=*)"acceptance ratio = ",acc/n
!   WRITE(unit=*,fmt=*)"# Results (simulation vs. exact results):"
!   WRITE(unit=*,fmt=format1)"etot = ",etot/n,1.0_dp/(8.0_dp*sigma**2)&
!       +0.5_dp*sigma**2
!   WRITE(unit=*,fmt=format1)"ekin = ",ekin/n,1.0_dp/(8.0_dp*sigma**2)
!   WRITE(unit=*,fmt=format1)"epot = ",epot/n,0.5_dp*sigma**2
!   WRITE(unit=*,fmt=format1)"<x>  = ",x1/n,0.0_dp
!   WRITE(unit=*,fmt=format1)"<x^2>= ",x2/n,sigma**2
  
  
   OPEN (unit=1,file="Metropolis_exp.dat", position="append")  
   WRITE(unit=1,fmt=*)N, delta, 1._dp/(4*sigma**2), etot/N, sqrt(etot2/N - (etot/N)**2), acc/N
   
   CLOSE(1)
END SUBROUTINE metropolis_exp









!-------------------------------------------------------------------------
! METROPOLIS sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2),
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = B(a^2-x^2) for |x|<a; 0 elsewhere
!-------------------------------------------------------------------------

SUBROUTINE metropolis_par(a, N, x0, delta)
   IMPLICIT NONE
   INTEGER, PARAMETER        :: dp=selected_real_kind(13)
   INTEGER                   :: i
   INTEGER, INTENT(IN)       :: N
   REAL(kind=dp), INTENT(IN) :: a, x0, delta

   REAL(kind=dp):: etot,ekin,epot,rnd,ekinL,epotL,etot2,var,ekin_th,epot_th,etot_th
   REAL(kind=dp):: x,x1,x2,xp,psi,psip,p,acc,var_th
   CHARACTER(len=13), save :: format1 = "(a7,2x,1f9.5)"
   CHARACTER(len=17), save :: format2 = "(a17,2(2x,1f9.5))"

   INTEGER                             :: sizer, values(1:8)                !Seed variables
   INTEGER, DIMENSION(:), ALLOCATABLE  :: seed                              !Seed variables

   acc   = 0.0_dp
   x     = x0
   x1    = 0.0_dp
   x2    = 0.0_dp
   ekin  = 0.0_dp
   epot  = 0.0_dp
   etot2 = 0.0_dp
   
   IF((abs(x) >= a)) THEN
      PRINT*, "abs(x) < a"
      STOP
   END IF   

!--------------------------------------------
!Seed generation
!--------------------------------------------

   CALL random_seed(sizer)
   ALLOCATE(seed(sizer))
   CALL date_and_time(values=values)
   seed(:) = values(8)
   CALL random_seed(put=seed)



   DO i = 1, N
      ekinL = 1/(a**2-x**2)
      epotL = 0.5_dp * x**2
      ekin = ekin + ekinL
      epot = epot + epotL
      etot = ekin + epot
      etot2 = etot2 + (ekinL + epotL)**2
      x1 = x1 + x
      x2 = x2 + x**2
      !-------------------------------
      psi = a**2-x**2                !
      CALL random_number(rnd)        !
      xp = x + delta * (rnd-0.5_dp)  !
      IF (abs(xp) >= a) THEN         !
         psip = 0.0_dp !             !
      ELSE                           !
            psip = a**2-xp**2        !   metropolis
      END IF                         !
      p = (psip/psi)**2              !   algorithm
      CALL random_number(rnd)        !
      IF (p > rnd) then              !
         x = xp                      !
      !-------------------------------
         acc = acc + 1.0_dp                
      END IF
   END DO

   ekin_th = 5./4/a**2
   epot_th = a**2/14
   etot_th = ekin_th + epot_th
   var_th = (15./8/a**5) * ( a + 2*a**9/(7*45) + 2*a**5/15) - &
           (5./4/a**2+a**2/14)**2

!   WRITE(unit=*,fmt=*)"acceptance ratio = ",acc/n
!   WRITE(unit=*,fmt=*)"# Results:"
!   WRITE(unit=*,fmt=format2)"etot (num./th.)= ",etot/n,etot_th
!   WRITE(unit=*,fmt=format2)"ekin (num./th.)= ",ekin/n,ekin_th
!   WRITE(unit=*,fmt=format2)"epot (num./th.)= ",epot/n,epot_th
!   WRITE(unit=*,fmt=format2)"evar (num./th.)= ",etot2/n-(etot/n)**2,var_th
!   WRITE(unit=*,fmt=format1)"<x>  = ",x1/n
!   WRITE(unit=*,fmt=format1)"<x^2>= ",x2/n
   
   OPEN(unit=7,file='Metropolis_par.dat',position='append')
   WRITE(7,*)N, delta, a, etot/N, sqrt(etot2/N-(etot/N)**2), acc/N
   CLOSE(7)
   
END SUBROUTINE metropolis_par













!---------------------------------------------------------------------
! METROPOLIS sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2 + 1/8x^4),
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = exp(-x^2/(4\sigma^2))
! \sigma=1 => psi^2(x) = costant * standard gaussian
!  P(x) =  exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!----------------------------------------------------------------------
SUBROUTINE metropolis_aniso_exp(N, sigma, x0, delta)
   IMPLICIT NONE
   INTEGER, PARAMETER        :: dp=selected_real_kind(13)
   INTEGER                   :: i
   INTEGER, INTENT(IN)       :: N
   REAL(kind=dp), INTENT(IN) :: sigma, x0, delta
   REAL(kind=dp)             :: etot,etot2,ekin,ekinL,epot,epotL,rnd
   REAL(kind=dp)             :: x,x1,x2,xp,expx,expxp,p,acc
   CHARACTER(len=13), SAVE   :: format1 = "(a7,2x,2f9.5)"


   INTEGER                             :: sizer, values(1:8)                !Seed variables
   INTEGER, DIMENSION(:), ALLOCATABLE  :: seed                              !Seed variables


!--------------------------------------------
!Seed generation
!--------------------------------------------

   CALL random_seed(sizer)
   ALLOCATE(seed(sizer))
   CALL date_and_time(values=values)
   seed(:) = values(8)
   CALL random_seed(put=seed)


   acc   = 0.0_dp
   x     = x0
   x1    = 0.0_dp
   x2    = 0.0_dp
   ekin  = 0.0_dp
   epot  = 0.0_dp
   etot  = 0.0_dp
   etot2 = 0.0_dp   
   
   DO i = 1, N
      ekinL = - 0.5_dp * ((x/(2*sigma**2))**2 - 1/(2*sigma**2))
      epotL = 0.5_dp * x**2 + 1._dp/8._dp*x**4
      ekin  = ekin + ekinL
      epot  = epot + epotL
      etot  = etot + ekinL + epotL
      etot2 = etot2 + (ekinL + epotL)**2
      
      x1 = x1 + x
      x2 = x2 + x**2
      !-------------------------------
      expx = - x**2 /(2*sigma**2)    !
      CALL random_number(rnd)        !
      xp = x + delta * (rnd-0.5_dp)  !
      expxp = - xp**2 /(2*sigma**2)  !   metropolis
      p = exp (expxp-expx)           !   algorithm
      CALL random_number(rnd)        !
      IF (p > rnd) THEN              !
         x = xp                      !
      !-------------------------------
         acc = acc + 1.0_dp                
      END IF
   END DO
  
!   WRITE(unit=*,fmt=*)"acceptance ratio = ",acc/n
!   WRITE(unit=*,fmt=*)"# Results (simulation vs. exact results):"
!   WRITE(unit=*,fmt=format1)"etot = ",etot/n,1.0_dp/(8.0_dp*sigma**2)&
!       +0.5_dp*sigma**2
!   WRITE(unit=*,fmt=format1)"ekin = ",ekin/n,1.0_dp/(8.0_dp*sigma**2)
!   WRITE(unit=*,fmt=format1)"epot = ",epot/n,0.5_dp*sigma**2
!   WRITE(unit=*,fmt=format1)"<x>  = ",x1/n,0.0_dp
!   WRITE(unit=*,fmt=format1)"<x^2>= ",x2/n,sigma**2
  
  
   OPEN (unit=1,file="Metropolis_aniso_exp.dat", position="append")  
   WRITE(unit=1,fmt=*)N, delta, 1._dp/(4*sigma**2), etot/N, sqrt(etot2/N - (etot/N)**2), acc/N    
   
   
   CLOSE(1)
END SUBROUTINE metropolis_aniso_exp














!-------------------------------------------------------------------------
! METROPOLIS sampling of several physical observables for the
! hamiltonian:         h = -1/2 \nabla^2 + 1/2 x^2 + 1/8x^4,
! comparison exact expected results with numerical results
! on psi^2(x), with  psi(x) = B(a^2-x^2) for |x|<a; 0 elsewhere
!-------------------------------------------------------------------------

SUBROUTINE metropolis_aniso_par(a, N, x0, delta)
   IMPLICIT NONE
   INTEGER, PARAMETER        :: dp=selected_real_kind(13)
   INTEGER                   :: i
   INTEGER, INTENT(IN)       :: N
   REAL(kind=dp), INTENT(IN) :: a, x0, delta

   REAL(kind=dp):: etot,ekin,epot,rnd,ekinL,epotL,etot2,var,ekin_th,epot_th,etot_th
   REAL(kind=dp):: x,x1,x2,xp,psi,psip,p,acc,var_th
   CHARACTER(len=13), save :: format1 = "(a7,2x,1f9.5)"
   CHARACTER(len=17), save :: format2 = "(a17,2(2x,1f9.5))"

   INTEGER                             :: sizer, values(1:8)                !Seed variables
   INTEGER, DIMENSION(:), ALLOCATABLE  :: seed                              !Seed variables

   acc   = 0.0_dp
   x     = x0
   x1    = 0.0_dp
   x2    = 0.0_dp
   ekin  = 0.0_dp
   epot  = 0.0_dp
   etot2 = 0.0_dp
   
   IF((abs(x) >= a)) THEN
      PRINT*, "abs(x) < a"
      STOP
   END IF   

!--------------------------------------------
!Seed generation
!--------------------------------------------

   CALL random_seed(sizer)
   ALLOCATE(seed(sizer))
   CALL date_and_time(values=values)
   seed(:) = values(8)
   CALL random_seed(put=seed)



   DO i = 1, N
      ekinL = 1/(a**2-x**2)
      epotL = 0.5_dp * x**2 + 1._dp/8._dp*x**4
      ekin = ekin + ekinL
      epot = epot + epotL
      etot = ekin + epot
      etot2 = etot2 + (ekinL + epotL)**2
      x1 = x1 + x
      x2 = x2 + x**2
      !-------------------------------
      psi = a**2-x**2                !
      CALL random_number(rnd)        !
      xp = x + delta * (rnd-0.5_dp)  !
      IF (abs(xp) >= a) THEN         !
         psip = 0.0_dp !             !
      ELSE                           !
            psip = a**2-xp**2        !   metropolis
      END IF                         !
      p = (psip/psi)**2              !   algorithm
      CALL random_number(rnd)        !
      IF (p > rnd) then              !
         x = xp                      !
      !-------------------------------
         acc = acc + 1.0_dp                
      END IF
   END DO

   ekin_th = 5./4/a**2
   epot_th = a**2/14
   etot_th = ekin_th + epot_th
   var_th = (15./8/a**5) * ( a + 2*a**9/(7*45) + 2*a**5/15) - &
           (5./4/a**2+a**2/14)**2

!   WRITE(unit=*,fmt=*)"acceptance ratio = ",acc/n
!   WRITE(unit=*,fmt=*)"# Results:"
!   WRITE(unit=*,fmt=format2)"etot (num./th.)= ",etot/n,etot_th
!   WRITE(unit=*,fmt=format2)"ekin (num./th.)= ",ekin/n,ekin_th
!   WRITE(unit=*,fmt=format2)"epot (num./th.)= ",epot/n,epot_th
!   WRITE(unit=*,fmt=format2)"evar (num./th.)= ",etot2/n-(etot/n)**2,var_th
!   WRITE(unit=*,fmt=format1)"<x>  = ",x1/n
!   WRITE(unit=*,fmt=format1)"<x^2>= ",x2/n
   
   OPEN(unit=7,file='Metropolis_aniso_par.dat',position='append')
   WRITE(7,*)N, delta, a, etot/N, sqrt(etot2/N-(etot/N)**2), acc/N
   CLOSE(7)
   
END SUBROUTINE metropolis_aniso_par



END PROGRAM VMC









