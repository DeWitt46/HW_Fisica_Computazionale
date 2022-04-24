!------------------------------------------------------------------------------
! METROPOLIS generation of random numbers with a Gaussian distribution
! P(x) = exp(-x**2/(2*sigma**2))/sqrt(2*pi*sigma**2)
!------------------------------------------------------------------------------

PROGRAM es_1
IMPLICIT NONE
INTEGER, PARAMETER                             :: nkind = selected_real_kind(13)
INTEGER, PARAMETER                             :: maxN = 1, maxD = 3, Nmax = 100000, lmax=30
REAL(kind=nkind), PARAMETER                    :: Dmax = 8

INTEGER                                        :: i, j, k, l, h
INTEGER                                        :: ibin, maxbin, m, Nmin
INTEGER, DIMENSION(maxN)                       :: N

REAL(kind=nkind)                               :: sigma,rnd,Dmin,x0,deltaisto,y
REAL(kind=nkind), DIMENSION(maxD)              :: delta
REAL(kind=nkind)                               :: x,z,x1,x1_BM,x2,x3,x4,xp,expx,expxp,w,acc
  
REAL(kind=nkind), DIMENSION(:), ALLOCATABLE    :: istog
REAL(kind=nkind), DIMENSION(0:lmax-1)          :: C = 0._nkind, C_BM = 0._nkind
REAL(kind=nkind), DIMENSION(0:lmax-1)          :: corr = 0._nkind, corr_BM = 0._nkind

LOGICAL, SAVE                                  :: log_plot = .false.                    !Flag to choose N for log or linear plot


INTEGER                                     :: sizer, values(1:8)                !Seed variables
INTEGER, DIMENSION(:), ALLOCATABLE          :: seed, box                         !Seed variables


PRINT*, "insert Nmin >"
READ*, Nmin
PRINT*,"insert sigma >"
READ*, sigma
PRINT*,"insert x0 >"
READ*, x0
PRINT*,"insert delta_min >"
READ*, Dmin
PRINT*,"insert max number of bin in the histogram  >"
READ*, maxbin

ALLOCATE(istog(-maxbin/2:maxbin/2))
istog = 0._nkind
deltaisto = 10.*sigma/maxbin  !histogram over a range of 10*sigma

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




IF(log_plot) THEN
   DO j = 1, maxN
      N(j) = 10**(log10(dfloat(Nmin)) + (j-1))
      IF(N(j) > Nmax)N(j) = Nmax
   END DO
   DO k = 1, maxD
      delta(k) = 10**(k*log10(Dmax)/maxD)
      IF(maxD == 1)delta(k)=Dmin
   END DO   
ELSE
   DO j = 1, maxN
      IF(maxN /= 1) THEN
         N(j) = Nmin + (j-1)*(Nmax - Nmin)/(maxN-1)
      ELSE
         N(j) = Nmin
      END IF   
   END DO
   DO k = 1, maxD
      IF(maxD /= 1) THEN
         delta(k) = Dmin + (k-1)*(Dmax - Dmin)/(maxD-1)
      ELSE
         delta(k) = Dmin
      END IF   
   END DO
END IF


OPEN(7,file='gauss_metropolis.dat',status='replace')
OPEN(3,file='corr_metropolis.dat',status='replace')
OPEN(1,file='hist_metropolis.dat',status='replace')

WRITE(unit=7,fmt=*)"N, x0, delta, acc ratio, <x>, <x> th, <x^2>, <x^2> th, sigma, sigma th, <x^3>, <x^3> th,<x^4>, <x^4> th"
WRITE(unit=3,fmt=*)"N, x0, delta, l, corr(l), corr_BM(l), <x_BM>"

DO j = 1, maxN
   DO k = 1, maxD
!------------------------------  
!Set accumulators
!------------------------------
      acc     = 0.0_nkind
      x       = x0
      z       = x0
      x1      = 0.0_nkind
      x1_BM   = 0.0_nkind
      x2      = 0.0_nkind
      x3      = 0.0_nkind
      x4      = 0.0_nkind
      corr    = 0.0_nkind
      corr_BM = 0.0_nkind
      istog   = 0.

!------------------------------  
!ALGORITHM
!------------------------------
      DO i = 1, N(j)
         x1    = x1 + x
         x1_BM = x1_BM + z
         x2    = x2 + x**2
         x3    = x3 + x**3
         x4    = x4 + x**4

         expx  = - x**2 /(2*sigma**2)            !
         CALL random_number(rnd)                 !
         xp    = x + delta(k) * (rnd-0.5_nkind)  !
         expxp = - xp**2 /(2*sigma**2)           !   metropolis
         w     = exp (expxp-expx)                !   algorithm
         CALL random_number(rnd)                 !
         IF(w > rnd) THEN                        !
            x = xp                               !
      
         acc = acc + 1.0_nkind                
         END IF
         
         
         
         CALL gasdev(z)
         C(mod(i,lmax))    = x
         C_BM(mod(i,lmax)) = z
         
         DO l = 1, lmax
            IF(i >= l) THEN
               DO h = 0, l-1
                  corr(h)    = corr(h) + C(mod((i - l), l))*C(mod((i - l + h), l))
                  corr_BM(h) = corr_BM(h) + C_BM(mod((i - l), l))*C_BM(mod((i - l + h), l))
               END DO
            END IF   
         END DO   

         ibin = nint(x/deltaisto)                !Choose bin and fill It
         IF(abs(ibin) < maxbin/2) istog(ibin) = istog(ibin) + 1
      END DO
      corr    = (corr - x1**2)/corr(0)
      corr_BM = (corr_BM - x1_BM**2)/corr_BM(0)
      
!-------------------------------------------------------------
!WRITE RESULTS
!-------------------------------------------------------------
      WRITE(unit=7,fmt=*)"#INDEX", (j-1)*maxD + k-1
      WRITE(unit=7,fmt=*)N(j), x0, delta(k), acc/N(j), x1/N(j), 0.0_nkind,&
      x2/N(j), sigma**2, x2/N(j)-(x1/N(j))**2, sigma**2,&
      x3/N(j), 0.0_nkind, x4/N(j), 3.0_nkind*sigma**4

      WRITE(unit=3,fmt=*)"#INDEX", (j-1)*maxD + k-1
      DO l = 0, lmax-1
         WRITE(unit=3,fmt=*)N(j), x0, delta(k), l, corr(l), corr_BM(l), x1_BM 
      END DO
      
!-------------------------------------------------------------  
!WRITE HISTOGRAMS
!-------------------------------------------------------------
      WRITE(unit=1,fmt=*)"#INDEX", (j-1)*maxD + k-1
      WRITE(unit=1,fmt=*)"# N, x0, delta = ", N(j), x0, delta(k)
      DO ibin = -maxbin/2, maxbin/2
         CALL gauss(0._nkind, sigma**2, dfloat(ibin), y)
         WRITE(1,*)ibin*deltaisto, istog(ibin)/real(N(j))/deltaisto, y
      END DO
      
!WE WANT GNUPLOT INDECES      
      WRITE(unit=7,fmt=*)""
      WRITE(unit=7,fmt=*)""      
      WRITE(unit=3,fmt=*)""
      WRITE(unit=3,fmt=*)""
   END DO
   
   WRITE(unit=1,fmt=*)""   
   WRITE(unit=1,fmt=*)""
END DO
PRINT*,'seed data stored in "gauss_metropolis.dat"'
PRINT*,'hist data stored in "hist_metropolis.dat"'

CLOSE(1)
CLOSE(7)
DEALLOCATE(istog)


CONTAINS


!-------------------------------------------------------------------
!Subroutine to have point y from Gauss distribution re-normalized
!-------------------------------------------------------------------
SUBROUTINE gauss(mu, var, x, y)
IMPLICIT NONE
REAL(kind=nkind), INTENT(IN)    :: mu, var, x
REAL(kind=nkind), INTENT(OUT)   :: y
REAL(kind=nkind), PARAMETER     :: pi = 3.14159265359


y = ((1._nkind/(2._nkind*pi*var))**0.5)*exp(-(mu - x)**2/(2*var))

END SUBROUTINE gauss

!-----------------------------------------------------------------------------
!Subroutine to have random points following gauss with Box-Muller algorithm
!-----------------------------------------------------------------------------
SUBROUTINE gasdev(rnd)
   IMPLICIT NONE
   REAL(kind=nkind), INTENT(OUT) :: rnd
   REAL(kind=nkind)              :: r2,a,b
   REAL(kind=nkind), SAVE        :: g
   
   LOGICAL, SAVE :: gaus_stored=.false.
   IF (gaus_stored) THEN
      rnd=g
      gaus_stored=.false.
   ELSE
      DO
         CALL random_number(a)
         CALL random_number(b)
         a=2._nkind*a-1._nkind
         b=2.*b-1._nkind
         r2=a**2+b**2
         if (r2 > 0. .and. r2 < 1.) exit
      end do
      r2=sqrt(-2._nkind*log(r2)/r2)
      rnd=a*r2
      g=b*r2
      gaus_stored=.true.
   end if
END SUBROUTINE gasdev


END PROGRAM es_1
