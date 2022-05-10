!!-------------------------------------------------------------------------
!!Metropolis algorithm to calculate <E>, <M>, in the canonical ensemble
!!(fix T,N,V) with a 2D Ising model
!!
!!Here: k_B = 1
!!       J   = 1
!!2D Ising model has known T_c, we measure T in units of T_c 
!!-------------------------------------------------------------------------
MODULE COMMON
   IMPLICIT NONE
   PUBLIC :: initial,metropolis,DeltaE
   PUBLIC :: data, output
   
   INTEGER, PUBLIC, PARAMETER                   :: double = selected_real_kind(13)
   INTEGER, PUBLIC, PARAMETER                   :: totpoints = 100, N_T = 15, N_L = 1       !Total points in plot 
   REAL(kind = double), PUBLIC, PARAMETER       :: T_i = 1., T_f = 4., T_c = 2.269
   REAL(kind = double), PUBLIC                  :: T,E,M
   REAL(kind = double), PUBLIC, DIMENSION(N_T)  :: e_ave = 0.
   INTEGER, PUBLIC, DIMENSION(:,:), ALLOCATABLE :: spin
   REAL(kind = double), PUBLIC, DIMENSION(-8:8) :: w
   INTEGER, PUBLIC                              :: N,L,nmcs,nequil, npoints,k
   INTEGER, PUBLIC                              :: accept
   LOGICAL, PUBLIC                              :: PBC = .FALSE.

   INTEGER, PUBLIC                                     :: sizer, values(1:8)                !Seed variables
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE          :: seed, box                         !Seed variables



CONTAINS


!!--------------------------------------------
!!SUBROUTINE TO SET LATTICE PROPERTIES
!!--------------------------------------------

   SUBROUTINE setting(nequil, cum)
      INTEGER, INTENT(OUT)                           :: nequil
      REAL(kind = double), DIMENSION(5), INTENT(OUT) :: cum
      INTEGER                                        :: i
      
      PRINT*, "Linear dimension of lattice L ="
      READ*, L
      ALLOCATE(spin(L,L))
      N = L*L
      PRINT*, "# MC steps per spin for equilibrium ="
      READ*, nequil
      nequil = 1000
      PRINT*, "# MC steps per spin for averages ="
      READ*, nmcs
      nmcs = 10000
      npoints = int(nmcs/totpoints)        


!!Initialize key for data plot
      IF(PBC) THEN
         OPEN(unit=1,file="ist_data_pbc.dat", status="replace", action="write")
         OPEN(unit=3,file="eq_data_pbc.dat", status="replace", action="write")
         OPEN(unit=7,file="ist_cum_data_pbc.dat", status="replace", action="write")
         OPEN(unit=8,file='ising-up_pbc.dat', status="replace", action="write")
         OPEN(unit=9,file='ising-down_pbc.dat', status="replace", action="write")
         OPEN(unit=13,file="cdev_pbc.dat", status="replace", action="write")         
         WRITE(unit=1,fmt=*)"#MC_step, E, M, |M|, c, chi"
         WRITE(unit=3,fmt=*)"#T, E_avg, E^2_avg, M_avg, M^2_avg,|M|_avg, c, chi, MC_step"
         WRITE(unit=7,fmt=*)"#MC_step, E_cum, M_cum, |M|_cum"
         WRITE(unit=13,fmt=*)"#T, c_dev"
         CLOSE(1)
         CLOSE(3)
         CLOSE(7)
         CLOSE(8)
         CLOSE(9)
         CLOSE(13)
         
      ELSE  !!OBC

         OPEN(unit=1,file="ist_data_obc.dat", status="replace", action="write")
         OPEN(unit=3,file="eq_data_obc.dat", status="replace", action="write")
         OPEN(unit=7,file="ist_cum_data_obc.dat", status="replace", action="write")
         OPEN(unit=8,file='ising-up_obc.dat', status="replace", action="write")
         OPEN(unit=9,file='ising-down_obc.dat', status="replace", action="write")         
         OPEN(unit=13,file="cdev_obc.dat", status="replace", action="write")                 
         WRITE(unit=1,fmt=*)"#MC_step, E, M, |M|, c, chi"
         WRITE(unit=3,fmt=*)"#T, E_avg, E^2_avg, M_avg, M^2_avg,|M|_avg, c, chi, MC_step"
         WRITE(unit=7,fmt=*)"#MC_step, E_cum, M_cum, |M|_cum"
         WRITE(unit=13,fmt=*)"#T, c_dev"         
         CLOSE(1)
         CLOSE(3)
         CLOSE(7)
         CLOSE(8)
         CLOSE(9)
         CLOSE(13)         
      END IF
!!Seed generation and saving
      CALL random_seed(sizer)
      ALLOCATE(seed(sizer))
      ALLOCATE(box(sizer))       
      CALL date_and_time(values=values)
      seed(:) = values(8)

      PRINT*,'Here the seed has ',sizer,'components'  
      CALL random_seed(put=seed)
      CALL random_seed(get=box)
      OPEN(unit=1,file="seed.dat",status="replace",action="write")
      DO i = 1, sizer
         WRITE(unit=1,fmt=*)box(i)
      END DO
      CLOSE(1)
      PRINT*,'seed data stored in "seed.dat"'

   END SUBROUTINE setting


!!-----------------------------------------------
!!SUBROUTINE TO CONFIGURE INITIAL STATE
!!-----------------------------------------------
   SUBROUTINE initial(temp)      
      INTEGER             :: x,y,up,right,sums,i,dE
      REAL                :: rnd
      REAL(kind=double)   :: temp


      !temp = temp*T_c
      M = 0.0_double
            
!!random initial configuration
!!compute initial magnetization
      IF(temp == T_i) THEN
         DO y = 1, L
            DO x = 1, L
               CALL random_number(rnd)
               IF(rnd < 0.5) THEN
                  spin(x,y) = 1
               ELSE
                  spin(x,y) = -1
               END IF
            END DO
         END DO
      END IF   
!!compute initial energy
      E = 0.0_double
      
      IF(PBC) THEN
!!periodic boundary conditions
         DO y = 1, L
            IF(y == L) THEN
               up = 1
            ELSE
               up = y + 1
            END IF
            DO x = 1, L
               IF(x == L) THEN
                  right = 1
               ELSE
                  right = x + 1
               END IF
               sums = spin(x,up) + spin(right,y)
!!calculate the initial energy summing all over pairs
!!(For a given spin, consider only the up NN and the right NN
!!- NOT the down and the left NN - : each interaction is counted once
               E = E - spin(x,y)*sums
               M = M + spin(x,y)
            END DO
         END DO
!!calculate the transition probability according 
!!to the Boltzmann distribution (exp(-deltaE/KT).
!!Choosing the interaction parameter  J=1, ***IN CASE OF P.B.C.***
!!possible energy variations per spin flip are -8,-4,0,+4,+8: 
      DO dE = -8,8, 4
         w(dE) = exp(-dE/temp)
      END DO
      accept = 0

      ELSE
!!open boundary conditions      
         DO y = 1, L-1
            up = y + 1
            DO x = 1, L-1
               right = x + 1
               sums = spin(x,up) + spin(right,y)
!!calculate the initial energy summing all over pairs
!!(For a given spin, consider only the up NN and the right NN
!!- NOT the down and the left NN - : each interaction is counted once
               E = E - spin(x,y)*sums
               M = M + spin(x,y)
            END DO
         END DO
!!calculate the transition probability according 
!!to the Boltzmann distribution (exp(-deltaE/KT).
!!Choosing the interaction parameter  J=1, ***IN CASE OF O.B.C.***
!!possible energy variations per spin flip are -8,-6,-4,-2,0,+2,+4,+6,+8: 
      DO dE = -8,8, 2
         w(dE) = exp(-dE/temp)
      END DO
      accept = 0
         
      END IF 


END SUBROUTINE initial



!!------------------------------------------------------------------------
!!SUBROUTINE TO IMPLEMENT METROPOLIS ALGORITHM - SPIN-FLIP DYNAMICS
!!------------------------------------------------------------------------
   SUBROUTINE metropolis()
!!one Monte Carlo step per spin
      INTEGER :: ispin,x,y,dE
      REAL    :: rnd
      

      DO ispin = 1, N
!!random x and y coordinates for trial spin
         CALL random_number(rnd)
         x = int(L*rnd) + 1
         CALL random_number(rnd)
         y  = int(L*rnd) + 1
         dE = DeltaE(x,y)
         CALL random_number(rnd)
         IF(rnd <= w(dE)) THEN
            spin(x,y) = -spin(x,y)
            accept = accept + 1
            M = M + 2*spin(x,y)  !!factor 2 is to account for the variation:
            E = E + dE           !!(-(-)+(+))
         END IF
      END DO
      
   END SUBROUTINE metropolis



!!--------------------------------------------------------------
!!FUNCTION TO CALCULATE dE IMPLEMENTING PBC or OBC
!!--------------------------------------------------------------
   FUNCTION DeltaE(x,y) RESULT (DeltaE_result)
!!periodic boundary conditions
      INTEGER, INTENT(IN) :: x,y
      
      INTEGER :: DeltaE_result
      INTEGER :: left
      INTEGER :: right
      INTEGER :: up
      INTEGER :: down
      
      IF(PBC) THEN !!PBC
         IF(x == 1) THEN
            left = spin(L,y)
            right = spin(2,y)
         ELSE IF(x == L) THEN
            left = spin(L-1,y)
            right = spin(1,y)
         ELSE
            left = spin(x-1,y)
            right = spin(x+1,y)
         END IF
         IF(y == 1) THEN
            up = spin(x,2)
            down = spin(x,L)
         ELSE IF(y == L) THEN
            up = spin(x,1)
            down = spin(x,L-1)
         ELSE
            up = spin(x,y+1)
            down = spin(x,y-1)
         END IF

      ELSE !!OBC
         
         IF(x == 1) THEN
            left = 0
            right = spin(2,y)
         ELSE IF(x == L) THEN
            left = spin(L-1,y)
            right = 0
         ELSE
            left = spin(x-1,y)
            right = spin(x+1,y)
         END IF
         IF(y == 1) THEN
            up = spin(x,2)
            down = 0
         ELSE IF(y == L) THEN
            up = 0
            down = spin(x,L-1)
         ELSE
            up = spin(x,y+1)
            down = spin(x,y-1)
         END IF

      END IF
      
      DeltaE_result = 2*spin(x,y)*(left + right + up + down)
!!also here the factor 2 is to account for the variation

   END FUNCTION DeltaE



!!------------------------------------------------------------------------
!!SUBROUTINE TO ACCUMULATE DATA OF INTEREST AFTER EVERY MC STEP
!!------------------------------------------------------------------------
   SUBROUTINE data(cum, imcs)
      INTEGER, INTENT(IN)                              :: imcs
      REAL(kind = double), DIMENSION(5), INTENT(inout) :: cum
      
      INTEGER :: j
      
      cum(1) = cum(1) + E
      cum(2) = cum(2) + E*E
      cum(3) = cum(3) + M
      cum(4) = cum(4) + M*M
      cum(5) = cum(5) + abs(M)

      IF(PBC) THEN
         OPEN(unit=7,file="ist_cum_data_pbc.dat", position="append")  
         WRITE(unit=7,fmt=*)imcs, (cum(j)/real(N)/real(imcs), j=1,5, 2)
         CLOSE(7)
         
      ELSE  !!OBC

         OPEN(unit=7,file="ist_cum_data_obc.dat", position="append")  
         WRITE(unit=7,fmt=*)imcs, (cum(j)/real(N)/real(imcs), j=1,5, 2)
         CLOSE(7)

      END IF
END SUBROUTINE data



!!---------------------------------------------------------------------
!!SUBROUTINE TO OUTPUT THE CUM RESULTS
!!---------------------------------------------------------------------
   SUBROUTINE output(cum, mcs_tot)
      INTEGER, INTENT(IN)                              :: mcs_tot
      REAL(kind = double), DIMENSION(5), INTENT(inout) :: cum
      
      REAL(kind = double) :: eave,e2ave,mave,m2ave,abs_mave,c,chi
      REAL(kind = double) :: acceptance_prob, eave_old
      
      acceptance_prob = real(accept)/real(N)/real(mcs_tot+nequil)
      eave     = cum(1)/real(N)/real(mcs_tot)    !!to avoid interger overflow
      e2ave    = cum(2)/real(N*N)/real(mcs_tot)
      mave     = cum(3)/real(N)/real(mcs_tot)
      m2ave    = cum(4)/real(N*N)/real(mcs_tot)
      abs_mave = cum(5)/real(N)/real(mcs_tot)
      c        = dfloat(N)*(e2ave - eave**2)/(T**2)
      chi      = dfloat(N)*(m2ave - abs_mave**2)/T
      
      IF(PBC) THEN
         OPEN(unit=1,file="ist_data_pbc.dat", position="append")
         WRITE(unit=1,fmt=*)mcs_tot, E/N, M/N, abs(M)/N, c, chi
         CLOSE(1)
               
         IF(mcs_tot == nmcs) THEN
            OPEN(unit=3,file="eq_data_pbc.dat", position="append")
            WRITE(unit=3,fmt=*)T, eave, e2ave, mave, m2ave, abs_mave, c, chi, mcs_tot
            E_ave(k) = eave
            CLOSE(3)
         END IF
      
      ELSE !!OBC
      
         OPEN(unit=1,file="ist_data_obc.dat", position="append")
         WRITE(unit=1,fmt=*)mcs_tot, E/N, M/N, abs(M)/N, c, chi
         CLOSE(1)

         IF(mcs_tot == nmcs) THEN
            OPEN(unit=3,file="eq_data_obc.dat", position="append")
            WRITE(unit=3,fmt=*)T, eave, e2ave, mave, m2ave, abs_mave, c, chi, mcs_tot
            E_ave(k) = eave
            CLOSE(3)
         END IF

      END IF
      
      
   END SUBROUTINE output

END MODULE COMMON




PROGRAM 2D_ising
!!metropolis algorithm for the ising model on a square lattice
   USE COMMON
   REAL(kind = double), DIMENSION(5) :: cum
   REAL(kind = double)               :: c_dev
   
   INTEGER :: imcs,ispin,jspin,z


   OPEN(unit=46,file="obc_try.dat", status="replace", action="write")
   WRITE(46,*)"#L, E_ave/N"
   DO z = 1, N_L
      L = 4 + (z-1)*(26)/(N_L - 1)      !Tried up to L=30 for computational time
      
      CALL setting(nequil,cum)

      DO k = 1, N_T

         c_dev    = 0._double
         cum = 0.0_double
         T   = T_i + (k-1)*(T_f - T_i)/(N_T-1)
         PRINT*,"Now simulating system at temperature T=",T

         CALL initial(T) 
         DO imcs = 1, nequil
            CALL metropolis()
         END DO
   
         IF(PBC) THEN
!!accumulate data while updating spins   
            OPEN(unit=1,file="ist_data_pbc.dat", position="append")  
            OPEN(unit=7,file="ist_cum_data_pbc.dat", position="append")
            WRITE(unit=1,fmt=*)"#EQUILIBRIUM ENDED"
            WRITE(unit=7,fmt=*)"#EQUILIBRIUM ENDED"
            CLOSE(1)
            CLOSE(7)

         ELSE  !OBC

!!accumulate data while updating spins   
            OPEN(unit=1,file="ist_data_obc.dat", position="append")  
            OPEN(unit=7,file="ist_cum_data_obc.dat", position="append")    
            WRITE(unit=1,fmt=*)"#EQUILIBRIUM ENDED"
            WRITE(unit=7,fmt=*)"#EQUILIBRIUM ENDED"
            CLOSE(1)
            CLOSE(7)

         END IF
      
         DO imcs = 1, nmcs
            CALL metropolis()
            CALL data(cum, imcs)
            IF(MOD(imcs, npoints) == 0)CALL output(cum, imcs)      
         END DO
         IF(MOD(nmcs, npoints) /= 0)CALL output(cum, nmcs)


         IF(PBC) THEN
!!SPACES TO INDEX FOR GNUPLOT
            OPEN(unit=1,file="ist_data_pbc.dat", position="append")  
            OPEN(unit=7,file="ist_cum_data_pbc.dat", position="append")    
            WRITE(unit=1,fmt=*)""
            WRITE(unit=1,fmt=*)""
            WRITE(unit=7,fmt=*)""
            WRITE(unit=7,fmt=*)""
            CLOSE(1)
            CLOSE(7)


!!write the coordinates of spins up and down on files for plotting
            OPEN(unit=8,file='ising-up_pbc.dat',position='append')
            OPEN(unit=9,file='ising-down_pbc.dat',position='append')
            WRITE(8,*)"#x, y"
            WRITE(9,*)"#x, y"         
            DO jspin = 1, L
               DO ispin = 1, L
                  IF(spin(ispin,jspin)==1)WRITE(8,*)ispin,jspin
                  IF(spin(ispin,jspin)==-1)WRITE(9,*)ispin,jspin         
               END DO
            END DO
            WRITE(8,*)""
            WRITE(8,*)""
            WRITE(9,*)""
            WRITE(9,*)""         
            CLOSE(8)
            CLOSE(9)
         
         ELSE
!!SPACES TO INDEX FOR GNUPLOT
            OPEN(unit=1,file="ist_data_obc.dat", position="append")  
            OPEN(unit=7,file="ist_cum_data_obc.dat", position="append")    
            WRITE(unit=1,fmt=*)""
            WRITE(unit=1,fmt=*)""
            WRITE(unit=7,fmt=*)""
            WRITE(unit=7,fmt=*)""
            CLOSE(1)
            CLOSE(7)


!!write the coordinates of spins up and down on files for plotting
            OPEN(unit=8,file='ising-up_obc.dat',position='append')
            OPEN(unit=9,file='ising-down_obc.dat',position='append')
            WRITE(8,*)"#x, y"
            WRITE(9,*)"#x, y"
            DO jspin = 1, L
               DO ispin = 1, L         
                  IF(spin(ispin,jspin)==1)WRITE(8,*)ispin,jspin
                  IF(spin(ispin,jspin)==-1)WRITE(9,*)ispin,jspin
               END DO
            END DO
            WRITE(8,*)""
            WRITE(8,*)""
            WRITE(9,*)""
            WRITE(9,*)""         
            CLOSE(8)
            CLOSE(9)
         END IF
      END DO
   
      IF(PBC) THEN
         OPEN(unit=13,file="cdev_pbc.dat", position="append")
      ELSE
         OPEN(unit=13,file="cdev_obc.dat", position="append")         
      END IF
   
      IF(N_T /= 1) THEN
         DO i = 1, N_T
           ! dT = (T_f - T_i)/(N_T-1)
            IF(i==1) THEN
               c_dev = (E_ave(i+1) - E_ave(i))/dT
            ELSE IF(i==N_T) THEN
               c_dev = (E_ave(i) - E_ave(i-1))/dT
            ELSE
               c_dev = (E_ave(i+1) - E_ave(i-1))/(2*dT)
            END IF
            T   = T_i + (i-1)*(T_f - T_i)/(N_T-1)
            WRITE(13,*)T,c_dev
         END DO
      END IF
      
   WRITE(46,*)L, E_ave(1)
   DEALLOCATE(spin, seed, box)
   END DO

   CLOSE(13)
   CLOSE(46)
END PROGRAM 2D_ising
