PROGRAM lattice

    IMPLICIT NONE

    LOGICAL, ALLOCATABLE          :: lattice(:,:) 
    INTEGER, ALLOCATABLE          :: x(:),y(:)
    DOUBLE PRECISION, ALLOCATABLE :: dx(:),dy(:)

    INTEGER  :: Nsteps,Np,L
    INTEGER  :: istep,isubstep, irun  
    INTEGER  :: dir,i,j,nfail,njumps
    INTEGER  :: k, dk
    

    INTEGER, PARAMETER                 :: MAXINT = 30                               !!VARIABLES FOR COUNTING
    INTEGER, PARAMETER                 :: RUNNER = 5, k_max=4
    INTEGER, PARAMETER                 :: N_i = 200, N_f = 1750, N_max = 40
    DOUBLE PRECISION, PARAMETER        :: b = 2
    LOGICAL, PARAMETER                 :: log_plot = .FALSE.
!!ALLOWED DIRECTIONS 
    INTEGER                            :: free(4),nfree, i_N   
    INTEGER                            :: dxtrial(4),dytrial(4) 
    INTEGER                            :: xnew(4),ynew(4)

    REAL, DIMENSION(2)                 :: rnd(2)
    REAL                               :: rnd1
    DOUBLE PRECISION                   :: dxsum,dysum,dxsqsum,dysqsum 
    DOUBLE PRECISION                   :: t,deltat,drsqave,D,a,help
    DOUBLE PRECISION                   :: Dsum, D2sum, Dtot=0.0d0, errtot=0.0d0, errtot_T = 0.0d0,old_step = 1.0d0
    DOUBLE PRECISION, DIMENSION(k_max) :: Avg, Avg2    


    INTEGER                                     :: sizer, values(1:8)               !!SEED VARIABLES
    INTEGER, DIMENSION(:), ALLOCATABLE          :: seed, box                        !!SEED VARIABLES


!!SET AVERAGE TIME BETWEEN JUMPS AND JUMP LENGTH UNITS IS SEC AND CM
!!ALTHOUGH ACTUALLY IS NOT NEEDED FOR THE SIMULATION
    deltat = 1.d0 !!OR 1e-9 FOR 1 ns 
    a      = 1.d0 !!OR 2e-8 FOR 2 Ang

    PRINT*," Nsteps>"
    READ*, Nsteps  
    !PRINT*," Number of particles>"
    !READ*, Np 
    PRINT*," L>"
    READ*, L
    
    dk = Nsteps/dfloat(k_max) !!Length of intervals for the D block average
    
    OPEN(unit=1,file="disp_data.dat", status="replace", action="write")
    WRITE(1,*)'#EACH INDEX IS A DIFFERENT RUN WITH THE SAME SETTINGS'
    WRITE(1,*)'#TOT step:',Nsteps,'Number of particles:',Np,' dimension lattice:',L,'x',L
    WRITE(1,*)'#density:', real(Np)/L**2
    WRITE(1,*)'#t, <R^2(t)>, D(t), <D>_t, var, block average D, var'

    OPEN(unit=3,file="seed.dat",status="replace",action="write")
    WRITE(3,*)'#EACH INDEX IS A DIFFERENT RUN'

    OPEN(unit=5,file="disp_data_tot.dat", status="replace", action="write")
    WRITE(5,*)'#EACH INDEX IS A DIFFERENT RUN WITH THE SAME SETTINGS, HERE WE TAKE ALL MC STEPS'
    WRITE(5,*)'#TOT step:',Nsteps,'Number of particles:',Np,' dimension lattice:',L,'x',L
    WRITE(5,*)'#density:', real(Np)/L**2
    WRITE(5,*)'#t, <R^2(t)>, D(t), <D>_t, var, block average D, var'

    OPEN(unit=46,file="D.dat",status="replace",action="write")
    WRITE(46,*)'#rho, D, var'



    DO i_N = 1, N_max
        Np = N_i + dfloat(i_N-1)*(N_f - N_i)/dfloat(N_max-1)
    
        Dtot     = 0.0d0
        errtot   = 0.0d0
        errtot_T = 0.0d0
        DO irun = 1, RUNNER
    !!SEED GENERATION AND SAVING
            CALL random_seed(sizer)
            ALLOCATE(seed(sizer))
            ALLOCATE(box(sizer))       
            CALL date_and_time(values = values)
            seed(:) = values(8)*(3.4 + 2.**irun)
            CALL random_seed(put=seed)
            CALL random_seed(get=box)
            DO i = 1, sizer
                WRITE(3,*)box(i)
            END DO
    !!SAVED
            PRINT*,'Doing lattice gas walk to',Nsteps,'MC steps'
            PRINT*,'using',Np,' particles  on a',L,'x',L,'square lattice'
            IF(Np >= L*L) THEN
                PRINT*,'Number of particles > number of sites !!!' 
                STOP 'Too small lattice' 
            END IF
        
            ALLOCATE(lattice(0:L-1,0:L-1))
            ALLOCATE(x(Np),y(Np)) 
            ALLOCATE(dx(Np),dy(Np))

        !!MARK ALL POSITIONS AS EMPTY 
            DO i = 0, L-1
                DO j = 0, L-1
                    lattice(i,j) = .FALSE. 
                END DO
            END DO

        !!ENUMERATION OF DIRECTIONS: 1->L 2->R 3->U 4->D 
            dxtrial(1) = +1; dytrial(1)= 0;   
            dxtrial(2) = -1; dytrial(2)= 0;   
            dxtrial(3) =  0; dytrial(3)=+1; 
            dxtrial(4) =  0; dytrial(4)=-1;

            nfail=0; njumps=0; 
        !!GENERATE PARTICLES ON LATTICE
            DO i=1,Np
                DO    !!LOOP UNTIL EMPTY POSITION FOUNDED
                      !!TO BE ON SAFE SIDE, CHECK THAT UPPER LIMIT NOT REACHED
                    CALL random_number(rnd)
                    x(i) = INT(rnd(1)*L);  IF(x(i)>=L) x(i) = L-1;  
                    y(i) = INT(rnd(2)*L);  IF(y(i)>=L) y(i) = L-1; 
                    IF(lattice(x(i),y(i))) THEN
                !!POSITION ALREADY FILLED, LOOP TO FIND NEW TRIAL 
                        CYCLE 
                    ELSE
                    lattice(x(i),y(i)) = .TRUE. 
                !!SUCCESS, GO TO NEXT PARTICLE
                        EXIT 
                    END IF
                END DO
                dx(i) = 0.0d0; dy(i) = 0.0d0; 
            END DO




            t        = 0.0
            Dsum     = 0.0d0
            D2sum    = 0.0d0
            Avg      = 0.0d0
            Avg2     = 0.0d0
            old_step = 1.0d0
            DO istep = 0, Nsteps-1 !!LOOP OVER MC STEPS
                DO isubstep = 1, Np !!DO ALL PARTICLES ON AVERAGE ONCE EVERY MC STEP
            !!PICK ONE PARTICLE RANDOM 
                    CALL random_number(rnd1)
                    i = INT(rnd1*Np) + 1; IF(i>Np) i = Np;

            !!FIND POSSIBLE DIRECTIONS, STORE IT IN FREE()
                    nfree=0 
                    DO j = 1, 4
                        xnew(j) = x(i) + dxtrial(j); 
                        IF(xnew(j) >= L) xnew(j) = 0; IF(xnew(j)<0) xnew(j) = L-1; 
                        ynew(j) = y(i) + dytrial(j); 
                        IF(ynew(j) >= L) ynew(j) = 0; IF(ynew(j)<0) ynew(j) = L-1; 
                        IF(.NOT. lattice(xnew(j),ynew(j))) THEN
                !!SUCCESS:POSITION FREE 
                            nfree = nfree + 1 
                            free(nfree) = j 
                        END IF
                    END DO

            !!IF NO POSSIBLE DIRECTIONS, GET NEW PARTICLE
                    IF(nfree == 0) THEN
                        nfail = nfail + 1 
                        CYCLE 
                    END IF
                    njumps = njumps + 1

            !!PICK ONE OF THE POSSIBLE DIRECTIONS RANDOMLY  
            !!NOTE THAT THE dir > nfree CHECK HERE IS REALLY NEEDED  
                    CALL random_number(rnd1)
                    dir = INT(rnd1*nfree) + 1; IF(dir>nfree) dir = nfree 
                    j   = free(dir)

            !!NOW x(i), y(i) IS OLD POSITION AND xnew(j), ynew(j) THE NEW ONES 
            !!DOUBLE CHECK THAT NEW SITE REALLY IS FREE 
                    IF(lattice(xnew(j),ynew(j))) THEN
                        PRINT*,'ERROR:    THIS   SHOULD   BE   IMPOSSIBLE'   
                        PRINT*,i,j,dir,nfree  
                        PRINT*,free  
                        PRINT*,x(i),y(i),xnew(j),ynew(j) 
                        STOP 'ERROR  new  site  bug'  
                    END IF
            !!EMPTY OLD POSITION AND FILL NEW
                    lattice(x(i),y(i))       = .FALSE. 
                    lattice(xnew(j),ynew(j)) = .TRUE.

                    x(i)  = xnew(j); y(i) = ynew(j);          
                    dx(i) = dx(i) + dxtrial(j); dy(i) = dy(i) + dytrial(j); 
                END DO




            !!CALCULATE AND PRINT INTERMEDIATE RESULTS
            !!GET TOTAL DISPLACEMENT FROM dx,dy 
                dxsum   = 0.0d0; dysum   = 0.0d0; 
                dxsqsum = 0.0d0; dysqsum = 0.0d0; 
                DO i = 1, Np
                    dxsum   = dxsum + dx(i);   dysum = dysum + dy(i);   
                    dxsqsum = dxsqsum+dx(i)*dx(i);
                    dysqsum = dysqsum+dy(i)*dy(i); 
                END DO
                    drsqave = (dxsqsum+dysqsum)/(1.0*Np)

                IF(t > 0.0) THEN
            !!GET DIFFUSION COEFFICIENT BY PROPER SCALING  
                    D       = drsqave*a*a/(4*t)
                    Dsum    = Dsum + D
                    D2sum   = D2sum + D*D
                    k       = int(istep/dk) + 1
                    Avg(k)  = Avg(k) + D
                    Avg2(k) = Avg2(k) + D*D
                    WRITE(5,*)t, drsqave, D, Dsum/(istep+1), D2sum/(istep+1) - (Dsum/(istep+1))**2,&
                                Avg(k)/((istep+1) -(k-1)*dk), Avg2(k)/((istep+1) -(k-1)*dk) - (Avg(k)/((istep+1) -(k-1)*dk))**2
                            
                    IF(log_plot) THEN
                        IF(dfloat(istep)/old_step == b) THEN
                            WRITE(1,*)t, drsqave, D, Dsum/(istep+1), D2sum/(istep+1) - (Dsum/(istep+1))**2,&
                                Avg(k)/((istep+1) -(k-1)*dk), Avg2(k)/((istep+1) -(k-1)*dk) - (Avg(k)/((istep+1) -(k-1)*dk))**2
                            old_step = istep
                        END IF    
                    ELSE
                        IF(MOD(istep, MAXINT) == 0)WRITE(1,*)t, drsqave, D, Dsum/(istep+1), D2sum/(istep+1) - (Dsum/(istep+1))**2,&
                                Avg(k)/((istep+1) -(k-1)*dk), Avg2(k)/((istep+1) -(k-1)*dk) - (Avg(k)/((istep+1) -(k-1)*dk))**2
                    END IF
                    
                    END IF


                t = t + deltat
            END DO
    
        !!GET TOTAL DISPLACEMENT FROM dx,dy
            dxsum = 0.0d0; dysum = 0.0d0; 
            dxsqsum = 0.0d0; dysqsum = 0.0d0; 
            DO i = 1, Np
                dxsum   = dxsum + dx(i);         dysum   = dysum + dy(i);   
                dxsqsum = dxsqsum + dx(i)*dx(i); dysqsum = dysqsum + dy(i)*dy(i);     
            END DO
            PRINT*,'dxsum',dxsum,'   dysum',dysum 
            PRINT*,'dxsqsum',dxsqsum,' dysqsum',dysqsum
        
            drsqave = (dxsqsum+dysqsum)/dfloat(Np)
        
            PRINT*,'drsqave',drsqave 
            PRINT*,'Number of  failed jumps',nfail,' number of  successes',njumps 
        !!GET DIFFUSION COEFFICIENT BY PROPER SCALING
            D = drsqave*a*a/(4*t) 
            Dsum    = Dsum + D
            D2sum   = D2sum + D*D
            k       = int(Nsteps/dk)
            Avg(k)  = Avg(k) + D
            Avg2(k) = Avg2(k) + D*D
            WRITE(5,*)t, drsqave, D, Dsum/(Nsteps+1), D2sum/(Nsteps+1) - (Dsum/(Nsteps+1))**2,&
                                Avg(k)/(Nsteps+1-(k-1)*dk), Avg2(k)/(Nsteps+1-(k-1)*dk) - (Avg(k)/(Nsteps+1-(k-1)*dk))**2
            IF(MOD(Nsteps, MAXINT) == 0)WRITE(1,*)t, drsqave, D, Dsum/(Nsteps+1), D2sum/(Nsteps+1) - (Dsum/(Nsteps+1))**2,&
                                Avg(k)/(Nsteps+1-(k-1)*dk), Avg2(k)/(Nsteps+1-(k-1)*dk) - (Avg(k)/(Nsteps+1-(k-1)*dk))**2
            PRINT*,'RUNNER:',irun
            PRINT*,'Density Np/L^2=',real(Np)/L**2,' t=',t,'; <D>_T=',Dsum/Nsteps,'\pm',(D2sum/Nsteps - (Dsum/Nsteps)**2)**0.5,&
            '; Block average D=',SUM(Avg)/Nsteps,'\pm',((dfloat(dk)/Nsteps)*SUM(Avg2/dk - (Avg/dk)**2))**0.5
        
            Dtot     = Dtot + SUM(Avg)/Nsteps
            errtot   = errtot + ((dfloat(dk)/Nsteps)*SUM(Avg2/dk - (Avg/dk)**2))**0.5
            errtot_T = errtot_T + (D2sum/Nsteps - (Dsum/Nsteps)**2)**0.5
            DEALLOCATE(seed, box, lattice, x, y, dx, dy)
            WRITE(1,*)""
            WRITE(1,*)""
            WRITE(3,*)""
            WRITE(3,*)""
            WRITE(5,*)""
            WRITE(5,*)""
        END DO

        PRINT*,'--------------------------------------------------------------------------------------------------------'
        PRINT*,'Density Np/L^2=',real(Np)/L**2,' t=',t,'; <D>_T=',Dtot/RUNNER,'\pm',errtot_T/RUNNER,&
               '; Block average D=',Dtot/RUNNER,'\pm',errtot/RUNNER
        PRINT*,'--------------------------------------------------------------------------------------------------------'
        WRITE(46,*)dfloat(Np)/(L**2), Dtot/RUNNER, errtot/RUNNER
    END DO
    
        CLOSE(1)
        CLOSE(3)
        CLOSE(5)
        PRINT*,'seed data stored in "seed.dat"'
        PRINT*,'Displacement data for SOME steps stored in "disp_data.dat"'
        PRINT*,'Displacement data for ALL steps stored in "disp_data_tot.dat"'
    
END PROGRAM lattice
