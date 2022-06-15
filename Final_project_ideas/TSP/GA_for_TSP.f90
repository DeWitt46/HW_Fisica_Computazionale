!!---------------------------------------------------------------------------------------------------
!!Final project for Master course: Laboratorio di fisica computazionale, UniTS 2021/22
!!
!!Code using genetic algorithm for TSP: we are interested in crossover operators comparison
!!
!!Also scaling time comparison with brute-force method
!!
!!Reference: https://doi.org/10.1155/2017/7430125 but elitism selection strategy implemented.
!!---------------------------------------------------------------------------------------------------
MODULE common
    IMPLICIT NONE
    !!---------------------------------------------------------------------------------
    !! RUN PARAMETERS
    INTEGER, PARAMETER, PUBLIC :: nkind = 16, N_sim = 30
    LOGICAL, PUBLIC :: RND_PLACE_FLAG=.FALSE. !! CHECK N IF THIS IS ENABLED
    LOGICAL, PUBLIC :: BF_FLAG = .TRUE. !! Automatically shuts down for N geq 10
    
    INTEGER, PUBLIC                                     :: sizer, values(1:8)                !Seed variables
    INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE          :: seed, box                         !Seed variables
    !!---------------------------------------------------------------------------------
    !!---------------------------------------------------------------------------------
    !! GA PARAMETERS
    INTEGER, PUBLIC :: N_gen = 30 !! Maximum number of generations
    INTEGER, PUBLIC :: N_pop = 30 !! Number of inidividuals in a certain generation, NOTE: must be even
    REAL(KIND=nkind), PUBLIC :: p_cross = 0.95 !! Crossover probability
    REAL(KIND=nkind), PUBLIC :: p_mut = 0.1 !! Mutation probability
    !!---------------------------------------------------------------------------------
    !!---------------------------------------------------------------------------------
    !! TSP PARAMETERS
    INTEGER, PUBLIC :: N = 7  !! Number of cities
    INTEGER, PUBLIC :: L = 30 !! Dimension of squared lattice
    !!---------------------------------------------------------------------------------
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: place, pos !! Matrix representing the place cities are
                                                                !! and their coordinates (x,y)
    
    REAL(KIND=nkind), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: distances !! Matrix of distances between each pair of cities
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: parents, offspring !! Array of parents and offspring
                                                               !! a row is a sequence, columns are list of indeces
    
    REAL(KIND=nkind), DIMENSION(:), ALLOCATABLE, PUBLIC :: fit_p, fit_o, dummy_fit_p, dummy_fit_o !! Arrays of fitnesses of the parents and offspring
    
    INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: best_GA_seq !! Array representing best GA sequence among all simulations
    
    
    CONTAINS
    
    
    
    !!Maybe useful later
    SUBROUTINE city_gen(L, x, y)
    !!-----------------------------------------------------------------
    !! Generate a city position given dimension of squared lattice
    !!-----------------------------------------------------------------
        INTEGER, INTENT(IN) :: L
        INTEGER, INTENT(OUT) :: x, y
        
        REAL(KIND=nkind), DIMENSION(2) :: rnd
        
        
        
        CALL random_number(rnd)
        x = int(rnd(1)*L)
        y = int(rnd(2)*L)
        IF(x==0)x=1
        IF(y==0)y=1
        
    END SUBROUTINE city_gen



    SUBROUTINE INITIAL()
    !!----------------------------------------------
    !! Initial settings for the run
    !!----------------------------------------------
        INTEGER :: i
        
        
        
        !PRINT*, "Number of cities N ="
        !READ*,N
        !PRINT*, "Linear dimension of lattice L ="
        !READ*,L
        !IF(L<N*N)PRINT*,"INCOMPATIBLE DIMENSIONS"
    
        IF(N >= 10)BF_FLAG=.FALSE.
        !! -----------------------------------
        !! Seed generation and saving
        !! -----------------------------------
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
        !! Generated and saved
        !! --------------------------------------
        
        !! ------------------------------------------------------------------
        !! PLACE SETTINGS: YOU CAN GENERATE IT OR MAKE THE SCRIPT READ ONE
        !! JUST CHANGE THIS FLAG AND CHECK THAT DIMENSIONS ARE PROPERLY SET
        !! ------------------------------------------------------------------
        
        ALLOCATE(pos(N,2),distances(N,N))
        ALLOCATE(parents(N_pop,N),offspring(N_pop,N))
        ALLOCATE(dummy_fit_p(N_pop),dummy_fit_o(N_pop))
        ALLOCATE(fit_p(N_pop), fit_o(N_pop))
        
        
        
        IF(RND_PLACE_FLAG) THEN !! Generate place
            PRINT*,"GENERATING RANDOM CONFIGURATION"
            PRINT*,"IF YOU WANTED AN EXISTING ONE PLEASE MODIFY"
            PRINT*,"RND_PLACE_FLAG INSIDE THE CODE AND RE-COMPILE"

            CALL place_gen(N,L,place)
            
            CALL positions(place, pos) !! write cities position
        
            !! If you want to see (x,y) plane you have to rotate the matrix
            !! In the notation used 1 argument is row and 2 columns
            !! If you want that 1 stands for x and 2 for y, you have to rotate
            OPEN(unit=3,file="place.dat", status="replace", action="write")
            WRITE(3,*)'#PLACE OF CURRENT RUN'
            WRITE(3,*)'#N =',N,"L=",L
            DO i = L, 1, -1
                WRITE(3,*)place(:,i)
            END DO
            CLOSE(3)

            OPEN(unit=1,file="pos.dat", status="replace", action="write")
            WRITE(1,*)'#CITY POSITIONS OF CURRENT RUN'
            WRITE(1,*)'#N =',N,"L=",L
            DO i = 1, N
                WRITE(1,*)pos(i,:)
                WRITE(1,*)"" !! For GNUplot index
                WRITE(1,*)""
            END DO
            CLOSE(1)

        ELSE !! TAKE EXISTING PLACE FROM A FILE OR KEY
    
            PRINT*,"YOU ONLY NEED DISTANCES MATRIX"
            
        END IF
        


        !!--------------------------------------------------------------------------
        !!Open and reset data files
        !!--------------------------------------------------------------------------
        OPEN(unit=1,file="data_gen.dat", status="replace", action="write")
        WRITE(1,*)'#FITNESS DATA OF CURRENT RUN PER GENERATION'
        WRITE(1,*)'#N=',N,"L=",L, "N_sim=",N_sim
        WRITE(1,*)"#EACH INDEX IS A SIMULATION"
        WRITE(1,*)"#gen, avg_fit, dev, best_fit, best_sequence"
        CLOSE(1)
        
        OPEN(unit=1,file="data_sim.dat", status="replace", action="write")
        WRITE(1,*)'#FITNESS DATA OF CURRENT RUN'
        WRITE(1,*)'#N=',N,"L=",L, "N_sim=",N_sim
        WRITE(1,*)"#sim, avg_fit, dev, best_fit, best_sequence"
        CLOSE(1)
        
    END SUBROUTINE



    SUBROUTINE place_gen(N, L, place)
    !!-----------------------------------------------------------------
    !! Generate a place given dimension of squared lattice and 
    !! number of cities
    !!-----------------------------------------------------------------
        INTEGER, INTENT(IN) :: N, L
        INTEGER             :: i, x, y
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: place
        REAL(KIND=nkind), DIMENSION(2) :: rnd
        ALLOCATE(place(L,L))
 
        place = 0
        DO i = 1, N
            DO
                CALL random_number(rnd)
                x = int(rnd(1)*L)
                y = int(rnd(2)*L)
                IF(x==0)x=1
                IF(y==0)y=1
                
                IF(place(x,y) == 0) THEN
                    place(x,y) = i
                    EXIT
                END IF    
            END DO        
        END DO
    END SUBROUTINE place_gen

    
    
    SUBROUTINE positions(place, pos)
    !!------------------------------------------------
    !! Gives array of positions of the cities given the 
    !! place where they are located
    !!------------------------------------------------
        INTEGER, DIMENSION(L,L), INTENT(IN)  :: place
        INTEGER, DIMENSION(N,2), INTENT(OUT) :: pos
        
        INTEGER :: i
        
        
        
        IF (MAXVAL(place) /= SIZE(pos,1))PRINT*,"NUMBER CITIES NOT CONSISTENT" !Check number of cities
                
        DO i = 1, SIZE(pos, 1)
            pos(i,:) = FINDLOC(place, i)
        END DO
        
    END SUBROUTINE positions 
    
    
    
    SUBROUTINE distance(x1, y1, x2, y2, d)
    !!-------------------------------------------------------
    !! Compute distance given the position of two cities
    !!-------------------------------------------------------
        INTEGER, INTENT(IN) :: x1, y1, x2, y2
        REAL(KIND = nkind), INTENT(OUT) :: d
        REAL(KIND=nkind) :: d_x, d_y
        
        d_x = x2 - x1
        d_y = y2 - y1
        d = sqrt(d_x*d_x + d_y*d_y)
    END SUBROUTINE distance    



    SUBROUTINE fit(sequence, distances, fitness)
    !!-------------------------------------------------------
    !! Compute fitness of a sequence given
    !! the place and distances matrix
    !!-------------------------------------------------------
        INTEGER, DIMENSION(N), INTENT(IN) :: sequence
        REAL(KIND=nkind), DIMENSION(N,N), INTENT(IN) :: distances
        REAL(KIND=nkind), INTENT(OUT) :: fitness
        
        REAL(KIND=nkind) :: d_total
        INTEGER :: i, j
        INTEGER, DIMENSION(2) :: pos1, pos2
        INTEGER, DIMENSION(SIZE(sequence)+1) :: dummy_sequence
        
        
        
        DO i = 1, N 
            dummy_sequence(i) = sequence(i)
        END DO
        dummy_sequence(N+1) = sequence(1)
        
        
        j = 0
        DO i = 1, SIZE(sequence)
            j = i + 1
            d_total = d_total + distances(dummy_sequence(i), dummy_sequence(j))
        END DO
        fitness = (L*L)/d_total
    END SUBROUTINE
    
    
    
    SUBROUTINE seq_gen(N, sequence)
    !!------------------------------------------------------
    !! Generates a random sequence given N number of cities
    !!------------------------------------------------------
        INTEGER, INTENT(IN) :: N
        INTEGER, DIMENSION(N), INTENT(OUT) :: sequence
        
        INTEGER :: i, j, candidate
        REAL(KIND=nkind) :: rnd
        
        
        sequence = 0
        DO i = 1, N
            DO
                CALL random_number(rnd)
                candidate = 1 + FLOOR(N*rnd)
                IF(ANY(candidate == sequence) .EQV. .FALSE.) EXIT
            END DO
            sequence(i) = candidate
        END DO                
    END SUBROUTINE    
    
    
    
    SUBROUTINE roulette(fit, ind)
    !!----------------------------------------------------------
    !! Roulette wheel method to choose parents to breed
    !! Takes fitness array and gives the indices of the parents
    !!----------------------------------------------------------
        REAL(KIND=nkind), DIMENSION(N_pop), INTENT(IN) :: fit
        INTEGER, DIMENSION(2), INTENT(OUT) :: ind

        LOGICAL :: ind_taken
        INTEGER :: i, j
        REAL(KIND=nkind) :: rnd
        REAL(KIND=nkind), DIMENSION(N_pop) :: dummy_fit
        REAL(KIND=nkind), DIMENSION(SIZE(fit) + 1) :: wheel

        !! Take the probabilities proportional to fitness
        !! Fix the slices of the wheel according to probabilities        

        dummy_fit = 0
        wheel(1) = 0
        DO i = 1, N_pop
        !! Normalize fitness to probabilites
            dummy_fit(i) = fit(i)/SUM(fit)
            wheel(i+1) = wheel(i) + dummy_fit(i)
        END DO
        
        !! See in which slice of the wheel the ball stops
        ind = 0
        DO i = 1, 2
            DO
                ind_taken = .FALSE.
                CALL random_number(rnd)
                DO j = 1, N_pop-1
                    IF(rnd >= wheel(j)) THEN
                        IF(rnd < wheel(j+1)) THEN
                            ind(i) = j
                            EXIT
                        END IF
                    END IF
                END DO
                IF(ind(i) == 0)ind(i)=N_pop  !! If it didn't stop at the last but one, then It's the last
                IF(ind(2) == ind(1)) ind_taken=.TRUE.
                IF(ind_taken .EQV. .FALSE.) EXIT
            END DO
        END DO

    END SUBROUTINE
    
    
    
    SUBROUTINE MUT(old_seq, new_seq)
    !!-------------------------------------------------------------------
    !! Takes a sequence and mute It into a new one
    !! Basic two-"bits" swap
    !! Should be the direct transpose of one-bit flip for this problem
    !!-------------------------------------------------------------------
        INTEGER, DIMENSION(N), INTENT(IN) :: old_seq
        INTEGER, DIMENSION(SIZE(old_seq)), INTENT(OUT) :: new_seq    

        LOGICAL :: bit_taken
        REAL(KIND=nkind) :: rnd
        INTEGER :: bit1, bit2, pos_bit1, pos_bit2
        
        
        
        new_seq = old_seq
        !! Random search a "bit" to swap
        CALL random_number(rnd)
        bit1 = 1 + FLOOR(rnd*N)
        pos_bit1 = FINDLOC(old_seq, bit1, 1)
        !! Random search of the second "bit" to swap with
        DO
            bit_taken = .FALSE.
            CALL random_number(rnd)
            bit2 = 1 + FLOOR(rnd*N)
            IF(bit2 == bit1)bit_taken = .TRUE.
            IF(bit_taken .EQV. .FALSE.)EXIT
        END DO    
        pos_bit2 = FINDLOC(old_seq, bit2, 1)

        new_seq(pos_bit1) = old_seq(pos_bit2)
        new_seq(pos_bit2) = old_seq(pos_bit1)

        !PRINT*,"MUT CHECK" !! Debug
        !CALL CHECK(new_seq)
    END SUBROUTINE
    
    
    
    SUBROUTINE CX_CYCLE(par1, par2, off1, off2)
    !!------------------------------------------------------
    !! Main process of the CX2 crossover operator:
    !! Search for a candidate in par1-par2 correspondence.
    !! Returns offsprings with only no-conflict "bits" filled.
    !! The remaining ones are 0.
    !!------------------------------------------------------
        INTEGER, DIMENSION(N), INTENT(IN) :: par1, par2
        INTEGER, DIMENSION(N), INTENT(OUT) :: off1, off2

        INTEGER :: position, candidate
        INTEGER :: i, iterator1, iterator2

        
        
        off1 = 0
        off2 = 0
        
        iterator1 = 1
        iterator2 = 1
        position  = 1
        candidate = par2(position)
        DO i = 1, N  !! Max number of tries
            IF(ANY(candidate == off1) .EQV. .FALSE.) THEN
                off1(iterator1) = candidate
                iterator1 = iterator1 + 1
            END IF    
            !! 2 moves for true off2 candidate
            position = FINDLOC(par1,candidate,1)
            candidate = par2(position)
            position = FINDLOC(par1,candidate,1)
            candidate = par2(position)
            !! 2 moves done
            IF(ANY(candidate == off2) .EQV. .FALSE.) THEN
                off2(iterator2) = candidate
                iterator2 = iterator2 + 1
            END IF
            !! 1 move for off1 candidate
            position = FINDLOC(par1,candidate,1)
            candidate = par2(position)
            !! 1 move done
       END DO                                 


    END SUBROUTINE



    SUBROUTINE CROSS(par1, par2, off1, off2, op)
    !!-------------------------------------------------------
    !! Cross two parents to have an offspring;
    !! Three crossover operators are implemented and can be
    !! chosen: PMX, OX, CX2
    !!-------------------------------------------------------
        INTEGER, INTENT(IN) :: op
        INTEGER, DIMENSION(N), INTENT(IN) :: par1, par2
        INTEGER, DIMENSION(N), INTENT(OUT) :: off1, off2
 
        REAL(KIND=nkind) :: rnd
        INTEGER :: cut_start, cut_end, cut_dummy, PROVA
        INTEGER :: i, j, k, candidate, iterator1, iterator2, position
        INTEGER, DIMENSION(2) :: pos
        INTEGER, DIMENSION(SIZE(par1)) :: dummy_par1, dummy_par2, dummy_off1, dummy_off2
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: MAP
        
        
        
        IF(op == 1) THEN !! Partially Mixed crossover operator
        !!----------------------------------------------------------------------------------
        !! Choose two random cut points: the portion of both sequences between these
        !! are transposed in the other parent's offspring, the rest is exchanged.
        !! The map is between each "bit" of the two portions. The "bits" not in this
        !! map are filled like in the parent's string, the others are filled
        !! according to the map.
        !!----------------------------------------------------------------------------------
            !! Choose two random cut points
            CALL random_number(rnd)
            cut_start = 1 + FLOOR(rnd*N)
            CALL random_number(rnd)
            cut_end = 1 + FLOOR(rnd*N) !! Portion can even be a single "bit"
            IF(cut_end < cut_start) THEN !! Order start and end of the portion
                cut_dummy = cut_end
                cut_end = cut_start
                cut_start = cut_dummy
            END IF    

            !! Transpose portion in offspring
            off1 = 0
            off2 = 0
            DO i = cut_start, cut_end
                off1(i) = par2(i)
                off2(i) = par1(i)
            END DO
        
            !! Fill the "bits" with no conflicts
            DO i = 1, N
                IF(off1(i) == 0) THEN
                    IF(ANY(par1(i) == off1) .EQV. .FALSE.)off1(i)=par1(i) !! Fill if no conflicts are found
                END IF    
                IF(off2(i) == 0) THEN
                    IF(ANY(par2(i) == off2) .EQV. .FALSE.)off2(i)=par2(i) !! Fill if no conflicts are found
                END IF
            END DO     
                
            !! Exchange remaining "bits" informations
            ALLOCATE(MAP((1+cut_end-cut_start),2)) !! MAP stores PORTION information like this
                                                   !! (parent 1 first bit, parent 2 first bit) ...
                                                   !! (parent 1 last bit, parent 2 last bit)
                                                   
            DO i = 1, (cut_end - cut_start+1)      !! MAP starts from 1, but the portion not necessarily
                MAP(i,1) = par1(cut_start + i -1)
                MAP(i,2) = par2(cut_start + i -1)
            END DO    
            
            DO i = 1, N
                IF(off1(i) == 0) THEN
                    candidate = par1(i)
                    DO j = 1, (cut_end-cut_start+1) !! Max time one can fail in finding no conflicts "bit"
                        pos = FINDLOC(MAP,candidate,1) !! Find the location in maap of the value P1(i)
                                                        !! if It's the first try, otherwise of the candidate
                        candidate = MAP(pos(2),1)      !! Candidate is found using portion MAP
                        IF(ANY(candidate == off1) .EQV. .FALSE.) EXIT !!Check If candidate is already in offspring
                    END DO    
                    off1(i) = candidate
                END IF
                
                IF(off2(i) == 0) THEN !! Same but with offspring of the parents 2
                    candidate = par2(i)
                    DO j = 1, (cut_end-cut_start+1)
                        pos = FINDLOC(MAP,candidate,1)
                        candidate = MAP(pos(1),2)
                        IF(ANY(candidate == off2) .EQV. .FALSE.) EXIT
                    END DO    
                    off2(i) = candidate
                END IF
            END DO
            
            DEALLOCATE(MAP)
            
        ELSE IF(op == 2) THEN !! Order crossover operator
        !!-----------------------------------------------------------------------------
        !! Choose two random cut points: the portion of both sequences between these
        !! remains at the same place in the same order in his own offspring.
        !! Remaining "bits" starting from the last cut are transposed from a parent
        !! into the other one's offspring preserving the order; already existing "bits"
        !! are discarded.
        !!-----------------------------------------------------------------------------
            !! Choose two random cut points
            CALL random_number(rnd)
            cut_start = 1 + FLOOR(rnd*N)
            CALL random_number(rnd)
            cut_end = 1 + FLOOR(rnd*N) !! Portion can even be a single "bit"
            IF(cut_end < cut_start) THEN !! Order start and end of the portion
                cut_dummy = cut_end
                cut_end = cut_start
                cut_start = cut_dummy
            END IF 

            !! Transpose portion in offspring
            off1 = 0
            off2 = 0
            DO i = cut_start, cut_end
                off1(i) = par1(i)
                off2(i) = par2(i)
            END DO
            
            !! Exchange remaining "bits" informations
            ALLOCATE(MAP(2,(N-(cut_end-cut_start+1)))) !! MAP stores REMAINING information like this:
                                                        !! (parent 1 first bit, parent 1 second bit,...)
                                                       !! (parent 2 first bit, parent 2 second bit,...)
            iterator1 = 1
            iterator2 = 1
            DO i = cut_end+1, N      !! Divide in two cycles to preserve the order
                candidate = par1(i)
                IF(ANY(candidate == off2) .EQV. .FALSE.) THEN !! Take candidate only if no conflicts are present
                    MAP(1,iterator1) = candidate
                    iterator1 = iterator1 +1
                END IF
                candidate = par2(i)
                IF(ANY(candidate == off1) .EQV. .FALSE.) THEN
                    MAP(2,iterator2) = candidate
                    iterator2 = iterator2 +1
                END IF
            END DO
            DO i = 1, cut_end      !! Divide in two cycles to preserve the order
                candidate = par1(i)
                IF(ANY(candidate == off2) .EQV. .FALSE.) THEN !! Take candidate only If no conflicts are present
                    MAP(1,iterator1) = candidate
                    iterator1 = iterator1 +1
                END IF
                candidate = par2(i)
                IF(ANY(candidate == off1) .EQV. .FALSE.) THEN
                    MAP(2,iterator2) = candidate
                    iterator2 = iterator2 +1
                END IF
            END DO
        
            DO i = 1, (N - cut_end)
                off1(cut_end+i) = MAP(2,i)
                off2(cut_end+i) = MAP(1,i)
            END DO    
            DO i = N-cut_end+1, (N - (cut_end - cut_start + 1))
                off1(cut_end-N+i) = MAP(2,i)
                off2(cut_end-N+i) = MAP(1,i)
            END DO
            
            DEALLOCATE(MAP)
            
        ELSE IF(op == 3) THEN !! Cycle crossover operator 2 (Their proposal)
        !!--------------------------------------------------------------------------------
        !! Choose as first offspring first "bit", the one of the parent 2.
        !! We call "move" that search process of OX; for second offspring you do
        !! two moves, for next first offspring "bits" you do one.
        !! If all the "bits" are filled the job is done;
        !! If some "bits" are in conflict, you discard all of these from both parents
        !! and re-start the process from the first no-conflict "bit" of the second parent
        !!--------------------------------------------------------------------------------
        
            off1 = 0
            off2 = 0
            dummy_off1 = 0
            dummy_off2 = 0
            dummy_par1 = par1  !! First process with full parents
            dummy_par2 = par2
            
            iterator1 = 1
            iterator2 = 1
            DO PROVA = 1, N
                CALL CX_CYCLE(dummy_par1, dummy_par2, dummy_off1, dummy_off2)
                DO i = 1, N
                    IF(dummy_off1(i) /= 0) THEN !! "bits" zero must not overlap
                        candidate = dummy_off1(i)
                        IF(ANY(candidate == off1) .EQV. .FALSE.) THEN !! Check for conflicts
                            off1(iterator1) = candidate
                            iterator1 = iterator1 + 1
                        END IF
                    END IF
                    IF(dummy_off2(i) /= 0) THEN !! "bits" zero must not overlap
                        candidate = dummy_off2(i)
                        IF(ANY(candidate == off2) .EQV. .FALSE.) THEN !! Check for conflicts 
                            off2(iterator2) = candidate
                            iterator2 = iterator2 + 1
                        END IF
                    END IF    
                END DO
                DO i = 1, N !! If a cycle is such that not all no-conflict "bits" in It
                            !! are filled, then do It
                    IF(ANY(i == off1)) THEN
                        IF(ANY(i == off2) .EQV. .FALSE.) THEN
                            off2(iterator2) = i
                            iterator2 = iterator2 + 1
                        END IF    
                    END IF
                    IF(ANY(i == off2)) THEN
                        IF(ANY(i == off1) .EQV. .FALSE.) THEN
                            off1(iterator1) = i
                            iterator1 = iterator1 + 1
                        END IF    
                    END IF
                END DO    
                IF(ANY(0 == off1)) THEN !! Here the process re-start ONLY IF offspring 1 is not full
                    dummy_par1 = 0
                    dummy_par2 = 0
                    j = 1 !! dummy_par1 iterator
                    k = 1 !! dummy_par2 iterator
                    DO i = 1, N !! Prepare parents with no conflicts and re-start the procedure
                        IF(ANY(par1(i) == off1) .EQV. .FALSE.) THEN
                            IF(ANY(par1(i) == off2) .EQV. .FALSE.) THEN
                                dummy_par1(j) = par1(i)
                                j = j + 1
                            END IF
                        END IF
                        IF(ANY(par2(i) == off1).EQV. .FALSE.) THEN
                            IF(ANY(par2(i) == off2).EQV. .FALSE.) THEN
                                dummy_par2(k) = par2(i)
                                k = k + 1
                            END IF
                        END IF                        
                    END DO
                
                ELSE 
                    
                    EXIT
                    
                END IF
            END DO    
                
        ELSE
            PRINT*,"Operator not recognized: to choose among:"
            PRINT*,"1 for 'PMX' - 2 for 'OX' - 3 for 'CX2'"
        END IF
        
        !PRINT*,"CROSS CHECK" !! Debug
        !CALL CHECK(off1)
        !CALL CHECK(off2)
        
    END SUBROUTINE


    
    SUBROUTINE CHECK(seq)
    !!-------------------------------------------------------------------
    !! Checks If a sequence is healthy --- Debug
    !!-------------------------------------------------------------------
        INTEGER, DIMENSION(N), INTENT(IN) :: seq
        
        INTEGER, DIMENSION(SIZE(seq)) :: dummy_seq
        INTEGER :: i
        
        
        
        DO i = 1, N
            dummy_seq = seq
            dummy_seq(i) = 0
            IF(ANY(seq(i) == dummy_seq))PRINT*,"CONFLICT FOUND IN SEQUENCE"
            IF(seq(i) == 0)PRINT*,"ZERO FOUND IN SEQUENCE"
        END DO    

    END SUBROUTINE



    SUBROUTINE OUTPUT_GEN(generation, population, fitness)
    !!-------------------------------------------------------
    !! To save output for each generation in the run
    !!-------------------------------------------------------
        INTEGER, INTENT(IN) :: generation
        REAL(KIND=nkind), DIMENSION(N_pop), INTENT(IN) :: fitness
        INTEGER, DIMENSION(N_pop,N), INTENT(IN) :: population
    
        REAL(KIND=nkind) :: best, fit, fit2, dev
        INTEGER :: i, best_pos
    
    
        
        OPEN(unit=1,file="data_gen.dat", position="append")
        !! Takes observables
        fit  = 0
        fit2 = 0
        DO i = 1, N_pop
            fit  = fit  + fitness(i)
            fit2 = fit2 + fitness(i)*fitness(i)
        END DO
        dev = sqrt(fit2/N_pop - (fit/N_pop)*(fit/N_pop))
        best = MAXVAL(fitness)
        best_pos = MAXLOC(fitness, 1)
        WRITE(1,*)generation, fit/N_pop, dev, best, population(best_pos, :)
        CLOSE(1)

    END SUBROUTINE



    SUBROUTINE OUTPUT_SIM(simulation, population, fitness)
    !!-------------------------------------------------------
    !! To save output for each generation in the run
    !!-------------------------------------------------------
        INTEGER, INTENT(IN) :: simulation
        REAL(KIND=nkind), DIMENSION(N_gen), INTENT(IN) :: fitness
        INTEGER, DIMENSION(N_gen,N), INTENT(IN) :: population
    
        REAL(KIND=nkind) :: best, fit, fit2, dev
        INTEGER :: i, best_pos
    
    
        
        OPEN(unit=1,file="data_sim.dat", position="append")
        !! Takes observables
        fit  = 0
        fit2 = 0
        DO i = 1, N_gen
            fit  = fit  + fitness(i)
            fit2 = fit2 + fitness(i)*fitness(i)
        END DO
        dev = sqrt(fit2/N_gen - (fit/N_gen)*(fit/N_gen))
        best = MAXVAL(fitness)
        best_pos = MAXLOC(fitness, 1)
        WRITE(1,*)simulation, fit/N_gen, dev, best, population(best_pos, :)
        CLOSE(1)

    END SUBROUTINE



    RECURSIVE SUBROUTINE permutate(E, P)
    !!---------------------------------------------------------------
    !! Takes a sequence and returns a matrix with all the possible
    !! permutation of the "bits" in the sequence
    !!---------------------------------------------------------------
        INTEGER, INTENT(IN)  :: E(:)       ! array of objects
        INTEGER, INTENT(OUT) :: P(:,:)     ! permutations of E
        INTEGER  :: N, Nfac, i, k, S(SIZE(P,1)/SIZE(E), SIZE(E)-1)
    
    
    
        N = SIZE(E)
        Nfac = SIZE(P,1)
    
        DO i = 1, N                           ! cases with E(i) in front
              IF( N>1 ) CALL permutate((/E(:i-1), E(i+1:)/), S)
              FORALL(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/)
        END DO
        
    END SUBROUTINE permutate



END MODULE common




!!------------------------------------------------------------------------------------------------

!! HERE STARTS THE RUN: 

!! 1) BUILD THE PLACE RANDOMLY OR IMPORT A CONFIGURATION 

!! 2) N_sim SIMULATIONS OF GA WITH CHOSEN CROSSOVER OPERATORS AND FIXED PARAMETERS

!! 3) BRUTE-FORCE SIMULATION TO COMPARE THE RESULTS

!! LOOKING FOR AVERAGE FITNESS AND TIMES GA GOT A GOOD RESULT FOR EACH CROSSOVER OPERATOR

!!------------------------------------------------------------------------------------------------


PROGRAM GS_for_TSP
    USE common
    INTEGER :: i, x, y, gen, op, N_elite, elite_iterator
    INTEGER, DIMENSION(2) :: ind !! Indices of the parents chosen for breeding
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: best_list !! array of best sequence for each generation
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: sim_best_list !! array of best sequence for each simulation
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: P !! Matrix of all possible sequences
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: best_seq !! best seq of a generation
    
    REAL(KIND=nkind), DIMENSION(:), ALLOCATABLE :: fit_best !! Array of best fitness for each generation
    REAL(KIND=nkind), DIMENSION(:), ALLOCATABLE :: fit_sim !! Array of best fitness for each simulation       
    INTEGER, DIMENSION(:), ALLOCATABLE :: worst_positions
    INTEGER, DIMENSION(:), ALLOCATABLE :: best_positions, best_fitnesses
   
    INTEGER :: sim, best_pos, worst_pos,sim_elite, elite_counter !! sim_elite=first simulation in which emerged elite
    REAL(KIND=nkind) :: best, worst, dummy_best, dummy_worst !! best and worst fitness of the simulation 
    REAL(KIND=nkind) :: t_start, t_end, elitism
    LOGICAL :: ELITE_FLAG = .TRUE.


    elitism = 0.15 
    N_elite = int(elitism*N_pop)
    ALLOCATE(best_positions(N_elite), best_fitnesses(N_elite), best_seq(N_elite,N), worst_positions(N_elite))
    IF(BF_FLAG)ALLOCATE(P(PRODUCT((/(i, i=1,N)/)),  N))
    ALLOCATE(best_list(N_gen,N),sim_best_list(N_sim,N))
    ALLOCATE(fit_best(N_gen), best_GA_seq(N), fit_sim(N_sim))
    
    
    
    CALL INITIAL() !! Set initial configuration
    
    IF(RND_PLACE_FLAG) THEN
        distances = -1
        !!--------------------------------------------------------------------------------------------
        !! Build distances matrix to avoid a lot of calculations
        !! Optimization to write in the TeX
        !!--------------------------------------------------------------------------------------------
        DO i = 1, N
            DO j = i, N
                CALL distance(pos(i,1), pos(i,2), pos(j,1), pos(j,2), distances(i,j))
            END DO
        END DO

        DO i = 1, N
            DO j = 1, N
                IF(distances(i,j) < 0)distances(i,j)=distances(j,i)
            END DO
        END DO
    ELSE
        PRINT*,"INSERT DISTANCES MATRIX"
        READ*,distances
    END IF
    ! PRINT*,"DISTANCES" !! DEBUG
    ! DO i = 1, N
    !     PRINT*,distances(:,i)
    ! END DO
    !!--------------------------------------------------------------------------------------------





    !! Choose crossover operator - 1 FOR PMX, 2 FOR OX, 3 FOR CX2
    op = 3
    
    
    
    !! Set values
    best_list = 0
    fit_best = 0
    best = 0
    worst = L*L
    
    CALL CPU_TIME(t_start)
    DO sim = 1, N_sim
    
        OPEN(unit=1,file="data_gen.dat", position="append")
        WRITE(1,*)"#INDEX",sim-1
        CLOSE(1)
        
        !! Generate initial parents population and compute their fitness
        DO i = 1, N_pop 
            CALL seq_gen(N, parents(i,:))
            CALL fit(parents(i,:), distances, fit_p(i))
        END DO   
        CALL OUTPUT_GEN(0,parents,fit_p)
        
        DO gen = 1, N_gen
        !! Starts evolution - a step is when you CAN replace all the parents 
            DO i = 1, N_pop, 2
                offspring(i,:) = parents(i,:)     !! Set the new generation to be like the old one
                offspring(i+1,:) = parents(i+1,:) !! Not all of the times a parents choice
                                                  !! lead to an offspring
                                              
                CALL roulette(fit_p, ind) !! Choose two parents to breed
                !! Elite possible improvement, or different weight to major fitness pressure
                !! To write in the TeX file
                
                
                !! Crossover: choose between 1 for "PMX", 2 for "OX", 3 for "CX2"
                CALL random_number(rnd)                
                IF(rnd < p_cross)CALL CROSS(parents(ind(1),:), parents(ind(2),:), offspring(i,:), offspring(i+1,:), op=op)

                
                !! Mutation: basic two-"bits" swap
                CALL random_number(rnd)
                IF(rnd < p_mut)CALL MUT(offspring(i,:), offspring(i,:))
                CALL random_number(rnd)
                IF(rnd < p_mut)CALL MUT(offspring(i+1,:), offspring(i+1,:))
                
                !! Compute offspring fitness
                CALL fit(offspring(i,:), distances, fit_o(i))
                CALL fit(offspring(i+1,:), distances, fit_o(i+1))
            END DO

            IF(ELITE_FLAG .EQV. .FALSE.) THEN
                DO j = 1, N_pop
                    parents(j,:) = offspring(j,:) !! Replace the individuals for a new generation step
                    fit_p(j) = fit_o(j)
                END DO
            END IF
            
            IF(ELITE_FLAG) THEN
                !! The worst of a generation is replaced with
                !! the best of the old one
                dummy_fit_p = fit_p
                dummy_fit_o = fit_o
                DO i = 1, N_elite
                    best_positions(i) = MAXLOC(dummy_fit_p, 1)
                    best_fitnesses(i) = fit_p(best_positions(i))
                    best_seq(i,:) = parents(best_positions(i),:)

                    worst_positions(i) = MINLOC(dummy_fit_o, 1)
       
                    dummy_fit_p(best_positions(i)) = 0
                    dummy_fit_o(worst_positions(i)) = L*L
                END DO
                
                elite_iterator = 1
                DO j = 1, N_pop
                    IF(ANY(j == worst_positions)) THEN
                        parents(j,:) = best_seq(elite_iterator, :) !! Replace the individuals for a new generation step
                        fit_p(j) = best_fitnesses(elite_iterator)
                        elite_iterator = elite_iterator + 1
                    ELSE    
                        parents(j,:) = offspring(j,:) !! Replace the individuals for a new generation step
                        fit_p(j) = fit_o(j)
                    END IF    
                END DO
            END IF
            
            !! Take the best of each generation
            fit_best(gen) = MAXVAL(fit_p, 1)
            best_pos = MAXLOC(fit_p, 1)
            best_list(gen,:) = parents(best_pos,:)
            !! Taken
            
            CALL OUTPUT_GEN(gen,parents,fit_p)
        END DO
        !!-----------------------------------------------------------------------
        !! Evolution ended
        !!-----------------------------------------------------------------------

        !! Takes the best of all simulations choosing
        !! from the best of each generation
        !! for all simulations
        DO i = 1, N_gen
            IF(fit_best(i) > best) THEN
                best = fit_best(i)
                best_GA_seq = best_list(i,:)
            END IF
        END DO
        best_pos = MAXLOC(fit_best, 1)
        fit_sim(sim)=fit_best(best_pos)
        sim_best_list(sim,:) = best_list(best_pos,:)
        !! Taken If worth
        
        !! Takes the worst from the last generation
        !! of all the simulations
        dummy_worst = MINVAL(fit_p, 1)
        IF(worst > dummy_worst)worst=dummy_worst
        !! Taken
        
        !! Each index in this file is a simulation
        OPEN(unit=1,file="data_gen.dat", position="append")
        WRITE(1,*)"" !! For GNUplot index
        WRITE(1,*)""
        CLOSE(1)
        
        !! Write simulation output to see the best of each
        !! generation in comparison
        CALL OUTPUT_SIM(sim,best_list,fit_best)
    END DO
    CALL CPU_TIME(t_end)
    !!-----------------------------------------------------------------------
    !! Simulation ended
    !!-----------------------------------------------------------------------
    PRINT*,"Data of population per generation saved in 'data_gen.dat'"
    PRINT*,"Data of the best individuals per generation saved in 'data_sim.dat'"
    
    elite_counter = 0 !! Times elite emerged from evolution
    DO sim=1, N_sim
        IF(fit_sim(sim) == best)elite_counter=elite_counter+1
    END DO
    
    sim_elite = 0 !! First simulation in which elite emerged
    DO sim=1, N_sim
        IF(fit_sim(sim) == best) THEN
            sim_elite=sim
            EXIT
        END IF    
    END DO
    
    
    
    
    
    
    
    OPEN(unit=17,file="analysis.dat", position="append")
    !!------------------------------------------------
    !! GA results are printed for a faster comparison
    !! with Brute-force method
    !!------------------------------------------------
    PRINT*,"--------------------------------------"
    PRINT*,"N =",N,"PROBLEM WITH N_sim=",N_sim
    PRINT*,"GA BEST SEQUENCE",best_GA_seq
    PRINT*,"WITH FITNESS",best
    PRINT*,"WORST FITNESS",worst
    PRINT*,"IN",t_end-t_start,"s"
    PRINT*,"FIRST APPEARED IN SIMULATION N", sim_elite
    PRINT*,"APPEARED", elite_counter, "TIMES"
    PRINT*,""
    WRITE(17,*)N,t_end-t_start,elite_counter
    !!------------------------------------------------   
    !! Brute-force computation of the best sequence
    !! and fitness
    !!------------------------------------------------
    IF(BF_FLAG) THEN
        CALL CPU_TIME(t_start)
        CALL permutate( (/(i, i=1,n)/),  P ) !! Compute all permutations

        best = 0
        best_pos = 1
        DO i = 1, SIZE(P,1)
            CALL fit(P(i,:), distances, dummy_best)
            IF(dummy_best>best) THEN
                best = dummy_best
                best_pos = i
            END IF
        END DO
        CALL CPU_TIME(t_end)
        PRINT*,"BRUTE FORCE BEST SEQUENCE", P(best_pos,:)
        PRINT*,"WITH FITNESS", best
        PRINT*,"IN",t_end-t_start,"s"
        WRITE(17,*)N,t_end-t_start, best
    END IF


    DEALLOCATE(seed, box)
    IF(RND_PLACE_FLAG)DEALLOCATE(place, pos)
    DEALLOCATE(distances, best_list)
    DEALLOCATE(parents, offspring, fit_p, fit_o, best_GA_seq)
    IF(BF_FLAG)DEALLOCATE(P)
    
END PROGRAM    
