
module resources_montecarlo
    use std_types
    use helpers
    use nakajima_zwanzig_shared
    use module_goft

    implicit none

    ! declarations
    integer(i4b)    ::    TRAJECTORIES

    real(dp), dimension(:), allocatable    :: global_lambda, global_temp, global_tau
    real(dp), dimension(:,:), allocatable  :: Ham
    real(dp)                               :: rwa
    logical ::  origin = .false., g_functions = .true.,                                        &
                                       modified_unraveling2 = .false.,                         &
                                       vynechani_G_ifu = .false., use_exciton_basis = .false.

    character(len=64), parameter :: external_dir = "external", config_filename = "config.prm", directory = "."

    integer(i4b)    ::  STEPS = -1
    integer(i4b)    ::  Nl                                   ! number of one-excitons
    real(dp)        ::  timeStep                             ! timestep of the method
    real(dp)        ::  jumps_in_one_run                     ! jumps in average in one run
    real(dp)        ::  p_of_no_jump                         ! probability of no jump (important for restarting)
    !complex(dpc), dimension(:,:,:), allocatable       :: rho, rho_coherent

    complex(dpc), dimension(:,:,:,:,:), allocatable   :: Ueg, Uee, Ugg

    logical, parameter :: only_coherences = .true.

    interface generate_trajectory
        module procedure generate_trajectory_re
        module procedure generate_trajectory_cmplx
    end interface

    interface jump_probability
        module procedure jump_probability_144
        module procedure jump_probability_444
    end interface
    interface jump_probability_total
        module procedure jump_probability_total_
        module procedure jump_probability_total_IJ
    end interface

    contains

    !
    ! Do all the simulations within the module
    !
    subroutine do_montecarlo_work()
        integer(i4b)            :: r,s,k,l,i
        complex(dpc), dimension(:,:,:), allocatable :: rho, rho_coherent
        character(len=128)      :: cbuff
        character :: type

        if(only_coherences) then
            type = 'O'
        else
            type = 'E'
        end if

        call init_monte_carlo()

        allocate(rho(Nl,Nl,STEPS))
        allocate(rho_coherent(Nl,Nl,STEPS))

        if(vynechani_G_ifu) then
            write(cbuff, '(A F6.4 A)') "WARNING! FAST_G METHOD CONTAINS SERIOUS ERROR, USE ONLY FOR TESTING!"
            call print_warning_message(cbuff, -1)
        end if

        write(cbuff, '(A L)') "G-functions ", g_functions
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "Exciton basis ", use_exciton_basis
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "Modified Unraveling 2 ", modified_unraveling2
        call print_log_message(cbuff, 5)
        write(cbuff, '(A L)') "FastG ", vynechani_G_ifu
        call print_log_message(cbuff, 5)

        !*************************************************************
        ! Calculation of evolution superops
        !*************************************************************
        do r=N1_from_type(type),1,-1
        do s=1,N2_from_type(type)
            rho = 0.0_dp
            rho_coherent = 0.0_dp

            call perform_montecarlo_parallel(r,s,rho,rho_coherent)

            do k=1,N1_from_type(type)
            do l=1,N2_from_type(type)

            rho(k,l,1) = 0.0_dp
            rho(r,s,1) = 1.0_dp

            rho_coherent(k,l,1) = 0.0_dp
            rho_coherent(r,s,1) = 1.0_dp

            ! saving result into evolution superoperator
            do i=1,STEPS

                if(only_coherences) then
                    Ueg(k,l, r,s, i) = rho(k,l,i)
                else
                    Uee(k,l, r,s, i) = rho(k,l,i)
                endif

            end do

            if(only_coherences) then
                write(cbuff,'(a,i1,i1,i1,i1,a)') 'Ueg(',k,l,r,s,') filled'
            else
                write(cbuff,'(a,i1,i1,i1,i1,a)') 'Uee(',k,l,r,s,') filled'
            end if

            call print_log_message(trim(cbuff),5)

            end do
            end do

        end do
        end do

        ! transformation into excitonic picture
        if(use_exciton_basis) then
            do i=1,STEPS
                if(only_coherences) then
                    call superops_to_exc(Ueg(:,:, :,:, i), 'O')
                else
                    call superops_to_exc(Uee(:,:, :,:, i), 'E')
                endif

            end do
        end if

        if(only_coherences) then
            call open_files('O')
            call write_evolution_operators('O',Ueg)
            call close_files('O')
        end if

    end subroutine do_montecarlo_work



    subroutine init_monte_carlo()
        use_exciton_basis = .false.

        STEPS = 2048
        rwa = 0
        Nl = 2

        allocate(global_lambda(Nl))
        allocate(global_temp(Nl))
        allocate(global_tau(Nl))
        allocate(Ham(Nl,Nl))
        !allocate(Ueg(N1_from_type('O'), N2_from_type('O'), N1_from_type('O'), N2_from_type('O'), Nt(1)  ))
        !allocate(Uee(N1_from_type('E'), N2_from_type('E'), N1_from_type('E'), N2_from_type('E'), Nt(1)  ))
        !allocate(Ugg(1,1,1,1, Nt(1)  ))

        global_temp = 0.0
        global_tau = 0.0
        global_lambda = 0.0
        Ham = 0.0

        !global_en(1) = -400/Energy_internal_to_cm
        !global_en(2) = 0/Energy_internal_to_cm
        !J_coupl(1,2) = 300/Energy_internal_to_cm
        !J_coupl(2,1) = 300/Energy_internal_to_cm
        !Ham(1,1) = global_en(1)
        !Ham(2,2) = global_en(2)
        !Ham(1,2) = J_coupl(1,2)
        !Ham(2,1) = J_coupl(1,2)

        !global_temp(1) = 300
        !global_temp(2) = 300
        !global_lambda(1) = 100/Energy_internal_to_cm
        !global_lambda(2) = 1000/Energy_internal_to_cm
        !global_tau(1) = 100
        !global_tau(2) = 100

        call read_config_file(Nl)

        write(*,*) ';lambda',global_lambda
        write(*,*) ';tau',global_tau
        write(*,*) ';temp',global_temp
        write(*,*) ';Ham', Ham

        call init_random_seed()
        call init_nakajima_zwanzig_shared(Ham)

        TRAJECTORIES = 10000*250
        !jumps_in_one_run = 18

        if(jump_probability_total() < 1e-5) then
            timeStep = 0.25 !(dt*gt(1))/RUNS
        else
            write(*,*) jumps_in_one_run,',',jump_probability_total(),',',STEPS
            timeStep = jumps_in_one_run/jump_probability_total()/STEPS
            write(*,*) timeStep
        end if

        jumps_in_one_run = jump_probability_total()*timeStep*STEPS

        p_of_no_jump = (1.0_dp - jump_probability_total()*timeStep)**STEPS

        call init_goft(Nl, STEPS*2, global_lambda, global_tau, global_temp, timeStep)

        allocate(Ueg(Nl, 1,Nl, 1,STEPS) )
        allocate(Uee(Nl,Nl,Nl,Nl,STEPS) )
        allocate(Ugg(1,1,1,1,STEPS) )

        Ueg = 0.0_dp
        Uee = 0.0_dp
        Ugg = 1.0_dp
    end subroutine init_monte_carlo

    subroutine clean_montecarlo()
    end subroutine clean_montecarlo

    subroutine perform_montecarlo_parallel(i0, j0, rho, rho_coherent)
        integer(i4b), intent(in)                     :: i0, j0
        complex(dpc), dimension(:,:,:), intent(out)  :: rho, rho_coherent
        complex(dpc), dimension(Nl,Nl,STEPS)  :: rho_local, rho_coherent_local

        integer(i4b), dimension(3) :: time
        character(len=256) :: buff
        integer(i4b) :: i

        !call omp_set_num_threads(2)

        write(buff,'(f12.3)') jump_probability_total()*timeStep
        buff = 'Averagely ' // trim(buff) // ' jumps in one run'
        call print_log_message(trim(buff),5)

        rho = 0.0_dp
        rho_coherent = 0.0_dp

        ! $ OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) PRIVATE(i, rho_local, rho_coherent_local) REDUCTION(+:rho,rho_coherent)
        do i=1,360
          call perform_montecarlo(i0, j0, rho_local, rho_coherent_local, TRAJECTORIES/360)

          rho = rho + rho_local
          rho_coherent = rho_coherent + rho_coherent_local

          call itime(time)
          write(buff, '(A I2 A I2.2 A I2.2 A I6 )') '  time ',time(1), ':', time(2), ':', time(3) ,  &
                 ', STEPS/360: ', i
          call print_log_message(trim(buff), 5)
        end do

    end subroutine perform_montecarlo_parallel

    subroutine perform_montecarlo(i0, j0, rho, rho_coherent, traj)
        integer(i4b), intent(in)                     :: i0, j0
        complex(dpc), dimension(:,:,:), intent(out)  :: rho, rho_coherent
        integer(i4b), intent(in) :: traj

        integer(i4b) :: i, j, a, b, jj
        real(dp)     :: r, norm

        integer(i1b), dimension(2,STEPS) :: draha
        complex(dpc), dimension(STEPS)   :: factor, Gfactor, Ifactor
        complex(dpc) :: factor_in
        complex(dpc), dimension(Nl, Nl) :: rho_init

        factor_in = 1.0_dp
        Gfactor = cmplx(1,0,dpc)

        do i=1, traj
            call random_number(r)
            a = INT(r * Nl + 1.0_dp)
            call random_number(r)
            b = INT(r * Nl + 1.0_dp)

            if(only_coherences) then
                b = 1
            end if

            rho_init(:,:) = 0.0_dp
            rho_init(i0,j0) = 1.0_dp

            call generate_trajectory(draha, factor, 1.0_dp, i0, j0, i0, j0)

            if(g_functions) then
                call calculate_Gfactor_from_trajectory_history(draha,cmplx(1,0,dp),Gfactor)
            end if

            call calculate_Ifactor_from_trajectory_history(draha,cmplx(1,0,dp),Ifactor)
            !Ifactor = Ifactor/Ifactor(1)    ! coherent at the restart point
            !Ifactor(1:STEPS*(run-1)) = 0.0_dp             ! zero before it

            do j=1, STEPS
                if(maxval(draha(:,j)) > 0 .and. minval(draha(:,j)) <= Nl) then
                    rho(draha(1,j),draha(2,j),j) =             &
                    rho(draha(1,j),draha(2,j),j)               &
                    + factor(j)*conjg(Gfactor(j))*Ifactor(j)

                    rho_coherent(draha(1,j),draha(2,j),j) =    &
                    rho_coherent(draha(1,j),draha(2,j),j)      &
                    + factor(j)*Ifactor(j)
                end if
            end do

        end do

        !norm = abs(maxval(abs(rho_coherent(:,1,1:STEPS*run))))
!        norm = abs(maxval(abs(rho_coherent(:,1,1:1 + 5)))) ! do not include the end
!        do j=STEPS, 1, -1
!            if(.not. only_coherences) then
!                rho(:,:,j) = rho(:,:,j)/trace(rho(:,:,j))
!                rho_coherent(:,:,j) = rho_coherent(:,:,j)/trace(rho_coherent(:,:,j))
!            else
!                do jj=1, Nl
!                  rho(jj,1,j) = rho(jj,1,j)/norm
!                  rho_coherent(jj,1,j) = rho_coherent(jj,1,j)/norm
!                end do
!            end if
!        end do

    end subroutine perform_montecarlo

    recursive function goft_site(m,n,t) result(dgoft)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: dgoft
        integer    (i4b)            :: t_index

        if(t < 0) then
            dgoft = conjg(goft_site(n,m,-t))
            return
        end if

        !t_index = INT(t/dt)+1
        t_index = INT(t/timeStep)+1

        dgoft = 0.0_dp

        if(m == n) then

            dgoft = goft(m,t_index) !@!

        end if

    end function goft_site

    recursive function hoft_site(m,n,t1,t2) result(dgoft)
    ! hoft(a,a,t1,t2) = hoft(a,a,-t2,-t1)
        real(dp), intent(in)        :: t1,t2
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: dgoft

        dgoft = 0.0_dp

        if(m == n) then
            dgoft = goft_site(m,m,t1) - goft_site(m,m,t1-t2) + goft_site(m,m,-t2)
        end if

    end function hoft_site

    !
    ! Defines the jump probability at given circumstances
    !
    function jump_probability_444(state_where_I_am,state_where_to_go,side) result(res)
        integer(i4b), intent(in)   :: state_where_I_am
        integer(i4b), intent(in)   :: side, state_where_to_go
        real(dp)                   :: res

        res = Ham(state_where_I_am,state_where_to_go)

        if(state_where_I_am == state_where_to_go) then
          res = 0.0
        end if

    end function jump_probability_444

    ! just for interface-purposes
    function jump_probability_144(state_where_I_am,state_where_to_go,side) result(res)
        integer(i1b), intent(in)                 :: state_where_I_am
        integer(i4b), intent(in)                 :: side, state_where_to_go
        real(dp)                                 :: res

        integer(i4b)                             :: a

        a = state_where_I_am

        res = jump_probability_444(a,state_where_to_go,side)
    end function jump_probability_144

    !
    ! Total jump probability at given circumstances
    !
    function jump_probability_total_IJ(i0, j0) result(cumulative_probability)
        integer(i4b), intent(in)   :: i0, j0
        real(dp)                   :: cumulative_probability

        integer(i4b)               :: a

        cumulative_probability = 0.0_dp
        do a=1,Nl
            cumulative_probability = cumulative_probability + jump_probability(i0, a, 1)
        end do
        if(.not. only_coherences) then
            do a=1,Nl
                cumulative_probability = cumulative_probability + jump_probability(j0, a, 2)
            end do
        end if
    end function

    function jump_probability_total_() result(cumulative_probability)
        real(dp)                   :: cumulative_probability

        integer(i4b)               :: a

        cumulative_probability = 0.0_dp
        do a=1,Nl
            cumulative_probability = cumulative_probability + jump_probability_total_IJ(a,a)
        end do

    end function

    !
    ! Defines the jump phase change at given circumstances
    !
    function jump_phase_change(trajectory,side,i,a) result(res)
        integer(i1b), intent(in), dimension(:,:) :: trajectory
        integer(i4b), intent(in)                 :: side, i,a
        complex(dpc)                             :: res

        res = 0.0_dp

        if(side == 1) then
            res = cmplx(0,-1,dpc)
        else if(side == 2) then
            res = cmplx(0,+1,dpc)
        end if
    end function jump_phase_change    !

    !
    ! Generates trajectory with STEPS steps and related pure-coupling-factors
    !
    subroutine generate_trajectory_re(trajectory, factor_out, factor_in, i0, j0, i00, j00)
        integer(i1b), intent(out), dimension(:,:) :: trajectory
        integer(i4b), intent(in)    :: i0, j0, i00, j00
        real(dp), intent(in)        :: factor_in
        complex(dpc), intent(out), dimension(:)    :: factor_out

        complex(dpc) :: ffactor

        ffactor = factor_in

        call generate_trajectory_cmplx(trajectory, factor_out, ffactor, i0, j0, i00, j00)
    end subroutine generate_trajectory_re

    recursive subroutine generate_trajectory_cmplx(trajectory, factor_out, factor_in, i0, j0, i00, j00)
        integer(i1b), intent(out), dimension(:,:) :: trajectory
        integer(i4b), intent(in)    :: i0, j0, i00, j00 ! beginning of run-th run and beginning of 1-st run
        complex(dpc), intent(in)    :: factor_in
        complex(dpc), intent(out), dimension(:)    :: factor_out

        integer(i4b) :: i, a, b, side, inside_run, u, v
        real(dp) :: cumulative_random, cumulative_probability

        if(size(trajectory,2) /= STEPS .or. size(trajectory,1) /= 2 .or. &
            size(factor_out,1) /= STEPS) then
            call print_error_message(-1, "dimension error in generate_trajectory()")
        end if

        if((.not.(i0 == -1 .and. j0 == -1))) then
            ! if depository is on or we are in first recursion, we zero the trajectory and factor
            factor_out = 0.0_dp
            trajectory = 0

            trajectory(1, 1) = i0
            trajectory(2, 1) = j0
        end if

        factor_out(1) = factor_in
        trajectory(1,1) = i0
        trajectory(2,1) = j0

        do i=1, STEPS-1
            factor_out(i+1) = factor_out(i)
            trajectory(1,i+1) = trajectory(1,i)
            trajectory(2,i+1) = trajectory(2,i)

            call random_number(cumulative_random)
            cumulative_probability = 0.0_dp

            do a=1,Nl
            do side=1,2
                if(only_coherences .and. side == 2) then
                    exit
                end if

                cumulative_probability = cumulative_probability + jump_probability(trajectory(side,i),a,side)*timeStep

                if(cumulative_random < cumulative_probability) then
                    trajectory(side,i+1) = a
                    factor_out(i+1) = factor_out(i)*jump_phase_change(trajectory,side,i,a)
                    if(.not. modified_unraveling2) then
                        factor_out(1:i) = 0.0_dp
                    end if
                    goto 42
                end if

            end do
            end do

42          continue
        end do

        ! normalization according to the jump probability
        if(modified_unraveling2) then
            cumulative_probability = jump_probability_total(i0,j0)*timeStep

            do i=1, STEPS-1
                factor_out(i) = factor_out(i) / (1.0 - cumulative_probability)**(i)
            end do
        end if

    end subroutine generate_trajectory_cmplx

    subroutine calculate_Ifactor_from_trajectory_history(trajectory,initial_factor,factor_out)
        integer(i1b), dimension(:,:), intent(in)     :: trajectory
        complex(dpc), intent(in)                     :: initial_factor
        complex(dpc), dimension(:), intent(out)      :: factor_out

        integer(i4b)  :: i

        if(size(trajectory,1) /= 2 .or. size(factor_out) /= size(trajectory,2)) then
            call print_error_message(-1, "dimension error in calculate_Ifactor_from_trajectory_history()")
        end if

        if(maxval(trajectory) > Nl) then
            call print_error_message(-1, "value error in calculate_Ifactor_from_trajectory_history()")
        end if

        factor_out = 0.0_dp
        factor_out(1) = initial_factor

        do i=2,size(factor_out)
            if(minval(trajectory(:,i)) < 1) then
                exit
            else
                if(.not. only_coherences) then
                    factor_out(i) = factor_out(i-1)*exp(cmplx(0,-1,dpc)*timeStep*&
                      (Ham(trajectory(1,i),trajectory(1,i)) - Ham(trajectory(2,i),trajectory(2,i))) )
                else
                    factor_out(i) = factor_out(i-1)*exp(cmplx(0,-1,dpc)*timeStep*&
                      (Ham(trajectory(1,i),trajectory(1,i)) - rwa) )
                end if
            end if
        end do

    end subroutine calculate_Ifactor_from_trajectory_history

     ! Let we have U^+_A(t1)U^+_B(t2)U^+_C(t3)...U_c(t3)U_b(t2)U_a(t1)W, then the cummulant expansion is derived
     ! from it starts as (1-i h_A(t1)-g_A(t1))^*(1-i h_B(t1)-g_B(t1))(1-i h_B(t1+t2)-g_B(t1+t2))^*
     ! (1-i h_C(t1+t2)-g_C(t1+t2)) ... (1-i h_Z(t1+..+tn)-g_Z(t1+..+tn))^* (1-i h_z(t1+..+tn)-g_z(t1+..+tn))
     ! (1-i h_z(t1+..+t(n-1))-g_z(t1+..+t(n-1)))^* (1-i h_y(t1+..+t(n-1))-g_y(t1+..+t(n-1)))
    subroutine calculate_Gfactor_from_trajectory_history(trajectory,initial_factor,factor_out)
        integer(i1b), dimension(:,:), intent(in)  :: trajectory
        complex(dpc), intent(in)                  :: initial_factor
        complex(dpc), dimension(:), intent(out)   :: factor_out

        ! we overallocate jump-fields not to have to always calculate their dimensions
        integer(i4b), parameter :: overallocation = 10000
        integer(i4b), dimension(0:overallocation)    ::    &
                            jumps_index, jumps_tovalue_ket, jumps_tovalue_bra
        real(dp), dimension(0:overallocation)        ::    &
                            jumps_time
        logical, dimension(0:overallocation)         ::    &
                            jumps_ket

        integer(i4b) ::     i, j, k, number_of_jumps, tovalue,        &
                            fromvalue, ket_i, bra_i, ket_j, bra_j

        real(dp)     :: time_i, time_j, time
        logical      :: last_round_i
        complex(dpc) :: additive_factor

        complex(dpc) :: additive_factor_not_last_turn_z_minuleho
        integer(i4b) :: max_i_minuleho
        logical      :: lze_pouzit_AFNLTZM

        logical, parameter :: debug_G = .false.

        if(size(trajectory,1) /= 2 .or. size(factor_out) /= size(trajectory,2)) then
            call print_error_message(-1, "dimension error in calculate_factor_from_trajectory_history()")
        end if

        ! we prepare fields representing the jumps
        jumps_index(0)          = 1
        jumps_time(0)           = 0.0_dp
        jumps_tovalue_ket(0)    = trajectory(1,1)
        jumps_tovalue_bra(0)    = trajectory(2,1)

        tovalue = 1
        fromvalue = 1
        i = 1
        do while(.true.)
            if(tovalue == -1 .or. fromvalue == -1) then
                exit
            end if
            if(i > overallocation) then
                call print_error_message(-1,'overallocation insufficient in calculate_Gfactor_from_trajectory_history()')
            end if

            call next_step_of_trajectory(trajectory, jumps_index(i-1), jumps_index(i), fromvalue, tovalue, jumps_ket(i))

            if(jumps_ket(i)) then
                jumps_tovalue_ket(i) = tovalue
                jumps_tovalue_bra(i) = jumps_tovalue_bra(i-1)
            else
                jumps_tovalue_ket(i) = jumps_tovalue_ket(i-1)
                jumps_tovalue_bra(i) = tovalue
            end if

            jumps_time(i) = (jumps_index(i)-1)*timeStep
            i = i + 1
        end do

        number_of_jumps = i-2
        max_i_minuleho = -1


        ! now we perform evaluation of brackets of g and h-functions

        ! g-functions come first
        factor_out = 0.0_dp

        do k=1, size(factor_out)

            time = (k-1)*timeStep

            if(trajectory(1,k) < 1 .or. trajectory(1,k) > Nl .or. trajectory(2,k) < 1 .or. trajectory(2,k) > Nl) then
                return
            end if

            factor_out(k) = initial_factor
            additive_factor = 0.0_dp

            if(debug_G) then
                write(*,*)
                write(*,*) 'jumps_tovalue_ket', jumps_tovalue_ket(:5)
                write(*,*) 'jumps_tovalue_bra', jumps_tovalue_bra(:5)
                write(*,*) 'jumps_time', jumps_time(:5)
            end if

            ! last operator has time t instead of time of last jump
            do i=1, number_of_jumps + 1
                if(jumps_time(i) >= time .or. i > number_of_jumps) then
                    last_round_i = .true.
                    time_i = time
                    ket_i = trajectory(1,k)
                    bra_i = trajectory(2,k)

                    if(max_i_minuleho == i) then
                        lze_pouzit_AFNLTZM = vynechani_G_ifu
                    else
                        lze_pouzit_AFNLTZM = .false.
                        additive_factor_not_last_turn_z_minuleho = -9999
                        ! if not correctly set later, it will screw up everything and make the mistake obvious
                    end if
                    max_i_minuleho = i

                else
                    last_round_i = .false.
                    time_i = jumps_time(i)
                    ket_i = jumps_tovalue_ket(i)
                    bra_i = jumps_tovalue_bra(i)

                    if(lze_pouzit_AFNLTZM) then
                        additive_factor = additive_factor_not_last_turn_z_minuleho
                        cycle ! we wait until we are in the last round
                    end if
                end if

                ! there are always two 1hg-brackets for both bra and ket for each jump-time, with site i and (i-1).
                ! exception is the last round, where the bra and ket meet into two brackets instead of four

                if(.not. last_round_i) then
                    ! two ket-brackets
                    additive_factor = additive_factor+(-goft_site(ket_i,ket_i,time_i))
                        if(debug_G) then
                            write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if
                    additive_factor = additive_factor+&
                            conjg((-goft_site(jumps_tovalue_ket(i-1),jumps_tovalue_ket(i-1),time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i, jumps_tovalue_ket(i-1)
                        end if

                    ! two bra-brackets
                    if(.not. only_coherences) then
                        additive_factor = additive_factor+conjg((-goft_site(bra_i,bra_i,time_i)))
                            if(debug_G) then
                                write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                        additive_factor = additive_factor+&
                                (-goft_site(jumps_tovalue_bra(i-1),jumps_tovalue_bra(i-1),time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i,jumps_tovalue_bra(i-1)
                            end if
                    end if


                    ! hh between bra and ket brackets of the same i (4 terms => 6 pairs)
                    if(ket_i == bra_i .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(ket_i,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H1', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, bra_i
                            end if
                    endif
                    if(ket_i == jumps_tovalue_bra(i-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(ket_i,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H2', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, jumps_tovalue_bra(i-1)
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == bra_i .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(bra_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H3', real(additive_factor), k, i, time, time_i, last_round_i, jumps_tovalue_ket(i-1), bra_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == jumps_tovalue_bra(i-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_bra(i-1),jumps_tovalue_bra(i-1),&
                                                                    time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H4', real(additive_factor), k, i, time, time_i, last_round_i, jumps_tovalue_ket(i-1), jumps_tovalue_bra(i-1)
                            end if
                    endif
                    if(ket_i == jumps_tovalue_ket(i-1)) then
                        additive_factor = additive_factor+(hoft_site(ket_i,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            end if
                    endif
                    if(bra_i == jumps_tovalue_bra(i-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(bra_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H6', real(additive_factor), k, i, time, time_i, last_round_i, bra_i, bra_i
                            end if
                    endif

                ! one additional factor where bra and kets meet
                else
                    additive_factor = additive_factor+conjg((-goft_site(ket_i,ket_i,time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_g--', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if

                    if(.not. only_coherences) then
                        additive_factor = additive_factor+(-goft_site(bra_i,bra_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg--', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                    end if

                    if(bra_i == ket_i .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(bra_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_h (last)', real(additive_factor), k, i, time, time_i, last_round_i
                            end if
                    end if
                end if


                do j=1, i-1
!                    if(last_round_i .and. j == i-1) then
!                        exit
!                    end if

                    time_j = jumps_time(j)
                    ket_j = jumps_tovalue_ket(j)
                    bra_j = jumps_tovalue_bra(j)

                    ! hh between all combinations of i-j pairs
                    if(ket_i == ket_j .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor1', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == ket_j) then
                        additive_factor = additive_factor+(hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor2', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == ket_j .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor3', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == ket_j .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(ket_j,ket_j,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor4', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(ket_i == jumps_tovalue_ket(j-1) .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,+time_i))
                            if(debug_G) then
                                write(*,*) 'factor5', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == jumps_tovalue_ket(j-1)) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor6', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == jumps_tovalue_ket(j-1) .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor7', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == jumps_tovalue_ket(j-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_ket(j-1),jumps_tovalue_ket(j-1),time_j,+time_i))
                            if(debug_G) then
                                write(*,*) 'factor8', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(ket_i == jumps_tovalue_bra(j-1) .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor9', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == jumps_tovalue_bra(j-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor10', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == jumps_tovalue_bra(j-1) .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_j,time_i))            !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor11', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == jumps_tovalue_bra(j-1) .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(jumps_tovalue_bra(j-1),jumps_tovalue_bra(j-1),time_j,+time_i))        !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor12', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(ket_i == bra_j .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(hoft_site(bra_j,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor13', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_ket(i-1) == bra_j .and. .not. only_coherences) then
                        additive_factor = additive_factor+(-hoft_site(bra_j,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor14', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(bra_i == bra_j .and. .not. only_coherences .and. .not. last_round_i) then
                        additive_factor = additive_factor+(-hoft_site(bra_j,bra_j,time_j,time_i))                          !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor15', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(jumps_tovalue_bra(i-1) == bra_j .and. .not. only_coherences) then
                        additive_factor = additive_factor+(hoft_site(bra_j,bra_j,time_j,time_i))                          !! MOZNA CHYBA V CASECH
                            if(debug_G) then
                                write(*,*) 'factor16', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                end do

                ! At the time of jump, code is not correct. It is safer to exclude these cases than fix it.
                if(abs(time-jumps_time(i)) < timeStep/2.0_dp) then
                    additive_factor = -9999999.0

                    if(debug_G) then
                        write(*,*) 'case excluded', real(additive_factor)
                    end if
                end if

                if(real(additive_factor) > 0 .and. debug_G) then
                    write(*,*) 'zbyly-faktor', real(additive_factor)
                end if

                if(last_round_i) then
                    exit
                else
                    if(.not. lze_pouzit_AFNLTZM) then
                        additive_factor_not_last_turn_z_minuleho = additive_factor
                    end if
                end if
            end do

            factor_out(k) = factor_out(k)*exp(additive_factor)

        end do

        !write(*,*) 'factor', factor_out

    end subroutine calculate_Gfactor_from_trajectory_history

    pure subroutine next_step_of_trajectory(trajectory, from_index, tindex, fromvalue, tovalue, ket)
        integer(i1b), dimension(:,:), intent(in) :: trajectory
        integer(i4b), intent(in)                 :: from_index
        integer(i4b), intent(out)                :: tindex, fromvalue, tovalue
        logical, intent(out)                     :: ket

        integer(i4b)    :: i

        if(from_index < 0) then
            return
        end if

        do i=from_index+1, size(trajectory,2)
            if(trajectory(1,i) /= trajectory(1,i-1)) then
                ket = .true.
                fromvalue = trajectory(1,i-1)
                tovalue = trajectory(1,i)
                tindex = i
                return
            elseif(trajectory(2,i) /= trajectory(2,i-1)) then
                ket = .false.
                fromvalue = trajectory(2,i-1)
                tovalue = trajectory(2,i)
                tindex = i
                return
            end if
        end do

        fromvalue = -1
        tovalue   = -1
        tindex = size(trajectory,2)
    end subroutine next_step_of_trajectory

    pure function number_of_trajectory_steps(trajectory) result(nav)
        integer(i4b) :: nav

        integer(i1b), dimension(:,:), intent(in)    :: trajectory
        integer(i4b)                                :: tindex, fromvalue, tovalue, from_index
        logical                                     :: ket

        tovalue = 1
        fromvalue = 1
        from_index = 1
        nav = -1

        do while(tovalue /= -1)
            nav = nav + 1
            call next_step_of_trajectory(trajectory, from_index, tindex, fromvalue, tovalue, ket)
            from_index = tindex
        end do
    end function number_of_trajectory_steps



    subroutine get_coherent_dynamics(rho,t)
        complex(dpc), dimension(:,:), intent(inout)    :: rho
        real(dp), intent(in)                           :: t

        integer(i4b) :: i,j
        complex(dpc), dimension(size(rho,1),1)         :: rho_O

        if(only_coherences) then
            rho_O(:,1) = rho(:,1)
            call operator_to_exc(rho_O,'O')
            rho(:,1) = rho_O(:,1)
        else
            call operator_to_exc(rho,'E')
        end if

        do i=1,size(rho,1)
        do j=1,size(rho,2)

            if(.not. only_coherences) then
                    rho(i,j) = rho(i,j)*exp(cmplx(0,-1,dpc)*t*&
                      (Ham(i,i) - Ham(j,j)) )
                else
                    rho(i,j) = rho(i,j)*exp(cmplx(0,-1,dpc)*t*&
                      (Ham(i,i) - rwa))
            end if

        end do
        end do

        if(only_coherences) then
            rho_O(:,1) = rho(:,1)
            call operator_from_exc(rho_O,'O')
            rho(:,1) = rho_O(:,1)
        else
            call operator_from_exc(rho,'E')
        end if

    end subroutine get_coherent_dynamics

    pure function ind(i,j,k,l,read) result(res)
      integer(i4b) :: res
      integer(i4b), intent(in) :: i,j,k,l
      logical, intent(in) :: read

      res = 22 + i+Nl*(j-1)+Nl*Nl*(k-1)+Nl*Nl*Nl*(l-1)

      if(read) then
        res = res + 10000
      end if

    end function ind

    subroutine write_evolution_operators(type, actual_U)
      character, intent(in)      :: type
      complex(dpc), dimension(:,:,:,:,:), intent(in)     :: actual_U

      integer (i4b)       :: i
      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      do i=1,size(actual_U,5)
      do Uelement=1,N1_from_type(type)
      do Uelement2=1,N2_from_type(type)
      do Utnemele=1,N1_from_type(type)
      do Utnemele2=1,N2_from_type(type)


        write(ind(Uelement,Uelement2,Utnemele,Utnemele2,.false.) &
        ,*) timeStep*i,' ',real(actual_U(Uelement,Uelement2,Utnemele,Utnemele2,i)),' ',&
                             aimag(actual_U(Uelement,Uelement2,Utnemele,Utnemele2,i))

      end do
      end do
      end do
      end do
      end do
    end subroutine write_evolution_operators

    subroutine open_files(type)
      character, intent(in) :: type

      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2, Ublock
      character(len=4)    :: no1,no2,no3,no4
      character(len=100)  :: name
      character(len=50)   :: prefix


      Ublock = 1

      ! We set indices range according to block we evaluate. Because rho0 is
      ! whole density matrix, while evolution operators are only from particular
      ! block, offset is set between these indices.
      if (type == '2') then
        !actual_U => evops(Ublock,Ublock)%Ufe
        prefix = 'Evops_fe'
      else if (type == 'E') then
        !actual_U => evops(Ublock,Ublock)%Uee
        prefix = 'Evops_ee'
      else if (type == 'O') then
        !actual_U => evops(Ublock,Ublock)%Ueg
        prefix = 'Evops_eg'
      end if

      do Uelement=1,N1_from_type(type)
      do Uelement2=1,N2_from_type(type)
      do Utnemele=1,N1_from_type(type)
      do Utnemele2=1,N2_from_type(type)


      if(Uelement < 10) then
        write(no1,'(i1)')   Uelement
      else if (Uelement < 100) then
        write(no1,'(i2)')   Uelement
      else
        write(no1,'(i3)')   Uelement
      endif
      if(Uelement2 < 10) then
        write(no2,'(i1)')   Uelement2
      else if (Uelement2 < 100) then
        write(no2,'(i2)')   Uelement2
      else
        write(no2,'(i3)')   Uelement2
      endif
      if(Utnemele < 10) then
        write(no3,'(i1)')   Utnemele
      else if (Uelement2 < 100) then
        write(no3,'(i2)')   Utnemele
      else
        write(no3,'(i3)')   Utnemele
      endif
      if(Utnemele2 < 10) then
        write(no4,'(i1)')   Utnemele2
      else if (Uelement2 < 100) then
        write(no4,'(i2)')   Utnemele2
      else
        write(no4,'(i3)')   Utnemele2
      endif

      name = trim(prefix) // trim(no1) // '-'//trim(no2)//'--'// trim(no3) // '-'//trim(no4)//'X.dat'
      open(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,.false.), FILE = trim(name))

      end do
      end do
      end do
      end do
    end subroutine open_files

    subroutine close_files(type)
      character, intent(in) :: type
      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      do Uelement=1,N1_from_type(type)
      do Uelement2=1,N2_from_type(type)
      do Utnemele=1,N1_from_type(type)
      do Utnemele2=1,N2_from_type(type)

      close(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,.false.))

      end do
      end do
      end do
      end do
    end subroutine close_files

    subroutine random_test()
      real(dp)     :: r, rr, rrr
      integer(i4b) :: i,j

      write(*,*) 'Montecarlo'

      call init_random_seed()

      call random_number(rr)
      call random_number(rrr)

      i = 0
      j = 0
      do while (rr < 2)
         i = i + 1
         call random_number(r)

         if (r == rr) then
           write(*,*) "perioda nalezena po ", i
           stop
         end if

         if (mod(i,100000000) == 0) then
           i = 0
           j = j + 1
           write(*,*) j
         end if
      end do

      stop
    end subroutine random_test

    subroutine read_config_file(Nsys)
        integer(i4b), intent(in)     :: Nsys
        character(len=256)           :: buff = ""
        real(dp), dimension(Nsys)    :: value
        real(dp)                     :: svalue
        integer(i4b)                 :: i = 0, j

        open(unit=32,file=trim(trim(directory)//'/'//trim(config_filename) ) , err=32, status='old')

        jumps_in_one_run = 18 !default
        rwa = 0.0             !default

        do while(i == 0)
          read(32, *, iostat=i) buff

        !global_en, global_lambda, global_temp, global_tau

          if(trim(adjustl(buff)) == 'lambda') then
            read(32, *, iostat=i) value
            write(*,*) buff, value
            global_lambda(1:Nsys) = value(1:Nsys)/Energy_internal_to_cm

          elseif(trim(adjustl(buff)) == 'jumps') then
            read(32, *, iostat=i) svalue
            write(*,*) buff, int(svalue)
            jumps_in_one_run = int(svalue)

          elseif(trim(adjustl(buff)) == 'rwa') then
            read(32, *, iostat=i) svalue
            write(*,*) buff, svalue
            rwa = svalue

          elseif(trim(adjustl(buff)) == 'tauC') then
            read(32, *, iostat=i) value
            write(*,*) buff, value
            global_tau(1:Nsys) = value(1:Nsys)

           elseif(trim(adjustl(buff)) == 'temperature') then
            read(32, *, iostat=i) value
            write(*,*) buff, value
            global_temp(1:Nsys) = value(1:Nsys)

          elseif(trim(adjustl(buff)) == 'hamiltonian') then
            do j=1,Nsys
              read(32, *, iostat=i) value
              write(*,*) buff, j, value
              Ham(j,1:Nsys) = value(1:Nsys)/Energy_internal_to_cm
            end do

            buff = "-"

          else
            !write(*,*) 'unrecognised:', buff, value
          end if
        end do

        close(32)

        return

32      write(*,*) "couldn't read the supplementary config file"
        stop
    end subroutine read_config_file

end module resources_montecarlo

