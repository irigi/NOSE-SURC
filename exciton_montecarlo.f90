
module exciton_montecarlo
    use std_types
    use helpers
    use nakajima_zwanzig_shared
    use module_goft
    use basic_montecarlo

    implicit none

    complex(dpc), dimension(:,:,:,:,:), allocatable   :: gg_MC, hh_MC, cc_MC

    logical, private ::  origin = .false., g_functions = .true.,                               &
                                       modified_unraveling2 = .false.,                         &
                                       vynechani_G_ifu = .false.

    interface goft_general
        module procedure goft_general_2
        module procedure goft_general_4
    end interface
    interface dgoft_general
        module procedure dgoft_general_2
        module procedure dgoft_general_4
    end interface

    contains

    !
    ! Do all the simulations within the module
    !
    subroutine do_montecarlo_work_exciton()
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

            call perform_montecarlo_parallel_exciton(r,s,rho,rho_coherent)

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

        if(only_coherences) then
            call open_files('O')
            call write_evolution_operators('O',Ueg)
            call close_files('O')
        end if

    end subroutine do_montecarlo_work_exciton

    subroutine perform_montecarlo_parallel_exciton(i0, j0, rho, rho_coherent)
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
          call perform_montecarlo_exciton(i0, j0, rho_local, rho_coherent_local, TRAJECTORIES/360)

          rho = rho + rho_local
          rho_coherent = rho_coherent + rho_coherent_local

          call itime(time)
          write(buff, '(A I2 A I2.2 A I2.2 A I6 )') '  time ',time(1), ':', time(2), ':', time(3) ,  &
                 ', STEPS/360: ', i
          call print_log_message(trim(buff), 5)
        end do

    end subroutine perform_montecarlo_parallel_exciton

    subroutine perform_montecarlo_exciton(i0, j0, rho, rho_coherent, traj)
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
                call calculate_Gfactor_from_trajectory_history_general_basis(draha,cmplx(1,0,dp),Gfactor)
            end if

            !call calculate_Ifactor_from_trajectory_history(draha,cmplx(1,0,dp),Ifactor)
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

    end subroutine perform_montecarlo_exciton

     ! Let we have U^+_A(t1) V_AB U^+_B(t2) V_BC U^+_C(t3)...U_c(t3)  V_cb  U_b(t2) V_ba U_a(t1)W,
     ! then the cummulant expansion is derived from it starts as
     ! (1-i h_A(t1)-g_A(t1))^* V_AB(t1) (1-i h_B(t1)-g_B(t1))(1-i h_B(t1+t2)-g_B(t1+t2))^*  V_BC(t1+t2)
     ! (1-i h_C(t1+t2)-g_C(t1+t2)) ... (1-i h_Z(t1+..+tn)-g_Z(t1+..+tn))^* (1-i h_z(t1+..+tn)-g_z(t1+..+tn))
     ! (1-i h_z(t1+..+t(n-1))-g_z(t1+..+t(n-1)))^* V_zy(t1+..+t(n-1) ) (1-i h_y(t1+..+t(n-1))-g_y(t1+..+t(n-1)))
    subroutine calculate_Gfactor_from_trajectory_history_general_basis(trajectory,initial_factor,factor_out)
        integer(i1b), dimension(:,:), intent(in) :: trajectory
        complex(dpc), intent(in)                 :: initial_factor
        complex(dpc), dimension(:), intent(out)  :: factor_out

        ! we overallocate jump-fields not to have to always calculate their dimensions
        integer(i4b), parameter :: overallocation = 10000
        integer(i4b), dimension(0:overallocation)   ::  &
                            jumps_index, jumps_tovalue_ket, jumps_tovalue_bra
        real(dp), dimension(0:overallocation)       ::  &
                            jumps_time
        logical, dimension(0:overallocation)        ::  &
                            jumps_ket

        integer(i4b) ::     i, j, k, number_of_jumps, last_index, tovalue,      &
                            fromvalue, ket_i, bra_i, ket_j, bra_j,  ket_i_, bra_i_, ket_j_, bra_j_

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
        jumps_index(0)              = 1
        jumps_time(0)            = 0.0_dp
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
                    ket_i_ = jumps_tovalue_ket(i-1)
                    bra_i_ = jumps_tovalue_bra(i-1)

                    if(lze_pouzit_AFNLTZM) then
                        additive_factor = additive_factor_not_last_turn_z_minuleho
                        cycle ! we wait until we are in the last round
                    end if
                end if

                ! there are always two 1hg-brackets for both bra and ket for each jump-time, with site i and (i-1).
                ! exception is the last round, where the bra and ket meet into two brackets instead of four
                !
                ! The general ordering is ket_, ket, bra, bra_
                if(.not. last_round_i) then
                    ! two ket-brackets
                    additive_factor = additive_factor+(-goft_general(ket_i,ket_i,time_i))
                        if(debug_G) then
                            write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if
                    additive_factor = additive_factor+conjg((-goft_general(ket_i_,ket_i_,time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i, ket_i_
                        end if

                    ! two bra-brackets
                    if(.not. only_coherences) then
                        additive_factor = additive_factor+conjg((-goft_general(bra_i,bra_i,time_i)))
                            if(debug_G) then
                                write(*,*) 'factor_g', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                        additive_factor = additive_factor+(-goft_general(bra_i_,bra_i_,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg', real(additive_factor), k, i, time, time_i, last_round_i,jumps_tovalue_bra(i-1)
                            end if
                    end if

                    ! hh between bra and ket brackets of the same i (4 terms: C(4,2) => 6 pairs)
                        ! ket_i == ket_i_
                        additive_factor = additive_factor + (hoft_general(ket_i_,ket_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            end if

                    if(.not. only_coherences) then
                        ! ket_i == bra_i
                        additive_factor = additive_factor + (hoft_general(ket_i,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H1', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, bra_i
                            end if

                        ! ket_i == bra_i_
                        additive_factor = additive_factor + (-hoft_general(ket_i,bra_i_,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H2', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, bra_i_
                            end if

                        ! ket_i_ == bra_i
                        additive_factor = additive_factor + (-hoft_general(ket_i_,bra_i,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H3', real(additive_factor), k, i, time, time_i, last_round_i, ket_i_, bra_i
                            end if

                        ! ket_i_ == bra_i_
                        additive_factor = additive_factor + (hoft_general(ket_i_,bra_i_,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H4', real(additive_factor), k, i, time, time_i, last_round_i, ket_i_, bra_i_
                            end if

                        ! bra_i == bra_i_
                        additive_factor = additive_factor + (hoft_general(bra_i,bra_i_,time_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_H6', real(additive_factor), k, i, time, time_i, last_round_i, bra_i, bra_i
                            end if
                    endif

                    ! hV, Vh between bra and ket brackets of the same i (4 h-terms, 2 V-terms => 8 pairs)
                        ! ket_i == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   ket_i,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (hV_general(ket_i_,   ket_i_,ket_i,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    if(.not. only_coherences) then
                        ! bra_i == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_i,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(bra_i,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (Vh_general(bra_i,bra_i_,   bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_i,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_i_,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if



                    ! VV term
                        additive_factor = additive_factor + (VV_general(ket_i_,ket_i,   bra_i,bra_i_,time_i,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    end if


                ! one additional factor where bra and kets meet
                else ! not last round
                    additive_factor = additive_factor+conjg((-goft_general(ket_i,ket_i,time_i)))
                        if(debug_G) then
                            write(*,*) 'factor_g--', real(additive_factor), k, i, time, time_i, last_round_i, ket_i
                        end if

                    if(.not. only_coherences) then
                        additive_factor = additive_factor+(-goft_general(bra_i,bra_i,time_i))
                            if(debug_G) then
                                write(*,*) 'factor_gg--', real(additive_factor), k, i, time, time_i, last_round_i, bra_i
                            end if
                    end if

                    ! bra_i == ket_i
                    if(.not. only_coherences) then
                        additive_factor = additive_factor + (hoft_general(ket_i,bra_i,time_i,time_i))

                            if(debug_G) then
                                write(*,*) 'factor_h (last)', real(additive_factor), k, i, time, time_i, last_round_i
                            end if
                    end if

                    ! no V-terms here, I guess
                end if


                do j=1, i-1
                    time_j = jumps_time(j)
                    ket_j = jumps_tovalue_ket(j)
                    bra_j = jumps_tovalue_bra(j)
                    ket_j_ = jumps_tovalue_ket(j-1)
                    bra_j_ = jumps_tovalue_bra(j-1)

                    ! hh between all combinations of i-j pairs (8 terms: C(8,2) - 2*C(4,2)  => 16 pairs)
                    ! the correct ordering is:  ket_j, ket_j_, ket_i, ket_i_, bra_i_, bra_i, bra_j_, bra_j
                        ! ket_i_ == ket_j
                        additive_factor = additive_factor + (hoft_general(ket_j,ket_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor2', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                        ! ket_i_ == ket_j_
                        additive_factor = additive_factor + (-hoft_general(ket_j_,ket_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor6', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                    if(.not. only_coherences) then
                        ! bra_i_ == ket_j_
                        additive_factor = additive_factor + (hoft_general(ket_j_,bra_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor8', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i_ == ket_j
                        additive_factor = additive_factor + (-hoft_general(ket_j,bra_i_,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor4', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i_ == bra_j_
                        additive_factor = additive_factor + (hoft_general(ket_i_,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor10', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i_ == bra_j_
                        additive_factor = additive_factor + (-hoft_general(bra_i_,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor12', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i_ == bra_j
                        additive_factor = additive_factor + (-hoft_general(ket_i_,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor14', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i_ == bra_j
                        additive_factor = additive_factor + (hoft_general(bra_i_,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor16', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                    endif
                    if(.not. last_round_i) then
                        ! ket_i == ket_j
                        additive_factor = additive_factor + (-hoft_general(ket_j,ket_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor1', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i == ket_j_
                        additive_factor = additive_factor + (hoft_general(ket_j_,ket_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor5', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif
                    if(.not. only_coherences .and. .not. last_round_i) then
                        ! bra_i == ket_j
                        additive_factor = additive_factor + (hoft_general(ket_j,bra_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor3', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i == ket_j_
                        additive_factor = additive_factor + (-hoft_general(ket_j_,bra_i,time_j,time_i))
                            if(debug_G) then
                                write(*,*) 'factor7', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i == bra_j_
                        additive_factor = additive_factor + (-hoft_general(ket_i,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor9', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i == bra_j_
                        additive_factor = additive_factor + (hoft_general(bra_i,bra_j_,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor11', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! ket_i == bra_j
                        additive_factor = additive_factor + (hoft_general(ket_i,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor13', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if

                        ! bra_i == bra_j
                        additive_factor = additive_factor + (-hoft_general(bra_i,bra_j,time_i,time_j))
                            if(debug_G) then
                                write(*,*) 'factor15', real(additive_factor), k, i, j, time, time_i, time_j, last_round_i
                            end if
                    endif


                    ! hV, Vh between bra and ket brackets of the same i (8 h-terms, 4 V-terms, 8 x 4 - 2*8 => 16 pairs)
                        ! ket_i == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   ket_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (hV_general(ket_j,   ket_i_,ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (hV_general(ket_j_,   ket_i_,ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    if(.not. only_coherences) then
                        ! bra_j == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_j,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_j_ == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (Vh_general(ket_i_,ket_i,   bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_j == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (Vh_general(bra_i,bra_i_,   bra_j,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_j_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (Vh_general(bra_i,bra_i_,   bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_j,   bra_i,bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_j_ == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (hV_general(ket_j_,   bra_i,bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   bra_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(ket_j_,ket_j)
                        additive_factor = additive_factor + (Vh_general(ket_j_,ket_j,   bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(bra_i,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! bra_i_ == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(bra_i_,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(ket_i,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! ket_i_ == V(bra_j,bra_j_)
                        additive_factor = additive_factor + (hV_general(ket_i_,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if



                    ! VV terms C(4,2) - 2 = 4
                        ! V(ket_j_,ket_j) == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (VV_general(ket_j_,ket_j,   ket_i_,ket_i,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! V(bra_j,bra_j_) == V(ket_i_,ket_i)
                        additive_factor = additive_factor + (VV_general(ket_i_,ket_i,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! V(ket_j_,ket_j) == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (VV_general(ket_j_,ket_j,   bra_i,bra_i_,time_j,time_i))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                        ! V(bra_j,bra_j_) == V(bra_i,bra_i_)
                        additive_factor = additive_factor + (VV_general(bra_i,bra_i_,   bra_j,bra_j_,time_i,time_j))
                            !if(debug_G) then
                            !    write(*,*) 'factor_H5', real(additive_factor), k, i, time, time_i, last_round_i, ket_i, ket_i
                            !end if

                    end if

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

    end subroutine calculate_Gfactor_from_trajectory_history_general_basis


    ! this is basicaly ordinary correlation function in general basis, so
    ! it has four idices, since it is quantity of second order in SBC
    recursive function VV_general(m,n,   i,j,t1,t2) result(res)
        real(dp), intent(in)        :: t1,t2
        integer(i4b), intent(in)    :: m,n,i,j
        complex(dpc)                :: res

        integer(i4b)                :: t_index
        real(dp)                    :: t

        t = t1 - t2

        if(t < 0) then
            res = conjg(VV_general(i,j,m,n,t2,t1))
            return
        end if

        t_index = INT(t/timeStep)+1

        !!! what about terms ii-jj, are they in this term too, or are they in g-fctions of excitons?

        res = cc_MC(m,n,i,j,t_index)

    end function VV_general

    function goft_general_2(m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: res

        res = goft_general_4(m,m,n,n,t)
    end function goft_general_2

    recursive function goft_general_4(i,j,m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: i,j,m,n
        complex(dpc)                :: res
        integer (i4b)               :: a, t_index
        character(len=300)          :: buff

        if(t < 0) then
            res = conjg(goft_general_4(m,n,i,j,-t))
            return
        end if

        t_index = INT(t/timeStep)+1

        res = 0.0_dp

!        if(m < 1 .or. m > size(iblocks(1,1)%sblock%gindex) .or. n < 1 .or. n > size(iblocks(1,1)%sblock%gindex)) then
!            write(buff,'(i2)') m
!            buff = "m="//trim(buff)//" exceeds Nl1 size in goft_site()"
!            call print_error_message(-1,buff)
!            stop
!        end if
!
!        if(t_index < 1 .or. t_index > size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)) then
!            write(*,*) t, t_index, m,  size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)
!            call print_error_message(-1,"tau exceeds goft size in dgoft_site()")
!            stop
!        end if

        res = gg_MC(i,j,m,n,t_index)

    end function goft_general_4

    function dgoft_general_2(m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: res

        res = dgoft_general_4(m,m,n,n,t)
    end function dgoft_general_2

    recursive function dgoft_general_4(i,j,m,n,t) result(res)
        real(dp), intent(in)        :: t
        integer(i4b), intent(in)    :: m,n,i,j
        complex(dpc)                :: res
        integer (i4b)               :: a, t_index
        character(len=300)          :: buff

        if(t < 0) then
            res = -conjg(dgoft_general_4(m,n,i,j,-t))
            return
        end if

        t_index = INT(t/timeStep)+1

        res = 0.0_dp

!        if(m < 1 .or. m > size(iblocks(1,1)%sblock%gindex) .or. n < 1 .or. n > size(iblocks(1,1)%sblock%gindex)) then
!            write(buff,'(i2)') m
!            buff = "m="//trim(buff)//" exceeds Nl1 size in goft_site()"
!            call print_error_message(-1,buff)
!            stop
!        end if
!
!        if(t_index < 1 .or. t_index > size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)) then
!            write(*,*) t, t_index, m,  size(all_goft(iblocks(1,1)%sblock%gindex(m))%gg,1)
!            call print_error_message(-1,"tau exceeds goft size in dgoft_site()")
!            stop
!        end if

        res = hh_MC(i,j,m,n,t_index)

    end function dgoft_general_4

    function hoft_general(m,n,t1,t2) result(res)
    ! hoft(a,a,t1,t2) = hoft(a,a,-t2,-t1)
        real(dp), intent(in)        :: t1,t2
        integer(i4b), intent(in)    :: m,n
        complex(dpc)                :: res

        res = 0.0_dp

        res = goft_general(m,n,t1) - goft_general(m,n,t1-t2) + goft_general(m,n,-t2)

    end function hoft_general

    recursive function hV_general(i,   m,n,th,tV) result(res)
        real(dp), intent(in)        :: th,tV
        integer(i4b), intent(in)    :: m,n,i
        complex(dpc)                :: res

        res = 0.0_dp

        res = goft_general(i,i,m,n,th) - dgoft_general(i,i,m,n,-tV) + dgoft_general(i,i,m,n,th-tV)

    end function hV_general

    recursive function Vh_general(m,n,   i,tV,th) result(res)
        real(dp), intent(in)        :: th,tV
        integer(i4b), intent(in)    :: m,n,i
        complex(dpc)                :: res

        res = 0.0_dp

        res = dgoft_general(m,n,i,i,tV) + goft_general(m,n,i,i,-th) - dgoft_general(m,n,i,i,tV-th)

    end function Vh_general

    subroutine deinit_goft_general()

        deallocate(gg_MC)
        deallocate(hh_MC)
        deallocate(cc_MC)

    end subroutine deinit_goft_general

    subroutine init_goft_general()
        integer(i4b) :: i,j,k,l,r,s,t_index

        call init_nakajima_zwanzig_shared (Ham)

        allocate(gg_MC(Nl, Nl, Nl, Nl, size(goft,2)) )
        allocate(hh_MC(Nl, Nl, Nl, Nl, size(goft,2)) )
        allocate(cc_MC(Nl, Nl, Nl, Nl, size(goft,2)) )

        gg_MC = 0.0_dp
        hh_MC = 0.0_dp
        cc_MC = 0.0_dp

        ! init of exciton goft and hoft
        do t_index=1,size(goft,2)

        do i=1,Nl
        do j=1,Nl
        do k=1,Nl
        do l=1,Nl

        do r=1,Nl
        do s=1,Nl
            gg_MC(i,j,k,l,t_index) =     S1(i,r) * goft(r, t_index) * SS(r,j) &
                                            * S1(k,s) * goft(s, t_index) * SS(s,l)

            hh_MC(i,j,k,l,t_index) =     S1(i,r) * hoft(r, t_index) * SS(r,j) &
                                            * S1(k,s) * hoft(s, t_index) * SS(s,l)

            cc_MC(i,j,k,l,t_index) =     S1(i,r) * coft(r, t_index) * SS(r,j) &
                                            * S1(k,s) * coft(s, t_index) * SS(s,l)
        end do
        end do

        end do
        end do
        end do
        end do

        end do

    end subroutine init_goft_general

end module exciton_montecarlo

