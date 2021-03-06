!
!
! Fortran implementation of certain convenient
! matrix operation functions (inspired by the Scilab)
!
!
!
!
module numer_matrix
	use std_types
	use std_lapack

	implicit none

    interface inv
    	module procedure inv_real
    	module procedure inv_cmplx
    end interface

    interface spec
    	module procedure spec_real
    	module procedure spec_cmplx
    end interface

    interface spec2
    	module procedure spec2_real
    	module procedure spec2_cmplx
    end interface

    interface eigsrt
    	module procedure eigsrt_real
    	module procedure eigsrt_cmplx
    	module procedure eigsrt_real_mat
    	module procedure eigsrt_cmplx_mat
    end interface

    interface eigconvention
    	module procedure eigconvention_real
    	module procedure eigconvention_cmplx
    end interface

    interface matrix_power
    	module procedure matrix_power_real
    	module procedure matrix_power_real2
    	module procedure matrix_power_real3
    	module procedure matrix_power_cmplx
    	module procedure matrix_power_cmplx2
    	module procedure matrix_power_cmplx3
    end interface

    interface matrix_exp
    	module procedure matrix_exp_real
    	module procedure matrix_exp_cmplx
    end interface

    interface spec_generalized
    	module procedure spec_generalized_real
    	module procedure spec_generalized_cmplx
    end interface

    interface entropy
        module procedure entropy_real
        module procedure entropy_cmplx
    end interface

    interface cholesky
        module procedure cholesky_real
        module procedure cholesky_cmplx
    end interface

    public :: spec, spec2, inv, matrix_power, matrix_exp

	private :: eigsrt, eigconvention, iminloc
	private :: inv_real, inv_cmplx, spec_real, spec_cmplx, spec2_real, spec2_cmplx
	private :: eigsrt_real, eigsrt_cmplx, eigconvention_real, eigconvention_cmplx
	private :: matrix_power_real, matrix_power_real2, matrix_power_real3
	private :: matrix_power_cmplx, matrix_power_cmplx2, matrix_power_cmplx3
	private :: matrix_exp_real !, matrix_exp_real2, matrix_exp_real3
	private :: matrix_exp_cmplx !, matrix_exp_cmplx2, matrix_exp_cmplx3
	private :: spec_generalized_real, spec_generalized_cmplx, eigsrt_real_mat, eigsrt_cmplx_mat
	private :: entropy_real, entropy_cmplx


contains

  !
  ! Eigenvalue problem
  !
  subroutine spec_real(A,S1,WM)

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)), intent(out) :: S1
    real(dp), dimension(size(A,1),size(A,2)), intent(out) :: WM

    ! local
    real(dp), dimension(size(A,1)) :: W,W2
    real(dp), dimension(size(A,1),size(A,2)):: AC
    integer(i4b) :: N,i
	integer :: info

    N = size(A,1)
    AC = A
	WM = 0.0_dp

    if (N.ne.size(A,2)) then
       stop "error in spec: square matrix required on input"
    end if

    call std_lapack_geev(AC,W,W2,VR=S1,INFO=info)

    if (info.ne.0) then
       write(*,*) "error in spec: lapack error", info
       stop
    end if

    call eigsrt(W,S1)
    call eigconvention(S1)

    do i = 1, N
	WM(i,i) = W(i)
    end do


  end subroutine spec_real

  subroutine spec_cmplx(A,S1,WM)

    complex(dpc), dimension(:,:), intent(in) :: A
    complex(dpc), dimension(size(A,1),size(A,2)), intent(out) :: S1
    complex(dpc), dimension(size(A,1),size(A,2)), intent(out) :: WM

    ! local
    complex(dpc), dimension(size(A,1)) :: W,W2
    complex(dpc), dimension(size(A,1),size(A,2)):: AC
    integer(i4b) :: N,i
	integer :: info

    N = size(A,1)
    AC = A
	WM = 0.0_dp

    if (N.ne.size(A,2)) then
       stop "error in spec: square matrix required on input"
    end if

	! W2 not used - only for the same number af args
    call std_lapack_geev(AC,W,W2,VR=S1,INFO=info)

    if (info.ne.0) then
       write(*,*) "error in spec: lapack error", info
       stop
    end if

    call eigsrt(W,S1)
!    call eigconvention(S1)

    do i = 1, N
	WM(i,i) = W(i)
    end do


  end subroutine spec_cmplx

  !
  ! Eigenvalues without sorting the eigenvalues
  !
  subroutine spec2_real(A,S1,WM)

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)), intent(out) :: S1
    real(dp), dimension(size(A,1),size(A,2)), intent(out) :: WM

    ! local
    real(dp), dimension(size(A,1)) :: W,W2
    real(dp), dimension(size(A,1),size(A,2)):: AC
    integer(i4b) :: N, i
    integer :: info

    N = size(A,1)
    AC = A

    if (N.ne.size(A,2)) then
       stop "error in spec: square matrix required on input"
    end if

    call std_lapack_geev(AC,W,W2,VR=S1,INFO=info)

    if (info.ne.0) then
       write(*,*) "error in spec: lapack error", info
       stop
    end if

    do i = 1, N
	WM(i,i) = W(i)
    end do


  end subroutine

  subroutine spec2_cmplx(A,S1,WM)

    complex(dpc), dimension(:,:), intent(in) :: A
    complex(dpc), dimension(size(A,1),size(A,2)), intent(out) :: S1
    complex(dpc), dimension(size(A,1),size(A,2)), intent(out) :: WM

    ! local
    complex(dpc), dimension(size(A,1)) :: W,W2
    complex(dpc), dimension(size(A,1),size(A,2)):: AC
    integer(i4b) :: N, i
    integer :: info

    N = size(A,1)
    AC = A

    if (N.ne.size(A,2)) then
       stop "error in spec: square matrix required on input"
    end if

    call std_lapack_geev(AC,W,W2,VR=S1,INFO=info)

    if (info.ne.0) then
       write(*,*) "error in spec: lapack error", info
       stop
    end if

    do i = 1, N
	WM(i,i) = W(i)
    end do


  end subroutine



  !
  ! Matrix inversion
  !
  subroutine inv_real(A,Ai)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)), intent(out) :: Ai

    ! local
    integer(i4b), dimension(size(A,1)) :: ipiv
    real(dp), dimension(size(A,1),size(A,2)) :: B
    integer(i4b) :: i, N
    integer :: info

    N = size(A,1)
    if (N.ne.size(A,2)) then
       stop "error in inv: square matrix required on input"
    end if

    Ai = A
    B = 0.0_dp
    do i=1,N
	B(i,i) = 1.0_dp
    end do

    call std_lapack_gesv(Ai,B,INFO=info)

    if (info.ne.0) then
       stop "error in inv: lapack error"
    end if

    Ai = B

  end subroutine inv_real

  !
  ! Matrix inversion (complex)
  !
  subroutine inv_cmplx(A,Ai)
    complex(dpc), dimension(:,:), intent(in) :: A
    complex(dpc), dimension(size(A,1),size(A,2)), intent(out) :: Ai

    ! local
    integer(i4b), dimension(size(A,1)) :: ipiv
    complex(dpc), dimension(size(A,1),size(A,2)) :: B
    integer(i4b) :: i, N
    integer :: info

    N = size(A,1)
    if (N.ne.size(A,2)) then
       stop "error in inv: square matrix required on input"
    end if

    Ai = A
    B = 0.0_dp
    do i=1,N
	B(i,i) = 1.0_dp
    end do

    call std_lapack_gesv(Ai,B,INFO=info)

    if (info.ne.0) then
       stop "error in inv: lapack error"
    end if

    Ai = B

  end subroutine inv_cmplx

  !
  ! Sorting eigenvalues
  !
  subroutine eigsrt_real(d,v)
  	real(dp), dimension(:), intent(inout)   :: d
  	real(dp), dimension(:,:), intent(inout) :: v
  	real(dp) :: dsw
  	real(dp), dimension(size(d,1)) :: vsw

  	integer(i4b) :: i,j,k,n
  	n = size(d,1)
  	do i =  1,n-1
  		j = iminloc(d(i:n))+i-1
  		if (j /= i) then
  			dsw = d(i)
  			d(i) = d(j)
  			d(j) = dsw
   			vsw(:) = v(:,i)
   			v(:,i) = v(:,j)
   			v(:,j) = vsw(:)
  		end if
  	end do
  end subroutine eigsrt_real

  subroutine eigsrt_cmplx(d,v)
  	complex(dp), dimension(:), intent(inout)   :: d
  	complex(dpc), dimension(:,:), intent(inout) :: v
  	complex(dp) :: dsw
  	complex(dpc), dimension(size(d,1)) :: vsw

  	integer(i4b) :: i,j,k,n
  	n = size(d,1)
  	do i =  1,n-1
  		j = iminloc(aimag(d(i:n)))+i-1
  		if (j /= i) then
  			dsw = d(i)
  			d(i) = d(j)
  			d(j) = dsw
   			vsw(:) = v(:,i)
   			v(:,i) = v(:,j)
   			v(:,j) = vsw(:)
  		end if
  	end do


  	do i =  1,n-1
  		j = iminloc(real(d(i:n)))+i-1
  		if (j /= i) then
  			dsw = d(i)
  			d(i) = d(j)
  			d(j) = dsw
   			vsw(:) = v(:,i)
   			v(:,i) = v(:,j)
   			v(:,j) = vsw(:)
  		end if
  	end do
  end subroutine eigsrt_cmplx

  subroutine eigsrt_real_mat(d,v)
  	real(dp), dimension(:,:), intent(inout)   :: d
  	real(dp), dimension(:,:), intent(inout) 	:: v

  	real(dp), dimension(size(d,1)) :: dd
  	integer(i4b) :: i

  	do i=1,size(d,1)
  		dd(i) = d(i,i)
  	end do

  	call eigsrt(dd,v)

	d = 0.0_dp
  	do i=1,size(d,1)
  		d(i,i) = dd(i)
  	end do

  end subroutine eigsrt_real_mat

  subroutine eigsrt_cmplx_mat(d,v)
  	complex(dpc), dimension(:,:), intent(inout)   :: d
  	complex(dpc), dimension(:,:), intent(inout) 	:: v

  	complex(dpc), dimension(size(d,1)) :: dd
  	integer(i4b) :: i

  	do i=1,size(d,1)
  		dd(i) = d(i,i)
  	end do

  	call eigsrt(dd,v)

	d = 0.0_dp
  	do i=1,size(d,1)
  		d(i,i) = dd(i)
  	end do

  end subroutine eigsrt_cmplx_mat

  !
  ! Gives unique sign convention to normalized eigenvectors,
  ! i.e. real part of the maximum abs() component must be positive
  !
  subroutine eigconvention_real(v)
  	real(dp), dimension(:,:), intent(inout) 	:: v
  	real(dp), dimension(size(v(:,1),1))		:: vect
  	real(dp)									:: maximum
  	integer										:: i,k,j

	j = 1
	do while (j <= size(v(1,:)))
		vect =v(:,j)
  		i = 1
  		k = 0
  		maximum = -1
  		do while (i <= size(vect))
  	  		if (abs(vect(i)) > maximum ) then
  	  			maximum = abs(vect(i))
  	  			k = i
  	  		end if
  	  		i = i + 1
  		end do
  	  	v(:,j) = vect * vect(k) / maximum

  		j = j + 1
  	end do

  end subroutine eigconvention_real

    subroutine eigconvention_cmplx(v)
  	complex(dpc), dimension(:,:), intent(inout) 	:: v
  	complex(dpc), dimension(size(v(:,1),1))		:: vect
  	real(dp)									:: maximum
  	integer										:: i,k,j

	j = 1
	do while (j <= size(v(1,:)))
		vect =v(:,j)
  		i = 1
  		k = 0
  		maximum = -1
  		do while (i <= size(vect))
  	  		if (abs(vect(i)) > maximum ) then
  	  			maximum = abs(vect(i))
  	  			k = i
  	  		end if
  	  		i = i + 1
  		end do
  	  	v(:,j) = vect * vect(k) / maximum

  		j = j + 1
  	end do

  end subroutine eigconvention_cmplx


  !
  ! location of the minimum
  !
  function iminloc(ar) result(k)
  	integer(i4b) :: k
  	real(dp), dimension(:), intent(in) :: ar
  	integer(i4b), dimension(1) :: imax
  	imax = minloc(ar)
  	k = imax(1)
  end function iminloc

  !
  ! general power of the matrix
  !
  subroutine matrix_power_real(A,r)
  	real(dp), dimension(:,:), intent(inout) 	:: A
  	real(dp), intent(in)						:: r

  	real(dp), dimension(size(A,1),size(A,2))	:: S, invS, AA
  	integer(i4b)								:: i

	call spec(A,S,AA)
	call inv(S,invS)

	A = matmul(matmul(invS,A),S)
	do i=1,size(A,1)
		A(i,i) = A(i,i)**r
	end do
	A = matmul(matmul(S,A),invS)

  end subroutine matrix_power_real

  subroutine matrix_power_real2(A,r)
  	real(dp), dimension(:,:), intent(inout) 	:: A
  	real(sp), intent(in)						:: r

	call matrix_power_real(A, real(r,dp))

  end subroutine matrix_power_real2

  subroutine matrix_power_real3(A,r)
  	real(dp), dimension(:,:), intent(inout) 	:: A
  	integer(i4b), intent(in)					:: r

	call matrix_power_real(A, real(r,dp))

  end subroutine matrix_power_real3

  subroutine matrix_power_cmplx(A,r)
  	complex(dpc), dimension(:,:), intent(inout) 	:: A
  	real(dp), intent(in)							:: r

  	complex(dpc), dimension(size(A,1),size(A,2))	:: S, invS, AA
  	integer(i4b)								:: i

	call spec(A,S,AA)
	call inv(S,invS)

	A = matmul(matmul(invS,A),S)
	do i=1,size(A,1)
		A(i,i) = A(i,i)**r
	end do
	A = matmul(matmul(S,A),invS)

  end subroutine matrix_power_cmplx

  subroutine matrix_power_cmplx2(A,r)
  	complex(dpc), dimension(:,:), intent(inout) 	:: A
  	real(sp), intent(in)							:: r

	call matrix_power_cmplx(A, real(r,dp))

  end subroutine matrix_power_cmplx2

  subroutine matrix_power_cmplx3(A,r)
  	complex(dpc), dimension(:,:), intent(inout) 	:: A
  	integer(i4b), intent(in)						:: r

	call matrix_power_cmplx(A, real(r,dp))

  end subroutine matrix_power_cmplx3

  !
  ! matrix exponential
  !
  subroutine matrix_exp_real(A)
  	real(dp), dimension(:,:), intent(inout) 	:: A

  	real(dp), dimension(size(A,1),size(A,2))	:: S, invS, AA
  	integer(i4b)								:: i

	call spec(A,S,AA)
	call inv(S,invS)

	A = matmul(matmul(invS,A),S)
	do i=1,size(A,1)
		A(i,i) = exp(A(i,i))
	end do
	A = matmul(matmul(S,A),invS)

  end subroutine matrix_exp_real

  subroutine matrix_exp_cmplx(A)
  	complex(dpc), dimension(:,:), intent(inout) 	:: A

  	complex(dpc), dimension(size(A,1),size(A,2))	:: S, invS, AA
  	integer(i4b)								:: i

	call spec(A,S,AA)
	call inv(S,invS)

	A = matmul(matmul(invS,A),S)
	do i=1,size(A,1)
		A(i,i) = exp(A(i,i))
	end do
	A = matmul(matmul(S,A),invS)

  end subroutine matrix_exp_cmplx

  !
  ! eigenvalue problem with overlap matrix
  !
  subroutine spec_generalized_real(Ain,Sin,S1,WM)
    real(dp), dimension(:,:), intent(in) :: Ain,Sin
    real(dp), dimension(size(Ain,1),size(Ain,2)), intent(out) :: S1
    real(dp), dimension(size(Ain,1),size(Ain,2)), intent(out) :: WM

    real(dp), dimension(size(Ain,1),size(Ain,2)) 	:: Shalf,Shalf_inv, F
    integer(i4b)									:: i

	Shalf = Sin
    call matrix_power(Shalf,0.5_dp)
    call inv(Shalf,Shalf_inv)

    F = matmul(matmul(Shalf_inv,Ain),Shalf_inv)
    call spec(F,S1,WM)

	do i=1,size(Ain,1)
		S1(:,i) = matmul(Shalf_inv,S1(:,i))
	end do

	call eigsrt(WM,S1)

	! normalize
	do i=1,size(Ain,1)
		S1(:,i) = S1(:,i)/sqrt( dot_product(S1(:,i),matmul(Sin,S1(:,i))) )
	end do

  end subroutine spec_generalized_real

  subroutine spec_generalized_cmplx(Ain,Sin,S1,WM)
    complex(dpc), dimension(:,:), intent(in) :: Ain,Sin
    complex(dpc), dimension(size(Ain,1),size(Ain,2)), intent(out) :: S1
    complex(dpc), dimension(size(Ain,1),size(Ain,2)), intent(out) :: WM

    complex(dpc), dimension(size(Ain,1),size(Ain,2)) 	:: Shalf,Shalf_inv, F
    integer(i4b)									:: i

	Shalf = Sin
    call matrix_power(Shalf,0.5_dp)
    call inv(Shalf,Shalf_inv)

    F = matmul(matmul(Shalf,Ain),Shalf_inv)
    call spec(F,S1,WM)

	do i=1,size(Ain,1)
		S1(:,i) = matmul(Shalf_inv,S1(:,i))
	end do

	call eigsrt(WM,S1)

	! normalize
	do i=1,size(Ain,1)
		S1(:,i) = S1(:,i)/sqrt( dot_product(S1(:,i),matmul(Sin,S1(:,i))) )
	end do

  end subroutine spec_generalized_cmplx

  !
  ! matrix exponential
  !
  real(dp) function entropy_real(AAA) result(res)
    real(dp), dimension(:,:), intent(in)          :: AAA

    real(dp), dimension(size(AAA,1),size(AAA,2))  :: S, invS, AA, A
    integer(i4b)                                  :: i

    A = AAA
    res = 0.0_dp

    call spec(AAA,S,AA)
    call inv(S,invS)

    A = matmul(matmul(invS,A),S)
    do i=1,size(A,1)
        if(A(i,i) > 0 .and. A(i,i) < 1) then
            res = res - A(i,i)*log(A(i,i))
        end if
    end do

  end function entropy_real

  real(dp) function entropy_cmplx(AAA) result(res)
    complex(dpc), dimension(:,:), intent(in)          :: AAA

    complex(dpc), dimension(size(AAA,1),size(AAA,2))  :: S, invS, AA, A
    integer(i4b)                                      :: i

    A = AAA
    res = 0.0_dp

    call spec(AAA,S,AA)
    call inv(S,invS)

    A = matmul(matmul(invS,A),S)
    do i=1,size(A,1)
        if(real(A(i,i)) > 0 .and. real(A(i,i)) < 1) then
            res = res - real(A(i,i))*log(real(A(i,i)))
        end if
    end do

  end function entropy_cmplx

  subroutine cholesky_real(A,U)
    real(dp), dimension(:,:), intent(in)              :: A
    real(dp), dimension(:,:), intent(out)             :: U

    integer(i4b)                                      :: info, i, j, kd
    real(dp), dimension(size(A,1),size(A,1))          :: X

    kd = size(A,1)-1

    X = 0.0_dp

    do j=1,size(A,1)
    do i=max(1,j-kd),j
      X(kd+1+i-j,j) = A(i,j)
    end do
    end do

    !write(*,*) X
    !write(*,*)

    call DPBTRF( 'U', size(A,1), kd, X, size(A,1), info )

    !write(*,*) X

    U = 0.0_dp

    do j=1,size(A,1)
    do i=max(1,j-kd),j
      U(i,j) = X(kd+1+i-j,j)
    end do
    end do

    write(*,*) 'cholesky info', info
  end subroutine cholesky_real

  subroutine cholesky_cmplx(A,U)
    complex(dpc), dimension(:,:), intent(in)          :: A
    complex(dpc), dimension(:,:), intent(out)         :: U

    integer(i4b)                                      :: info, i, j, kd
    complex(dpc), dimension(size(A,1),size(A,1))      :: X

    kd = size(A,1)-1

    X = 0.0_dp

    do j=1,size(A,1)
    do i=max(1,j-kd),j
      X(kd+1+i-j,j) = A(i,j)
    end do
    end do

    !write(*,*) X
    !write(*,*)

    call ZPBTRF( 'U', size(A,1), kd, X, size(A,1), info )

    !write(*,*) X

    U = 0.0_dp

    do j=1,size(A,1)
    do i=max(1,j-kd),j
      U(i,j) = X(kd+1+i-j,j)
    end do
    end do

    write(*,*) 'cholesky info', info

  end subroutine cholesky_cmplx

end module numer_matrix
