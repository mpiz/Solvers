	subroutine umnog(diag,ggu,ggl,ig,jg,X,Rez,N,kr) ! A*X=Rez
	implicit none
	real*8  :: diag(*), ggu(*), ggl(*), X(*), Rez(*)
      integer :: ig(*), jg(*)
      integer :: i, N, j, kr
      real*8  :: s
      
	do i = 1, N
		Rez(i) = diag(i)*X(i)
	end do

	do i = 1, N
	    s = 0.D0
		do j = ig(i), ig(i + 1) - 1
			s          = s + ggl(j)*X(jg(j))
			Rez(jg(j)) = Rez(jg(j)) + ggu(j)*X(i)
		end do
	    Rez(i) = Rez(i) + s
	end do
	end
	

	real*8 function scalar3(X,Y,N)
	implicit real*8 (a-h, o-z)	
	complex*16 :: X(*), Y(*)
	real*8     :: s, a, b
	integer    :: i, N
      s = 0.D0
	a = 0.D0
	b = 0.D0
	do i = 1, N
		a = REAL(X(i))
		b = REAL(Y(i))
		s = s + a*b
		a = AIMAG(X(i))
		b = AIMAG(Y(i))
		s = s + a*b
	enddo
	scalar3 = s
	end


	subroutine bcgstab(ig,jg,ggl,ggu,di,right,x,N,kr,eps)
	implicit none
	real*8  :: eps, enorma, fnorma, eps3, f_norma, a_norma
	real*8  :: di(*), ggu(*), ggl(*)
	real*8  :: right(*),x(*)
      integer :: ig(*), jg(*)
      integer :: N, kr, i, ierr, j, maxiter
      real*8  :: scalar2, a1, a2, alfa, beta, a3, w
      
	real*8, allocatable :: r(:)
	real*8, allocatable :: rt(:)
	real*8, allocatable :: p(:)
	real*8, allocatable :: s(:)
	real*8, allocatable :: q(:)
	real*8, allocatable :: h(:)
	real*8, allocatable :: x0(:)
      
	eps3    = 1e-15
	maxiter = 10000000
	allocate(r(N),stat = ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(rt(N),stat = ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(p(N),stat = ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(s(N),stat = ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(q(N),stat = ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(h(N),stat = ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(x0(N),stat = ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'
	! нулевая итерация
	rt = 0.D0
	r  = 0.D0
	s  = 0.D0
	p  = 0.D0
	q  = 0.D0
	h  = 0.D0
	x0 = 0.D0

!	call umnog(di,ggu,ggl,ig,jg,x0,q,N,kr) ! A*x0=q
	do i = 1, N
	   r(i)  = right(i) - q(i)
	   rt(i) = r(i) ! 1d0
	   p(i)  = r(i)
	end do

	f_norma = scalar2(right,right,N)
	f_norma = sqrt(f_norma)
	! основной цикл	
	do i = 1, maxiter
		a1 = scalar2(r,rt,N)
	    call umnog(di,ggu,ggl,ig,jg,p,q,N,kr) ! A*p=>q
		a2   = scalar2(q,rt,N)    
	    alfa = a1/a2
	    do j = 1, N
	       s(j) = r(j) - alfa*q(j)
	    end do
	    call umnog(di,ggu,ggl,ig,jg,s,h,N,kr) ! A*s=>h
		a2 = scalar2(s,h,N)
		a3 = scalar2(h,h,N)
		w  = a2/a3	    
	    do j = 1, N
	       x0(j) = x(j)
	       x(j)  = x(j) + alfa*p(j) + w*s(j)
	       r(j)  = s(j) - w*h(j)
	    end do
		a2   = scalar2(r,rt,N)
	    beta = a2/a1*alfa/w
	    
	    do j = 1, N
		   p(j) = r(j) + beta*(p(j) - w*q(j))
	    end do
	    a_norma = scalar2(r,r,N)
	    a_norma = sqrt(a_norma)/f_norma
		
	    if (a_norma.lt.eps) then
	       goto 10
	    end if
	    if (abs(beta).le.eps3) then
	        print *, i, "NO SOLUTION"
	       goto 11
	    end if

	print * , i, a_norma
	end do

10    print * , i, a_norma
	open  (11, file = 'vector')
	do i = 1,N
		write(11,*), i, x(i)
	end do 
	close(11)
	deallocate(r)
	deallocate(rt)
	deallocate(p)
	deallocate(s)
	deallocate(q)
	deallocate(h)
	deallocate(x0)	
11	end


	real*8 function scalar2(a1,a2,N)
	implicit none
	real*8  :: a1(N), a2(N)
      real*8  :: s
      integer :: i, N
      
	s = 0d0
	do i = 1, N
		s = s + a1(i)*a2(i)
	end do
	scalar2 = s
	end