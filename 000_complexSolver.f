	subroutine umnogComplexReal(diag,ggu,ggl,ig,jg,X,Rez,N,kr) ! A*X=Rez
	implicit complex*16 (a-h, o-z)

	real*8     :: diag(N),ggu(kr),ggl(kr)
      complex*16 :: X(N),Rez(N)
      integer    :: ig(N+1),jg(kr)
	Rez=(0.D0,0.D0)  
	do i=1, N
		Rez(i)=diag(i)*X(i)
	end do

	do i=1, N
	    s=(0.D0,0.D0)
		do j=ig(i),ig(i+1)-1
			s=s+ggl(j)*X(jg(j))
			Rez(jg(j))=Rez(jg(j))+ggu(j)*X(i)
		end do
	    Rez(i)=Rez(i)+s
	end do
	end
      
      subroutine umnogComplex(diag,ggu,ggl,ig,jg,X,Rez,N,kr) ! A*X=Rez
	implicit complex*16 (a-h, o-z)

	dimension diag(N),ggu(kr),ggl(kr), ig(N+1),jg(kr),X(N),Rez(N)
	Rez=(0.D0,0.D0)  
	do i=1, N
		Rez(i)=diag(i)*X(i)
	end do

	do i=1, N
	    s=(0.D0,0.D0)
		do j=ig(i),ig(i+1)-1
			s=s+ggl(j)*X(jg(j))
			Rez(jg(j))=Rez(jg(j))+ggu(j)*X(i)
		end do
	    Rez(i)=Rez(i)+s
	end do
	end

	subroutine CoCG(ig,jg,ggl,ggu,di,right,x,N,kr,eps)
	implicit complex*16  (a-h, o-z)
	real*8 eps, enorma, fnorma
	integer*2 t1,t2,t3,t4,t5,t6,t7,t8
	dimension ig(N+1),jg(kr),ggl(kr),ggu(kr),di(N),right(N),x(N)

	complex*16, allocatable :: r(:)
	complex*16, allocatable :: p(:)
	complex*16, allocatable :: w(:)
	complex*16, allocatable :: u(:)
	complex*16, allocatable :: q(:)
						
	allocate(r(N), stat=ierr1)
	if(ierr1/=0)stop 'Allocation error in r'

	allocate(w(N), stat=ierr1)
	if(ierr1/=0)stop 'Allocation error in w'
	
	allocate(p(N), stat=ierr1)
	if(ierr1/=0)stop 'Allocation error in p'		

	allocate(u(N), stat=ierr1)
	if(ierr1/=0)stop 'Allocation error in x'

	allocate(q(N), stat=ierr1)
	if(ierr1/=0)stop 'Allocation error in q'

	w=(0.D0,0.D0)
	r=(0.D0,0.D0)
	p=(0.D0,0.D0)
	u=(0.D0,0.D0)
	x=(0.D0,0.D0)   !x=x0
	q=(0.D0,0.D0)
	beta=(0.D0,0.D0)
	alfa=(0.D0,0.D0)
	alfa1=(0.D0,0.D0)
	alfa2=(0.D0,0.D0)
	maxiter=10000000

	call umnogComplex(di,ggu,ggl,ig,jg,x,q,N,kr) ! A*X=Rez
	do i=1, N
		r(i)=right(i)-q(i)
	    w(i)=r(i)/di(i)
	enddo
		alfa1=CoScalar(r,w,N)
	    fnorma=Scalar(right,right,N)
	do ix=1, maxiter
		do i=1, N
			p(i)=w(i)+beta*p(i)
	    enddo
	    call umnogComplex(di,ggu,ggl,ig,jg,p,u,N,kr)
	    alfa2=CoScalar(u,p,N)
	    alfa=alfa1/alfa2
		do i=1, N
			x(i)=x(i)+alfa*p(i)
			r(i)=r(i)-alfa*u(i)
		enddo
		enorma=Scalar(r,r,N)
	    enorma=sqrt(enorma/fnorma)
	    if(enorma.lt.eps)then
	       goto 1
	    else
	       do i=1, N
	          w(i)=r(i)/di(i)
	       enddo
		   alfa2=alfa1
	       alfa1=CoScalar(r,w,N)
	       beta=alfa1/alfa2
	       print *, ix, enorma
	    endif
	enddo

1	print *, enorma
	w=0d0
	call umnogComplex(di,ggu,ggl,ig,jg,x,w,N,kr)  ! q<=Ax0
	do i=1, N
	   r(i)=right(i)-w(i)
	enddo
	enorma=Scalar(r,r,N)
	print *, ix,enorma
	write(1,*),ix,enorma

	enorma=sqrt(enorma/fnorma)		
	print *, enorma

	deallocate(p)
	deallocate(w)
	deallocate(r)
	deallocate(u)
	deallocate(q)
	end

	function CoScalar(X,Y,N)
	implicit complex*16 (a-h, o-z)
	dimension X(N), Y(N)

	s=(0.D0,0.D0)

	do i=1, N
		s=s+X(i)*Y(i)
	enddo
	CoScalar=s
	end


	subroutine CoBcgstab(ig,jg,ggl,ggu,di,right,x,N,kr,eps)
	implicit complex*16  (a-h, o-z)
	real*8 eps, enorma, f_norma,eps3,a_norma 
	dimension di(*), ggu(*), ggl(*), ig(*), jg(*)
	dimension right(*),x(*)

	complex*16, allocatable :: r(:)
	complex*16, allocatable :: rt(:)
	complex*16, allocatable :: p(:)
	complex*16, allocatable :: s(:)
	complex*16, allocatable :: q(:)
	complex*16, allocatable :: h(:)
	complex*16, allocatable :: x0(:)
	eps3=1e-15
	maxiter=10000000
	allocate(r(N),stat=ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(rt(N),stat=ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(p(N),stat=ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(s(N),stat=ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(q(N),stat=ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(h(N),stat=ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'

	allocate(x0(N),stat=ierr)
	if(ierr.ne.0) stop 'Allocation error in bcg'
	! нулевая итерация
	rt=(0.d0,0.d0)
	r=(0.d0,0.d0)
	s=(0.d0,0.d0)
	p=(0.d0,0.d0)
	q=(0.d0,0.d0)
	h=(0.d0,0.d0)
	x0=(0.d0,0.d0)   !!!!ОПАСНОЕ МЕСТО
!	x0(1)=(1.d0,0.d0)
      beta=(0.d0,0.d0)
      a1=(0.d0,0.d0)
      a2=(0.d0,0.d0)
      alfa=(0.d0,0.d0)
      w=(0.d0,0.d0)
!	call umnogComplex(di,ggu,ggl,ig,jg,x0,q,N,kr) ! A*x0=q
	do i=1, N
	   r(i)=right(i)-q(i)
	   rt(i)=(1.d0,0.d0) !r(i)! 1d0
	   p(i)=r(i)
	end do 
	f_norma=Scalar(right,right,N)
	f_norma=sqrt(f_norma)
	! основной цикл	
	do i=1,maxiter
		a1=CoScalar2(r,rt,N)
	  call umnogComplex(di,ggu,ggl,ig,jg,p,q,N,kr) ! A*p=>q
		a2=CoScalar2(q,rt,N)    
	    alfa=a1/a2
	    do j=1, N
	       s(j)=r(j)-alfa*q(j)
	    end do
	    call umnogComplex(di,ggu,ggl,ig,jg,s,h,N,kr) ! A*s=>h
		a2=CoScalar2(s,h,N)
		a3=CoScalar2(h,h,N)
		w=a2/a3	         
	    do j=1, N
	       x0(j)=x(j)
	       x(j)=x(j)+alfa*p(j)+w*s(j)
	       r(j)=s(j)-w*h(j)
	    end do
		a2=CoScalar2(r,rt,N)
	    beta=a2/a1*alfa/w
	    
	    do j=1, N
		   p(j)=r(j)+beta*(p(j)-w*q(j))
	    end do
	    a_norma=Scalar(r,r,N)
	    a_norma=sqrt(a_norma)/f_norma
		
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
	open  (11,file='vector')
	do i=1,N
		write(11,*), i, x(i)
	end do 
	close(11)
      do i=1, N
         print *, x(i)
      enddo
      read(*,*), ch
	deallocate(r)
	deallocate(rt)
	deallocate(p)
	deallocate(s)
	deallocate(q)
	deallocate(h)
	deallocate(x0)	
11	end

	function CoScalar2(X,Y,N)
	implicit complex*16 (a-h, o-z)
	dimension X(N), Y(N)
      real*8 :: a,b, a1,a2,b1,b2
	s=(0.D0,0.D0)
	do i=1, N
        a1=Real(X(i))
        b1=AIMAG(X(i))

        a2=Real(Y(i))
        b2=AIMAG(Y(i))

        a=a1*a2-b1*b2
        b=a2*b1+a1*b2
        s=s+cmplx(a,b)
	enddo
	CoScalar2=s
	end