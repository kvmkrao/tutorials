	implicit real(a-h,o-z)
	common i,j,k,l,m,n,n1,h,h1,ajacob,nel,npts,barlt,xi,nn
	common alpha1,alpha2,dampcof,maxgqpts,ngqpts,igaus1,igaus2,maxgaus
	common ibtype1,ibtype2,bdata1(5),bdata2(5)
	common f1,f2,f3,eps,alpha,maxits,sorelax
	dimension ea(1000,5),f(1000,5),weight(10),gqpt(10)
	dimension zeta(10),sh(10,10),sh1(10,10),stiffmt(1000,1000)
	dimension x(1000),elementk(10,10),forcemt(1000)
	dimension finalk(1000,1000),finalf(1000),soln(1000)
	character*30 key
	common input,twicese
	input =5
c
	read(input,*) key
	write(6,*)key
102     read(input,*) key
!        write(6,*)key
        if (key.eq.'done') then
          goto 101
        else if (key.eq.'order') then
          read(input,*)npts
	else if(key.eq.'numel') then
	  read(input,*) nel
	else if(key.eq.'barlength') then
          read(input,*) barlt
	else if(key.eq.'initialPoint') then
          read(input,*) xi
	else if(key.eq.'EAconst1') then
	  read(input,*) alpha1
	else if(key.eq.'EAconst2') then
	 read(input,*) alpha2
	else if(key.eq.'Fconst1') then
	 read(input,*) f1
	else if(key.eq.'Fconst2') then
	 read(input,*) f2
	else if(key.eq.'Fconst3') then
	 read(input,*) f3
	else if(key.eq.'dampCoeff') then
	 read(input,*) dampcof
	else if (key.eq.'RelaxFact') then
	 read(input,*) sorelax
	else if (key.eq.'epsilon') then
	 read(input,*) eps
	else if(key.eq.'maxiter') then
	 read(input,*) maxits
	else if(key.eq.'bctype1') then
	 read(input,*) ibtype1
	else if(key.eq.'bcdata01') then
	 read(input,*) bdata1(1)
	else if(key.eq.'bcdata02') then
	 read(input,*) bdata1(2)
	else if(key.eq.'bctype2') then
	 read(input,*) ibtype2
	else if(key.eq.'bcdatal1') then
	 read(input,*) bdata2(1)
	else if(key.eq.'bcdatal2') then
	 read(input,*) bdata2(2)
	 end if
	go to 102
101	continue
c
	igaus1 = npts
!	calculating (Na1, Na2) term's ploynomial order
	if(npts.gt.1) igaus2 = (npts*2)-2
!	evaluating the (Na,f) term's ploynomial order
	if(f2.ne.0) igaus1 = npts+1
	if(f3.ne.0) igaus1 = npts+2
!	calculating Max order of ploynomial for Max number of 
!	gauss quadrature points
	if(igaus1.ge.igaus2) then
	maxgaus = igaus1
	else 
	maxgaus = igaus2
	end if
	write(6,*)" max order of the polymonial",maxgaus
c
	ngqpts = (maxgaus+1)/2
!	roundoff the max gauss quadrature points for
!	maximum  odd order of the polynomial( e.g 2.5 = 3)
	if(mod(maxgaus,2) .eq.0) ngqpts = 1 + (maxgaus+1)/2
c
	write(6,*) " number of G-Q points",ngqpts
c
	h1   = barlt/nel
	nn     =  nel * npts + 1  	!total number of nodes
c
	write(6,*) "total number of nodes ",nn
c
	call shapefun(npts,nel,sh,sh1,barlt,ajacob)
	open(unit =10, file ='ShapeCoeff.dat', status='unknown')
c
	do i=1,npts+1
	write(10,*)"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
		write(10,*) "shape function",i, "is"
		write(10,*) (sh(i,k),'x**',npts-k+1,k=1,npts+1)
	 	write(10,*)"derivative of shape function",i,"is"
		write(10,*) (sh1(i,j),'x**',npts-j,j=1,npts)
	end do
	close(10)
c
	call mesh(nn,xi,barlt,x)
	open(unit=11,file='mesh.dat',status='unknown')
	write(11,*)"number of nodes",nn
	write(11,*)"number of elements",nel
	write(11,*) "node no location"
	write(11,111)(i,x(i),i=1,nn)
111	format(i4,4x,1f9.5)
	close(11)
c
!	converting  EA as function of X to zeta
	call matdata(npts,nel,alpha1,alpha2,h1,x,ea)
c
	call forcedata(f1,f2,f3,h1,nel,npts,x,f)
!	calculating Gauss quadrature points and weights corresponding to 
!	gauss quadrature points  
	call gauss(ngqpts,npts,weight,gqpt,maxgqpts)
c
!	 Global stiffness matrix and load vector formation
	call elemstiff(ea,sh,sh1,f,ajacob,weight,gqpt,dampcof,npts,nn,
     $				nel,ngqpts,h1,stiffmt,forcemt)
c
	write(6,*)"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
	write(6,*) (forcemt(i),i=1,npts+1)
	do i=1,nn
		write(6,131) (stiffmt(i,j),j=1,nn)
	end do
	write(6,*)"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
c
!	Implementing boundary conditions
!	modified Global stiffness matrix and load vector
	call boundary(nn,stiffmt,forcemt,ibtype1,ibtype2,bdata1,
     $			bdata2,finalk,finalf)
c
	open(unit=12,file='forceVector.dat',status='unknown')
	write(12,*)" force vector is"
	open(unit=13,file='Stiffmatrix.dat',status='unknown')
	write(13,*)"stiffness matrix of order",nn,'x',nn
	write(13,*) "stiffness matrix is"
c
	do i=1,nn
		write(13,131)(finalk(i,j),j=1,nn)
		write(12,121) finalf(i)
	end do
121	format(1e12.5)
131	format(4(2x,1f9.5))
	close(12)
	close(13)
c
	call solver(finalk,finalf,soln,nn,maxits,eps,sorelax)
c
	write(6,*)"solution"
	do i=1,nn
	write(6,*) soln(i)
	end do

	twicese =0.0d0
	do i=1,nn
	twicese = twicese +finalf(i) * soln(i)
	end do
	write(6,*) "twice strain energy",twicese
c
	open(unit =14,file='GSpoints.dat',status='unknown')
	write(14,*) "number of  gauss quadrature pts",ngqpts
	write(14,*) "Gauss Quadrature points"
	write(14,*)" pointNo. location, weight " 
	 do i=1, ngqpts
	write(14,*) i,gqpt(i),weight(i)
	end do
	close(14)
c
	write(6,*) "max gauss quadrature pts",ngqpts
c
100	format(a30)
	write(6,*)" jacobian ",ajacob
c
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine boundary(nn,stiffmt,forcemt,ibtype1,ibtype2,bdata1,
     $				bdata2,condK,condF)

	dimension condK(1000,1000), condF(1000)
	dimension stiffmt(1000,1000), forcemt(1000)
	dimension bdata1(5),bdata2(5)

!******************** bytype1 is specified at X=0 ********************
!	ibctype1 =1 for displacement boundary condition at x=0
!			        Uo  o-----------
!	 Uo = bdata1(1)
	if (ibtype1.eq.1) then
	condF(1) = bdata1(1)

	do i=2,nn
		do j=1,nn
		condK(i,j) = stiffmt(i,j)
		end do
	condF(i) = forcemt(i)
	end do
	do i=2,nn
		condF(i) = forcemt(i) -stiffmt(i,1)*bdata1(1)
		condK(i,1) = 0.0d0
		condK(1,j) = 0.0d0
	end do
	condK(1,1) = 1.0d0
!	ibctype1 =2 for force boundary condition at X=0
!				P1 <------ ||---------
!	P1 = bdata1(1)
	else if(ibtype1.eq.2) then 
		condF(1) = forcemt(1) - bdata1(1)
	do i=1,nn
		do j=1,nn
		condK(i,j) = stiffmt(i,j)
		end do
	if (i.ge.2)  condF(i) = forcemt(i)
	end do
!	ibctype1 =3 for spring load at X =0   
!				           ||----^^^^---
!					         k1, delta1		
!	k1       = bdata1(1)
!	delata1  = bdata1(2) 					
	else if (ibtype1.eq.3) then
		condK(1,1) = stiffmt(1,1) + bdata1(1)
		condF(1)   = forcemt(1) + bdata1(1)*bdata1(2)
		do i=1,nn
			do j=1,nn
			if(i.ne.1.and.j.ne.1) condK(i,j) = stiffmt(i,j)
			end do
		if (i.ge.2)  condF(i) = forcemt(i)
		end do
! ************	 Implementation of Bctype1 is completed **************	
! ************    Implementing  bytype2 at X=L          **************
!	 ibctype2 =1 for displacement boundary condition at x=L
!				  -----------o U_L
!	U_L = bdtata2(1)
	else if (ibtype2.eq.1) then   
		condF(nn) = bdata2(1)
		do i=1,nn-1
		condF(i) = forcemt(i) - stiffmt(i,nn) * bdata2(1)
		condK(i,nn) = 0.0d0
		condK(nn,i) = 0.0d0
		end do
!       ibctype2 =2 for force boundary condition at X=0
!                                -----------|| -----> P2
!	P2 = bdata2(1)
	else if (ibtype2.eq.2) then
		condF(nn) = forcemt(nn) + bdata2(1)
		do i=1,nn
			do j=1,nn
			condK(i,j) = stiffmt(i,j)
			end do
		if (i.ne.1) condF(i) = forcemt(i)
		end do
!       ibctype2 =3 for spring load at X =0   
!                                ----^^^^---||
!                                   k2, delta2 
!	k2       = bdata2(1)
!	delata2  = bdata2(2)
	else if (ibtype2.eq.3) then
		condK(nn,nn) = stiffmt(nn,nn) + bdata2(1)
		condF(nn) = forcemt(nn) - bdata2(1)*bdata2(2)
		do i=1,nn
			do j=1,nn
			if(i.ne.nn.and.j.ne.nn) condK(i,j) = stiffmt(i,j)
			enddo
		if(i.ne.nn) condF(i) = forcemt(i)
		end do
	end if
! ************  Implementing  bytype2 at X=L is completed  **************
c
	return
	end
c
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine elemstiff(ea,sh,sh1,f,ajacob,weight,gqpt,dampcof,
     $			npts,nn,nel,ngqpts,h1,globalk,globalf)
	implicit real(a-h,o-z)
	dimension sh(10,10),sh1(10,10),f(1000,5),weight(10)
	dimension eleK(10,10),eleF(10),gqpt(10),ea(1000,5)
	dimension globalk(1000,1000), globalf(1000) 
	real*8 temp(1000,10),temp1(1000,10),ym(1000,10),ftemp(1000,10)
	integer ij,jk,itemp
c
!*******************	INTIALISATION	**********************************
	do  i=1,npts+1
		do j=1, ngqpts
			temp(i,j) = 0.0d0
			if (i.ge.2) then
				temp1(i-1,j) = 0.0
			end if
			ym(i,j)      = 1.0d0
			ftemp(i,j)   = 1.0d0
		end do
	end do
c
	do i = 1,nn
		do j=1,nn
			globalk(i,j) = 0.0d0
		end do
	globalf(i) = 0.0d0
	end do
!*************************************************************************
!	local shape functions and its derivatives at Gauss 
!	Quadrature  points"
	do k=1,ngqpts
		do i=1,npts+1
			do j=1,npts+1
				temp(i,k) = temp(i,k) + sh(i,j)*  
     $				(gqpt(k)**(npts+1-j))
			end do
			do j=1,npts
				temp1(i,k) = temp1(i,k) + (sh1(i,j)*
     $				gqpt(k)**(npts-j))
			end do
		end do
	end do
!	EA and f are as a function of X. so, these are not constant at each node
!	EA  and f is evaluating at zeta for each node in an element 
	l = 1            	!starting node number of an element
	m = npts+1       	!ending node number of an element
	do itemp = 1,nel   			!number of elements
c
	do k=1,ngqpts
		do i=l,m
			ym(i,k) = ea(i,1) + ea(i,2) * gqpt(k)
			ftemp(i,k) = f(i,1) + f(i,2)*gqpt(k)+ f(i,3) * 
     $			gqpt(k)**2
!	global to local transformation for EA and f in an element   
			do j=1,npts+1
				ym(j,k)    = ym(i,k)
				ftemp(j,k) = ftemp(i,k) 
			end do
		end do
	end do
!******* LOCAL STIFFNESS MATRIX and LOAD VECTOR for an element ***********
	do k=1,ngqpts
		do i=1,npts+1
			do j=1,npts+1
			eleK(i,j) = ym(i,k)*temp1(i,k)*temp1(j,k)
     $					*(2.0/h1) *weight(k) + 
     $					dampcof*temp(i,k)*
     $					temp(j,k)*ajacob*weight(k)
			end do
			eleF(i) = ftemp(i,k)*temp(i,k)*weight(k)*ajacob
		end do
	end do
!*************************************************************************
!*****************************  ASSEMBLY **********************************
!	Adding the local element stiffness matrix and load vector to 
!	global stiffness matrix and load vector for an element
	do i=1, npts+1
		j=itemp 
		ij = (j-1)*npts + i
		globalf(ij) = globalf(ij) + eleF(i)
			do k = 1,npts+1
				jk = (j-1)*npts + k
				globalk(ij,jk) = globalk(ij,jk)+ eleK(i,k) 
			end do
	end do
!!**************************************************************************
	l = l+ npts 
	m = m+ npts 
	end do
	return
	end 
c
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine gauss(ngqpts,npts,wght,qpt,maxn)
	dimension wght(10),qpt(10)
!	ONE POINT GAUSS QUADRATURE RULE
	if(ngqpts.eq.1) then
                qpt(ngqpts) =  0.0d0
                wght(ngqpts) = 2.0d0
!       TWO POINT GAUSS QUADRATURE RULE
        else if (ngqpts.eq.2) then
                qpt(ngqpts-1) =   -sqrt(1.0/3)
                wght(ngqpts-1)  =   1.0d0
                qpt(ngqpts)   =    -qpt(ngqpts-1)
                wght(ngqpts)   =    wght(ngqpts-1)
!       THREE POINT GAUSS QUADRATURE RULE
        else if (ngqpts.eq.3) then
                qpt(ngqpts -2) =   - sqrt(3.0d0/5)
                wght(ngqpts-2)  =   5.0d0/9
                qpt(ngqpts-1)  =   0.0d0
                wght(ngqpts-1)  =  8.0d0/9
                qpt(ngqpts)   =   sqrt(3.0d0/5)
                wght(ngqpts)    =  wght(ngqpts-2)
!       FOUR POINT GAUSS QUADRATURE RULE
        else if(ngqpts.eq.4) then
                qpt(ngqpts-3)  =    sqrt((3.0 + (2.0*sqrt(6.0/5)))/7.0)
                wght(ngqpts-3)  =    (18.0 - sqrt(30.0))/36
                qpt(ngqpts-2)  =    sqrt((3.0 - (2.0*sqrt(6.0/5)))/7)
                wght(ngqpts-2)  =    (18.0 + sqrt (30.0))/36
                qpt(ngqpts-1)  =   -sqrt((3.0 - (2.0*sqrt(6.0/5)))/7)
                wght(ngqpts-1)  =    wght(ngqpts-2)
                qpt(ngqpts)    =   -sqrt((3.0 + (2.0*sqrt(6.0/5)))/7)
                wght(ngqpts)    =    wght(ngqpts-3)
!       FIVE POINT GAUSS QUADRATURE RULE
        else if(ngqpts.eq.5) then
                qpt(ngqpts-4)  =   (sqrt(5.0+(2*sqrt(10.0/7))))/3
                wght(ngqpts-4)  =   (322.0 - (13.0 * sqrt(70.0)))/900
                qpt(ngqpts-3)  =   (sqrt(5.0-(2*sqrt(10.0/7))))/3
                wght(ngqpts-3)  =   (322.0 + (13.0 * sqrt(70.0)))/900
                qpt(ngqpts-2)  =   0.0d0
                wght(ngqpts-2)  =   128.0/225
                qpt(ngqpts-1)  =   -(sqrt(5.0-(2*sqrt(10.0/7))))/3
                wght(ngqpts-1)  =    wght(ngqpts-3)
                qpt(ngqpts)    =   -(sqrt(5.0+(2*sqrt(10.0/7))))/3
                wght(ngqpts)   =     wght(ngqpts-4)
        end if
c
	return
	end 
c
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine forcedata(f1,f2,f3,h1,nel,npts,x,force)
c
	implicit real(a-h,o-z)
	dimension force(1000,5),x(1000)
c
!	 f     = f1 + f2 *x + f3 * x**2
!        f(i)  = f(i,1) + f(i,2) *zeta + f(i,3) * zeta**2
        do i=1,npts*nel+1
        force(i,1) = f1 + (f2*h1/2.0d0)+(f3*h1*h1/4.0d0)+
     $	x(i)*(f2+f3*x(i)+f3*h1)
        force(i,2) = (f2*h1/2.0d0)+(f3*h1*x(i))+(f3*h1*h1/2.0d0)
	force(i,3) = f3*h1**2/4.0d0 
        end do
c
	return
	end

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	subroutine matdata (npts,nel,alpha1,alpha2,h1,x,young)
c
	implicit real (a-h,o-z)
	dimension young(1000,5) ,x(1000)
!	EA     = alpha1 + alpha2 * x
!	EA(i)  = young(i,1) + young(i,2) * zeta
	do i=1,npts*nel+1
	young(i,1) = (alpha1 +(alpha2*h1)/2.0d0) +alpha2 *x(i)
	young(i,2) = alpha2*h1/2.0d0
	end do
c
	return
	end
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
       subroutine shapefun(npts,nel,b,b1,barlt,ajacob)
        implicit real(a-h,o-z)
c
	dimension zeta(10),b(10,10),b1(10,10)
	common i,j,k,l,m,n,n1,h
	 real*8 const(10),c0(10),c1(10)
	 real*8  c3(10),c4(10),c5(10),c6(10),c2(10)
c
        do i=1,npts+1
                c0(i) = 0.0d0
                c1(i) = 0.0d0
                c2(i) = 0.0d0
                c3(i) = 0.0d0
                c4(i) = 0.0d0
                c5(i) = 0.0d0
                c6(i) = 0.0d0
                const(i) = 1.0d0
		do j=1,npts+1
		b(i,j)  = 1.0d0
		b1(i,j) = 1.0d0
		end do
        end do
c
!       INITIAL AND FINAL POINT
        do i=1,npts+1
        zeta(i) = 1.0 - 2.0*(i-1)/npts
        end do
        write(6,*)"number of grid points in a element", npts+1
        do i=1,npts+1
        write(6,*)"location of point ",i
        write(6,*) -zeta(i)
        end do
c       
	ajacob = barlt/(nel*2.0d0)
c
        do i=1,npts+1
        c0(i) = c0(i)+1.0d0

          do j=1,npts+1
                if(j.ne.i) then
                 c1(i) = c1(i)+zeta(j)
                   do k=j+1,npts+1
                      if (k.ne.i) then
                        c2(i) = c2(i) + zeta(j)*zeta(k)
                       do l=k+1,npts+1
                          if (l.ne.i) then
                            c3(i) = c3(i) + zeta(j)*zeta(k)*zeta(l)
                             do m=l+1,npts+1
                              if (m.ne.i) then
                               c4(i) = c4(i) + zeta(j)*zeta(k)*zeta(l)*
     $                                  zeta(m)
                                do n = m+1,npts+1
                                  if (n.ne.i) then
                                     c5(i) = c5(i) + zeta(j)*zeta(k)
     $                                       *zeta(l)* zeta(m)*zeta(n)
                                      do n1 =n+1,npts+1
                                        if(n1.ne.i) then
                                        c6(i) = c6(i) + zeta(j)*zeta(k)*
     $                                  zeta(l)*zeta(m)*zeta(n)*
     $                                  zeta(n1)
                                        end if
                                      end do
                                  end if
                                end do
                              end if
                             end do
                          end if
                       end do
                      end if
                   end do
                end if
           end do

       do j=1,npts+1
          if (j.ne.i) then
            const(i) = const(i)*(-zeta(i) +zeta(j))
          end if
        end do
c
        b(i,1) = c0(i)/const(i)
        b(i,2) = c1(i)/const(i)
        b(i,3) = c2(i)/const(i)
        b(i,4) = c3(i)/const(i)
        b(i,5) = c4(i)/const(i)
        b(i,6) = c5(i)/const(i)
        b(i,7) = c6(i)/const(i)
c
        write(6,*)"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
        write(6,*) "shape function",i, "is"
        write(6,*)( b(i,k),'x**',npts-k+1,k=1,npts+1)
c
        do k=1,npts
        b1(i,k) = b(i,k) * (npts-k+1)
        end do
        write(6,*)"derivative of shape function",i,"is"
        write(6,*) (b1(i,k),'x**',npts-k,k=1,npts)
        end do
c        
	return
        end

