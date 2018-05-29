MODULE cheb

  USE prms

  IMPLICIT none

  real, allocatable, dimension(:,:,:)  :: Gr,Gt
  real, allocatable, dimension(:,:,:)  :: Dccr,Dddr,Dcdr
  real, allocatable, dimension(:,:,:)  :: Dcct,Dddt,Dcdt
  real, allocatable, dimension(:,:)    :: Tccr,Tddr,Tcdr,Thddr,Thccr,Thdcr
  real, allocatable, dimension(:,:)    :: Tcct,Tddt,Tcdt,Thddt,Thcct,Thdct
  real, allocatable, dimension(:,:)    :: Tfdr,Tfdt
  real, allocatable, dimension(:,:,:)  :: Hdr
  real, allocatable, dimension(:)      :: Xcr,Xdr,Xct,Xdt,Xfr,Xft
  real, allocatable, dimension(:)      :: Rc,Tc,Rd,Td,Rf,Tf,Rdd,Tdd

  real,allocatable,dimension(:,:)     :: XXdr,XX2dr,XXdt,XX2dt,XX3dr,XX3dt,XX4dr,XX4dt

  integer, parameter :: NDmax = 6

CONTAINS

  SUBROUTINE initcheb

    call getcol
    call makeTs
    call makeDs
    call makeGs
    call makeHs
    call makeRmat

  END SUBROUTINE initcheb

  SUBROUTINE makeRmat

    integer   :: i,j

    allocate(XXdr(0:Ndr,0:Ndr))
    allocate(XX2dr(0:Ndr,0:Ndr))
    allocate(XXdt(0:Ndr,0:Ndr))
    allocate(XX2dt(0:Ndr,0:Ndr))
    allocate(XX3dr(0:Ndr,0:Ndr))
    allocate(XX3dt(0:Ndr,0:Ndr))
    allocate(XX4dr(0:Ndr,0:Ndr))
    allocate(XX4dt(0:Ndr,0:Ndr))

    Do i = 0,Ndr
      Do j = 0,Ndr
        If (i == j) then
          XXdr(i,j) = Rdd(i)
          XXdt(i,j) = Tdd(i)
        else
          XXdr(i,j) = 0
          XXdt(i,j) = 0
        End If
      End Do
    End Do

    XX2dr(:,:) = Matmul(XXdr(:,:),XXdr(:,:))
    XX2dt(:,:) = Matmul(XXdt(:,:),XXdt(:,:))

    XX3dr(:,:) = Matmul(XX2dr(:,:),XXdr(:,:))
    XX3dt(:,:) = Matmul(XX2dt(:,:),XXdt(:,:))

    XX4dr(:,:) = Matmul(XX2dr(:,:),XX2dr(:,:))
    XX4dt(:,:) = Matmul(XX2dt(:,:),XX2dt(:,:))

    END SUBROUTINE makeRmat


  SUBROUTINE getcol

    integer  :: i,j

    allocate(Xcr(0:Ncr),Xct(0:Nct))
    allocate(Xdr(0:Ndr),Xdt(0:Ndt))
    allocate(Xfr(0:Nfr),Xft(0:Nft))
    allocate(Rc(0:Ncr),Tc(0:Nct))
    allocate(Rd(0:Ndr),Td(0:Ndt))
    allocate(Rf(0:Nfr),Tf(0:Nft))
    allocate(Rdd(0:Ndr),Tdd(0:Ndr))


    call makeX(Ncr,Xcr)
    call makeX(Nct,Xct)
    call makeX(Ndr,Xdr)
    call makeX(Ndt,Xdt)
    call makeX(Nfr,Xfr)
    call makeX(Nft,Xft)

    Rc = 0.5*Rmax*(Xcr+1.)
    Tc = 0.5*Tmax*(Xct+1.)+Rmax

    Rd = 0.5*Rmax*(Xdr+1.)
    Td = 0.5*Tmax*(Xdt+1.)+Rmax

    Rf = 0.5*Rmax*(Xfr+1.)
    Tf = 0.5*Tmax*(Xft+1.)+Rmax

    !inverse order
    Do i=0,Ndr
    Rdd(i) = Rd(Ndr-i)
    Tdd(i) = Td(Ndr-i)
    End Do

    END SUBROUTINE getcol

  SUBROUTINE makeX(N,X)
    integer                   :: N
    real, dimension(0:N)  :: X

    integer :: i

    do i = 0,N
       X(i) = COS(REAL(i)*Pi/N)
    end do

    END SUBROUTINE makeX

  SUBROUTINE makeTs

    allocate(Tccr(0:Ncr,0:Ncr))
    allocate(Tddr(0:Ndr,0:Ndr))
    allocate(Tcdr(0:Ncr,0:Ndr))

    allocate(Tcct(0:Nct,0:Nct))
    allocate(Tddt(0:Ndt,0:Ndt))
    allocate(Tcdt(0:Nct,0:Ndt))

    allocate(Thccr(0:Ncr,0:Ncr))
    allocate(Thddr(0:Ndr,0:Ndr))
    allocate(Thdcr(0:Ndr,0:Ncr))

    allocate(Thcct(0:Nct,0:Nct))
    allocate(Thddt(0:Ndt,0:Ndt))
    allocate(Thdct(0:Ndt,0:Nct))

    allocate(Tfdr(0:Nfr,0:Ndr))
    allocate(Tfdt(0:Nft,0:Ndt))

    call makeT(Ncr,Ncr,Tccr)
    call makeT(Ndr,Ndr,Tddr)
    call makeT(Ncr,Ndr,Tcdr)

    call makeTh(Ncr,Ncr,Thccr)
    call makeTh(Ndr,Ndr,Thddr)
    call makeTh(Ndr,Ncr,Thdcr)

    call makeT(Nct,Nct,Tcct)
    call makeT(Ndt,Ndt,Tddt)
    call makeT(Nct,Ndt,Tcdt)

    call makeTh(Nct,Nct,Thcct)
    call makeTh(Ndt,Ndt,Thddt)
    call makeTh(Ndt,Nct,Thdct)

    call makeT(Nft,Ndt,Tfdt)
    call makeT(Nfr,Ndr,Tfdr)

   END SUBROUTINE makeTs

  SUBROUTINE makeGs

    integer :: m
    !deallocate(Gr)
    allocate(Gr(0:Ndr,0:Ndr,0:NDmax))
    allocate(Gt(0:Ndt,0:Ndt,0:NDmax))

    call makeG(Ndr,Gr(0,0,1))
    call makeG(Ndt,Gt(0,0,1))

    call identity(Ndr,Gr(0,0,0))
    call identity(Ndt,Gt(0,0,0))

    Gr(:,:,1) = Gr(:,:,1)*(2./Rmax)
    Gt(:,:,1) = Gt(:,:,1)*(2./Tmax)

    do m = 2,NDmax
       Gr(:,:,m) = MATMUL(Gr(:,:,1),Gr(:,:,m-1))
       Gt(:,:,m) = MATMUL(Gt(:,:,1),Gt(:,:,m-1))
    end do

    END SUBROUTINE makeGs

  SUBROUTINE makeGG(p,q,r,GG)
    integer                                  :: p,q,r
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: GG

    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: GGr,GGt
    integer                                  :: i,j

    do j = 0,Ndt
       GGr(:,j,:,j) = MATMUL(Hdr(:,:,r),Gr(:,:,p))
    end do
    do i = 0,Ndr
       GGt(i,:,i,:) = Gt(:,:,q)
    end do
    call MATMUL2(GGr,GGt,GG)

    END SUBROUTINE makeGG

  SUBROUTINE MATMUL2(A,B,C)
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt)  :: A,B,C

    integer               :: i,j,ii,jj

    C = 0.
    do jj = 0,Ndt
       do ii = 0,Ndr
          do j = 0,Nct
             do i = 0,Ncr
                C(i,j,ii,jj) = SUM(A(i,j,:,:)*B(:,:,ii,jj))
             end do
          end do
       end do
    end do

    END SUBROUTINE MATMUL2

  SUBROUTINE makeDs

    integer :: m

    allocate(Dccr(0:Ncr,0:Ncr,0:NDmax))
    allocate(Dddr(0:Ndr,0:Ndr,0:NDmax))
    allocate(Dcdr(0:Ncr,0:Ndr,0:NDmax))

    allocate(Dcct(0:Nct,0:Nct,0:NDmax))
    allocate(Dddt(0:Ndt,0:Ndt,0:NDmax))
    allocate(Dcdt(0:Nct,0:Ndt,0:NDmax))

    call makeD(Ncr,Ncr,Dccr)
    call makeD(Ndr,Ndr,Dddr)
    call makeD(Ncr,Ndr,Dcdr)

    do m = 1,NDmax
       Dccr(:,:,m) = Dccr(:,:,m)*(2./Rmax)**m
       Dddr(:,:,m) = Dddr(:,:,m)*(2./Rmax)**m
       Dcdr(:,:,m) = Dcdr(:,:,m)*(2./Rmax)**m
    end do

    call makeD(Nct,Nct,Dcct)
    call makeD(Ndt,Ndt,Dddt)
    call makeD(Nct,Ndt,Dcdt)


    do m = 1,NDmax
       Dcct(:,:,m) = Dcct(:,:,m)*(2./Tmax)**m
       Dddt(:,:,m) = Dddt(:,:,m)*(2./Tmax)**m
       Dcdt(:,:,m) = Dcdt(:,:,m)*(2./Tmax)**m
    end do


    END SUBROUTINE makeDs

  SUBROUTINE makeD(Nc,Nd,D)
    integer                     :: Nc,Nd
    real, dimension(0:Nc,0:Nd,0:NDmax)  :: D

    real, dimension(0:Nc,0:Nd)      :: Tcd
    real, dimension(0:Nd,0:Nd)      :: Thdd, Gdd, Gdd1

    !!$    real, dimension(Nd+1) :: Work
    !!$    integer, dimension(0:Nd)  :: ipiv
    !!$    integer                   :: info

    real                 :: dum
    integer                  :: m,ibc

    call makeT(Nc,Nd,Tcd)
    call makeTh(Nd,Nd,Thdd)
    call makeG(Nd,Gdd)

    D(:,:,0) = MATMUL(Tcd,Thdd)

    Gdd1 = Gdd

    do m = 1,NDmax
       D(:,:,m) = MATMUL(Tcd,MATMUL(Gdd,Thdd))
       Gdd = MATMUL(Gdd,Gdd1)
    end do

    END SUBROUTINE makeD

  SUBROUTINE makeTpDerVec(NN,p,x,bv)
    integer          :: NN
    integer          :: p  ! derivative order
    real             :: x  ! x = +/- 1
    real, dimension(0:NN) :: bv  ! coef of a_n

    integer          :: n

    do n = 0,NN
       call makeTpDer(NN,p,n,x,bv(n))
    end do

    END SUBROUTINE makeTpDerVec

  SUBROUTINE makeTpDer(NN,p,n,x,b)
    integer          :: NN
    integer          :: p  ! derivative order
    integer          :: n  ! as in T_n
    real             :: x  ! x = +/- 1
    real             :: b  ! coef of a_n

    integer          :: k

    b = (x)**(n+p)

    do k = 0,p-1
       b = b * (n**2 - k**2)/(2.*k+1)
    end do

    END SUBROUTINE makeTpDer

  SUBROUTINE identity(Np,A)
    integer                    :: Np
    real, dimension(0:Np,0:Np) :: A

    integer :: j

    A = 0.
    do j = 0,Np
       A(j,j) = 1.
    end do

    END SUBROUTINE identity

  SUBROUTINE makeG(Np,G)
    integer                        :: Np
    real, dimension(0:Np,0:Np) :: G

    integer  :: p,n

    do p = 0,Np
       do n = 0,Np
          if (p.ge.n.or.is_even(p+n)) then
             G(p,n) = 0.
          else
             G(p,n) = 2.*REAL(n)/c(Np,p)
          end if
       end do
    end do

    END SUBROUTINE makeG

  SUBROUTINE makeHs

    allocate(Hdr(0:Ndr,0:Ndr,0:4))

    call identity(Ndr,Hdr(0,0,0))
    call makeH(Ndr,Hdr(0,0,1))
    Hdr(:,:,2) = MATMUL(Hdr(:,:,1),Hdr(:,:,1))
    Hdr(:,:,3) = MATMUL(Hdr(:,:,1),Hdr(:,:,2))
    Hdr(:,:,4) = MATMUL(Hdr(:,:,2),Hdr(:,:,2))

    END SUBROUTINE makeHs

  SUBROUTINE makeH(Np,H)
    integer                        :: Np
    real, dimension(0:Np,0:Np) :: H

    integer  :: p,n

    H = 0.
    do p = 0,Np
       do n = 0,Np
          if (ABS(p-n).eq.1) then
             H(p,n) = 0.5
          else if (p.eq.n) then
             H(p,n) = 1.0
          end if
       end do
    end do
    H(1,0) = 1.
    H = H*(Rmax/2.)

    END SUBROUTINE makeH

  SUBROUTINE makeHinv(Np,H)
    integer                        :: Np
    real, dimension(0:Np,0:Np) :: H

    integer  :: p,n

    do p = 0,Np
       do n = 0,Np
          if (is_even(p+n)) then
             H(p,n) = 0.
          else if (n.lt.p .and. is_even(p+1)) then
             H(p,n) = 2*(-1)**((p-n-1)/2)
          else if (n.gt.p .and. is_even(p)) then
             H(p,n) = 2*(-1)**((p-n)/2)/c(Np,p)
          end if
       end do
    end do

    END SUBROUTINE makeHinv

  SUBROUTINE fromcheb2(MM,NN,u,a)
    integer                    :: MM,NN
    real, dimension(0:MM,0:NN) :: u,a

    real, dimension(0:MM,0:NN) :: ua
    real, dimension(0:NN,0:MM) :: uat,ut

    integer                    :: i,j,m,n

    do j = 0,NN
       call fromcheb1(MM,ua(0,j),a(0,j))
    end do
    uat = TRANSPOSE(ua)
    do i = 0,MM
       call fromcheb1(NN,ut(0,i),uat(0,i))
    end do
    u = TRANSPOSE(ut)

    END SUBROUTINE fromcheb2

  SUBROUTINE tocheb2(MM,NN,u,a)
    integer                    :: MM,NN
    real, dimension(0:MM,0:NN) :: u,a

    real, dimension(0:MM,0:NN) :: ua
    real, dimension(0:NN,0:MM) :: uat,at

    integer                    :: i,j,m,n

    do j = 0,NN
       call tocheb1(MM,u(0,j),ua(0,j))
    end do
    uat = TRANSPOSE(ua)
    do i = 0,MM
       call tocheb1(NN,uat(0,i),at(0,i))
    end do
    a = TRANSPOSE(at)

   END SUBROUTINE tocheb2

  SUBROUTINE tocheb1(NN,u,a)
    integer                :: NN
    real, dimension(0:NN)  :: u,a

    real, dimension(0:NN)  :: X
    integer                :: i,n

    call makeX(NN,X)

    a = 0.
    do i = 0,NN
       do n = 0,NN
          a(n) = a(n) + 2./REAL(NN)/c(NN,n)/c(NN,i)*u(i)*chebT(n,X(i))
       end do
    end do

    END SUBROUTINE tocheb1

  SUBROUTINE fromcheb1(NN,u,a)
    integer                :: NN
    real, dimension(0:NN)  :: u,a

    real, dimension(0:NN)  :: X
    integer                :: i,n

    call makeX(NN,X)

    u = 0.
    do i = 0,NN
       do n = 0,NN
          u(i) = u(i) + a(n)*chebT(n,X(i))
       end do
    end do

    END SUBROUTINE fromcheb1

  SUBROUTINE makeT(Nc,Nd,T)
    integer                   :: Nc,Nd

    real, dimension(0:Nc,0:Nd) :: T

    real, dimension(0:Nc)  :: X
    integer :: n,j

    call makeX(Nc,X)

    do n = 0,Nd
       do j = 0,Nc
          T(j,n) = chebT(n,X(j))
       end do
    end do

    END SUBROUTINE makeT

  SUBROUTINE makeTh(Nc,Nd,Th)
    integer                   :: Nc,Nd

    real, dimension(0:Nc,0:Nd) :: Th

    real, dimension(0:Nd) :: X
    integer :: n,j

    call makeX(Nd,X)

    do n = 0,Nc
       do j = 0,Nd
          Th(n,j) = chebT(n,X(j))/c(Nc,n)/c(Nd,j)
       end do
    end do
    Th = 2.*Th/REAL(Nd)

    END SUBROUTINE makeTh



  FUNCTION c(Np,j) RESULT (cval)
    integer :: Np,j
    real    :: cval

    if (j.eq.0 .or. j.eq.Np) then
       cval = 2.
    else
       cval = 1.
    end if

    END FUNCTION c

  FUNCTION chebT(n,xx) RESULT (T)
    integer  :: n
    real     :: xx
    real     :: T

    !  T = COS(Pi*REAL(n*j)/REAL(Nc))
    T = COS(REAL(n)*ACOS(xx))

    END FUNCTION chebT

  FUNCTION is_even (i) RESULT (y)
    integer :: i
    logical :: y

    if (NINT(REAL(i)/2.0) .eq. i/2) then
       y = .TRUE.
    else
       y = .FALSE.
    end if

    END FUNCTION is_even

END MODULE cheb
