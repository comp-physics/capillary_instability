PROGRAM film

  USE prms
  USE cheb
  USE out
  USE ev

  IMPLICIT none

  real :: dum
  integer :: i,j

  call init
!!$  call testint; stop
!!$  call testTs!; stop
!!$  call testPoiMat; stop
!!$  call testders; stop
!!$  call testdersRT; stop
!!$  call testrhs; stop
!!$  call testChebXforms; stop
!!$  call testChebXforms2; stop
!!$  call testH; print *,"TEST H DONE"; stop
!!$  call testHR; print *,"TEST HR DONE"; stop
!!$  call testDeraR; print *,"TEST DER a R DONE"; stop

  call solve
  call finalize

CONTAINS

  SUBROUTINE init

    call initprms
    call initcheb

  END SUBROUTINE init

  SUBROUTINE solve

!    call poisson1D
!    call biharm1D
!    call poisson2D
!    call biharm2D
!    call Xpoisson1D
!    call poisson2Dcyl
!    call biharm2Dcyl !RTTEST
    call evsolve

  END SUBROUTINE solve

  SUBROUTINE finalize

    call closeout

  END SUBROUTINE finalize

  SUBROUTINE biharm2Dcyl

    real, dimension(0:Ndr,0:Ndt) :: u,a,rx,ra

    real, dimension(0:Ndr)       :: bvr,abr,ubr,bvr1
    real, dimension(0:Ndt)       :: bvt,abt,ubt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: BB
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrrGGrrrr,GGtttt,rrGGrrtt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrGGrrr
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrGGr,rGGrtt,rrGGrr,GGtt,rGGr

    integer, dimension((Ndr+1)*(Ndt+1))  :: iepiv
    integer :: info

    integer                      :: Np
    integer                      :: i,j,ii,jj

    Np = (Ndr+1)*(Ndt+1)

    ra = 0.


    call makeGG(4,0,4,rrrrGGrrrr)
    call makeGG(2,2,2,rrGGrrtt)
    call makeGG(0,4,0,GGtttt)
    call makeGG(3,0,3,rrrGGrrr)
    call makeGG(1,2,1,rGGrtt)
    call makeGG(2,0,2,rrGGrr)
    call makeGG(0,2,0,GGtt)
    call makeGG(1,0,1,rGGr)


    BB = rrrrGGrrrr + 2.*rrGGrrtt + GGtttt &
       + 2.*rrrGGrrr - 2.*rGGrtt - rrGGrr  + 4.*GGtt + rGGr


    call makeTpDerVec(Ndr,0,-1.,bvr)
    ubt = 0.
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+1,j,:,j) = bvr
       ra(Ncr+1,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,0,1.,bvr)
    ubt = 0.
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+2,j,:,j) = bvr
       ra(Ncr+2,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,1,-1.,bvr)
    ubt = 0.
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+3,j,:,j) = bvr*(2./Rmax)
       ra(Ncr+3,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,2,1.,bvr)
    call makeTpDerVec(Ndr,1,1.,bvr1)
    ubt = SIN(2.*Pi*Td/Tmax)    !!! SURFACE TRACTION
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+4,j,:,j) = bvr*(2./Rmax)**2 - bvr1*(2./Rmax)/Rmax
       ra(Ncr+4,j) = abt(j)
    end do




    call makeTpDerVec(Ndt,0,-1.,bvt)
    ubr = 0.
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+1,i,:) = bvt
       ra(i,Nct+1) = abr(i)
    end do

    call makeTpDerVec(Ndt,0,1.,bvt)
    ubr =  0.
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+2,i,:) = bvt
       ra(i,Nct+2) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,-1.,bvt)
    ubr = 0.
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+3,i,:) = bvt*(2./Tmax)
       ra(i,Nct+3) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,1.,bvt)
    ubr =  0.
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+4,i,:) = bvt*(2./Tmax)
       ra(i,Nct+4) = abr(i)
    end do

    call DGETRF(Np,Np,BB(0,0,0,0),Np,iepiv,info)
    call DGETRS('N',Np,1,BB(0,0,0,0),Np,iepiv,ra(0,0),Np,info)

    call fromCheb2(Ndr,Ndt,u,ra)

    !!$    print "('         ',20F7.2)",Xdt
    !!$    print *,"NUM"
    !!$    do i = 0,Ndr
    !!$       print "(F7.2,'  ',20F7.2)",Xdr(i),u(i,:)
    !!$    end do

    call writefield(u)
    call writefieldF(u)
    call writevel(u)

    END SUBROUTINE biharm2Dcyl

  SUBROUTINE biharm2DcylRTTEST

    real, dimension(0:Ndr,0:Ndt) :: u,a,rx,ra,aa,ys

    real, dimension(0:Ndr)       :: bvr,abr,ubr,bvr1
    real, dimension(0:Ndt)       :: bvt,abt,ubt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: BB
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrrGGrrrr,GGtttt,rrGGrrtt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrGGrrr
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrGGr,rGGrtt,rrGGrr,GGtt,rGGr

    real, dimension((Ndr+1)*(Ndt+1)) :: yD
    real, dimension((Ndr+1)*(Ndt+1),(Ndr+1)*(Ndt+1)) :: BD

    integer, dimension((Ndr+1)*(Ndt+1))  :: iepiv
    integer :: info

    integer                      :: Np
    integer                      :: i,j,ii,jj

    Np = (Ndr+1)*(Ndt+1)

    do j = 0,Ndt
       do i = 0,Ndr
          rx(i,j) = &
               Rd(i)**4*(120.*Rd(i)*Td(j) &
               + 680.*Rd(i)*Td(j)**3 &
               + 225.*Rd(i)*Td(j)**5)
       end do
    end do

    call tocheb2(Ndr,Ndt,rx,ra)

    call makeGG(4,0,4,rrrrGGrrrr)
    call makeGG(2,2,2,rrGGrrtt)
    call makeGG(0,4,0,GGtttt)
    call makeGG(3,0,3,rrrGGrrr)
    call makeGG(1,2,1,rGGrtt)
    call makeGG(2,0,2,rrGGrr)
    call makeGG(0,2,0,GGtt)
    call makeGG(1,0,1,rGGr)


    BB = rrrrGGrrrr + 2.*rrGGrrtt + GGtttt &
       + 2.*rrrGGrrr - 2.*rGGrtt - rrGGrr  + 4.*GGtt + rGGr

    call makeTpDerVec(Ndr,0,-1.,bvr)
    ubt = Rd(Ndr)**5*Td**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+1,j,:,j) = bvr
       ra(Ncr+1,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,0,1.,bvr)
    ubt = Rd(0)**5*Td**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+2,j,:,j) = bvr
       ra(Ncr+2,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,1,-1.,bvr)
    ubt = 5.*Rd(Ndr)**4*Td**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+3,j,:,j) = bvr*(2./Rmax)
       ra(Ncr+3,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,2,1.,bvr)
    call makeTpDerVec(Ndr,1,1.,bvr1)
    ubt = 5.*4.*Rd(0)**3*Td**5 - 5.*Rd(0)**4*Td**5/Rmax
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+4,j,:,j) = bvr*(2./Rmax)**2 - bvr1*(2./Rmax)/Rmax
       ra(Ncr+4,j) = abt(j)
    end do




    call makeTpDerVec(Ndt,0,-1.,bvt)
    ubr = Rd**5*Td(Ndt)**5
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+1,i,:) = bvt
       ra(i,Nct+1) = abr(i)
    end do

    call makeTpDerVec(Ndt,0,1.,bvt)
    ubr =  Rd**5*Td(0)**5
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+2,i,:) = bvt
       ra(i,Nct+2) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,-1.,bvt)
    ubr = 5.*Rd**5*Td(Ndt)**4
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+3,i,:) = bvt*(2./Tmax)
       ra(i,Nct+3) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,1.,bvt)
    ubr =  5.*Rd**5*Td(0)**4
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+4,i,:) = bvt*(2./Tmax)
       ra(i,Nct+4) = abr(i)
    end do

    do j = 0,Ndt
       do i = 0,Ndr
          print "(7(5F7.3,'  '))",BB(i,j,:,:),ra(i,j)
       end do
       print *
    end do



    call DGETRF(Np,Np,BB(0,0,0,0),Np,iepiv,info)
    call DGETRS('N',Np,1,BB(0,0,0,0),Np,iepiv,ra(0,0),Np,info)



    call fromCheb2(Ndr,Ndt,u,ra)

    print "('         ',20F7.2)",Td
    print *,"NUM"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Rd(i),u(i,:)
    end do
    print *,"EX"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Rd(i),Rd(i)**5*Td(:)**5
    end do
    print *
    print *,"ERR"
    do i = 0,Ndr
       print "(F7.2,'  ',20E9.1)",Rd(i),u(i,:)-Rd(i)**5*Td(:)**5
    end do
    do j = 0,Ndt
       do i = 0,Ndr
          ys(i,j) = Rd(i)**5*Td(j)**5
       end do
    end do
    print *,"MAXERR = ",MAXVAL(ABS(u-ys))

    END SUBROUTINE biharm2DcylRTTEST

  SUBROUTINE biharm2DcylTEST

    real, dimension(0:Ndr,0:Ndt) :: u,a,rx,ra,aa,ys
    !!$    real, dimension(0:Ndr,0:Ndr,4) :: Gr
    !!$    real, dimension(0:Ndt,0:Ndt,4) :: Gt

    real, dimension(0:Ndr)       :: bvr,abr,ubr
    real, dimension(0:Ndt)       :: bvt,abt,ubt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: BB
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrrGGrrrr,GGtttt,rrGGrrtt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrGGrrr
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: rrrGGr,rGGrtt,rrGGrr,GGtt,rGGr

    real, dimension((Ndr+1)*(Ndt+1)) :: yD
    real, dimension((Ndr+1)*(Ndt+1),(Ndr+1)*(Ndt+1)) :: BD

    integer, dimension((Ndr+1)*(Ndt+1))  :: iepiv
    integer :: info

    integer                      :: Np
    integer                      :: i,j,ii,jj

    Np = (Ndr+1)*(Ndt+1)

    do j = 0,Ndt
       do i = 0,Ndr
          rx(i,j) = &
               Xdr(i)**4*(120.*Xdr(i)*Xdt(j) &
               + 680.*Xdr(i)*Xdt(j)**3 &
               + 225.*Xdr(i)*Xdt(j)**5)
       end do
    end do

    call tocheb2(Ndr,Ndt,rx,ra)

    call makeGG(4,0,4,rrrrGGrrrr)
    call makeGG(2,2,2,rrGGrrtt)
    call makeGG(0,4,0,GGtttt)
    call makeGG(3,0,3,rrrGGrrr)
    call makeGG(1,2,1,rGGrtt)
    call makeGG(2,0,2,rrGGrr)
    call makeGG(0,2,0,GGtt)
    call makeGG(1,0,1,rGGr)


    BB = rrrrGGrrrr + 2.*rrGGrrtt + GGtttt &
       + 2.*rrrGGrrr - 2.*rGGrtt - rrGGrr  + 4.*GGtt + rGGr


    call makeTpDerVec(Ndr,0,-1.,bvr)
    ubt = Xdr(Ndr)**5*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+1,j,:,j) = bvr
       ra(Ncr+1,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,0,1.,bvr)
    ubt = Xdr(0)**5*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+2,j,:,j) = bvr
       ra(Ncr+2,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,1,-1.,bvr)
    ubt = 5.*Xdr(Ndr)**4*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+3,j,:,j) = bvr
       ra(Ncr+3,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,2,1.,bvr)
    ubt = 5.*4.*Xdr(0)**3*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+4,j,:,j) = bvr
       ra(Ncr+4,j) = abt(j)
    end do




    call makeTpDerVec(Ndt,0,-1.,bvt)
    ubr = Xdr**5*Xdt(Ndt)**5
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+1,i,:) = bvt
       ra(i,Nct+1) = abr(i)
    end do

    call makeTpDerVec(Ndt,0,1.,bvt)
    ubr =  Xdr**5*Xdt(0)**5
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+2,i,:) = bvt
       ra(i,Nct+2) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,-1.,bvt)
    ubr = 5.*Xdr**5*Xdt(Ndt)**4
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+3,i,:) = bvt
       ra(i,Nct+3) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,1.,bvt)
    ubr =  5.*Xdr**5*Xdt(0)**4
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+4,i,:) = bvt
       ra(i,Nct+4) = abr(i)
    end do

    do j = 0,Ndt
       do i = 0,Ndr
          print "(7(5F7.3,'  '))",BB(i,j,:,:),ra(i,j)
       end do
       print *
    end do


    do jj = 0,Ndt
       do ii = 0,Ndt
          do j = 0,Ndr
             do i = 0,Ndr
                BD(1+i+ii*(Ndr+1),1+j+jj*(Ndr+1)) = BB(i,ii,j,jj)
             end do
          end do
       end do
    end do


    call DGETRF(Np,Np,BB(0,0,0,0),Np,iepiv,info)
    call DGETRS('N',Np,1,BB(0,0,0,0),Np,iepiv,ra(0,0),Np,info)



    call fromCheb2(Ndr,Ndt,u,ra)

    print "('         ',20F7.2)",Xdt
    print *,"NUM"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),u(i,:)
    end do
    print *,"EX"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),Xdr(i)**5*Xdt(:)**5
    end do
    print *
    print *,"ERR"
    do i = 0,Ndr
       print "(F7.2,'  ',20E9.1)",Xdr(i),u(i,:)-Xdr(i)**5*Xdt(:)**5
    end do
    do j = 0,Ndt
       do i = 0,Ndr
          ys(i,j) = Xdr(i)**5*Xdt(j)**5
       end do
    end do
    print *,"MAXERR = ",MAXVAL(u-ys)

    END SUBROUTINE biharm2DcylTEST

  SUBROUTINE biharm2D

    real, dimension(0:Ndr,0:Ndt) :: u,a,rx,ra,aa,ys
    !!$    real, dimension(0:Ndr,0:Ndr,4) :: Gr
    !!$    real, dimension(0:Ndt,0:Ndt,4) :: Gt

    real, dimension(0:Ndr)       :: bvr,abr,ubr
    real, dimension(0:Ndt)       :: bvt,abt,ubt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: BB
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: GGrrrr,GGtttt,GGrrtt

    real, dimension((Ndr+1)*(Ndt+1)) :: yD
    real, dimension((Ndr+1)*(Ndt+1),(Ndr+1)*(Ndt+1)) :: BD

    integer, dimension((Ndr+1)*(Ndt+1))  :: iepiv
    integer :: info

    integer                      :: Np
    integer                      :: i,j,ii,jj

    Np = (Ndr+1)*(Ndt+1)

    do j = 0,Ndt
       do i = 0,Ndr
          rx(i,j) = 120.*Xdr(i)**5*Xdt(j) &
               + 800.*Xdr(i)**3*Xdt(j)**3 &
               + 120.*Xdr(i)*Xdt(j)**5
       end do
    end do

    call tocheb2(Ndr,Ndt,rx,ra)


    call makeGG(4,0,0,GGrrrr)
    call makeGG(2,2,0,GGrrtt)
    call makeGG(0,4,0,GGtttt)

    BB = GGrrrr+2.*GGrrtt+GGtttt


    call makeTpDerVec(Ndr,0,-1.,bvr)
    ubt = Xdr(Ndr)**5*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+1,j,:,j) = bvr
       ra(Ncr+1,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,0,1.,bvr)
    ubt = Xdr(0)**5*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+2,j,:,j) = bvr
       ra(Ncr+2,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,1,-1.,bvr)
    ubt = 5.*Xdr(Ndr)**4*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+3,j,:,j) = bvr
       ra(Ncr+3,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,1,1.,bvr)
    ubt = 5.*Xdr(0)**4*Xdt**5
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+4,j,:,j) = bvr
       ra(Ncr+4,j) = abt(j)
    end do




    call makeTpDerVec(Ndt,0,-1.,bvt)
    ubr = Xdr**5*Xdt(Ndt)**5
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+1,i,:) = bvt
       ra(i,Nct+1) = abr(i)
    end do

    call makeTpDerVec(Ndt,0,1.,bvt)
    ubr =  Xdr**5*Xdt(0)**5
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+2,i,:) = bvt
       ra(i,Nct+2) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,-1.,bvt)
    ubr = 5.*Xdr**5*Xdt(Ndt)**4
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+3,i,:) = bvt
       ra(i,Nct+3) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,1.,bvt)
    ubr =  5.*Xdr**5*Xdt(0)**4
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+4,i,:) = bvt
       ra(i,Nct+4) = abr(i)
    end do

    do j = 0,Ndt
       do i = 0,Ndr
          print "(7(5F7.3,'  '))",BB(i,j,:,:),ra(i,j)
       end do
       print *
    end do


    do jj = 0,Ndt
       do ii = 0,Ndt
          do j = 0,Ndr
             do i = 0,Ndr
                BD(1+i+ii*(Ndr+1),1+j+jj*(Ndr+1)) = BB(i,ii,j,jj)
             end do
          end do
       end do
    end do


    call DGETRF(Np,Np,BB(0,0,0,0),Np,iepiv,info)
    call DGETRS('N',Np,1,BB(0,0,0,0),Np,iepiv,ra(0,0),Np,info)



    call fromCheb2(Ndr,Ndt,u,ra)

    print "('         ',20F7.2)",Xdt
    print *,"NUM"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),u(i,:)
    end do
    print *,"EX"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),Xdr(i)**5*Xdt(:)**5
    end do
    print *
    print *,"ERR"
    do i = 0,Ndr
       print "(F7.2,'  ',20E9.1)",Xdr(i),u(i,:)-Xdr(i)**5*Xdt(:)**5
    end do
    do j = 0,Ndt
       do i = 0,Ndr
          ys(i,j) = Xdr(i)**5*Xdt(j)**5
       end do
    end do
    print *,"MAXERR = ",MAXVAL(u-ys)

    END SUBROUTINE biharm2D

  SUBROUTINE poisson2D
    real, dimension(0:Ndr,0:Ndt) :: u,a,rx,ra,aa,ys

    real, dimension(0:Ndr)       :: bvr,abr,ubr
    real, dimension(0:Ndt)       :: bvt,abt,ubt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: BB,GGrr,GGtt

    real, dimension((Ndr+1)*(Ndt+1)) :: yD
   real, dimension((Ndr+1)*(Ndt+1),(Ndr+1)*(Ndt+1)) :: BD

    integer, dimension((Ndr+1)*(Ndt+1))  :: iepiv
    integer :: info
    integer                      :: i,j,ii,jj

    do j = 0,Ndt
       do i = 0,Ndr
          rx(i,j) = 2.*(Xdr(i)**2 + Xdt(j)**2)
       end do
    end do

    call tocheb2(Ndr,Ndt,rx,ra)

    call makeGG(2,0,0,GGrr)
    call makeGG(0,2,0,GGtt)
    BB = GGrr + GGtt

    call makeTpDerVec(Ndr,0,-1.,bvr)
    ubt = Xdr(Ndr)**2*Xdt**2
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+1,j,:,j) = bvr
       ra(Ncr+1,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,1,1.,bvr)
    ubt = 2.*Xdr(0)*Xdt**2
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+2,j,:,j) = bvr
       ra(Ncr+2,j) = abt(j)
    end do


    call makeTpDerVec(Ndt,0,-1.,bvt)
    ubr = Xdr**2*Xdt(Ndt)**2
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+1,i,:) = bvt
       ra(i,Nct+1) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,1.,bvt)
    ubr =  2.*Xdr**2*Xdt(0)
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+2,i,:) = bvt
       ra(i,Nct+2) = abr(i)
    end do

    do j = 0,Ndt
       do i = 0,Ndr
          print "(7(5F7.3,'  '))",BB(i,j,:,:),ra(i,j)
       end do
       print *
    end do




    call DGETRF((Ndr+1)*(Ndt+1),(Ndr+1)*(Ndt+1),BB(0,0,0,0),(Ndr+1)*(Ndt+1),iepiv,info)
    call DGETRS('N',(Ndr+1)*(Ndt+1),1,BB(0,0,0,0),(Ndr+1)*(Ndt+1),iepiv,ra(0,0),(Ndr+1)*(Ndt+1),info)




    call fromCheb2(Ndr,Ndt,u,ra)

    print "('         ',20F7.2)",Xdt
    print *,"NUM"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),u(i,:)
    end do
    print *,"EX,"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),Xdr(i)**2*Xdt(:)**2
    end do
    print *
    print *,"ERR"
    do i = 0,Ndr
       print "(F7.2,'  ',20E9.1)",Xdr(i),u(i,:)-Xdr(i)**2*Xdt(:)**2
    end do
    do j = 0,Ndt
       do i = 0,Ndr
          ys(i,j) = Xdr(i)**2*Xdt(j)**2
       end do
    end do
    print *,"MAXERR = ",MAXVAL(u-ys)



    END SUBROUTINE poisson2D

  SUBROUTINE poisson2Dcyl

    real, dimension(0:Ndr,0:Ndt) :: u,a,rx,ra,aa,ys

    real, dimension(0:Ndr)       :: bvr,abr,ubr
    real, dimension(0:Ndt)       :: bvt,abt,ubt
    real, dimension(0:Ndr,0:Ndt,0:Ndr,0:Ndt) :: BB,rrGGrr,rGGr,GGtt

    real, dimension((Ndr+1)*(Ndt+1)) :: yD
   real, dimension((Ndr+1)*(Ndt+1),(Ndr+1)*(Ndt+1)) :: BD

    integer, dimension((Ndr+1)*(Ndt+1))  :: iepiv
    integer :: info
    integer                      :: i,j,ii,jj

    do j = 0,Ndt
       do i = 0,Ndr
    !          rx(i,j) = 2.*(Xdr(i)**2 + Xdt(j)**2)

          rx(i,j) = Xdr(i)**2*(4.*Xdt(j)**2 + 2.)
       end do
    end do

    call tocheb2(Ndr,Ndt,rx,ra)

    call makeGG(2,0,2,rrGGrr)
    call makeGG(1,0,1,rGGr)
    call makeGG(0,2,0,GGtt)
    BB = rrGGrr + rGGr + GGtt

    call makeTpDerVec(Ndr,0,-1.,bvr)
    ubt = Xdr(Ndr)**2*Xdt**2
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+1,j,:,j) = bvr
       ra(Ncr+1,j) = abt(j)
    end do

    call makeTpDerVec(Ndr,1,1.,bvr)
    ubt = 2.*Xdr(0)*Xdt**2
    call tocheb1(Ndt,ubt,abt)
    do j = 0,Nct
       BB(Ncr+2,j,:,j) = bvr
       ra(Ncr+2,j) = abt(j)
    end do


    call makeTpDerVec(Ndt,0,-1.,bvt)
    ubr = Xdr**2*Xdt(Ndt)**2
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+1,i,:) = bvt
       ra(i,Nct+1) = abr(i)
    end do

    call makeTpDerVec(Ndt,1,1.,bvt)
    ubr =  2.*Xdr**2*Xdt(0)
    call tocheb1(Ndt,ubr,abr)
    do i = 0,Ndr
       BB(i,Nct+2,i,:) = bvt
       ra(i,Nct+2) = abr(i)
    end do

    do j = 0,Ndt
       do i = 0,Ndr
          print "(7(5F7.3,'  '))",BB(i,j,:,:),ra(i,j)
       end do
       print *
    end do


    call DGETRF((Ndr+1)*(Ndt+1),(Ndr+1)*(Ndt+1),BB(0,0,0,0),(Ndr+1)*(Ndt+1),iepiv,info)
    call DGETRS('N',(Ndr+1)*(Ndt+1),1,BB(0,0,0,0),(Ndr+1)*(Ndt+1),iepiv,ra(0,0),(Ndr+1)*(Ndt+1),info)


    call fromCheb2(Ndr,Ndt,u,ra)

    print "('         ',20F7.2)",Xdt
    print *,"NUM"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),u(i,:)
    end do
    print *,"EX,"
    do i = 0,Ndr
       print "(F7.2,'  ',20F7.2)",Xdr(i),Xdr(i)**2*Xdt(:)**2
    end do
    print *
    print *,"ERR"
    do i = 0,Ndr
       print "(F7.2,'  ',20E9.1)",Xdr(i),u(i,:)-Xdr(i)**2*Xdt(:)**2
    end do
    do j = 0,Ndt
       do i = 0,Ndr
          ys(i,j) = Xdr(i)**2*Xdt(j)**2
       end do
    end do
    print *,"MAXERR = ",MAXVAL(u-ys)


    END SUBROUTINE poisson2Dcyl

  SUBROUTINE Xpoisson1D

    real, dimension(0:Ndr) :: u,a,rx,ra

    real, dimension(0:Ndr)       :: bv
    real, dimension(0:Ndr,0:Ndr) :: BB,H

    integer, dimension(0:Ndr)  :: iepiv
    integer :: info
    integer                      :: i
    ! x^2 d^2 u / d x^2 == 12 x^4
    rx = (12.+4.)*Xdr**4

    call tocheb1(Ndr,rx,ra)

    call makeH(Ndr,H)
    BB = MATMUL(H,MATMUL(H,Gr(:,:,2))) + MATMUL(H,Gr(:,:,1))

    BB(Ncr+1:Ndr,:) = 0.

    call makeTpDerVec(Ndr,0,-1.,bv)
    BB(Ncr+1,:) = bv
    !  print *,bv

    call makeTpDerVec(Ndr,1,1.,bv)
    BB(Ncr+2,:) = bv
    !  print *,bv

    ra(Ncr+1) = 1.
    ra(Ncr+2) = 4.

    call DGETRF(Ndr+1,Ndr+1,BB(0,0),Ndr+1,iepiv,info)
    call DGETRS('N',Ndr+1,1,BB(0,0),Ndr+1,iepiv,ra(0),Ndr+1,info)

    call fromCheb1(Ndr,u,ra)

    print *,"X, NUM, EX, DIFF"
    do i = 0,Ndr
       print *,Xdr(i),u(i),Xdr(i)**4,u(i)-Xdr(i)**4
    end do

    END SUBROUTINE Xpoisson1D

  SUBROUTINE poisson1D

    real, dimension(0:Ndr) :: u,a,rx,ra
    real, dimension(0:Ndr,0:Ndr,0:4) :: G

    real, dimension(0:Ndr)       :: bv
    real, dimension(0:Ndr,0:Ndr) :: BB

    integer, dimension(0:Ndr)  :: iepiv
    integer :: info
    integer                      :: i

    rx = 30.*Xdr**4

    call tocheb1(Ndr,rx,ra)
    print *,ra

    call makeG(Ndr,G(0,0,1))
    G(:,:,2) = MATMUL(G(:,:,1),G(:,:,1))
    G(:,:,3) = MATMUL(G(:,:,2),G(:,:,1))
    G(:,:,4) = MATMUL(G(:,:,2),G(:,:,2))

    BB(0:Ncr,:) = G(0:Ncr,:,2)

    call makeTpDerVec(Ndr,0,-1.,bv)
    BB(Ncr+1,:) = bv
    print *,bv

    call makeTpDerVec(Ndr,1,1.,bv)
    BB(Ncr+2,:) = bv
    print *,bv

    ra(Ncr+1) = 1.
    ra(Ncr+2) = 6.

    call DGETRF(Ndr+1,Ndr+1,BB(0,0),Ndr+1,iepiv,info)
    call DGETRS('N',Ndr+1,1,BB(0,0),Ndr+1,iepiv,ra(0),Ndr+1,info)

    call fromCheb1(Ndr,u,ra)

    print *,"X, NUM, EX, DIFF"
    do i = 0,Ndr
       print *,Xdr(i),u(i),Xdr(i)**6,u(i)-Xdr(i)**6
    end do

    END SUBROUTINE poisson1D

  SUBROUTINE biharm1D

    real, dimension(0:Ndr) :: u,a,rx,ra
    real, dimension(0:Ndr,0:Ndr,0:4) :: G

    real, dimension(0:Ndr)       :: bv
    real, dimension(0:Ndr,0:Ndr) :: BB

    integer, dimension(0:Ndr)  :: iepiv
    integer :: info
    integer                      :: i

    rx = 6.*5.*4.*3.*Xdr**2

    call tocheb1(Ndr,rx,ra)
    print *,ra

    call makeG(Ndr,G(0,0,1))
    G(:,:,2) = MATMUL(G(:,:,1),G(:,:,1))
    G(:,:,3) = MATMUL(G(:,:,2),G(:,:,1))
    G(:,:,4) = MATMUL(G(:,:,2),G(:,:,2))

    BB(0:Ncr,:) = G(0:Ncr,:,4)

    call makeTpDerVec(Ndr,0,-1.,bv)
    BB(Ncr+1,:) = bv
    ra(Ncr+1) = 1.
    print *,bv

    call makeTpDerVec(Ndr,1,1.,bv)
    BB(Ncr+2,:) = bv
    ra(Ncr+2) = 6.
    print *,bv

    call makeTpDerVec(Ndr,0,1.,bv)
    BB(Ncr+3,:) = bv
    ra(Ncr+3) = 1.
    print *,bv

    call makeTpDerVec(Ndr,2,-1.,bv)
    BB(Ncr+4,:) = bv
    ra(Ncr+4) = 30.
    print *,bv


    call DGETRF(Ndr+1,Ndr+1,BB(0,0),Ndr+1,iepiv,info)
    call DGETRS('N',Ndr+1,1,BB(0,0),Ndr+1,iepiv,ra(0),Ndr+1,info)

    call fromCheb1(Ndr,u,ra)

    print *,"X, NUM, EX, DIFF"
    do i = 0,Ndr
       print *,Xdr(i),u(i),Xdr(i)**6,u(i)-Xdr(i)**6
    end do

    END SUBROUTINE biharm1D

  SUBROUTINE testdersRT

    integer  :: j
    real, dimension(0:Ndr) :: f,df,dfe
    real, dimension(0:Ncr) :: dfc,dfce

    real, dimension(0:Ndt) :: ft,dft,dfet
    real, dimension(0:Nct) :: dfct,dfcet


    f = Rd**3

    df = MATMUL(Dddr(:,:,2),f)
    dfe(0:Ndr) = 3.*Rd
    print *,"DD-R"
    do j = 0,Ndr
       write(*,"(6E35.25)")Rd(j),df(j),dfe(j),df(j)-dfe(j)
    end do

    dfc = MATMUL(Dcdr(:,:,2),f)
    dfce(0:Ncr) = 3.*Rc
    print *,"CD-R"
    do j = 0,Ncr
       write(*,"(6E35.25)")Rc(j),dfc(j),dfce(j),dfc(j)-dfce(j)
    end do


    ft = Td**3

    dft = MATMUL(Dddt(:,:,2),ft)
    dfet(0:Ndt) = 6.*Td
    print *,"DD-T"
    do j = 0,Ndt
       write(*,"(6E35.25)")Td(j),dft(j),dfet(j),dft(j)-dfet(j)
    end do

    dfct = MATMUL(Dcdt(:,:,2),ft)
    dfcet(0:Nct) = 6.*Tc
    print *,"CD-T"
    do j = 0,Nct
       write(*,"(6E35.25)")Tc(j),dfct(j),dfcet(j),dfct(j)-dfcet(j)
    end do

    END SUBROUTINE testdersRT

  SUBROUTINE testders

    integer  :: j
    real, dimension(0:Ndr) :: f,df,dfe
    real, dimension(0:Ncr) :: dfc,dfce


    f = Xdr**2


    df = MATMUL(Dddr(:,:,1),f)
    dfe(0:Ndr) = 2.*Xdr
    print *,"DD-R"
    do j = 0,Ndr
       write(*,"(6E35.25)")Xdr(j),df(j),dfe(j),df(j)-dfe(j)
    end do

    dfc = MATMUL(Dcdr(:,:,1),f)
    dfce(0:Ncr) = 2.*Xcr
    print *,"CD-R"
    do j = 0,Ncr
       write(*,"(6E35.25)")Xcr(j),dfc(j),dfce(j),dfc(j)-dfce(j)
    end do


    f = Xdt**2

    df = MATMUL(Dddt(:,:,1),f)
    dfe(0:Ndr) = 2.*Xdt
    print *,"DD-T"
    do j = 0,Ndr
       write(*,"(6E35.25)")Xdt(j),df(j),dfe(j),df(j)-dfe(j)
    end do

    dfc = MATMUL(Dcdt(:,:,1),f)
    dfce(0:Nct) = 2.*Xct
    print *,"CD-T"
    do j = 0,Nct
       write(*,"(6E35.25)")Xct(j),dfc(j),dfce(j),dfc(j)-dfce(j)
    end do

    END SUBROUTINE testders

  SUBROUTINE testTs

    real, dimension(Ndr,Ndr)  :: Fdd
    real, dimension(Ncr,Ncr)  :: Fcc
    integer                   :: i

    Fdd = MATMUL(Tddr,Thddr)
    do i = 1,Ndr
       print "(20F10.2)", Fdd(i,:)
    end do

    print *

    Fcc = MATMUL(Tccr,Thccr)
    do i = 1,Ncr
       print "(20F10.2)", Fcc(i,:)
    end do

    print *

    Fcc = MATMUL(Tcdr,Thdcr)
    do i = 1,Ncr
       print "(20F10.2)", Fcc(i,:)
    end do


    END SUBROUTINE testTs

  SUBROUTINE TestChebXforms

    real, dimension(0:Ndr)  :: u,uo,a

    integer :: n,i

    u = Xdr**5
    uo = u

    call tocheb1(Ndr,u,a)

    print *,"n,a(n)"
    do n = 0,Ndr
       print *,n,a(n)
    end do

    call fromcheb1(Ndr,u,a)

    print *
    print *,"i,u,uo,uo-u"
    do i = 0,Ndr
       print *,i,u(i),uo(i),uo(i)-u(i)
    end do

    END SUBROUTINE TestChebXforms

  SUBROUTINE TestHR

    real, dimension(0:Ndr)  :: u,uo,a
    real, dimension(0:Ndr,0:Ndr) :: H,Hinv,HHinv

    integer, dimension(0:Ndr)  :: iepiv
    integer :: n,i,info

    u = Rd**6
    uo = u

    call tocheb1(Ndr,u,a)

    print *,"n,a(n)"
    do n = 0,Ndr
       print *,n,a(n)
    end do

    call makeH(Ndr,H)
    call makeHinv(Ndr,Hinv)

    HHinv = MATMUL(H,Hinv)

    do i = 0,Ndr
       write(*,"(20F6.2)") HHinv(i,:)
    end do

    print *

    !!$    do i = Ndr,1,-1
    !!$       H(i,:) = H(i-1,:)
    !!$    end do
    !!$    H(0,:) = 0.
    do i = 0,Ndr
       write(*,"(20F6.2)") H(i,:)
    end do

    a = MATMUL(H,a)
    a = MATMUL(H,a)

    !!$    call DGETRF(Ndr+1,Ndr+1,H(0,0),Ndr+1,iepiv,info)
    !!$    call DGETRS('N',Ndr+1,1,H(0,0),Ndr+1,iepiv,a(0),Ndr+1,info)
    !!$

    print *
    print *,"n,b(n)"
    do n = 0,Ndr
       print *,n,a(n)
    end do


    call fromcheb1(Ndr,u,a)

    print *
    print *,"i,x u,x uo,x uo-u"
    do i = 0,Ndr
       print *,i,u(i),uo(i)*Rd(i)**2,uo(i)*Rd(i)**2-u(i)
    end do

    END SUBROUTINE TestHR

  SUBROUTINE TestDeraR

    real, dimension(0:Ndr)  :: u,a
    real, dimension(0:Ndr)  :: ut,at
    real, dimension(0:Ndr,0:Ndr) :: H,Hinv,HHinv

    integer, dimension(0:Ndr)  :: iepiv
    integer :: n,i,info

    u = Rd**6

    call tocheb1(Ndr,u,a)

    print *,"n,a(n)"
    do n = 0,Ndr
       print *,n,a(n)
    end do

    a = MATMUL(Gr(:,:,4),a)

    print *

    print *
    print *,"n,b(n)"
    do n = 0,Ndr
       print *,n,a(n)
    end do

    call fromcheb1(Ndr,u,a)

    print *
    print *,"i,u'',ex, diff"
    do i = 0,Ndr
       print *,i,u(i),6.*5.*4.*3.*Rd(i)**2,6.*5.*4.*3.*Rd(i)**2-u(i)
    end do


    !!!!!!!!!!

    ut = Td**6

    call tocheb1(Ndr,ut,at)

    print *,"n,a(n)"
    do n = 0,Ndr
       print *,n,at(n)
    end do

    at = MATMUL(Gt(:,:,2),at)

    print *

    print *
    print *,"n,b(n)"
    do n = 0,Ndt
       print *,n,at(n)
    end do

    call fromcheb1(Ndt,ut,at)

    print *
    print *,"i,u'',ex, diff"
    do i = 0,Ndr
       print *,i,ut(i),6.*5.*Td(i)**4,6.*5.*Td(i)**4-ut(i)
    end do

    END SUBROUTINE TestDeraR

  SUBROUTINE TestH

    real, dimension(0:Ndr)  :: u,uo,a
    real, dimension(0:Ndr,0:Ndr) :: H,Hinv,HHinv

    integer, dimension(0:Ndr)  :: iepiv
    integer :: n,i,info

    u = Xdr**6
    uo = u

    call tocheb1(Ndr,u,a)

    print *,"n,a(n)"
    do n = 0,Ndr
       print *,n,a(n)
    end do

    call makeH(Ndr,H)
    call makeHinv(Ndr,Hinv)

    HHinv = MATMUL(H,Hinv)

    do i = 0,Ndr
       write(*,"(20F6.2)") HHinv(i,:)
    end do

    print *

    do i = 0,Ndr
       write(*,"(20F6.2)") H(i,:)
    end do

    a = MATMUL(H,a)
    a = MATMUL(H,a)

    !!$    call DGETRF(Ndr+1,Ndr+1,H(0,0),Ndr+1,iepiv,info)
    !!$    call DGETRS('N',Ndr+1,1,H(0,0),Ndr+1,iepiv,a(0),Ndr+1,info)
    !!$

    print *
    print *,"n,b(n)"
    do n = 0,Ndr
       print *,n,a(n)
    end do


    call fromcheb1(Ndr,u,a)

    print *
    print *,"i,x u,x uo,x uo-u"
    do i = 0,Ndr
       print *,i,u(i),uo(i)*Xdr(i)**2,uo(i)*Xdr(i)**2-u(i)
    end do

    END SUBROUTINE TestH

  SUBROUTINE TestChebXforms2

    real, dimension(0:Ndr,0:Ndt)  :: u,uo,a

    integer :: n,i

    do i = 0,Ndr
       u(i,:) = Xdr(i)**5*Xdt**4
    end do
    uo = u

    call tocheb2(Ndr,Ndt,u,a)

    print *,"n,a(n)"
    do n = 0,Ndr
       print *,n,a(n,:)
    end do

    call fromcheb2(Ndr,Ndt,u,a)

    print *
    print *,"i,u,uo,uo-u"
    do i = 0,Ndr
       print *,i,u(i,:)
       print *,i,uo(i,:)
       print *,i,uo(i,:)-u(i,:)
       print *
    end do

    END SUBROUTINE TestChebXforms2

END PROGRAM film

