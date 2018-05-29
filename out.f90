MODULE out

  USE prms
  USE data
  USE cheb

  IMPLICIT none


CONTAINS

  SUBROUTINE initout

  END SUBROUTINE initout

  SUBROUTINE closeout
    

  END SUBROUTINE closeout

  SUBROUTINE writefield(y)
    real, dimension(0:Ndr,0:Ndt) :: y

    integer                    :: i,j


    open(1,file='D/y.plt',form='FORMATTED')
    write(1, '(A)') 'VARIABLES = X, Y, Dil'    
    write(1, '(A,I9,A,I9,A)') 'ZONE I =', Ndr+1,' J=', Ndt+1, ' F=POINT'  
    do j = 0,Ndt
       do i = 0,Ndr
          write(1,"(E20.5)") Rd(i)*COS(Td(j)),Rd(i)*SIN(Td(j)),y(i,j)
       end do
    end do
    close(1)

  END SUBROUTINE writefield

  SUBROUTINE writevel(y)
    real, dimension(0:Ndr,0:Ndt) :: y

    real, dimension(0:Ndr,0:Ndt) :: u,v
    real, dimension(0:Nfr,0:Nft) :: uF,vF

    integer                    :: i,j

    do j = 0,Ndt
       v(:,j) = -MATMUL(Dddr(:,:,1),y(:,j))
    end do
    do i = 0,Ndr
       u(i,:) = MATMUL(Dddt(:,:,1),y(i,:))
    end do
!!$
!!$    uF = u; vF = v
    call makefine(u,uF)
    call makefine(v,vF)
    do j = 0,Nft
       uF(:,j) = uF(:,j)/Rf(:)
    end do

    open(1,file='D/uv.plt',form='FORMATTED')
    write(1, '(A)') 'VARIABLES = X, Y, U, V'    
    write(1, '(A,I9,A,I9,A)') 'ZONE I =', Nfr,' J=', Nft+1, ' F=POINT'  
    do j = 0,Nft
       do i = 0,Nfr-1
          write(1,"(E20.5)") Rf(i)*COS(Tf(j)), &
               Rf(i)*SIN(Tf(j)),uF(i,j)*COS(Tf(j))-vF(i,j)*SIN(Tf(j)), &
               uF(i,j)*SIN(Tf(j))+vF(i,j)*COS(Tf(j))
       end do
    end do
    close(1)


  END SUBROUTINE writevel

  SUBROUTINE writefieldF(y)
    real, dimension(0:Ndr,0:Ndt) :: y

    real, dimension(0:Nfr,0:Nft) :: yF

    integer                    :: i,j

    call makefine(y,yF)

    open(1,file='D/yF.plt',form='FORMATTED')
    write(1, '(A)') 'VARIABLES = X, Y, Dil'    
    write(1, '(A,I9,A,I9,A)') 'ZONE I =', Nfr+1,' J=', Nft+1, ' F=POINT'  
    do j = 0,Nft
       do i = 0,Nfr
          write(1,"(E20.5)") Rf(i)*COS(Tf(j)),Rf(i)*SIN(Tf(j)),yF(i,j)
       end do
    end do
    close(1)


  END SUBROUTINE writefieldF

  SUBROUTINE makefine(y,yF)
    real, dimension(0:Ndr,0:Ndt) :: y
    real, dimension(0:Nfr,0:Nft) :: yF

    real, dimension(0:Nfr,0:Ndt) :: yi
    integer                          :: i,j

    do j = 0,Ndt
       yi(:,j) = MATMUL(Tfdr,MATMUL(Thddr,y(:,j)))
    end do
    do i = 0,Nfr
       yF(i,:) = MATMUL(Tfdt,MATMUL(Thddt,yi(i,:)))
    end do

  END SUBROUTINE makefine

  SUBROUTINE writeFDresid(y)
    real, dimension(0:Ndr,0:Ndt) :: y
    real, dimension(0:Nfr,0:Nft) :: yF

    real, dimension(0:Nfr,0:Nft) :: ytttt,yrrrr,yttrr,yrrr,yrr,yr,yrtt,ytt,res
    integer :: i,j

    call makefine(y,yF)

    call der(0,4,yF,ytttt)
    call der(4,0,yF,yrrrr)
    call der(2,2,yF,yttrr)
    call der(3,0,yF,yrrr)
    call der(2,0,yF,yrr)
    call der(1,0,yF,yr)
    call der(1,2,yF,yrtt)
    call der(0,2,yF,ytt)

    do i = 0,Nfr-1
       res(i,:) = yrrrr(i,:) + 2./Rf(i)**2*yttrr(i,:) &
                + 1./Rf(i)**4*ytttt(i,:) &
                + 2./Rf(i)*yrrr(i,:) - 2./Rf(i)**3*yrtt(i,:) &
                - 1./Rf(i)**2*yrr(i,:) + 4./Rf(i)**4*ytt(i,:) &
                + 1./Rf(i)**3*yr(i,:)
    end do

    open(1,file='D/res.plt',form='FORMATTED')
    write(1, '(A)') 'VARIABLES = X, Y, Dil'    
    write(1, '(A,I9,A,I9,A)') 'ZONE I =', Nfr,' J=', Nft+1, ' F=POINT'  
    do j = 0,Nft
       do i = 0,Nfr-1
          write(1,"(E20.5)") Rf(i)*COS(Tf(j)),Rf(i)*SIN(Tf(j)),res(i,j)
       end do
    end do
    close(1)

  END SUBROUTINE writeFDresid

  SUBROUTINE der(nr,nt,f,df)
    integer                       :: nr,nt
    real, dimension(0:Nfr,0:Nft)  :: f,df    

    integer :: l

    df = f
    print *,"nr = ",nr,"  nt = ",nt
    print *,"r"
    do l = 1,nr
       print *,l
       call derR(df)
    end do
    print *,"t"
    do l = 1,nt
       print *,l
       call derT(df)
    end do
    print *
    print *

  END SUBROUTINE der

  SUBROUTINE derR(df)
    real, dimension(0:Nfr,0:Nft)  :: df
    real, dimension(0:Nfr,0:Nft)  :: w

    integer  :: i

    w(0,:) = (df(1,:)-df(0,:))/(Rf(1)-Rf(0))
    do i = 1,Nfr-1
       w(i,:) = (df(i+1,:)-df(i-1,:))/(Rf(i+1)-Rf(i-1))
    end do
    w(Nfr,:) = (df(Nfr,:)-df(Nfr-1,:))/(Rf(Nfr)-Rf(Nfr-1))

    df = w

  END SUBROUTINE derR


  SUBROUTINE derT(df)
    real, dimension(0:Nfr,0:Nft)  :: df
    real, dimension(0:Nfr,0:Nft)  :: w

    integer  :: i

    w(:,0) = (df(:,1)-df(:,0))/(Tf(1)-Tf(0))
    do i = 1,Nft-1
       w(:,i) = (df(:,i+1)-df(:,i-1))/(Tf(i+1)-Tf(i-1))
    end do
    w(:,Nft) = (df(:,Nft)-df(:,Nft-1))/(Tf(Nft)-Tf(Nft-1))

    df = w

  END SUBROUTINE derT


END MODULE out
