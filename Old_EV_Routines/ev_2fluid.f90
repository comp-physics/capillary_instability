!3d two fluid newtonian pipe flow solver!

MODULE ev

  USE prms
  USE data
  USE cheb

  IMPLICIT none

INTEGER          NOUT
PARAMETER        (NOUT=6)
INTEGER          NB, NMAX, lt
PARAMETER        (NB=4000,NMAX=4000)
INTEGER          LDA, LDB, LDVR, LWORK
PARAMETER        (LDA=NMAX,LDB=NMAX,LDVR=NMAX,LWORK=NMAX+NMAX*NB)
COMPLEX *16      ii
PARAMETER        (ii=(0,1.0))

!FLOW PARAMETERS
    Real             GAM,RHO1,RHO2,MU1,MU2,NU1,NU2,aa,HH,AA1,AA2,KAP

    !take these out later!!!
    real              mu, rho, nu

    PARAMETER        (GAM=1.0,RHO1=1.0,RHO2=1.0,MU1=0.5,MU2=0.5,NU1=MU1/RHO1,NU2=MU2,KAP=-1.0)

    REAL          K
    integer       iter,iterMax

CONTAINS

SUBROUTINE evsolve

  INTEGER          INFO,LWKOPT, N, U, W
  COMPLEX *16      A(LDA,NMAX), ALPHA(NMAX), B(LDB,NMAX),BETA(NMAX), DUMMY(1,1), VR(LDVR,NMAX),WORK(LWORK)
  DOUBLE PRECISION RWORK(8*NMAX), SMALL, DLAMCH

  integer :: i,j
  character(30)          :: fn

  Real, allocatable, dimension(:,:,:)  :: xDddr,x2Dddr,xDddt,x2Dddt,x3Dddr,x3Dddt,x4Dddr,x4Dddt

  !stop if the number of collocation points in each domain is not the same
  If (Ndt .ne. Ndr) Stop

  allocate(xDddr(0:Ndr,0:Ndr,0:3))
  allocate(x2Dddr(0:Ndr,0:Ndr,0:3))
  allocate(x3Dddr(0:Ndr,0:Ndr,0:3))
  allocate(x4Dddr(0:Ndr,0:Ndr,0:3))

  allocate(xDddt(0:Ndr,0:Ndr,0:3))
  allocate(x2Dddt(0:Ndr,0:Ndr,0:3))
  allocate(x3Dddt(0:Ndr,0:Ndr,0:3))
  allocate(x4Dddt(0:Ndr,0:Ndr,0:3))

  lt=1
  N=6*(Ndr+1)

  aa=RMAX
  HH=RMAX+TMAX

  AA1 = (kap/4)*((aa**2)*(1/mu1-1/mu2)-(HH**2)*(1/mu2))
  AA2 = (-kap/4)*(HH**2)/mu2


  iterMAx=100

  !fix the sign error for odd derivatives
  Dddr(:,:,1) = -1*Dddr(:,:,1)
  Dddr(:,:,3) = -1*Dddr(:,:,3)

  Dddt(:,:,1) = -1*Dddt(:,:,1)
  Dddt(:,:,3) = -1*Dddt(:,:,3)

  !make the xDddr matrix that premultiplies Dddr(:) by XXdr
  Do i=0,3
  xDddr(:,:,i) = Matmul(XXdr(:,:),Dddr(:,:,i))
  xDddt(:,:,i) = Matmul(XXdt(:,:),Dddt(:,:,i))
  End Do

  !make the x2Dddr matrix that premultiplies Dddr(:) by x^2 vector matrix
  Do i=0,3
  x2Dddr(:,:,i) = Matmul(XX2dr(:,:),Dddr(:,:,i))
  x2Dddt(:,:,i) = Matmul(XX2dt(:,:),Dddt(:,:,i))
  End Do

  Do i=0,3
  x3Dddr(:,:,i) = Matmul(XX3dr(:,:),Dddr(:,:,i))
  x3Dddt(:,:,i) = Matmul(XX3dt(:,:),Dddt(:,:,i))
  End Do

  Do i=0,3
  x4Dddr(:,:,i) = Matmul(XX4dr(:,:),Dddr(:,:,i))
  x4Dddt(:,:,i) = Matmul(XX4dt(:,:),Dddt(:,:,i))
  End Do

  !write the xDddr matrix
  open(unit=13, file='graphxDddr.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(13, '(1600F9.2)')( Real(xDddr(i,j,1)) ,j=0,Ndr)
  end do
  close(13)

  !write the x2Dddr matrix
  open(unit=157, file='graphx2Dddr.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(157, '(1600F9.2)')( Real(x2Dddr(i,j,1)) ,j=0,Ndr)
  end do
  close(157)

  !looping over K
  Do iter=1,itermax

  k=Real(iter)/Real(itermax)
  print*, 'K =', k

  call buildA(A,xDddr,x2Dddr,xDddt,x2Dddt,x3Dddr,x3Dddt,x4Dddr,x4Dddt)
  call buildB(B,xDddr,x2Dddr,xDddt,x2Dddt,x3Dddr,x3Dddt,x4Dddr,x4Dddt)

  !lapack gevp solver

  IF (N.LE.NMAX) THEN
         CALL ZGGEV('No left vectors','Vectors (right)',N,A,LDA,B,LDB,&
     &              ALPHA,BETA,DUMMY,1,VR,LDVR,WORK,LWORK,RWORK,INFO)

            !write eigenvalues to files in directory D
            Do J=1,N
               write(fn,"('D/eval',I9.9)")iter
                open(1,file=fn)
                If (Real(BETA(j))>0.00001) Then
                  write(1,"(3F30.20)")  ALPHA(J)/BETA(J)
                End if
            End Do
               close(1)

            !write eigenvectors to files in directory V
            Do lt=1,N
              Do I=1,N
                If (Real(Beta(lt))>0.00001) Then
                  write(fn,"('V/x.',I9.9)")lt
                    open(2,file=fn)
                      write(2,"(3F20.10)")  VR(I,lt)
                End if
              End Do
            End Do

            close(2)

         IF (INFO.GT.0) THEN
            WRITE (NOUT,99999) 'Failure in ZGGEV. INFO =', INFO
         ELSE
            SMALL = DLAMCH('Sfmin')
            DO 20 J = 1,N
               IF ((ABS(ALPHA(J)))*SMALL.GE.ABS(BETA(J))) THEN
                  WRITE (NOUT,99998) 'Eigenvalue(', J, ') is numerically infinite or undetermined'
               ELSE
               !print eigenvalues to the shell

               WRITE (NOUT,99997) 'Eigenvalue(', J, ') = ',ALPHA(J)/BETA(J)
               END IF
               !print eigenvectors to the shell

              !WRITE (NOUT,*)
              !WRITE (NOUT,99996) 'Eigenvector(', J, ')',&
     !&           (VR(I,J),I=1,N)
   20       CONTINUE

            LWKOPT = WORK(1)
            IF (LWORK.LT.LWKOPT) THEN
               WRITE (NOUT,99995) 'Optimum workspace required = ',&
     &           LWKOPT, 'Workspace provided         = ', LWORK
            END IF
         END IF
      ELSE
         WRITE (NOUT,*) 'NMAX too small'
      END IF

      99999 FORMAT (1X,A,I4)
      99998 FORMAT (1X,A,I2,2A,/1X,2(A,I2,A,'(',1P,E11.4,',',1P,E11.4,')'))
      99997 FORMAT (1X,A,I2,A,'(',1P,E11.4,',',1P,E11.4,')')
      99996 FORMAT (1X,A,I2,A,/3(1X,'(',1P,E11.4,',',1P,E11.4,')',:))
      99995 FORMAT (1X,A,I5,/1X,A,I5)

  End Do

  END SUBROUTINE evsolve

SUBROUTINE buildA(A,xDddr,x2Dddr,xDddt,x2Dddt,x3Dddr,x3Dddt,x4Dddr,x4Dddt)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: A11, A12, A13, A14, A15, A16, A21, A22, A23, A24, A25, A26, A31, A32, A33, A34, A35, A36
  COMPLEX, allocatable, dimension(:,:)  :: A41, A42, A43, A44, A45, A46, A51, A52, A53, A54, A55, A56, A61, A62, A63, A64, A65, A66
  Real, allocatable, dimension(:,:,:)   :: xDddr,x2Dddr,xDddt,x2Dddt,x3Dddr,x3Dddt,x4Dddr,x4Dddt
  COMPLEX *16                           :: A(LDA,NMAX)

  allocate(A11(0:Ndr,0:Ndr))
  allocate(A12(0:Ndr,0:Ndr))
  allocate(A13(0:Ndr,0:Ndr))
  allocate(A14(0:Ndr,0:Ndr))
  allocate(A15(0:Ndr,0:Ndr))
  allocate(A16(0:Ndr,0:Ndr))

  allocate(A21(0:Ndr,0:Ndr))
  allocate(A22(0:Ndr,0:Ndr))
  allocate(A23(0:Ndr,0:Ndr))
  allocate(A24(0:Ndr,0:Ndr))
  allocate(A25(0:Ndr,0:Ndr))
  allocate(A26(0:Ndr,0:Ndr))

  allocate(A31(0:Ndr,0:Ndr))
  allocate(A32(0:Ndr,0:Ndr))
  allocate(A33(0:Ndr,0:Ndr))
  allocate(A34(0:Ndr,0:Ndr))
  allocate(A35(0:Ndr,0:Ndr))
  allocate(A36(0:Ndr,0:Ndr))

  allocate(A41(0:Ndr,0:Ndr))
  allocate(A42(0:Ndr,0:Ndr))
  allocate(A43(0:Ndr,0:Ndr))
  allocate(A44(0:Ndr,0:Ndr))
  allocate(A45(0:Ndr,0:Ndr))
  allocate(A46(0:Ndr,0:Ndr))

  allocate(A51(0:Ndr,0:Ndr))
  allocate(A52(0:Ndr,0:Ndr))
  allocate(A53(0:Ndr,0:Ndr))
  allocate(A54(0:Ndr,0:Ndr))
  allocate(A55(0:Ndr,0:Ndr))
  allocate(A56(0:Ndr,0:Ndr))

  allocate(A61(0:Ndr,0:Ndr))
  allocate(A62(0:Ndr,0:Ndr))
  allocate(A63(0:Ndr,0:Ndr))
  allocate(A64(0:Ndr,0:Ndr))
  allocate(A65(0:Ndr,0:Ndr))
  allocate(A66(0:Ndr,0:Ndr))



  !FLUID 1 EQUATIONS

  !R-MOMENTUM

  !A11
  A11(:,:) = MU1 * x2Dddr(:,:,2) + MU1 * xDddr(:,:,1) - MU1*(K**2)*x2Dddr(:,:,0) &
   & - MU1*Dddr(:,:,0) - RHO1*ii*K*(KAP/(4*MU1))*Dddr(:,:,0) - RHO1*ii*K*AA1*x2Dddr(:,:,0)

      !BC for singularity condition on Ur at r=0
      A11(0,:) = Dddr(0,:,0)

      !BC for stress condition at r=aa
      !!!! check the sign here !!!! this reflects page 15 of notes
      A11(Ndr,:) = (-1)*Dddr(Ndr,:,0)*GAM*((K**2)-1./(aa**2))

  !A12
  A12(:,:) = 0

      !BC for singularity on Ur
      A12(0,:) = 0

      !BC for stress at r=aa
      A12(Ndr,:) = 0

  !A13
  A13(:,:) = (-1)*x2Dddr(:,:,1)

      !BC singularity condition on Ur
      A13(0,:) = 0

      !BC for stress at r=aa
      A13(Ndr,:) = 0

  !A14
  A14(:,:) = 0

  !A15
  A15(:,:) = 0

  !A16
  A16(:,:) = 0

  !X-MOMENTUM

  !A21
  A21(:,:) = -RHO1*(KAP/(MU1*2))*x2Dddr(:,:,0)

      !BC for singularity on Ux
      A21(0,:) = 0

      !tang stress BC at r=a
      A21(Ndr,:) = (-1)*MU1*ii*K*Dddr(Ndr,:,0)

  !A22
  A22(:,:) = MU1*xDddr(:,:,2) + MU1*xDddr(:,:,1) - K*MU1*xDddr(:,:,0) - &
   & RHO1*ii*K*(KAP/(MU1*4))*x3Dddr(:,:,0) - RHO1*ii*K*AA1*xDddr(:,:,0)

      !BC singularity condition on Ux
      A22(0,:) = Dddr(0,:,1)

      !BC for tang stress at r=a
      A22(Ndr,:) = (-1)*MU1*Dddr(Ndr,:,1)


  !A23
  A23(:,:) = (-1)*ii*K*xDddr(:,:,0)

      !BC singularity condition on Ux
      A23(0,:) = 0

      !BC for tang stress at r=a
      A23(Ndr,:) = 0

  !A24
  A24(:,:) = 0

      !BC for tang stress at r=a
      A24(Ndr,:) = MU2*ii*K*Dddt(0,:,0)

  !A25
  A25(:,:) = 0

      !BC for tang stress at r=a
      A25(Ndr,:) = MU2*Dddt(0,:,0)

  !A26
  A26(:,:) = 0

  !CONTINUITY

  !A31
  A31(:,:) = Dddr(:,:,0) + xDddr(:,:,1)

      !BC singularity condition on P
      A31(0,:) = 0

  !A32
  A32(:,:) = ii*K*xDddr(:,:,0)

      !BC singularity condition on P
      A32(0,:) = 0

  !A33
  A33(:,:) = 0

      !bc on singularity condition on P
      A33(0,:) = Dddr(0,:,1)

  !A34
  A34(:,:) = 0

  !A35
  A35(:,:) = 0

  !A36
  A36(:,:) = 0


  !FLUID 2 EQUATIONS

  !R-MOMENTUM

  !A41
  A41(:,:) = 0

    !bc for continuity at interface
    A41(0,:) = Dddr(Ndr,:,0)

  !A42
  A42(:,:) = 0

  !A43
  A43(:,:) = 0

  !A44
  A44(:,:) = MU2*xDddt(:,:,1) + MU2*x2Dddt(:,:,2) - MU2*(K**2)*x2Dddt(:,:,0) &
    & - MU2*Dddt(:,:,0) - RHO2*ii*K*(KAP/(MU2*4))*x4Dddt(:,:,0) - RHO2*ii*K*AA2*x2Dddt(:,:,0)

    !bc for no penetration at r=H
    A44(Ndr,:) = Dddt(Ndr,:,0)

    !bc for continuity at interface
    A44(0,:) = (-1)*Dddt(0,:,0)

  !A45
  A45(:,:) = 0

  !A46
  A46(:,:) = (-1)*x2Dddt(:,:,1)

    !bc for no penetration at r=H
    A46(Ndr,:) = 0

    !bc for continuity at interface
    A46(0,:) = 0

  !X-MOMENTUM

  !A51
  A51(:,:) = 0

  !A52
  A52(:,:) = 0

    !bc for no slip at interface
    A52(0,:) = Dddr(Ndr,:,0)

  !A53
  A53(:,:) = 0

  !A54
  A54(:,:) = -RHO2*(KAP/(2*MU2))*x2Dddt(:,:,1)

    !bc no slip at r=H
    A54(Ndr,:) = 0

    !bc no slip at interface
    A54(0,:) = 0

  !A55
  A55(:,:) = MU2*xDddt(:,:,2) + MU2*Dddt(:,:,1) - (K**2)*MU2*xDddt(:,:,0) - RHO2*ii*K*(KAP/(MU2*4))*x3Dddt(:,:,0) &
    & - RHO2*ii*K*AA2*xDddt(:,:,0)

    !bc no slip at r=H
    A55(Ndr,:) = Dddt(Ndr,:,0)

    !bc no slip at interface
    A55(0,:) = (-1)*Dddt(0,:,0)

  !A56
  A56(:,:) = -1*ii*K*xDddt(:,:,0)

    !bc for no slip at r=H
    A56(Ndr,:) = 0

    !bc no slip at interface
    A56(0,:) = 0

  !CONTINUITY

  !A61
  A61(:,:) = 0

  !A62
  A62(:,:) = 0

  !A63
  A63(:,:) = 0

  !A64
  A64(:,:) = Dddt(:,:,0) + xDddt(:,:,1)

  !A65
  A65(:,:) = ii*K*xDddt(:,:,0)

  !A66
  A66(:,:) = 0



  !put together the LHS

  !first row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1) = A11(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+ 1*(Ndr+1)) = A12(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+ 2*(Ndr+1)) = A13(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+ 3*(Ndr+1)) = A14(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+ 4*(Ndr+1)) = A15(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+ 5*(Ndr+1)) = A16(i,j)
  End Do
  End Do

  !second row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1) = A21(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1+ 1*(Ndr+1)) = A22(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1+ 2*(Ndr+1)) = A23(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1+ 3*(Ndr+1)) = A24(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1+ 4*(Ndr+1)) = A25(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1+ 5*(Ndr+1)) = A26(i,j)
  End Do
  End Do


  !third row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 2*(Ndr+1),j+1) = A31(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 2*(Ndr+1),j+1+ 1*(Ndr+1)) = A32(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 2*(Ndr+1),j+1+ 2*(Ndr+1)) = A33(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 2*(Ndr+1),j+1+ 3*(Ndr+1)) = A34(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 2*(Ndr+1),j+1+ 4*(Ndr+1)) = A35(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 2*(Ndr+1),j+1+ 5*(Ndr+1)) = A36(i,j)
  End Do
  End Do

  !fourth row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 3*(Ndr+1),j+1) = A41(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 3*(Ndr+1),j+1+ 1*(Ndr+1)) = A42(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 3*(Ndr+1),j+1+ 2*(Ndr+1)) = A43(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 3*(Ndr+1),j+1+ 3*(Ndr+1)) = A44(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 3*(Ndr+1),j+1+ 4*(Ndr+1)) = A45(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 3*(Ndr+1),j+1+ 5*(Ndr+1)) = A46(i,j)
  End Do
  End Do

  !fifth row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 4*(Ndr+1),j+1) = A51(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 4*(Ndr+1),j+1+ 1*(Ndr+1)) = A52(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 4*(Ndr+1),j+1+ 2*(Ndr+1)) = A53(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 4*(Ndr+1),j+1+ 3*(Ndr+1)) = A54(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 4*(Ndr+1),j+1+ 4*(Ndr+1)) = A55(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 4*(Ndr+1),j+1+ 5*(Ndr+1)) = A56(i,j)
  End Do
  End Do

  !sixth row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 5*(Ndr+1),j+1) = A61(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 5*(Ndr+1),j+1+ 1*(Ndr+1)) = A62(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 5*(Ndr+1),j+1+ 2*(Ndr+1)) = A63(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 5*(Ndr+1),j+1+ 3*(Ndr+1)) = A64(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 5*(Ndr+1),j+1+ 4*(Ndr+1)) = A65(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+ 5*(Ndr+1),j+1+ 5*(Ndr+1)) = A66(i,j)
  End Do
  End Do


  call matrixAPrint(A,A11,A12,A13,A21,A22,A23,A31,A32,A33)

  END SUBROUTINE buildA

SUBROUTINE buildB(B,xDddr,x2Dddr,xDddt,x2Dddt,x3Dddr,x3Dddt,x4Dddr,x4Dddt)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: B11, B12, B13, B14, B15, B16, B21, B22, B23, B24, B25, B26, B31, B32, B33, B34, B35, B36
  COMPLEX, allocatable, dimension(:,:)  :: B41, B42, B43, B44, B45, B46, B51, B52, B53, B54, B55, B56, B61, B62, B63, B64, B65, B66

  Real, allocatable, dimension(:,:,:)   :: xDddr,x2Dddr,xDddt,x2Dddt,x3Dddr,x3Dddt,x4Dddr,x4Dddt
  COMPLEX *16                           :: B(LDB,NMAX)

  allocate(B11(0:Ndr,0:Ndr))
  allocate(B12(0:Ndr,0:Ndr))
  allocate(B13(0:Ndr,0:Ndr))
  allocate(B14(0:Ndr,0:Ndr))
  allocate(B15(0:Ndr,0:Ndr))
  allocate(B16(0:Ndr,0:Ndr))

  allocate(B21(0:Ndr,0:Ndr))
  allocate(B22(0:Ndr,0:Ndr))
  allocate(B23(0:Ndr,0:Ndr))
  allocate(B24(0:Ndr,0:Ndr))
  allocate(B25(0:Ndr,0:Ndr))
  allocate(B26(0:Ndr,0:Ndr))

  allocate(B31(0:Ndr,0:Ndr))
  allocate(B32(0:Ndr,0:Ndr))
  allocate(B33(0:Ndr,0:Ndr))
  allocate(B34(0:Ndr,0:Ndr))
  allocate(B35(0:Ndr,0:Ndr))
  allocate(B36(0:Ndr,0:Ndr))

  allocate(B41(0:Ndr,0:Ndr))
  allocate(B42(0:Ndr,0:Ndr))
  allocate(B43(0:Ndr,0:Ndr))
  allocate(B44(0:Ndr,0:Ndr))
  allocate(B45(0:Ndr,0:Ndr))
  allocate(B46(0:Ndr,0:Ndr))

  allocate(B51(0:Ndr,0:Ndr))
  allocate(B52(0:Ndr,0:Ndr))
  allocate(B53(0:Ndr,0:Ndr))
  allocate(B54(0:Ndr,0:Ndr))
  allocate(B55(0:Ndr,0:Ndr))
  allocate(B56(0:Ndr,0:Ndr))

  allocate(B61(0:Ndr,0:Ndr))
  allocate(B62(0:Ndr,0:Ndr))
  allocate(B63(0:Ndr,0:Ndr))
  allocate(B64(0:Ndr,0:Ndr))
  allocate(B65(0:Ndr,0:Ndr))
  allocate(B66(0:Ndr,0:Ndr))

  !FLUID 1 EQUATIONS

  !R-MOMENTUM

  !B11
  B11(:,:) = RHO1*x2Dddr(:,:,0)

      !singularity bc on Ur
      B11(0,:) = 0

      !normal stress at r=a
      B11(Ndr,:) = -2.*MU1*Dddr(Ndr,:,1)

  !B12
  B12(:,:) = 0

  !B13
  B13(:,:) = 0

      !bc for normal stress at r=a
      B13(Ndr,:) = 1.*Dddr(Ndr,:,0)

  !B14
  B14(:,:) = 0

      !bc for normal stress at r=a
      B14(Ndr,:) = 2*MU2*Dddt(0,:,1)

  !B15
  B15(:,:) = 0

  !B16
  B16(:,:) = 0

      !bc for normal stress at r=a
      B16(Ndr,:) = (-1)*Dddt(0,:,0)

  !X-MOMENTUM

  !B21
  B21(:,:) = 0

  !B22
  B22(:,:) = RHO1*xDddr(:,:,0)

      !singularity BC on Ux
      B22(0,:) = 0

      !tang stress at r=a
      B22(Ndr,:) = 0

  !B23
  B23(:,:) = 0

  !B24
  B24(:,:) = 0

  !B25
  B25(:,:) = 0

  !B26
  B26(:,:) = 0

  !CONTINUITY

  !B31
  B31(:,:) = 0

  !B32
  B32(:,:) = 0

  !B33
  B33(:,:) = 0

  !B34
  B34(:,:) = 0

  !B35
  B35(:,:) = 0

  !B36
  B36(:,:) = 0

  !FLUID 2 EQUATIONS

  !R-MOMENTUM

  !B41
  B41(:,:) = 0

  !B42
  B42(:,:) = 0

  !B43
  B43(:,:) = 0

  !B44
  B44(:,:) = RHO2*x2Dddt(:,:,0)

    !bc for no penetration at r=H
    B44(Ndr,:) = 0

    !bc for continuity at interface
    B44(0,:) = 0

  !B45
  B45(:,:) = 0

  !B46
  B46(:,:) = 0

  !X-MOMENTUM

  !B51
  B51(:,:) = 0

  !B52
  B52(:,:) = 0

  !B53
  B53(:,:) = 0

  !B54
  B54(:,:) = 0

  !B55
  B55(:,:) = RHO2*xDddt(:,:,0)

    !bc for no slip at r=H
    B55(Ndr,:) = 0

    !bc no slip at interface
    B55(0,:) = 0

  !B56
  B56(:,:) = 0

  !CONTINUITY

  !B61
  B61(:,:) = 0

  !B62
  B62(:,:) = 0

  !B63
  B63(:,:) = 0

  !B64
  B64(:,:) = 0

  !B65
  B65(:,:) = 0

  !B66
  B66(:,:) = 0


  !put together the RHS

  !first row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1) = B11(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+ 1*(Ndr+1)) = B12(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+ 2*(Ndr+1)) = B13(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+ 3*(Ndr+1)) = B14(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+ 4*(Ndr+1)) = B15(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+ 5*(Ndr+1)) = B16(i,j)
  End Do
  End Do

  !second row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1) = B21(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1+ 1*(Ndr+1)) = B22(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1+ 2*(Ndr+1)) = B23(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1+ 3*(Ndr+1)) = B24(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1+ 4*(Ndr+1)) = B25(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1+ 5*(Ndr+1)) = B26(i,j)
  End Do
  End Do


  !third row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 2*(Ndr+1),j+1) = B31(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 2*(Ndr+1),j+1+ 1*(Ndr+1)) = B32(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 2*(Ndr+1),j+1+ 2*(Ndr+1)) = B33(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 2*(Ndr+1),j+1+ 3*(Ndr+1)) = B34(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 2*(Ndr+1),j+1+ 4*(Ndr+1)) = B35(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 2*(Ndr+1),j+1+ 5*(Ndr+1)) = B36(i,j)
  End Do
  End Do

  !fourth row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 3*(Ndr+1),j+1) = B41(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 3*(Ndr+1),j+1+ 1*(Ndr+1)) = B42(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 3*(Ndr+1),j+1+ 2*(Ndr+1)) = B43(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 3*(Ndr+1),j+1+ 3*(Ndr+1)) = B44(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 3*(Ndr+1),j+1+ 4*(Ndr+1)) = B45(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 3*(Ndr+1),j+1+ 5*(Ndr+1)) = B46(i,j)
  End Do
  End Do

  !fifth row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 4*(Ndr+1),j+1) = B51(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 4*(Ndr+1),j+1+ 1*(Ndr+1)) = B52(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 4*(Ndr+1),j+1+ 2*(Ndr+1)) = B53(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 4*(Ndr+1),j+1+ 3*(Ndr+1)) = B54(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 4*(Ndr+1),j+1+ 4*(Ndr+1)) = B55(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 4*(Ndr+1),j+1+ 5*(Ndr+1)) = B56(i,j)
  End Do
  End Do

  !sixth row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 5*(Ndr+1),j+1) = B61(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 5*(Ndr+1),j+1+ 1*(Ndr+1)) = B62(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 5*(Ndr+1),j+1+ 2*(Ndr+1)) = B63(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 5*(Ndr+1),j+1+ 3*(Ndr+1)) = B64(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 5*(Ndr+1),j+1+ 4*(Ndr+1)) = B65(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+ 5*(Ndr+1),j+1+ 5*(Ndr+1)) = B66(i,j)
  End Do
  End Do

  call matrixBPrint(B,B11,B12,B13,B21,B22,B23,B31,B32,B33)

  END SUBROUTINE buildB

SUBROUTINE matrixAPrint(A,A11,A12,A13,A21,A22,A23,A31,A32,A33)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: A11, A12, A13, A21, A22, A23, A31, A32, A33
  COMPLEX *16      A(LDA,NMAX)

  !print each individual A matrix separately

  open(unit=1, file='M/graphA11.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(1, '(1600F9.2)')( Real(A11(i,j)) ,j=0,Ndr)
  end do
  close(1)

  open(unit=2, file='M/graphA12.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(2, '(1600F9.2)')( Real(A12(i,j)) ,j=0,Ndr)
  end do
  close(2)

  open(unit=3, file='M/graphA13.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(3, '(1600F9.2)')( Real(A13(i,j)) ,j=0,Ndr)
  end do
  close(3)

  open(unit=4, file='M/graphA21.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(4, '(1600F9.2)')( Real(A21(i,j)) ,j=0,Ndr)
  end do
  close(4)

  open(unit=111, file='M/graphA22.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(111, '(1600F9.2)')( Real(A22(i,j)) ,j=0,Ndr)
  end do
  close(111)

  open(unit=112, file='M/graphA23.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(112, '(1600F9.2)')( aImag(A23(i,j)) ,j=0,Ndr)
  end do
  close(112)

  open(unit=113, file='M/graphA31.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(113, '(1600F9.2)')( Real(A31(i,j)) ,j=0,Ndr)
  end do
  close(113)

  open(unit=114, file='M/graphA32.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(114, '(1600F9.2)')( aIMag(A32(i,j)) ,j=0,Ndr)
  end do
  close(114)

  open(unit=115, file='M/graphA33.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(115, '(1600F9.2)')( Real(A33(i,j)) ,j=0,Ndr)
  end do
  close(115)

  !write the output LHS matrix

  open(unit=11, file='M/graphLHS.txt', ACTION="write", STATUS="replace")
  do i=0,(3*Ndr)+2
    write(11, '(1600F7.2)')( Real(A(i+1,j+1)) ,j=0,(3*Ndr)+2)
  end do
  close(11)

  open(unit=24, file='M/graphLHSImag.txt', ACTION="write", STATUS="replace")
  do i=0,(3*Ndr)+2
    write(24, '(1600F7.2)')( aImag(A(i+1,j+1)) ,j=0,(3*Ndr)+2)
  end do
  close(24)

  END SUBROUTINE matrixAPrint

SUBROUTINE matrixBPrint(B,B11,B12,B13,B21,B22,B23,B31,B32,B33)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: B11, B12, B13, B21, B22, B23, B31, B32, B33
  COMPLEX *16      B(LDA,NMAX)

  open(unit=116, file='M/graphB11.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(116, '(1600F9.2)')( Real(B11(i,j)) ,j=0,Ndr)
  end do
  close(116)

  open(unit=117, file='M/graphB12.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(117, '(1600F9.2)')( Real(B12(i,j)) ,j=0,Ndr)
  end do
  close(117)

  open(unit=118, file='M/graphB13.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(118, '(1600F9.2)')( Real(B13(i,j)) ,j=0,Ndr)
  end do
  close(118)

  open(unit=119, file='M/graphB21.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(119, '(1600F9.2)')( Real(B21(i,j)) ,j=0,Ndr)
  end do
  close(119)

  open(unit=120, file='M/graphB22.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(120, '(1600F9.2)')( Real(B22(i,j)) ,j=0,Ndr)
  end do
  close(120)

  open(unit=121, file='M/graphB23.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(121, '(1600F9.2)')( Real(B23(i,j)) ,j=0,Ndr)
  end do
  close(121)

  open(unit=122, file='M/graphB31.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(122, '(1600F9.2)')( Real(B31(i,j)) ,j=0,Ndr)
  end do
  close(122)

  open(unit=122, file='M/graphB32.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(122, '(1600F9.2)')( Real(B32(i,j)) ,j=0,Ndr)
  end do
  close(122)

  open(unit=123, file='M/graphB33.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(123, '(1600F9.2)')( Real(B33(i,j)) ,j=0,Ndr)
  end do
  close(123)

  !write the output RHS matrix

  open(unit=12, file='M/graphRHS.txt', ACTION="write", STATUS="replace")
  do i=0,(3*Ndr)+2
    write(12, '(1600F7.2)')( Real(B(i+1,j+1)) ,j=0,(3*Ndr)+2)
  end do
  close(12)

  open(unit=142, file='M/graphImagRHS.txt', ACTION="write", STATUS="replace")
  do i=0,(3*Ndr)+2
    write(142, '(1600F7.2)')( aImag(B(i+1,j+1)) ,j=0,(3*Ndr)+2)
  end do
  close(142)

  END SUBROUTINE matrixBPrint

END MODULE ev

