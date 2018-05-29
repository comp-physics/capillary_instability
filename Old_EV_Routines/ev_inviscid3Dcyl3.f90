!3d inviscid solver!

MODULE ev

  USE prms
  USE data
  USE cheb
  
  IMPLICIT none 

INTEGER          NOUT
PARAMETER        (NOUT=6)
INTEGER          NB, NMAX, lt
PARAMETER        (NB=100,NMAX=100)
INTEGER          LDA, LDB, LDVR, LWORK
PARAMETER        (LDA=NMAX,LDB=NMAX,LDVR=NMAX,LWORK=NMAX+NMAX*NB)
COMPLEX *16      ii
PARAMETER        (ii=(0,1.0))

!FLOW PARAMETERS
    Real             GAM,RHO,MU,aa
    PARAMETER        (GAM=1.0,RHO=1.0,MU=1.0)

    REAL          K
    integer       iter,iterMax

CONTAINS

SUBROUTINE evsolve

  INTEGER          INFO,LWKOPT, N, U, W
  COMPLEX *16      A(LDA,NMAX), ALPHA(NMAX), B(LDB,NMAX),BETA(NMAX), DUMMY(1,1), VR(LDVR,NMAX),WORK(LWORK)
  DOUBLE PRECISION RWORK(8*NMAX), SMALL, DLAMCH

  integer :: i,j
  character(30)          :: fn

  Real, allocatable, dimension(:,:,:)  :: xDddr

  allocate(xDddr(0:Ndr,0:Ndr,0:3))

  lt=1
  N=3*(Ndr+1)

  iterMAx=100
  aa=RMAX

  !fix the sign error for odd derivatives
  Dddr(:,:,1) = -1*Dddr(:,:,1)
  Dddr(:,:,3) = -1*Dddr(:,:,3)

  !Dddr(:,:,2) = -1*Dddr(:,:,2)
  !Dddr(:,:,4) = -1*Dddr(:,:,4)  

  !make the xDddr matrix that premultiplies Dddr(:) by XXdr
  Do i=0,3
  xDddr(:,:,i) = Matmul(XXdr(:,:),Dddr(:,:,i))
  End Do

  !write the xDddr matrix 
  open(unit=13, file='graphxDddr.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(13, '(1600F9.2)')( Real(xDddr(i,j,1)) ,j=0,Ndr)
  end do
  close(13)

  !write the xDddr matrix 
  open(unit=133, file='graphxxdr.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(133, '(1600F9.2)')( Real(xxdr(i,j)) ,j=0,Ndr)
  end do
  close(133)

  !looping over K
  Do iter=1,itermax

  k=Real(iter)/Real(itermax)
  print*, k

  call buildA(A,xDddr)
  call buildB(B,xDddr)


  !lapack gevp solver
  print*, N, "N"
  IF (N.LE.NMAX) THEN
         CALL ZGGEV('No left vectors','Vectors (right)',N,A,LDA,B,LDB,&
     &              ALPHA,BETA,DUMMY,1,VR,LDVR,WORK,LWORK,RWORK,INFO)

            !write eigenvalues to files in directory D 
            Do J=1,N
               write(fn,"('D/eval',I9.9)")iter
                open(1,file=fn)
                If (Real(BETA(j))>0.00001) Then
                  write(1,"(3F35.25)")  ALPHA(J)/BETA(J)
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
      !STOP

      99999 FORMAT (1X,A,I4)
      99998 FORMAT (1X,A,I2,2A,/1X,2(A,I2,A,'(',1P,E11.4,',',1P,E11.4,')'))
      99997 FORMAT (1X,A,I2,A,'(',1P,E11.4,',',1P,E11.4,')')
      99996 FORMAT (1X,A,I2,A,/3(1X,'(',1P,E11.4,',',1P,E11.4,')',:))
      99995 FORMAT (1X,A,I5,/1X,A,I5)

  End Do

  END SUBROUTINE evsolve

SUBROUTINE buildA(A,xDddr)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: A11, A12, A13, A21, A22, A23, A31, A32, A33
  Real, allocatable, dimension(:,:,:)   :: xDddr
  COMPLEX *16                           :: A(LDA,NMAX)

  allocate(A11(0:Ndr,0:Ndr))  
  allocate(A12(0:Ndr,0:Ndr))  
  allocate(A13(0:Ndr,0:Ndr))  

  allocate(A21(0:Ndr,0:Ndr))
  allocate(A22(0:Ndr,0:Ndr))  
  allocate(A23(0:Ndr,0:Ndr))  

  allocate(A31(0:Ndr,0:Ndr))  
  allocate(A32(0:Ndr,0:Ndr))
  allocate(A33(0:Ndr,0:Ndr))


  !construct A11 matrix

  A11(:,:) = 0

      !BC for stress condition at r=aa

      A11(Ndr,:) = Dddr(Ndr,:,0)*GAM*(1./(aa**2)-K**2.)

      !BC for singularity condition on Ur at r=0

      A11(0,:) = Dddr(0,:,0)

  !construct A12 matrix

  A12(:,:) = 0

      !BC for stress at r=aa

      A12(0,:) = 0

      !BC for singularity on Ur

      A12(0,:) = 0

  !construct A13 matrix

  A13(:,:) = (-1/RHO)*Dddr(:,:,1)

      !BC for stress at r=aa

      A13(Ndr,:) = 0

      !BC singularity condition on Ur

      A13(0,:) = 0

  !construct A21 matrix

  A21(:,:) = 0

      !BC for singularity on Ux

      A21(0,:) = 0

  !construct A22 matrix

  A22(:,:) = 0

      !BC singularity condition on Ux

      A22(0,:) = Dddr(0,:,1)

  !construct A23 matrix

  A23(:,:) = (-1/RHO)*ii*K*Dddr(:,:,0)

      !BC singularity condition on Ux

      A23(0,:) = 0

  !construct A31 matrix
  
  A31(:,:) = Dddr(:,:,0) + xDddr(:,:,1)

      !BC singularity condition on P

      A31(0,:) = 0

  !construct A32 matrix

  A32(:,:) = ii*K*xDddr(:,:,0)

      !BC singularity condition on P

      A32(0,:) = 0

  !construct A33 matrix

  A33(:,:) = 0

      !bc on singularity condition on P

      A33(0,:) = Dddr(0,:,1)


  !put together the LHS

  !first row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1) = A11(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+Ndr+1) = A12(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+Ndr+1+Ndr+1) = A13(i,j)
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
    A(i+1+Ndr+1,j+1+Ndr+1) = A22(i,j)
  End Do
  End Do  

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1+Ndr+1+Ndr+1) = A23(i,j)
  End Do
  End Do

  !third row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1+Ndr+1,j+1) = A31(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1+Ndr+1,j+1+Ndr+1) = A32(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1+Ndr+1,j+1+Ndr+1+Ndr+1) = A33(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Print*, 'A11', A11(Ndr,i)
  end do


  call matrixAPrint(A,A11,A12,A13,A21,A22,A23,A31,A32,A33)

  END SUBROUTINE buildA

SUBROUTINE buildB(B,xDddr)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: B11, B12, B13, B21, B22, B23, B31, B32, B33
  Real, allocatable, dimension(:,:,:)   :: xDddr
  COMPLEX *16                           :: B(LDB,NMAX)

  allocate(B11(0:Ndr,0:Ndr))  
  allocate(B12(0:Ndr,0:Ndr))  
  allocate(B13(0:Ndr,0:Ndr))  

  allocate(B21(0:Ndr,0:Ndr))
  allocate(B22(0:Ndr,0:Ndr))  
  allocate(B23(0:Ndr,0:Ndr))  

  allocate(B31(0:Ndr,0:Ndr))  
  allocate(B32(0:Ndr,0:Ndr))
  allocate(B33(0:Ndr,0:Ndr))

  !construct B11 matrix

  B11(:,:) = Dddr(:,:,0)

      !stress bc on r = a

      B11(Ndr,:) = 0

      !singularity bc on Ur

      B11(0,:) = 0

  !construct B12 matrix

  B12(:,:) = 0

  !construct B13 matrix

  B13(:,:) = 0

      !bc for stress at r=a

      B13(Ndr,:) = -1.*Dddr(Ndr,:,0)

  !construct B21 matrix

  B21(:,:) = 0

  !construct B22 matrix

  B22(:,:) = Dddr(:,:,0)

      !singularity BC on Ux

      B22(0,:) = 0

  !construct B23 matrix

  B23(:,:) = 0

  !construct B31 matrix
  
  B31(:,:) = 0

  !construct B32 matrix

  B32(:,:) = 0

  !construct B33 matrix

  B33(:,:) = 0


  !put together the RHS

  !first row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1) = B11(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+Ndr+1) = B12(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+Ndr+1+Ndr+1) = B13(i,j)
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
    B(i+1+Ndr+1,j+1+Ndr+1) = B22(i,j)
  End Do
  End Do  

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1+Ndr+1+Ndr+1) = B23(i,j)
  End Do
  End Do

  !third row of blocks
  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1+Ndr+1,j+1) = B31(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1+Ndr+1,j+1+Ndr+1) = B32(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1+Ndr+1,j+1+Ndr+1+Ndr+1) = B33(i,j)
  End Do
  End Do

  !do i=0,Ndr
  !print*, 'B13', B13(Ndr,i)
  !end do

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



