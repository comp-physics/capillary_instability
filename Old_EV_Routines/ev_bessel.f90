!bessel eq solver!

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
  !1D
  N=Ndr+1

  iterMAx=1
  aa=RMAX

  !fix the sign error for odd derivatives
  Dddr(:,:,1) = -1*Dddr(:,:,1)
  Dddr(:,:,3) = -1*Dddr(:,:,3)

  !make the xDddr matrix that premultiplies Dddr(:) by XXdr
  Do i=0,3
  xDddr(:,:,i) = Matmul(XXdr(:,:),Dddr(:,:,i))
  End Do

  !write the xDddr matrix 
  open(unit=13, file='graphxDddr.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(13, '(1600F9.2)')( Real(xDddr(i,j,2)) ,j=0,Ndr)
  end do
  close(13)

  !write the xDddr matrix 
  open(unit=133, file='graphxxdr.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(133, '(1600F9.2)')( Real(xxdr(i,j)) ,j=0,Ndr)
  end do
  close(133)

  open(unit=135, file='graphDddr1.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(135, '(1600F9.2)')( Real(Dddr(i,j,1)) ,j=0,Ndr)
  end do
  close(135)

 open(unit=136, file='graphDddr2.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(136, '(1600F9.2)')( Real(Dddr(i,j,2)) ,j=0,Ndr)
  end do
  close(136)

  !looping over K
  Do iter=1,itermax

  k=1.0

  Do i=0,Ndr
  Do j=0,Ndr
  A(i+1,j+1) = xDddr(i,j,2) + Dddr(i,j,1)
  End Do
  End Do

  Do i=0,Ndr
  A(1,i+1) = Dddr(0,i,1)
  End Do

  Do i=0,Ndr
  A(Ndr+1,i+1) = Dddr(Ndr,i,0)
  End Do
  
  Do i=0,Ndr
  Do j=0,Ndr
  B(i+1,j+1) = xDddr(i,j,0)
  End Do
  End Do

  Do i=0,Ndr
  B(Ndr+1,i+1) = 0
  End Do

  Do i=0,Ndr
  B(1,i+1) = 0
  End Do

  call matrixAprint(A)
  call matrixBprint(B)
 
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
                  write(1,"(3F20.10)")  ALPHA(J)/BETA(J)
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

SUBROUTINE matrixAPrint(A)

  integer :: i,j
  COMPLEX *16      A(LDA,NMAX)

  open(unit=11, file='M/graphLHS.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(11, '(1600F7.2)')( Real(A(i+1,j+1)) ,j=0,Ndr)
  end do
  close(11)

  open(unit=24, file='M/graphLHSImag.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(24, '(1600F7.2)')( aImag(A(i+1,j+1)) ,j=0,Ndr)
  end do
  close(24)

  END SUBROUTINE matrixAPrint

SUBROUTINE matrixBPrint(B)

  integer :: i,j
  COMPLEX *16      B(LDA,NMAX)

  !write the output RHS matrix

  open(unit=12, file='M/graphRHS.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(12, '(1600F7.2)')( Real(B(i+1,j+1)) ,j=0,Ndr)
  end do
  close(12)

  open(unit=142, file='M/graphImagRHS.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(142, '(1600F7.2)')( aImag(B(i+1,j+1)) ,j=0,Ndr)
  end do
  close(142)

  END SUBROUTINE matrixBPrint

END MODULE ev



