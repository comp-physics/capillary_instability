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


CONTAINS

  SUBROUTINE evsolve

      INTEGER          INFO,LWKOPT, N, U, W
      COMPLEX *16      A(LDA,NMAX), ALPHA(NMAX), B(LDB,NMAX),BETA(NMAX), DUMMY(1,1), VR(LDVR,NMAX),WORK(LWORK)
      DOUBLE PRECISION RWORK(8*NMAX), SMALL, DLAMCH

	  integer :: i,j
    character(30)          :: fn

call buildA(A)
call buildB(B)

!lapack gevp solver

lt=1
N=2*(Ndr+1)

IF (N.LE.NMAX) THEN
         CALL ZGGEV('No left vectors','Vectors (right)',N,A,LDA,B,LDB,&
     &              ALPHA,BETA,DUMMY,1,VR,LDVR,WORK,LWORK,RWORK,INFO)

!write eigenvalues to files in directory D 
            Do J=1,N
               write(fn,"('D/x.',I9.9)")lt
                open(1,file=fn)
                If (Real(BETA(j))>0.00001) Then
                  write(1,"(3F20.10)")  ALPHA(J)/BETA(J)
                End if

            End Do
               close(1)

!write eigenvectors to files in directory V
            Do lt=1,N
              Do I=1,N
               write(fn,"('V/x.',I9.9)")lt
                open(2,file=fn)
               write(2,"(3F20.10)")  VR(I,lt)
              End Do
            End Do
               close(2)

         IF (INFO.GT.0) THEN
!            WRITE (NOUT,*)
            WRITE (NOUT,99999) 'Failure in ZGGEV. INFO =', INFO
         ELSE
            SMALL = DLAMCH('Sfmin')


            DO 20 J = 1,N
!               WRITE (NOUT,*)
               IF ((ABS(ALPHA(J)))*SMALL.GE.ABS(BETA(J))) THEN
                  WRITE (NOUT,99998) 'Eigenvalue(', J, ') is numerically infinite or undetermined'
               ELSE
               !print eigenvalues to the shell

                  WRITE (NOUT,99997) 'Eigenvalue(', J, ') = ',ALPHA(J)/BETA(J)
               END IF
               !print eigenvectors to the shell 

!               WRITE (NOUT,*)
!              WRITE (NOUT,99996) 'Eigenvector(', J, ')',&
!     &           (VR(I,J),I=1,N)
   20       CONTINUE

            LWKOPT = WORK(1)
            IF (LWORK.LT.LWKOPT) THEN
               WRITE (NOUT,*)
               WRITE (NOUT,99995) 'Optimum workspace required = ',&
     &           LWKOPT, 'Workspace provided         = ', LWORK
            END IF
         END IF
      ELSE
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'NMAX too small'
      END IF
      STOP

99999 FORMAT (1X,A,I4)
99998 FORMAT (1X,A,I2,2A,/1X,2(A,I2,A,'(',1P,E11.4,',',1P,E11.4,')'))
99997 FORMAT (1X,A,I2,A,'(',1P,E11.4,',',1P,E11.4,')')
99996 FORMAT (1X,A,I2,A,/3(1X,'(',1P,E11.4,',',1P,E11.4,')',:))
99995 FORMAT (1X,A,I5,/1X,A,I5)

  END SUBROUTINE evsolve


  SUBROUTINE buildA(A)

    integer :: i,j
    real, allocatable, dimension(:,:)  :: TR, BL
    COMPLEX *16      A(LDA,NMAX)

  allocate(TR(0:Ndr,0:Ndr))  
  allocate(BL(0:Ndr,0:Ndr))


!construct the TR matrix
TR(:,:) = Dddr(:,:,1)

!construct the BL matrix
BL(:,:) = Dddr(:,:,1)

BL(0,:) = Dddr(0,:,0)
BL(Ndr,:) = Dddr(Ndr,:,0)

!LHS

!insert top right part

Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+Ndr+1) = TR(i,j)
  End Do
End Do

!insert bottom left part

Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1) = BL(i,j)
  End Do
End Do




!write the first derivative matrix

open(unit=1, file='graphTR.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(1, '(1600F7.2)')( TR(i,j) ,j=0,Ndr)
  end do
close(1)

!write the second derivative matrix

open(unit=2, file='graphBL.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(2, '(1600F7.2)')( BL(i,j) ,j=0,Ndr)
  end do
close(2)

!write the output LHS matrix

open(unit=11, file='graphLHS.txt', ACTION="write", STATUS="replace")
  do i=0,(2*Ndr)+1
    write(11, '(1600F7.2)')( Real(A(i+1,j+1)) ,j=0,(2*Ndr)+1)
  end do
close(11)


  END SUBROUTINE buildA


  SUBROUTINE buildB(B)

    integer :: i,j
    real, allocatable, dimension(:,:)  :: RHS
    COMPLEX *16      A(LDA,NMAX), B(LDB,NMAX)

    allocate(RHS(0:2*Ndr,0:2*Ndr))


  !construct the RHS matrix

Do i=0,2*Ndr
  Do j=0,2*Ndr
    If (i == j) Then
      RHS(i,j) = 1
    Else
      RHS(i,j) = 0
    End If
  End Do
End Do

RHS(0,0) = 0
RHS(Ndr,Ndr) = 0
RHS(Ndr+1,Ndr+1) = 0
RHS(2*Ndr+1,2*Ndr+1) = 0

!RHS

Do i=0,2*Ndr
  Do j=0,2*Ndr
    B(i+1,j+1) = RHS(i,j)
  End Do
End Do

!write the RHS matrix 

  open(unit=10, file='graphRHS.txt', ACTION="write", STATUS="replace")
    do i=0,2*Ndr+1
      write(10, '(1600F8.2)')( Real(B(i+1,j+1)) ,j=0,2*Ndr+1)
    end do
  close(10)

  END SUBROUTINE buildB

 END MODULE ev

