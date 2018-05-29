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
    COMPLEX *16 ii

!FLOW PARAMETERS
    Real             GAM,K,RHO,MU

CONTAINS

SUBROUTINE evsolve

  INTEGER          INFO,LWKOPT, N, U, W
  COMPLEX *16      A(LDA,NMAX), ALPHA(NMAX), B(LDB,NMAX),BETA(NMAX), DUMMY(1,1), VR(LDVR,NMAX),WORK(LWORK)
  DOUBLE PRECISION RWORK(8*NMAX), SMALL, DLAMCH

  integer :: i,j
  character(30)          :: fn

  ii=(0,1)
  K=1
  
  GAM=1
  RHO=1
  MU=1

  !fix the sign error for odd derivatives
  Dddr(:,:,1) = -1*Dddr(:,:,1)
  Dddr(:,:,3) = -1*Dddr(:,:,3)

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
      STOP

      99999 FORMAT (1X,A,I4)
      99998 FORMAT (1X,A,I2,2A,/1X,2(A,I2,A,'(',1P,E11.4,',',1P,E11.4,')'))
      99997 FORMAT (1X,A,I2,A,'(',1P,E11.4,',',1P,E11.4,')')
      99996 FORMAT (1X,A,I2,A,/3(1X,'(',1P,E11.4,',',1P,E11.4,')',:))
      99995 FORMAT (1X,A,I5,/1X,A,I5)

  END SUBROUTINE evsolve

SUBROUTINE buildA(A)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: TL, TR, BL, BR
  COMPLEX *16      A(LDA,NMAX)

  allocate(TL(0:Ndr,0:Ndr))  
  allocate(TR(0:Ndr,0:Ndr))  
  allocate(BL(0:Ndr,0:Ndr))  
  allocate(BR(0:Ndr,0:Ndr))

  !construct the TL matrix

  TL(:,:) = MU*ii*(K*Dddr(:,:,3)-(K**3)*Dddr(:,:,1))
  
    !add on boundary condition tan stress at +/- a

    TL(0,:) = Dddr(0,:,2) + Dddr(0,:,0)
    TL(Ndr,:) = Dddr(Ndr,:,2) + Dddr(Ndr,:,0)

  !construct the TR matrix

  TR(:,:) = ii*K*Dddr(:,:,0)

    !add on boundary condition tan stress at +/- a

    TR(0,:) = 0
    TR(Ndr,:) = 0

  !construct the BL matrix

  BL(:,:) = MU*(Dddr(:,:,2) - (K**2)*Dddr(:,:,0))

    !add on boundary condition normal stress at +/- a

    BL(0,:) = GAM*(K**2)*Dddr(0,:,0)
    BL(Ndr,:) = GAM*(K**2)*Dddr(Ndr,:,0)

  !construct the BR matrix

  BR(:,:) = Dddr(:,:,1)

    !add on boundary condition normal stress at +/- a

    BR(0,:) = 0
    BR(Ndr,:) = 0

  !put together the LHS

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1) = TL(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1,j+1+Ndr+1) = TR(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1) = BL(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    A(i+1+Ndr+1,j+1+Ndr+1) = BR(i,j)
  End Do
  End Do

  call matrixAPrint(BL,BR,TL,TR,A)

  END SUBROUTINE buildA

SUBROUTINE buildB(B)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: TL2, TR2, BL2, BR2
  COMPLEX *16      B(LDB,NMAX)

  allocate(TL2(0:Ndr,0:Ndr))  
  allocate(TR2(0:Ndr,0:Ndr))  
  allocate(BL2(0:Ndr,0:Ndr))  
  allocate(BR2(0:Ndr,0:Ndr))

  !construct the TL matrix

  TL2(:,:) = ii*RHO*K*Dddr(:,:,1)
  
    !add on boundary condition tan stress at +/- a

    TL2(0,:) = 0
    TL2(Ndr,:) = 0

  !construct the TR matrix

  TR2(:,:) = 0

    !add on boundary condition tan stress at +/- a

    TR2(0,:) = 0
    TR2(Ndr,:) = 0

  !construct the BL matrix

  BL2(:,:) = RHO*Dddr(:,:,0)

    !add on boundary condition normal stress at +/- a

    BL2(0,:) = MU*Dddr(0,:,0)
    BL2(Ndr,:) = MU*Dddr(Ndr,:,0)

  !construct the BR matrix

  BR2(:,:) = 0

    !add on boundary condition normal stress at +/- a

    BR2(0,:) = -1*Dddr(0,:,0)
    BR2(Ndr,:) = -1*Dddr(Ndr,:,0)

  !put together the RHS

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1) = TL2(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1,j+1+Ndr+1) = TR2(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1) = BL2(i,j)
  End Do
  End Do

  Do i=0,Ndr
  Do j=0,Ndr
    B(i+1+Ndr+1,j+1+Ndr+1) = BR2(i,j)
  End Do
  End Do

  call matrixBPrint(BL2,BR2,TL2,TR2,B)

  END SUBROUTINE buildB

SUBROUTINE matrixAPrint(BL,BR,TL,TR,A)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: TL, TR, BL, BR
  COMPLEX *16      A(LDA,NMAX)

  !write TL

  open(unit=1, file='graphTL.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(1, '(1600F9.2)')( Real(TL(i,j)) ,j=0,Ndr)
  end do
  close(1)

  !write TR

  open(unit=2, file='graphTR.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(2, '(1600F9.2)')( Real(TR(i,j)) ,j=0,Ndr)
  end do
  close(2)

  !write BL

  open(unit=3, file='graphBL.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(3, '(1600F9.2)')( Real(BL(i,j)) ,j=0,Ndr)
  end do
  close(3)

  !write BR

  open(unit=4, file='graphBR.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(4, '(1600F9.2)')( Real(BR(i,j)) ,j=0,Ndr)
  end do
  close(4)

  !write the output LHS matrix

  open(unit=11, file='graphLHS.txt', ACTION="write", STATUS="replace")
  do i=0,(2*Ndr)+1
    write(11, '(1600F7.2)')( Real(A(i+1,j+1)) ,j=0,(2*Ndr)+1)
  end do
  close(11)

  END SUBROUTINE matrixAPrint

SUBROUTINE matrixBPrint(BL2,BR2,TL2,TR2,B)

  integer :: i,j
  COMPLEX, allocatable, dimension(:,:)  :: TL2, TR2, BL2, BR2
  COMPLEX *16      B(LDA,NMAX)

  !write TL

  open(unit=5, file='graphTL2.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(5, '(1600F9.2)')( Real(TL2(i,j)) ,j=0,Ndr)
  end do
  close(5)

  !write TR

  open(unit=66, file='graphTR2.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(66, '(1600F9.2)')( Real(TR2(i,j)) ,j=0,Ndr)
  end do
  close(66)

  !write BL

  open(unit=7, file='graphBL2.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(7, '(1600F9.2)')( Real(BL2(i,j)) ,j=0,Ndr)
  end do
  close(7)

  !write BR

  open(unit=8, file='graphBR2.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(8, '(1600F9.2)')( Real(BR2(i,j)) ,j=0,Ndr)
  end do
  close(8)

  !write the output LHS matrix

  open(unit=12, file='graphRHS.txt', ACTION="write", STATUS="replace")
  do i=0,(2*Ndr)+1
    write(12, '(1600F7.2)')( Real(B(i+1,j+1)) ,j=0,(2*Ndr)+1)
  end do
  close(12)

  !write the XXdr matrix 
  open(unit=13, file='graphXXdr.txt', ACTION="write", STATUS="replace")
  do i=0,Ndr
    write(13, '(1600F9.2)')( Real(XXdr(i,j)) ,j=0,Ndr)
  end do
  close(13)

  END SUBROUTINE matrixBPrint

END MODULE ev



