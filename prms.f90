MODULE prms

  IMPLICIT none

  integer         :: Ncr
  integer         :: Nbr
  integer         :: Ndr

  integer         :: Nct
  integer         :: Nbt
  integer         :: Ndt

  integer         :: Nfac
  integer         :: Nfr
  integer         :: Nft

  real(16)        :: Rmax
  real(16)        :: Tmax

  real(16)            :: Pi

  Real(16) ::            GAM,RHO1,RHO2,MU1,MU2,aa,HH,KAP,Reynolds,K
  !PARAMETER        (GAM=0.0,RHO1=1.0,RHO2=1.0)

  CONTAINS

  SUBROUTINE initprms

    Pi = 4.*ATAN(1.)

    call readprms

  END SUBROUTINE initprms

  SUBROUTINE readprms

    open(1,file='circle.in')

    read(1,*) Ncr
    read(1,*) Nbr;  Ndr = Ncr + Nbr
    read(1,*) Nct
    read(1,*) Nbt;  Ndt = Nct + Nbt

    read(1,*) Nfac; Nfr = Nfac*Ndr; Nft = Nfac*Ndt

    read(1,*) Tmax
    read(1,*) Rmax

    close(1)

   !Tmax = Tmax/10. !hickox
   Tmax = 0.1 !joseph

   GAM = 0.0
   RHO1 = 1.0
   RHO2 = 1.0
   K = 1.
   aa=RMAX
   HH=RMAX+TMAX
   Reynolds = 100.
   MU2 = HH/Reynolds !joseph
   kap = -4.*MU2/(HH**2) !joseph

   !mu1 = 10.0
   !mu2 = 200.0
   !kap = -1.0


  END SUBROUTINE readprms

END MODULE prms
