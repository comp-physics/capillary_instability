MODULE data

  USE prms
  USE cheb

  IMPLICIT none

  real(16), allocatable, dimension(:,:)  :: y

CONTAINS

  SUBROUTINE initdata

    allocate(y(0:Ndr,0:Ndt))

  END SUBROUTINE initdata

  

  SUBROUTINE closedata

    deallocate(y)

  END SUBROUTINE closedata

END MODULE data
