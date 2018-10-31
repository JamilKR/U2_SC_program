program test
  !
  use def_par
  use mix_inter
  use f95_lapack, only: la_syevr
  !
  implicit none
  Double precision,allocatable::H(:,:),E(:),V(:,:)
  integer:: i,j
  !
  N1=4
  N2=2
  !
  allocate(V(0:N1,0:N2))
  !
  V = V_mix_op(2)
  !
  do j = 0,N2
     write(*,*) (V(i,j), i=0,N1)
  end do
  !
end program test
