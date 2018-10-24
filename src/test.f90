program test
  use def_par
  use f95_lapack, only: la_syevr
  implicit none
  Double precision,allocatable::H(:,:),E(:)
  integer:: i,j
  !
  alp1=0.0
  bet1=0.0
  del1=1.0
  N1=4
  !
  allocate(H(0:N1,0:N1),E(0:N1))
  !
  call U2_ham(alp1,bet1,del1,0.0d0,N1,H)
  !
  do i = 0,N1
     write(*,*) (H(i,j),j=0,N1)
  end do
  !
  call la_syevr(H,E,'V','U')
  !
  do i = 0,N1
     write(*,*) E(i)
  end do
  !
end program test
