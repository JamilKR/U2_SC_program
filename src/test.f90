program test
  !
  use def_par
  use mix_inter
  use f95_lapack, only: la_syevr
  !
  implicit none
  Double precision,allocatable::H1(:,:),H2(:,:),E(:),V(:,:),H(:,:)
  integer:: i,j
  !
  N1=3
  N2=2
  !
  allocate(V(0:N1,0:N2),H1(0:N1,0:N1),H2(0:N2,0:N2))
  allocate(H(0:N1+N2+1,0:N1+N2+1))
  !
  call U2_ham (1.0d0,1.0d0,1.0d0,0.0d0,N1,H1)
  call U2_ham (1.0d0,1.0d0,1.0d0,1.0d0,N2,H2)
  V = V_mix_op(0) + V_mix_op(N1-N2)
  !
  H = Global_matrix(H1,H2,V)
  !
  do i=0,N1+N2+1
     write(*,*) (H(i,j), j=0,N1+N2+1)
  end do
  !
end program test
