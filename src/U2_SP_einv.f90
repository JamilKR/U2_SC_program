Program U2_SP_einv
  !
  use def_par
  use mix_inter
  use f95_lapack, only: la_syevr
  !
  implicit none
  !
  !This program compute the U(2) Shape Coexistence Hamiltonian
  !for a given interaction
  !
  !H=alp1 n + bet1 n**2 + del1 Jz**2 +
  !  alp2 n + bet2 n**2 + del2 Jz**2 +
  !  V_{int}
  !
  !V_{int} are different operator that are built with the funtion V_mix_op(p)
  !in the mix_iter molule.
  !
  Integer, parameter:: ird=5, iwe=6, isav=7 !read, write and save units
  Integer:: MAX !Maximun number of interactions
  Integer,allocatable::INT_orders(:) !Interactions orders
  Double precision,allocatable::INT_coef(:) !Interaction parameters
  Character (len=65)::output_file, namelist_inp
  Logical::Eigenvec
  Double precision,allocatable:: H1(:,:), H2(:,:), V_int(:,:),V_aux(:,:)
  Double precision,allocatable:: H_global(:,:),E_lev(:)
  Double precision:: Energy_zero !Zero energy
  Integer i
  !
  Namelist/OUT/ output_file,Eigenvec
  Namelist/U2dim/ N1, N2
  Namelist/U2par1/ alp1, bet1, del1
  Namelist/U2par2/ alp2, bet2, del2, delta
  Namelist/INT_dim/ MAX
  Namelist/INT_par/ INT_orders, INT_coef
  !
  read(*,*) namelist_inp
  !
  !Read data and allocate arrays
  !
  open(unit=11,file=trim(namelist_inp),status='old',action='read')
  !
  read(11,OUT)
  read(11,U2dim)
  read(11,U2par1)
  read(11,U2par2)
  read(11,INT_dim)
  !
  allocate(INT_orders(1:MAX),INT_coef(1:MAX))
  !
  read(11,INT_par)
  !
  close(11)
  !
  Allocate(H1(0:N1,0:N1),H2(0:N2,0:N2),V_int(0:N1,0:N2),V_aux(0:N1,0:N2))
  !
  ! H1 and H2 building
  !
  call U2_ham(alp1,bet1,del1,0.0d0,N1,H1)
  call U2_ham(alp2,bet2,del2,delta,N2,H2)
  !
  ! Potential boulding
  !
  V_int=0.0d0
  !
  do i=1,MAX
     V_aux=0.0d0
     V_aux= V_mix_op(INT_orders(i))
     V_int=V_int+INT_coef(i)*V_aux
  enddo
  !
  allocate(H_global(0:N1+N2+1,0:N1+N2+1),E_lev(0:N1+N2+1))
  !
  H_global = Global_matrix(H1,H2,V_int)
  !
  !Diagonalization
  !
  CALL LA_SYEVR(A=H_global, W=E_lev, JOBZ='V', UPLO='U')
  !
  Energy_zero=minval(E_lev)
  !
  E_lev=E_lev-Energy_zero
  !
  write(*,*) 'The zero energy is: ', Energy_zero
  write(*,*)
  !
  do i=0,N1+N2+1
     write(*,*) 'Level ',i,' Energy= ', E_lev(i)
  enddo
  !
end Program U2_SP_einv

