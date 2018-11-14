Program U2_SP_min
  !
  use def_par
  use mix_inter
  use obs_def
  use f95_lapack, only: la_syevr
  use min_funct_mod
  !
  implicit none
  !
  !This program fit a set of given energies with their sigmas and quantum numbers
  !to the U(2) Shape-Coexistence Hamiltonian:
  !
  !   H  = P1 H1 P1  +  P2 H2 P2  +  Vint
  !
  !   H1   = alp1 n + bet1 n**2 + del1 Jz**2
  !   H2   = alp2 n + bet2 n**2 + del2 Jz**2 + delta
  !   Vint = P1[\sum_{p} Ap (s^\dagger)**p (t^\dagger)^{N1-N2-p}]P2 + h.c.
  !
  !
  !
  Integer, parameter:: ird=5, iwe=6, isav=7 !read, write and save units
  Character (len=65)::output_file, namelist_inp,exp_file
  Character (len=200)::fixed_par
  Integer i
  Double precision,allocatable:: IPR1(:), IPR2(:), IPRT(:)
  Double precision, allocatable::par_fit(:) !Minimization parameters
  !
  Namelist/OUT/ output_file,exp_file
  Namelist/U2dim/ N1, N2
  Namelist/U2par1/ alp1, bet1, del1
  Namelist/U2par2/ alp2, bet2, del2, delta
  Namelist/INT_dim/ MAX
  Namelist/INT_par/ INT_orders, INT_coef
  Namelist/Mnt/ fixed_par
  !
  read(*,*) namelist_inp
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
  read(11,Mnt)
  !
  close(11)
  !
  !ERRORS DEBUGGER
  if (N2 .ge. N1) then
     write(*,*) "----------->ERROR<----------"
     write(*,*) " N2 must be smaller than N1 "
     stop
  else if (MAX .gt. int(N1-N2+1)) then
     write(*,*) "--------------------------->ERROR<--------------------------"
     write(*,*) " The total number of interactions cant be higher than N1-N2 "
     stop
  else if (maxval(INT_orders) .gt. int(N1-N2)) then
     write(*,*) "--------------------->ERROR<--------------------"
     write(*,*) " The highest order interaction allowed is N1-N2 "
     stop
  endif
  !
  Allocate(H1(0:N1,0:N1),H2(0:N2,0:N2),V_int(0:N1,0:N2),V_aux(0:N1,0:N2),)
  Allocate(IPR1(0:N1+N2+1), IPR2(0:N1+N2+1), IPRT(0:N1+N2+1))
  Allocate(H_global(0:N1+N2+1,0:N1+N2+1),E_lev(0:N1+N2+1))
  Allocate(nop_1(0:N1,0:N1),n2op_1(0:N1,0:N1),Jz2_1(0:N1,0:N1))
  Allocate(nop_2(0:N1,0:N1),n2op_2(0:N1,0:N1),Jz2_2(0:N1,0:N1))
  Allocate(V_int_mtr(0:N1,0:N2,1:MAX)) !Interaction matrixes V(:,:,i)
  Allocate(par_fit(1:7+MAX))
  !
  write(*,*) "Input values read correctly"
  write(*,*)
  write(*,*) "Matrices allocated"
  write(*,*)
  !
  !First we compute operators matrices (only one time)
  !
  nop_1=n_op(N1)
  nop_2=n_op(N2)
  n2op_1=n2_op(N1)
  n2op_2=n2_op(N2)
  Jz2_1=Jz2_op(N1)
  Jz2_2=Jz2_op(N2)
  !
  do i=1,max
     V_int_mtr(:,:,i)=V_mix_op(INT_orders(i))
  enddo
  !
  write(*,*) "Operators matrices computed"
  write(*,*)
  !
  !Minuit input file:
  !
  open(unit=12,file='minuit_input.inp',status='replace')
  inquire(12)
  !
  write(12,'(A9)')  "SET TITLE"
  write(12,'(A21)') "'Minuit minimization'"
  write(12,'(A10)') "PARAMETERS"
  write(12,20) alp1
  write(12,21) bet1
  write(12,22) del1
  write(12,23) alp2
  write(12,24) bet2
  write(12,25) del2
  write(12,26) delta
  !
20 format("1        'alp1 ' ",D20.10,"          0.1D-02")
21 format("2        'bet1 ' ",D20.10,"          0.1D-02")
22 format("3        'del1 ' ",D20.10,"          0.1D-02")
23 format("4        'alp2 ' ",D20.10,"          0.1D-02")
24 format("5        'bet2 ' ",D20.10,"          0.1D-02")
25 format("6        'del2 ' ",D20.10,"          0.1D-02")
26 format("7        'delta ' ",D20.10,"          0.1D-02")
  !
  do i=1,MAX
     if (i+7 .lt. 10) then
        write(12,27) i+7,INT_orders(i),INT_coef(i)
     else
        write(12,28) i+7,INT_orders(i),INT_coef(i)
     endif
  enddo
  !
27 format(I1,"       'INT_p",I3," ' ",D20.10,"          0.1D-02")
28 format(I2,"      'INT_p",I3," ' ",D20.10,"          0.1D-02")
  !
  write(12,*)
  write(12,'(a)') adjustl(fixed_par) 
  write(12,'(A10)') 'set stra 2'
  write(12,'(A14)') 'minimize 10000'
  write(12,'(A6)') 'call 3'
  write(12,'(A4)') 'exit'
  write(12,*)
  write(12,*)
  close(12)
  !
  write(*,*) "Minuit input file created"
  write(*,*)
  !
  !Experimental data read:
  !
  write(*,*) "Input energies should be sorted starting by the lowest one"
  write(*,*)
  !
  call exp_read(exp_file)
  !
  write(*,*) "Experimental data read"
  write(*,*)
  !
  !Parameters inizializations
  !
  par_fit(1)=alp1
  par_fit(2)=bet1
  par_fit(3)=del1
  par_fit(4)=alp2
  par_fit(5)=bet2
  par_fit(6)=del2
  par_fit(7)=delta
  do i=8,7+max
     par_fit(i)=INT_coef(i)
  enddo
  !
  write(*,*) "Parameters inicialized"
  write(*,*)
  !
  open(unit=ird,file="minuit_input.inp",status='old')
  open(unit=iwe,file=trim(output_file),status='replace')
  !
  !Call Minuit:
  call mintio(ird,iwe,isav)
  call minuit(FCN,chi2)
  !
  close(ird,status='delete')
  close(iwe)
  !
end Program U2_SP_min
