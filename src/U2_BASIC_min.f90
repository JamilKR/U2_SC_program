Program U2_BASIC_min
  !
  use def_par
  use f95_lapack, only: la_syevr
  use min_U2_BASIC_mod
  !
  implicit none
  !
  !This program fit a U(2) Hamiltonian to a data set
  !
  !   H = alp1 n + bet1 n**2 + del1 Jz**2
  !
  !   Basis= {|N1 n> n=0,...N1}
  !
  Integer, parameter:: ird=5, iwe=6, isav=7 !read, write and save units
  Character (len=65)::output_file, namelist_inp,exp_file
  Character (len=200)::fixed_par
  !
  Namelist/OUT/ output_file, exp_file
  Namelist/U2/ N1, alp1, bet1, del1
  Namelist/Mnt/fixed_par
  !
  read(*,*) namelist_inp
  !
  open(unit=11,file=trim(namelist_inp),status='old',action='read')
  !
  read(11,OUT)
  read(11,U2)
  read(11,Mnt)
  !
  close(11)
  !
  write(*,*) 'Namelist read'
  write(*,*)
  !
  Allocate(n_mt(0:N1,0:N1),n2_mt(0:N1,0:N1))
  Allocate(H(0:N1,0:N1),En(0:N1))
  !
  write(*,*) 'Matrices allocated'
  write(*,*)
  !
  n_mt=n_op(N1)
  n2_mt=n2_op(N1)
  Jz2_mt=Jz2_op(N1)
  !
  write(*,*) "Operators matrices computed"
  write(*,*)
  !
  !Minuit input file:
  !
  open(unit=12,file='minuit_input_U2.inp',status='replace')
  inquire(12)
  !
  write(12,'(A9)')  "SET TITLE"
  write(12,'(A21)') "'Minuit minimization'"
  write(12,'(A10)') "PARAMETERS"
  write(12,20) alp1
  write(12,21) bet1
  write(12,22) del1
  !
20 format("1        'alp1 ' ",D20.10,"          0.1D-02")
21 format("2        'bet1 ' ",D20.10,"          0.1D-02")
22 format("3        'del1 ' ",D20.10,"          0.1D-02")
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
  write(*,*) output_file
  call exp_read(exp_file)
  !
  write(*,*) "Experimental data read"
  write(*,*)
  !
  
  open(unit=ird,file="minuit_input_U2.inp",status='old')
  open(unit=iwe,file=trim(output_file),status='replace')
  !
  !Call Minuit:
  call mintio(ird,iwe,isav)
  call minuit(FCN,chi2)
  !
  close(ird,status='delete')
  close(iwe)
  !
end Program U2_BASIC_min
