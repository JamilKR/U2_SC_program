Module min_funct_mod
  !
  use def_par
  use mix_inter
  use obs_def
  use f95_lapack, only: la_syevr
  !
  implicit none
  Double precision, allocatable::E_exp(:),sig_exp(:) !Energies and errors exp
  Integer,allocatable::QN_exp(:)!Quantum number exp
  Integer:: long_exp !Total number of experimental data
  Integer:: QN_max !Maximun QN U(2)
  Double precision,allocatable:: H1(:,:), H2(:,:), V_int(:,:),V_aux(:,:)
  Double precision,allocatable::V_int_mtr(:,:,:)
  Double precision,allocatable:: H_global(:,:),E_lev(:)
  Double precision:: Energy_zero !Zero energy
  Double precision, allocatable:: nop_1(:,:),nop_2(:,:),n2op_1(:,:),n2op_2(:,:),&
       Jz2_1(:,:),Jz2_2(:,:) !Operator matrices, compute only one time
  Integer:: MAX !Maximun number of interactions
  Integer,allocatable::INT_orders(:) !Interactions orders
  Double precision,allocatable::INT_coef(:) !Interaction parameters
  !
Contains
  Subroutine exp_read(input_file)
    !*************************************************************************
    !
    !Read experimental data
    !
    !Input:  input_file--> path + name of the energies file
    !
    !Output: All of them defined in global module
    !
    !*************************************************************************
    implicit none
    Character(len=65):: input_file
    Integer:: i
    !
    !Firstly we open the file to count the total lines number
    !
    open(unit=90,file=trim(input_file),status='old',action='read')
    !
    read(90,*) long_exp
    !
    !Allocate arrays
    !
    Allocate(E_exp(1:long_exp),sig_exp(1:long_exp),QN_exp(1:long_exp))
    !
    E_exp=0.0d0
    sig_exp=0.0d0
    QN_exp=0.0d0
    !
    do i=1,long_exp
       read(90,*) E_exp(i),sig_exp(i),QN_exp(i)
    enddo
    !
    QN_max=maxval(QN_exp)
    !
    if(QN_max > N1+N2+1) then
       write(*,*) "------------------------>ERROR<-------------------------"
       write(*,*) "The higher quantum number of the experimental data file "
       write(*,*) "is bigger than the total Hilbert space dimension N1+N2+1"
       stop
    endif    
    !
  end Subroutine exp_read
  !
  !***************************************************************************
  !
  Function chi2(par_fit)
    !*************************************************************************
    !
    !Compute the chi squre
    !
    !*************************************************************************
    implicit none
    Double precision,intent(in)::par_fit(*)
    Double precision::chi2, E0
    integer::i,j
    !
    H1=0.0d0
    H2=0.0d0
    V_int=0.0d0
    !
    H1=par_fit(1)*nop_1 + par_fit(2)*n2op_1 + par_fit(3)*Jz2_1
    H2=par_fit(4)*nop_2 + par_fit(5)*n2op_2 + par_fit(6)*Jz2_2
    do i = 0,N2
       H2(i,i)=H2(i,i)+par_fit(7)
    enddo
    !
    do i=1,max
       V_int=V_int+par_fit(i+7)*V_int_mtr(:,:,i)
    enddo
    !
    H_global=0.0d0
    E_lev=0.0d0
    !
    H_global=Global_matrix(H1,H2,V_int)
    !
    CALL LA_SYEVR(A=H_global, W=E_lev, JOBZ='V', UPLO='U')
    !
    E0=minval(E_lev)
    E_lev=E_lev-E0
    !
    chi2=0.0d0
    !
    do j=1,long_exp !Fix the exp energy
       chi2=chi2+(((E_lev(QN_exp(j))-E_exp(j))/sig_exp(j))**2.0d0)
    enddo
    !
  end Function chi2
  !
  !***************************************************************************
  !
  Subroutine FCN(npar,grad,fval,xval,iflag,chi2)
    !*************************************************************************
    !
    !Minuit function
    !
    !*************************************************************************
    implicit none
    Double Precision:: grad(*),xval(*),fval
    Integer:: iflag,npar
    Double Precision:: chi2,rms
    !
    ! if (iflag .eq. 1)  then
    !    write(*,*) " FCN iflag = ", iflag
    ! endif
    fval=chi2(xval)
    write(*,*) 'CHI2 = ',fval
    ! rms=0.0d0
    ! rms=sqrt(fval/(dble(totdat-npar)))
    ! write(*,*)" ___________________________________ "
    ! write(*,*)"| rms = ",rms," |"
    ! write(*,*)"|___________________________________|"
    if (iflag .eq. 3) then
       write(*,*)
       write(*,*)
       write(*,*)
       call pretty_out()
       write(*,*)
       write(*,*)
       write(*,*)
       rms=0.0d0
       rms=sqrt(fval/(dble(long_exp-npar)))
       write(*,*)" ___________________________________ "
       write(*,*)"| rms = ",rms," |"
       write(*,*)"|___________________________________|"
    endif
 
  end subroutine FCN
  !
  !***************************************************************************
  !
  Subroutine pretty_out()
    !*************************************************************************
    !
    !
    !
    !*************************************************************************
    implicit none
    integer::j
    !
    write(*,*) "***************************************************"
    write(*,*) "********************RESIDUALS**********************"
    write(*,*) "***************************************************"
    do j=1,long_exp !Fix the exp energy
       if (E_exp(j) .ne. 0.0) then
          write(*,101) QN_exp(j),(E_lev(QN_exp(j))-E_exp(j)),(E_lev(QN_exp(j))-E_exp(j))*100.0d0/E_exp(j)
       endif
    enddo
101 format("Estado ",I4," Residual ",D15.5,'(',F7.2,'%)')
    !
  end Subroutine pretty_out
  
end Module min_funct_mod

    
