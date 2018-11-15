Module min_U2_BASIC_mod
  !
  use def_par
  use f95_lapack, only: la_syevr
  !
  implicit none
  Double precision,allocatable:: n_mt(:,:), n2_mt(:,:),Jz2_mt(:,:)
  Integer,allocatable::QN_exp(:)!Quantum number exp
  Double precision,allocatable:: E_exp(:),  sig_exp(:) !Experimental energies
  Integer:: long_exp !Total number of experimental data
  Integer:: QN_max !Maximun QN U(2)
  Double precision,allocatable:: H(:,:),En(:)
  Double precision:: Energy_zero !Zero energy
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
    if(QN_max > N1) then
       write(*,*) "------------------------>ERROR<-------------------------"
       write(*,*) "The higher quantum number of the experimental data file "
       write(*,*) "is bigger than the total Hilbert space dimension N1"
       stop
    endif    
    !
  end Subroutine exp_read
  !
  !***************************************************************************
  !
    Function chi2(a,b,c)
    !*************************************************************************
    !
    !Compute the chi squre
    !
    !*************************************************************************
    implicit none
    Double precision::chi2, E0
    Double precision:: a,b,c
    integer::i,j
    !
    H=0.0d0
    !
    H=a*n_mt + b*n2_mt + c*Jz2_mt
    !
    CALL LA_SYEVR(A=H, W=En, JOBZ='V', UPLO='U')
    !
    E0=minval(En)
    En=En-E0
    !
    chi2=0.0d0
    !
    do j=1,long_exp !Fix the exp energy
       chi2=chi2+(((En(QN_exp(j))-E_exp(j))/sig_exp(j))**2.0d0)
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
    fval=chi2(xval(1),xval(2),xval(3))
    write(*,*) 'CHI2 = ',fval
    !
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
          write(*,101) QN_exp(j),(En(QN_exp(j))-E_exp(j)),(En(QN_exp(j))-E_exp(j))*100.0d0/E_exp(j)
       endif
    enddo
101 format("Estado ",I4," Residual ",D15.5,'(',F7.2,'%)')
    !
  end Subroutine pretty_out
  !
end Module min_U2_BASIC_mod

