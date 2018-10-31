Module def_par
  !
  !Module with global U(2) variables, functions and subroutines
  !
  implicit none
  !
  Integer:: N1, N2, dim
  !Ni: U(2) of i-representation
  !dim: global matrix dimension dim=N1+1+N2+1
  Double precision:: alp1, alp2, bet1, bet2, del1, del2, delta
  !Hamiltonian terms: H1(alp1,bet1,del1)
  !                   H2(alp2,bet2,del2,delta)
  !                   Vint(A,B)
  !
Contains
    !
  Subroutine U2_ham(alp,bet,del,cte,N,H_mx) !------------------------>[Tested]
    !*************************************************************************
    !
    !Build the U(2) Hamiltonian matrix with dimension N+1
    !
    ! H=cte + alp*n + bet*n**2 + del*Jz**2
    !
    !Inputs: alp, bet, del, cte, N
    !
    !Outputs: H
    !
    !*************************************************************************
    implicit none
    Double precision:: alp, bet, del, cte
    Integer:: N
    Double precision, allocatable:: n_mx(:,:), n2_mx(:,:), Jz2_mx(:,:)
    Double precision::H_mx(0:N,0:N)
    Integer:: i,j
    !
    allocate(n_mx(0:N,0:N), n2_mx(0:N,0:N), Jz2_mx(0:N,0:N))
    !
    n_mx = n_op(N)
    n2_mx = n2_op(N)
    Jz2_mx =  Jz2_op(N)
    !Operators contributions
    H_mx = alp*n_mx + bet*n2_mx + del*Jz2_mx
    !
    do i=0,N !Constant contribution
       H_mx(i,i) = H_mx(i,i) + cte
    end do
    !
  end Subroutine U2_ham
  !
  !***************************************************************************
  !
  Function n_op(N) !------------------------------------------------->[Tested]
    !*************************************************************************
    !
    !n operator matrix in U(1) basis
    !
    !*************************************************************************
    implicit none
    Integer, intent(in):: N
    Double precision:: n_op(0:N,0:N)
    Integer:: i
    !
    n_op=0.0d0
    !
    do i=0,N
       n_op(i,i)=dble(i)
    enddo
      !  
  end Function n_op
  !
  !***************************************************************************
  !
  Function n2_op(N) !------------------------------------------------>[Tested]
    !*************************************************************************
    !
    !n**2 operator matrix in U(1) basis
    !
    !*************************************************************************
    implicit none
    Integer, intent(in):: N
    Double precision:: n2_op(0:N,0:N)
    Integer:: i
    !
    n2_op=0.0d0
    !
    do i=0,N
       n2_op(i,i)=dble(i)*dble(i)
    enddo
    !  
  end Function n2_op
  !
  !***************************************************************************
  !
  Function Jz2_op(N) !----------------------------------------------->[Tested]
    !*************************************************************************
    !
    !Jz**2 operator matrix in U(1) basis
    !
    !*************************************************************************
    implicit none
    Integer, intent(in):: N
    Double precision:: Jz2_op(0:N,0:N)
    Integer i
    !
    Jz2_op=0.0d0
    !
    do i = 0, N !Mean diagonal
       Jz2_op(i,i) = 0.25d0 * ( dble(N) + 2.0d0*(dble(i)*dble(N-i)) )
    enddo
    !
    do i = 0,N-2 !First diagonal
       Jz2_op(i,i+2) = -0.25d0 * sqrt(dble(N-i-1)*dble(N-i)*dble(i+1)*dble(i+2))
       Jz2_op(i+2,i) = Jz2_op(i,i+2)
    enddo
    !
  end Function Jz2_op
  !
  !***************************************************************************
  !
end Module def_par
  
