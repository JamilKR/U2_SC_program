Module mix_inter
  !
  use def_par
  !
  implicit none
contains
    Function V_mix_op(p)
      !*************************************************************************
      !
      !Potential: V= P_1(s^{\dagger})^p (t^{\dagger})^{N_1-N_2-p} P_2
      !
      !Base: {|N1 0>, ..., ..., |N1 N1>, |N2 0>, ..., |N2 N2>} N2<N1
      !
      !Inputs: N1--------> First representation boson number
      !        N2--------> Second representation boson number
      !        p --------> Related with the interaction order
      !                    0 <= p <= N1-N2
      !
      !Output: V_mix_op--> Interaction matrix: N_1+1 rows and N_2+1 columns
      !                    Defined: dim(0:N1,0:N2)
      !
      !*************************************************************************
      implicit none
      Integer,intent(in):: p
      Double precision:: V_mix_op(0:N1,0:N2)
      Double precision:: a,b
      Integer:: n,l
      !
      V_mix_op=0.0d0
      !
      do l=0,N2 !Fix the column
         !
         n = l+N1-N2-p
         a = prod_1(p,l)
         b = prod_2(p,l)
         !
         V_mix_op(n,l) = a*b
      enddo
      !
    end Function V_mix_op
    !
    !***************************************************************************
    !
    Function prod_1(p,l)
      !*************************************************************************
      !
      !Inputs: N1--------> First representation boson number
      !        N2--------> Second representation boson number
      !        p --------> Related with the interaction order
      !                    0 <= p <= N1-N2
      !        l --------> N2-Hilbert Space quantum number
      !                    0 <= l <= N2
      !
      !Output: prod_1= \prod_{i=1}^{N1-N2-p} \sqrt{(l+i)}
      !
      !*************************************************************************
      Integer,intent(in):: p,l
      Double precision:: prod_1
      Integer:: i
      !
      prod_1=1.0d0
      !
      do i=1,N1-N2-p
         prod_1=prod_1*sqrt(dble(l+i))
      enddo
      !
    end Function prod_1
    !
    !***************************************************************************
    !
    Function prod_2(p,l)
      !*************************************************************************
      !
      !Inputs: N2--------> Second representation boson number
      !        p --------> Related with the interaction order
      !                    0 <= p <= N1-N2
      !        l --------> N2-Hilbert Space quantum number
      !                    0 <= l <= N2
      !
      !Output: prod_2= \prod_{j=1}^{p} \sqrt{(N2-l+j)}
      !
      !*************************************************************************
      Integer,intent(in):: p,l
      Double precision:: prod_2
      Integer:: j
      !
      prod_2=1.0d0
      !
      do j=1,p
         prod_2=prod_2*sqrt(dble(N2-l+j))
      enddo
      !
    end Function prod_2
    !
  end Module mix_inter
  
