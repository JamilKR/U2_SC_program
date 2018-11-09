Module obs_def
  !
  use def_par
  use mix_inter
  !
  implicit none
contains
  Function total_U2_IPR(V)
    !*************************************************************************
    !
    !Compute the IPR in the U(2) basis V[N1] \oplut V[N2]
    !
    !Input:  V-------------> Eigenstate coordinates
    !
    !Output: total_U2_IPR--> Inverse Participation Ratio
    !
    !*************************************************************************
    implicit none
    Double precision, intent(in):: V(0:N1+N2+1)
    Double precision:: total_U2_IPR
    integer:: i
    !
    total_U2_IPR=0.0d0
    !
    do i=0,N1+N2+1
       total_U2_IPR = total_U2_IPR + dble(V(i)**4)
    enddo
    !
    total_U2_IPR=1.0d0/total_U2_IPR
    !
  end Function total_U2_IPR
  !
  !***************************************************************************
  !
  Function partial_U2_IPR(V,bs)
    !*************************************************************************
    !
    !Compute the partial IPR in the U(2) basis V[N1] if bs==0 or in V[N2] if
    !bs==1
    !
    !Input:  V-------------> Eigenstate coordinates
    !
    !Output: partial_U2_IPR--> Partial Inverse Participation Ratio
    !
    !*************************************************************************
    implicit none
    Integer,intent(in)::bs
    Double precision, intent(in):: V(0:N1+N2+1)
    Double precision:: partial_U2_IPR
    integer:: i
    !
    partial_U2_IPR=0.0d0
    !
    if (bs==0) then
       do i=0,N1
          partial_U2_IPR = partial_U2_IPR + dble(V(i)**4)
       enddo
       partial_U2_IPR = 1.0d0/partial_U2_IPR
    elseif (bs==1) then
       do i=N1+1,N1+N2+1
          partial_U2_IPR = partial_U2_IPR + dble(V(i)**4)
       end do
       partial_U2_IPR = 1.0d0/partial_U2_IPR
    else
       write(*,*) '-----> ERROR <-----'
       write(*,*) ' bs must be 0 or 1 '
    end if
  end Function partial_U2_IPR

end Module obs_def
