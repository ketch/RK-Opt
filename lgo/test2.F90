module globals
    implicit none
    integer, parameter :: s=2
    integer, parameter :: p=2
    character, parameter :: RKclass*3 = 'erk'
end module globals

program main
    use globals
    real*8 obj
    real*8, allocatable :: con(:), x(:)
    integer n1, n2, n3, ncons, nvars

    call set_ncons(n1,n2,n3)
    ncons = n1+n2+n3
    allocate(con(ncons))

    call set_nvars(nvars)
    allocate(x(nvars))
    x = [0.d0,1.d0,0.5d0,0.5d0,1.d0,1.d0]

    !construct some dummy data and call obj_and_constraints()

    call obj_and_constraints(x,obj,con)
    call obj_and_constraints(x,obj,con)

endprogram main


subroutine set_ncons(n_order_conditions,n_consistency_conditions, n_am_conditions)

  use globals
  implicit none
  integer n_order_conditions, n_consistency_conditions, n_am_conditions

  ! Number of order conditions
  select case(p)
    case(1)
        n_order_conditions = 1
    case(2)
        n_order_conditions = 2
    case(3)
        n_order_conditions = 4
    case(4)
        n_order_conditions = 8
    case(5)
        n_order_conditions = 17
    case(6)
        n_order_conditions = 37
    case(7)
        n_order_conditions = 85
  end select

  n_consistency_conditions = s

  n_am_conditions = (s+1)*(s+2)

end subroutine set_ncons


subroutine obj_and_constraints(x, obj, con)
    use globals
    implicit none

    real*8 x(*)
    real*8 obj
    real*8 con(*)
    real*8 A(s,s), b(s), c(s)
    real*8, allocatable :: order_conditions(:)
    real*8 K(s+1,s+1), kiprkinv(s+1,s+1), Kpow(s+1,s+1)
    real*8 stage_consistency(s)
    real*8 r
    real*8 am_conditions_1((s+1)**2), am_conditions_2(s+1)
    integer ncons
    
    integer n_consistency_conditions, n_order_conditions
    integer n_am_conditions

    integer nvars, i

    call set_nvars(nvars)
    call set_ncons(n_order_conditions,n_consistency_conditions,n_am_conditions)

    
    ncons = n_order_conditions + n_consistency_conditions + n_am_conditions

    call unpack_rk(x,A,b,c)

    r   = x(nvars)
    obj = -r

    allocate(order_conditions(n_order_conditions))

    ! Compute order conditions
    call oc_butcher(A,b,c,order_conditions)

    ! Compute consistency conditions
    stage_consistency = c - sum(A,2)

    ! Compute abs. mon. conditions
    K(:,:) = 0.d0
    K(1:s,1:s) = A
    K(s+1,1:s) = b

    kiprkinv(:,:) = 0.d0
    Kpow(:,:) = 0.d0
    ! (I + rK)^(-1) = I - rK + (rK)^2 - (rK^3) + ...
    do i = 1,s+1
        kiprkinv(i,i) = 1.d0
        Kpow(i,i)    = 1.d0
    end do
    do i = 1, s
        Kpow = matmul(K,Kpow)
        kiprkinv = kiprkinv + (-r)**i*Kpow
    end do
    kiprkinv = matmul(K,kiprkinv)

    ! K(I+rK)^{-1} >= 0
    am_conditions_1 = -reshape(kiprkinv, [(s+1)**2])
    ! rK(I+rK)^{-1}e <= 1
    am_conditions_2 = r*sum(kiprkinv,2)-1.d0

    ! Concatenate all constraints, equalities first
    con(1:ncons) = [order_conditions,stage_consistency,am_conditions_1, &
            am_conditions_2]

end subroutine obj_and_constraints


subroutine set_nvars(n)
    ! Set total number of decision variables
    ! Note that (unlike the MATLAB code) we store all s
    ! abscissae for explicit methods

    use globals
    implicit none
    integer n

    select case(RKclass)
        !=====================
        ! RK classes
        !=====================
        case ('erk') ! Explicit (A is strictly lower triangular)
          n=2*s+s*(s-1)/2 + 1;
    end select
        !case 'irk'  %Fully implicit
        !  n=2*s+s^2+1;
        !case 'irk5'  %Implicit, p>=5 (so first row of A is zero)
        !  n=2*s+s*(s-1);
        !case 'dirk' %Diagonally Implicit (A is lower triangular)
        !  n=2*s+s*(s+1)/2+1;
        !case 'dirk5' %Diagonally Implicit, p>=5 (lower tri. and first row of A is zero)
        !  n=2*s+s*(s+1)/2-1;
        !case 'sdirk' %Singly Diagonally Implicit
        !  n=2*s+s*(s+1)/2+1-s+1;
    return
end subroutine set_nvars

subroutine unpack_rk(X,A,b,c)
! Extracts the coefficient arrays from the vector of decision variables.
!
! The coefficients are tored in a single vector x as::
!
!       x=[c' b' A]
!
! A is stored row-by-row.
! Note that we always store all of c (unlike the MATLAB code).

    use globals
    implicit none

    real*8 A(s,s)
    real*8 b(s,1)
    real*8 c(s,1)
    real*8 X(*)
    integer i
    real*8 K(s+1,s+1)

    A(:,:) = 0.d0

    ! Generate Butcher coefficients (for RK) or Shu-Osher coefficients (for lsRK)
    select case(RKclass)
        case ('erk')
            c = reshape( [X(1:s)], [s,1])
            b = reshape( X(s+1:2*s) , [s,1])
            do i = 1, s
                A(i,1:i-1)=X(2*s+1+(i-2)*(i-1)/2:2*s+i*(i-1)/2);
            end do
    end select

            K(:,:) = 0.d0
            K(1:s,1:s) = A
            K(s+1,1:s) = b(:,1)
    return
end
!    case ('irk')
!        c=X(1:s)'; b=X(s+1:2*s)';
!        A=reshape(X(2*s+1:2*s+s^2),s,s)';
!    case ('irk5')
!        c=[0 X(1:s-1)]'; b=X(s:2*s-1)';
!        for i=2:s
!            A(i,:)=X(i*s:(i+1)*s-1);
!        end
!    case ('dirk')
!        c=X(1:s)'; b=X(s+1:2*s)';
!        for i=1:s
!            A(i,1:i)=X(2*s+1+i*(i-1)/2:2*s+i*(i+1)/2);
!        end
!    case ('dirk5')
!        c=[0 X(1:s-1)]'; b=X(s:2*s-1)'; 
!        for i=2:s
!            A(i,1:i)=X(2*s-1+i*(i-1)/2:2*s-2+i*(i+1)/2);
!        end
!    case ('sdirk')
!        c=X(1:s)'; b=X(s+1:2*s)';
!        A(1,1)=X(2*s+1);
!        for i=2:s
!            A(i,1:i-1)=X(2*s+2+(i-2)*(i-1)/2:2*s+1+i*(i-1)/2);
!            A(i,i)=X(2*s+1);
!        end


subroutine oc_butcher(A,b,c,coneq)

    use globals
    real*8 A(s,s), b(s), c(s)
    real*8 coneq(*)

    if (p>=1) then
      ! order 1 conditions:
        coneq(1)=sum(b)-1d0
    endif

    if (p>=2) then
      ! order 2 conditions:
        coneq(2)=dot_product(b,c)-1/2d0
    endif

    if (p>=3) then
      ! order 3 conditions:
        coneq(3)=dot_product(b,matmul(A,c))-1/6d0
        coneq(4)=dot_product(b,c**2)-1/3d0
    endif

    if (p>=4) then
      ! order 4 conditions:
        coneq(5)=dot_product(b,matmul(A,c**2))-1/12d0
        coneq(6)=dot_product(b,matmul(A,matmul(A,c)))-1/24d0
        coneq(7)=dot_product(b,matmul(A,c)*c)-1/8d0
        coneq(8)=dot_product(b,c**3)-1/4d0
    endif

    if (p>=5) then
      ! order 5 conditions:
        coneq(9)=dot_product(b,matmul(A,c**3))-1/20d0
        coneq(10)=dot_product(b,matmul(A,matmul(A,c**2)))-1/60d0
        coneq(11)=dot_product(b,matmul(A,matmul(A,matmul(A,c))))-1/120d0
        coneq(12)=dot_product(b,matmul(A,c*matmul(A,c)))-1/40d0
        coneq(13)=dot_product(b,matmul(A,c**2)*c)-1/15d0
        coneq(14)=dot_product(b,matmul(A,matmul(A,c))*c)-1/30d0
        coneq(15)=dot_product(b,matmul(A,c)*c**2)-1/10d0
        coneq(16)=dot_product(b,matmul(A,c)*matmul(A,c))-1/20d0
        coneq(17)=dot_product(b,c**4)-1/5d0
    endif

    if (p>=6) then
      ! order 6 conditions:
        coneq(18)=dot_product(b,matmul(A,c**4))-1/30d0
        coneq(19)=dot_product(b,matmul(A,matmul(A,c**3)))-1/120d0
        coneq(20)=dot_product(b,matmul(A,matmul(A,matmul(A,c**2))))-1/360d0
        coneq(21)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,c)))))-1/720d0
        coneq(22)=dot_product(b,matmul(A,matmul(A,c*matmul(A,c))))-1/240d0
        coneq(23)=dot_product(b,matmul(A,c*matmul(A,c**2)))-1/90d0
        coneq(24)=dot_product(b,matmul(A,c*matmul(A,matmul(A,c))))-1/180d0
        coneq(25)=dot_product(b,matmul(A,c**2*matmul(A,c)))-1/60d0
        coneq(26)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,c)))-1/120d0
        coneq(27)=dot_product(b,matmul(A,c**3)*c)-1/24d0
        coneq(28)=dot_product(b,matmul(A,matmul(A,c**2))*c)-1/72d0
        coneq(29)=dot_product(b,matmul(A,matmul(A,matmul(A,c)))*c)-1/144d0
        coneq(30)=dot_product(b,matmul(A,c*matmul(A,c))*c)-1/48d0
        coneq(31)=dot_product(b,matmul(A,c**2)*c**2)-1/18d0
        coneq(32)=dot_product(b,matmul(A,matmul(A,c))*c**2)-1/36d0
        coneq(33)=dot_product(b,matmul(A,c)*c**3)-1/12d0
        coneq(34)=dot_product(b,matmul(A,c**2)*matmul(A,c))-1/36d0
        coneq(35)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,c))-1/72d0
        coneq(36)=dot_product(b,matmul(A,c)*matmul(A,c)*c)-1/24d0
        coneq(37)=dot_product(b,c**5)-1/6d0
    endif

    if (p>=7) then
      ! order 7 conditions:
        coneq(38)=dot_product(b,matmul(A,c**5))-1/42d0
        coneq(39)=dot_product(b,matmul(A,matmul(A,c**4)))-1/210d0
        coneq(40)=dot_product(b,matmul(A,matmul(A,matmul(A,c**3))))-1/840d0
        coneq(41)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,c**2)))))-1/2520d0
        coneq(42)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,matmul(A,c))))))-1/5040d0
        coneq(43)=dot_product(b,matmul(A,matmul(A,matmul(A,c*matmul(A,c)))))-1/1680d0
        coneq(44)=dot_product(b,matmul(A,matmul(A,c*matmul(A,c**2))))-1/630d0
        coneq(45)=dot_product(b,matmul(A,matmul(A,c*matmul(A,matmul(A,c)))))-1/1260d0
        coneq(46)=dot_product(b,matmul(A,matmul(A,c**2*matmul(A,c))))-1/420d0
        coneq(47)=dot_product(b,matmul(A,matmul(A,matmul(A,c)*matmul(A,c))))-1/840d0
        coneq(48)=dot_product(b,matmul(A,c*matmul(A,c**3)))-1/168d0
        coneq(49)=dot_product(b,matmul(A,c*matmul(A,matmul(A,c**2))))-1/504d0
        coneq(50)=dot_product(b,matmul(A,c*matmul(A,matmul(A,matmul(A,c)))))-1/1008d0
        coneq(51)=dot_product(b,matmul(A,c*matmul(A,c*matmul(A,c))))-1/336d0
        coneq(52)=dot_product(b,matmul(A,c**2*matmul(A,c**2)))-1/126d0
        coneq(53)=dot_product(b,matmul(A,c**2*matmul(A,matmul(A,c))))-1/252d0
        coneq(54)=dot_product(b,matmul(A,c**3*matmul(A,c)))-1/84d0
        coneq(55)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,c**2)))-1/252d0
        coneq(56)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,matmul(A,c))))-1/504d0
        coneq(57)=dot_product(b,matmul(A,c*matmul(A,c)*matmul(A,c)))-1/168d0
        coneq(58)=dot_product(b,matmul(A,c**4)*c)-1/35d0
        coneq(59)=dot_product(b,matmul(A,matmul(A,c**3))*c)-1/140d0
        coneq(60)=dot_product(b,matmul(A,matmul(A,matmul(A,c**2)))*c)-1/420d0
        coneq(61)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,c))))*c)-1/840d0
        coneq(62)=dot_product(b,matmul(A,matmul(A,c*matmul(A,c)))*c)-1/280d0
        coneq(63)=dot_product(b,matmul(A,c*matmul(A,c**2))*c)-1/105d0
        coneq(64)=dot_product(b,matmul(A,c*matmul(A,matmul(A,c)))*c)-1/210d0
        coneq(65)=dot_product(b,matmul(A,c**2*matmul(A,c))*c)-1/70d0
        coneq(66)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,c))*c)-1/140d0
        coneq(67)=dot_product(b,matmul(A,c**3)*c**2)-1/28d0
        coneq(68)=dot_product(b,matmul(A,matmul(A,c**2))*c**2)-1/84d0
        coneq(69)=dot_product(b,matmul(A,matmul(A,matmul(A,c)))*c**2)-1/168d0
        coneq(70)=dot_product(b,matmul(A,c*matmul(A,c))*c**2)-1/56d0
        coneq(71)=dot_product(b,matmul(A,c**2)*c**3)-1/21d0
        coneq(72)=dot_product(b,matmul(A,matmul(A,c))*c**3)-1/42d0
        coneq(73)=dot_product(b,matmul(A,c)*c**4)-1/14d0
        coneq(74)=dot_product(b,matmul(A,c**2)*matmul(A,c**2))-1/63d0
        coneq(75)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,c**2))-1/126d0
        coneq(76)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,matmul(A,c)))-1/252d0
        coneq(77)=dot_product(b,matmul(A,c**3)*matmul(A,c))-1/56d0
        coneq(78)=dot_product(b,matmul(A,matmul(A,c**2))*matmul(A,c))-1/168d0
        coneq(79)=dot_product(b,matmul(A,matmul(A,matmul(A,c)))*matmul(A,c))-1/336d0
        coneq(80)=dot_product(b,matmul(A,c*matmul(A,c))*matmul(A,c))-1/112d0
        coneq(81)=dot_product(b,matmul(A,c**2)*matmul(A,c)*c)-1/42d0
        coneq(82)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,c)*c)-1/84d0
        coneq(83)=dot_product(b,matmul(A,c)*matmul(A,c)*c**2)-1/28d0
        coneq(84)=dot_product(b,matmul(A,c)*matmul(A,c)*matmul(A,c))-1/56d0
        coneq(85)=dot_product(b,c**6)-1/7d0
    endif

return
end subroutine
