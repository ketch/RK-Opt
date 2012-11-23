
! SSP model formulation received from David Ketcheson, KAUST, October 16, 2012 

       module globals
            implicit none
            integer, parameter :: s=5
            integer, parameter :: p=4
            character, parameter :: RKclass*3 = 'erk'
            real*8, allocatable :: order_conditions(:)           

            contains
                subroutine alloc_oc(n_order_conditions)
                    integer n_order_conditions

                    allocate(order_conditions(n_order_conditions)) 
                end subroutine

        end module globals
    
        program main
!     LGOMAIN.F90 --- Main program file for SSP optimization
  
!     Basic program structure:
!     LGOMAIN calls LGO_RUN; the latter calls USER_FCT (in LGOFCT.FOR) 
!     iteratively.
!     ------------------------------------------------------------------

!     This main program includes the explicit definition of all LGO I/O 
!     arguments. The results are optionally written also to output files.
 
!     The detailed example below illustrates the input information setup.

! --- Declarations ---
  
        use globals

      IMPLICIT REAL*8 (a-h,o-z)
!     Maximal currently handled model size: maxdim variables, maxcon con
      PARAMETER  (maxdim=5000,maxcon=2000)
      CHARACTER*20 modname, varname(maxdim), conname(maxcon), objname
      DIMENSION varlb(maxdim), varnom(maxdim), varub(maxdim)
      DIMENSION varopt(maxdim), conopt(maxcon)
      INTEGER opmode, g_maxfct, tlimit, sstat, ctype(maxcon)
      REAL*8 kt_tol, l_target
      common / modsize/ nvars, ncons
!     CHARACTER*40 sst(100),mst(100)

      integer ncons, nvars

      integer n_consistency_conditions, n_order_conditions
      integer n_am_conditions

! --- Model descriptors ---

! Model name (character)
      modname='SSP RK'


! Number of variables (positive integer)
      call set_nvars(nvars)

! Number of general constraints (non-negative integer)
      call set_ncons(n_order_conditions,n_consistency_conditions,n_am_conditions)
      neqcons = n_order_conditions + n_consistency_conditions
      ncons=neqcons + n_am_conditions


      call alloc_oc(n_order_conditions)
! Default variable bound, set for variables without explicit bounds
      defbnd=5.d0 

      varlb(:) = 0.d0
      varub(:) = 0.d0
      varnom(:) = 0.d0
      varopt(:) = 0.d0
      conopt(:) = 0.d0
      ctype(:) = 0.d0
      varname(:) = ''
      conname(:) = ''

! Variable lower bounds (real)
      do i=1,nvars
        varlb(i)=0.d0
      end do  

! Variable upper bounds (real)
      do i=1,nvars-1
        varub(i)=1.d0
      end do     
      varub(nvars) = s-p+1.d0

! Variable nominal values (real)
      do i=1,nvars
        varnom(i)=0.5d0*(varlb(i)+varub(i))
!       varnom(i)=dble(i)
!     if (i.eq.nvars) varnom(i)=1.        
      end do
     
! Objective function name (character)
      objname='r'

      do i=1,nvars
        varname(i)=''
      end do  

! Constraint names (character)
      do i=1,ncons
        conname(i)=''
      end do  

! Constraint types: equations 0, <= inequalities -1
      do i=1,neqcons
         ctype(i)=0
      end do

      do i=neqcons+1, ncons
         ctype(i)=-1
      end do      

!      write(*,*) n_order_conditions, n_consistency_conditions, n_am_conditions
!      write(*,*) nvars, ncons, ctype
!      stop

! --- Solver options and parameters ---
  
! Operational mode selected: 0: LS; 1: B&B+LS; 2: GARS+LS; 3: MS+LS
! Suggested default for global search: 3; for local search only: 0
!     opmode=0
!     opmode=1
!     opmode=2        
      opmode=3
  
! Maximal no. of model function evaluations in global search phase (integer)
! Suggested optional default settings
      !g_maxfct=100
      g_maxfct=10000*(nvars+ncons+1)**2
  
! Maximal no. of model function evaluations in each executed local search 
! phase (integer)
! Suggested optional default setting
      l_maxfct=10000*(nvars+ncons+1)**2        
  
! Maximal no. of model function evaluations in global search phase, w/o improvement (integer)
! Suggested default settings: g_maxfct/10 <= max_nosuc <= g_maxfct        
! Suggested optional default setting
      max_nosuc=5000*(nvars+ncons+1)**2    
  
! Constraint penalty multiplier (positive real, used in merit function evaluation)
! Suggested default setting: 1.; one can incresase it to enforce feasibility       
      pen_mult=100.d0
  
! Target objective function value in global search phase (real)
! Suggested default setting: "close" to optimum value when known
! Otherwise can be set to an "unrealistic" value such as e.g. -1.d30
      g_target=-8.
  
! Target objective function value in local search phase (real)
! Suggested default setting: "even closer" to optimum value when known
! l_target <= g_target settings are expected
! Otherwise can be set to an "unrealistic" value such as e.g. -1.d20
      l_target=-9.
  
! Required local search precision parameter (positive real)
! Suggested default setting: 1.d-8
      fi_tol=1.d-8
  
! Constraint violation tolerance in final solution (positive real)
! Suggested default setting: 1.d-8
!     con_tol=1.d-8
! in this model we want high precision to avoid the many pseudo-solutions
      con_tol=1.d-8 
  
! Kuhn-Tucker local optimality conditions tolerance in local search phase
! (positive real)
! Suggested default setting: 1.d-8
      kt_tol=1.d-8
  
! Seed value for built-in random number generator (non-negative integer)
! Suggested default setting: 0
      irngs=0
  
! Program execution time limit (seconds, positive integer)
! Suggested default setting: 300, but large-size and/or difficult models 
! may require longer runs
      tlimit=3000
  
! Optional result and logfile printout level 
! 0: no printout of results 
! 1: summary printout to file LGO.SUM (unit 10) 
! 2: additional more detailed printout to file LGO.OUT (unit 9)
! 3: additional printout in all iterations: input argument and return 
!    function values, to file LGO.LOG (unit 11)
! Suggested default setting: 1; use 2 for more information, or 3 in 
! checking and troubleshooting sessions
      iprl=1
 
! --- LGO solver call ---
  
      write(*,*) nvars, ncons
      CALL lgo_run (modname, nvars, ncons, varname, objname, conname, &
     ctype, defbnd, varlb, varnom, varub, opmode, g_maxfct, l_maxfct, &
     max_nosuc, pen_mult, g_target, l_target, fi_tol, con_tol, kt_tol, &
     irngs, tlimit, iprl, varopt, conopt, objopt, resmax, nevals, &
     runtime, sstat, mstat)
  
! --- LGO return values and status indicators ---
 
!     varopt: optimized variable values (real array)  
!     conopt: constraint function values at varopt (real array)  
!     objopt: objective function value at varopt (real)  
!     resmax: maximal constraint violation at varopt (real)  
!     nevals: total no. of model function evaluations (integer)  
!     runtime: LGO runtime, given to 0.01 sec precision (real)

!     sstat: solver status (integer)
!     Solver status return values
!     sst(1)='Normal completion'        ! Solver terminated with (expectedly) good result    
!     sst(2)='Iteration interrupt'      ! Interrupt by user
!     sst(3)='Resource interrupt'       ! Preset no. of FCT evals or time limit (both set by user exceeded)
!     sst(4)='Run terminated by solver' ! Due to set solver parameters or solver limitations     
!     sst(7)='Licensing error'          ! Built-in size limitations are exceeded by the user's model 

!     mstat: model status (integer)     
!     Model status return values
!     mst(1)='Global search based feasible solution' ! Solution is generated by global+local search 
!     mst(2)='Local search based feasible solution'  ! Solution is generated by local search only 
!     mst(3)='Unbounded solution (error)'   ! In presence of variable bounds, this implies model error
!     mst(4)='Infeasible solution'          ! Solution found does not meet the set constraint violation tolerance
!     mst(7)='Intermediate solution'        ! This value is returned when sstat=2, 3, 4, or 7

!     Verification printout to file 
!     This can (and typically will) be omitted in LGO solver versions
!     linked to other modeling environments: AIMMS, GAMS, AMPL, MPL, ...
      OPEN(unit=12, file="LGOresults", status="unknown")
      n=nvars
      m=ncons
      WRITE (12,250) modname,nevals,objname,objopt, &
      (i,varname(i),varopt(i),i=1,n)
250    FORMAT (/,5x,'--- LGO Solver Results Summary ---',//,5x, &
      'Model name ',a20,//,5x,'Number of function evaluations', &
      15x,i15,//,5x,'Objective function name ',a20,//,5x, &
      'Objective function value',16x,f20.10,//,5x, &
      'Variables',/,(5x,i9,3x,a20,8x,f20.10))
      IF (m.gt.0) THEN
        WRITE (12,270) (j,conname(j),conopt(j),j=1,m)
270     FORMAT (/,5x,'Constraints', &
       /,(5x,i9,3x,a20,8x,f20.10))         
        WRITE (12,290) resmax
290     FORMAT (/,5x,'Maximal constraint violation',12x,f20.10)        
      END IF
!     WRITE( 9,103) sstat, sst(sstat), mstat, mst(mstat)          
!103   FORMAT(/,5x,'Solver status indicator value ',i1, 1x, a40, &
!        //,5x,'Model status indicator value  ',i1, 1x, a40)
      WRITE(12,300) sstat, mstat          
300   FORMAT(/,5x,'Solver status indicator value ',i1, &
        //,5x,'Model status indicator value  ',i1)
      WRITE (12,310) runtime           
310   FORMAT (/,5x,'LGO solver system execution time (seconds)',3x, &
       f15.2)      
      CLOSE (12) 

      STOP
      end program main

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


subroutine USER_FCT(x, obj, con)
    use globals
    implicit none

    real*8 x(*)
    real*8 obj
    real*8 con(*)
    real*8 A(s,s), b(s), c(s)
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
!    am_conditions_1 = -reshape(kiprkinv, [(s+1)**2])
    am_conditions_1 = -reshape( kiprkinv, (/ (s+1)**2 /) )
    ! rK(I+rK)^{-1}e <= 1
    am_conditions_2 = r*sum(kiprkinv,2)-1.d0

    ! Concatenate all constraints, equalities first
!    con(1:ncons) = [order_conditions,stage_consistency,am_conditions_1, &
!            am_conditions_2]
    con(1:ncons) = (/ order_conditions,stage_consistency,am_conditions_1, &
             am_conditions_2 /)
end subroutine ! objective_and_constraints


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
!            c = reshape( [X(1:s)], [s,1])
!            b = reshape( X(s+1:2*s) , [s,1])
            c = reshape( (/ X(1:s) /), (/ s,1 /) )
            b = reshape( X(s+1:2*s) , (/ s,1 /) )
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
        coneq(1)=sum(b)-1.d0
    endif

    if (p>=2) then
      ! order 2 conditions:
        coneq(2)=dot_product(b,c)-1/2.d0
    endif

    if (p>=3) then
      ! order 3 conditions:
        coneq(3)=dot_product(b,matmul(A,c))-1/6.d0
        coneq(4)=dot_product(b,c**2)-1/3.d0
    endif

    if (p>=4) then
      ! order 4 conditions:
        coneq(5)=dot_product(b,matmul(A,c**2))-1/12.d0
        coneq(6)=dot_product(b,matmul(A,matmul(A,c)))-1/24.d0
        coneq(7)=dot_product(b,matmul(A,c)*c)-1/8.d0
        coneq(8)=dot_product(b,c**3)-1/4.d0
    endif

    if (p>=5) then
      ! order 5 conditions:
        coneq(9)=dot_product(b,matmul(A,c**3))-1/20.d0
        coneq(10)=dot_product(b,matmul(A,matmul(A,c**2)))-1/60.d0
        coneq(11)=dot_product(b,matmul(A,matmul(A,matmul(A,c))))-1/120.d0
        coneq(12)=dot_product(b,matmul(A,c*matmul(A,c)))-1/40.d0
        coneq(13)=dot_product(b,matmul(A,c**2)*c)-1/15.d0
        coneq(14)=dot_product(b,matmul(A,matmul(A,c))*c)-1/30.d0
        coneq(15)=dot_product(b,matmul(A,c)*c**2)-1/10.d0
        coneq(16)=dot_product(b,matmul(A,c)*matmul(A,c))-1/20.d0
        coneq(17)=dot_product(b,c**4)-1/5.d0
    endif

    if (p>=6) then
      ! order 6 conditions:
        coneq(18)=dot_product(b,matmul(A,c**4))-1/30.d0
        coneq(19)=dot_product(b,matmul(A,matmul(A,c**3)))-1/120.d0
        coneq(20)=dot_product(b,matmul(A,matmul(A,matmul(A,c**2))))-1/360.d0
        coneq(21)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,c)))))-1/720.d0
        coneq(22)=dot_product(b,matmul(A,matmul(A,c*matmul(A,c))))-1/240.d0
        coneq(23)=dot_product(b,matmul(A,c*matmul(A,c**2)))-1/90.d0
        coneq(24)=dot_product(b,matmul(A,c*matmul(A,matmul(A,c))))-1/180.d0
        coneq(25)=dot_product(b,matmul(A,c**2*matmul(A,c)))-1/60.d0
        coneq(26)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,c)))-1/120.d0
        coneq(27)=dot_product(b,matmul(A,c**3)*c)-1/24.d0
        coneq(28)=dot_product(b,matmul(A,matmul(A,c**2))*c)-1/72.d0
        coneq(29)=dot_product(b,matmul(A,matmul(A,matmul(A,c)))*c)-1/144.d0
        coneq(30)=dot_product(b,matmul(A,c*matmul(A,c))*c)-1/48.d0
        coneq(31)=dot_product(b,matmul(A,c**2)*c**2)-1/18.d0
        coneq(32)=dot_product(b,matmul(A,matmul(A,c))*c**2)-1/36.d0
        coneq(33)=dot_product(b,matmul(A,c)*c**3)-1/12.d0
        coneq(34)=dot_product(b,matmul(A,c**2)*matmul(A,c))-1/36.d0
        coneq(35)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,c))-1/72.d0
        coneq(36)=dot_product(b,matmul(A,c)*matmul(A,c)*c)-1/24.d0
        coneq(37)=dot_product(b,c**5)-1/6.d0
    endif

    if (p>=7) then
      ! order 7 conditions:
        coneq(38)=dot_product(b,matmul(A,c**5))-1/42.d0
        coneq(39)=dot_product(b,matmul(A,matmul(A,c**4)))-1/210.d0
        coneq(40)=dot_product(b,matmul(A,matmul(A,matmul(A,c**3))))-1/840.d0
        coneq(41)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,c**2)))))-1/2520.d0
        coneq(42)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,matmul(A,c))))))-1/5040.d0
        coneq(43)=dot_product(b,matmul(A,matmul(A,matmul(A,c*matmul(A,c)))))-1/1680.d0
        coneq(44)=dot_product(b,matmul(A,matmul(A,c*matmul(A,c**2))))-1/630.d0
        coneq(45)=dot_product(b,matmul(A,matmul(A,c*matmul(A,matmul(A,c)))))-1/1260.d0
        coneq(46)=dot_product(b,matmul(A,matmul(A,c**2*matmul(A,c))))-1/420.d0
        coneq(47)=dot_product(b,matmul(A,matmul(A,matmul(A,c)*matmul(A,c))))-1/840.d0
        coneq(48)=dot_product(b,matmul(A,c*matmul(A,c**3)))-1/168.d0
        coneq(49)=dot_product(b,matmul(A,c*matmul(A,matmul(A,c**2))))-1/504.d0
        coneq(50)=dot_product(b,matmul(A,c*matmul(A,matmul(A,matmul(A,c)))))-1/1008.d0
        coneq(51)=dot_product(b,matmul(A,c*matmul(A,c*matmul(A,c))))-1/336.d0
        coneq(52)=dot_product(b,matmul(A,c**2*matmul(A,c**2)))-1/126.d0
        coneq(53)=dot_product(b,matmul(A,c**2*matmul(A,matmul(A,c))))-1/252.d0
        coneq(54)=dot_product(b,matmul(A,c**3*matmul(A,c)))-1/84.d0
        coneq(55)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,c**2)))-1/252.d0
        coneq(56)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,matmul(A,c))))-1/504.d0
        coneq(57)=dot_product(b,matmul(A,c*matmul(A,c)*matmul(A,c)))-1/168.d0
        coneq(58)=dot_product(b,matmul(A,c**4)*c)-1/35.d0
        coneq(59)=dot_product(b,matmul(A,matmul(A,c**3))*c)-1/140.d0
        coneq(60)=dot_product(b,matmul(A,matmul(A,matmul(A,c**2)))*c)-1/420.d0
        coneq(61)=dot_product(b,matmul(A,matmul(A,matmul(A,matmul(A,c))))*c)-1/840.d0
        coneq(62)=dot_product(b,matmul(A,matmul(A,c*matmul(A,c)))*c)-1/280.d0
        coneq(63)=dot_product(b,matmul(A,c*matmul(A,c**2))*c)-1/105.d0
        coneq(64)=dot_product(b,matmul(A,c*matmul(A,matmul(A,c)))*c)-1/210.d0
        coneq(65)=dot_product(b,matmul(A,c**2*matmul(A,c))*c)-1/70.d0
        coneq(66)=dot_product(b,matmul(A,matmul(A,c)*matmul(A,c))*c)-1/140.d0
        coneq(67)=dot_product(b,matmul(A,c**3)*c**2)-1/28.d0
        coneq(68)=dot_product(b,matmul(A,matmul(A,c**2))*c**2)-1/84.d0
        coneq(69)=dot_product(b,matmul(A,matmul(A,matmul(A,c)))*c**2)-1/168.d0
        coneq(70)=dot_product(b,matmul(A,c*matmul(A,c))*c**2)-1/56.d0
        coneq(71)=dot_product(b,matmul(A,c**2)*c**3)-1/21.d0
        coneq(72)=dot_product(b,matmul(A,matmul(A,c))*c**3)-1/42.d0
        coneq(73)=dot_product(b,matmul(A,c)*c**4)-1/14.d0
        coneq(74)=dot_product(b,matmul(A,c**2)*matmul(A,c**2))-1/63.d0
        coneq(75)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,c**2))-1/126.d0
        coneq(76)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,matmul(A,c)))-1/252.d0
        coneq(77)=dot_product(b,matmul(A,c**3)*matmul(A,c))-1/56.d0
        coneq(78)=dot_product(b,matmul(A,matmul(A,c**2))*matmul(A,c))-1/168.d0
        coneq(79)=dot_product(b,matmul(A,matmul(A,matmul(A,c)))*matmul(A,c))-1/336.d0
        coneq(80)=dot_product(b,matmul(A,c*matmul(A,c))*matmul(A,c))-1/112.d0
        coneq(81)=dot_product(b,matmul(A,c**2)*matmul(A,c)*c)-1/42.d0
        coneq(82)=dot_product(b,matmul(A,matmul(A,c))*matmul(A,c)*c)-1/84.d0
        coneq(83)=dot_product(b,matmul(A,c)*matmul(A,c)*c**2)-1/28.d0
        coneq(84)=dot_product(b,matmul(A,c)*matmul(A,c)*matmul(A,c))-1/56.d0
        coneq(85)=dot_product(b,c**6)-1/7.d0
    endif

return
end subroutine

