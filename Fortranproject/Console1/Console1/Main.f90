    !  Console1.f90
    !
    !  FUNCTIONS:
    !  Console1 - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: Console1
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    program Console1

    use Header
    use routines

    implicit none

    ! Variables
    type (structparamstype) :: params
    type (gridsType) :: grids
    real (kind=rk) :: minInc(Tperiods), maxInc(Tperiods), maxA(Tperiods+1)
    real (kind=rk) :: tol, minCons, r, startA, span, loggrid(numPointsA)
    integer :: ixT, i, temp(numPointsA)

    r=0.01
    startA =0
    params%beta = 0.970 !                ! Discount factor
    params%gamma = 31!             ! Coefficient of relative risk aversion
    params%mu = 0 !;                      ! mean of initial log income
    params%sigma = 0.0922 !mean([0.007 0.005 0.010 0.012])^0.5;! variance of innovations to log income
    params%rho = 0.96 !                ! persistency of log income
    params%nu = 0.4637 ! mean([0.422,0.417,0.5167,0.499]); !lesuire
    params%hrsWrk = 0.3159 !1840/5824;
    params%delta(3) = 8.2270 !mean([9.083298	8.5256042	7.906901	7.392151]);
    params%delta(2)=0.0767 !mean([0.029333	0.061582208	0.094386	0.121426]);
    params%delta(1)=-7.2503e-04 !mean([-0.00023	-0.0005901	-0.00091	-0.00117]);
    params%spouseInc = 6.2318e+03!mean([5121, 5282, 6840, 7684]);

    tol = 1e-10 !                ! max allowed error
    params%minCons = 1e-5 !              ! min allowed consumption
    
    call getIncomeGrid(params, grids%Ygrid, grids%incTransitionMrx, minInc, maxInc, grids%AIMEgrid, grids%benefit)
    
    !Set maximum assets
    maxA(1) = startA;
    do ixt = 2, Tperiods+1
        maxA(ixt) = (maxA(ixt - 1) + maxInc(ixt-1) ) * (1+r)
    end do
    
    !Create asset grid
    do ixt = 1, Tperiods+1
        span =  log(1.0+log(1.0+log(1+maxA(ixt)))) - log(1.0+log(1.0+log(1.0)))/(numPointsA-1)
        loggrid = log(1.0+log(1.0+log(1.0))) + span*(/(i,i=0,numPointsA-1)/) 
        !temp = (/(i,i=1,numPointsA)/)
        !Agrid(ixt, :) = real(temp, rk)
        grids%Agrid(ixt, :) = (/(exp(exp(exp(loggrid(i))-1.0)-1.0)-1.0,i=1,numPointsA)/) !exp(exp(exp(loggrid)-1)-1)-1
        !getGrid(borrowCon(ixt), maxAss(ixt), numPointsA, gridMethod)
    end do
    
     print *, "Hello World!"
    
    end program Console1

