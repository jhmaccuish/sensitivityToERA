    module routines
    use Header
    use routines_generic

    implicit none

    contains

    subroutine getIncomeGrid(params, Ygrid, incTransitionMrx, minInc, maxInc, AIMEgrid, benefit)
    implicit none

    !inputs
    type (structparamstype), intent(inout) :: params

    !outputs
    real (kind=rk) :: Ygrid(:,:), incTransitionMrx(:,:), minInc(:), maxInc(:), AIMEgrid(:,:), benefit(:)

    !local
    real (kind=rk) :: sig_inc, ly(numPointsY), Q(numPointsY,numPointsY), upper(Tperiods+1), a
    integer :: t, i

    sig_inc = params%sigma/((1-params%rho**2)**0.5)

    !Why 3
    call tauchen(numPointsY,params%mu,params%rho,params%sigma,3,ly,Q)

    !Ygrid = exp(repmat(ly', T, 1)+repmat(polyval(delta,1:T)'-fc(1:T),1,numPointsY));
    !minInc = exp(repmat((-normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((-normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));
    !maxInc = exp(repmat((normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));

    params%pension = 107.45*52
    upper(1) = 0
    do t=1 , Tperiods
        Ygrid(t,:)= exp(ly+params%delta(1)*t**2+params%delta(2)*t+params%delta(3))
        minInc(t) = exp((-normBnd * sig_inc)+params%delta(1)*t**2+params%delta(2)*t+params%delta(3))
        maxInc(t) = exp((normBnd * sig_inc)+params%delta(1)*t**2+params%delta(2)*t+params%delta(3))
        upper(t+1) = upper(t) + maxInc(t)
        if (t <=45) then
            a = (upper(t+1)/t)/(numAIME-1)
            AIMEgrid(t,:) = a*(/(i,i=0,numAIME-1)/) !linspace(0,upper(t+1)/t,numAIME);
        else
            AIMEgrid(t,:) = AIMEgrid(t-1,:)
        end if
        if (t <=5) then
            benefit(t) = 57.90*52
        else
            benefit(t) = 73.10*52
        end if
    end do
    AIMEgrid(Tperiods+1,:) = AIMEgrid(Tperiods,:)

    !benefit = [57.90*52*ones(5,1);73.10*52*ones(length(minInc)-5,1)]

    end subroutine

    subroutine tauchen(N,mu,rho,sigma,m,Z,Zprob)

    implicit none
    !Function TAUCHEN
    !
    !Purpose:    Finds a Markov chain whose sample paths
    !            approximate those of the AR(1) process
    !                z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
    !            where eps are normal with stddev sigma
    !
    !Format:     {Z, Zprob} = Tauchen(N,mu,rho,sigma,m)
    !
    !Input:      N       scalar, number of nodes for Z
    !            mu      scalar, unconditional mean of process
    !            rho     scalar
    !            sigma   scalar, std. dev. of epsilons
    !            m       max +- std. devs.
    !
    !Output:     Z       N*1 vector, nodes for Z
    !            Zprob   N*N matrix, transition probabilities
    !
    !    Martin Floden
    !    Fall 1996
    !
    !    This procedure is an implementation of George Tauchen's algorithm
    !    described in Ec. Letters 20 (1986) 177-181.
    !

    !inputs
    integer, intent(in) :: N, m
    real (kind=rk), intent(in):: mu, rho, sigma

    !outputs
    real (kind=rk), intent(out) :: Z(N), Zprob(N,N)

    !local
    real (kind=rk) :: a, zstep
    integer :: i, j, k

    a     = (1-rho)*mu;

    Z(N)  = m * sqrt(sigma**2 / (1 - rho**2))
    Z(1)  = -Z(N)
    zstep = (Z(N) - Z(1)) / (N - 1)

    do i = 2, (N-1)
        Z(i) = Z(1) + zstep * (i - 1)
    end do

    Z = Z + a / (1-rho);

    do j = 1, N
        do k = 1, N
            if (k == 1) then
                Zprob(j,k) = cdf_normal((Z(1) - a - rho * Z(j) + zstep / 2) / sigma)
            elseif (k == N) then
                Zprob(j,k) = 1 - cdf_normal((Z(N) - a - rho * Z(j) - zstep / 2) / sigma)
            else
                Zprob(j,k) = cdf_normal((Z(k) - a - rho * Z(j) + zstep / 2) / sigma) - &
                    cdf_normal((Z(k) - a - rho * Z(j) - zstep / 2) / sigma);
            end if
        end do
    end do

    end subroutine

    function cdf_normal(x)
    REAL(KIND=rk), INTENT(in)  :: X
    REAL(KIND=rk) 		:: cdf_normal

    cdf_normal =  0.5d0*( 1+erf(x/SQRT(2.d0)) )

    end function

    subroutine solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, EdU )
    implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids

    !outputs
    real (kind=rk), intent(out) :: V(Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(out) :: policyA1(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(out) :: policyC(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(out) :: policyL(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(out) :: EV(Tperiods+1, numPointsA, numPointsY, numAIME);
    real (kind=rk), intent(out) :: EdU(Tperiods,   numPointsA, numPointsY, numAIME);

    !local
    integer :: ixt, ixAIME, ixA, ixY, ixL
    real (kind=rk) :: negV, A, Y, lbA1, ubA1, AIME, EV1(numPointsA ,numAIME), Agrid1(numPointsA)
    real (kind=rk) :: AIME1grid(numAIME), policyA1temp, negVtemp, realisedV(numPointsA)

    !Set the terminal value function and expected value function to 0
    EV(Tperiods + 1, :,:,:)  = 0;          ! continuation value at T-1
    V(Tperiods + 1,:,:,:) = 0;

    do ixt=Tperiods,1, -1                               ! Loop from time T-1 to 1
        !Agrid1 = grids%Agrid(ixt + 1, :);               ! The grid on assets tomorrow
        AIME1grid = grids%AIMEgrid(ixt + 1, :);
        do ixAIME = 1, numAIME
            do ixA = 1, numPointsA                   ! points on asset grid
                !Although doesn't recieve income still need to loop around
                !hypothetical income because participation doesn't effect
                !earning potential
                ! STEP 1. solve problem at grid points in assets, income + labour choices
                ! ---------------------------------------------------------
                do ixY = 1, numPointsY               ! points on income grid
                    negV = -log(0.0) !inf
                    do ixL = 0,(numPointsL-1),1           ! points of income choice
                        ! Value of income and information for optimisation
                        A    = grids%Agrid(ixt, ixA);            ! assets today
                        Y    = ixL*grids%Ygrid(ixt, ixY)+ (1-ixL)*grids%benefit(ixt);
                        if (ixL==1 .AND.  ixt < 65) then
                            AIME =  grids%Ygrid(ixt, ixY)/ixt + AIME * (ixt-1)/ixt;
                        elseif (ixL==0 .AND. ixt >= 65) then
                            AIME = grids%AIMEgrid(ixt,ixAIME);
                            !Y =   Y + dbPension(AIME);
                        end if
                        !AIME only depends on earned income so add spousal
                        !and pensions after calculating it
                        Y = Y + (ixt >= Tretire)*params%pension+params%spouseInc + (ixt >= 45)*params%pension;
                        lbA1 = grids%Agrid(ixt + 1, 1);          ! lower bound: assets tomorrow
                        ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
                        EV1  = EV(ixt + 1,:, ixY,:);  ! relevant section of EV matrix (in assets tomorrow)

                        ! Compute solution
                        if (ubA1 - lbA1 < params%minCons) then         ! if liquidity constrained
                            !negVtemp = objectivefunc(lbA1, A, Y,ixL,ixt, AIME,EV1);
                            policyA1temp = lbA1;
                        else                               ! if interior solution
                            ![policyA1temp, negVtemp] = ...
                            !    fminbnd(@(A1) objectivefunc(A1, A, Y,ixL,ixt, AIME), lbA1, ubA1, optimset('TolX',tol));
                        end if! if (ubA1 - lbA1 < minCons)
                        if (negVtemp < negV) then
                            negV = negVtemp;
                            policyA1(ixt,ixA,ixY,ixAIME)=policyA1temp;
                            policyL(ixt,ixA,ixY,ixAIME)=ixL;
                        end if
                    end do
                    !end
                    ! Store solution and its value
                    policyC(ixt, ixA, ixY,ixAIME) = A + Y - policyA1(ixt, ixA, ixY,ixAIME)/(1+params%r);
                    V(ixt, ixA, ixY,ixAIME)       = -negV;
                    !dU(ixt, ixA, ixY)      = getmargutility(policyC(ixt, ixA, ixY),policyL(ixt,ixA,ixY));
                end do

                ! STEP 2. integrate out income today conditional on income
                ! yesterday to get EV and EdU
                ! --------------------------------------------------------
                realisedV(:) = V(ixt, ixA, :, ixAIME);
                !realiseddU(:,:) = dU(ixt, ixA, :);
                do ixY = 1,numPointsY,1
                    EV(ixt, ixA, ixY, ixAIME)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV);
                    !EdU(ixt, ixA, ixY) = incTransitionMrx(ixY,:)*realiseddU;
                end do !ixY

            end do!ixA
        end do!ixAIME
        !fprintf('Passed period %d of %d.\n',ixt, T)
    end do!ixt

    end subroutine

    function objectivefunc(params, grids, A1, A0, Y,L,ixP, AIME, EV1)
    implicit none
    !!
    !-------------------------------------------------------------------------------!
    ! This function returns the following quantity:
    ! - (u(c) +  b V( A1))
    ! where c is calculated from today's assets and tomorrow's assets
    !-------------------------------------------------------------------------------!
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    real (kind=rk), intent(in) :: A1, A0, Y, AIME, EV1(:,:)
    integer, intent(in) :: ixP, L
    !ouptut
    real (kind=rk) :: objectivefunc
    !local
    real (kind=rk) :: cons, VA1, mortal(Tperiods+1)
    !! ------------------------------------------------------------------------
    ! Declare global we need this file have access to
    !mortal = (/1,2,&
    !        3,4/)
    mortal = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011071, 0.011907, 0.012807, &
        0.013676, 0.01475, 0.015818, 0.016846, 0.018079, 0.019343, 0.020659, 0.0218, 0.023505, 0.025202, 0.02696, &
        0.028831, 0.031017, 0.033496, 0.036024, 0.038912, 0.042054, 0.045689, 0.049653, 0.054036, 0.05886, 0.064093, &
        0.069636, 0.07533, 0.081069, 0.086912, 0.093067, 0.099807, 0.107483, 0.116125, 0.125196, 0.134361, 0.143881, &
        0.1542, 0.165675, 0.17842, 0.192363, 1.0/);

    !! ------------------------------------------------------------------------
    !Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    !total value (u(c) + b * VA1
    cons = A0  + Y - (A1)/(1+params%r);
    !VA1 = interp1(Agrid1,EV1,A1, interpMethod, 'extrap');
    ![X,Y] = meshgrid(Agrid1,AIME1grid);
    !VA1 = interp2(Agrid1, AIME1grid , EV1, A1, AIME);
    call linearinterp2_withextrap(grids%Agrid(ixP + 1, :), grids%AIMEgrid(ixP + 1, :), numPointsA, numAIME, A1, AIME, VA1, 1, 1, EV1)
    !interp2(EV1, A1/Agrid1(20), AIME/AIME1grid(10));
    objectivefunc = utility(cons,L) + params%beta * (1- mortal(ixP))* VA1;

    !! ------------------------------------------------------------------------
    !The optimisation routine that we will use searches for the minimum of the
    !function. We want the maximum. So we multiply out function here by -1 so
    !that the optimiser will fill the minimum of the negative of our function,
    !i.e. the maximum of our functino

    objectivefunc = - objectivefunc;

    end

    function utility(params,cons,L)
    !This function takes consumption as an argument and returns utility. The
    !utility functin is CRRA except that we add a very small number (eps) to
    !consumption so that the computer can deal wi
    !inputs
    type (structparamstype), intent(in) :: params
    integer, intent(in) :: L
    !real (kind=rk), intent(in) :: A1
    !outpus
    real (kind=rk) :: utility, les

    if (cons<=0) then
        error('Error in utility. Consumption is <=0');
    end
    !10/112 comuniting time -10/112
    les=(L)*(1-params%hrsWrk -10/112)+(~L);
    if params%gamma == 1
    utility = log(cons^params%nu*les**(1-params%nu));
    else
        utility= ((cons**params%nu*les**(1-params%nu))^(1-params%gamma)  )/(1-params%gamma);
    end

    end
    end module routines