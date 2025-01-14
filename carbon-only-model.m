% Resource Allocation in Annual Plants MATLAB Script

function [t, S, R, F, U, L] = growthpath(S_Final, R_Final)
  %Terminal Conditions
  p.S_Final = S_Final;  %Final value of shoots
  p.R_Final = R_Final;  %Final value of roots

  % Primary parameters
  p.nu_R = 1; %C:N ratio in roots
  beta = 3;  %Ratio of C:N ratio in roots to C:N ratio in shoots
  p.n = 2^10;  %Number of grid points per unit time
  p.T = 10; %Length of growing season
  p.BGPI_tol = 10^-10;  %Tolerance for BG-PI Hamiltonian condition
  p.IBG_tol = 10^-10; %Tolerance for initial phase to BG controls condition

  % Secondary parameters
  p.C_Final = C(p.S_Final); %Final rate of carbon fixation
  p.N_Final = N(p.R_Final); %Final rate of nitrogen uptake
  p.dCdS_Final = dCdS(p.S_Final); %Final rate of change of carbon fixation rate wrt shoots
  p.dNdR_Final = dNdR(p.R_Final); %Final rate of change of nitrogen uptake rate wrt roots
  p.nu_S = p.nu_R/beta; %C:N ratio in roots
  p.t_star = p.T - 1/p.dCdS_Final; %Compute t_star

  % Create structure w/ terminal conditions for PI solver
    % Note that since we solve forward in z, terminal/initial are in reference to z, not t
  p.PI.args.step = 1;  %Used to compute RK4 step size via h = p.PI.args.step/p.n;
  p.PI.args.z_start = 0; %Starting z1C value for the PI solver
  p.PI.args.z_end = 100;  %End z1C for PI solver, 100 is an arbitrary initial guess
  p.PI.args.L2_init = 0;  %PI initial condition for L2
  p.PI.args.S_init = p.S_Final; %PI initial condition for S
  p.PI.args.t_init = p.t_star;  %PI initial condition for t

  % Penultimate interval solver - first run
  [tP, zP, SP, RP, UP, LP] = PI(p);

  p.PI.BGPI_zwindow_begin_index = find(LP(2,:)<0.99,1,'last');  %Find last index where L2 < 1

  % BG-PI transition point finder - first run
  [tP, zP, SP, RP, UP, LP, p] = BGPI_Find(tP, zP, SP, UP, LP, p);


  % Reorder by time - note that we don't need to reorder R b/c it is constant here
  [~,tPsort] = sort(tP);
  tP = tP(tPsort);
  SP = SP(tPsort);
  LP = LP(:,tPsort);
  UP = UP(:,tPsort);
  zP = zP(tPsort);

  % Refine integration mesh and iterate until BG-PI transition located to within p.BGPI_tol
  while p.BGPI.H > p.BGPI_tol
    [tP2, zP2, SP2, RP2, UP2, LP2] = PI(p);
    p.PI.BGPI_zwindow_begin_index = 1;
    [tP2, zP2, SP2, RP2, UP2, LP2, p] = BGPI_Find(tP2, zP2, SP2, UP2, LP2, p);
  end

  % Append value of variables at BG-PI transition point
  tP = [p.PI.args.t_init, tP];
  SP = [p.PI.args.S_init, SP];
  RP = [p.R_Final RP];
  LP = [[1/g_prime(p.PI.args.z_start);p.PI.args.L2_init], LP];
  zP = [p.PI.args.z_start, zP];
  UP = [[1-UP2(2,end); UP2(2,end); 0; 1; 0], UP];

  clear SP2 RP2 LP2 UP2 tP2 zP2 %Clear unnecessary variables

  % Create struture w/ terminal conditions for balanced growth solver
  p.BG.args.step = 1;
  p.BG.args.t_start = 0;
  p.BG.args.t_end = tP(1);
  p.BG.args.S_end = SP(1);
  p.BG.args.R_end = p.R_Final;
  p.BG.args.L1_end = LP(1,1);
  p.BG.args.L2_end = LP(2,1);
  p.BG.args.z1_end = zP(1);
  [p.BG.args.z2_end,~,~] = fsolve(@(Z)(ZLeft(Z,zP(1),p)),1,optimoptions('fsolve','MaxFunEvals',10000,'Display','off','OptimalityTolerance',1e-20)); %Solve for z2C using necessary conditions


  % Solve for fruits during PI
  FP = FruitsPI(tP, SP, UP(1,:));

  % Final interval
  [tF, SF, RF, FF, UF, LF] = final(FP(end), p);

  % Balanced growth
  [tB, zB, SB, RB, UB, LB] = BG(p);

  % Find beginning of balanced growth
  if(tB(1) < 2/p.n)   %If U < 1 and t = 0 then there is no intial phase
    FB = zeros(1,length(tB));
    t = [tB, tP, tF];
    S = [SB, SP, SF];
    R = [RB, RP, RF];
    F = [FB, FP, FF];
    L = [LB, LP, LF];
    U = [UB, UP, UF];
  else
    [tB, zB, SB, RB, UB, LB, p] = CBG_Find(tB, zB, SB, RB, UB, LB, p);  %Search for beginning of BG - first run
    CBG_ind = 0;
    while p.CBG.Ucondition > p.IBG_tol  %Refine integration mesh and iterate until beginning of BG found to within p.IBG_tol
      [tB2, zB2, SB2, RB2, UB2, LB2] = BG(p);
      CBG_ind = 1;
      [tB2, zB2, SB2, RB2, UB2, LB2, p] = CBG_Find(tB2, zB2, SB2, RB2, UB2, LB2, p);
    end
    if CBG_ind == 0 %Move on if no refinement is possible
      SB2 = nan;
      RB2 = nan;
      LB2 = [nan;nan];
      UB2 = [nan;nan;nan;nan;nan];
      tB2 = nan;
      zB2 = [nan;nan];
    end

    % Append value of variables at beginning of balanced growth
    SB = [SB2, SB];
    RB = [RB2, RB];
    LB = [LB2, LB];
    UB = [UB2, UB];
    tB = [tB2, tB];
    zB = [zB2, zB];
    FB = zeros(1,length(tB));
    clear SB2 RB2 LB2 UB2 tB2 zB2 %Clear unnecessary variables

    % Compute conditions required for S-only or R-only initial phase
    [S_only_condition, R_only_condition] = convergence_conditions(p);

    % Initial phase
    if S_only_condition < R_only_condition
      [tC, SC, RC, FC, UC, LC] = shootonly(p);
    elseif R_only_condition < S_only_condition
      [tC, SC, RC, FC, UC, LC] = rootonly(p);
    else
      error('could not determind C-BG transition, both shoot-only and root-only growth are possible')
    end

    % Combine Vectors From 4 Stages
    t = [tC, tB, tP, tF];
    S = [SC, SB, SP, SF];
    R = [RC, RB, RP, RF];
    F = [FC, FB, FP, FF];
    L = [LC, LB, LP, LF];
    U = [UC, UB, UP, UF];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Functions called in growthpath
% Growth stage functions

% Penultimate interval function
function [t, z, S, R, U, L] = PI(p)
  %RK4 Parameters
  h = p.PI.args.step/p.n; %Step size
  h2 = h/2; %Half step size for RK4
  h6 = h/6; %h/6 for RK4 update

  %Initialize vectors for penultimate interval
  z = p.PI.args.z_start:h:p.PI.args.z_end;  %z1C
  length_z = length(z);
  S = p.PI.args.S_init*ones(1,length_z);
  lambda_1 = 1./g_prime(z);
  lambda_2 = p.PI.args.L2_init*ones(1,length_z);  %L1 given by necessary conditions
  t = p.PI.args.t_init*ones(1,length_z);

  %Solve forward in z1C via RK4
  for j = 1:(length_z-1)

    [k11, k12, k13] = PIRK4(S(j)         , z(j)     , p);
    [k21, k22, k23] = PIRK4(S(j) + h2*k11, z(j) + h2, p);
    [k31, k32, k33] = PIRK4(S(j) + h2*k21, z(j) + h2, p);
    [k41, k42, k43] = PIRK4(S(j) + h*k31 , z(j) + h , p);

    S(j+1) = S(j) + h6*(k11 + 2*(k21 + k31) + k41);
    lambda_2(j+1) = lambda_2(j) + h6*(k12 + 2*(k22 + k32) + k42);
    t(j+1) = t(j) + h6*(k13 + 2*(k23 + k33) + k43);

    u1C = p.nu_S*p.N_Final*z(j+1)/C(S(j+1)); %Update u1c

    % Stop and truncate if controls are unbounded
    if u1C > 1  || u1C < 0
      if j == 1
        warning('No PI - controls unbounded immediately')
        t = nan;
        S = nan;
        R = nan;
        z = nan;
        L = nan;
        U = nan;
        break
      else
        ind = j + 1;
        t = t(1:ind);
        S = S(1:ind);
        R = p.R_Final*ones(1,ind);
        z = z(1:ind);
        L(1,:) = lambda_1(1:ind);
        L(2,:) = lambda_2(1:ind);
        u1C = p.nu_S*p.N_Final*z./C(S);
        U(1,:) = 1 - u1C;  %u0C
        U(2,:) = u1C;  %u1C
        U(3,:) = zeros(1,ind); %u2C
        U(4,:) = ones(1,ind);  %u1N
        U(5,:) = zeros(1,ind); %u2N
        break
      end
    elseif j == length_z - 1
      ind = j + 1;
      t = t(1:ind);
      S = S(1:ind);
      R = p.R_Final*ones(1,ind);
      z = z(1:ind);
      L(1,:) = lambda_1(1:ind);
      L(2,:) = lambda_2(1:ind);
      u1C = p.nu_S*p.N_Final*z./C(S);
      U(1,:) = 1 - u1C;  %u0C
      U(2,:) = u1C;  %u1C
      U(3,:) = zeros(1,ind); %u2C
      U(4,:) = ones(1,ind);  %u1N
      U(5,:) = zeros(1,ind); %u2N
      break
    end
  end
end

% BG-PI transition point finder
function [tP, zP, SP, RP, UP, LP, p] = BGPI_Find(t, z, S, U, L, p)
  left_index = p.PI.BGPI_zwindow_begin_index;
  right_index = length(z);
  [ZL,~,~] = fsolve(@(Z)(ZLeft(Z,z(left_index),p)),1,optimoptions('fsolve','MaxFunEvals',10000,'Display','off','OptimalityTolerance',1e-20));
  [ZR,~,~] = fsolve(@(Z)(ZLeft(Z,z(right_index),p)),1,optimoptions('fsolve','MaxFunEvals',10000,'Display','off','OptimalityTolerance',1e-20));
  HL = L(2,left_index) *(p.nu_R*p.N_Final*g_prime_recip(ZL) + C(S(left_index))*g_prime(ZL));
  HR = L(2,right_index) *(p.nu_R*p.N_Final*g_prime_recip(ZR) + C(S(right_index))*g_prime(ZR));
  HDiffsignL = sign(HL - p.C_Final);
  HDiffsignR = sign(HR - p.C_Final);

  if HDiffsignL*HDiffsignR > 0
    warning('sgn(H - H*) same on both sides of potential BG-PI transition point ')
    stop = 1;
    SP = nan;
    RP = nan;
    LP = nan;
    UP = nan;
    tP = nan;
    zP = nan;
    return
  end

  while abs(right_index - left_index) > 4
    mid_index = floor((right_index + left_index)/2);
    [ZM,~,~] = fsolve(@(Z)(ZLeft(Z,z(mid_index),p)),1,optimoptions('fsolve','MaxFunEvals',10000,'Display','off','OptimalityTolerance',1e-20));
    HM = L(2,mid_index) *(p.nu_R*p.N_Final*g_prime_recip(ZM) + C(S(mid_index))*g_prime(ZM));
    if sign(HM - p.C_Final) == HDiffsignR
      right_index = mid_index;
    else
      left_index = mid_index;
    end
  end
  U1NM = ((C(S(mid_index))./p.N_Final) - p.nu_R*ZM)./(p.nu_S*z(mid_index) - p.nu_R*ZM);
  U1CM = z(mid_index)*p.nu_S*U1NM*p.N_Final/C(S(mid_index));
  if U1NM > 1 || U1NM < 0 || U1CM > 1 || U1CM < 0
    warning('Conrols Not Bounded At Potential BG-PI Boundary')
  end

  % Update variables and arguments for PI solver
  tP = t(1:left_index);
  SP = S(1:left_index);
  RP = p.R_Final*ones(1,left_index);
  LP = L(:,1:left_index);
  UP = U(:,1:left_index);
  zP = z(1:left_index);

  p.PI.args.step = p.PI.args.step/100; %Reduce step size for next iteration
  p.PI.args.z_start = z(left_index);
  p.PI.args.z_end = z(right_index);
  p.PI.args.L2_init = L(2,left_index);
  p.PI.args.S_init = S(left_index);
  p.PI.args.t_init = t(left_index);
  p.BGPI.H = abs(HM - p.C_Final); %Hamiltonian condition
end

% Fruits during the penultimate interval
function F = FruitsPI(t,S,u0C)
  F = zeros(1,length(t));
  hF = diff(t); %Account for non-uniform time stepping
  hF6 = hF/6; %For RK4

  % Solve forwared in time via RK4, using midpoints for half time-steps when necessary
  for i = 1:length(t)-1
    k1 = (2/3)*u0C(i)*C(S(i));
    k2 = (2/3)*0.5*(u0C(i) + u0C(i+1))*(C(0.5*(S(i) + S(i+1))));
    k4 = (2/3)*u0C(i+1)*C(S(i+1));
    F(i+1) = F(i) + hF6(i)*(k1 + 4*k2 + k4);
  end
end

% Final interval function
function [t, S, R, F, U, L] = final(F_star,p)
  h = 1/p.n;
  t = p.t_star:h:p.T;
  length_t = length(t);
  F = (2/3)*C(p.S_Final)*(t - p.t_star) + F_star;
  S = p.S_Final*ones(1,length_t);
  R = p.R_Final*ones(1,length_t);
  L(1,:) = dCdS(p.S_Final)*(p.T-t);
  L(2,:) = zeros(1,length_t);
  U = zeros(5,length_t);
  U(1,:) = ones(1,length_t);
  U(4,:) = ones(1,length_t);
end

% Balanced growth function
function [t, Z, S, R, U, L] = BG(p)
  %Parameters
  h = p.BG.args.step/p.n; %Step size
  h2 = h/2; %Half step size for RK4
  h6 = h/6; %h/6 for RK4 update

  %  Verify that there is a balanced growth phase
  u1N_hat = (C(p.BG.args.S_end) - (p.nu_R*N(p.BG.args.R_end)*p.BG.args.z2_end))/(N(p.BG.args.R_end)*(p.nu_S*p.BG.args.z1_end - p.nu_R*p.BG.args.z2_end)); %Compute u1N at the end of BG
  u1C_hat = (p.nu_S*N(p.BG.args.R_end)*p.BG.args.z1_end*u1N_hat)/C(p.BG.args.S_end); %Compute u1C at the end of BG

  if u1N_hat > 1 || u1N_hat < 0 || u1C_hat > 1 || u1C_hat < 0  %Check that controls are bounded
    S = p.BG.args.S_end;
    R = p.BG.args.R_end;
    L = [p.BG.args.L1_end; p.BG.args.L2_end];
    U = [0; 1; 0; 1; 0];
    t = p.BG.args.t_end;
    Z = [p.BG.args.z1_end; p.BG.args.z2_end];
    warning('No BG - controls unbounded immediately')
    return
  end

  %Initialize vectors for balanced growth
  t = p.BG.args.t_start:h:p.BG.args.t_end;
  length_t = length(t);
  S = p.BG.args.S_end*ones(1,length_t);
  R = p.BG.args.R_end*ones(1,length_t);
  lambda_1 = p.BG.args.L1_end*ones(1,length_t);
  lambda_2 = p.BG.args.L2_end*ones(1,length_t);
  z1 = p.BG.args.z1_end*ones(1,length_t);
  z2 = p.BG.args.z2_end*ones(1,length_t);

  % Solve backwards using RK4
  for i = 1:length_t-1
    j = length_t + 1 - i;

    [k11, k12, k13, k14, k15, k16] = BGRK4(S(j)         , R(j)         , lambda_1(j)         , z1(j)         , z2(j)         , p);
    [k21, k22, k23, k24, k25, k26] = BGRK4(S(j) - h2*k11, R(j) - h2*k12, lambda_1(j) - h2*k13, z1(j) - h2*k15, z2(j) - h2*k16, p);
    [k31, k32, k33, k34, k35, k36] = BGRK4(S(j) - h2*k21, R(j) - h2*k22, lambda_1(j) - h2*k23, z1(j) - h2*k25, z2(j) - h2*k26, p);
    [k41, k42, k43, k44, k45, k46] = BGRK4(S(j) - h*k31 , R(j) - h*k32 , lambda_1(j) - h*k33 , z1(j) - h*k35 , z2(j) - h*k36 , p);

    S(j-1) = S(j) - h6*(k11 + 2*(k21 + k31) + k41);
    R(j-1) = R(j) - h6*(k12 + 2*(k22 + k32) + k42);
    lambda_1(j-1) = lambda_1(j) - h6*(k13 + 2*(k23 + k33) + k43);
    lambda_2(j-1) = lambda_2(j) - h6*(k14 + 2*(k24 + k34) + k44);
    z1(j-1) = z1(j) - h6*(k15 + 2*(k25 + k35) + k45);
    z2(j-1) = z2(j) - h6*(k16 + 2*(k26 + k36) + k46);

 % Compute controls at current time step
    u1N = (C(S(j-1)) - (p.nu_R*N(R(j-1)))*z2(j-1))/(N(R(j-1))*(p.nu_S*z1(j-1) - p.nu_R*z2(j-1)));
    u1C = (p.nu_S*N(R(j-1))*z1(j-1)*u1N)/C(S(j-1));

    % Stop and update variables if controls become unbounded or if t = 0 reached
    if u1N > 1  || u1N < 0 || u1C > 1  || u1C < 0 || j==2
      balgrowthindex = j-1;
      u1N = (C(S(balgrowthindex:end)) - (p.nu_R*N(R(balgrowthindex:end)).*z2(balgrowthindex:end)))./(N(R(balgrowthindex:end)).*(p.nu_S*z1(balgrowthindex:end) - p.nu_R*z2(balgrowthindex:end)));
      u1C = (p.nu_S*N(R(balgrowthindex:end)).*z1(balgrowthindex:end).*u1N)./C(S(balgrowthindex:end));
      S = S(balgrowthindex:end);
      R = R(balgrowthindex:end);
      L(1,:) = lambda_1(balgrowthindex:end);
      L(2,:) = lambda_2(balgrowthindex:end);
      U(1,:) = zeros(1,length(S)); %u0C
      U(2,:) = u1C;  %u1C
      U(3,:) = 1 - U(2,:);  %u2C
      U(4,:) = u1N;  %u1N
      U(5,:) = 1 - U(4,:);  %u2N
      Z(1,:) = z1(balgrowthindex:end);
      Z(2,:) = z2(balgrowthindex:end);
      t = t(balgrowthindex:end);
      break
    end
  end
end

% Balanced Growth Constraints on Z - used for computing left limit of z2C from BG
function X = ZLeft(Z2,Z1,p)
  X = (g_prime(Z2)/g_prime_recip(Z2)) - (p.nu_R/p.nu_S)*(g_prime(Z1)/g_prime_recip(Z1));
end

% Convergence stage to balanced growth phase transition point finder
function [t, z, S, R, U, L, p] = CBG_Find(t, z, S, R, U, L, p)
  % Check to see if transition is between first two time steps
  if sign(U(2,1) - U(4,1))*sign(U(2,2) - U(4,2)) < 0
    left_index = 1;
    right_index = 2;
  else  %If not, do full binary search over entire BG interval
    left_index = 1;
    right_index = length(t);
    ULsign = sign(U(2,left_index) - U(4,left_index));
    URsign = sign(U(2,right_index) - U(4,right_index));

    if ULsign*URsign > 0
      error('sgn(u1c - u1n) same on both sides of potential C-BG transition point')
    end

    while abs(right_index - left_index) > 4
      mid_index = floor((right_index + left_index)/2);
      UMsign = sign(U(2,mid_index) - U(4,mid_index));
      if UMsign == URsign
        right_index = mid_index;
      else
        left_index = mid_index;
      end
    end
  end

  % Update variables and terminal conditions for next refinement
  p.BG.args.t_start = t(left_index);
  S = S(right_index:end);
  R = R(right_index:end);
  L = L(:,right_index:end);
  U = U(:,right_index:end);
  t = t(right_index:end);
  z = z(:,right_index:end);
  p.BG.args.step = p.BG.args.step/1000;
  p.BG.args.t_end = t(1);
  p.BG.args.S_end = S(1);
  p.BG.args.R_end = R(1);
  p.BG.args.L1_end = L(1,1);
  p.BG.args.L2_end = L(2,1);
  p.BG.args.z1_end = z(1,1);
  p.BG.args.z2_end = z(2,1);
  p.CBG.Ucondition = abs(U(2,1) - U(4,1));
end

% Function to compute conditions for each type of initial phase
function [S_only_condition, R_only_condition] = convergence_conditions(p)
  S_only_condition = abs(p.nu_S*p.BG.args.L1_end*N(p.BG.args.R_end)*g(C(p.BG.args.S_end)/(p.nu_S*N(p.BG.args.R_end))) - p.C_Final);
  R_only_condition = abs(p.nu_R*p.BG.args.L2_end*N(p.BG.args.R_end)*g(C(p.BG.args.S_end)/(p.nu_R*N(p.BG.args.R_end))) - p.C_Final);
end

% Initial (convergence) stage
% Shoot-only growth
function [t, S, R, F, U, L] = shootonly(p)
  h = 1/p.n;
  h2 = h/2;
  h6 = h/6;

  % Define constants
  p.S.N_end = N(p.BG.args.R_end); %Initial nitrogen
  p.S.dNdR_end = dNdR(p.BG.args.R_end); %Initial nitrogen derivative

  %Initialize vectors for shoot-only growth
  t = 0:h:p.BG.args.t_end;
  length_t = length(t);
  S = p.BG.args.S_end*ones(1,length_t);
  R = p.BG.args.R_end*ones(1,length_t);
  lambda_1 = p.BG.args.L1_end*ones(1,length_t);
  lambda_2 = p.BG.args.L2_end*ones(1,length_t);
  F = zeros(1,length_t);

  % Solve backwards in time using RK4
  for i = 1:length_t-1
    j = length_t + 1 - i;

    [k11, k12, k13] = ShootonlyRK4(S(j)         ,  lambda_1(j)         , p);
    [k21, k22, k23] = ShootonlyRK4(S(j) - h2*k11,  lambda_1(j) - h2*k12, p);
    [k31, k32, k33] = ShootonlyRK4(S(j) - h2*k21,  lambda_1(j) - h2*k22, p);
    [k41, k42, k43] = ShootonlyRK4(S(j) - h*k31 ,  lambda_1(j) - h*k32 , p);

    S(j-1) = S(j) - h6*(k11 + 2*(k21 + k31) + k41);
    lambda_1(j-1) = lambda_1(j) - h6*(k12 + 2*(k22 + k32) + k42);
    lambda_2(j-1) = lambda_2(j) - h6*(k13 + 2*(k23 + k33) + k43);
  end
  U = zeros(5,length_t);
  U(2,:) = ones(1,length_t);
  U(4,:) = ones(1,length_t);
  L = [lambda_1; lambda_2];
end

% Root-only growth
function [t, S, R, F, U, L] = rootonly(p)
  h = 1/p.n;
  h2 = h/2;
  h6 = h/6;

  % Define constants
  p.R.C_end = C(p.BG.args.S_end); %Initial carbon
  p.R.dCdS_end = dCdS(p.BG.args.S_end); %Initial carbon derivative

  %Initialize vectors for root-only growth
  t = 0:h:p.BG.args.t_end;
  length_t = length(t);
  S = p.BG.args.S_end*ones(1,length_t);
  R = p.BG.args.R_end*ones(1,length_t);
  lambda_1 = p.BG.args.L1_end*ones(1,length_t);
  lambda_2 = p.BG.args.L2_end*ones(1,length_t);
  F = zeros(1,length_t);

  % Solve backwards in time using RK4
  for i = 1:length_t-1
    j = length_t + 1 - i;

    [k11, k12, k13] = RootonlyRK4(R(j)         ,  lambda_2(j)         , p);
    [k21, k22, k23] = RootonlyRK4(R(j) - h2*k11,  lambda_2(j) - h2*k13, p);
    [k31, k32, k33] = RootonlyRK4(R(j) - h2*k21,  lambda_2(j) - h2*k23, p);
    [k41, k42, k43] = RootonlyRK4(R(j) - h*k31 ,  lambda_2(j) - h*k33 , p);

    R(j-1) = R(j) - h6*(k11 + 2*(k21 + k31) + k41);
    lambda_1(j-1) = lambda_1(j) - h6*(k12 + 2*(k22 + k32) + k42);
    lambda_2(j-1) = lambda_2(j) - h6*(k13 + 2*(k23 + k33) + k43);
  end
  U = zeros(5,length_t);
  U(3,:) = ones(1,length_t);
  U(5,:) = ones(1,length_t);
  L = [lambda_1; lambda_2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RK4 Functions

% Penultimate interval RK4 function
function [k1, k2, k3] = PIRK4(S, z, p)
  k1 = (p.nu_S*p.N_Final*g(z)*g_2prime(z))/(dCdS(S)*(g_prime(z))^2); %S
  k2 = -(p.dNdR_Final*p.nu_S*g_prime_recip(z)*g_2prime(z))/(dCdS(S)*(g_prime(z))^3); %L2
  k3 = (g_2prime(z))/(dCdS(S)*(g_prime(z))^2);  %t
end

% Balanced growth RK4 function
function [k1, k2, k3, k4, k5, k6] = BGRK4(S, R, L1, z1, z2, p)
  u = ((C(S)/N(R)) - p.nu_R*z2)/(p.nu_S*z1 - p.nu_R*z2);  %u_1N

  k1 = p.nu_S*u*N(R)*g(z1);  %S
  k2 = p.nu_R*(1 - u)*N(R)*g(z2);  %R
  k3 = -dCdS(S)*L1*g_prime(z1);  %l1
  k4 = -dNdR(R)*L1*p.nu_S*g_prime_recip(z1);  %l2
  k5 = (dNdR(R)*p.nu_S*p.nu_R*g(z2)*g_prime_recip(z1) - dCdS(S)*g_prime(z1)*(p.nu_S*g_prime_recip(z1) + z2*p.nu_R*g_prime(z1)))/(g_2prime(z1)*(p.nu_S*z1 - p.nu_R*z2)); %z1C
  k6 = (dNdR(R)*p.nu_S*g_prime_recip(z1)*(g_prime(z2))^2 - dCdS(S)*(g_prime(z1))^2*g_prime(z2) + g_2prime(z1)*g_prime(z2)*k5)/(g_prime(z1)*g_2prime(z2)); %z2C
end

% Shoots-only growth RK4 function
function [k1, k2, k3] = ShootonlyRK4(S, L1, p)
  z = C(S)/(p.nu_S*p.S.N_end);  %z1C
  k1 = p.nu_S*p.S.N_end*g(z); %S
  k2 = -dCdS(S)*L1*g_prime(z);  %L1
  k3 = -p.S.dNdR_end*p.nu_S*L1*g_prime_recip(z);  %L2
end

% Root-only Growth RK4 Function
function [k1, k2, k3] = RootonlyRK4(R, L2, p)
  z = p.R.C_end/(p.nu_R*N(R));  %z2C
  k1 = p.nu_R*N(R)*g(z);  %R
  k2 = -p.R.dCdS_end*L2*g_prime(z); %L1
  k3 = -dNdR(R)*p.nu_R*L2*g_prime_recip(z); %L2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Functions - G/C/N and Derivatives

% Carbon functions
% Note: Modify both C and dCdS together
function Carbon = C(S) %Rate of carbon fixation for a given amount of shoot biomass
 a = 1; %Proportionality constant for carbon production
 Carbon = a*S; %Rate of carbon fixation is proportional to shoot biomass in carbon
end

function Cs = dCdS(S) %Derivative of carbon wrt shoot
 Cs = 1;
end

% Nitrogen functions
% Note: Modify both N and dNdR together
function Nitrogen = N(R) %Rate of nitrogen assimilation for a given amount of root biomass
 b = 1; %Proportionality constant
 Nitrogen = b*R; %Rate of nitrogen assimilation is proportional to root biomass in carbon
end

function Nr = dNdR(R) %Derivative of nitrogen wrt root
 Nr = 1;
end

% G functions and derivatives

% G(z)
function G = g(z)
 G = z.*(1+z)./(1+z+z.^2);
end

% G'(z)
function gprime = g_prime(z)
 gprime = (1+2*z)./(1+z+z.^2).^2;
end

% G''(z)
function g2prime = g_2prime(z)
  g2prime = -(6*z*(z+1))./(z.^2 + z + 1).^3;
end

% G'(1/z)
function x = g_prime_recip(z)
  x = ((z+2).*(z.^3))./(1+z+z.^2).^2;
end

% d/dz(G'(1/z))
function x = g_prime_recip_prime(z)
  x = (6*(z.^2).*(z+1))./(1+z+z.^2).^3;
end
