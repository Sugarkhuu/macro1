% --------------------------------------------
% ps_mod.m
% sets up the fiscal and monetary model,decentralized version.
% Compute derivatives symbolically.
% --------------------------------------------


% -----------------------------------
% define symbolic variables
% -----------------------------------

var_par = {'beta', 'gamma', 'alpha', 'pi_star', 'R_star', 's_star', 'b_star', 'rho_e','epsilon_par'};

% Variables
% today
var_st  = {'bt','st'};
var_ct  = {'Rtb','Rt','btb','pit'};
% tomorrow 
var_stp = strcat(var_st, 'p');
var_ctp = strcat(var_ct, 'p');

% -----------------------------------
% Define the vector of controls, y, and states, x; and names_
% -----------------------------------

% syms
cellfun(@syms,{var_par, var_st,var_ct,var_stp,var_ctp});

x  = [cellfun(@eval,var_st)];
y  = [cellfun(@eval,var_ct)];
xp = [cellfun(@eval,var_stp)];
yp = [cellfun(@eval,var_ctp)];

x_ = strvcat(string((var_st)));
y_ = strvcat(string((var_ct)));
paramsym = [cellfun(@eval,var_par)]; % parameters

% -----------------------------------
% model equations; create function f
% -----------------------------------
jkl      = 0;
% - Euler equation
jkl      = jkl+1; f(jkl,1) = 1/Rt - beta/pitp;
% - Monetary policy rule
% jkl      = jkl+1; f(jkl,1) = 1/Rt - 1/R_star - alpha*(1/pit-1/pi_star); % - epsilon_Mt
jkl      = jkl+1; f(jkl,1) = 1/Rt - 1/R_star - alpha*(1/pit-1/pi_star); % - epsilon_Mt
% - Fiscal policy rule
jkl      = jkl+1; f(jkl,1) = st - s_star - gamma*(btb/Rtb - b_star/R_star);
% - Budget constraint
jkl      = jkl+1; f(jkl,1) = st + bt/Rt - btb/pit;
% - identity
jkl      = jkl+1; f(jkl,1) = btbp - bt;
% - identity
jkl      = jkl+1; f(jkl,1) = Rtbp - Rt;
% - identity
% jkl      = jkl+1; f(jkl,1) = pitbp - pit;
% - identity
% jkl      = jkl+1; f(jkl,1) = epsilon_Mtp - epsilon_par - rho_e * (epsilon_Mt - epsilon_par);
% - identity
% jkl      = jkl+1; f(jkl,1) = bt - Bt/pt;
% - identity
% jkl      = jkl+1; f(jkl,1) = pitp - ptp/pt;


nx = length(x);
ny = length(y);
% -----------------------------------
% approximate for log-deviations
% -----------------------------------
%f = subs(f, [x,y,xp,yp], (exp([x,y,xp,yp])));
 
% -----------------------------------
% compute symbolic derivatives 
% -----------------------------------
fx  = jacobian(f,x);
fxp = jacobian(f,xp);
fy  = jacobian(f,y);
fyp = jacobian(f,yp);

