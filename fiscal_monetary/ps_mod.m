% --------------------------------------------
% ps_mod.m
% sets up the fiscal and monetary model,decentralized version.
% Compute derivatives symbolically.
% --------------------------------------------


% -----------------------------------
% define symbolic variables
% -----------------------------------

var_par = {'beta', 'gamma', 'alpha', 'pi_star', 'R_star', 's_star', 'b_star', 'y_star'};

% Variables
% today
var_st  = {'Rt','st'};
var_ct  = {'ct','pit','bt'};
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
jkl      = jkl+1; f(jkl,1) = 1/Rt - 1/R_star - alpha*(1/pit-1/pi_star);
% - Fiscal policy rule
jkl      = jkl+1; f(jkl,1) = stp - s_star - gamma*(1/Rt*bt - b_star/R_star);
% - consumer budget
jkl      = jkl+1; f(jkl,1) = ctp + stp + btp/Rtp - y_star - bt/pitp;
% - market clearing
jkl      = jkl+1; f(jkl,1) = ct - y_star;


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

