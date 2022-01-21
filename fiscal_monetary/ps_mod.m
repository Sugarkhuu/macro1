% --------------------------------------------
% ps_mod.m
% sets up the fiscal and monetary model,decentralized version.
% Compute derivatives symbolically.
% --------------------------------------------


% -----------------------------------
% define symbolic variables
% -----------------------------------
syms beta gamma alpha PIstar Rstar sstar bstar ypar zstar;     % parameters
syms ct PIt Pt Rt st taut zt Bt Btb bt;                  % variables: today 
syms ctp PItp Ptp Rtp stp tautp ztp Btp Btbp btp;         % variables: tomorrow 

% -----------------------------------
% Define the vector of controls, y, and states, x; and names_
% -----------------------------------

x = [Pt Bt];
y = [ct PIt Rt st taut zt Btb bt];
xp = [Ptp Btp];
yp = [ctp PItp Rtp stp tautp ztp Btbp btp];

x_ = strvcat('Pt','Bt');
y_ = strvcat('ct', 'PIt', 'Rt', 'st', 'taut', 'zt', 'Btb', 'bt');

paramsym = [beta gamma alpha PIstar Rstar sstar bstar ypar zstar]; % parameters

% -----------------------------------
% model equations; create function f
% -----------------------------------
jkl      = 0;
% - Euler equation
jkl      = jkl+1; f(jkl,1) = 1/(Pt*Rt) - beta/Ptp;
% - Inflation
jkl      = jkl+1; f(jkl,1) = PItp - Ptp/Pt;
% - Primary surplus
jkl      = jkl+1; f(jkl,1) = st - taut + zt;
% - Monetary policy rule
jkl      = jkl+1; f(jkl,1) = 1/Rt - 1/Rstar - alpha*(1/PIt-1/PIstar);
% - Fiscal policy rule
jkl      = jkl+1; f(jkl,1) = stp - sstar - gamma*(1/Rt*Bt/Pt - bstar/Rstar);
% - 
jkl      = jkl+1; f(jkl,1) = bt - Bt/Pt;
% - consumer budget
jkl      = jkl+1; f(jkl,1) = Pt*ct + Pt*taut + Bt/Rt - Pt*ypar - Pt*zt - Btb;
% - market clearing
jkl      = jkl+1; f(jkl,1) = ct - ypar;
% - Backward p is current period
jkl      = jkl+1; f(jkl,1) = Btbp - Bt;
% 
jkl      = jkl+1; f(jkl,1) = zt - zstar;


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

