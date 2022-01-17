% --------------------------------------------
% RBC_model.m
% sets up the RBC model,decentralized version.
% Compute derivatives symbolically.
% --------------------------------------------


% -----------------------------------
% define symbolic variables
% -----------------------------------
syms zetapar Apar rhopar betapar nupar etapar pipar xipar chipar kappapar sigmaApar;            % parameters
syms ct et wt Psit lambdat stoch_betat yt At lt xt ut ft Upsilont Jt vt mt qt w_firmt;                  % variables: today 
syms ctp etp wtp Psitp lambdatp stoch_betatp ytp Atp ltp xtp utp ftp Upsilontp Jtp vtp mtp qtp w_firmtp; % variables: tomorrow 
syms taut;
syms tautp;

% -----------------------------------
% model equations; create function f
% -----------------------------------
jkl      = 0;
% - budget constraint
jkl      = jkl+1; f(jkl,1) = ct - et*wt - Psit; % consumers pay taxes for the wage subsidy
% - marginal utility of consumption
jkl      = jkl+1; f(jkl,1) = lambdat - 1/ct;
% - output
jkl      = jkl+1; f(jkl,1) = yt - At*lt;
% - productivity process
jkl      = jkl+1; f(jkl,1) = Atp - Apar - rhopar*(At - Apar);
% - labor good production
jkl      = jkl+1; f(jkl,1) = lt - yt/xt;
% - labor market clearing
jkl      = jkl+1; f(jkl,1) = lt - et;
% - unemployment
jkl      = jkl+1; f(jkl,1) = et - 1 + ut;
% - employment dynamic
jkl      = jkl+1; f(jkl,1) = etp - (1-nupar)*et - ft*ut;
% - labor firm profit
jkl      = jkl+1; f(jkl,1) = Upsilont - xt + wt*(1-taut);
% - labor firm value
jkl      = jkl+1; f(jkl,1) = Jt - Upsilont - (1-nupar)*stoch_betat*Jtp;
% - wage equation
jkl      = jkl+1; f(jkl,1) = wt - etapar*xt - (1-etapar)*pipar;
% - match
jkl      = jkl+1; f(jkl,1) = mt - chipar*(ut^xipar)*(vt^(1-xipar));
% - vacancy fill prob.
jkl      = jkl+1; f(jkl,1) = qt - mt/vt;
% - get employed prob.
jkl      = jkl+1; f(jkl,1) = ft - mt/ut;
% - firm zero profit cond.
jkl      = jkl+1; f(jkl,1) = kappapar - qt*stoch_betat*Jtp;
% - stoch_betat
jkl      = jkl+1; f(jkl,1) = stoch_betat - betapar*lambdatp/lambdat;
% - firm profit
jkl      = jkl+1; f(jkl,1) = Psit - Upsilont*et + kappapar*vt;
% - fixed unemployment
jkl      = jkl+1; f(jkl,1) = ut - utp;
% - firm wage
jkl      = jkl+1; f(jkl,1) = w_firmt - wt*(1-taut);


% -----------------------------------
% Define the vector of controls, y, and states, x; and names_
% -----------------------------------

x = [At];
y = [et ct wt Psit lambdat stoch_betat yt lt xt ut ft Upsilont Jt vt mt qt taut w_firmt];
xp = [Atp];
yp = [etp ctp wtp Psitp lambdatp stoch_betatp ytp ltp xtp utp ftp Upsilontp Jtp vtp mtp qtp tautp w_firmtp];

x_ = strvcat('At');
y_ = strvcat('et','ct', 'wt', 'Psit', 'lambdat', 'stoch_betat', 'yt', 'lt', 'xt', 'ut', 'ft', 'Upsilont', 'Jt', 'vt', 'mt', 'qt', 'taut','w_firmt');

paramsym = [zetapar Apar rhopar betapar nupar etapar pipar xipar chipar kappapar sigmaApar]; % parameters

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

