% --------------------------------------------
% RBC_model.m
% sets up the RBC model,decentralized version.
% Compute derivatives symbolically.
% --------------------------------------------


% -----------------------------------
% define symbolic variables
% -----------------------------------
syms zetapar Apar rhopar betapar thetapar etapar pipar xipar chipar kappapar sigmaApar;         % parameters
syms ct et wt Psit lambdat stoch_betat yt At lt xt ut ft Upsilont Jt vt mt qt;                  % variables: today 
syms ctp etp wtp Psitp lambdatp stoch_betatp ytp Atp ltp xtp utp ftp Upsilontp Jtp vtp mtp qtp; % variables: tomorrow 

% -----------------------------------
% model equations; create function f
% -----------------------------------
jkl      = 0;
% - budget constraint
jkl      = jkl+1; f(jkl,1) = ct - et*wt - Psit;
% - marginal utility of consumption
jkl      = jkl+1; f(jkl,1) = lambdat - 1/ct;
% - labor demand
jkl      = jkl+1; f(jkl,1) = yt - At*lt;
% - productivity
jkl      = jkl+1; f(jkl,1) = Atp - Apar - rhopar*(At - Apar);
% - labor good production
jkl      = jkl+1; f(jkl,1) = lt - yt/xt;
% - labor market 
jkl      = jkl+1; f(jkl,1) = lt - et;
% - unemployment
jkl      = jkl+1; f(jkl,1) = et - 1 + ut;
% - employment dynamic
jkl      = jkl+1; f(jkl,1) = etp - (1-thetapar)*et - ft*ut;
% - labor firm profit
jkl      = jkl+1; f(jkl,1) = Upsilont - xt + wt;
% - labor firm profit
jkl      = jkl+1; f(jkl,1) = Jt - Upsilont - (1-thetapar)*stoch_betat*Jtp;
% - wage equation
jkl      = jkl+1; f(jkl,1) = wt - etapar*xt - (1-etapar)*pipar;
% - match
jkl      = jkl+1; f(jkl,1) = mt - chipar*(ut^xipar)*(vt^(1-xipar));
% - vacancy fill
jkl      = jkl+1; f(jkl,1) = qt - mt/vt;
% - get employed
jkl      = jkl+1; f(jkl,1) = ft - mt/ut;
% - get employed
jkl      = jkl+1; f(jkl,1) = kappapar - qt*stoch_betat*Jtp;
% - stoch_betat
jkl      = jkl+1; f(jkl,1) = stoch_betat - betapar*lambdatp/lambdat;
% - stoch_betat
jkl      = jkl+1; f(jkl,1) = Psit - Upsilont*et + kappapar*vt;

% double counting
% % - output market clearing
% jkl      = jkl+1; f(jkl,1) = yt - ct + kapppar*vt;

% -----------------------------------
% Define the vector of controls, y, and states, x; and names_
% -----------------------------------

x = [et At];
y = [ct wt Psit lambdat stoch_betat yt lt xt ut ft Upsilont Jt vt mt qt];
xp = [etp Atp];
yp = [ctp wtp Psitp lambdatp stoch_betatp ytp ltp xtp utp ftp Upsilontp Jtp vtp mtp qtp];

x_ = strvcat('et','At');
y_ = strvcat('ctp', 'wtp', 'Psitp', 'lambdatp', 'stoch_betatp', 'ytp', 'ltp', 'xtp', 'utp', 'ftp', 'Upsilontp', 'Jtp', 'vtp', 'mtp', 'qtp');

paramsym = [zetapar Apar rhopar betapar thetapar etapar pipar xipar chipar kappapar sigmaApar]; % parameters

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

