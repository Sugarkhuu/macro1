% --------------------------------------------
% RBC_model.m
% sets up the RBC model,decentralized version.
% Compute derivatives symbolically.
% --------------------------------------------


% -----------------------------------
% define symbolic variables
% -----------------------------------
syms zetapar Apar rhopar betapar nupar etapar pipar xipar chipar kappapar sigmaApar taupar;            % parameters
syms ct et wt Psit lambdat stoch_betat yt At lt xt ut ft Upsilont Jt vt mt qt w_firmt subst;                  % variables: today 
syms ctp etp wtp Psitp lambdatp stoch_betatp ytp Atp ltp xtp utp ftp Upsilontp Jtp vtp mtp qtp w_firmtp substp; % variables: tomorrow 

% -----------------------------------
% model equations; create function f
% -----------------------------------
jkl      = 0;
% - budget constraint
jkl      = jkl+1; f(jkl,1) = ct - et*wt*(1-taupar) - Psit; %et*wt*(1-taut) - Psit;
% - marginal utility of consumption
jkl      = jkl+1; f(jkl,1) = lambdat - 1/ct;
% - labor demand
jkl      = jkl+1; f(jkl,1) = yt - At*lt;
% - productivity
jkl      = jkl+1; f(jkl,1) = Atp - Apar - rhopar*(At - Apar);
% - labor good production
jkl      = jkl+1; f(jkl,1) = lt - yt/xt;
% - labor market clearing
jkl      = jkl+1; f(jkl,1) = lt - et;%et
% - unemployment
jkl      = jkl+1; f(jkl,1) = et - 1 + ut;
% - employment dynamic
jkl      = jkl+1; f(jkl,1) = et - (1-nupar)*et - ft*ut; % epar
% - labor firm profit
jkl      = jkl+1; f(jkl,1) = Upsilont - xt + w_firmt;
% - value
jkl      = jkl+1; f(jkl,1) = Jt - Upsilont - (1-nupar)*stoch_betat*Jtp;
% - wage equation
jkl      = jkl+1; f(jkl,1) = wt - etapar*xt - (1-etapar)*pipar;
% - match
jkl      = jkl+1; f(jkl,1) = mt - chipar*(ut^xipar)*(vt^(1-xipar));
% - vacancy fill prob.
jkl      = jkl+1; f(jkl,1) = qt - mt/vt;
% - get employed prob.
jkl      = jkl+1; f(jkl,1) = ft - mt/ut;
% - search profit zero cond.
jkl      = jkl+1; f(jkl,1) = kappapar - qt*stoch_betat*Jtp;
% - stoch_betat
jkl      = jkl+1; f(jkl,1) = stoch_betat - betapar*lambdatp/lambdat;
% - firm profit
jkl      = jkl+1; f(jkl,1) = Psit - Upsilont*et + kappapar*vt;
% % - subsidy
% jkl      = jkl+1; f(jkl,1) = taut - (wt - (etapar*Apar+(1-etapar)*pipar) - (xt-Apar))/wt;
% - wage of firm
jkl      = jkl+1; f(jkl,1) = w_firmt - wt*(1-taupar);
% - subsidy
jkl      = jkl+1; f(jkl,1) = subst - wt + w_firmt;
%
% jkl      = jkl+1; f(jkl,1) = Jt - Jtp;
% jkl      = jkl+1; f(jkl,1) = taut - 0;

% % output and consumption identity 
% jkl      = jkl+1; f(jkl,1) = yt - ct - kappapar*vt;

% - output market clearing - double counting
% jkl      = jkl+1; f(jkl,1) = yt - ct + kapppar*vt;

% -----------------------------------
% Define the vector of controls, y, and states, x; and names_
% -----------------------------------

x = [et At];
y = [ct wt Psit lambdat stoch_betat yt lt xt ut ft Upsilont Jt vt mt qt w_firmt subst];
xp = [etp Atp];
yp = [ctp wtp Psitp lambdatp stoch_betatp ytp ltp xtp utp ftp Upsilontp Jtp vtp mtp qtp w_firmtp substp];

x_ = strvcat('et','At');
y_ = strvcat('ct', 'wt', 'Psit', 'lambdat', 'stoch_betat', 'yt', 'lt', 'xt', 'ut', 'ft', 'Upsilont', 'Jt', 'vt', 'mt', 'qt', 'w_firmt', 'subst');

paramsym = [zetapar Apar rhopar betapar nupar etapar pipar xipar chipar kappapar sigmaApar taupar]; % parameters

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

