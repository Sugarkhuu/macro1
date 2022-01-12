% --------------------------------------------
% RBC_model.m
% sets up the RBC model,decentralized version.
% Compute derivatives symbolically.
% --------------------------------------------


% -----------------------------------
% define symbolic variables
% -----------------------------------
syms alphapar betapar deltapar psipar rhoapar sigmaapar chipar;    % parameters
syms at ct ht invt kt rkt wt yt;                                   % variables: today 
syms atp ctp htp invtp ktp rktp wtp ytp;                            % variables: tomorrow 

% -----------------------------------
% model equations; create function f
% -----------------------------------
jkl      = 0;
% - labor supply
jkl      = jkl+1; f(jkl,1) = wt - chipar*ht^psipar*ct;
% - consumption Euler equation 
jkl      = jkl+1; f(jkl,1) = 1 - ( betapar*ct/ctp*(rktp + (1-deltapar)) );
% - labor demand
jkl      = jkl+1; f(jkl,1) = wt -  (1-alphapar)*at*kt^alphapar*(ht)^(-alphapar);
% - capital demand
jkl      = jkl+1; f(jkl,1) = rkt - alphapar*at*kt^(alphapar-1)*(ht)^(1-alphapar);
% - goods-market clearing
jkl      = jkl+1; f(jkl,1) = ct + invt - yt;
% - output
jkl      = jkl+1; f(jkl,1) = yt - (at*kt^alphapar*ht^(1-alphapar));
% - capital accumulation
jkl      = jkl+1; f(jkl,1) = ktp - ( (1-deltapar)*kt + invt );
% - shock
jkl      = jkl+1; f(jkl,1) = log(atp) - rhoapar*log(at);

% -----------------------------------
% Define the vector of controls, y, and states, x; and names_
% -----------------------------------
x = [kt at];
y = [ct ht invt rkt wt yt];
xp = [ktp atp];
yp = [ctp htp invtp rktp wtp ytp];

x_ = strvcat('kt','at');
y_ = strvcat('ct','ht','invt','rkt','wt','yt');

paramsym = [alphapar betapar deltapar psipar rhoapar sigmaapar chipar];

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

