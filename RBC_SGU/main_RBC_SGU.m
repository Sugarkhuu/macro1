% ----------------------------------------------------------------------------
% main_RBC_SGU.m
% 
% codes the RBC model in PhD Macro 1.
% Computes second moments and impulse responses using first-order
% perturbation.
% --------------------------------------------------------------------------

% - housekeeping
clear all, close all, clc

% - search path
thispath_ = cd;
addpath([thispath_,'\auxfiles\'] )

% -----------------------------------------
% Section 1. Set parameters and targets
% -----------------------------------------

% - set parameters
deltaparval   = (1.10)^0.25-1; % 10 percent depreciation per year
psiparval     = 2;             % Frisch elasticity of labor supply of .5
rhoaparval    = 0.98;          % persistence TFP shock
sigmaaparval  = 0.007;         % std deviation TFP shock
alphaparval   = 0.33;          % capital share of 1/3
betaparval    = 0.99;          % time discount factor (targets 4 percent real return after depreciation)

% - set targets
hbar       = 1/3;           % st-st target for hours worked

% ----------------------------
% Section 2. Compute steady state and implied params
% ----------------------------
rkbar  = 1/betaparval - (1-deltaparval);            % rental rate of capital (from investment FOC)
abar   = 1;                                         % st-st productivity
ykbar  = rkbar/alphaparval;                         % output-capital ratio (from capital-demand FOC)
khbar  = (ykbar/abar)^(-1/(1-alphaparval));  % capital-labor ratio (from production function)
wbar   = (1-alphaparval)/alphaparval*khbar*rkbar;   % wage (from labor-demand FOC)
ybar   = abar*khbar^alphaparval*hbar;               % output (from production function)
kbar   = 1/ykbar*ybar;                              % from definition of capital-output ratio
invbar = deltaparval*kbar;                          % investment (from capital LOM)
cbar   = ybar - invbar;                             % consumption (from resource constraint)
chiparval = wbar/cbar/hbar^psiparval;               % scaling parameter disutility of work (from labor-supply FOC)
    
xstst = [kbar,abar];
ystst = [cbar, hbar, invbar, rkbar, wbar, ybar];
paramvals = [alphaparval betaparval deltaparval psiparval rhoaparval sigmaaparval chiparval];

% --------------------------
% Section 2. run model file and get derivatives
% --------------------------
RBC_model;

% --------------------------
% Section 3. Check steady state, and evaluate derivatives at steady state
% --------------------------
nf  = eval(subs(f,[x,xp,y,yp,paramsym],[xstst,xstst,ystst,ystst,paramvals]));
if max(abs(nf))>1.e-8
    error('did not solve for steady state')
end
nfx  = eval(subs(fx, [x,xp,y,yp,paramsym],[xstst,xstst,ystst,ystst, paramvals]));
nfxp = eval(subs(fxp,[x,xp,y,yp,paramsym],[xstst,xstst,ystst,ystst, paramvals]));
nfy  = eval(subs(fy, [x,xp,y,yp,paramsym],[xstst,xstst,ystst,ystst, paramvals]));
nfyp = eval(subs(fyp,[x,xp,y,yp,paramsym],[xstst,xstst,ystst,ystst, paramvals]));

% -------------------------
% solve first-order approximation
% -------------------------
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

% -------------------------
% impulse response to TFP shock
% -------------------------
tperiods = 25;
IRF      = zeros(ny+nx,tperiods);

%x        = [0; sigmaaparval];  % one standard deviation shock
x        = [0; 0.01];           % one percent shock

for t = 1:tperiods
    IRF(:,t) = [gx*x; x];
    x = hx*x;
end

% -- plot (in percent deviation)
yxst = [ystst,xstst]';
yx_  = strvcat(y_,x_);

nrow = 3; 
ncol = 3;

figure(1)
for jkl=1:(nx+ny);
    subplot(nrow,ncol,jkl)
    scalepar = yxst(jkl);
    plot(1:tperiods, IRF(jkl,:)/scalepar*100, 'k');
    title(yx_(jkl,:))
    ylabel('percent')
    xlabel('quarters')
    axis tight
end

orient landscape
print -dpdf rbcfig_2021.pdf

% --------------------------------------
% simulate series
% --------------------------------------
nsimul = 10000;
nburn  = 1000;

Inno = randn(1,nsimul+nburn);
x    = zeros(nx,1); 
YXsimul = zeros(ny+nx, nsimul+nburn);
eta = [0; sigmaaparval];

for jkl = 1:(nsimul+nburn);
    x = hx*x+eta*Inno(1,jkl);
    y = gx*x;
    YXsimul(:,jkl) = [y;x]./yxst;
end

% - drop burnin 
YXsimul = YXsimul(:, nburn+1:end);

% - hpfilter variables
[yt,yd] = hp(YXsimul',1600);
YXfilter = yd*100;

% - std deviations
ypos = 6;
VCOV = cov(YXfilter);
stds = sqrt(diag(VCOV));                % standard deviations
CORR = corr(YXfilter);
corrys = CORR(:,ypos);     % correlation with output
ar1 = [];
for jkl=1:(ny+nx);
    CORR = corr(YXfilter(1:end-1,jkl), YXfilter(2:end,jkl));
    ar1(jkl,1) = CORR(1,1);
end

% - moment table
info.rnames = strvcat(' ',yx_);
info.cnames = strvcat('std', 'relstd','ar1', 'corry');
mprint([stds, stds./stds(ypos) ar1, corrys], info)

 