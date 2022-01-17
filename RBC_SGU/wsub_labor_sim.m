% ----------------------------------------------------------------------------
% main_RBC_SGU.m
% 
% codes the RBC model in PhD Macro 1.
% Computes second moments and impulse responses using first-order
% perturbation.
% --------------------------------------------------------------------------

% - housekeeping
clear all, close all, clc
fontSize = 20;
% - search path
thispath_ = cd;
addpath([thispath_,'\auxfiles\'] )

% -----------------------------------------
% Section 1. Set parameters and targets
% -----------------------------------------

% - set parameters
zetaparval   = 0;
rhoparval    = 0.95;
betaparval   = 0.997;
nuparval     = 0.0265;
etaparval    = 0.4;
piparval     = 1.02; 
xiparval     = 0.5; 
chiparval    = 0.38; 
kappaparval  = 0.24;
sigmaAparval = 0.5/100;
Aparval      = 1.0576; % not given Should implement root finding method here!!!

% ----------------------------
% Section 2. Compute steady state and implied params
% ----------------------------
Abar = Aparval;
xbar = Abar;
wbar = etaparval*xbar + (1-etaparval)*piparval;
Upsilonbar = xbar - wbar;
Jbar = Upsilonbar/(1-(1-nuparval)*betaparval);
stoch_betabar = betaparval;
qbar = kappaparval/(stoch_betabar*Jbar);
fbar = (chiparval*qbar^(xiparval-1))^(1/xiparval); 
ebar = fbar/(fbar+nuparval);
ubar = 1 - ebar;
mbar = fbar*ubar;
vbar = (mbar/(chiparval*ubar^xiparval))^(1/(1-xiparval));
lbar = ebar;
ybar = Abar*lbar;
cbar = ybar-kappaparval*vbar;
lambdabar = 1/cbar;
Psibar = Upsilonbar*ebar-kappaparval*vbar;
taubar = 0;
w_firmbar = wbar;
%% checks
% mbar - chiparval*ubar^xiparval*vbar^(1-xiparval)
% % vbar = mbar/qbar;
% cbar-ebar*wbar - Psibar

xstst = [Abar];
ystst = [ebar, cbar, wbar, Psibar, lambdabar, stoch_betabar, ybar, lbar, xbar, ubar, fbar, Upsilonbar, Jbar, vbar, mbar, qbar, taubar, w_firmbar];
paramvals = [zetaparval Aparval rhoparval betaparval nuparval etaparval piparval xiparval chiparval kappaparval sigmaAparval];

% --------------------------
% Section 2. run model file and get derivatives
% --------------------------
wsub_labor_model;

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
tperiods = 40;25;
IRF      = zeros(ny+nx,tperiods);

% x        = [0; sigmaAparval];     % one standard deviation shock
x        = [-0.01];[0; -0.01];           % one percent shock

for t = 1:tperiods
    IRF(:,t) = [gx*x; x];
    x = hx*x;
end

% -- plot (in percent deviation)
yxst = [ystst,xstst]';
yx_  = strvcat(y_,x_);

nrow = 4; 
ncol = 5;

figure(1)
for jkl=1:(nx+ny);
    subplot(nrow,ncol,jkl)
    scalepar = yxst(jkl);
    plot(1:tperiods, IRF(jkl,:)/scalepar*100, 'k');
    title(yx_(jkl,:), 'Interpreter','None', 'FontSize', fontSize)
    ylabel('percent', 'FontSize', fontSize)
    xlabel('quarters', 'FontSize', fontSize)
    axis tight
end

vars = {'yt','ut'};
figure(2)
for i=1:numel(vars);
    jkl = find(strcmp(cellstr(yx_),vars{i}));
    subplot(2,1,i)
    scalepar = yxst(jkl);
    plot(1:tperiods, IRF(jkl,:)/scalepar*100, 'k');
    title(yx_(jkl,:), 'Interpreter','None', 'FontSize', fontSize)
    ylabel('percent', 'FontSize', fontSize)
    xlabel('quarters', 'FontSize', fontSize)
    axis tight
end

vars = {'yt','wt','w_firmt'};
subplot(2,1,1)
jkl = find(strcmp(cellstr(yx_),vars{1}));
scalepar = yxst(jkl);
plot(1:tperiods, IRF(jkl,:)/scalepar*100, 'k');
title(yx_(jkl,:), 'Interpreter','None','FontSize', fontSize)
ylabel('percent','FontSize', fontSize)
xlabel('quarters','FontSize', fontSize)
axis tight
subplot(2,1,2)
jkl = find(strcmp(cellstr(yx_),vars{2}));
jkl_ = find(strcmp(cellstr(yx_),vars{3}));
scalepar = yxst(jkl);
plot(1:tperiods, (IRF(jkl,:)-IRF(jkl_,:))/scalepar*100, 'k');
title('Subsidy (% of steady state wage)', 'Interpreter','None','FontSize', fontSize)
ylabel('percent','FontSize', fontSize)
xlabel('quarters','FontSize', fontSize)
axis tight

% --------------------------------------
% simulate series
% --------------------------------------
nsimul = 10000;
nburn  = 1000;

Inno = randn(1,nsimul+nburn);
x    = zeros(nx,1); 
YXsimul = zeros(ny+nx, nsimul+nburn);
eta = [sigmaAparval];[0; sigmaAparval];
% keyboard
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
ypos = find(strcmp(cellstr(yx_),'yt'));6;
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

 