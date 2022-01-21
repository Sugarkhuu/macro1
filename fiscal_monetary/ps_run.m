% ----------------------------------------------------------------------------
% ps_run.m
% 
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
beta   = 0.99;
gamma  = 0.5;
alpha  = 1.5;
PIstar = 1;
ypar = 1;
bstar  = ypar/2;
Rstar  = PIstar/beta;
sstar  = (1/beta-1)*bstar/Rstar;
zstar = 0;
% ----------------------------
% Section 2. Compute steady state and implied params
% ----------------------------
cbar = ypar;
PIbar = PIstar;
Pbar = 1;
Rbar = Rstar;
sbar = sstar;
bbar = bstar;
Bbar = bbar;
Bbbar = Bbar;
zbar = zstar;
taubar = zbar + sbar;


%% checks
% mbar - chiparval*ubar^xiparval*vbar^(1-xiparval)
% % vbar = mbar/qbar;
% cbar-ebar*wbar - Psibar

xstst = [Pbar, Bbar];
ystst = [cbar, PIbar, Rbar, sbar, taubar, zbar, Bbbar, bbar];
paramvals = [beta gamma alpha PIstar Rstar sstar bstar ypar zstar];

% --------------------------
% Section 2. run model file and get derivatives
% --------------------------
ps_mod;

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
x        = [0; -0.01];           % one percent shock

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

% variables to plot
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

% --------------------------------------
% simulate series
% --------------------------------------
nsimul = 10000;
nburn  = 1000;

Inno = randn(1,nsimul+nburn);
x    = zeros(nx,1); 
YXsimul = zeros(ny+nx, nsimul+nburn);
eta = [0; sigmaAparval];

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

% save baseline data in main_data.mat
main_data = struct();
main_data.yx_ = yx_;
main_data.IRF = IRF;
main_data.YXsimul  = YXsimul';
main_data.YXfilter = YXfilter;

save main_data.mat main_data;
