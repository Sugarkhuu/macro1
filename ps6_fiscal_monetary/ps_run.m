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
betaparval   = 0.99;
gammaparval  = 0.5;
alphaparval  = 1.5;
pi_starparval = 1;
b_starparval  = 1/2;
R_starparval  = pi_starparval/betaparval;
s_starparval  = (1/betaparval-1)*b_starparval/R_starparval;
% p_starparval = 1;
rho_eparval = 0.9;
eparparval = 0;

shock_aux = 0;
shock_F   = 0.01;0.01;
shock_M   = 0;0.01;
% ----------------------------
% Section 2. Compute steady state and implied params
% ----------------------------
pitbar = pi_starparval;
Rtbar = R_starparval;
stbar = s_starparval;
btbar = b_starparval;
% ptbar = p_starparval;
% Btbar = btbar*ptbar;

xstst = [Rtbar,stbar,eparparval]; 
ystst = [btbar,btbar,pitbar,Rtbar];
paramvals = [betaparval gammaparval alphaparval pi_starparval R_starparval s_starparval b_starparval rho_eparval eparparval];

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
tperiods = 25;
IRF      = zeros(ny+nx,tperiods);

% x        = [0; sigmaAparval];     % one standard deviation shock
x        = [shock_aux;shock_F;shock_M]; %; ; shock_M          % one percent shock

for t = 1:tperiods
    IRF(:,t) = [gx*x; x];
    x = hx*x;
end

% -- plot (in percent deviation)
yxst = [ystst,xstst]';
yx_  = strvcat(y_,x_);

nrow = ceil((nx+ny)/3); 
ncol = 3;

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

% % variables to plot
% vars = {'yt','ut'};
% figure(2)
% for i=1:numel(vars);
%     jkl = find(strcmp(cellstr(yx_),vars{i}));
%     subplot(2,1,i)
%     scalepar = yxst(jkl);
%     plot(1:tperiods, IRF(jkl,:)/scalepar*100, 'k');
%     title(yx_(jkl,:), 'Interpreter','None', 'FontSize', fontSize)
%     ylabel('percent', 'FontSize', fontSize)
%     xlabel('quarters', 'FontSize', fontSize)
%     axis tight
% end

% --------------------------------------
% simulate series
% --------------------------------------
nsimul = 10000;
nburn  = 1000;

Inno = randn(1,nsimul+nburn);
x    = zeros(nx,1); 
YXsimul = zeros(ny+nx, nsimul+nburn);
eta = [shock_aux;shock_F;shock_M]; %[0; sigmaAparval];

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
ypos = find(strcmp(cellstr(yx_),'Rt'));6;
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

% % save baseline data in main_data.mat
% main_data = struct();
% main_data.yx_ = yx_;
% main_data.IRF = IRF;
% main_data.YXsimul  = YXsimul';
% main_data.YXfilter = YXfilter;

% save main_data.mat main_data;
