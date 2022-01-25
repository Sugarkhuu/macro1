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
pi_star = 1;
b_star  = 1/2;
R_star  = pi_star/beta;
s_star  = (1/beta-1)*b_star/R_star;
p_star = 1;

shock_F = 0.01;0.01;
shock_M = 0;-0.01;

% ----------------------------
% Section 2. Compute steady state and implied params
% ----------------------------
pitbar = pi_star;
Rtbar = R_star;
stbar = s_star;
btbar = b_star;
ptbar = p_star;
Btbar = btbar*ptbar;
rho_eval = 0;
epsilon_par = 0;

xstst = [btbar,stbar]; 
ystst = [Rtbar, Rtbar, btbar, pitbar];
paramvals = [beta gamma alpha pi_star R_star s_star b_star rho_eval epsilon_par];

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
x        = [shock_F;shock_M]; %; ; shock_M          % one percent shock

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
eta = [shock_F; shock_M]; %[0; sigmaAparval];

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