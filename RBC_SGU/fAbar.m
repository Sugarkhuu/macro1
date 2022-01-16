function [y] = fAbar(Aparval, etaparval,piparval, nuparval, betaparval, kappaparval, xiparval, chiparval)

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
y=ybar-1;
end

% wbar = etaparval*xbar + (1-etaparval)*piparval;
% Upsilonbar = xbar - wbar;
% Jbar = Upsilonbar/(1-(1-nuparval)*betaparval);
% stoch_betabar = betaparval;
% qbar = kappaparval/(stoch_betabar*Jbar);
% fbar = (chiparval*qbar^(xiparval-1))^(1/xiparval); 
% lbar = fbar/(fbar+nuparval);
% y = xbar*lbar - 1;

% A = fsolve(@(xbar)fAbar(xbar, etaparval,piparval, nuparval, betaparval, kappaparval, xiparval, chiparval),0.001)