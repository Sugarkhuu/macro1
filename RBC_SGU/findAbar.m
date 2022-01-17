[A, fval] = fminbnd(@(xbar)fAbar(xbar, etaparval,piparval, nuparval, betaparval, kappaparval, xiparval, chiparval),1, 10)

disp(["Found steady state value A = " A]);

function [y] = fAbar(xbar, etaparval,piparval, nuparval, betaparval, kappaparval, xiparval, chiparval)
wbar = etaparval*xbar + (1-etaparval)*piparval;
Upsilonbar = xbar - wbar;
Jbar = Upsilonbar/(1-(1-nuparval)*betaparval);
stoch_betabar = betaparval;
qbar = kappaparval/(stoch_betabar*Jbar);
fbar = (chiparval*qbar^(xiparval-1))^(1/xiparval); 
lbar = fbar/(fbar+nuparval);
y = abs(xbar*lbar - 1);
end