% given parameters
rhoparval    = 0.95;
etaparval    = 0.4;
piparval     = 1.02; 
nuparval     = 0.0265;
betaparval   = 0.997;
kappaparval  = 0.24;
xiparval     = 0.5; 
chiparval    = 0.38; 

[A, fval] = fminbnd(@(Abar)fAbar(Abar, etaparval,piparval, nuparval, betaparval, kappaparval, xiparval, chiparval),1, 10)
disp(["Found steady state value A = " A]);


function [y] = fAbar(Abar, etaparval,piparval, nuparval, betaparval, kappaparval, xiparval, chiparval)
    % steady state values
    xbar = Abar;
    wbar = etaparval*xbar + (1-etaparval)*piparval;
    Upsilonbar = xbar - wbar;
    Jbar = Upsilonbar/(1-(1-nuparval)*betaparval);
    stoch_betabar = betaparval;
    qbar = kappaparval/(stoch_betabar*Jbar);
    fbar = (chiparval*qbar^(xiparval-1))^(1/xiparval); 
    lbar = fbar/(fbar+nuparval);

    % error
    y = abs(Abar*lbar - 1);
end