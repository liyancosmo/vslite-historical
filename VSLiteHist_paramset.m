function [paramset] = VSLiteHist_paramset(paramstruct)
T1 = getparam(paramstruct, 'T1', []);
T2 = getparam(paramstruct, 'T2', []);
M1 = getparam(paramstruct, 'M1', []);
M2 = getparam(paramstruct, 'M2', []);
D1 = getparam(paramstruct, 'D1', []);
D2 = getparam(paramstruct, 'D2', []);
taui = getparam(paramstruct, 'taui', []);
taue = getparam(paramstruct, 'taue', []);
eoi = getparam(paramstruct, 'eoi', []);
dampth = getparam(paramstruct, 'dampth', exp(1));
T = getparam(paramstruct, 'T', []);
P = getparam(paramstruct, 'P', []);
M = getparam(paramstruct, 'M', []);
phi = getparam(paramstruct, 'phi', 0);
D = getparam(paramstruct, 'D', []);
gE = getparam(paramstruct, 'gE', []);
gT = getparam(paramstruct, 'gT', []);
gM = getparam(paramstruct, 'gM', []);
g0 = getparam(paramstruct, 'g0', []);
eD = getparam(paramstruct, 'eD', []);
gD = getparam(paramstruct, 'gD', []);
[Mmax, Mmin, alph, m_th, mu_th, rootd, M0, substep] = ...
    dealarray(getparam(paramstruct, 'lbparams', [0.76, 0.01, 0.093, 4.886, 5.80, 1000, 0.2, 0]));
[I_0,I_f] = dealarray(getparam(paramstruct, 'intwindow', [1,12]));
outsel = getparam(paramstruct, 'outsel', 'trw');
paramset = {T1,T2,M1,M2,D1,D2,taui,taue,eoi,dampth,T,P,M,phi,D,gE,gT,gM,g0,eD,gD,Mmax,Mmin,alph,m_th,mu_th,rootd,M0,substep,I_0,I_f,outsel};
end

function [ret] = getparam(parameters, paramname, defaultval)
if isfield(parameters, paramname); ret = parameters.(paramname);
else ret = defaultval;
end
end
