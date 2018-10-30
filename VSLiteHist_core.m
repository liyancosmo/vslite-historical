function [result,details] = core_VSLiteHist(iyears,T1,T2,M1,M2,D1,D2,taui,taue,eoi,dampth,T,P,M,phi,D,gE,gT,gM,g0,eD,gD,Mmax,Mmin,alph,m_th,mu_th,rootd,M0,substep,I_0,I_f,outsel)
% Core implementation of VSLiteHist.
% Refer to VSLiteHist.m for argument descriptions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [T1,T2,M1,M2,D1,D2,taui,taue,eoi,dampth,T,P,M,phi,D,gE,gT,gM,g0,eD,gD,Mmax,Mmin,alph,m_th,mu_th,rootd,M0,substep,I_0,I_f,outsel] = deal(paramset{:});
nyrs = length(iyears);
%%% check input params %%%
if isempty(T); throw(MException('VSLiteHist:VSLiteHist', 'T is not set')); end
if isempty(phi); throw(MException('VSLiteHist:VSLiteHist', 'phi is not set')); end
if isempty(D); throw(MException('VSLiteHist:VSLiteHist', 'D is not set')); end
if isempty(P) && isempty(M); throw(MException('VSLiteHist:VSLiteHist', 'neither P and M is set')); end
%%% convert taui&taue to natural exponential time scale
taui = -taui/log(dampth);
taue = -taue/log(dampth);
dampth = exp(-1);
%%% convert outsel to boolean flags
outtrw = strcmpi(outsel, 'trw');
outg0 = strcmpi(outsel, 'g0');
outgET = strcmpi(outsel, 'gET');
outgTE = strcmpi(outsel, 'gTE');
outgEM = strcmpi(outsel, 'gEM');
outgME = strcmpi(outsel, 'gME');
outgE = strcmpi(outsel, 'gE');
outgT = strcmpi(outsel, 'gT');
outgM = strcmpi(outsel, 'gM');
outgD = strcmpi(outsel, 'gD');
%%% decide the parts required to compute
needtrw = outtrw;
needg0 = needtrw || outg0;
needgET = outgET;
needgTE = outgTE;
needgEM = outgEM;
needgME = outgME;
needgE = needg0 || needgTE || needgET || needgME || needgMT || outgE;
needgT = needg0 || needgTE || needgET || outgT;
needgM = needg0 || needgME || needgEM || outgM;
needM = needgM;
needgD = needtrw || outgD;
needeD = needgD;
%%% specify a negative I_0 means we use previous year's data. as we skipped
%%% 0 here (I_0=1 for first month of current year, while I_0=-1 for last
%%% month of previous year), we increase negative I_0 by 1.
if I_0<0; I_0=I_0+1; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the required parts
% compute gE
if needgE && isempty(gE)
    gE = Compute_gE(phi);
    diaggE = diag(gE);
end
% compute gT
if needgT && isempty(gT)
    gT = NaN(12,nyrs);
    gT(T<T1) = 0;
    gT(T>T2) = 1;
    gT((T>=T1)&(T<=T2)) = (T((T>=T1)&(T<=T2))-T1)/(T2-T1);
end
% compute M
if needM && isempty(M)
    M = NaN(12,nyrs);
    if substep == 1
        M = leakybucket_submonthly(iyears,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
    elseif substep == 0
        M = leakybucket_monthly(iyears,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
    elseif substep ~=1 && substep ~= 0
        disp('''substep'' must either be set to 1 or 0.');
        return
    end
end
% compute gM
if needgM && isempty(gM)
    gM = NaN(12,nyrs);
    gM(M<M1) = 0;
    gM(M>M2) = 1;
    gM((M>=M1)&(M<=M2)) = (M((M>=M1)&(M<=M2))-M1)/(M2-M1);
end
% compute gEM and gME
if needgEM || needgME
    gEM = diaggE * gM;
    gME = gEM';
end
% compute gET and gTE
if needgET || needgTE
    gET = diaggE * gT;
    gTE = gET';
end
% compute g0 %%% consume: 0.0001
if needg0 && isempty(g0)
    gr = diag(gE)*min(gT,gM);
    g0 = integrate_months(gr,phi,I_0,I_f);
end
% compute eD
if needeD && isempty(eD)
    eD = disturbance_effects(D,taui,taue,eoi);
end
% compute gD
if needgD && isempty(gD)
    gD = NaN(1,nyrs);
    gD(eD<D1) = 0;
    gD(eD>D2) = 1;
    gD((eD>=D1)&(eD<=D2)) = (eD((eD>=D1)&(eD<=D2))-D1)/(D2-D1);
end
% compute trw
if needtrw
    trw = g0.*gD;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% handle the result component
% get the selected output
if outtrw; result = trw; end
if outg0;  result = g0; end
if outgET; result = integrate_months(gET,phi,I_0,I_f); end
if outgTE; result = integrate_months(gTE',phi,I_0,I_f); end
if outgEM; result = integrate_months(gEM,phi,I_0,I_f); end
if outgME; result = integrate_months(gME',phi,I_0,I_f); end
if outgE;  result = integrate_months(repmat(gE,[1,nyrs]),phi,I_0,I_f); end
if outgT;  result = integrate_months(gT,phi,I_0,I_f); end
if outgM;  result = integrate_months(gM,phi,I_0,I_f); end
if outgD;  result = gD; end
% stardarize the output by Z score %%% consume 0.00007~0.00010
result = (result-mean(result))/std(result);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% store the detailed components into "details"
if ~isempty(T);  details.T  = T;  end
if ~isempty(P);  details.P  = P;  end
if ~isempty(M);  details.M  = M;  end
if ~isempty(D);  details.D  = D;  end
if ~isempty(gE); details.gE = gE; end
if ~isempty(gT); details.gT = gT; end
if ~isempty(gM); details.gM = gM; end
if ~isempty(g0); details.g0 = g0; end
if ~isempty(eD); details.eD = eD; end
if ~isempty(gD); details.gD = gD; end
end % END VSLiteHist

function [gout] = integrate_months(gin, phi, I_0, I_f)
%%% integrate 12xN monthly data into 1xN annual data
%%% according to the latitute and integrating window
if phi>0 % if site is in the Northern Hemisphere:
    grstack = [NaN(12,1),gin(:,1:end-1);gin];
    grstack(isnan(grstack(:,1)),1) = mean( grstack(isnan(grstack(:,1)),2:end), 2 );
    startmo = 12+I_0;
    endmo = 12+I_f;
elseif phi<0 % if site is in the Southern Hemisphere:
    grstack = [gin;gin(:,2:end),NaN(12,1)];
    grstack(isnan(grstack(:,1)),end) = mean( grstack(isnan(grstack(:,1)),1:end-1), 2 );
    % (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo = 6+I_0; % (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo = 6+I_f; % (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
end
gout = sum(grstack(startmo:endmo,:));
end