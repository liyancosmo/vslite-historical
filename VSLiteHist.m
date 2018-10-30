function [result,details] = VSLiteHist(iyears,varargin)
% VSLite_v2_3.m - Simulate tree ring width index given monthly climate inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Usage:
%    trw = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P)
%    gives just simulated tree ring as ouput.
%
%   [trw,gT,gM,gE,M] = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P))
%    also includes growth response to temperature, growth response to soil
%    moisture, scaled insolation index, and soil moisture estimate in outputs.
%
% Basic Inputs:
%   syear = start year of simulation.
%   eyear = end year of simulation.
%
% Advanced Inputs (must be specified as property/value pairs):
%   'T1':        scalar temperature threshold below which temp. growth response is zero (in deg. C)
%   'T2':        scalar temperature threshold above which temp. growth response is one (in deg. C)
%   'M1':        scalar soil moisture threshold below which moist. growth response is zero (in v/v)
%   'M2':        scalar soil moisture threshold above which moist. growth response is one (in v/v)
%   'D1':        scalar historical disturbance threshold below which historical disturbance growth response is zero (in times)
%   'D2':        scalar historical disturbance threshold above which historical disturbance growth response is zero (in times)
%   'taui':      time constant of inhibition effect of historical disturbance (in years)
%   'taue':      time constant of promotion effect of historical disturbance (in years)
%   'eoi':       historical disturbance effect ratio (strength of promotion effect versus strength of inhibition effect)
%                (the above parameters are all required. Note that these growth re-
%                 sponse parameters may be estimated using code estimate_vslite_pa-
%                 rams_v2_3.m)
%   'dampth'     scaler strength threshold below which historical disturba-
%                nce effect is considered as vanished. default is the Euler
%                number (e = 2.718...).
%   'T':         (12 x nyrs) matrix of ordered mean monthly temperatures 
%                (in degEes C) (required)
%   'P':         (12 x nyrs) matrix of ordered accumulated monthly precipi-
%                tation (in mm) (required if 'M' is absent)
%   'M':         (12 x nyrs) matrix of ordered calculated monthly moisture
%                (required if 'P' is absent)
%   'phi':       latitude of site (in degrees N) (required)
%   'D':         historical disturbance strenghs (required)
%   'gE':        calculated 12x1 growth response to insolation
%   'gT':        calculated 12xN growth response to temperature
%   'gM':        calculated 12xN growth response to moisture
%   'g0':        calculated 1xN growth response to present conditions
%   'eD':        calculated 1xN present effects of historical disturbances
%   'gD':        calculated 1xN growth response to historical disturbances
%   'lbparams':  Parameters of the Leaky Bucket model of soil moisture.
%                These may be specified in an 8 x 1 vector in the following
%                order (otherwise the default values are read in):
%                   Mmax: scalar maximum soil moisture content (in v/v),
%                     default value is 0.76
%                   Mmin: scalar minimum soil moisture (in v/v), default
%                     value is 0.01
%                   alph: scalar runoff parameter 1 (in inverse months),
%                     default value is 0.093
%                   m_th: scalar runoff parameter 3 (unitless), default
%                     value is 4.886
%                   mu_th: scalar runoff parameter 2 (unitless), default
%                     value is 5.80
%                   rootd: scalar root/"bucket" depth (in mm), default
%                     value is 1000
%                   M0: initial value for previous month's soil moisture at
%                     t = 1 (in v/v), default value is 0.2
%                   substep: logical 1 or 0; perform monthly substepping in
%                     leaky bucket (1) or not (0)? Default value is 0.
%   'intwindow': Integration window. Which months' growth responses should
%                be integrated to compute the annual ring-width index?
%                Specified as a 2 x 1 vector of integer values. Both
%                elements are given in integer number of months since January
%                (July) 1st of the current year in the Northern (Southern)
%                hemisphere, and specify the beginning and end of the integration
%                window, respectively. Defaults is [1 ; 12] (eg. integrate
%                response to climate over the corresponding calendar year,
%                assuming location is in the northern hemisphere).
%   'outsel':    Vhe output growth selection. VSLiteHist will calculate only
%                the necessary parts for "outsel". And will return the value
%                of "outsel" as the first output. All of the output are
%                standarized by Z score. Selections are:
%                    'trw':  overall tree-ring output (default option)
%                    'gE':   12x1 growth response to insolation
%                    'gT':   1xN integrated growth response to temperature
%                    'gM':   1xN integrated growth response to moisture
%                    'gET':  1xN integrated growth response to temperature multiplied by gE
%                    'gTE':  transposed gET
%                    'gEM':  1xN integrated growth response to moisture multiplied by gE
%                    'gME':  transposed gEM
%                    'g0':   1xN growth response to present conditions
%                    'gD':   1xN growth response to historical disturbances
%
%
%
% Outputs:
%   'result':   The output growth component specified by 'outsel'. If 'outsel' 
%               is absent, then the overall predicted tree-ring is output.
%   'details':  A struct that stores the calculated components. Available co-
%               mponents are: T, P, M, D, gE, gT, gM, g0, eD, gD. These comp-
%               onents are the same as the corresponding advanced inputs.
%               (refer to the corresponding parts in the "Advanced Inputs" for
%                their descriptions.)
%
%
% For more detailed documentation, see:
% 1) Tolwinski-Ward et al., An efficient forward model of the climate
% controls on interannual variation in tree-ring width, Climate Dynamics (2011)
% DOI: 10.1007/s00382-010-0945-5
%
% 2) Tolwinski-Ward et al., Erratum to: An efficient forward model of the climate
% controls on interannual variation in tree-ring width, Climate Dynamics (2011)
% DOI: 10.1007/s00382-011-1062-9
%
% 3) Tolwinski-Ward et al., Bayesian parameter estimation and
% interpretation for an intermediate model of tree-ring width, Clim. Past
% (2013), DOI: 10.5194/cp-9-1-2013
%
% 4) Documentation available with the model at http://www.ncdc.noaa.gov/paleo/softlib/softlib.html
%
% Revision History
% v0.1 - Original coding at monthly timestep from full daily timestep model (SETW, 4/09)
% v1.0 - Changed soil moisture module to the CPC Leaky Bucket model (SETW, 5/09)
% v1.1 - No upper parametric bounds for gT, gW as in full model; no density module (SETW, 9/09)
% v1.2 - Added adjustable integration window parameters (SETW, 1/10)
% v2.0 - Minor debugging for Octave compatibility, final version for publication (SETW, 10/10)
% v2.1 - Error in evapotranspiration calculation corrected (SETW, 7/11)
% v2.2 - Add switch to allow for monthly sub-stepping in soil moisture computation (SETW, N.Graham, K.Georgakaos, 9/11)
% v2.3 - Add switch to allow moisture M to be given as input rather than estimated
%        from T and P; add variable input options and improve commenting (SETW, 7/13)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nyrs = length(iyears);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in advanced inputs if user-specified; else read in parameter defaults:
varargin = VarArgs(varargin);
T1 = varargin.get('T1', []);
T2 = varargin.get('T2', []);
M1 = varargin.get('M1', []);
M2 = varargin.get('M2', []);
D1 = varargin.get('D1', []);
D2 = varargin.get('D2', []);
taui = varargin.get('taui', []);
taue = varargin.get('taue', []);
eoi = varargin.get('eoi', []);
dampth = varargin.get('dampth', exp(1));
T = varargin.get('T', []);
P = varargin.get('P', []);
M = varargin.get('M', []);
phi = varargin.get('phi', 0);
D = varargin.get('D', []);
gE = varargin.get('gE', []);
gT = varargin.get('gT', []);
gM = varargin.get('gM', []);
g0 = varargin.get('g0', []);
eD = varargin.get('eD', []);
gD = varargin.get('gD', []);
[Mmax, Mmin, alph, m_th, mu_th, rootd, M0, substep] = ...
    dealarray(varargin.get('lbparams', [0.76, 0.01, 0.093, 4.886, 5.80, 1000, 0.2, 0]));
[I_0,I_f] = dealarray(varargin.get('intwindow', [1,12]));
outsel = varargin.get('outsel', 'trw');
%%% check input params %%%
if isempty(T); throw(MException('VSLiteHist:VSLiteHist', 'T is not set')); end
if isempty(phi); throw(MException('VSLiteHist:VSLiteHist', 'phi is not set')); end
if isempty(D); throw(MException('VSLiteHist:VSLiteHist', 'D is not set')); end
if isempty(P) && isempty(M); throw(MException('VSLiteHist:VSLiteHist', 'neither P and M is set')); end
%%% convert taui&taue to natural exponential time scale
taui = -taui/log(dampth);
taue = -taue/log(dampth);
dampth = exp(1);
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
% compute g0
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
% stardarize the output by Z score
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
