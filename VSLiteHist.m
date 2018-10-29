function [trw,varargout] = VSLiteHist(iyears,varargin)
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
%   'T':         (12 x Nyrs) matrix of ordered mean monthly temperatures 
%                (in degEes C) (required)
%   'P':         (12 x Nyrs) matrix of ordered accumulated monthly precipi-
%                tation (in mm) (required if 'M' is absent)
%   'M':         (12 x Nyrs) matrix of ordered calculated monthly moisture
%                (required if 'P' is absent)
%   'phi':       latitude of site (in degrees N) (required)
%   'gE':        calculated monthly growth response of insolation
%   'D':         historical disturbance strenghs (required)
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
%                be intregrated to compute the annual ring-width index?
%                Specified as a 2 x 1 vector of integer values. Both
%                elements are given in integer number of months since January
%                (July) 1st of the current year in the Northern (Southern)
%                hemisphere, and specify the beginning and end of the integration
%                window, respectively. Defaults is [1 ; 12] (eg. integrate
%                response to climate over the corresponding calendar year,
%                assuming location is in the northern hemisphere).
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
gE = varargin.get('gE', []);
D = varargin.get('D', []);
[Mmax, Mmin, alph, m_th, mu_th, rootd, M0, substep] = ...
    dealarray(varargin.get('lbparams', [0.76, 0.01, 0.093, 4.886, 5.80, 1000, 0.2, 0]));
[I_0,I_f] = dealarray(varargin.get('intwindow', [1,12]));
%%% check input params %%%
if isempty(T); throw(MException('VSLiteHist:VSLiteHist', 'T is not set')); end
if isempty(phi); throw(MException('VSLiteHist:VSLiteHist', 'phi is not set')); end
if isempty(D); throw(MException('VSLiteHist:VSLiteHist', 'D is not set')); end
if isempty(P) && isempty(M); throw(MException('VSLiteHist:VSLiteHist', 'neither P and M is set')); end
%%% convert taui&taue to natural exponential time scale
taui = -taui/log(dampth);
taue = -taue/log(dampth);
damph = exp(1);
%%% specify a negative I_0 means we use previous year's data. as we skipped
%%% 0 here (I_0=1 for first month of current year, while I_0=-1 for last
%%% month of previous year), we increase negative I_0 by 1.
if I_0<0; I_0=I_0+1; end
%%% Pre-allocate storage for outputs: %%%%
gT = NaN(12,nyrs);
gM = NaN(12,nyrs);
gD = NaN(1,nyrs);
eD = NaN(1,nyrs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in or estimate soil moisture:
if isempty(M)
    M = NaN(12,nyrs);
    % Compute soil moisture:
    if substep == 1
        M = leakybucket_submonthly(iyears,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
    elseif substep == 0
        M = leakybucket_monthly(iyears,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
    elseif substep ~=1 && substep ~= 0
        disp('''substep'' must either be set to 1 or 0.');
        return
    end
end
% Compute gE, the scaled monthly proxy for insolation:
if isempty(gE)
    gE = Compute_gE(phi);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now compute growth responses to climate, and simulate proxy:
%%%%%%%%%%%%%%%%
% syear = start (first) year of simulation
% eyear = end (last) year of simulation
% Compute monthly growth response to T & M, and overall growth response g0:
gT(T<T1) = 0;
gT(T>T2) = 1;
gT((T>=T1)&(T<=T2)) = (T((T>=T1)&(T<=T2))-T1)/(T2-T1);
gM(M<M1) = 0;
gM(M>M2) = 1;
gM((M>=M1)&(M<=M2)) = (M((M>=M1)&(M<=M2))-M1)/(M2-M1);
gr = diag(gE)*min(gT,gM);
%%%%%%%%%%%%%% Compute proxy quantity from growth responses %%%%%%%%%%%%%%%
if phi>0 % if site is in the Northern Hemisphere:
    grstack = [NaN(12,1),gr(:,1:end-1);gr];
    grstack(isnan(grstack(:,1)),1) = mean( grstack(isnan(grstack(:,1)),2:end), 2 );
    startmo = 12+I_0;
    endmo = 12+I_f;
elseif phi<0 % if site is in the Southern Hemisphere:
    grstack = [gr;gr(:,2:end),NaN(12,1)];
    grstack(isnan(grstack(:,1)),end) = mean( grstack(isnan(grstack(:,1)),1:end-1), 2 );
    % (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo = 6+I_0; % (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo = 6+I_f; % (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
end
g0 = sum(grstack(startmo:endmo,:));
% compute growth response to historical disturbance
for i = 1:nyrs
    Di = D(1:i);
    tau = i - (1:i);
    eD(i) = disturbance_effect(Di,tau,taue,taui,eoi);
end
gD(eD<D1) = 0;
gD(eD>D2) = 1;
gD((eD>=D1)&(eD<=D2)) = (eD((eD>=D1)&(eD<=D2))-D1)/(D2-D1);
% compute final response trwidth
trwidth = g0.*gD;

%
trw = ((trwidth-mean(trwidth))/std(trwidth))'; % proxy series is standardized width.
%
if nargout >=1
    varargout(1) = {gT};
    varargout(2) = {gM};
    varargout(3) = {gE};
    varargout(4) = {gD};
    varargout(5) = {M};
    varargout(6) = {trwidth};
    varargout(7) = {mean(trwidth)};
    varargout(8) = {std(trwidth)};
end
%
end
