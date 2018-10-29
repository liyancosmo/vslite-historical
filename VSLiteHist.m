function [trw,varargout] = VSLiteHist(syear,eyear,varargin)
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
%                (Note that optimal growth response parameters T1, T2, M1, M2 may be estimated
%                 using code estimate_vslite_params_v2_3.m also freely available at
%                 the NOAA NCDC Paleoclimatology software library.)
%   'T':         (12 x Nyrs) matrix of ordered mean monthly temperatures 
%                (in degEes C) (required)
%   'P':         (12 x Nyrs) matrix of ordered accumulated monthly precipi-
%                tation (in mm) (required if 'M' is absent)
%   'M':         (12 x Nyrs) matrix of ordered calculated monthly moisture
%                (required if 'P' is absent)
%   'phi':       latitude of site (in degrees N) (required)
%   'gE':        calculated monthly growth response of insolation
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
iyear = syear:eyear;
nyrs = length(syear:eyear);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in advanced inputs if user-specified; else read in parameter defaults:
if nargin > 6
    %%%% First fill parameter values in with defaults: %%%%%
    % Parameters of the Leaky Bucket model:
    T1 = [];
    T2 = [];
    M1 = [];
    M2 = [];
    T = [];
    P = [];
    M = [];
    phi = [];
    gE = [];
    Mmax = 0.76;
    Mmin = 0.01;
    alph = 0.093;
    m_th = 4.886;
    mu_th = 5.80;
    rootd = 1000;
    M0 = 0.2;
    substep = 0;
    % Integration window parameters:
    I_0 = 1;
    I_f = 12;
    % then over-write defaults if user-specified:
    Nvararg = length(varargin);
    for i = 1:Nvararg/2
        namein = varargin{2*(i-1)+1};
        valin = varargin{2*i};
        switch namein
            case 'T1'
                T1 = valin;
            case 'T2'
                T2 = valin;
            case 'M1'
                M1 = valin;
            case 'M2'
                M2 = valin;
            case 'T'
                T = valin;
            case 'P'
                P = valin;
            case 'M'
                M = valin;
            case 'phi'
                phi = valin;
            case 'gE'
                gE = valin;
            case 'lbparams'
                Mmax = valin(1);
                Mmin = valin(2);
                alph = valin(3);
                m_th = valin(4);
                mu_th = valin(5);
                rootd = valin(6);
                M0 = valin(7);
                substep = valin(8);
            case 'intwindow'
                I_0 = valin(1);
                I_f = valin(2);
        end
    end
else
    throw(MException('VSLiteHist:VSLiteHist', 'no enough params received'));
end
%%% check input params %%%
if isempty(P) && isempty(M); throw(MException('VSLiteHist:VSLiteHist', 'neither P and M is set')); end
if isempty(phi) && isempty(gE); throw(MException('VSLiteHist:VSLiteHist', 'neither phi and gE is set')); end
%%% specify a negative I_0 means we use previous year's data. as we skipped
%%% 0 here (I_0=1 for first month of current year, while I_0=-1 for last
%%% month of previous year), we increase negative I_0 by 1.
if I_0<0; I_0=I_0+1; end
%%% Pre-allocate storage for outputs: %%%%
gT = NaN(12,nyrs);
gM = NaN(12,nyrs);
potEv = NaN(12,nyrs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in or estimate soil moisture:
if isempty(M)
    M = NaN(12,nyrs);
    % Compute soil moisture:
    if substep == 1;
        M = leakybucket_submonthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
    elseif substep == 0
        M = leakybucket_monthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
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
% Compute monthly growth response to T & M, and overall growth response G:
gT(T<T1) = 0;
gT(T>T2) = 1;
gT((T>=T1)&(T<=T2)) = (T((T>=T1)&(T<=T2))-T1)/(T2-T1);
gM(M<M1) = 0;
gM(M>M2) = 1;
gM((M>=M1)&(M<=M2)) = (M((M>=M1)&(M<=M2))-M1)/(M2-M1);
Gr = diag(gE)*min(gT,gM);
%%%%%%%%%%%%%% Compute proxy quantity from growth responses %%%%%%%%%%%%%%%
if phi>0 % if site is in the Northern Hemisphere:
    Grstack = [NaN(12,1),Gr(:,1:end-1);Gr];
    Grstack(isnan(Grstack(:,1)),1) = mean( Grstack(isnan(Grstack(:,1)),2:end), 2 );
    startmo = 12+I_0;
    endmo = 12+I_f;
elseif phi<0 % if site is in the Southern Hemisphere:
    Grstack = [Gr;Gr(:,2:end),NaN(12,1)];
    Grstack(isnan(Grstack(:,1)),end) = mean( Grstack(isnan(Grstack(:,1)),1:end-1), 2 );
    % (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo = 6+I_0; % (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo = 6+I_f; % (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
end
trwidth = sum(Grstack(startmo:endmo,:));
%
trw = ((trwidth-mean(trwidth))/std(trwidth))'; % proxy series is standardized width.
%
if nargout >=1
    varargout(1) = {gT};
    varargout(2) = {gM};
    varargout(3) = {gE};
    varargout(4) = {M};
    varargout(5) = {potEv};
    varargout{6} = {trwidth};
    varargout{7} = {mean(trwidth)};
    varargout{8} = {std(trwidth)};
end
%
end
