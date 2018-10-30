function [result,details] = VSLiteHist(iyears,parameters)
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
% Advanced Inputs (must be specified by a single struct):
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
T1 = getparam(parameters, 'T1', []);
T2 = getparam(parameters, 'T2', []);
M1 = getparam(parameters, 'M1', []);
M2 = getparam(parameters, 'M2', []);
D1 = getparam(parameters, 'D1', []);
D2 = getparam(parameters, 'D2', []);
taui = getparam(parameters, 'taui', []);
taue = getparam(parameters, 'taue', []);
eoi = getparam(parameters, 'eoi', []);
dampth = getparam(parameters, 'dampth', exp(-1));
T = getparam(parameters, 'T', []);
P = getparam(parameters, 'P', []);
M = getparam(parameters, 'M', []);
phi = getparam(parameters, 'phi', 0);
D = getparam(parameters, 'D', []);
gE = getparam(parameters, 'gE', []);
gT = getparam(parameters, 'gT', []);
gM = getparam(parameters, 'gM', []);
g0 = getparam(parameters, 'g0', []);
eD = getparam(parameters, 'eD', []);
gD = getparam(parameters, 'gD', []);
[Mmax, Mmin, alph, m_th, mu_th, rootd, M0, substep] = ...
    dealarray(getparam(parameters, 'lbparams', [0.76, 0.01, 0.093, 4.886, 5.80, 1000, 0.2, 0]));
[I_0,I_f] = dealarray(getparam(parameters, 'intwindow', [1,12]));
outsel = getparam(parameters, 'outsel', 'trw');

% [result,details] = VSLiteHist_core(iyears, VSLiteHist_paramset(parameters));
[result,details] = VSLiteHist_core(iyears,T1,T2,M1,M2,D1,D2,taui,taue,eoi,dampth,T,P,M,phi,D,gE,gT,gM,g0,eD,gD,Mmax,Mmin,alph,m_th,mu_th,rootd,M0,substep,I_0,I_f,outsel);
end

function [ret] = getparam(parameters, paramname, defaultval)
if isfield(parameters, paramname); ret = parameters.(paramname);
else ret = defaultval;
end
end
