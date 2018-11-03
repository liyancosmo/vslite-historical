function [T1,T2,M1,M2,D1,D2,taui,taue,eoi,varargout] = estimate_vslitehist_params(RW,varargin)
% Given calibration-interval temperature, precipitation, and ring-width data,
% and site latitude, estimate_vslite_params_v2_3.m performes a Bayesian parameter
% estimation of the growth response parameters T1, T2, M1 and M2. The
% default prior distributions are based on a survey of current literature
% pertaining to biological growth thresholds in trees; however uniform or
% user-defined four-parameter-beta distributions may also optionally be
% used as priors.  The scheme supports an assumption of either independent, 
% Gaussian model errors, or AR(1) error structure, and in either case the
% parameters of the error model may also be estimated.
% 
% For more details, see Tolwinski-Ward et al., 'Bayesian parameter estimation
% and interpretation for an intermediate model of tree-ring width', Clim. Past, 
% 9, 1-13, 3013, doi: 10.5194/cp-9-1-2013
%
% Basic Usage: [T1,T2,M1,M2] = estimate_vslite_params_v2_3(RW,'T',T,'P',P,'phi',phi)
%
% Basic Usage Inputs:
% T: monthly calibration-interval temperature, dimension 12 x number of calibration years.
% P: monthly calibration-interval precipitation, dimension 12 x number of cal. years.
% phi: site latitude in degrees N.
% RW: standardized annual calibration-interval ring-width index.
%
% Basic Usage Ouptuts:
% T1, T2, M1, M2: point estimates given by the median of their respective
%                 posterior distributions.
%
% Advanced Usage: [T1,T2,M1,M2,varargout] = vslite_bayes_param_cal(T,P,phi,RW,varargin)
%
% Advanced Optional Inputs:
%     Must be specified as property/value pairs.  Valid property/value pairs are:
%     'T':            (12 x Nyrs) matrix of ordered mean monthly temperatures 
%                     (in degEes C) (required)
%     'P':            (12 x Nyrs) matrix of ordered accumulated monthly precipi-
%                     tation (in mm) (required if 'M' is absent)
%     'M':            (12 x Nyrs) matrix of ordered calculated monthly moisture
%                     (required if 'P' is absent)
%     'phi':          latitude of site (in degrees N) (required)
%     'gE':           calculated monthly growth response of insolation
%     'D':            historical disturbance strenghs (required)
%     'errormod'      Error model. Options are [0], [1], and [2] for white Gaussian
%                     noise, AR(1) model, or AR(2) model.  Default is [0].
%     'gparscalint'   Indices of years to use to estimate the growth response
%                     parameters T1, T2, M1, M2. Default is all years.
%     'eparscalint'   Indices of years to use to estimate the parameters of the
%                     error model if these parameters are to be estimated.
%                     Must be contiguous if using AR(1) error model. Default
%                     is all years. (Note: may underestimate error if not disjoint
%                     from interval used to fit growth response parameters
%                     as specified in 'gparscalint'.)
%     'errorpars'     Vector holding values of error model parameters is user
%                     wishes to fix their values rather than estimate them.
%                     For errormod == 0 (white noise model), values is a scalar
%                     with fixed value for sigma2w; for errormod == 1 (AR(1)),
%                     errorpars = [phi1 tau^2]; for errormod == 2 (AR(2)),
%                     errorpars = [phi1 phi2 tau^2]. No default (since default is
%                     to estimate parameters of a white noise error model).
%     'pt_ests'       Choices are ['mle'] or ['med'] to return either the
%                     posterior draw that maximizes the likelihood of the data
%                     or the marginal posterior medians as the parameter point
%                     estimates. Default is 'mle'.
%     'substep'       If hydroclim == 'P', then 'substep' is logical 0/1
%                     depending on whether leaky bucket model without/with
%                     substepping is preferred.  Default is [0].
%     'intwindow'     VS-Lite integration window, specified as vector [I_0 I_f]
%                     Default is [0 12].
%     'nsamp'         200<=integer<=10,000 fixing number of MCMC iterations.
%                     Default is [1000].
%     'nbi'           Number of burn-in samples. Default is [200].
%     'nchain'        Integer number of comp threads for the computation
%                     of Rhat. Default is [3].
%     'T1priorsupp'   2x1 vector with elements giving lower and upper bounds
%                     for support of uniform T1 prior. If not included in input
%                     argument list, default used is [0.0 8.5]; If set to
%                     a single scaler, it is not estimated and used as a
%                     certain parameter fed into the model.
%     'T2priorsupp'   T2 prior. Default is [9.0 20.0]
%     'M1priorsupp'   M1 prior. Default is [0.01 0.03]
%     'M2priorsupp'   M2 prior. Default is [0.1 0.5]
%     'D1priorsupp'   D1 prior. Default is [-5 -0.5]
%     'D2priorsupp'   D2 prior. Default is [0.2 2]
%     'tauipriorsupp' taui prior. Default is [0 10]
%     'tauepriorsupp' taue prior. Default is [100 500]
%     'eoipriorsupp'  eoi prior. Default is [0 1]
%     'dampth'        the scaler strength threshold below which historical 
%                     disturbance effect is considered as vanished.
%                     default is the Euler number (e = 2.718...)
%     'convthresh'    Scalar value greater than 0.  Threshold for MCMC
%                     convergence; warning is displayed if abs(Rhat-1)>convthresh.
%                     Default value is [0.1].
%     'verbose'       Logical [0] or [1]; print progress to screen? Default 
%                     is [1].
%
% Advanced Optional Ouptuts (must be specified in the following order):
% T1dist, T2dist, M1dist, M2dist: Returns the entire numerical posterior distributions
%                 of the growth response parameters if the user wants to check for
%                 convergence, autocorrelation, multi-modality, etc., or to
%                 use the full distributions of the parameters to quantify
%                 uncertainty.
% Rhats:          Returns the convergence statistics associated with T1, T2, M1, M2,
%                 and sigma2rw if it was estimated.
% convwarning:    Logical [0] or [1] depending on whether any Rhat values were
%                 outside of the threshold distance from 1.
% -- Next ordered outputs for white noise error model (errormod==0): -- %%%
% sig2rw:         Point estimate of model error variance
% sigma2rwdist:   Returns the entire numerical posterior distribution
% Gdist:          Returns the numerical posterior distribution of monthly growth
%                 responses to the input climate for the corresponding posterior
%                 parameter sets; has dimension 12 x Nyrs x Nsamps.
%                 
% -- Next ordered outputs for AR(1) error model (errormod==1): -- %%%%%%%%%
% phi1:           Point estimate of AR(1) coefficient
% phi1dist:       Numerical posterior distribution of AR(1) coefficient 
% tau2:           Point estimate of error model innovation variance
% tau2dist:       Numerical distribution of error model innovation variance
% Gdist:          Returns the numerical posterior distribution of monthly growth
%                 responses to the input climate for the corresponding posterior
%                 parameter sets; has dimension 12 x Nyrs x Nsamps.
%
% SETW 9/20/2011:  version 1.0. Estimates T1, T2, M1, M2 for fixed sigma2w, assuming
%                  normally- and independent and identically distributed model residuals.
% SETW 4/15/2012:  version 2.0. Added simultaneous estimation of sigma2w under assumption
%                  of inverse-gamma prior, literature-based priors, Rhat convergence
%                  metric, and output plots.
% SETW 12/15/2012: version 2.1 for publication; added options for user-defined
%                  prior distributions.
% SETW 5/10/2013:  version 2.2: Revised data-level model and added options for white or
%                  AR(1) error models following reviewer comments in
%                  Climate of the Past Discussions; options for flexible
%                  calibration-intervals for growth parameters and error
%                  model parameters also added; MLE added as option for 
%                  point estimates;
%                  version 2.3: additional commenting added; option to
%                  condition on user-supplied input soil moisture data included 
%                  as opposed to necessarily estimating M from T & P via
%                  Leaky Bucket.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargin = VarArgs(varargin);
T = varargin.get('T', []);
P = varargin.get('P', []);
M = varargin.get('M', []);
phi = varargin.get('phi', []);
gE = varargin.get('gE', []);
D = varargin.get('D', []);
errormod = varargin.get('errormod', 0);
gparscalint = varargin.get('gparscalint', 1:length(RW));
eparscalint = varargin.get('eparscalint', 1:length(RW));
errorpars = varargin.get('errorpars', 0.5);
pt_ests = varargin.get('pt_ests', 'mle');
substep = varargin.get('substep', 0);
intwindow = varargin.get('intwindow', [0 12]);
nsamp = varargin.get('nsamp', 1000);
nbi = varargin.get('nbi', 200);
nchain = varargin.get('nchain', 3);
[aT1,   bT1,   cT1]   = dealprior(varargin.get('T1priorsupp',   [0, 9]));
[aT2,   bT2,   cT2]   = dealprior(varargin.get('T2priorsupp',   [10, 24]));
[aM1,   bM1,   cM1]   = dealprior(varargin.get('M1priorsupp',   [0, .1]));
[aM2,   bM2,   cM2]   = dealprior(varargin.get('M2priorsupp',   [.1, .5]));
[aD1,   bD1,   cD1]   = dealprior(varargin.get('D1priorsupp',   [-10, -5]));
[aD2,   bD2,   cD2]   = dealprior(varargin.get('D2priorsupp',   [0, 5]));
[ataui, btaui, ctaui] = dealprior(varargin.get('tauipriorsupp', [0, 10]));
[ataue, btaue, ctaue] = dealprior(varargin.get('tauepriorsupp', [100, 500]));
[aeoi,  beoi,  ceoi]  = dealprior(varargin.get('eoipriorsupp',  [0, 1]));
dampth = varargin.get('dampth', exp(-1));
convthresh = varargin.get('convthresh', .1);
verbose = varargin.get('verbose', 1);
%%% check input params %%%
if isempty(T); throw(MException('VSLiteHist:estimate_params', 'T is not set')); end
if isempty(phi); throw(MException('VSLiteHist:estimate_params', 'phi is not set')); end
if isempty(D); throw(MException('VSLiteHist:estimate_params', 'D is not set')); end
if isempty(P) && isempty(M); throw(MException('VSLiteHist:estimate_params', 'neither P and M is set')); end
%%% convert taui&taue to natural exponential time scale
ataui = -ataui/log(dampth); btaui = -btaui/log(dampth);
ataue = -ataue/log(dampth); btaue = -btaue/log(dampth);
dampth = exp(-1);
%%% Take zscore of RW data to fulfill assumptions of model error/noise structure
RW = zscore(RW);
Nyrs = size(T,2);
%%%% Now do the MCMC sampling: %%%%%%%%%%%%%
if ~cT1; Ttchains = NaN(nsamp+nbi, nchain); Ttensemb = NaN(1, nsamp*nchain); end
if ~cT2; Tochains = NaN(nsamp+nbi, nchain); Toensemb = NaN(1, nsamp*nchain); end
if ~cM1; Mtchains = NaN(nsamp+nbi, nchain); Mtensemb = NaN(1, nsamp*nchain); end
if ~cM2; Mochains = NaN(nsamp+nbi, nchain); Moensemb = NaN(1, nsamp*nchain); end
if ~cD1; Dtchains = NaN(nsamp+nbi, nchain); Dtensemb = NaN(1, nsamp*nchain); end
if ~cD2; Dochains = NaN(nsamp+nbi, nchain); Doensemb = NaN(1, nsamp*nchain); end
if ~ctaui; tauichains = NaN(nsamp+nbi, nchain); tauiensemb = NaN(1, nsamp*nchain); end
if ~ctaue; tauechains = NaN(nsamp+nbi, nchain); taueensemb = NaN(1, nsamp*nchain); end
if ~ceoi; eoichains = NaN(nsamp+nbi, nchain); eoiensemb = NaN(1, nsamp*nchain); end
if errormod == 0
    sig2rwchains = NaN(nsamp+nbi, nchain); sig2rwensemb = NaN(1, nsamp*nchain);
elseif errormod == 1
    phi1chains = NaN(nsamp+nbi, nchain); phi1ensemb = NaN(1, nsamp*nchain);
    tau2chains = NaN(nsamp+nbi, nchain); tau2ensemb = NaN(1, nsamp*nchain);
end
logLchains = NaN(nsamp+nbi, nchain); logLensemb = NaN(1, nsamp*nchain);
for chain = 1:nchain
    % Storage space for realizations of parameters
    Tt = NaN(1,nsamp+nbi);
    To = NaN(1,nsamp+nbi);
    Mt = NaN(1,nsamp+nbi);
    Mo = NaN(1,nsamp+nbi);
    Dt = NaN(1,nsamp+nbi);
    Do = NaN(1,nsamp+nbi);
    taui = NaN(1,nsamp+nbi);
    taue = NaN(1,nsamp+nbi);
    eoi = NaN(1,nsamp+nbi);
    logLdata = NaN(1,nsamp+nbi);
    %
    if verbose; disp(['Working on chain ' num2str(chain)...
            ' out of ' num2str(nchain) '...']); end
    %
    % Initialize the MCMC:
    sim = 1;
    %
    % Initialize growth response parameters:
    paramscurr = struct('T',T,'P',P,'phi',phi,'D',D,'intwindow',intwindow,'dampth',dampth,'substep',substep);
    Tt(sim)   = unifrnd(aT1,bT1);     paramscurr.T1   = Tt(sim);     if cT1;   Tt(:)   = Tt(sim);   end
    To(sim)   = unifrnd(aT2,bT2);     paramscurr.T2   = To(sim);     if cT2;   To(:)   = To(sim);   end
    Mt(sim)   = unifrnd(aM1,bM1);     paramscurr.M1   = Mt(sim);     if cM1;   Mt(:)   = Mt(sim);   end
    Mo(sim)   = unifrnd(aM2,bM2);     paramscurr.M2   = Mo(sim);     if cM2;   Mo(:)   = Mo(sim);   end
    Dt(sim)   = unifrnd(aD1,bD1);     paramscurr.D1   = Dt(sim);     if cD1;   Dt(:)   = Dt(sim);   end
    Do(sim)   = unifrnd(aD2,bD2);     paramscurr.D2   = Do(sim);     if cD2;   Do(:)   = Do(sim);   end
    taui(sim) = unifrnd(ataui,btaui); paramscurr.taui = taui(sim);   if ctaui; taui(:) = taui(sim); end
    taue(sim) = unifrnd(ataue,btaue); paramscurr.taue = taue(sim);   if ctaue; taue(:) = taue(sim); end
    eoi(sim)  = unifrnd(aeoi,beoi);   paramscurr.eoi  = eoi(sim);    if ceoi;  eoi(:)  = eoi(sim);  end
    %
    sim = sim+1;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% The main sampling MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create storage space to hold current values of the parameters to pass
    % to auxiliary sampling functions, and storage for MCMC realizations of these
    % parameters. Then initialize values of parameters.
    if errormod == 0
        sigma2w = NaN(nsamp+nbi,1);
        sigma2w(1) = unifrnd(0,1); % initialize estimate from the prior.
        errorpars = sigma2w(1); % errorpars holds current value of error model parameters.
    elseif errormod == 1
        tau2 = NaN(nsamp+nbi,1);
        phi1 = NaN(nsamp+nbi,1);
        % initialize estimates from the joint prior:
        phi1(1) = unifrnd(0,1);
        tau2(1) = unifrnd(0,1);
        while tau2(1) > 1-phi1(1)^2
            phi1(1) = unifrnd(0,1);
            tau2(1) = unifrnd(0,1);
        end
        % hold current values of error model parameters:
        errorpars(1) = phi1(1); errorpars(2) = tau2(1);
    end
    %%% Calculate the first simulation. Save the middle parameters for reuse.
    paramscurr.M = [];
    paramscurr.gE = [];
    paramscurr.g0 = [];
    paramscurr.gT = [];
    paramscurr.gM = [];
    paramscurr.eD = [];
    paramscurr.gD = [];
    [paramscurr.gcurr,details] = VSLiteHist(1:Nyrs,paramscurr);
    paramscurr.M = details.M;
    paramscurr.gE = details.gE;
    paramscurr.g0 = details.g0;
    paramscurr.gT = details.gT;
    paramscurr.gM = details.gM;
    paramscurr.eD = details.eD;
    paramscurr.gD = details.gD;
    %
    while sim < nsamp+nbi+1
        % estimate each parameter using Gibbs method
        if ~cT1;   [Tt(sim),   paramscurr] = param_U_aux('T1',aT1,bT1,paramscurr,errorpars,RW,gparscalint,{'g0','gT'});       end;
        if ~cT2;   [To(sim),   paramscurr] = param_U_aux('T2',aT2,bT2,paramscurr,errorpars,RW,gparscalint,{'g0','gT'});       end;
        if ~cM1;   [Mt(sim),   paramscurr] = param_U_aux('M1',aM1,bM1,paramscurr,errorpars,RW,gparscalint,{'g0','gM'});       end;
        if ~cM2;   [Mo(sim),   paramscurr] = param_U_aux('M2',aM2,bM2,paramscurr,errorpars,RW,gparscalint,{'g0','gM'});       end;
        if ~cD1;   [Dt(sim),   paramscurr] = param_U_aux('D1',aD1,bD1,paramscurr,errorpars,RW,gparscalint,{'gD'});            end;
        if ~cD2;   [Do(sim),   paramscurr] = param_U_aux('D2',aD2,bD2,paramscurr,errorpars,RW,gparscalint,{'gD'});            end;
        if ~ctaui; [taui(sim), paramscurr] = param_U_aux('taui',ataui,btaui,paramscurr,errorpars,RW,gparscalint,{'gD','eD'}); end;
        if ~ctaue; [taue(sim), paramscurr] = param_U_aux('taue',ataue,btaue,paramscurr,errorpars,RW,gparscalint,{'gD','eD'}); end;
        if ~ceoi;  [eoi(sim),  paramscurr] = param_U_aux('eoi',aeoi,beoi,paramscurr,errorpars,RW,gparscalint,{'gD','eD'});    end;
        % Now draw from error model parameters:
        if errormod == 0
            [errorpars,logLdata(sim)] = ...
                errormodel0_aux(errorpars,paramscurr,RW,eparscalint);
            sigma2w(sim) = errorpars;
        elseif errormod == 1
            [errorpars,logLdata(sim)] = ...
                errormodel1_aux(errorpars,paramscurr,RW,eparscalint);
            phi1(sim) = errorpars(1);
            tau2(sim) = errorpars(2);
        end
        %
        sim = sim+1;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~cT1;   Ttchains   (:,chain) = Tt;   Ttensemb   (1,(chain-1)*nsamp+(1:nsamp)) = Tt   (nbi+1:end); end
    if ~cT2;   Tochains   (:,chain) = To;   Toensemb   (1,(chain-1)*nsamp+(1:nsamp)) = To   (nbi+1:end); end
    if ~cM1;   Mtchains   (:,chain) = Mt;   Mtensemb   (1,(chain-1)*nsamp+(1:nsamp)) = Mt   (nbi+1:end); end
    if ~cM2;   Mochains   (:,chain) = Mo;   Moensemb   (1,(chain-1)*nsamp+(1:nsamp)) = Mo   (nbi+1:end); end
    if ~cD1;   Dtchains   (:,chain) = Dt;   Dtensemb   (1,(chain-1)*nsamp+(1:nsamp)) = Dt   (nbi+1:end); end
    if ~cD2;   Dochains   (:,chain) = Do;   Doensemb   (1,(chain-1)*nsamp+(1:nsamp)) = Do   (nbi+1:end); end
    if ~ctaui; tauichains (:,chain) = taui; tauiensemb (1,(chain-1)*nsamp+(1:nsamp)) = taui (nbi+1:end); end
    if ~ctaue; tauechains (:,chain) = taue; taueensemb (1,(chain-1)*nsamp+(1:nsamp)) = taue (nbi+1:end); end
    if ~ceoi;  eoichains  (:,chain) = eoi;  eoiensemb  (1,(chain-1)*nsamp+(1:nsamp)) = eoi  (nbi+1:end); end
    if errormod == 0
        sig2rwchains(:,chain) = sigma2w; sig2rwensemb(1,(chain-1)*nsamp+(1:nsamp)) = sigma2w(nbi+1:end);
    elseif errormod == 1
        phi1chains(:,chain) = phi1; phi1ensemb(1,(chain-1)*nsamp+(1:nsamp)) = phi1(nbi+1:end);
        tau2chains(:,chain) = tau2; tau2ensemb(1,(chain-1)*nsamp+(1:nsamp)) = tau2(nbi+1:end);
    end
    logLchains(:,chain) = logLdata; logLensemb(1,(chain-1)*nsamp+(1:nsamp)) = logLdata(nbi+1:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESS SAMPLES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assess convergence:
if strcmp(pt_ests,'med')
    fconvgergent = @(ensemb) median(ensemb);
elseif strcmp(pt_ests,'mle')
    mle_ind = find(logLensemb==max(logLensemb));
    if length(mle_ind)>1; mle_ind = mle_ind(1); end
    fconvgergent = @(ensemb) ensemb(mle_ind);
end
if ~cT1;   T1   = fconvgergent(Ttensemb);   else; T1   = aT1;   end
if ~cT2;   T2   = fconvgergent(Toensemb);   else; T2   = aT2;   end
if ~cM1;   M1   = fconvgergent(Mtensemb);   else; M1   = aM1;   end
if ~cM2;   M2   = fconvgergent(Moensemb);   else; M2   = aM2;   end
if ~cD1;   D1   = fconvgergent(Dtensemb);   else; D1   = aD1;   end
if ~cD2;   D2   = fconvgergent(Doensemb);   else; D2   = aD2;   end
if ~ctaui; taui = fconvgergent(tauiensemb); else; taui = ataui; end
if ~ctaue; taue = fconvgergent(taueensemb); else; taue = ataue; end
if ~ceoi;  eoi  = fconvgergent(eoiensemb);  else; eoi  = aeoi;  end
% if errormod == 0
%     sig2rw = fconvgergent(sig2rwensemb);
% elseif errormod == 1
%     phi1hat = fconvgergent(phi1ensemb);
%     tau2hat = fconvgergent(tau2ensemb);
% end
%
Rhats = struct;
if ~cT1;   RhatT1   = gelmanrubin92(nsamp,nbi,Ttchains);   Rhats.T1   = RhatT1;   end
if ~cT2;   RhatT2   = gelmanrubin92(nsamp,nbi,Tochains);   Rhats.T2   = RhatT2;   end
if ~cM1;   RhatM1   = gelmanrubin92(nsamp,nbi,Mtchains);   Rhats.M1   = RhatM1;   end
if ~cM2;   RhatM2   = gelmanrubin92(nsamp,nbi,Mochains);   Rhats.M2   = RhatM2;   end
if ~cD1;   RhatD1   = gelmanrubin92(nsamp,nbi,Dtchains);   Rhats.D1   = RhatD1;   end
if ~cD2;   RhatD2   = gelmanrubin92(nsamp,nbi,Dochains);   Rhats.D2   = RhatD2;   end
if ~ctaui; Rhattaui = gelmanrubin92(nsamp,nbi,tauichains); Rhats.taui = Rhattaui; end
if ~ctaue; Rhattaue = gelmanrubin92(nsamp,nbi,tauechains); Rhats.taue = Rhattaue; end
if ~ceoi;  Rhateoi  = gelmanrubin92(nsamp,nbi,eoichains);  Rhats.eoi  = Rhateoi;  end
if errormod == 0
    Rhatsig2rw = gelmanrubin92(nsamp,nbi,sig2rwchains); Rhats.sig2rw = Rhatsig2rw;
elseif errormod == 1
    Rhatphi1 = gelmanrubin92(nsamp,nbi,phi1chains); Rhats.phi1 = Rhatphi1;
    Rhattau2 = gelmanrubin92(nsamp,nbi,tau2chains); Rhats.tau2 = Rhattau2;
end
%
Rhatfields = fieldnames(Rhats);
nRhats = length(Rhatfields);
Rhatvals = NaN(1,nRhats);
for i = 1:nRhats; Rhatvals(i) = Rhats.(Rhatfields{i}); end
if verbose == 1
    texthead = '    Rhat for ';
    textvals = '             ';
    for i = 1:nRhats
        headitem = sprintf('% 6s', Rhatfields{i});
        len = length(headitem);
        valitem = sprintf(['% ', num2str(len), 's'], sprintf('%1.4f', Rhatvals(i)));
        if i > 1
            headitem = sprintf(' | %s', headitem);
            valitem = sprintf(' | %s', valitem);
        end
        texthead = sprintf('%s%s', texthead, headitem);
        textvals = sprintf('%s%s', textvals, valitem);
    end
    disp(texthead);
    disp(textvals);
end
convwarning = any(abs(Rhatvals-1)>convthresh);
if convwarning
    disp('Gelman and Rubin metric suggests MCMC has not yet converged to within desired threshold;')
    disp('Parameter estimation code should be re-run using a greater number of MCMC iterations.')
    disp('(See ''nsamp'' advanced input option.)')
end
%
if nargout > 0
    varargout{1} = Rhats; varargout{2} = convwarning;
end
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CONDITIONAL PARAMETER SAMPLING SUBROUTINES %%%%%%%%%%%
function [retvalue,retparam] = param_U_aux(paramname,boundlower,boundupper,paramscurr,errorpars,RW,cyrs,recalcparams)
% param_U_aux: Do a single step of Gibbs sampling.
% Inputs:
%   recalcparams: names of params that should be re-calculated when sampling from the proposal distributing
% Outputs:
%   retvalue: the sampled param value
%   retparam: the new parameters struct that is acceptable by "VSLiteHist"
if nargin < 8; recalcparams = {}; end
if isfield(paramscurr, 'gcurr'); gcurr = paramscurr.gcurr;
else; gcurr = VSLiteHist(cyrs,paramscurr);
end
paramsprop = paramscurr;
for i = 1:length(recalcparams); paramsprop.(recalcparams{i}) = []; end
paramsprop.(paramname) = unifrnd(boundlower,boundupper);
[gprop,detailsprop] = VSLiteHist(cyrs,paramsprop);
%
if length(errorpars) == 1 % White noise error model:
    sigma2rw = errorpars;
    expcurr = sum((RW(cyrs)-sqrt(1-sigma2rw)*gcurr).^2);
    expprop = sum((RW(cyrs)-sqrt(1-sigma2rw)*gprop).^2);
    HR = exp(-.5*(expprop-expcurr)/sigma2rw);
elseif length(errorpars) == 2 % AR(1) error model:
    phi1 = errorpars(1); tau2 = errorpars(2);
    sigma2rw = tau2/(1-phi1^2);
    %
    [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
    %
    logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*gprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*gprop);
    logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*gcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*gcurr);
    HR = exp(logLprop-logLcurr);
end

% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    % accept
    retvalue = paramsprop.(paramname);
    for i = 1:length(recalcparams); paramsprop.(recalcparams{i}) = detailsprop.(recalcparams{i}); end
    paramsprop.gcurr = gprop;
    retparam = paramsprop;
else
    % reject
    retvalue = paramscurr.(paramname);
    paramsprop.gcurr = gcurr;
    retparam = paramscurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma2rw,logLdata] = errormodel0_aux(sigma2rwcurr,paramscurr,RW,cyrs)
% RW = vector of observed annual ring widths
% Gterms = vector of terms that sum together to give the simulated raw ring with index for all
% months (rows) and all years (columns)
% SETW 4/5/2013
%
%%%%%%%%%% account for variable integration window:
if isfield(paramscurr, 'gcurr'); gcurr = paramscurr.gcurr;
else; gcurr = VSLiteHist(cyrs,paramscurr);
end
%%%%%%%%%%%%
% % sample proposal from the prior:
% sigma2rwprop = unifrnd(0,1); % restricted to the unit interval since sigma^2_rw = (1)/(1+SNR^2) in this model

% try in terms of a uniform prior on the std dev. of the model error,
% rather than a uniform prior on the variance....
sigma2rwprop = (unifrnd(0,1))^2;

%
% accept or reject?
Nyrs = length(cyrs);
%
logprop = -.5*sum((RW(cyrs)-sqrt(1-sigma2rwprop)*gcurr).^2)/sigma2rwprop;
logcurr = -.5*sum((RW(cyrs)-sqrt(1-sigma2rwcurr)*gcurr).^2)/sigma2rwcurr;
HR = ((sigma2rwcurr/sigma2rwprop)^(Nyrs/2))*exp(logprop-logcurr);
if binornd(1,min(HR,1))==1
    sigma2rw = sigma2rwprop;
    logLdata = logprop-Nyrs/2*log(sigma2rwprop);
else
    sigma2rw = sigma2rwcurr;
    logLdata = logcurr-Nyrs/2*log(sigma2rwcurr);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pars,logLdata] = errormodel1_aux(currpars,paramscurr,RW,cyrs)
% RW = vector of observed annual ring widths
% Gterms = vector of terms that sum together to give the simulated raw ring with index for all
% months (rows) and all years (columns)
% SETW 4/5/2013
%
%%%%%%%%%% account for variable integration window:
if isfield(paramscurr, 'gcurr'); gcurr = paramscurr.gcurr;
else; gcurr = VSLiteHist(cyrs,paramscurr);
end
%%%%%%%%%%%%
% read current values of parameters:
phi1curr = currpars(1);
tau2curr = currpars(2);
% if 0 % sample proposal from the prior:
phi1prop = unifrnd(0,1);
tau2prop = unifrnd(0,1);
while tau2prop > 1-phi1prop^2
    % satisfy conditions for stationarity, causality, and also
    % sigma2_w = tau2/(1-phi1^2) <= 1 since sigma2_w = 1/(1+SNR^2) in this model
    phi1prop = unifrnd(0,1);
    tau2prop = unifrnd(0,1);
end

% try in terms of a uniform prior on the std dev. of the model error,
% rather than a uniform prior on the variance....
% sig_w = unifrnd(0,1);
% phi1prop = unifrnd(0,1);
% tau2prop = (1-phi1prop^2)*sig_w^2;



% else % sample proposal using random walk step from current location:
%     phi1prop = phi1curr + .25*randn;
%     tau2prop = tau2curr + .1*randn;
%     while tau2prop > 1-phi1prop^2 || tau2prop <0 || phi1prop < 0
%         % satisfy conditions for stationarity, causality, and also
%         % sigma2_w = tau2/(1-phi1^2) <= 1 since sigma2_w = 1/(1+SNR^2) in this model
%         phi1prop = phi1curr + .25*randn;
%         tau2prop = tau2curr + .1*randn;
%     end
% end
%
% accept or reject?
Ny = length(cyrs);
%
[iSigprop,detSigprop] = makeAR1covmat(phi1prop,tau2prop,Ny);
[iSigcurr,detSigcurr] = makeAR1covmat(phi1curr,tau2curr,Ny);
alphaprop = sqrt(1-tau2prop/(1-phi1prop^2));
alphacurr = sqrt(1-tau2curr/(1-phi1curr^2));
%
logLprop = -.5*(RW(cyrs)-alphaprop*gcurr)'*iSigprop*(RW(cyrs)-alphaprop*gcurr);
logLcurr = -.5*(RW(cyrs)-alphacurr*gcurr)'*iSigcurr*(RW(cyrs)-alphacurr*gcurr);
HR = sqrt(detSigcurr/detSigprop)*exp(logLprop-logLcurr);

if binornd(1,min(HR,1))==1
    pars(1) = phi1prop;
    pars(2) = tau2prop;
    logLdata = logLprop-log(detSigprop);
else
    pars(1) = phi1curr;
    pars(2) = tau2curr;
    logLdata = logLcurr-log(detSigcurr);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [invSigma,detSigma] = makeAR1covmat(phi1,tau2,N)
%%% [Sigma,invSigma,detSigma] = makeAR1covmat(phi1,tau2,N)
% Make approximate joint covariance matrix, inverse covariace matrix, and covariance matrix
% determinant for N sequential observations that follow the AR(1) model
% X_t = phi1 X_t-1 + eps_t, where eps_t ~ N(0,tau^2)

A = -phi1*eye(N);
superdiag = sub2ind([N N],(1:(N-1))',(2:N)');
A(superdiag) = 1;
% Now Var(A* e) \approx tau2*I, so
% Sigma \approx tau2* inv(A)*inv(A')
% and
% invSigma \approx (1/tau2)* A'*A

%Sigma = tau2*(A\eye(N))*(A'\eye(N));
invSigma = (1/tau2)*(A')*A;
detSigma = (tau2/phi1^2)^N;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rhat] = gelmanrubin92(Nsamp,Nbi,chains)
% Usage: Rhat = gelmanrubin92(Nsamp,Nbi,chains)
% Nsamp = number of iterations in each sample
% Nbi = number to consider "burn-in".
% chains must have dimension Nsamp x N, where N is the number of chains.
% SETW 1/26/2011

% Number of chains:
m = size(chains,2);
% number of non-burn-in iterations:
n = Nsamp-Nbi;
%
Xbar = NaN(1,m);
Xs2 = NaN(1,m);
allX = NaN(1,n*m);
for i = 1:m
    X = chains(:,i);
    % within-chain means of X:
    Xbar(i) = nanmean(X(Nbi+1:length(X)));
    % within-chain variances of X:
    Xs2(i) = nanvar(X(Nbi+1:length(X)));
    allX((i-1)*n+1:i*n) = X(Nbi+1:Nsamp);
end
% mean across chains of mean X in each month:
Xbarbar = mean(Xbar);
%
BX = n*(sum(((Xbar-repmat(Xbarbar,1,m)).^2),2))/(m-1);
%
WX = nanmean(Xs2,2);
%
sig2hatX = (n-1)*WX/n + BX/n; % G&R92 eqn. 3
%
VhatX = sig2hatX + BX/(m*n);
varhatXs2 = var(Xs2,0,2);
%
covhatXs2Xbar2 = sum((Xs2-nanmean(Xs2)).*(Xbar.^2-nanmean(Xs2.^2)))/m; 
covhatXs2Xbar = sum((Xs2-nanmean(Xs2)).*(Xbar-nanmean(Xs2)))/m;
%
covhatXs2Xbar2 = covhatXs2Xbar2';covhatXs2Xbar = covhatXs2Xbar';
%
varhatVhatX = (((n-1)/n)^2)*varhatXs2/m + (((m+1)/(m*n))^2)*2*BX.^2/(m-1) +...
    2*((m+1)*(n-1)/(m*n^2))*n*(covhatXs2Xbar2-2*Xbarbar.*covhatXs2Xbar)/m;
dfX = 2*(VhatX.^2)./varhatVhatX;
%
Rhat = (VhatX./WX).*(dfX./(dfX-2));
end

function [av,bv,cv] = dealprior(priorarr)
% resolve the prior argument
% if prior is a 2D vector, it is used as the proposal sampling range. in
% this case, av and bv are respectively resolved as the lower and upper
% bound. and cv is set to FALSE.
% if prior is a scalar, it is used as the certain parameter and not estima-
% ted. in this case, av and bv are both resolved as the certain parameter
% value. and cv is set to TRUE.
% 
n = length(priorarr);
if n == 1
    av = priorarr;
    bv = priorarr;
    cv = true;
elseif n == 2
    av = priorarr(1);
    bv = priorarr(2);
    cv = false;
else
    throw(MException('VSLiteHist:estimate_prior', 'prior must be a scalar or a 2D vector'));
end
    
end
