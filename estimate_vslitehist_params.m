function [T1,T2,M1,M2,varargout] = estimate_vslitehist_params(RW,varargin)
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
%     'phi':          latitude of site (in degrees N) (required if 'gE' is absent)
%     'gE':           calculated monthly growth response of insolation (required
%                     set if phi is absent)
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
%                     argument list, default used is [0.0 8.5]
%     'T2priorsupp'   " T2 prior. Default is [9.0 20.0]
%     'M1priorsupp'   " M1 prior. Default is [0.01 0.03]
%     'M2priorsupp'   " M2 prior. Default is [0.1 0.5]
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
if nargin > 4 % read in advanced options if user-specified:
    % first fill values in with defaults:
    T = [];
    P = [];
    M = [];
    phi = [];
    gE = [];
    errormod = 0;
    gparscalint = 1:length(RW);
    eparscalint = 1:length(RW);
    pt_ests = 'mle';
    substep = 0;
    intwindow = [0 12];
    nsamp = 1000;
    nbi = 200;
    nchain = 3;
    aT1 = 0; bT1 = 9;
    aT2 = 10; bT2 = 24;
    aM1 = 0; bM1 = .1;
    aM2 = .1; bM2 = .5;
    convthresh = .1;
    verbose = 1;
    % then over-write defaults if user-specified:
    Nvararg = length(varargin);
    for i = 1:Nvararg/2
        namein = varargin{2*(i-1)+1};
        valin = varargin{2*i};
        switch namein
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
            case 'errormod'
                errormod = valin;
            case 'gparscalint'
                gparscalint = valin;
            case 'eparscalint'
                eparscalint = valin;
            case 'errorpars'
                errorpars = valin;
            case 'pt_ests'
                pt_ests = valin;
            case 'substep'
                substep = valin;
            case 'intwindow'
                intwindow = valin;
            case 'nsamp'
                nsamp = valin;
            case 'nbi'
                nbi = valin;
            case 'nchain'
                nchain = valin;
            case 'T1priorsupp'
                aT1 = valin(1); bT1 = valin(2);
            case 'T2priorsupp'
                aT2 = valin(1); bT2 = valin(2);
            case 'M1priorsupp'
                aM1 = valin(1); bM1 = valin(2);
            case 'M2priorsupp'
                aM2 = valin(1); bM2 = valin(2);
            case 'convthresh'
                convthresh = valin;
            case 'verbose'
                verbose = valin;
        end
    end
else % otherwise, read in defaults:
    throw(MException('VSLiteHist:estimate_params', 'no enough params received'));
end
%%% check input params %%%
if isempty(P) && isempty(M); throw(MException('VSLiteHist:estimate_params', 'neither P and M is set')); end
if isempty(phi) && isempty(gE); throw(MException('VSLiteHist:estimate_params', 'neither phi and gE is set')); end
%
% Take zscore of RW data to fulfill assumptions of model error/noise structure
RW = zscore(RW);
Gterms_dist = NaN(size(T));
%
% Compute soil moisture:
Mmax =.76; % maximum soil moisture; v/v
Mmin =.01; % minimum soil moisture; v/v
muth = 5.8; % mu from thornthwaite's Ep scheme
mth = 4.886; % m from thornthwaite's Ep scheme
alpha = .093;
Minit = 200; % initial value for soil moisture; v/v
dr = 1000; % root depth
%
Nyrs = size(T,2);
%
% Read in or compute estimate of soil moisture M:
if isempty(M)
    % then estimate soil moisture from T and P inputs via Leaky Bucket:
    if substep == 1;
        M = leakybucket_submonthly(1,Nyrs,phi,T,P,Mmax,Mmin,alpha,mth,muth,dr,Minit/dr);
    elseif substep == 0
        M = leakybucket_monthly(1,Nyrs,phi,T,P,Mmax,Mmin,alpha,mth,muth,dr,Minit/dr);
    end
end
% Compute monthly growth response to insolation, gE:
if isempty(gE)
    gE = Compute_gE(phi);
end
%
%%%% Now do the MCMC sampling: %%%%%%%%%%%%%
Ttchains = NaN(nsamp+nbi, nchain); Ttensemb = NaN(1, nsamp*nchain);
Tochains = NaN(nsamp+nbi, nchain); Toensemb = NaN(1, nsamp*nchain);
Mtchains = NaN(nsamp+nbi, nchain); Mtensemb = NaN(1, nsamp*nchain);
Mochains = NaN(nsamp+nbi, nchain); Moensemb = NaN(1, nsamp*nchain);
Gterms_distchains = NaN([size(T), nsamp+nbi, nchain]); Gterms_distensemb = NaN([size(T), nsamp*nchain]);
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
    logLdata = NaN(1,nsamp+nbi);
    %
    if verbose; disp(['Working on chain ' num2str(chain)...
            ' out of ' num2str(nchain) '...']); end
    %
    % Initialize the MCMC:
    gT = NaN(size(T));
    gM = NaN(size(M));
    sim = 1;
    %
    % Initialize growth response parameters:
    % Initialize Tt and To with draws from priors:
    Tt(sim) = unifrnd(aT1,bT1);
    To(sim) = unifrnd(aT2,bT2);
    % Initialize Mt and Mo with draws from priors:
    Mt(sim) = unifrnd(aM1,bM1);
    Mo(sim) = unifrnd(aM2,bM2);
    %
    gT(T<Tt(sim)) = 0;
    gT(T>To(sim)) = 1;
    gT(T>Tt(sim)&T<To(sim)) = (T(T>Tt(sim)&T<To(sim))-Tt(sim))/(To(sim)-Tt(sim));
    %
    gM(M<Mt(sim)) = 0;
    gM(M>Mo(sim)) = 1;
    gM(M>Mt(sim)&M<Mo(sim)) = (M(M>Mt(sim)&M<Mo(sim))-Mt(sim))/(Mo(sim)-Mt(sim));
    %
    Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
    %
    sim = sim+1;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% The main sampling MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create storage space to hold current values of the parameters to pass
    % to auxiliary sampling functions, and storage for MCMC realizations of these
    % parameters. Then initialize values of parameters.
    if errormod == 0;
        sigma2w = NaN(nsamp+nbi,1);
        sigma2w(1) = unifrnd(0,1); % initialize estimate from the prior.
        errorpars = sigma2w(1); % errorpars holds current value of error model parameters.
    elseif errormod == 1;
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
    %
    while sim < nsamp+nbi+1
        %
        Tt(sim) = Tt_U_aux(Tt(sim-1),T,To(sim-1),gM,RW',errorpars,...
            gE,Gterms,aT1,bT1,intwindow,gparscalint);
        gT(T<Tt(sim)) = 0;
        gT(T>To(sim-1)) = 1;
        gT(T>Tt(sim)&T<To(sim-1)) = (T(T>Tt(sim)&T<To(sim-1))-Tt(sim))/(To(sim-1)-Tt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        To(sim) = To_U_aux(To(sim-1),T,Tt(sim),gM,RW',errorpars,...
            gE,Gterms,aT2,bT2,intwindow,gparscalint);
        gT(T<Tt(sim)) = 0;
        gT(T>To(sim)) = 1;
        gT(T>Tt(sim)&T<To(sim)) = (T(T>Tt(sim)&T<To(sim))-Tt(sim))/(To(sim)-Tt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        Mt(sim) = Mt_U_aux(Mt(sim-1),M,Mo(sim-1),gT,RW',errorpars,...
            gE,Gterms,aM1,bM1,intwindow,gparscalint);
        gM(M<Mt(sim)) = 0;
        gM(M>Mo(sim-1)) = 1;
        gM(M>Mt(sim)&M<Mo(sim-1)) = (M(M>Mt(sim)&M<Mo(sim-1))-Mt(sim))/(Mo(sim-1)-Mt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        Mo(sim) = Mo_U_aux(Mo(sim-1),M,Mt(sim),gT,RW',errorpars,...
            gE,Gterms,aM2,bM2,intwindow,gparscalint);
        gM(M<Mt(sim)) = 0;
        gM(M>Mo(sim)) = 1;
        gM(M>Mt(sim)&M<Mo(sim)) = (M(M>Mt(sim)&M<Mo(sim))-Mt(sim))/(Mo(sim)-Mt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        Gterms_dist(:,:,sim) = Gterms; 
        %
        % Now draw from error model parameters:
        if errormod == 0
            [errorpars,logLdata(sim)] = ...
                errormodel0_aux(errorpars,RW,Gterms,intwindow,eparscalint);
            sigma2w(sim) = errorpars;
        elseif errormod == 1
            [errorpars,logLdata(sim)] = ...
                errormodel1_aux(errorpars,RW,Gterms,intwindow,eparscalint);
            phi1(sim) = errorpars(1);
            tau2(sim) = errorpars(2);
        end
        %
        sim = sim+1;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    Ttchains(:,chain) = Tt; Ttensemb(1,(chain-1)*nsamp+(1:nsamp)) = Tt(nbi+1:end);
    Tochains(:,chain) = To; Toensemb(1,(chain-1)*nsamp+(1:nsamp)) = To(nbi+1:end);
    Mtchains(:,chain) = Mt; Mtensemb(1,(chain-1)*nsamp+(1:nsamp)) = Mt(nbi+1:end);
    Mochains(:,chain) = Mo; Moensemb(1,(chain-1)*nsamp+(1:nsamp)) = Mo(nbi+1:end);
    Gterms_distchains(:,:,:,chain) = Gterms_dist; Gterms_distensemb(:,:,(chain-1)*nsamp+(1:nsamp)) = Gterms_dist(:,:,nbi+1:end);
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
if strcmp(pt_ests,'med');
    T1 = median(Ttensemb); T2 = median(Toensemb);
    M1 = median(Mtensemb); M2 = median(Moensemb);
    if errormod == 0
        sig2rw = median(sig2rwensemb);
    elseif errormod == 1
        phi1hat = median(phi1ensemb);
        tau2hat = median(tau2ensemb);
    end
elseif strcmp(pt_ests,'mle')
    mle_ind = find(logLensemb==max(logLensemb));
    if length(mle_ind)>1; mle_ind = mle_ind(1); end
    T1 = Ttensemb(mle_ind); T2 = Toensemb(mle_ind);
    M1 = Mtensemb(mle_ind); M2 = Moensemb(mle_ind);
    if errormod == 0
        sig2rw = sig2rwensemb(mle_ind);
    elseif errormod == 1
        phi1hat = phi1ensemb(mle_ind);
        tau2hat = tau2ensemb(mle_ind);
    end
end
%
RhatT1 = gelmanrubin92(nsamp,nbi,Ttchains);
RhatT2 = gelmanrubin92(nsamp,nbi,Tochains);
RhatM1 = gelmanrubin92(nsamp,nbi,Mtchains);
RhatM2 = gelmanrubin92(nsamp,nbi,Mochains);
if errormod == 0
    Rhatsig2rw = gelmanrubin92(nsamp,nbi,sig2rwchains);
elseif errormod == 1
    Rhatphi1 = gelmanrubin92(nsamp,nbi,phi1chains);
    Rhattau2 = gelmanrubin92(nsamp,nbi,tau2chains);
end
%
Rhats = [RhatT1 RhatT2 RhatM1 RhatM2];
if verbose == 1
    if errormod == 0
        Rhats = [Rhats Rhatsig2rw];
        disp('    Rhat for T1, T2, M1, M2, sigma2rw:');
        disp([RhatT1 RhatT2 RhatM1 RhatM2 Rhatsig2rw]);
    elseif errormod == 1
        Rhats = [Rhats Rhatphi1 Rhattau2];
        disp('    Rhat for T1, T2, M1, M2, phi1, tau2:');
        disp([RhatT1 RhatT2 RhatM1 RhatM2 Rhatphi1 Rhattau2]);
    end
end
if any(abs(Rhats-1)>convthresh)
    disp('Gelman and Rubin metric suggests MCMC has not yet converged to within desired threshold;')
    disp('Parameter estimation code should be re-run using a greater number of MCMC iterations.')
    disp('(See ''nsamp'' advanced input option.)')
    convwarning = 1;
else
    convwarning = 0;
end
%
if nargout > 0
    varargout{1} = Ttensemb; varargout{2} = Toensemb;
    varargout{3} = Mtensemb; varargout{4} = Moensemb;
    varargout{5} = Rhats; varargout{6} = convwarning;
    if errormod == 0
        varargout{7} = sig2rw; varargout{8} = sig2rwensemb;
        varargout{9} = Gterms_distensemb;
    elseif errormod == 1
        varargout{7} = phi1hat; varargout{8} = tau2hat;
        varargout{9} = phi1ensemb; varargout{10} = tau2ensemb;
        varargout{11} = Gterms_distensemb;
    end
end
%
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SOIL MOISTURE SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SCALED DAYLENGTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CONDITIONAL PARAMETER SAMPLING SUBROUTINES %%%%%%%%%%%
function [Tt] = Tt_U_aux(Ttcurr,T,To,gM,RW,errorpars,gE,Gterms,att,btt,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% T is the matrix of temperature for all years and all months of the previous simulation.
% att is the lower bound on the support of the uniform prior distribution for Tt
% btt is the upper bound on the support of the uniform prior distribution for Tt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
%
if 1 % Sample from prior as proposal distribution!
    Ttprop = unifrnd(att,btt);
    gTprop = NaN*ones(12,Ny);
    gTprop(T<Ttprop) = 0;
    gTprop(T>To) = 1;
    gTprop(T<To&T>Ttprop) = (T(T<To&T>Ttprop)-Ttprop)/(To-Ttprop);
    gprop = diag(gE)*min(gM,gTprop);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Tt = Ttprop;
else
    Tt = Ttcurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [To] = To_U_aux(Tocurr,T,Tt,gM,RW,errorpars,gE,Gterms,ato,bto,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% T is the matrix of temperature for all years and all months of the previous simulation.
% att is the lower bound on the support of the uniform prior distribution for Tt
% btt is the upper bound on the support of the uniform prior distribution for Tt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Toprop = unifrnd(ato,bto);
    gTprop = NaN*ones(12,Ny);
    gTprop(T<Tt) = 0;
    gTprop(T>Toprop) = 1;
    gTprop(T<Toprop&T>Tt) = (T(T<Toprop&T>Tt)-Tt)/(Toprop-Tt);
    gprop = diag(gE)*min(gM,gTprop);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    To = Toprop;
else
    To = Tocurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mt] = Mt_U_aux(Mtcurr,M,Mo,gT,RW,errorpars,gE,Gterms,amt,bmt,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% M is the matrix of soil moisture for all years and all months of the previous simulation.
% amt is the lower bound on the support of the uniform prior distribution for Mt
% bmt is the upper bound on the support of the uniform prior distribution for Mt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Mtprop = unifrnd(amt,bmt);
    gWprop = NaN*ones(12,Ny);
    gWprop(M<Mtprop) = 0;
    gWprop(M>Mo) = 1;
    gWprop(M<Mo&M>Mtprop) = (M(M<Mo&M>Mtprop)-Mtprop)/(Mo-Mtprop);
    gprop = diag(gE)*min(gWprop,gT);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Mt = Mtprop;
else
    Mt = Mtcurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mo] = Mo_U_aux(Mocurr,M,Mt,gT,RW,errorpars,gE,Gterms,amo,bmo,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% M is the matrix of soil moisture for all years and all months of the previous simulation.
% amt is the lower bound on the support of the uniform prior distribution for Mt
% bmt is the upper bound on the support of the uniform prior distribution for Mt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Moprop = unifrnd(amo,bmo);
    gWprop = NaN*ones(12,Ny);
    gWprop(M<Mt) = 0;
    gWprop(M>Moprop) = 1;
    gWprop(M<Moprop&M>Mt) = (M(M<Moprop&M>Mt)-Mt)/(Moprop-Mt);
    gprop = diag(gE)*min(gWprop,gT);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Mo = Moprop;
else
    Mo = Mocurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma2rw,logLdata] = errormodel0_aux(sigma2rwcurr,RW,Gterms,intwindow,cyrs)
% RW = vector of observed annual ring widths
% Gterms = vector of terms that sum together to give the simulated raw ring with index for all
% months (rows) and all years (columns)
% SETW 4/5/2013
%
%%%%%%%%%% account for variable integration window:
I_0 = intwindow(1); I_f = intwindow(2);
if I_0<0; % if we include part of the previous year in each year's modeled growth:
    startmo = 13+I_0;
    endmo = I_f;
    prevseas = [mean(Gterms(startmo:12,:),2) Gterms(startmo:12,1:end-1)];
    Gterms = Gterms(1:endmo,:);
    Gterms = [prevseas; Gterms];
else % no inclusion of last year's growth conditions in estimates of this year's growth:
    startmo = I_0+1;
    endmo = I_f;
    Gterms = Gterms(startmo:endmo,:);
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
Gamma = squeeze(sum(Gterms));
Gammabar = mean(Gamma);
siggamma = std(Gamma);
%
logprop = -.5*sum((RW(cyrs)-sqrt(1-sigma2rwprop)*(Gamma(cyrs)-Gammabar)/siggamma).^2)/sigma2rwprop;
logcurr = -.5*sum((RW(cyrs)-sqrt(1-sigma2rwcurr)*(Gamma(cyrs)-Gammabar)/siggamma).^2)/sigma2rwcurr;
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
function [pars,logLdata] = errormodel1_aux(currpars,RW,Gterms,intwindow,cyrs)
% RW = vector of observed annual ring widths
% Gterms = vector of terms that sum together to give the simulated raw ring with index for all
% months (rows) and all years (columns)
% SETW 4/5/2013
%
%%%%%%%%%% account for variable integration window:
I_0 = intwindow(1); I_f = intwindow(2);
if I_0<0; % if we include part of the previous year in each year's modeled growth:
    startmo = 13+I_0;
    endmo = I_f;
    prevseas = [mean(Gterms(startmo:12,:),2) Gterms(startmo:12,1:end-1)];
    Gterms = Gterms(1:endmo,:);
    Gterms = [prevseas; Gterms];
else % no inclusion of last year's growth conditions in estimates of this year's growth:
    startmo = I_0+1;
    endmo = I_f;
    Gterms = Gterms(startmo:endmo,:);
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
Gamma = squeeze(sum(Gterms));
Gammabar = mean(Gamma);
siggamma = std(Gamma);
% VS-lite estimate of TRW at current parameter values:
What = ((Gamma(cyrs)-Gammabar)/siggamma)';
%
[iSigprop,detSigprop] = makeAR1covmat(phi1prop,tau2prop,Ny);
[iSigcurr,detSigcurr] = makeAR1covmat(phi1curr,tau2curr,Ny);
alphaprop = sqrt(1-tau2prop/(1-phi1prop^2));
alphacurr = sqrt(1-tau2curr/(1-phi1curr^2));
%
logLprop = -.5*(RW(cyrs)'-alphaprop*What)'*iSigprop*(RW(cyrs)'-alphaprop*What);
logLcurr = -.5*(RW(cyrs)'-alphacurr*What)'*iSigcurr*(RW(cyrs)'-alphacurr*What);
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
muhatX = nanmean(allX,2);
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

