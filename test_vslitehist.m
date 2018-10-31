function varargout = test_vslitehist(site)
% Basic usage:  test_vslite(site)
% Input: 'site' is either 1, 2, 3, or 4, runs VS-Lite at one of four test sites.
%   Site 1: 36.45N, -118.22E, 'ca530'
%   Site 2: 34.17N, -117.12E, 'ca544'
%   Site 3: 39.02N, -122.82E, 'ca615'
%   Site 4: 39.32N, -106.08E, 'co523'
% 
% Or can specify outputs: [trw,gT,gM,gE,M] = test_vslite_v2_3
%
% Monthly input climate data is from the OSU PRISM data product for the years 1895 through 1984. 
% Observed chronologies are from the International Tree Ring Data Bank, but have been standardized
% for the period 1895-1984.
% Parameters used are the modes of those found to be 'optimal' over 100 calibration/verification
% tests in Tolwinski-Ward et al., 2010
% SETW 11/2010
% SETW: Updated 11/2011 to run with updated VSLite version 2.2 
% SETW: Updated 7/2013 to run with VS-Lite v2.3, and to additionally 
% estimate parameters of the model at each site using Bayesian scheme 
% detailed in Tolwinski-Ward et al 2013 Climate of the Past paper.

% If user doesn't specify site, just run the first one.
if nargin == 0; site = 1; end

% % load the test data:
% load vslite_testdata.mat
% 
% % select climate and location data for the chosen test site:
% T = m08clim(site).T;
% P = m08clim(site).P;
% phi = sitecoords(site,1);
% D = zeros(1,size(T,2));

% get data
filename = uigetfile('*.xlsx;*.xls');
[RW,T,P,D] = read_data(filename);
phi = inputdlg('Input the latitute:');
phi = str2double(phi{1});

% estimate the climate response parameters:
disp('Performing Bayesian estimation of VS-Lite parameters for chosen site.')
tic;
% [T1,T2,M1,M2,D1,D2,taui,taue,eoi] = estimate_vslitehist_params(RW,'T',T,'P',P,'D',D,'phi',phi,'nbi',200,'nsamp',2000,'gparpriors','uniform',...
%     'D1priorsupp',-3,'D2priorsupp',1,'tauipriorsupp',5,'tauepriorsupp',100,'eoipriorsupp',1);
[T1,T2,M1,M2,D1,D2,taui,taue,eoi] = estimate_vslitehist_params(RW,'T',T,'P',P,'D',D,'phi',phi,'nbi',200,'nsamp',2000,'gparpriors','uniform');
toc;
% save('vslitehist_params.mat', 'T1','T2','M1','M2','D1','D2','taui','taue','eoi');
% load('vslitehist_params.mat');

% Run VS-Lite.
[trw,details] = VSLiteHist(syear:eyear,struct('phi',phi,'T',T,'P',P,'D',D,...
    'T1',T1,'T2',T2,'M1',M1,'M2',M2,'D1',D1,'D2',D2,'taui',taui,'taue',taue,'eoi',eoi));
gM = details.gM;
gT = details.gT;
% Draw some output.
figure;
set(gcf,'units','normalized','position',[.25 .25 .5 .4])

if exist('worldmap')
    subplot(2,2,2);
    plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
    plot(mean(gT,2),'r');
    title('Mean gT (red) and gM (blue)')
    subplot(2,2,1)
    worldmap([27 50],[-127 -65])
    load coast
    plotm(lat,long,'k')
    hold on; plotm(sitecoords(site,1),sitecoords(site,2),'*r','markersize',8)
    title('Site location')
    subplot(2,1,2);
    plot(syear:eyear,trw,'rx-'); hold on
    plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
    eval(['title(''Simulated RW (red) and observed (black)'')']);
    ylim([min(min(trw),min(trw_obs(:,site))) max(max(trw),max(trw_obs(:,site)))]);
elseif exist('m_proj')
    subplot(2,2,2);
    plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
    plot(mean(gT,2),'r');
    title('Mean gT (red) and gM (blue)')
    subplot(2,2,1)
    m_proj('Equidistant Cylindrical','longitude',[-127 -65],'latitude',[27 50])
    m_coast('color','k')
    m_grid
    m_line(sitecoords(site,2),sitecoords(site,1),'marker','*','color','r','markersize',8)
    title('Site location')
    subplot(2,1,2);
    plot(syear:eyear,trw,'rx-'); hold on
    plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
    eval(['title(''Simulated RW (red) and observed (black)'')']);
    ylim([min(min(trw),min(trw_obs(:,site))) max(max(trw),max(trw_obs(:,site)))]);
else
    subplot(2,1,1);
    plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
    plot(mean(gT,2),'r');
    title('Mean gT (red) and gM (blue)')
    subplot(2,1,2);
    plot(syear:eyear,trw,'rx-'); hold on
    plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
    eval(['title(''Simulated RW (red) and observed (black)'')']);
end

if nargout > 0 
   varargout{1} = trw;
   varargout{2} = gT;
   varargout{3} = gM;
   varargout{4} = gE;
   varargout{5} = M;
end
    