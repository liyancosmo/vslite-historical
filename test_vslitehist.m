function varargout = test_vslitehist(varargin)
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
%
% Input Arguments (passed in name-value pair)
%   'inputfile':    The input file path. If absent, will prompt to open a file.
%   'outputfile':   The output file path If absent, will prompt to open a file.
%   'choice':       The choice of estimating parameters. If absent, will
%                   prompt to choose. Available values are:
%                     1: Presents only
%                     2: Presents & Historical disturbances
%   'doplot':       Whether to plot the tree-ring. Default is TRUE.


varargin = VarArgs(varargin);
inputfile = varargin.get('inputfile','');
outputfile = varargin.get('outputfile','');
choice = varargin.get('choice',-1);
doplot = varargin.get('doplot',true);


%%% get data and choice
if isempty(inputfile)
    [inputfile,inputfilepath] = uigetfile('*.xlsx;*.xls', 'Select the input data');
    inputfile = fullfile(inputfilepath,inputfile);
    [years,RW,phi,T,P,D] = read_data(inputfile);
    % phi = inputdlg('Input the latitute:');
    % phi = str2double(phi{1});
end

if choice <= 0
    CHOICE_PRESENTSONLY = 'Presents only';
    CHOICE_PRESENTHIST = 'Presents & Historical disturbances';
    choicestr = questdlg('Choose a mode','VSLiteHist',CHOICE_PRESENTSONLY,CHOICE_PRESENTHIST,'');
    switch choicestr
        case CHOICE_PRESENTSONLY
            choice = 1;
        case CHOICE_PRESENTHIST
            choice = 2;
    end
end

%%% check the arguments
if choice <= 0; throw(MException('VSLiteHist:test', 'choice is not set')); end

%%% standarize RW
RW = zscore(RW);

% estimate the climate response parameters:
disp('Performing Bayesian estimation of VS-Lite parameters for chosen site.')
tic;
if choice == 1
    D = zeros(1,length(RW));
    [T1,T2,M1,M2,D1,D2,taui,taue,eoi] = estimate_vslitehist_params(RW,'T',T,'P',P,'D',D,'phi',phi,'nbi',200,'nsamp',2000,'gparpriors','uniform',...
        'D1priorsupp',-3,'D2priorsupp',1,'tauipriorsupp',5,'tauepriorsupp',100,'eoipriorsupp',1);
elseif choice == 2
    [T1,T2,M1,M2,D1,D2,taui,taue,eoi] = estimate_vslitehist_params(RW,'T',T,'P',P,'D',D,'phi',phi,'nbi',200,'nsamp',2000,'gparpriors','uniform');
end
toc;
% save('vslitehist_params.mat', 'T1','T2','M1','M2','D1','D2','taui','taue','eoi');
% load('vslitehist_params.mat');

% Run VS-Lite.
[trw,details] = VSLiteHist(years,struct('phi',phi,'T',T,'P',P,'D',D,...
    'T1',T1,'T2',T2,'M1',M1,'M2',M2,'D1',D1,'D2',D2,'taui',taui,'taue',taue,'eoi',eoi));

% % Draw some output.
% gM = details.gM;
% gT = details.gT;
% figure;
% set(gcf,'units','normalized','position',[.25 .25 .5 .4])
% if exist('worldmap')
%     subplot(2,2,2);
%     plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
%     plot(mean(gT,2),'r');
%     title('Mean gT (red) and gM (blue)')
%     subplot(2,2,1)
%     worldmap([27 50],[-127 -65])
%     load coast
%     plotm(lat,long,'k')
%     hold on; plotm(sitecoords(site,1),sitecoords(site,2),'*r','markersize',8)
%     title('Site location')
%     subplot(2,1,2);
%     plot(syear:eyear,trw,'rx-'); hold on
%     plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
%     eval(['title(''Simulated RW (red) and observed (black)'')']);
%     ylim([min(min(trw),min(trw_obs(:,site))) max(max(trw),max(trw_obs(:,site)))]);
% elseif exist('m_proj')
%     subplot(2,2,2);
%     plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
%     plot(mean(gT,2),'r');
%     title('Mean gT (red) and gM (blue)')
%     subplot(2,2,1)
%     m_proj('Equidistant Cylindrical','longitude',[-127 -65],'latitude',[27 50])
%     m_coast('color','k')
%     m_grid
%     m_line(sitecoords(site,2),sitecoords(site,1),'marker','*','color','r','markersize',8)
%     title('Site location')
%     subplot(2,1,2);
%     plot(syear:eyear,trw,'rx-'); hold on
%     plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
%     eval(['title(''Simulated RW (red) and observed (black)'')']);
%     ylim([min(min(trw),min(trw_obs(:,site))) max(max(trw),max(trw_obs(:,site)))]);
% else
%     subplot(2,1,1);
%     plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
%     plot(mean(gT,2),'r');
%     title('Mean gT (red) and gM (blue)')
%     subplot(2,1,2);
%     plot(syear:eyear,trw,'rx-'); hold on
%     plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
%     eval(['title(''Simulated RW (red) and observed (black)'')']);
% end

if isempty(outputfile);
    [outputfile,outputfilepath] = uiputfile('*.xlsx;*.xls', 'Save the output data');
    outputfile = fullfile(outputfilepath,outputfile);
end
if ~isempty(outputfile); write_data(outputfile,years,RW,trw,phi,T,P,D); end

% Draw only plots
if doplot
    figure;
    set(gcf, 'Unit', 'inches');
    set(gcf, 'Position', [3 3 9 3]);
    hold on;
    plot(years, trw, 'r-x');
    plot(years, RW, 'k-*');
    hold off;
    legend('Predicted','Ground truth');
end

if nargout > 0 
   varargout{1} = trw;
   varargout{2} = details;
end
    