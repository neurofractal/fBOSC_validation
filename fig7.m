%%
% Let's analyse some openly available iEEG data with BOSC, eBOSC and fBOSC.
% Specifically, calculate theta abundance for all electrodes
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Trust Centre for Neuroimaging

% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

%% Start fBOSC
root = fileparts(matlab.desktop.editor.getActiveFilename);
cd(fullfile(root,'..'));

start_fBOSC
cd(fullfile(root,'..','validation'));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load some Resting State ECOG Data from 10 patients
%
% This dataset is fixation_PAC.zip, which can be downloaded from:
% https://searchworks.stanford.edu/view/zk881ps0522
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Change this to the fixation_PAX directory
data_dir = '/Volumes/Robert T5/fixation_PAC';

subject_list = {'bp','cc','hl','jc','jm','jp','ug','wc','wm','zt'};
Fs           = 1000;


for sub = 5:length(subject_list)
    disp(subject_list{sub});
    cd(fullfile(data_dir,'data',subject_list{sub}));
    load([subject_list{sub} '_base.mat']);
    
    % Organise into a Fieldtrip structure
    data_ft         = [];
    data_ft.fsample = 1000;
    data_ft.Fs      = 1000;
    
    % Channel names
    for chan = 1:size(data,2)
        data_ft.label{chan,1} = ['electrode_' num2str(chan)];
    end
    
    data_ft.trial{1} = data';
    
    data_ft.time{1} = [1/Fs:1/Fs:length(data)/Fs];
    
    %% Resample to 200Hz
    cfg = [];
    cfg.resamplefs = 200;
    data_ft = ft_resampledata(cfg,data_ft);
    
    %% Filter the data
    % Spectral Interpolation
    cfg                     = [];
    cfg.channel             = 'all';
    cfg.bsfilter           = 'yes';
    cfg.bsfreq             = [59 61];
    data_ft              = ft_preprocessing(cfg,data_ft);
    
    % High-pass filter at 2Hz to remove low-frequency artefacts
    cfg                     = [];
    cfg.hpfilter            = 'yes';
    cfg.hpfreq              = 0.5;
    data_ft           = ft_preprocessing(cfg,data_ft);
    
    % Low Pass Filter
    cfg                     = [];
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = 98;
    data_ft        = ft_preprocessing(cfg,data_ft);

%     % Plot data
%     cfg = [];
%     cfg.blocksize = 10;
%     cfg.viewmode = 'butterfly';
%     ft_databrowser(cfg,data_ft);
    
    % Calculate power spectra with Welch's method
    [psds, freqs] = pwelch(data_ft.trial{1}', 100, [], [], data_ft.fsample);
    
    figure;
    title(subject_list{sub});
    semilogy(log10(freqs),psds,'k');
    xlim([0.1 1.9]);
    set(gca,'FontSize',20)
    xlabel('Log Frequency (Hz)','FontSize',28);
    ylabel('Log Power (au.)','FontSize',28);
    drawnow;
%     print('zt','-dpng','-r400');
    
%%
    %%%%%%%%%%
    % fBOSC
    %%%%%%%%%%
    %% Set-up fBOSC parameters
    
    % general setup
    cfg.fBOSC.F                 = 2.^[1:.125:6];
    cfg.fBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
    cfg.fBOSC.fsample           = data_ft.fsample;         % current sampling frequency of EEG data
    
    % padding
    cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
    cfg.fBOSC.pad.detection_s   = .1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
    cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)
    
    % fooof parameters - fit with fixed line or allow a knee
    cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
    cfg.fBOSC.fooof.version           = 'python';
    
    % threshold settings
    cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
    cfg.fBOSC.threshold.percentile  = .95;                              % percentile of background fit for power threshold
    
    % episode post-processing
    cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now
    
    % general processing settings
    cfg.fBOSC.channel           = []; % select posterior channels (default: all)
    cfg.fBOSC.trial             = []; % select trials (default: all)
    cfg.fBOSC.trial_background  = []; % select trials for background (default: all)
    
    %% Run fBOSC
    clear fBOSC
    [fBOSC, cfg] = fBOSC_wrapper(cfg, data_ft);
%     
%     %% Plot the Results of the 1/f fit
%     cfg.log_freqs = 1;
%     cfg.plot_old = 0;
%     fBOSC_fooof_plot(cfg,fBOSC)
    
    %% Get mean SNR of theta bursts
    SNR_mean = [];
    for chan = 1:length(data_ft.label)
        T = fBOSC.episodes(fBOSC.episodes.Channel == chan,:);
        T = T(T.FrequencyMean >= 3 & T.FrequencyMean < 7,:);
        SNR_mean(chan) = mean(T.SNRMean);
        clear T
    end
    
    save SNR_mean SNR_mean
    
    %% Get abundance in theta
    freqs       = 2.^[1:.125:5];
    theta_freqs = find(freqs < 7 & freqs > 3);

    theta_abundance = mean(squeeze(fBOSC.abundance_ep(:,:,theta_freqs)),2);
    
    save theta_abundance theta_abundance
    
    %% Clear variables
    clear theta_abundance SNR_mean T fBOSC cfg
    
     %%%%%%%%%%%%%%%%%
    % Now BOSC
    %%%%%%%%%%%%%%%%%
    %% Set-up BOSC parameters
    cfg                         = [];
    % general setup
    cfg.fBOSC.F                 = 2.^[1:.125:6];
    cfg.fBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
    cfg.fBOSC.fsample           = data_ft.fsample;         % current sampling frequency of data
    
    % padding
    cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
    cfg.fBOSC.pad.detection_s   = .1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
    cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)
    
    % fooof parameters - fit with fixed line or allow a knee
    cfg.fBOSC.fooof.aperiodic_mode    = 'old';
    cfg.fBOSC.fooof.version           = 'python';
    
    % threshold settings
    cfg.fBOSC.threshold.duration	= repmat(3, 1, numel(cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
    cfg.fBOSC.threshold.percentile  = .95;                              % percentile of background fit for power threshold
    
    % episode post-processing
    cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now
    
    % general processing settings
    cfg.fBOSC.channel           = []; % select posterior channels (default: all)
    cfg.fBOSC.trial             = []; % select trials (default: all)
    cfg.fBOSC.trial_background  = []; % select trials for background (default: all)
    
    %% Run fBOSC
    
    clear fBOSC
    [fBOSC, cfg] = fBOSC_wrapper(cfg, data_ft);
       
%     %% Plot the Results of the 1/f fit
%     cfg.log_freqs = 1;
%     cfg.plot_old = 0;
%     fBOSC_fooof_plot(cfg,fBOSC)
    
    %% Get mean SNR of theta bursts
    SNR_mean = [];
    for chan = 1:length(data_ft.label)
        T = fBOSC.episodes(fBOSC.episodes.Channel == chan,:);
        T = T(T.FrequencyMean >= 3 & T.FrequencyMean < 7,:);
        SNR_mean(chan) = mean(T.SNRMean);
        clear T
    end
    
    save SNR_mean_BOSC SNR_mean
    
    %% Get abundance in theta
    freqs       = 2.^[1:.125:6];
    theta_freqs = find(freqs < 7 & freqs > 3);

    theta_abundance = mean(squeeze(fBOSC.abundance_ep(:,:,theta_freqs)),2);
    
    save theta_abundance_BOSC theta_abundance
    
    %% Clear variables
    clear theta_abundance SNR_mean T fBOSC cfg
    
    %%%%%%%%%%%%%%%%%
    % Now eBOSC
    %%%%%%%%%%%%%%%%%
    %% Set-up eBOSC parameters
    cfg                         = [];
    % general setup
    cfg.eBOSC.F                 = 2.^[1:.125:6];
    cfg.eBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
    cfg.eBOSC.fsample           = data_ft.fsample;         % current sampling frequency of data
    
    % Don't use alpha freqs for 1/f fit
    cfg.eBOSC.F_fit = cfg.eBOSC.F(~(cfg.eBOSC.F >= 8 & cfg.eBOSC.F <= 12));
    cfg.eBOSC.threshold.excludePeak = [];

    % padding
    cfg.eBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
    cfg.eBOSC.pad.detection_s   = .1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
    cfg.eBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)
    
    % threshold settings
    cfg.eBOSC.threshold.duration	= repmat(3, 1, numel(cfg.eBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
    cfg.eBOSC.threshold.percentile  = .95;                              % percentile of background fit for power threshold
    
    % episode post-processing
    cfg.eBOSC.postproc.use      = 'no';        % Post-processing turned off for now
    
    % general processing settings
    cfg.eBOSC.channel           = []; % select posterior channels (default: all)
    cfg.eBOSC.trial             = []; % select trials (default: all)
    cfg.eBOSC.trial_background  = []; % select trials for background (default: all)
    
    %% Run eBOSC
    
    [eBOSC, cfg] = eBOSC_wrapper(cfg, data_ft);
       
%     %% Plot the Results of the 1/f fit
%     cfg.log_freqs = 1;
%     cfg.plot_old = 0;
%     fBOSC_fooof_plot(cfg,fBOSC)
    
    %% Get mean SNR of theta bursts
    SNR_mean = [];
    for chan = 1:length(data_ft.label)
        T = eBOSC.episodes(eBOSC.episodes.Channel == chan,:);
        T = T(T.FrequencyMean >= 3 & T.FrequencyMean < 7,:);
        SNR_mean(chan) = mean(T.SNRMean);
        clear T
    end
    
    save SNR_mean_eBOSC SNR_mean
    
    %% Get abundance in theta
    freqs       = 2.^[1:.125:6];
    theta_freqs = find(freqs < 7 & freqs > 3);

    theta_abundance = mean(squeeze(eBOSC.abundance_ep(:,:,theta_freqs)),2);
    
    save theta_abundance_eBOSC theta_abundance
    
    %% Clear variables
    clear theta_abundance SNR_mean T eBOSC cfg data data_ft el_codes
    
end

%% 
%%%%%%%%%%%%%%%%%%%%
% Load all data
%%%%%%%%%%%%%%%%%%%%

theta_abundance_all         = [];
theta_abundance_all_BOSC    = [];
theta_abundance_all_eBOSC   = [];

for sub = 1:length(subject_list)
    disp(['Subject: ' num2str(subject_list{sub})]);
    
    cd(fullfile(data_dir,'data',subject_list{sub}));
    
    % fBOSC
    load('theta_abundance.mat');
    theta_abundance_all = vertcat(theta_abundance_all,...
                            theta_abundance);
    clear theta_abundance
    
    % BOSC
    load('theta_abundance_BOSC.mat');
    theta_abundance_all_BOSC = vertcat(theta_abundance_all_BOSC,...
                            theta_abundance);
    clear theta_abundance
    
    % eBOSC
    load('theta_abundance_eBOSC.mat');
    theta_abundance_all_eBOSC = vertcat(theta_abundance_all_eBOSC,...
                            theta_abundance);
    clear theta_abundance    
    
end


%% Plot Figure
figure;
boxplot([theta_abundance_all_BOSC ...
    theta_abundance_all_eBOSC ...
    theta_abundance_all],{'BOSC','eBOSC','fBOSC'});

%% Export for plotting in python
condition = [repmat({'BOSC'},1,length(theta_abundance_all)) ...
    repmat({'eBOSC'},1,length(theta_abundance_all)) ...
    repmat({'FBOSC'},1,length(theta_abundance_all))]';

abundance = [theta_abundance_all_BOSC',theta_abundance_all_eBOSC',theta_abundance_all']';
                          
cd(fullfile(root,'..','validation'));
t = table(abundance,condition);
writetable(t,'abundance_theta_ECOG_resting.csv');   






