% This code will reproduce all results and figures from the manuscript
% "Neurons in human motor thalamus encode reach kinematics and positional
% errors related to braking" by Tien et al, submitted to Nature 
% Communications in 2025. See Readme.txt for more information about the
% code

% Input is the directory where "SpikesKinData.mat" is stored
% These data are available at the following public repository:
% https://osf.io/3unka/?view_only=4e7e0d5cee0e4dd9a567fa637644b3b9

function Tien_et_al_2025_NatCom_Code(datadir)

% Load the data
disp('Loading SpikesKinData.mat');
load([datadir '\SpikesKinData.mat']);

ndat = length(Data);

% Set this to true to recalculate the processed data fields, otherwise, the
% precalculated data fields stored in "Data" will be used.
recalculate = false;

if recalculate
    % Suppress NaN warning
    warning('off', 'MATLAB:chckxy:IgnoreNaN');
    for dati = 1:ndat
        disp(['Processing Session ' num2str(dati) ' out of ' num2str(ndat)]);
    
        % get Stretch: mean firing rates (FRs) across reaches with time-stretching
        Stretch(dati) = get_Stretch(Settings, Data(dati));
    
        % get Stretch Boot: random mean firing rates to establish baseline for
        % significance testing
        StretchBoot(dati) = get_StretchBoot(Settings, Data(dati));
    
        % get Stretch Sig: significance testing for Stretch FRs
        StretchSig(dati) = get_StretchSig(Settings, Stretch(dati), StretchBoot(dati));
    
        % get Realtime: mean firing rates (FRs) across reaches without time-stretching
        Realtime(dati) = get_Realtime(Settings, Data(dati));
    
        % get Realtime Sig: significance testing for Realtime FRs
        RealtimeSig(dati) = get_RealtimeSig(Settings, Realtime(dati), StretchBoot(dati));
    
        % get Directionality: significance testing for directional tuning of
        % FRs
        Directionality(dati) = get_Directionality(Settings, Data(dati), Stretch(dati));
    
        % get Regression Buffer: kinematic and FR data for regressions
        RegressionBuffer(dati) = get_RegressionBuffer(Settings, Data(dati));
    
        % get Lag: lagged regression results
        Lag(dati) = get_Lag(Settings, RegressionBuffer(dati));
    
        % get LagBoot: regression results for randomly shifted FR sample times
        % to calculate significance of lagged regressions
        LagBoot(dati) = get_LagBoot(Settings, RegressionBuffer(dati));
    
        % get LagSig: significance of regressions
        LagSig(dati) = get_LagSig(Settings, Lag(dati), LagBoot(dati));
    
        % get ShapLag: Shapley decomposition of regressions at all lags
        Shap(dati).Reach = get_ShapLag(Settings, RegressionBuffer(dati), 'Reach');
        Shap(dati).ReachE = get_ShapLag(Settings, RegressionBuffer(dati), 'ReachE');
    
        % get FB Regression Buffer: data buffer for front (reach start) vs back
        % (reach end) windowed regressions
        FBRegressionBuffer(dati) = get_FBRegressionBuffer(Settings, Data(dati));
    
        % get FBL regressions: reach start vs. reach end lagged regressions
        FBL(dati) = get_FBL(Settings, FBRegressionBuffer(dati), RegressionBuffer(dati));
        
        % get FBL Boot: regression results for randomly shifted FR sample times
        % to calculate significance of lagged start vs. end regressions
        FBLBoot(dati) = get_FBLBoot(Settings, FBRegressionBuffer(dati), RegressionBuffer(dati));
    
        % get FBL Sig: significance of start vs. end regressions
        FBLSig(dati) = get_FBLSig(Settings, FBL(dati), FBLBoot(dati));
    end
else
    datnames = {'Stretch', 'StretchBoot', 'StretchSig', 'Realtime', 'RealtimeSig', 'Directionality', 'RegressionBuffer', 'Lag', 'LagBoot', 'LagSig', 'Shap', 'FBRegressionBuffer', 'FBL', 'FBLBoot', 'FBLSig'};
    ndatnames = length(datnames);
    for dati = 1:ndat
        for nami = 1:ndatnames
            eval([datnames{nami} '(dati) = Data(dati).(datnames{nami});']);
        end
    end
end

%% Report recording information
disp('');
disp('***Subjects and Sessions***');
% Number of subjects
disp(['# of subjects: ' num2str(length(SessionInfo.ismale))]);
% Number female
disp(['# of female subjects: ' num2str(sum(~SessionInfo.ismale))]);
% Age mean/sd
disp(['Age: ' num2str(mean(SessionInfo.ages)) ' p/m ' num2str(std(SessionInfo.ages))]);
% How many hemispheres / procedures?
nsurg = length(SessionInfo.procedure);
disp(['# of hemispheres/procedures: ' num2str(length(SessionInfo.procedure))]);

% Session info
nnospike = Settings.Global.nounitcount;
nlowreach = Settings.Global.lowreachcount;
nlowfr = Settings.Global.lowfrcount;
ntotsorted = Settings.Global.allnneu;
ngoodsess = Settings.Global.allnsess;

screenednreachpersess = [];
reachlens = [];
nunitspertrack = [];
mfrs = [];
ntargspersess = [];
for si = 1:ndat
    screenednreachpersess = [screenednreachpersess; length(Data(si).Time.reachstarts)];
    reachlens = [reachlens; Data(si).Time.reachstops - Data(si).Time.reachstarts];

    nneu = size(Data(si).N.SpkID,1);
    for neui = 1:nneu
        starttime = Data(si).Time.reachstarts(1) - Settings.Regression.reachpad;
        endtime = Data(si).Time.reachstops(end) + Settings.Regression.reachpad;
        mfrs = [mfrs; sum(Data(si).N.SpkTimes{neui} > starttime & Data(si).N.SpkTimes{neui} < endtime)/(endtime-starttime)];
    end

    % Get the reaches by targ
    rdirs = Data(si).Reach.reachdirs;
    rdirs(Data(si).Reach.outreach==0) = rdirs(Data(si).Reach.outreach==0)+180; % Flip inreaches
    rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
    rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
    rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360;
    rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360; % For good measure
    rdirs = round(rdirs/45)+1;
    % Fix for the new targdirs reflection
    rdirs = target_idx_swap(rdirs);
    udirs = unique(rdirs);
    ntargspersess = [ntargspersess; length(udirs)];
end

disp(' ');
disp('***Recordings***');

% How many no spikes?
disp(['# rejected for no spikes: ' num2str(nnospike)]);

% How many rejected due to low reaches?
disp(['# rejected for < ' num2str(Settings.reachmin) ' valid reaches: ' num2str(nlowreach)]);

% How many good sessions
disp(['# of good remaining sesssions: ' num2str(ngoodsess)]);

% How many rejected due to low FR?
disp(['# units rejected for < ' num2str(Settings.FRmin) ' Hz: ' num2str(nlowfr)]);

% How many remaining reaches?
disp(['# reaches total post-outlier: ' num2str(sum(screenednreachpersess))]);

% Reach duration mean/sd
disp(['Reach duration (s): ' num2str(mean(reachlens)) ' p/m ' num2str(std(reachlens))]);

% How many isolated units? mean/sd per microelectrode recording?
disp(['# of analyzed units: ' num2str(ntotsorted)]);

% Mean/SD of all FRs?
disp(['Mean FRs: ' num2str(mean(mfrs)) ' p/m ' num2str(std(mfrs))]);

% Lowest and highest depth of recording?
disp(['Lowest recording depth: ' num2str(min(Settings.Global.Locs.depths))])
disp(['Highest recording depth: ' num2str(max(Settings.Global.Locs.depths))])

%% Plot a segment of fingertip speed (Figure 1a)
plot_speed_segment(Settings, Data, 5, 4543);

%% Plot 3D kinematics for a session (Figure 1e)
plot_3DKin(Settings, Data, RegressionBuffer, 60);

%% Plot unit recording locations (Figure 2a)
plot_recording_locations(Settings);

%% Plot spike rasters (Figure 2b,c)
plot_Raster(Settings, Data, 51, 2);
plot_Raster(Settings, Data, 63, 1);

%% Plot mean firing rates by direction (Figure 2d,e)
plot_DirMeanFRs(Settings, Data, Stretch, 51, 2);
plot_DirMeanFRs(Settings, Data, Stretch, 63, 1);

%% Plot Time-Stretched FRs with significance (Figure 3a-e), report peri-reach modulation results
% Also plots recording locations of peri-reach modulated units (Extended Data Figure 5a,b)
RStretch = plot_report_Stretch(Settings, Stretch, StretchSig);

disp(' ');
disp('***Peri-reach Modulation***');

alphi = 2;
disp(['alpha: ' num2str(Settings.alph2do(alphi))]);
% How many units reach modulated, percentage?
disp(['# of peri-reach sig units: ' num2str(RStretch(alphi).nsig) '/' num2str(Settings.Global.allnneu) ' (' num2str(100*RStretch(alphi).nsig/Settings.Global.allnneu) '%)']);

% How many samples were positive out of total mod samples, percentage?
disp(['# of positive sig samples: ' num2str(RStretch(alphi).npossamples) '/' num2str(RStretch(alphi).nsigsamples) ' (' num2str(100*RStretch(alphi).npossamples/RStretch(alphi).nsigsamples) '%)']);

% How many mod units had only positive, only negative or both? percentage?
disp(['# only positive: ' num2str(sum(RStretch(alphi).onlypos)) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*sum(RStretch(alphi).onlypos)/RStretch(alphi).nsig) '%)']);
disp(['# only negative: ' num2str(sum(RStretch(alphi).onlyneg)) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*sum(RStretch(alphi).onlyneg)/RStretch(alphi).nsig) '%)']);
disp(['# positive and negative: ' num2str(sum(RStretch(alphi).posneg)) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*sum(RStretch(alphi).posneg)/RStretch(alphi).nsig) '%)']);

% How many units modulated before reach start, percentage?
disp(['# of units modulated before reach start: ' num2str(RStretch(alphi).segunits(1)) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*RStretch(alphi).segunits(1)/RStretch(alphi).nsig) '%)']);

% How many peri-reach samples were before reach start?
disp(['# of sig samples occuring before reach start: ' num2str(RStretch(alphi).segcount(1)) '/' num2str(RStretch(alphi).nsigsamples) ' (' num2str(100*RStretch(alphi).segcount(1)/RStretch(alphi).nsigsamples) '%)']);

% How many units modulated during reach, percentage?
disp(['# of units modulated during reach: ' num2str(RStretch(alphi).inunits) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*(RStretch(alphi).inunits)/RStretch(alphi).nsig) '%)']);

% How many peri-reach samples were during reach?
disp(['# of sig samples occuring during reach: ' num2str(RStretch(alphi).incount) '/' num2str(RStretch(alphi).nsigsamples) ' (' num2str(100*RStretch(alphi).incount/RStretch(alphi).nsigsamples) '%)']);

% How many units modulated 1st half of reach, percentage?
disp(['# of units modulated 1st half of reach: ' num2str(RStretch(alphi).segunits(2)) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*RStretch(alphi).segunits(2)/RStretch(alphi).nsig) '%)']);

% How many peri-reach samples were 1st half of reach?
disp(['# of sig samples occuring 1st half of reach: ' num2str(RStretch(alphi).segcount(2)) '/' num2str(RStretch(alphi).nsigsamples) ' (' num2str(100*RStretch(alphi).segcount(2)/RStretch(alphi).nsigsamples) '%)']);

% How many units modulated 2nd half of reach, percentage?
disp(['# of units modulated 2nd half of reach: ' num2str(RStretch(alphi).segunits(3)) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*RStretch(alphi).segunits(3)/RStretch(alphi).nsig) '%)']);

% How many peri-reach samples were 2nd half of reach?
disp(['# of sig samples occuring 2nd half of reach: ' num2str(RStretch(alphi).segcount(3)) '/' num2str(RStretch(alphi).nsigsamples) ' (' num2str(100*RStretch(alphi).segcount(3)/RStretch(alphi).nsigsamples) '%)']);

% How many units modulated after reach end, percentage?
disp(['# of units modulated after reach end: ' num2str(RStretch(alphi).segunits(4)) '/' num2str(RStretch(alphi).nsig) ' (' num2str(100*RStretch(alphi).segunits(4)/RStretch(alphi).nsig) '%)']);

% How many peri-reach samples were after reach end?
disp(['# of sig samples occuring after reach end: ' num2str(RStretch(alphi).segcount(4)) '/' num2str(RStretch(alphi).nsigsamples) ' (' num2str(100*RStretch(alphi).segcount(4)/RStretch(alphi).nsigsamples) '%)']);

% p-value reach vs not
disp(['ChiSq p-value inreach vs. outreach: ' num2str(RStretch(alphi).inoutchip)]);

% p-value first half vs. second half
disp(['ChiSq p-value 1st half vs. 2nd half of reach: ' num2str(RStretch(alphi).pairchip(RStretch(alphi).segpairs(:,1)==2 & RStretch(alphi).segpairs(:,2)==3))]);

% Within the reach, how many peri-reach significant samples were in the first and second half?
disp(['% of sig reach samples in first half of reach: ' num2str(100*RStretch(alphi).segcount(2)/RStretch(alphi).incount) '%']);
disp(['% of sig reach samples in second half of reach: ' num2str(100*RStretch(alphi).segcount(3)/RStretch(alphi).incount) '%']);

% p-value pre-reach vs. post-reach
disp(['ChiSq p-value pre-reach vs. post-reach: ' num2str(RStretch(alphi).pairchip(RStretch(alphi).segpairs(:,1)==1 & RStretch(alphi).segpairs(:,2)==4))]);

% time and how many sig mod at peak
disp(['Max concurrent modulated: ' num2str(RStretch(alphi).npeakmod) ' at sample ' num2str(RStretch(alphi).peakmodtime)])

%% Plot a spatial firing rate map (Figure 4a)
plot_PCRF(Settings, Data, RegressionBuffer, 6, 1, 15)

%% Plot Directional tuning (Figure 4b,c), report directional tuning results
% Also plots recording locations of directionally tuned units (Extended Data Figure 5c,d)
RDir = plot_report_Directionality(Settings, Directionality);

disp(' ');
disp('***DIRECTIONALITY***');
% How many targets per session mean/sd
disp(['# targets per session: ' num2str(mean(ntargspersess)) ' p/m ' num2str(std(ntargspersess))]);

alphi = 1;
disp(['alpha: ' num2str(Settings.alph2do(alphi))]);

% How many units reach dir tuned, percentage?
disp(['# of peri-reach sig units: ' num2str(RDir(alphi).nsig) '/' num2str(Settings.Global.allnneu) ' (' num2str(100*RDir(alphi).nsig/Settings.Global.allnneu) '%)']);

% How many units dir tuned before reach start, percentage?
disp(['# of units dir tuned before reach start: ' num2str(RDir(alphi).segunits(1)) '/' num2str(RDir(alphi).nsig) ' (' num2str(100*RDir(alphi).segunits(1)/RDir(alphi).nsig) '%)']);

% How many peri-reach samples were before reach start?
disp(['# of sig samples occuring before reach start: ' num2str(RDir(alphi).segcount(1)) '/' num2str(RDir(alphi).nsigsamples) ' (' num2str(100*RDir(alphi).segcount(1)/RDir(alphi).nsigsamples) '%)']);

% How many units dir tuned during reach, percentage?
disp(['# of units dir tuned during reach: ' num2str(RDir(alphi).inunits) '/' num2str(RDir(alphi).nsig) ' (' num2str(100*(RDir(alphi).inunits)/RDir(alphi).nsig) '%)']);

% How many peri-reach samples were during reach?
disp(['# of sig samples occuring during reach: ' num2str(RDir(alphi).incount) '/' num2str(RDir(alphi).nsigsamples) ' (' num2str(100*RDir(alphi).incount/RDir(alphi).nsigsamples) '%)']);

% How many units dir tuned 1st half of reach, percentage?
disp(['# of units dir tuned 1st half of reach: ' num2str(RDir(alphi).segunits(2)) '/' num2str(RDir(alphi).nsig) ' (' num2str(100*RDir(alphi).segunits(2)/RDir(alphi).nsig) '%)']);

% How many peri-reach samples were 1st half of reach?
disp(['# of sig samples occuring 1st half of reach: ' num2str(RDir(alphi).segcount(2)) '/' num2str(RDir(alphi).nsigsamples) ' (' num2str(100*RDir(alphi).segcount(2)/RDir(alphi).nsigsamples) '%)']);

% How many units dir tuned 2nd half of reach, percentage?
disp(['# of units dir tuned 2nd half of reach: ' num2str(RDir(alphi).segunits(3)) '/' num2str(RDir(alphi).nsig) ' (' num2str(100*RDir(alphi).segunits(3)/RDir(alphi).nsig) '%)']);

% How many peri-reach samples were 2nd half of reach?
disp(['# of sig samples occuring 2nd half of reach: ' num2str(RDir(alphi).segcount(3)) '/' num2str(RDir(alphi).nsigsamples) ' (' num2str(100*RDir(alphi).segcount(3)/RDir(alphi).nsigsamples) '%)']);

% How many units dir tuned after reach end, percentage?
disp(['# of units dir tuned after reach end: ' num2str(RDir(alphi).segunits(4)) '/' num2str(RDir(alphi).nsig) ' (' num2str(100*RDir(alphi).segunits(4)/RDir(alphi).nsig) '%)']);

% How many peri-reach samples were after reach end?
disp(['# of sig samples occuring after reach end: ' num2str(RDir(alphi).segcount(4)) '/' num2str(RDir(alphi).nsigsamples) ' (' num2str(100*RDir(alphi).segcount(4)/RDir(alphi).nsigsamples) '%)']);

% p-value reach vs not
disp(['ChiSq p-value inreach vs. outreach: ' num2str(RDir(alphi).inoutchip)]);

% p-value first half vs. second half
disp(['ChiSq p-value 1st half vs. 2nd half of reach: ' num2str(RDir(alphi).pairchip(RDir(alphi).segpairs(:,1)==2 & RDir(alphi).segpairs(:,2)==3))]);

% Within the reach, how many peri-reach significant samples were in the first and second half?
disp(['% of sig reach samples in first half of reach: ' num2str(100*RDir(alphi).segcount(2)/RDir(alphi).incount) '%)']);
disp(['% of sig reach samples in second half of reach: ' num2str(100*RDir(alphi).segcount(3)/RDir(alphi).incount) '%)']);

% p-value pre-reach vs. post-reach
disp(['ChiSq p-value pre-reach vs. post-reach: ' num2str(RDir(alphi).pairchip(RDir(alphi).segpairs(:,1)==1 & RDir(alphi).segpairs(:,2)==4))]);

%% Plot and report regression results (Figure 5a-h)
% Also plots recording locations of kinematic encoding units (Extended Data Figure 5e,f)
RReg = plot_report_ShapLag(Settings, RegressionBuffer, Lag, LagSig, Shap);

disp(' ');
disp('***Regression***');
segs = {'Reach', 'ReachE'};
alphi = 2;
for segi = 1:2
    seg = segs{segi};
    disp(['alpha: ' num2str(Settings.alph2do(alphi))]);
    disp(['segment: ' seg]);
    disp(['nsig: ' num2str(length(RReg(alphi).(seg).tops))])
    disp(['Mean Rsq at optimal: ' num2str(mean(RReg(alphi).(seg).optimalRsq)) ' p/m ' num2str(std(RReg(alphi).(seg).optimalRsq))]);
    disp('Tops: ');
    for parti = 1:length(RReg(alphi).(seg).parts)
        disp([RReg(alphi).(seg).parts{parti} ': ' num2str(sum(RReg(alphi).(seg).tops==parti)) '/' num2str(length(RReg(alphi).(seg).tops)) '(' num2str(sum(RReg(alphi).(seg).tops==parti)*100/length(RReg(alphi).(seg).tops)) '%)']);
    end
    disp('');
end
disp('');
disp(['# significantly encoding both models: ' num2str(sum(RReg(alphi).sigboth))])
disp(['Mean Adj Rsq Increase when adding error terms: ' num2str(mean(RReg(alphi).ReachE.sigbothAR - RReg(alphi).Reach.sigbothAR))]);

%% Plot peri-reach FRs without time-stretching (Extended Data Figures 1 and 2)
plot_Realtime(Settings, Realtime, RealtimeSig);

%% Plot results from reach start vs. reach end windowed regressions (Extended Data Figure 3)
plot_FBL(Settings, Data, FBRegressionBuffer, RegressionBuffer, FBL, FBLSig)

%% Plot duration dependency of peri-reach modulation and directional tuning (Extended Data Figure 4)
plot_Cutoffs_Reach(Settings, StretchSig);
plot_Cutoffs_Directionality(Settings, Directionality);