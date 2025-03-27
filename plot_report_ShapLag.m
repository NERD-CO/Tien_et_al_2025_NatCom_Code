function varargout = plot_report_ShapLag(Settings, RegressionBuffer, Lag, LagSig, Shap)

segs = {'Reach', 'ReachE'};

ndat = length(RegressionBuffer);

for segi = 1:length(segs)

    seg = segs{segi};

    if strcmp(seg,'ReachE')
        niceorder = [1 8 9 2 3 4 5 6 7];
    else
        niceorder = 1:7;
    end

    nalpha = length(Settings.alph2do);

    nlag = length(Lag(1).lagsamps);
    nall = Settings.Global.allnneu;
    nkin = size(Shap(1).(seg).shap,2);
    allshaps = nan(nall,nkin,nlag);
    allRsq = nan(nall,nlag);
    allAdjRsq = nan(nall,nlag);
    allRs = false(nall,nalpha,nlag);
    ticker = 1;
    for dati = 1:ndat
        nneu = size(Shap(dati).(seg).shap,1);
        allshaps(ticker:ticker+nneu-1,:,:) = Shap(dati).(seg).shap;
        allRsq(ticker:ticker+nneu-1,:) = Lag(dati).(seg).Rsq;

        for Ri = 1:size(Lag(dati).(seg).Rsq,1)
            n = length(RegressionBuffer(dati).Idx.(seg));
            if strcmp(seg,'ReachE')
                k = length(RegressionBuffer(dati).kinstrerr);
            else
                k = length(RegressionBuffer(dati).kinstr);
            end
            allAdjRsq(ticker+Ri-1,:) = 1 - ((1-Lag(dati).(seg).Rsq(Ri,:))*(n-1)/(n-k-1));
        end

        for alphi = 1:nalpha
            alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
            allRs(ticker:ticker+nneu-1,alphi,:) = LagSig(dati).(seg).(alphstr).sigs==1;
        end

        ticker = ticker+nneu;
    end

    alphi = 2;

    Rs = any(squeeze(allRs(:,alphi,:)),2);

    sigRsq = allRsq(Rs,:);
    msRsq = mean(sigRsq);

    optimalRsq = nan(sum(Rs),1);
    for sigi = 1:size(sigRsq,1)
        optimalRsq(sigi) = max(squeeze(sigRsq(sigi,:)));
    end
    R(alphi).(seg).optimalRsq = optimalRsq;

    sigAdjRsq = allAdjRsq(Rs,:);
    msAdjRsq = mean(sigAdjRsq);

    SigCombo{segi,alphi} = Rs;
    RsqCombo{segi,alphi} = allRsq;
    AdjRsqCombo{segi,alphi} = allAdjRsq;
    msRsqCombo{segi,alphi} = msRsq;
    msAdjRsqCombo{segi,alphi} = msAdjRsq;

    if strcmp(seg,'ReachE')
        kinstr = RegressionBuffer(1).kinstrerr;
        parts = Settings.Regression.parts_err;
        partslong = Settings.Regression.parts_err_long;
        pltcols = [.96 .74 .01; 0 0 .66; 0.33 0.33 1; 0.6 0.6 1; 0 0.5 0; 0 0.75 0; 0 1 0; 0.66 0 0; 1 0 0];
    else
        kinstr = RegressionBuffer(1).kinstr;
        parts = Settings.Regression.parts;
        partslong = Settings.Regression.parts_long;
        pltcols = [.96 .74 .01; 0 0 .66; 0.33 0.33 1; 0.6 0.6 1; 0 0.5 0; 0 0.75 0; 0 1 0];
    end
    nparts = length(parts);

    %     pltcols = [1 0 0; 0 0 1; 0 0 0.66; 0 0 0.33; 0 1 0; 0 0.66 0; 0 0.33 0; 0.66 0 0; 0.33 0 0];

    R(alphi).(seg).parts = parts;

    allshapspct = nan(size(allshaps));
    for neui = 1:nall
        for lagi = 1:nlag
            allshapspct(neui,:,lagi) = allshaps(neui,:,lagi)/sum(allshaps(neui,:,lagi));
        end
    end

    sumshaps = nan(nall, nparts, nlag);
    for i = 1:nparts
        for ii = 1:nlag
            sumshaps(:,i,ii) = sum(allshaps(:,contains(kinstr, parts{i}),ii),2);
        end
    end

    sigsumshapsabs = sumshaps(Rs,:,:);

    sumshapspct = nan(size(sumshaps));
    for neui = 1:nall
        for lagi = 1:nlag
            sumshapspct(neui,:,lagi) = sumshaps(neui,:,lagi)/sum(sumshaps(neui,:,lagi));
        end
    end
    sigsumshapspct = sumshapspct(Rs,:,:);

    pairalpha = 0.05;

    % SigsumPct but at optimal lags
    nsig = size(sigRsq,1);
    optisumshapspct = nan(nsig,nparts);
    for sigi = 1:nsig
        [~,rmaxi] = max(squeeze(sigRsq(sigi,:)));
        optisumshapspct(sigi,:) = sigsumshapspct(sigi,:,rmaxi);
    end

    [~, R(alphi).(seg).tops] = max(optisumshapspct,[],2);

    % Sort the kin groups by median
    medianshapspct = median(optisumshapspct,1);
    [medianshapspct, sorti] = sort(medianshapspct,'descend');
    sortpct = optisumshapspct(:,sorti);
    sortcols = pltcols(sorti,:);

    pairs = nchoosek(1:nparts,2);
    npairs = size(pairs,1);
    pairps = nan(npairs,1);
    for pairi = 1:npairs
        pairps(pairi) = signrank(sortpct(:,pairs(pairi,1)), sortpct(:,pairs(pairi,2)));
    end

    [pairs, pairps, pairps  < 0.01];
    pairsig = pairps < 0.01;

    % Find and plot total variance explained by each
    totvexpla = squeeze(sum(sigsumshapsabs,1));
    pctvexpla = nan(nparts,nlag);
    for lagi = 1:nlag
        pctvexpla(:,lagi) = totvexpla(:,lagi)/sum(totvexpla(:,lagi));
    end
    pctvexpla = pctvexpla(niceorder,:);

    scaledvexpla = flipud(repmat(msRsq,[nparts, 1]).*pctvexpla);

    mrsqaccounted = squeeze(mean(sigsumshapsabs(:,niceorder,:),1));
    scaledpctvexpla = pctvexpla.*repmat(msRsq,[size(pctvexpla,1) 1]);

    figure('Position', [360, 198, 840, 840]);
    hold on;
    for parti = 1:nparts
        swarmchart(parti*ones(size(sortpct,1),1),sortpct(:,parti),20,'o','MarkerEdgeColor',sortcols(parti,:),'XJitterWidth',0.7, 'LineWidth', 2);
        plot([parti-0.44 parti+0.44], [medianshapspct(parti), medianshapspct(parti)], 'LineWidth', 3, 'Color', sortcols(parti,:));
    end

    xticks(1:nparts)
    xticklabels([]);
    sortpartslong = partslong(sorti);
    for parti = 1:nparts
        text(parti,0, sortpartslong{parti}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Rotation', 30, 'FontSize', 20);
    end

    testis = nchoosek(1:nparts,2);

    set(gca,'TickDir','out');

    siglines = [1, .92, .84, .76, .68, .6, .52, .44]+0.08;

    % First draw the horizontal lines
    for parti = 1:nparts-1
        thesepairs = pairs(:,1)==parti;
        if ~any(thesepairs)
            continue;
        end
        sigpairs = find(thesepairs & pairsig);
        if isempty(sigpairs)
            continue;
        end
        lastone = max(pairs(sigpairs,2));
        plot([parti lastone], [siglines(parti), siglines(parti)],'-k', 'LineWidth', 3);
        plot([parti parti], [siglines(parti)+0.00395, siglines(parti)-0.0325], '-k', 'LineWidth', 3);
        for ppi = 1:length(sigpairs)
            plot([pairs(sigpairs(ppi),2), pairs(sigpairs(ppi),2)], [siglines(parti)+0.00395, siglines(parti)-0.0325], '-k', 'LineWidth', 3);
        end
        plot(median([parti lastone]),siglines(parti)+0.02,'*k');
    end

    set(gca, 'FontSize', 20)

    yticks([0 0.5 1]);
    ylim([0 1.1])
    xlim([0.4 nparts+0.6])

    set(gca,'Position',[0.18 0.21 0.775 0.63]);
    set(gca,'TickLength',[0.035 0.01]);
    

    figure('Renderer', 'painters');
    hold on;

    ax = gca;
    ax.YAxis(1).Color = 'k';
    thispltcols = pltcols(niceorder,:);
    for parti = 1:nparts
        plot(mrsqaccounted(parti,:), '-', 'Color', thispltcols(parti,:), 'LineWidth', 2);
    end
    xlim([1 nlag]);
    xticks(1:60:241);
    xticklabels(-1:0.5:1);
    ylabel('Mean R^2 accounted for');
    xlabel('Neural lag (s)');
    set(gca,'TickDir', 'out');
    set(gca,'TickLength',[0.035 0.01]);
    yticks(0:0.01:0.02);
    ylim([0 0.025]);
    yl = ylim;

    % Color-coded stems
    figure('Renderer','painters');
    hold on;
    lagtime = Lag(1).lagtime;
    nsig = size(sigRsq,1);
    maxrsq = nan(nsig,1);
    maxi = nan(nsig,1);
    maxp = nan(nsig,1);
    for ini = 1:nsig
        [maxrsq(ini),maxi(ini)] = max(sigRsq(ini,:));
        [~,maxp(ini)] = max(squeeze(sigsumshapsabs(ini,niceorder,maxi(ini))));
    end
    [sortmaxrsq,sorti] = sort(maxrsq,'descend');
    sortmaxi = maxi(sorti);
    sortmaxp = maxp(sorti);
    for ini = 1:nsig
        plot([lagtime(sortmaxi(ini)) lagtime(sortmaxi(ini))], [0 sortmaxrsq(ini)], '-', 'Color', thispltcols(sortmaxp(ini),:));
        plot(lagtime(sortmaxi(ini)), sortmaxrsq(ini), 'o', 'Color', thispltcols(sortmaxp(ini),:), 'LineWidth',1);
    end
    xlim([-1 1]);
    xticks([-1 -0.5 0 0.5 1]);
    ylim([0 0.4]);
    yticks([0 0.2 0.4]);
    set(gca,'TickDir', 'out');
    set(gca,'TickLength',[0.035 0.01]);

end

% 1 is Reach 2 is ReachE
alphi = 2;

maxAdjRsq1 = max(AdjRsqCombo{1,alphi},[],2);
maxAdjRsq2 = max(AdjRsqCombo{2,alphi},[],2);
sigAdjRsq1 = maxAdjRsq1(SigCombo{1,alphi});
sigAdjRsq2 = maxAdjRsq2(SigCombo{2,alphi});


sigboth = SigCombo{1,alphi} & SigCombo{2,alphi};
sum(sigboth);
R(alphi).sigboth = sigboth;

sigbothAR1 = maxAdjRsq1(sigboth);
sigbothAR2 = maxAdjRsq2(sigboth);
R(alphi).Reach.sigbothAR = sigbothAR1;
R(alphi).ReachE.sigbothAR = sigbothAR2;

figure('Renderer', 'painters');
set(gcf, 'Position', [10 10 420 420]);
hold on;
for bothi = 1:sum(sigboth)
    plot([1 2], [sigbothAR1(bothi) sigbothAR2(bothi)], '-', 'LineWidth', 1, 'Color', [0.75 0.75 0.75]);
end

% Plot medians
plot([0.8 1.2], [median(sigAdjRsq1) median(sigAdjRsq1)], '-', 'LineWidth', 5, 'Color', [0.5 0 0.5]);
plot([1.8 2.2], [median(sigAdjRsq2) median(sigAdjRsq2)], '-', 'LineWidth', 5, 'Color', [0.5 0 0.5]);

pdif = ranksum(sigAdjRsq1,sigAdjRsq2);
if pdif < 0.01
    plot([1 2], [0.44 0.44], '-k', 'LineWidth', 1);
    plot([1 1], [0.43 0.44], '-k', 'LineWidth', 1);
    plot([2 2], [0.43 0.44], '-k', 'LineWidth', 1);
    text(1.5, 0.44, '**', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

swarmchart(ones(size(sigAdjRsq1,1)), sigAdjRsq1, '.','MarkerEdgeColor',[0 0 0],'XJitterWidth',0.3, 'SizeData', 50);
swarmchart(2*ones(size(sigAdjRsq2,1)), sigAdjRsq2, '.','MarkerEdgeColor',[0 0 0],'XJitterWidth',0.3, 'SizeData', 50);

text(1, 0.3, ['n=' num2str(size(sigAdjRsq1,1))]);
text(2, 0.3, ['n=' num2str(size(sigAdjRsq2,1))]);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels(segs);
yticks(0:0.1:0.4);
ylim([0 0.5]);
set(gca,'TickDir', 'out');
set(gca,'TickLength',[0.035 0.01]);

% Now do mean Rsq
figure('Renderer', 'painters');
hold on;
[max1, maxi1] = max(msRsqCombo{1,alphi});
[max2, maxi2] = max(msRsqCombo{2,alphi});

plot([lagtime(maxi1), lagtime(maxi1)], [0, max1], 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
plot([lagtime(maxi2), lagtime(maxi2)], [0, max2], 'Color', [0.75 0.75 0.75], 'LineWidth', 1);

plot(lagtime, msRsqCombo{1,alphi}, '-k', 'LineWidth', 2);
plot(lagtime(maxi1), max1, '.k', 'MarkerSize', 20);
text(lagtime(maxi1), max1, num2str(lagtime(maxi1)));

plot(lagtime, msRsqCombo{2,alphi}, ':k', 'LineWidth', 2);
plot(lagtime(maxi2), max2, '.k', 'MarkerSize', 20);
text(lagtime(maxi2), max2, num2str(lagtime(maxi2)));

xticks(-1:0.5:1);
xlim([-1 1]);
set(gca,'TickDir', 'out');
set(gca,'TickLength',[0.035 0.01]);
ylim([0 0.1]);
yticks([0 0.05 0.1]);

sigi = SigCombo{1,2};

% do locations
usedepths = Settings.Global.Locs.depths;

locs.center = 2;
locs.anterior = 1;
locs.posterior = 3;
locs.medial = 4;
orderednames = {'anterior', 'center', 'posterior', 'medial'};

lwid = 0.2;
lgap = 0.05;

for typei = 1:Settings.Global.Locs.ntracktypes
    tracksigdepths{typei} = [];
    tracknotsigdepths{typei} = [];
end

allsigdepths = Settings.Global.Locs.depths(sigi);
allnotsigdepths = Settings.Global.Locs.depths(~sigi);

for alli = 1:length(Settings.Global.Locs.depths)
    xloc = locs.(Settings.Global.Locs.tracks{alli});
    if sigi(alli)
        tracksigdepths{xloc} = [tracksigdepths{xloc}; Settings.Global.Locs.depths(alli)];
    else
        tracknotsigdepths{xloc} = [tracknotsigdepths{xloc}; Settings.Global.Locs.depths(alli)];
    end
end

% do stats
trackdepthpvals = nan(Settings.Global.Locs.ntracktypes,1);
for typei = 1:Settings.Global.Locs.ntracktypes
    if ~(isempty(tracksigdepths{typei}) || isempty(tracknotsigdepths{typei}))
        trackdepthpvals(typei) = ranksum(tracksigdepths{typei}, tracknotsigdepths{typei});
    end
end
depthpvals = ranksum(allsigdepths,allnotsigdepths);

figure('Renderer','Painters');
hold on;
for typei = 1:Settings.Global.Locs.ntracktypes
    thisloc = typei;
    plot(thisloc*[1 1], [0 13.4], '-k');
    plot(thisloc*ones(length(tracksigdepths{typei}),1)-lwid/2, tracksigdepths{typei}, '.m')
    plot([thisloc-lwid thisloc], median(tracksigdepths{typei})*[1 1], '-m', 'LineWidth', 3)
    plot(thisloc*ones(length(tracknotsigdepths{typei}),1)+lwid/2, tracknotsigdepths{typei}, '.k')
    plot([thisloc thisloc+lwid], median(tracknotsigdepths{typei})*[1 1], '-k', 'LineWidth', 3)
    if trackdepthpvals(typei) < 0.05
        text(thisloc, 14, ['p = ' num2str(trackdepthpvals(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'r');
    else
        text(thisloc, 14, 'N.S.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'k');
    end
end

thisloc = Settings.Global.Locs.ntracktypes+2;
plot(thisloc*ones(length(allsigdepths),1)-lwid/2, allsigdepths, '.m')
plot([thisloc-lwid thisloc], median(allsigdepths)*[1 1], '-m', 'LineWidth', 3)
plot(thisloc*ones(length(allnotsigdepths),1)+lwid/2, allnotsigdepths, '.k')
plot([thisloc thisloc+lwid], median(allnotsigdepths)*[1 1], '-k', 'LineWidth', 3)
if depthpvals < 0.05
    text(thisloc, 14, ['p = ' num2str(depthpvals)], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'r');
else
    text(thisloc, 14, 'N.S.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'k');
end

xticks([1:Settings.Global.Locs.ntracktypes Settings.Global.Locs.ntracktypes+2]);
xticklabels({'ant', 'ctr', 'post', 'med', 'ALL'});
grid on;
xlim([0 Settings.Global.Locs.ntracktypes+3]);
yticks([0 5 10]);
if all(usedepths == Settings.Global.Locs.depths)
    ylabel('Depth (absolute, mm)');
else
    ylabel('Depth (relative, mm)');
end

% what about track proportions
totprop = sum(sigi)/length(sigi);
trackprops = nan(Settings.Global.Locs.ntracktypes,1);
trackns = nan(Settings.Global.Locs.ntracktypes,1);
for typei = 1:Settings.Global.Locs.ntracktypes
    trackns(typei) = length(tracksigdepths{typei}) + length(tracknotsigdepths{typei});
    trackprops(typei) = length(tracksigdepths{typei})/trackns(typei);
end

npair = Settings.Global.Locs.ntracktypes*(Settings.Global.Locs.ntracktypes-1)/2;
pairsigsz = [];
pairsigschi = [];
for typei1 = 1:Settings.Global.Locs.ntracktypes-1
    for typei2 = typei1+1:Settings.Global.Locs.ntracktypes
        jointprop = (trackns(typei1)*trackprops(typei1) + trackns(typei2)*trackprops(typei2))/(trackns(typei1)+trackns(typei2));

        % Do it with z test
        tspz = (trackprops(typei1)-trackprops(typei2))/sqrt(jointprop*(1-jointprop)*((1/trackns(typei1))+(1/trackns(typei2))));
        if tspz > 0
            tspp = 2*(1-normcdf(tspz));
        else
            tspp = 2*normcdf(tspz);
        end
        if tspp < 0.05
            pairsigsz = [pairsigs; [typei1 typei2]];
        end

        % Do it with chi-sq test
        n1 = length(tracksigdepths{typei1});
        n2 = length(tracksigdepths{typei2});
        N1 = trackns(typei1);
        N2 = trackns(typei2);
        x1 = [repmat('a',N1,1); repmat('b',N2,1)];
        x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
        [tbl,chi2stat,chipval] = crosstab(x1,x2);
        if chipval < 0.05
            pairsigschi = [pairsigschi; [typei1 typei2]];
        end
    end
end

ylabel([]);
yticklabels([]);
xticklabels([]);
title('');

figure('Renderer','Painters');
hold on;
for typei = 1:Settings.Global.Locs.ntracktypes
    bar(typei, trackprops(typei), 'EdgeColor', 'none', 'FaceColor', 'k')
    text(typei, 0, ['n=' num2str(trackns(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'w');
end
plot([0 Settings.Global.Locs.ntracktypes+1], [totprop totprop], '--m');
xlim([0 Settings.Global.Locs.ntracktypes+1]);
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
xticks(1:Settings.Global.Locs.ntracktypes);
xticklabels({'ant', 'ctr', 'post', 'med'});
text(Settings.Global.Locs.ntracktypes+1, totprop, {'Total';'Proportion'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'm');
if isempty(pairsigschi)
    text(mean([1 Settings.Global.Locs.ntracktypes]), 1, 'No Sig. Pairwise Diffs (Chi-sq test)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
else
    for pairi = 1:size(pairsigschi,1)
        plot([pairsigschi(1) pairsigschi(2)], [1 1]+pairi*0.05, '-k');
    end
end

text(7/2, 14, 'Green = Kinematic Encoding, Black = Not', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');

ylabel([]);
yticklabels([]);
xticklabels([]);
title('');

varargout{1} = R;