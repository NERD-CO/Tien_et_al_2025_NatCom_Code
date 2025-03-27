function varargout = plot_report_Stretch(Settings, Stretch, StretchSig)

nsess = length(Stretch);

npt = Settings.Stretch.nbef + Settings.Global.nbet + Settings.Stretch.naft;
nall = Settings.Global.allnneu;
nalpha = length(Settings.alph2do);

allmrates = nan(nall, npt);
allnormrates = nan(nall, npt);
allspeeds = [];
allsig = nan(nall,nalpha);
allallofem = nan(nall, npt, nalpha);
allsigofem = nan(nall, npt, nalpha);
allsighilo = nan(nall,npt,nalpha);

ticker = 1;
for si = 1:nsess
    nneu = size(Stretch(si).rates,3);

    allmrates(ticker:ticker+nneu-1,:) = squeeze(mean(Stretch(si).rates,2))';
    allspeeds = [allspeeds; Stretch(si).speeds'];

    for alphi = 1:nalpha
        alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];

        allsig(ticker:ticker+nneu-1,alphi) = StretchSig(si).(alphstr).sigi;
        allallofem(ticker:ticker+nneu-1,:,alphi) = StretchSig(si).(alphstr).ps < Settings.alph2do(alphi);
        allsigofem(ticker:ticker+nneu-1,:,alphi) = StretchSig(si).(alphstr).sigofem;
        allsighilo(ticker:ticker+nneu-1,:,alphi) = StretchSig(si).(alphstr).sighilo;
    end

    ticker = ticker+nneu;
end

allspeeds = allspeeds(~any(allspeeds>800,2),:);
meanspeed = nanmean(allspeeds,1);
stdspeed = nanstd(allspeeds,[],1);

plotx = 1:npt;
errdown = (meanspeed - stdspeed);
errup = (meanspeed + stdspeed);

for alli = 1:nall
    allnormrates(alli,:) = (allmrates(alli,:) - min(allmrates(alli,:)))/(range(allmrates(alli,:)));
end

% Peri-reach windows
periwin = [false(1,Settings.Stretch.nbef-Settings.Stretch.periwin) true(1,Settings.Stretch.periwin) true(1,Settings.Global.nbet) true(1,Settings.Stretch.periwin) false(1,Settings.Stretch.naft-Settings.Stretch.periwin)];
% These windows are for AFTER perireach window has been applied. Pre 400ms, first
% half of reach, 2nd half of reach, post 400ms
segwin = [true(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) true(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) true(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) true(1,Settings.Stretch.periwin)];

alphi = 2;
alphstr = ['alpha' num2str(100*Settings.alph2do(alphi), '%0.2i')];
sigi = any(allsigofem(:,periwin,alphi),2);
nsig = sum(sigi);
normrates = allnormrates(sigi,:);
sigofem = allsigofem(sigi,:,alphi);
sighilo = allsighilo(sigi,:,alphi);
[~,peakt] = max(normrates,[],2);
[~,sorti] = sort(peakt);
sorti = flipud(sorti);
normrates = normrates(sorti,:);
sigofem = sigofem(sorti,:);
sighilo = sighilo(sorti,:);

%% Plot onset times of significant modulation
posstarts = nan(nsig,1);
negstarts = nan(nsig,1);
for neui = 1:nsig
    if any(sighilo(neui,:)==1)
        posstarts(neui) = find(sighilo(neui,:)==1,1,'first');
    end
    if any(sighilo(neui,:)==-1)
        negstarts(neui) = find(sighilo(neui,:)==-1,1,'first');
    end
end

posstarts = posstarts(~isnan(posstarts));
negstarts = negstarts(~isnan(negstarts));

upos = unique(posstarts);
uneg = unique(negstarts);

pgroup = groupcounts(posstarts);
ngroup = groupcounts(negstarts);

nps = length(upos);
nns = length(uneg);

[posf,posxi] = ksdensity(posstarts,1:npt);
[negf,negxi] = ksdensity(negstarts,1:npt);
scaleposf = posf./max([posf negf]);
scalenegf = negf./max([posf negf]);

%% Plot offset times of significant modulation

nsig = sum(sigi);
posstops = nan(nsig,1);
negstops = nan(nsig,1);
for neui = 1:nsig
    if any(sighilo(neui,:)==1)
        posstops(neui) = find(sighilo(neui,:)==1,1,'last');
    end
    if any(sighilo(neui,:)==-1)
        negstops(neui) = find(sighilo(neui,:)==-1,1,'last');
    end
end

posstops = posstops(~isnan(posstops));
negstops = negstops(~isnan(negstops));

uposoff = unique(posstops);
unegoff = unique(negstops);

pgroupoff = groupcounts(posstops);
ngroupoff = groupcounts(negstops);

npsoff = length(uposoff);
nnsoff = length(unegoff);

[posfoff,posxioff] = ksdensity(posstops,1:npt);
[negfoff,negxioff] = ksdensity(negstops,1:npt);
scaleposfoff = posfoff./max([posfoff negfoff]);
scalenegfoff = negfoff./max([posfoff negfoff]);

% Do some stats!
winsigofem = sigofem(:,periwin);
winsighilo = sighilo(:,periwin);

% Count sig samples for chi square test
totcount = sum(sum(winsigofem));
segcount = nan(4,1);
segexpected = nan(4,1);
segunits = nan(4,1);
for segi = 1:4
    segcount(segi) = sum(sum(winsigofem(:,segwin(segi,:))));
    segexpected(segi) = totcount*(sum(segwin(segi,:))/size(segwin,2));
    segunits(segi) = sum(any(winsigofem(:,segwin(segi,:)),2));
end

inunits = sum(any(winsigofem(:,segwin(2,:)),2) | any(winsigofem(:,segwin(3,:)),2));

% Do chisq test for in reach vs out of reach
incount = sum(segcount(2:3));
outcount = sum(segcount([1 4]));
inexpected = sum(segexpected(2:3));
outexpected = sum(segexpected([1 4]));
inoutchi = (((incount-inexpected).^2)./inexpected) + (((outcount-outexpected).^2)./outexpected);
inoutchip = chi2cdf(inoutchi,1,'upper');

% Pairwise chisq tests
npairs = nchoosek(4,2);
pairchip = nan(npairs,1);
segpairs = nan(npairs,2);
ticker = 1;
for segi1 = 1:3
    for segi2 = (segi1+1):4
        count1 = segcount(segi1);
        count2 = segcount(segi2);
        expected1 = (count1+count2)*(sum(segwin(segi1,:)))/(sum(segwin(segi1,:))+sum(segwin(segi2,:)));
        expected2 = (count1+count2)*(sum(segwin(segi2,:)))/(sum(segwin(segi1,:))+sum(segwin(segi2,:)));
        pairchi = (((count1-expected1).^2)./expected1) + (((count2-expected2).^2)./expected2);
        pairchip(ticker) = chi2cdf(pairchi,1,'upper');
        segpairs(ticker,:) = [segi1 segi2];
        ticker = ticker+1;
    end
end

R(alphi).nsig = nsig;
R(alphi).sigi = sigi;
R(alphi).nsigsamples = totcount;
R(alphi).npossamples = sum(sum(winsighilo==1));
R(alphi).onlypos = all(winsighilo >= 0, 2);
R(alphi).onlyneg = all(winsighilo <= 0, 2);
R(alphi).posneg = any(winsighilo == -1, 2) & any(winsighilo == 1, 2);
R(alphi).segcount = segcount;
R(alphi).incount = incount;
R(alphi).inoutchip = inoutchip;
R(alphi).pairchip = pairchip;
R(alphi).segpairs = segpairs;
R(alphi).winsigofem = winsigofem;
R(alphi).segwin = segwin;
R(alphi).segunits = segunits;
R(alphi).inunits = inunits;
[R(alphi).npeakmod, R(alphi).peakmodtime] = max(squeeze(sum(sigofem,1)));

%% Make the fancy subplots
figdims = [0 0 700 1000];

% FR Bars plot
figure('Position', figdims); clf; hold on;
set(gcf, 'Renderer', 'painters');

colormap gray

imalpha = ones(size(normrates));
imalpha(isnan(normrates)) = 0;
imagesc(normrates, 'AlphaData', imalpha);
xticks([0.5, Settings.Stretch.nbef+0.5-Settings.Stretch.periwin, Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5+Settings.Stretch.periwin, npt+0.5])
xticklabels({});

xlim([0.5 plotx(end)+0.5]);
xlim1 = xlim(gca);
ylim([0.5 nsig+0.5]);
ax = gca;
ax.FontSize = 30;
set(gca, 'TickDir', 'out');
yticks([]);
yticklabels([]);
xticklabels([]);
yl = ylim;

% Plot colored boxes around the significantly modulated timepoints
for neui = 1:nsig
    thissig = sigofem(neui,:);
    thishilo = sighilo(neui,:);
    starts = find(diff(thissig)==1);
    if thissig(1)==1
        starts = [0, starts];
    end
    stops = find(diff(thissig)==-1);
    if thissig(end)==1
        stops = [stops, length(thissig)];
    end
    for boxi = 1:length(starts)
        ishilo = thishilo(starts(boxi)+1);
        if ishilo == 1
            lstyle = '-r';
        elseif ishilo == -1
            lstyle = '-c';
        end
        lwidth = 1;
        plot([starts(boxi)+0.5 starts(boxi)+0.5], [neui-0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
        plot([stops(boxi)+0.5 stops(boxi)+0.5], [neui-0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
        plot([starts(boxi)+0.5 stops(boxi)+0.5], [neui-0.5 neui-0.5], lstyle, 'LineWidth', lwidth);
        plot([starts(boxi)+0.5 stops(boxi)+0.5], [neui+0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
    end
end

set(gca,'TickLength',[0.035 0.01]);

%% Stacked plots
figure('Position', figdims); clf; hold on;

% Speed
sp1 = subplot(4,1,1); hold on;
patch([plotx, fliplr(plotx)], [errdown, fliplr(errup)], 'k', 'FaceAlpha', 0.25, 'LineStyle', 'none');
plot(meanspeed, '-k', 'LineWidth', 4);
xticks([]);
% ylim1 = ylim(gca);
% ylim([0, ylim1(2)]);
% STOPPP
% yl = ylim;
set(gca, 'TickDir', 'out');
ax = gca;
ax.FontSize = 30;
xlim([0.5, npt+0.5]);
ylim([0 500]);
yticks([0 250 500]);
yticklabels([]);

plot([npt-65 npt-15], [150 150], '-k', 'LineWidth', 8)
set(gca,'TickLength',[0.03 0.01]);
set(gca, 'YAxisLocation', 'right');

% Sigs
sp2 = subplot(4,1,2); hold on;
ntot = length(sigi);
plot(sum(sigofem,1)./ntot, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'black');
plot(sum(sighilo==-1,1)./ntot, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'cyan');
plot(sum(sighilo==1,1)./ntot, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'red');

xticks([0.5, Settings.Stretch.nbef+0.5-Settings.Stretch.periwin, Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5+Settings.Stretch.periwin, npt+0.5]);
set(gca, 'XTickLabels', []);

a = get(gca, 'XTickLabels');
set(gca, 'XTickLabels', a, 'FontSize', 30);
set(gca, 'TickDir', 'out');

xlim([0.5, npt+0.5]);
set(gca,'TickLength',[0.03 0.01]);
yticklabels([]);
set(gca,'YAxisLocation', 'right');

% Starts - NotViolin, stacked
tickwid = 0.4;
tickhi = 1;

sp4 = subplot(4,1,3); hold on;

for posi = 1:nps
    patch(upos(posi)*ones(1,4)+(tickwid*pgroup(posi))*[-0.5 -0.5 0.5 0.5], tickhi*[0 0.5 0.5 0], 'r', 'EdgeColor', 'none')
end
for negi = 1:nns
    patch(uneg(negi)*ones(1,4)+(tickwid*ngroup(negi))*[-0.5 -0.5 0.5 0.5], tickhi*[-0.5 0 0 -0.5], 'c', 'EdgeColor', 'none')
end

plot(posxi, scaleposf + tickhi/2, 'r', 'LineWidth', 3);
plot(posxi, scalenegf + tickhi/2, 'c', 'LineWidth', 3);

ylim([-tickhi/2, 1+tickhi/2]);
xlim([0.5, npt+0.5]);
xl = xlim;
yl = ylim;
xticks([]);

axis off

% Stops - NotViolin, stacked
sp5 = subplot(4,1,4); hold on;

for posi = 1:npsoff
    patch(uposoff(posi)*ones(1,4)+(tickwid*pgroupoff(posi))*[-0.5 -0.5 0.5 0.5], tickhi*[0 0.5 0.5 0], 'r', 'EdgeColor', 'none')
end
for negi = 1:nnsoff
    patch(unegoff(negi)*ones(1,4)+(tickwid*ngroupoff(negi))*[-0.5 -0.5 0.5 0.5], tickhi*[-0.5 0 0 -0.5], 'c', 'EdgeColor', 'none')
end

plot(posxioff, scaleposfoff + tickhi/2, 'r', 'LineWidth', 3);
plot(posxioff, scalenegfoff + tickhi/2, 'c', 'LineWidth', 3);

ylim([-tickhi/2, 1+tickhi/2]);
xlim([0.5, npt+0.5]);
xl = xlim;
yl = ylim;
xticks([]);

axis off

%% now do locations
% CHANGE THIS FOR RELATIVE OR ABSOLUTE DEPTH
usedepths = Settings.Global.Locs.depths;

locs.center = 2;
locs.anterior = 1;
locs.posterior = 3;
locs.medial = 4;
orderednames = {'anterior', 'center', 'posterior', 'medial'};

lwid = 0.2;
lgap = 0.05;

% Now get where the sigs are and where the nonsigs are
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

% plot those badboys

figure('Renderer','Painters');
hold on;
for typei = 1:Settings.Global.Locs.ntracktypes
    thisloc = typei;
    plot(thisloc*[1 1], [0 13.4], '-k');
    plot(thisloc*ones(length(tracksigdepths{typei}),1)-lwid/2, tracksigdepths{typei}, '.r')
    plot([thisloc-lwid thisloc], median(tracksigdepths{typei})*[1 1], '-r', 'LineWidth', 3)
    plot(thisloc*ones(length(tracknotsigdepths{typei}),1)+lwid/2, tracknotsigdepths{typei}, '.k')
    plot([thisloc thisloc+lwid], median(tracknotsigdepths{typei})*[1 1], '-k', 'LineWidth', 3)
    if trackdepthpvals(typei) < 0.05
        text(thisloc, 14, ['p = ' num2str(trackdepthpvals(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'r');
    else
        text(thisloc, 14, 'N.S.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'k');
    end
end

thisloc = Settings.Global.Locs.ntracktypes+2;
plot(thisloc*ones(length(allsigdepths),1)-lwid/2, allsigdepths, '.r')
plot([thisloc-lwid thisloc], median(allsigdepths)*[1 1], '-r', 'LineWidth', 3)
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
        if tspp < 0.05/npair
            pairsigsz = [pairsigsz; [typei1 typei2]];
        end

        % Do it with chi-sq test
        n1 = length(tracksigdepths{typei1});
        n2 = length(tracksigdepths{typei2});
        N1 = trackns(typei1);
        N2 = trackns(typei2);
        x1 = [repmat('a',N1,1); repmat('b',N2,1)];
        x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
        [tbl,chi2stat,chipval] = crosstab(x1,x2);
        if chipval < 0.05/npair
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
plot([0 Settings.Global.Locs.ntracktypes+1], [totprop totprop], '--r');
xlim([0 Settings.Global.Locs.ntracktypes+1]);
ylim([0 1]);
yticks([0 0.25 0.5 0.75 1]);
xticks(1:Settings.Global.Locs.ntracktypes);
xticklabels({'ant', 'ctr', 'post', 'med'});
text(Settings.Global.Locs.ntracktypes+1, totprop, {'Total';'Proportion'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'r');
if isempty(pairsigschi)
    text(mean([1 Settings.Global.Locs.ntracktypes]), 1, 'No Sig. Pairwise Diffs (Chi-sq test)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
else
    for pairi = 1:size(pairsigschi,1)
        plot([pairsigschi(1) pairsigschi(2)], [1 1]-pairi*0.05, '-k');
    end
end

ylabel([]);
yticklabels([]);
xticklabels([]);
title('');

varargout{1} = R;
