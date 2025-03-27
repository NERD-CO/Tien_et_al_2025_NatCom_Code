function plot_Realtime(Settings, Realtime, RealtimeSig)

npt = 4/Settings.Realtime.step;
nall = Settings.Global.allnneu;
nalpha = length(Settings.alph2do);

allmrates = nan(nall, npt);
allnormrates = nan(nall, npt);
allspeeds = [];
allsig = nan(nall,nalpha);
allsigofem = nan(nall, npt, nalpha);
allsighilo = nan(nall,npt,nalpha);

ticker = 1;
for dati = 1:length(Realtime)

    nneu = size(Realtime(dati).reach.rates,3);

    allmrates(ticker:ticker+nneu-1,:) = squeeze(mean(Realtime(dati).reach.rates,2))';
    allspeeds = [allspeeds; Realtime(dati).reach.speeds'];

    for alphi = 1:nalpha
        alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];

        allsig(ticker:ticker+nneu-1,alphi) = RealtimeSig(dati).reach.(alphstr).sigi;
        allsigofem(ticker:ticker+nneu-1,:,alphi) = RealtimeSig(dati).reach.(alphstr).sigofem;
        allsighilo(ticker:ticker+nneu-1,:,alphi) = RealtimeSig(dati).reach.(alphstr).sighilo;
    end

    ticker = ticker+nneu;
end

meanspeed = nanmean(allspeeds,1);
stdspeed = nanstd(allspeeds,[],1);

plotx = 1:npt;
errdown = (meanspeed - stdspeed);
errup = (meanspeed + stdspeed);

for alli = 1:nall
    allnormrates(alli,:) = (allmrates(alli,:) - min(allmrates(alli,:)))/(range(allmrates(alli,:)));
end

alpha = 2;
alphstr = ['alpha' num2str(100*Settings.alph2do(alphi), '%0.2i')];

sigi = any(allsigofem(:,((Settings.Realtime.nbef-Settings.Stretch.periwin):(Settings.Realtime.nbef+Settings.Realtime.nbet*2+Settings.Stretch.periwin))+1,alphi), 2);
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

nsig = sum(sigi);
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

% Do front and back separate
[frontposf,frontposxi] = ksdensity(posstarts(posstarts<=((npt/2))),1:((npt/2)));
[frontnegf,frontnegxi] = ksdensity(negstarts(negstarts<=((npt/2))),1:((npt/2)));
frontscaleposf = frontposf./max([frontposf frontnegf]);
frontscalenegf = frontnegf./max([frontposf frontnegf]);

[backposf,backposxi] = ksdensity(posstarts(posstarts>((npt/2))),((npt/2)+1):npt);
[backnegf,backnegxi] = ksdensity(negstarts(negstarts>((npt/2))),((npt/2)+1):npt);
backscaleposf = backposf./max([backposf backnegf]);
backscalenegf = backnegf./max([backposf backnegf]);


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

[frontposfoff,frontposxioff] = ksdensity(posstops(posstops<=((npt/2))),1:((npt/2)));
[frontnegfoff,frontnegxioff] = ksdensity(negstops(negstops<=((npt/2))),1:((npt/2)));
frontscaleposfoff = frontposfoff./max([frontposfoff frontnegfoff]);
frontscalenegfoff = frontnegfoff./max([frontposfoff frontnegfoff]);

[backposfoff,backposxioff] = ksdensity(posstops(posstops>((npt/2))),((npt/2)+1):npt);
[backnegfoff,backnegxioff] = ksdensity(negstops(negstops>((npt/2))),((npt/2)+1):npt);
backscaleposfoff = backposfoff./max([backposfoff backnegfoff]);
backscalenegfoff = backnegfoff./max([backposfoff backnegfoff]);

%% Make the fancy subplots

figdims = [0 0 700 1000];

% FR Bars plot
figure('Position', figdims); clf; hold on;
set(gcf, 'Renderer', 'painters');

colormap gray

nbreak = 10;
% Insert a couple nan columns to space out allnormrates

normrates = [normrates(:,1:((npt/2))), nan(size(normrates,1),nbreak), normrates(:,((npt/2)+1):end)];

imalpha = ones(size(normrates));
imalpha(isnan(normrates)) = 0;
imagesc(normrates, 'AlphaData', imalpha);
xticks([Settings.Realtime.nbef+0.5-150 Settings.Realtime.nbef+0.5-100 Settings.Realtime.nbef+0.5-50 Settings.Realtime.nbef+0.5 Settings.Realtime.nbef+0.5+50, (npt/2)+Settings.Realtime.nbet+0.5+nbreak-50 (npt/2)+Settings.Realtime.nbet+0.5+nbreak (npt/2)+Settings.Realtime.nbet+0.5+nbreak+50 (npt/2)+Settings.Realtime.nbet+0.5+nbreak+100 (npt/2)+Settings.Realtime.nbet+0.5+nbreak+150]);

xticklabels({});
xlim([0.5 plotx(end)+nbreak+0.5]);
xlim1 = xlim(gca);
ylim([0.5 nsig+0.5]);
ax = gca;
ax.FontSize = 30;
set(gca, 'TickDir', 'out');
yticks([]);
yticklabels([]);
xticklabels([]);
yl = ylim;

% Plot colored boxes around the significantly modulated timepoints. Do
% front and back separate
for neui = 1:nsig
    % do front
    fronti = 1:((npt/2));
    thissig = sigofem(neui,fronti);
    thishilo = sighilo(neui,fronti);
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

    % do back
    backi = ((npt/2)+1):((npt/2)+Settings.Realtime.nbet+Settings.Realtime.naft);
    thissig = sigofem(neui,backi);
    thishilo = sighilo(neui,backi);
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
        plot((npt/2)+nbreak+[starts(boxi)+0.5 starts(boxi)+0.5], [neui-0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
        plot((npt/2)+nbreak+[stops(boxi)+0.5 stops(boxi)+0.5], [neui-0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
        plot((npt/2)+nbreak+[starts(boxi)+0.5 stops(boxi)+0.5], [neui-0.5 neui-0.5], lstyle, 'LineWidth', lwidth);
        plot((npt/2)+nbreak+[starts(boxi)+0.5 stops(boxi)+0.5], [neui+0.5 neui+0.5], lstyle, 'LineWidth', lwidth);
    end
end

set(gca,'TickLength',[0.035 0.01]);
%% Stacked plots

fronti = 1:(Settings.Realtime.nbef+Settings.Realtime.nbet);
backi = (Settings.Realtime.nbef+Settings.Realtime.nbet+1):(Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+Settings.Realtime.naft);

plotfront = 1:(Settings.Realtime.nbef+Settings.Realtime.nbet);
plotback = (Settings.Realtime.nbef+Settings.Realtime.nbet+nbreak+1):(Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+Settings.Realtime.naft+nbreak);

figure('Position', figdims); clf; hold on;

% Speed
sp1 = subplot(4,1,1); hold on;

patch([plotfront, fliplr(plotfront)], [errdown(fronti), fliplr(errup(fronti))], 'k', 'FaceAlpha', 0.25, 'LineStyle', 'none');
patch([plotback, fliplr(plotback)], [errdown(backi), fliplr(errup(backi))], 'k', 'FaceAlpha', 0.25, 'LineStyle', 'none');

plot(plotfront, meanspeed(fronti), '-k', 'LineWidth', 4);
plot(plotback, meanspeed(backi), '-k', 'LineWidth', 4);
xticks([]);
%     ylim1 = ylim(gca);
ylim([0 500]);
%     yl = ylim;
yticks([0 250 500]);
set(gca, 'TickDir', 'out');
ax = gca;
ax.FontSize = 30;
nsamp = Settings.Realtime.nbef + Settings.Realtime.nbet*2 + Settings.Realtime.naft;
xlim([0.5, nsamp+nbreak+0.5]);
yticklabels([]);

plot(nbreak+[nsamp-50 nsamp], [150 150], '-k', 'LineWidth', 8)
set(gca,'TickLength',[0.03 0.01]);
set(gca, 'YAxisLocation', 'right');


% Sigs
sp2 = subplot(4,1,2); hold on;
plot(plotfront, sum(sigofem(:,fronti),1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'black');
plot(plotfront, sum(sighilo(:,fronti)==-1,1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'cyan');
plot(plotfront, sum(sighilo(:,fronti)==1,1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'red');

plot(plotback, sum(sigofem(:,backi),1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'black');
plot(plotback, sum(sighilo(:,backi)==-1,1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'cyan');
plot(plotback, sum(sighilo(:,backi)==1,1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'red');

% xticks([]);
xticks([Settings.Realtime.nbef+0.5-150 Settings.Realtime.nbef+0.5-100 Settings.Realtime.nbef+0.5-50 Settings.Realtime.nbef+0.5 Settings.Realtime.nbef+0.5+50, Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+0.5+nbreak-50 Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+0.5+nbreak Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+0.5+nbreak+50 Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+0.5+nbreak+100 Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+0.5+nbreak+150]);
set(gca, 'XTickLabels', []);

a = get(gca, 'XTickLabels');
set(gca, 'XTickLabels', a, 'FontSize', 30);
set(gca, 'TickDir', 'out');

xlim([0.5, nbreak+nsamp+0.5]);
yl = ylim;
set(gca,'TickLength',[0.03 0.01]);
yticklabels([]);
set(gca,'YAxisLocation', 'right');

% Starts - NotViolin, stacked
tickwid = 0.4;
tickhi = 1;

sp3 = subplot(4,1,3); hold on;

frontupos = upos(upos<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontpgroup = pgroup(upos<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontnps = length(frontupos);
for posi = 1:frontnps
    patch(frontupos(posi)*ones(1,4)+(tickwid*frontpgroup(posi))*[-0.5 -0.5 0.5 0.5], tickhi*[0 0.5 0.5 0], 'r', 'EdgeColor', 'none')
end
frontuneg = uneg(uneg<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontngroup = ngroup(uneg<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontnns = length(frontuneg);
for negi = 1:frontnns
    patch(frontuneg(negi)*ones(1,4)+(tickwid*frontngroup(negi))*[-0.5 -0.5 0.5 0.5], tickhi*[-0.5 0 0 -0.5], 'c', 'EdgeColor', 'none')
end

backupos = upos(upos>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backpgroup = pgroup(upos>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backnps = length(backupos);
for posi = 1:backnps
    patch(nbreak+backupos(posi)*ones(1,4)+(tickwid*backpgroup(posi))*[-0.5 -0.5 0.5 0.5], tickhi*[0 0.5 0.5 0], 'r', 'EdgeColor', 'none')
end
backuneg = uneg(uneg>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backngroup = ngroup(uneg>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backnns = length(backuneg);
for negi = 1:backnns
    patch(nbreak+backuneg(negi)*ones(1,4)+(tickwid*backngroup(negi))*[-0.5 -0.5 0.5 0.5], tickhi*[-0.5 0 0 -0.5], 'c', 'EdgeColor', 'none')
end

plot(frontposxi, frontscaleposf + tickhi/2, 'r', 'LineWidth', 3);
plot(frontnegxi, frontscalenegf + tickhi/2, 'c', 'LineWidth', 3);

plot(nbreak+backposxi, backscaleposf + tickhi/2, 'r', 'LineWidth', 3);
plot(nbreak+backnegxi, backscalenegf + tickhi/2, 'c', 'LineWidth', 3);

ylim([-tickhi/2, 1+tickhi/2]);
xlim([0.5, nbreak+nsamp+0.5]);
xl = xlim;
yl = ylim;
xticks([]);

axis off

% Stops - NotViolin, stacked
sp4 = subplot(4,1,4); hold on;

frontuposoff = uposoff(uposoff<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontpgroupoff = pgroupoff(uposoff<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontnps = length(frontuposoff);
for posi = 1:frontnps
    patch(frontuposoff(posi)*ones(1,4)+(tickwid*frontpgroupoff(posi))*[-0.5 -0.5 0.5 0.5], tickhi*[0 0.5 0.5 0], 'r', 'EdgeColor', 'none')
end
frontunegoff = unegoff(unegoff<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontngroupoff = ngroupoff(unegoff<=(Settings.Realtime.nbef+Settings.Realtime.nbet));
frontnns = length(frontunegoff);
for negi = 1:frontnns
    patch(frontunegoff(negi)*ones(1,4)+(tickwid*frontngroupoff(negi))*[-0.5 -0.5 0.5 0.5], tickhi*[-0.5 0 0 -0.5], 'c', 'EdgeColor', 'none')
end

backuposoff = uposoff(uposoff>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backpgroupoff = pgroupoff(uposoff>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backnps = length(backuposoff);
for posi = 1:backnps
    patch(nbreak+backuposoff(posi)*ones(1,4)+(tickwid*backpgroupoff(posi))*[-0.5 -0.5 0.5 0.5], tickhi*[0 0.5 0.5 0], 'r', 'EdgeColor', 'none')
end
backunegoff = unegoff(unegoff>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backngroupoff = ngroupoff(unegoff>(Settings.Realtime.nbef+Settings.Realtime.nbet));
backnns = length(backunegoff);
for negi = 1:backnns
    patch(nbreak+backunegoff(negi)*ones(1,4)+(tickwid*backngroupoff(negi))*[-0.5 -0.5 0.5 0.5], tickhi*[-0.5 0 0 -0.5], 'c', 'EdgeColor', 'none')
end

plot(frontposxioff, frontscaleposfoff + tickhi/2, 'r', 'LineWidth', 3);
plot(frontnegxioff, frontscalenegfoff + tickhi/2, 'c', 'LineWidth', 3);

plot(nbreak+backposxioff, backscaleposfoff + tickhi/2, 'r', 'LineWidth', 3);
plot(nbreak+backnegxioff, backscalenegfoff + tickhi/2, 'c', 'LineWidth', 3);

ylim([-tickhi/2, 1+tickhi/2]);
xlim([0.5, nbreak+nsamp+0.5]);
xl = xlim;
yl = ylim;
xticks([]);

axis off



%% And now do the multi mini

anames = fieldnames(Realtime);
anames = anames(~strcmp(anames,'reach'));

for ai = 1:length(anames)
    aistr = anames{ai};
    if any(strcmp(aistr,{'stop', 'start'}))
        continue;
    end

    npt = size(Realtime(1).(aistr).rates,1);
    nall = Settings.Global.allnneu;
    nalpha = length(Settings.alph2do);

    allmrates = nan(nall, npt);
    allnormrates = nan(nall, npt);
    allspeeds = [];
    allsig = nan(nall,nalpha);
    allsigofem = nan(nall, npt, nalpha);
    allsighilo = nan(nall,npt,nalpha);

    ticker = 1;
    for dati = 1:length(Realtime)

        nneu = size(Realtime(dati).(aistr).rates,3);

        allmrates(ticker:ticker+nneu-1,:) = squeeze(mean(Realtime(dati).(aistr).rates,2))';
        allspeeds = [allspeeds; Realtime(dati).(aistr).speeds'];

        for alphi = 1:nalpha
            alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];

            allsig(ticker:ticker+nneu-1,alphi) = RealtimeSig(dati).(aistr).(alphstr).sigi;
            allsigofem(ticker:ticker+nneu-1,:,alphi) = RealtimeSig(dati).(aistr).(alphstr).sigofem;
            allsighilo(ticker:ticker+nneu-1,:,alphi) = RealtimeSig(dati).(aistr).(alphstr).sighilo;
        end

        ticker = ticker+nneu;
    end

    meanspeed = nanmean(allspeeds,1);
    stdspeed = nanstd(allspeeds,[],1);

    plotx = 1:npt;
    errdown = (meanspeed - stdspeed);
    errup = (meanspeed + stdspeed);

    for alli = 1:nall
        allnormrates(alli,:) = (allmrates(alli,:) - min(allmrates(alli,:)))/(range(allmrates(alli,:)));
    end

    alphi = 2;
    alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
    %% Mini Stacked plots

    figure('Renderer', 'painters'); hold on;
    posit = get(gcf, 'Position');
    posit(3) = 500;
    set(gcf, 'Position', posit);

    % Speed
    sp1 = subplot(2,1,1); hold on;
    title(aistr)
    patch([plotx, fliplr(plotx)], [errdown, fliplr(errup)], 'k', 'FaceAlpha', 0.25, 'LineStyle', 'none');
    plot(meanspeed, '-k', 'LineWidth', 4);
    xticks([]);
    ylim([0 520]);
    set(gca, 'TickDir', 'out');
    ax = gca;
    ax.FontSize = 30;
    xlim([0.5, Settings.Realtime.naround+0.5]);
    yticks([0 250 500]);
    yticklabels([]);

    plot([5 55], [150 150], '-k', 'LineWidth', 8)
    set(gca,'TickLength',[0.03 0.01]);
    % Sigs
    sp2 = subplot(2,1,2); hold on;
    plot(sum(allsigofem(:,:,alphi),1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'black');
    plot(sum(allsighilo(:,:,alphi)==-1,1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'cyan');
    plot(sum(allsighilo(:,:,alphi)==1,1)./nall, 'LineStyle', '-', 'LineWidth', 3, 'Marker', 'none', 'Color', 'red');

    xticks([0.5 0.5+Settings.Realtime.naround/2 0.5+Settings.Realtime.naround])
    set(gca, 'XTickLabels', []);

    a = get(gca, 'XTickLabels');
    set(gca, 'XTickLabels', a, 'FontSize', 30);
    set(gca, 'TickDir', 'out');

    xlim([0.5, Settings.Realtime.naround+0.5]);

    if alphi == 1
        ylim([0 0.35]);
        yticks([0 .1 .2 .3]);
    else
        ylim([0 0.25]);
        yticks([0 0.1 0.2]);
    end
    yticklabels([]);
    set(gca,'TickLength',[0.03 0.01]);
    yticklabels([]);
    
end