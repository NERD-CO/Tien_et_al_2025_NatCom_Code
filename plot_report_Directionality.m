function varargout = plot_report_Directionality(Settings, Directionality)

nalph = length(Settings.alph2do);
npt = Settings.Stretch.nbef + Settings.Global.nbet + Settings.Stretch.naft;
nall = Settings.Global.allnneu;

sigs = nan(nall,npt,nalph);

ticker = 1;

for si = 1:length(Directionality)
    DirSig = Directionality(si);

    nneu = size(DirSig.pkw,1);

    for alphi = 1:nalph
        alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
        sigs(ticker:ticker+nneu-1,:,alphi) = DirSig.(alphstr).sigs;
    end

    ticker = ticker+nneu;
end

% Peri-reach windows
periwin = [false(1,Settings.Stretch.nbef-Settings.Stretch.periwin) true(1,Settings.Stretch.periwin) true(1,Settings.Global.nbet) true(1,Settings.Stretch.periwin) false(1,Settings.Stretch.naft-Settings.Stretch.periwin)];
% These windows are for AFTER perireach window has been applied. Pre 400ms, first
% half of reach, 2nd half of reach, post 400ms
segwin = [true(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) true(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) true(1,Settings.Global.nbet/2) false(1,Settings.Stretch.periwin); ...
    false(1,Settings.Stretch.periwin) false(1,Settings.Global.nbet/2) false(1,Settings.Global.nbet/2) true(1,Settings.Stretch.periwin)];

for alphi = 1:nalph
    alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
    sigi = any(sigs(:,periwin,alphi),2);

    thissigs = sigs(sigi,:,alphi);

    nsig = sum(sigi);

    sumsig = sum(thissigs,1);

    ssig{alphi} = sumsig/nall;

    % Do some stats!
    winsigofem = sigs(sigi,periwin,alphi);

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
    R(alphi).segcount = segcount;
    R(alphi).incount = incount;
    R(alphi).inoutchip = inoutchip;
    R(alphi).pairchip = pairchip;
    R(alphi).segpairs = segpairs;
    R(alphi).winsigofem = winsigofem;
    R(alphi).segwin = segwin;
    R(alphi).segunits = segunits;
    R(alphi).inunits = inunits;

    if alphi == 1
        %% now do locations
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

        figure('Renderer','Painters');
        hold on;
        for typei = 1:Settings.Global.Locs.ntracktypes
            thisloc = typei;
            plot(thisloc*[1 1], [0 13.4], '-k');
            plot(thisloc*ones(length(tracksigdepths{typei}),1)-lwid/2, tracksigdepths{typei}, '.b')
            plot([thisloc-lwid thisloc], median(tracksigdepths{typei})*[1 1], '-b', 'LineWidth', 3)
            plot(thisloc*ones(length(tracknotsigdepths{typei}),1)+lwid/2, tracknotsigdepths{typei}, '.k')
            plot([thisloc thisloc+lwid], median(tracknotsigdepths{typei})*[1 1], '-k', 'LineWidth', 3)
            if trackdepthpvals(typei) < 0.05
                text(thisloc, 14, ['p = ' num2str(trackdepthpvals(typei))], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'r');
            else
                text(thisloc, 14, 'N.S.', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', 'k');
            end
        end

        thisloc = Settings.Global.Locs.ntracktypes+2;
        plot(thisloc*ones(length(allsigdepths),1)-lwid/2, allsigdepths, '.b')
        plot([thisloc-lwid thisloc], median(allsigdepths)*[1 1], '-b', 'LineWidth', 3)
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
        plot([0 Settings.Global.Locs.ntracktypes+1], [totprop totprop], '--b');
        xlim([0 Settings.Global.Locs.ntracktypes+1]);
        ylim([0 1]);
        yticks([0 0.25 0.5 0.75 1]);
        xticks(1:Settings.Global.Locs.ntracktypes);
        xticklabels({'ant', 'ctr', 'post', 'med'});
        text(Settings.Global.Locs.ntracktypes+1, totprop, {'Total';'Proportion'}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'b');
        if isempty(pairsigschi)
            text(mean([1 Settings.Global.Locs.ntracktypes]), 1, 'No Sig. Pairwise Diffs (Chi-sq test)', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');
        else
            for pairi = 1:size(pairsigschi,1)
                plot([pairsigschi(1) pairsigschi(2)], [1 1]+pairi*0.05, '-k');
            end
        end

        text(7/2, 14, 'Blue = Directionally Tuned, Black = Not', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'k');

        ylabel([]);
        yticklabels([]);
        xticklabels([]);
        title('');
    end
end
%%
% Plot number signif at any timepoint - both alphas
figure;
hold on;

plot(ssig{1}, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 3);
plot(ssig{2}, '-k', 'LineWidth', 3);
yl = ylim;
plot([Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+0.5], yl, '--k', 'HandleVisibility', 'off');
plot([Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5], yl, '--k', 'HandleVisibility', 'off');
plot([npt npt-50], [0.03 0.03], '-', 'LineWidth', 5, 'Color', [0.494, .184, .557], 'HandleVisibility', 'off');
%     xticks([Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5]);
xticks([0.5, Settings.Stretch.nbef+0.5-Settings.Stretch.periwin, Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5+Settings.Stretch.periwin, npt+0.5]);
%     xticklabels({'Reach Start', 'Reach Stop'});
xticklabels([]);
yticks([0 0.08 0.16]);
xlim([0.5 npt+0.5]);
set(gca,'TickDir', 'out');
set(gca,'TickLength',[0.035 0.01]);
legend;

sigi = any(sigs(:,periwin,1),2);
sigs1 = sigs(sigi,:,1);
sigs2 = sigs(sigi,:,2);

% Plot out contigsigs segments
tunedsigs1 = sigs1;
tunedsigs2 = sigs2;
ntuned = size(tunedsigs1,1);
firsti = nan(ntuned,1);
for tuni = 1:ntuned
    firsti(tuni) = find(tunedsigs1(tuni,:) == 1, 1, 'first');
end

[~,sorti] = sort(firsti);
tunedsigs1 = tunedsigs1(sorti,:);
tunedsigs1 = flipud(tunedsigs1);

tunedsigs2 = tunedsigs2(sorti,:);
tunedsigs2 = flipud(tunedsigs2);

figure;
hold on;
imagesc(tunedsigs1+tunedsigs2);
% colormap gray
colormap([1 1 1; 0.8 0.8 0.8; 0 0 0])
yl = ylim;
plot([Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+0.5], yl, '--k');
plot([Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5], yl, '--k');
plot([0.5 npt+0.5], [0.5 0.5], '-k');
plot([0.5 npt+0.5], [ntuned+0.5 ntuned+0.5], '-k');
plot([0.5 0.5], [0.5 ntuned+0.5], '-k');
plot([npt+0.5 npt+0.5], [0.5 ntuned+0.5], '-k');
xlim([0.5 npt+0.5]);
ylim([0.5 ntuned+0.5]);
%     xticks([Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5]);
% xticks([0.5, Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5, npt+0.5])
xticks([0.5, Settings.Stretch.nbef+0.5-Settings.Stretch.periwin, Settings.Stretch.nbef+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5, Settings.Stretch.nbef+Settings.Global.nbet+0.5+Settings.Stretch.periwin, npt+0.5])
xticklabels({'', '', 'Reach Start', 'Reach Stop', '', ''});
set(gca, 'TickDir', 'out');
set(gca,'TickLength',[0.035 0.01]);
yticks([]);

plot([0.5, 50.5], [5 5], '-', 'LineWidth', 5, 'Color', [0.494, .184, .557])

varargout{1} = R;