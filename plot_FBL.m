function plot_FBL(Settings, Data, FBRegressionBuffer, RegressionBuffer, FBL, FBLSig)

% Do one example
doi = 22;
R = Data(doi);
FBRB = FBRegressionBuffer(doi);
RB = RegressionBuffer(doi);
gstr = 'SmoothG50';
figure;
hold on;
figpos = get(gcf,'Position');
figpos(3) = figpos(3)*2.5;
set(gcf,'Position',figpos);
doidx = (FBRB.Idx.FrontE(1)):(FBRB.Idx.BackE(end));
ptop = 1.1*max(R.Kin.(gstr).Speed(doidx));
urn = unique(FBRB.Reachnum.FrontE);
for rni = 1:length(urn)
    patch(R.Time.ktime([min(FBRB.Idx.FrontE(FBRB.Reachnum.FrontE==urn(rni))) max(FBRB.Idx.FrontE(FBRB.Reachnum.FrontE==urn(rni))) max(FBRB.Idx.FrontE(FBRB.Reachnum.FrontE==urn(rni))) min(FBRB.Idx.FrontE(FBRB.Reachnum.FrontE==urn(rni)))]), ...
        [0 0 ptop ptop], 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none');
    patch(R.Time.ktime([min(FBRB.Idx.BackE(FBRB.Reachnum.BackE==urn(rni))) max(FBRB.Idx.BackE(FBRB.Reachnum.BackE==urn(rni))) max(FBRB.Idx.BackE(FBRB.Reachnum.BackE==urn(rni))) min(FBRB.Idx.BackE(FBRB.Reachnum.BackE==urn(rni)))]), ...
        [0 0 ptop ptop], 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none');
end
plot(R.Time.ktime(doidx), R.Kin.(gstr).Speed(doidx),'-k', 'LineWidth', 2);
xlim([5171 5185]);
ylim([0 800]);
xticks([]);
yticks([0 400 800]);
yticklabels([]);
patch([5184 5185 5185 5184], [200 200 220 220], 'k', 'LineStyle', 'none');
set(gca,'TickDir', 'out');
set(gca,'TickLength',[0.035/2.5 0.01/2.5]);

nall = Settings.Global.allnneu;

segs = {'FrontE', 'BackE'};

nseg = length(segs);

nlag = length(FBL(1).lagtime);

nalpha = length(Settings.alph2do);

allRsq = nan(nall,nlag,nseg);
allRs = false(nall,nalpha,nlag,nseg);
ticker = 1;
anneu = 0;
for dati = 1:length(Data)

    Lag = FBL(dati);
    LagSig = FBLSig(dati);

    nneu = size(Lag.FrontE.Rsq,1);
    anneu = anneu+nneu;

    for segi = 1:nseg
        seg = segs{segi};
        allRsq(ticker:ticker+nneu-1,:,segi) = Lag.(seg).Rsq;    
        for alphi = 1:nalpha
            alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
            allRs(ticker:ticker+nneu-1,alphi,:,segi) = LagSig.(seg).(alphstr).sigs==1;
        end
    end

    ticker = ticker+nneu;
end

rng(0);

goodrows = find(~isnan(allRsq(:,1,1)));
allRs = allRs(goodrows,:,:,:);
allRsq = allRsq(goodrows,:,:);
zlagi = find(Lag.lagtime==0);
fronti = find(strcmp(segs, 'FrontE'));
backi = find(strcmp(segs, 'BackE'));
for alphi = 1:nalpha
    frontrs = allRs(:,alphi,zlagi,fronti);
    backrs = allRs(:,alphi,zlagi,backi);

    frontrsq = allRsq(frontrs,zlagi,fronti);
    backrsq = allRsq(backrs,zlagi,backi);

    zfronte = sum(allRs(:,alphi,zlagi,fronti));
    zbacke = sum(allRs(:,alphi,zlagi,backi));

    mfront = mean(frontrsq);
    mback = mean(backrsq);

    figure;
    hold on;
    swarmchart(ones(size(frontrsq,1),1), frontrsq, '.','MarkerEdgeColor',[1 0 0],'XJitterWidth',0.3, 'SizeData', 500);
    swarmchart(2*ones(size(backrsq,1),1), backrsq, '.','MarkerEdgeColor',[0 0 1],'XJitterWidth',0.3, 'SizeData', 500);
    xticks([1 2]);
    xlim([0.5 2.5]);
    ylim([0 0.5]);
    yticks([0 0.25 0.5]);
    yticklabels({});
    xticklabels({});
    set(gca,'TickDir', 'out');
    set(gca,'TickLength',[0.035 0.01]);

    % Do a quick permutation test
    nperm = 10000;
    labels = [ones(zfronte,1); 2*ones(zbacke,1)];
    nrs = zfronte+zbacke;
    rsqs = [frontrsq; backrsq];
    obsdiff = mback - mfront;
    testdiff = nan(nrs,1);
    for permi = 1:nperm
        shufflab = labels(randperm(nrs));
        testdiff(permi) = mean(rsqs(shufflab==2)) - mean(rsqs(shufflab==1));
    end
    meandiffpval = min([sum(testdiff > obsdiff), sum(testdiff < obsdiff)])/nperm;

    phat = (zfronte+zbacke)/nall;
    propz = ((zfronte/nall - zbacke/nall)-0)/sqrt(phat*(1-phat)*(2/nall));
    ppval = 2*(1-normcdf(abs(propz)));

end


