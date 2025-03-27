% get LagSig: significance of regressions
function LagSig = get_LagSig(Settings, Lag, LagBoot)

nseg = length(Settings.Regression.segstr);

for segi = 1:nseg
    seg = Settings.Regression.segstr{segi};

    nneu = size(Lag.(seg).Rsq,1);
    nlag = size(Lag.(seg).Rsq,2);
    nboot = size(LagBoot.(seg).Rsq,2);

    LagSig.(seg).ps = nan(nneu,nlag);

    for neui = 1:nneu
        for lagi = 1:nlag
            LagSig.(seg).ps(neui,lagi) = sum(LagBoot.(seg).Rsq(neui,:) > Lag.(seg).Rsq(neui,lagi))/nboot;
        end
    end

    for alphi = 1:length(Settings.alph2do)
        alphstr = ['alpha' num2str(Settings.alph2do(alphi)*100,'%0.2i')];
        LagSig.(seg).(alphstr).sigs = LagSig.(seg).ps < Settings.alph2do(alphi);
    end
end
end