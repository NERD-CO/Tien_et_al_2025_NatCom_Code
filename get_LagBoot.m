% get LagBoot: regression results for randomly shifted FR sample times
% to calculate significance of lagged regressions
function LagBoot = get_LagBoot(Settings, B)

nseg = length(Settings.Regression.segstr);

for segi = 1:nseg
    seg = Settings.Regression.segstr{segi};

    if strcmp(seg,'ReachE')
        dokin = [B.kin(B.Idx.(seg),:), B.errxy, B.emag];
    else
        dokin = B.kin(B.Idx.(seg),:);
    end

    npt = size(dokin,1);
    nkin = size(dokin,2);
    nneu = size(B.rates,2);
    nwhole = size(B.time,1);
    
    % Convert Idx numbers to idx 0s and 1s
    idx = false(nwhole,1);
    idx(B.Idx.(seg)) = true;

    % Pick neural max window as max lags away, adjust rates and idx
    % as such
    winstart = B.Idx.Full(1)-Settings.Regression.maxlagsamps;
    winstop = B.Idx.Full(end)+Settings.Regression.maxlagsamps;
    if winstart <= 0
        realwinstart = 1;
    else
        realwinstart = winstart;
    end
    if winstop > nwhole
        realwinstop = nwhole;
    else
        realwinstop = winstop;
    end
    winrates = B.rates(realwinstart:realwinstop,:);
    winidx = idx(realwinstart:realwinstop);
    nwin = size(winidx,1);
    % Seed rng for repeatability
    rng(0);
    shifti = randi(nwin, [Settings.nboot,1])-1;

    % Normalize - no constant term here!
    normkin = my_normalize(dokin);

    % Do the shift regressions
    Rsq = nan(nneu,Settings.nboot);
    Bs = nan(nneu, nkin+1, Settings.nboot);
    Bint = nan(nneu, nkin+1, 2, Settings.nboot);
    
    parfor bi = 1:Settings.nboot
        winshiftrates = circshift(winrates,shifti(bi));
        doshiftrates = winshiftrates(winidx,:);
        % re-normalized shifted rates
        winnormrates = my_normalize(doshiftrates);

        for ni = 1:nneu
            [bb, bbints, ~, ~, stats] = regress(winnormrates(:,ni), [ones(npt,1), normkin]);
            Bs(ni, :, bi) = bb;
            Bint(ni, :, :, bi) = bbints;
            Rsq(ni, bi) = stats(1);
        end
    end
    LagBoot.(seg) = make_save_struct(Rsq, Bs, Bint);
end
end