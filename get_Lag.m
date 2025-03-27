% get Lag: lagged regression results
function Lag = get_Lag(Settings, B)

nseg = length(Settings.Regression.segstr);

Lag.lagsamps = -Settings.Regression.maxlagsamps:Settings.Regression.maxlagsamps;
Lag.lagtime = Lag.lagsamps/Settings.vidframerate;
nlag = length(Lag.lagsamps);

for segi = 1:nseg
    thisseg = Settings.Regression.segstr{segi};

    if strcmp(thisseg,'ReachE')
        dokin = [B.kin(B.Idx.(thisseg),:), B.errxy, B.emag];
    else
        dokin = B.kin(B.Idx.(thisseg),:);
    end

    npt = size(dokin,1);
    nkin = size(dokin,2);
    nneu = size(B.rates,2);
    nwhole = size(B.time,1);

    % Pad rates if needed to accept max lags
    npadfront = max([Settings.Regression.maxlagsamps-B.Idx.(thisseg)(1)+1, 0]);
    npadend = max([B.Idx.(thisseg)(end)+Settings.Regression.maxlagsamps-nwhole+1, 0]);
    padrates = [nan(npadfront,nneu); B.rates; nan(npadend,nneu)];
    padidx = B.Idx.(thisseg) + npadfront; % Adjust index so we can still pick out the right ones

    % Normalize - no constant term here!
    normkin = my_normalize(dokin);

    % Do the lag regressions
    Rsq = nan(nneu,nlag);
    Bs = nan(nneu, nkin+1, nlag);
    Bint = nan(nneu, nkin+1, 2, nlag);

    for li = 1:nlag
        lagrates = padrates(padidx+Lag.lagsamps(li),:);

        % Re-normalize lagged rates
        lagnormrates = my_normalize(lagrates);

        for ni = 1:nneu
            [bb, bbints, ~, ~, stats] = regress(lagnormrates(:,ni), [ones(npt,1), normkin]);
            Bs(ni, :, li) = bb;
            Bint(ni, :, :, li) = bbints;
            Rsq(ni, li) = stats(1);
        end
    end
    Lag.(thisseg) = make_save_struct(Rsq, Bs, Bint);
end
end