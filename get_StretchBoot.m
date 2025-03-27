% get Stretch Boot: random mean firing rates to establish baseline for
% significance testing
function StretchBoot = get_StretchBoot(Settings, R)

startendpad = Settings.Stretch.naft*Settings.Stretch.step;

msedge = (R.Time.reachstarts(1)-Settings.msratepad:0.001:R.Time.reachstops(end)+Settings.msratepad)';
mstime = msedge(1:end-1)+0.0005;

nneu = size(R.N.SpkID,1);
nreach = length(R.Time.reachstarts);

rates = nan(nreach, nneu, Settings.nboot);

for neui = 1:nneu
    msrates = Spikes2FIR_Arbitrary( R.N.SpkTimes{neui}, msedge)';

    firsttime = R.Time.reachstarts(1)-startendpad;
    lasttime = R.Time.reachstops(end)+startendpad;
    boottimes = firsttime + (lasttime-firsttime)*rand(nreach,Settings.nboot);
    parfor booti = 1:Settings.nboot
        for reachi = 1:nreach
            rates(reachi, neui, booti) = GaussSmooth_Arbitrary( mstime, msrates, boottimes(reachi,booti), Settings.Stretch.gkern, 5);
        end
    end
end

StretchBoot = make_save_struct(rates);
end