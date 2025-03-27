% get Stretch: mean firing rates (FRs) across reaches with time-stretching
function Stretch = get_Stretch(Settings, R)

gkernstr = ['SmoothG' num2str(round(Settings.Stretch.gkern*1000))];

msedge = (R.Time.reachstarts(1)-Settings.msratepad:0.001:R.Time.reachstops(end)+Settings.msratepad)';
mstime = msedge(1:end-1)+0.0005;

nneu = size(R.N.SpkID,1);
nreach = length(R.Time.reachstarts);
npt = Settings.Stretch.nbef+Settings.Global.nbet+Settings.Stretch.naft;

Stretch.rates = nan(npt, nreach, nneu);
Stretch.speeds = nan(npt, nreach);

for neui = 1:nneu
    msrates = Spikes2FIR_Arbitrary( R.N.SpkTimes{neui}, msedge)';

    for reachi = 1:nreach
        thisbef = fliplr(R.Time.reachstarts(reachi)-(Settings.Stretch.step:Settings.Stretch.step:(Settings.Stretch.nbef*Settings.Stretch.step)))+Settings.Stretch.step/2;
        thisaft = R.Time.reachstops(reachi)+(Settings.Stretch.step:Settings.Stretch.step:(Settings.Stretch.naft*Settings.Stretch.step))-Settings.Stretch.step/2;
        betstep = ((R.Time.reachstops(reachi)-R.Time.reachstarts(reachi))/Settings.Global.nbet);
        thisbet = (R.Time.reachstarts(reachi):betstep:(R.Time.reachstops(reachi)-betstep))+betstep/2;
        thist = [thisbef, thisbet, thisaft]';
        Stretch.rates(:, reachi, neui) = GaussSmooth_Arbitrary( mstime, msrates, thist, Settings.Stretch.gkern, 5);
        if neui == 1
            % Resample speeds to the thist time
            Stretch.speeds(:, reachi) = interp1(R.Time.ktime, R.Kin.(gkernstr).Speed, thist);
            % Reinstate the nans
            for ti = 1:length(thist)
                % If there's a nan on either side of the thist, place a nan
                if any(R.Time.ktime == thist(ti))
                    if isnan(R.Kin.(gkernstr).Speed(R.Time.ktime == thist(ti)))
                        Stretch.speeds(ti,reachi) = nan;
                    end
                else
                    leftside = R.Kin.(gkernstr).Speed(find(R.Time.ktime < thist(ti),1,'last'));
                    rightside = R.Kin.(gkernstr).Speed(find(R.Time.ktime > thist(ti),1,'first'));
                    if isempty(leftside) || isempty(rightside) || isnan(leftside) || isnan(rightside)
                        Stretch.speeds(ti,reachi) = nan;
                    end
                end
            end
        end
    end
end
end