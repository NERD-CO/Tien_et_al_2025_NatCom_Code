% get Realtime: mean firing rates (FRs) across reaches without time-stretching
function Realtime = get_Realtime(Settings, R)

gkernstr = ['SmoothG' num2str(1000*Settings.Realtime.gkern, '%0.2i')];

% Gotta extend the time dangit
msedge = (R.Time.reachstarts(1)-Settings.msratepad:0.001:R.Time.reachstops(end)+Settings.msratepad)';
mstime = msedge(1:end-1)+0.0005;

nneu = size(R.N.SpkID,1);
nreach = length(R.Time.reachstarts);

Realtime.reach.rates = nan(Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+Settings.Realtime.naft, nreach, nneu);
Realtime.reach.speeds = nan(Settings.Realtime.nbef+Settings.Realtime.nbet+Settings.Realtime.nbet+Settings.Realtime.naft, nreach);

% Do other events... visual target cue, go cue, peak speed
alignnames = {'vis', 'go', 'start', 'peak', 'stop'};

for ai = 1:length(alignnames)
    if strcmp(alignnames{ai}, 'vis')
        Realtime.(alignnames{ai}).rates = nan(Settings.Realtime.naround, sum(R.Reach.outreach), nneu);
        Realtime.(alignnames{ai}).speeds = nan(Settings.Realtime.naround, sum(R.Reach.outreach));
    else
        Realtime.(alignnames{ai}).rates = nan(Settings.Realtime.naround, nreach, nneu);
        Realtime.(alignnames{ai}).speeds = nan(Settings.Realtime.naround, nreach);
    end
end


for neui = 1:nneu
    msrates = Spikes2FIR_Arbitrary( R.N.SpkTimes{neui}, msedge)';
    outcount = 0;
    for reachi = 1:nreach

        % And normal
        thisbef = fliplr(R.Time.reachstarts(reachi)-(Settings.Realtime.step:Settings.Realtime.step:(Settings.Realtime.nbef*Settings.Realtime.step)))+Settings.Realtime.step/2;
        thisaft = R.Time.reachstops(reachi)+(Settings.Realtime.step:Settings.Realtime.step:(Settings.Realtime.naft*Settings.Realtime.step))-Settings.Realtime.step/2;
        thisbet1 = R.Time.reachstarts(reachi)+(Settings.Realtime.step:Settings.Realtime.step:(Settings.Realtime.nbet*Settings.Realtime.step))-Settings.Realtime.step/2;
        thisbet2 = fliplr(R.Time.reachstops(reachi)-(Settings.Realtime.step:Settings.Realtime.step:(Settings.Realtime.nbet*Settings.Realtime.step)))+Settings.Realtime.step/2;

        thist = [thisbef, thisbet1, thisbet2, thisaft]';

        Realtime.reach.rates(:, reachi, neui) = GaussSmooth_Arbitrary( mstime, msrates, thist, Settings.Realtime.gkern, 5);

        if neui == 1
            Realtime.reach.speeds(:, reachi) = spline(R.Time.ktime, R.Kin.(gkernstr).Speed, thist);
            % Reinstate the nans
            for ti = 1:length(thist)
                % If there's a nan on either side of the thist, place a nan
                if any(R.Time.ktime == thist(ti))
                    if isnan(R.Kin.(gkernstr).Speed(R.Time.ktime == thist(ti)))
                        Realtime.reach.speeds(ti,reachi) = nan;
                    end
                else
                    leftside = R.Kin.(gkernstr).Speed(find(R.Time.ktime < thist(ti),1,'last'));
                    rightside = R.Kin.(gkernstr).Speed(find(R.Time.ktime > thist(ti),1,'first'));
                    if isempty(leftside) || isempty(rightside) || isnan(leftside) || isnan(rightside)
                        Realtime.reach.speeds(ti,reachi) = nan;
                    end
                end
            end
        end

        % Now do the alt reach periods
        for ai = 1:length(alignnames)
            if strcmp(alignnames{ai}, 'vis')
                if R.Reach.outreach(reachi)
                    outcount = outcount+1;
                    reftime = R.Time.dircue(reachi);
                else
                    continue;
                end
                bef = fliplr(reftime-(Settings.Realtime.step:Settings.Realtime.step:((Settings.Realtime.naround/2)*Settings.Realtime.step)))+Settings.Realtime.step/2;
                aft = reftime+(Settings.Realtime.step:Settings.Realtime.step:((Settings.Realtime.naround/2)*Settings.Realtime.step))-Settings.Realtime.step/2;
                thist = [bef aft];
                Realtime.(alignnames{ai}).rates(:,outcount,neui) = GaussSmooth_Arbitrary( mstime, msrates, thist, Settings.Realtime.gkern, 5);
                if neui == 1
                    Realtime.(alignnames{ai}).speeds(:,outcount) = spline(R.Time.ktime, R.Kin.(gkernstr).Speed, thist);
                    % Reinstate the nans
                    for ti = 1:length(thist)
                        % If there's a nan on either side of the thist, place a nan
                        if any(R.Time.ktime == thist(ti))
                            if isnan(R.Kin.(gkernstr).Speed(R.Time.ktime == thist(ti)))
                                Realtime.(alignnames{ai}).speeds(ti,outcount) = nan;
                            end
                        else
                            leftside = R.Kin.(gkernstr).Speed(find(R.Time.ktime < thist(ti),1,'last'));
                            rightside = R.Kin.(gkernstr).Speed(find(R.Time.ktime > thist(ti),1,'first'));
                            if isempty(leftside) || isempty(rightside) || isnan(leftside) || isnan(rightside)
                                Realtime.(alignnames{ai}).speeds(ti,outcount) = nan;
                            end
                        end
                    end
                end
            else
                if strcmp(alignnames{ai}, 'peak')
                    reftime = R.Time.reachpeaks(reachi);
                elseif strcmp(alignnames{ai}, 'start')
                    reftime = R.Time.reachstarts(reachi);
                elseif strcmp(alignnames{ai}, 'stop')
                    reftime = R.Time.reachstops(reachi);
                elseif strcmp(alignnames{ai}, 'go')
                    reftime = R.Time.gocue(reachi);
                end
                bef = fliplr(reftime-(Settings.Realtime.step:Settings.Realtime.step:((Settings.Realtime.naround/2)*Settings.Realtime.step)))+Settings.Realtime.step/2;
                aft = reftime+(Settings.Realtime.step:Settings.Realtime.step:((Settings.Realtime.naround/2)*Settings.Realtime.step))-Settings.Realtime.step/2;
                thist = [bef aft];
                Realtime.(alignnames{ai}).rates(:,reachi,neui) = GaussSmooth_Arbitrary( mstime, msrates, thist, Settings.Realtime.gkern, 5);
                if neui == 1
                    Realtime.(alignnames{ai}).speeds(:,reachi) = spline(R.Time.ktime, R.Kin.(gkernstr).Speed, thist);
                    % Reinstate the nans
                    for ti = 1:length(thist)
                        % If there's a nan on either side of the thist, place a nan
                        if any(R.Time.ktime == thist(ti))
                            if isnan(R.Kin.(gkernstr).Speed(R.Time.ktime == thist(ti)))
                                Realtime.(alignnames{ai}).speeds(ti,reachi) = nan;
                            end
                        else
                            leftside = R.Kin.(gkernstr).Speed(find(R.Time.ktime < thist(ti),1,'last'));
                            rightside = R.Kin.(gkernstr).Speed(find(R.Time.ktime > thist(ti),1,'first'));
                            if isempty(leftside) || isempty(rightside) || isnan(leftside) || isnan(rightside)
                                Realtime.(alignnames{ai}).speeds(ti,reachi) = nan;
                            end
                        end
                    end
                end
            end
        end
    end
end
end