function plot_Cutoffs_Reach(Settings, StretchSig)

cutoffs = 1:100;
nalpha = length(Settings.alph2do);
ncutoff = length(cutoffs);
nall = Settings.Global.allnneu;
npt = Settings.Stretch.nbef + Settings.Global.nbet + Settings.Stretch.naft;

cutsigs = nan(ncutoff, nalpha);

allps = nan(nall, npt, nalpha);

ticker = 1;
for dati = 1:length(StretchSig)
    
    for alphi = 1:nalpha
        alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];

        nneu = size(StretchSig(dati).(alphstr).ps,1);

        allps(ticker:ticker+nneu-1, :, alphi) = StretchSig(dati).(alphstr).ps;
    end

    ticker = ticker+nneu;
end

for alphi = 1:nalpha
    for cuti = 1:ncutoff
        cutoff = cutoffs(cuti);

        allofem = allps(:,:,alphi) < Settings.alph2do(alphi)/2;
        inarow = nan(nall,1);
        sigofem = nan(nall,npt);
        for ni = 1:nall
            bads = find(allofem(ni,:)==0);
            dbad = diff(bads);
            inarow(ni) = max([bads(1)-1, max(dbad)-1, npt-bads(end)]);

            % Trim < cutoff chunks
            sigofem(ni,:) = allofem(ni,:);
            incount = 0;
            firsti = 1;
            thislen = 0;
            for sampi = 1:npt
                if (sigofem(ni,sampi)==1 && incount==0)
                    firsti = sampi;
                    incount = 1;
                    thislen = 1;
                elseif (sigofem(ni,sampi)==1 && sampi < npt)
                    thislen = thislen + 1;
                elseif (sigofem(ni,sampi)==1)
                    thislen = thislen + 1;
                    if (thislen < cutoff)
                        sigofem(ni,firsti:sampi) = 0;
                    end
                else
                    incount = 0;
                    if (thislen < cutoff)
                        sigofem(ni,firsti:sampi-1) = 0;
                        thislen = 0;
                    end
                end
            end
        
            if sigofem(ni,end) == 1 && sigofem(ni,end-1) == 0
                sigofem(ni,end) = 0;
            end
        end
%         sigi = inarow>=cutoff;

        sigi = any(sigofem(:,(Settings.Stretch.nbef-Settings.Stretch.periwin+1):(Settings.Stretch.nbef+Settings.Global.nbet+Settings.Stretch.periwin+1)),2);      

        cutsigs(cuti,alphi) = sum(sigi);
    end
end

figure;
hold on;
lastcut = 79;
plot(cutoffs(1:lastcut+1)*Settings.Stretch.step*1000, cutsigs(1:lastcut+1,1)/length(sigi), '-k', 'LineWidth', 2);
plot(cutoffs(1:lastcut+1)*Settings.Stretch.step*1000, cutsigs(1:lastcut+1,2)/length(sigi), ':k', 'LineWidth', 2);
plot([5*Settings.Stretch.step*1000 5*Settings.Stretch.step*1000],[0 cutsigs(cutoffs==5,2)/length(sigi)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
plot([0 5*Settings.Stretch.step*1000],[cutsigs(cutoffs==5,2)/length(sigi) cutsigs(cutoffs==5,2)/length(sigi)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
plot(5*Settings.Stretch.step*1000,cutsigs(cutoffs==5,2)/length(sigi), '.k', 'MarkerSize', 20);

yticks([0 cutsigs(cutoffs==5,2)/length(sigi) 0.5 1]);
ylim([0 1]);
xticks([0 50 200 400 600 800]);
set(gca,'TickDir','out');
set(gca,'TickLength',[0.03 0.01]);
set(gca,'box', 'off')
legend({'\alpha = 0.05', '\alpha = 0.01'})