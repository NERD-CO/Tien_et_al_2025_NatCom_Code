% get Directionality: significance testing for directional tuning of
% FRs
function DirSig = get_Directionality(Settings, R, Stretch)

rates = Stretch.rates;
npt = size(Stretch.rates,1);
nneu = size(rates,3);

% Get the reaches by dir, not targ
rdirs = R.Reach.reachdirs;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360; % For good measure
rdirs = round(rdirs/45)+1;
% Fix for the new targdirs reflection
rdirs = target_idx_swap(rdirs);

pkw = nan(nneu, npt);

for neui = 1:nneu
    for pti = 1:npt
        % Get KW anova p value
        pkw(neui, pti) = kruskalwallis(squeeze(rates(pti,:,neui)), rdirs, 'off');
    end
end

DirSig = make_save_struct(pkw);

for alphi = 1:length(Settings.alph2do)
    alphstr = ['alpha' num2str(100*Settings.alph2do(alphi), '%0.2i')];
    DirSig.(alphstr).sigs = pkw < Settings.alph2do(alphi);

    for neui = 1:nneu
        % Now find inarow, only count if inarow is > cutoff
        allofem = DirSig.(alphstr).sigs(neui,:);
        bads = find(allofem == 0);
        dbad = diff(bads);
        inarow = max([bads(1)-1, max(dbad)-1, npt-bads(end)]);

        % Trim < cutoff chunks
        incount = 0;
        firsti = 1;
        thislen = 0;
        for sampi = 1:npt
            if (allofem(sampi)==1 && incount==0)
                firsti = sampi;
                incount = 1;
                thislen = 1;
            elseif (allofem(sampi)==1 && sampi < npt)
                thislen = thislen + 1;
            elseif (allofem(sampi)==1)
                thislen = thislen + 1;
                if (thislen < Settings.cutoff)
                    allofem(firsti:sampi) = 0;
                end
            else
                incount = 0;
                if (thislen < Settings.cutoff)
                    allofem(firsti:sampi-1) = 0;
                    thislen = 0;
                end
            end
        end
    
        if allofem(end) == 1 && allofem(end-1) == 0
            allofem(end) = 0;
        end
    
        DirSig.(alphstr).sigs(neui,:) = allofem;
    end
end
end