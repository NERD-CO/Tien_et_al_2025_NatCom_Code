% get Stretch Sig: significance testing for Stretch FRs
function StretchSig = get_StretchSig(Settings, Stretch, StretchBoot)

nneu = size(Stretch.rates,3);
nsamp = size(Stretch.rates,1);
nboot = size(StretchBoot.rates,3);

for alphi = 1:length(Settings.alph2do)
    alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];
    ps = nan(nneu, nsamp);
    hilo = zeros(nneu, nsamp);
    mrates = nan(nneu, nsamp);
    sigi = nan(nneu,1);
    sigofem = nan(nneu,nsamp);
    sighilo = nan(nneu,nsamp);
    
    for neui = 1:nneu
        mrates(neui,:) = squeeze(mean(Stretch.rates(:,:,neui),2));
        mbrates = squeeze(mean(StretchBoot.rates(:,neui,:),1));
    
        ps(neui,:) = min([sum(mbrates > mrates(neui,:))/nboot; sum(mbrates < mrates(neui,:))/nboot],[],1);
    
        hilo(neui,:) = sign(mrates(neui,:) - mean(mbrates));
    
        % Now find inarow, only count if inarow is > cutoff
        allofem = ps(neui,:) < Settings.alph2do(alphi)/2;
        bads = find(allofem == 0);
        dbad = diff(bads);
        inarow = max([bads(1)-1, max(dbad)-1, nsamp-bads(end)]);
        sigi(neui) = inarow >= Settings.cutoff;
    
        % Trim < cutoff chunks
        sigofem(neui,:) = allofem;
        incount = 0;
        firsti = 1;
        thislen = 0;
        for sampi = 1:nsamp
            if (sigofem(neui,sampi)==1 && incount==0)
                firsti = sampi;
                incount = 1;
                thislen = 1;
            elseif (sigofem(neui,sampi)==1 && sampi < nsamp)
                thislen = thislen + 1;
            elseif (sigofem(neui,sampi)==1)
                thislen = thislen + 1;
                if (thislen < Settings.cutoff)
                    sigofem(neui,firsti:sampi) = 0;
                end
            else
                incount = 0;
                if (thislen < Settings.cutoff)
                    sigofem(neui,firsti:sampi-1) = 0;
                    thislen = 0;
                end
            end
        end
    
        if sigofem(neui,end) == 1 && sigofem(neui,end-1) == 0
            sigofem(neui,end) = 0;
        end
    
        sighilo(neui,:) = hilo(neui,:).*sigofem(neui,:);
    end
    
    StretchSig.(alphstr) = make_save_struct(ps, hilo, mrates, sigi, sigofem, sighilo);
end
end