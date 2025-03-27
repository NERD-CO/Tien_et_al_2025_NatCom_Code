% get Realtime Sig: significance testing for Realtime FRs
function RealtimeSig = get_RealtimeSig(Settings, Realtime, StretchBoot)

for alphi = 1:length(Settings.alph2do)
    alphstr = ['alpha' num2str(100*Settings.alph2do(alphi),'%0.2i')];

    % Now do all the alt periods
    rtfn = fieldnames(Realtime);
    rtfn = rtfn(~strcmp(rtfn,'rates') & ~strcmp(rtfn,'speeds'));
    for ai = 1:length(rtfn)
        aistr = rtfn{ai};

        nneu = size(Realtime.(aistr).rates,3);
        npt = size(Realtime.(aistr).rates,1);
        nreach = size(Realtime.(aistr).rates,2);
        nboot = size(StretchBoot.rates,3);

        ps = nan(nneu, npt);
        hilo = zeros(nneu, npt);
        mrates = nan(nneu, npt);
        sigi = nan(nneu,1);
        sigofem = nan(nneu,npt);
        sighilo = nan(nneu,npt);

        for neui = 1:nneu
            mrates(neui,:) = squeeze(mean(Realtime.(aistr).rates(:,:,neui),2));
            mbrates = squeeze(mean(StretchBoot.rates(1:nreach,neui,:),1)); % Only take the nreach so that vis is matched

            ps(neui,:) = min([sum(mbrates > mrates(neui,:))/nboot; sum(mbrates < mrates(neui,:))/nboot],[],1);

            hilo(neui,:) = sign(mrates(neui,:) - mean(mbrates));

            % Now find inarow, only count if inarow is > cutoff. Do front
            % and back halves separately for reach
            if strcmp(aistr,'reach')
                halfnpt = npt/2;
                halfsamp{1} = 1:halfnpt;
                halfsamp{2} = (halfnpt+1):npt;

                for halfi = 1:2
                    allofem = ps(neui,halfsamp{halfi}) < Settings.alph2do(alphi)/2;

                    % Trim < cutoff chunks
                    sigofem(neui,halfsamp{halfi}) = allofem;
                    incount = 0;
                    firsti = 1;
                    thislen = 0;
                    for sampi = 1:halfnpt
                        if (sigofem(neui,halfsamp{halfi}(sampi))==1 && incount==0)
                            firsti = sampi;
                            incount = 1;
                            thislen = 1;
                        elseif (sigofem(neui,halfsamp{halfi}(sampi))==1 && sampi < halfnpt)
                            thislen = thislen + 1;
                        elseif (sigofem(neui,halfsamp{halfi}(sampi))==1)
                            thislen = thislen + 1;
                            if (thislen < Settings.cutoff)
                                sigofem(neui,halfsamp{halfi}(firsti:sampi)) = 0;
                            end
                        else
                            incount = 0;
                            if (thislen < Settings.cutoff)
                                sigofem(neui,halfsamp{halfi}(firsti:sampi-1)) = 0;
                                thislen = 0;
                            end
                        end
                    end
    
                    if sigofem(neui,halfsamp{halfi}(end)) == 1 && sigofem(neui,halfsamp{halfi}(end-1)) == 0
                        sigofem(neui,halfsamp{halfi}(end)) = 0;
                    end
    
                    sighilo(neui,halfsamp{halfi}) = hilo(neui,halfsamp{halfi}).*sigofem(neui,halfsamp{halfi});
                end

                sigi = any(sigofem,2);
            else
                allofem = ps(neui,:) < Settings.alph2do(alphi)/2;
                bads = find(allofem == 0);
                dbad = diff(bads);
                inarow = max([bads(1)-1, max(dbad)-1, npt-bads(end)]);
                sigi(neui) = inarow >= Settings.cutoff;

                % Trim < cutoff chunks
                sigofem(neui,:) = allofem;
                incount = 0;
                firsti = 1;
                thislen = 0;
                for sampi = 1:npt
                    if (sigofem(neui,sampi)==1 && incount==0)
                        firsti = sampi;
                        incount = 1;
                        thislen = 1;
                    elseif (sigofem(neui,sampi)==1 && sampi < npt)
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
        end

        RealtimeSig.(aistr).(alphstr) = make_save_struct(ps, hilo, mrates, sigi, sigofem, sighilo);
    end

end
end