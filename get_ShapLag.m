% get ShapLag: Shapley decomposition of regressions at all lags
function Shap = get_ShapLag(Settings, B, seg)    
    lagsamps = -Settings.Regression.maxlagsamps:Settings.Regression.maxlagsamps;
    nlag = length(lagsamps);
    
    parts = Settings.Regression.parts_err;
    dokin = [B.kin(B.Idx.(seg),:), B.errxy, B.emag];
    kinstr = B.kinstrerr;
    
    nkin = size(dokin,2);
    nneu = size(B.rates,2);
    nwhole = size(B.time,1);
    
    % Pad rates if needed to accept max lags
    npadfront = max([Settings.Regression.maxlagsamps-B.Idx.(seg)(1)+1, 0]);
    npadend = max([B.Idx.(seg)(end)+Settings.Regression.maxlagsamps-nwhole+1, 0]);
    padrates = [nan(npadfront,nneu); B.rates; nan(npadend,nneu)];
    padidx = B.Idx.(seg) + npadfront; % Adjust index so we can still pick out the right ones
    
    % Normalize - no constant term here!
    normkin = my_normalize(dokin);
    
    lagmap = true(nlag,1);
    shap = nan(nneu, nkin, nlag);
    parfor li = 1:nlag
        lagrates = padrates(padidx+lagsamps(li),:);
        lagnormrates = my_normalize(lagrates);
        for ni = 1:nneu
            shap(ni, :, li) = do_shapley_value_decomposition(lagnormrates(:,ni), normkin);
        end
    end
    
    Shap = make_save_struct(shap, seg, parts, kinstr, lagmap);
end

% Actually do it
function shap = do_shapley_value_decomposition(normrates, normkin)
    
    ndat = size(normrates,1);
    
    % Set up Shapley combination numbers
    nkin = size(normkin,2);
    shapk = nan(nkin,1);
    for kini = 1:nkin
        combs{kini} = nchoosek(1:nkin, kini);
        shapk(kini) = sum(any(combs{kini}==1,2));
    end
    
    addone = ones(ndat,1);
    
    % Now do Shapley!
    shap = zeros(nkin,1);
    for kini = 1:nkin
        ncomb = size(combs{kini},1);
        thisRsq = nan(ncomb,1);
        parfor combi = 1:ncomb
            [~,~,~,~,stats] = regress(normrates, [addone, normkin(:,combs{kini}(combi,:))]);
            thisRsq(combi) = stats(1);
        end
        for combi = 1:ncomb
            for kki = 1:nkin
                if any(combs{kini}(combi,:)==kki)
                    shap(kki) = shap(kki) + ((1/shapk(kini))*thisRsq(combi));
                else
                    shap(kki) = shap(kki) - ((1/shapk(kini+1))*thisRsq(combi));
                end
            end
        end
    end
    shap = shap/nkin;
end