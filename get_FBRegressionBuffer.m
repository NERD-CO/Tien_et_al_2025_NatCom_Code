% get FB Regression Buffer: data buffer for front (reach start) vs back
% (reach end) windowed regressions
function B = get_FBRegressionBuffer(Settings, R)

FBWin = 0.5;

gkernstr = ['SmoothG' num2str(1000*Settings.Regression.gkern)];

B.rates = R.FR.(gkernstr);
B.time = R.Time.ktime;

nneu = size(B.rates,2);
nratept = size(B.rates,1);
nk = length(R.Time.ktime);
ndim = size(R.Kin.(gkernstr).Pos,2);
nreach = length(R.Time.reachstarts);

% Create B.kin and B.kinstr, just put in all the B.kin fields
B.kin = [];
dims = {'X', 'Y', 'Z'};
kinfn = fieldnames(R.Kin.(gkernstr));
kinstrticker = 1;
for fni = 1:length(kinfn)
    thisndim = size(R.Kin.(gkernstr).(kinfn{fni}), 2);
    if thisndim == 1
        B.kin = [B.kin, R.Kin.(gkernstr).(kinfn{fni})];
        B.kinstr{kinstrticker} = kinfn{fni};
        kinstrticker = kinstrticker+1;
    else
        for dimi = 1:3
            B.kin = [B.kin, R.Kin.(gkernstr).(kinfn{fni})(:,dimi)];
            B.kinstr{kinstrticker} = [dims{dimi} ' ' kinfn{fni}];
            kinstrticker = kinstrticker+1;
        end
    end
end

B.Idx.FrontE = [];
B.Reachnum.FrontE = [];
for ri = 1:nreach
    usekt = find(B.time > R.Time.reachstarts(ri)-(FBWin/2) & B.time < R.Time.reachstarts(ri)+(FBWin/2));
    B.Idx.FrontE = [B.Idx.FrontE; usekt];
    B.Reachnum.FrontE =  [B.Reachnum.FrontE; ri*ones(size(usekt))];
end

B.Idx.BackE = [];
B.Reachnum.BackE = [];
for ri = 1:nreach
    usekt = find(B.time > R.Time.reachstops(ri)-(FBWin/2) & B.time < R.Time.reachstops(ri)+(FBWin/2));
    B.Idx.BackE = [B.Idx.BackE; usekt];
    B.Reachnum.BackE =  [B.Reachnum.BackE; ri*ones(size(usekt))];
end

% End positions for target error calculation. Take the average
% position from reachpad to errorwindow after reach stop
erri = cell(nreach,1);
for ri = 1:nreach
    erri{ri} = find(B.time > R.Time.reachstops(ri)+Settings.Regression.reachpad & B.time < R.Time.reachstops(ri)+Settings.Regression.errorwindow);
end

posi = find(contains(B.kinstr,'Pos'));
endpos = nan(nreach,ndim);
for ri = 1:nreach
    endpos(ri,:) = nanmean(B.kin(erri{ri},posi),1);
end

% Find x-y error term
B.Ferrxy = nan(length(B.Idx.FrontE),ndim);
for rii = 1:length(B.Idx.FrontE)
    B.Ferrxy(rii,:) = endpos(B.Reachnum.FrontE(rii),:) - B.kin(B.Idx.FrontE(rii),1:3);
end
B.Femag = vecnorm(B.Ferrxy')';

% Find x-y error term
B.Berrxy = nan(length(B.Idx.BackE),ndim);
for rii = 1:length(B.Idx.BackE)
    B.Berrxy(rii,:) = endpos(B.Reachnum.BackE(rii),:) - B.kin(B.Idx.BackE(rii),1:3);
end
B.Bemag = vecnorm(B.Berrxy')';

B.kinstrerr = [B.kinstr {'X Err', 'Y Err', 'Z Err', 'EMag'}];
end