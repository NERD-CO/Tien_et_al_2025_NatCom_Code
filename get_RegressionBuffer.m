% get Regression Buffer: kinematic and FR data for regressions
function B = get_RegressionBuffer(Settings, R)

gkernstr = ['SmoothG' num2str(1000*Settings.Regression.gkern)];

B.rates = R.FR.(gkernstr);
B.time = R.Time.ktime;

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

B.Idx.Reach = [];
B.Reachnum.Reach = [];
for ri = 1:nreach
    usekt = find(B.time > R.Time.reachstarts(ri)-Settings.Regression.reachpad & B.time < R.Time.reachstops(ri)+Settings.Regression.reachpad);
    B.Idx.Reach = [B.Idx.Reach; usekt];
    B.Reachnum.Reach =  [B.Reachnum.Reach; ri*ones(size(usekt))];
end
B.Idx.ReachE = B.Idx.Reach;
B.Reachnum.ReachE = B.Reachnum.Reach;

% Get Hold indices - only hold periods between two valid reaches. Hold
% reachnum reflects the number of the reach leading up to the hold.
B.Idx.Hold = [];
B.Reachnum.Hold = [];
for ri = 1:nreach-1
    % Only using holds between two good reaches
    if (R.Reach.outreach(ri) && R.Reach.reachnum(ri) == R.Reach.reachnum(ri+1) && ~R.Reach.outreach(ri+1)) || (~R.Reach.outreach(ri) && R.Reach.reachnum(ri) == R.Reach.reachnum(ri+1)-1 && R.Reach.outreach(ri+1))
        usekt = find(B.time >= R.Time.reachstops(ri)+Settings.Regression.reachpad & B.time <= R.Time.reachstarts(ri+1)-Settings.Regression.reachpad);
        B.Idx.Hold = [B.Idx.Hold; usekt];
        B.Reachnum.Hold = [B.Reachnum.Hold; ri*ones(size(usekt))];
    end
end

% Full is just good reaches and good holds combined
B.Idx.Full = sort([B.Idx.Reach; B.Idx.Hold]);

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
B.errxy = nan(length(B.Idx.Reach),ndim);
for rii = 1:length(B.Idx.Reach)
    B.errxy(rii,:) = endpos(B.Reachnum.Reach(rii),:) - B.kin(B.Idx.Reach(rii),1:3);
end
B.emag = vecnorm(B.errxy')';
B.kinstrerr = [B.kinstr {'X Err', 'Y Err', 'Z Err', 'EMag'}];
end