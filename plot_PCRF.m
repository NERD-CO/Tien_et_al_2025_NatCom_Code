function plot_PCRF(Settings, Data, RegressionBuffer, doi, neui, gridsq)

fg = figure;

gkernstr = ['SmoothG' num2str(1000*Settings.Plot.gkern,'%0.2i')];

R = Data(doi);
B = RegressionBuffer(doi);

usepos = R.Kin.(gkernstr).Pos;
userate = R.FR.(gkernstr);

fullholdi = find(ismember(B.Idx.Full, B.Idx.Hold));

posbuff = usepos(B.Idx.Full,:);
ratebuff = userate(B.Idx.Full,neui);

% Figure out top 2 PCA view
% Do PCA to find the plane of greatest variation (top two PCs)
[c,s,l,t,e] = pca(posbuff);

sco = s(:,1:2);

holdsco = sco(fullholdi,:);

% Rotate so that target 1 is left
rdirs = R.Reach.reachdirs;
rdirs(R.Reach.outreach == 0) = rdirs(R.Reach.outreach == 0) + 180;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360; % For good measure
rdirs = round(rdirs/45)+1;
% Fix for the new targdirs reflection
rdirs = target_idx_swap(rdirs);
udirs = unique(rdirs);
ndirs = length(udirs);
mhold = nan(ndirs,2);
for diri = 1:ndirs
    thisi = find(rdirs == udirs(diri) & R.Reach.outreach==1);
    holdi = ismember(B.Reachnum.Hold, thisi);
    if any(holdi)
        mhold(diri,:) = nanmean(holdsco(holdi,:),1);
    end
end
% Also find center mean hold pos
ctri = find(R.Reach.outreach==0);
holdi = ismember(B.Reachnum.Hold,ctri);
chold = nanmean(holdsco(holdi,:),1);

sco = sco - repmat(chold,[size(sco,1),1]); % subtract center from sco

% Elim those with no holds
udirs = udirs(~isnan(mhold(:,1)));
mhold = mhold(~isnan(mhold(:,1)),:);

% Subtract center and norm
for diri = 1:length(udirs)
    mhold(diri,:) = mhold(diri,:)-chold;
    mhold(diri,:) = mhold(diri,:)/norm(mhold(diri,:));
end

% Determine if we are going clockwise
pickdir = find(diff(udirs) < 4,1,'first');
ang1 = atan2d(mhold(pickdir,2), mhold(pickdir,1));
if ang1 < 0
    ang1 = ang1+360;
end
ang2 = atan2d(mhold(pickdir+1,2), mhold(pickdir+1,1));
if ang2 < 0
    ang2 = ang2+360;
end
adiff = ang2-ang1;
if adiff < 0
    adiff = adiff+360;
end
% Flip y if counterclockwise, otherwise do nothing
if sign(adiff) == 1
    sco(:,2) = -sco(:,2);
    mhold(:,2) = -mhold(:,2);
end

% Ideal target angles
ideal = [180 135 90 45 0 315 270 225]';
ideal = ideal(udirs);

% find average that we are off from ideal target angles
ahold = nan(size(mhold,1),1);
for diri = 1:length(udirs)
    ahold(diri) = atan2d(mhold(diri,2), mhold(diri,1));
    if ahold(diri) < 0
        ahold(diri) = ahold(diri)+360;
    end
end
doff = ideal-ahold;
doff(doff<0) = doff(doff<0) + 360;
if range(doff) > 180
    doff(doff<180) = doff(doff<180)+360;
end
moff = mean(doff);
rotmat = [cosd(moff) -sind(moff); sind(moff) cosd(moff)];
sco = (rotmat*sco')';
mhold = (rotmat*mhold')';

% Create mesh
gridlims = nan(2,2);
for i = 1:2
    gridlims(i,1) = floor(min(sco(:,i))/gridsq);
    gridlims(i,2) = ceil(max(sco(:,i))/gridsq);
end

grid1 = (gridlims(1,1)*gridsq):gridsq:(gridlims(1,2)*gridsq);
grid2 = (gridlims(2,1)*gridsq):gridsq:(gridlims(2,2)*gridsq);

ng1 = length(grid1)-1;
ng2 = length(grid2)-1;

frgrid = nan(ng1, ng2);
countgrid = nan(ng1, ng2);

for gi1 = 1:ng1
    for gi2 = 1:ng2
        insq = sco(:,1) > grid1(gi1) & sco(:,1) <= grid1(gi1+1) & sco(:,2) > grid2(gi2) & sco(:,2) <= grid2(gi2+1);
        countgrid(gi1,gi2) = sum(insq);
        frgrid(gi1,gi2) = mean(ratebuff(insq));
    end
end

set(0, 'CurrentFigure', fg);
clf;
hold on;

% Make custom hot colormap to yellow
topcut = 200;
cbits = 512;
hotty = hot;
hotty = hotty(1:topcut,:);
hotfull = nan(cbits,3);
hotfull(:,1) = interp1(1:topcut, hotty(:,1), linspace(1,topcut,cbits));
hotfull(:,2) = interp1(1:topcut, hotty(:,2), linspace(1,topcut,cbits));
hotfull(:,3) = interp1(1:topcut, hotty(:,3), linspace(1,topcut,cbits));

% Normalize to 1
frgrid = (frgrid - min(frgrid(isfinite(frgrid))))/range(frgrid(isfinite(frgrid)));
frgrid(~isfinite(frgrid)) = 0;
for gi1 = 1:ng1
    for gi2 = 1:ng2
        if countgrid(gi1,gi2) == 0
            alpha = 0;
        else
            alpha = 1;
        end

        patch([grid1(gi1) grid1(gi1) grid1(gi1+1) grid1(gi1+1)], [grid2(gi2) grid2(gi2+1) grid2(gi2+1) grid2(gi2)], hotfull(round(frgrid(gi1,gi2)*(cbits-1))+1,:), 'LineStyle', 'none', 'FaceAlpha', alpha);
    end
end

xlim([grid1(1)-gridsq grid1(end)+gridsq]);
ylim([grid2(1)-gridsq grid2(end)+gridsq]);

axis square

yl = ylim;
xl = xlim;
yr = diff(yl);
xr = diff(xl);
if yr > xr
    xlim([xl(1)-(yr-xr)/2 xl(2)+(yr-xr)/2]);
    ylim(yl);
elseif xr > yr
    xlim(xl);
    ylim([yl(1)-(xr-yr)/2 yl(2)+(xr-yr)/2]);
end

% Gray background
set(gca, 'Color', [0.7 0.7 0.7])

xticks([]);
yticks([]);
% Do scale bar instead of ticks
yl = ylim;
xl = xlim;
bwx = 0.015*diff(xl);
bwy = 0.015*diff(yl);
bl = 50;
patch([xl(1) xl(1) xl(1)+bwx xl(1)+bwx xl(1)+bl xl(1)+bl], [yl(2)-bl yl(2) yl(2) yl(2)-(bl-bwy) yl(2)-(bl-bwy) yl(2)-bl], 'k', 'LineStyle', 'none');
ylim(yl);
xlim(xl);
