function plot_3DKin(Settings, Data, RegressionBuffer, dati)

gkernstr = ['SmoothG' num2str(1000*Settings.Plot.gkern,'%0.2i')];

pltcols = jet(8)*0.9;

R = Data(dati);
B = RegressionBuffer(dati);

rdirs = R.Reach.reachdirs;
rdirs(R.Reach.outreach == 0) = rdirs(R.Reach.outreach == 0) + 180; % Flip the INREACHES silly!
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360; % For good measure
rdirs = round(rdirs/45)+1;
% Fix for the new targdirs reflection
rdirs = target_idx_swap(rdirs);

pos = R.Kin.(gkernstr).Pos;

posbuff = pos(B.Idx.Full,:);
npos = size(posbuff,1);
fullreachi = find(ismember(B.Idx.Full, B.Idx.Reach));
fullholdi = find(ismember(B.Idx.Full, B.Idx.Hold));

allreachtargi = nan(npos,1);
allreachtargi(fullreachi) = rdirs(B.Reachnum.Reach);
holdi = rdirs(B.Reachnum.Hold);
holdi(R.Reach.outreach(B.Reachnum.Hold) == 0) = 0;
allreachtargi(fullholdi) = holdi;

figure('Renderer', 'painters');
clf;
hold on;

% iterate through posbuff and plot contiguous segments
darti = diff(allreachtargi);
segstarts = find(darti~=0)+1;
for segi = 1:(length(segstarts)-1)
    thistarg = allreachtargi(segstarts(segi));
    if thistarg==0
        thiscol = [0 0 0];
    else
        thiscol = pltcols(thistarg,:);
    end
    thisi = segstarts(segi):segstarts(segi+1)-1;
    plot3(posbuff(thisi,1), posbuff(thisi,2), posbuff(thisi,3), 'Color', thiscol, 'LineWidth', 1.5);
end

axis equal
xlabel('X');
ylabel('Y');
zlabel('Z');

% Figure out top 2 PCA view
% Do PCA to find the plane of greatest variation (top two PCs)
c = pca(posbuff);
view(c(:,3));
camroll(170)
tickers = -250:100:250;
xticks(tickers)
yticks(tickers)
zticks(tickers)
xticklabels([]);
yticklabels([]);
zticklabels([]);
grid on;
box on;
