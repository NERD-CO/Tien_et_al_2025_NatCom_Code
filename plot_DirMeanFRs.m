function plot_DirMeanFRs(Settings, Data, Stretch, doi, neui)

R = Data(doi);

% Get the reaches by targ, not dir
rdirs = R.Reach.reachdirs;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360; % For good measure
rdirs = round(rdirs/45)+1;

% Fix for the new targdirs reflection
rdirs = target_idx_swap(rdirs);
udirs = unique(rdirs);
ndirs = length(udirs);

frontcut = Settings.Stretch.nbef-50;
backcut = Settings.Stretch.naft-50;
Stretch(doi).rates = Stretch(doi).rates(frontcut+1:end-backcut,:,:);
npt = size(Stretch(doi).rates,1);
cutSettings.Stretch.nbef = Settings.Stretch.nbef-frontcut;
cutSettings.Stretch.naft = Settings.Stretch.naft-backcut;

pltcols = 0.9*jet(8);

figure('Renderer', 'painters');
clf;
hold on;
mfrs = nan(npt,ndirs);
sfrs = nan(npt,ndirs);
for diri = 1:ndirs
    mfrs(:,diri) = nanmean(Stretch(doi).rates(:,rdirs==(udirs(diri)),neui),2);
    sfrs(:,diri) = nanstd(Stretch(doi).rates(:,rdirs==(udirs(diri)),neui),0,2);
    plot(mfrs(:,diri), 'Color', pltcols(udirs(diri),:), 'LineWidth', 3);
    tops = mfrs(:,diri)+sfrs(:,diri);
    bots = mfrs(:,diri)-sfrs(:,diri);
    bots(bots<0) = 0;
    patch([1:npt, fliplr(1:npt)], [tops; flipud(bots)], pltcols(udirs(diri),:), 'FaceAlpha', 0.1, 'LineStyle', 'none', 'HandleVisibility', 'off');
end
xlim([0, npt]);
yl = ylim;
plot([cutSettings.Stretch.nbef+0.5, cutSettings.Stretch.nbef+0.5], yl, '--k', 'HandleVisibility', 'off');
plot([cutSettings.Stretch.nbef+Settings.Global.nbet+0.5, cutSettings.Stretch.nbef+Settings.Global.nbet+0.5], yl, '--k', 'HandleVisibility', 'off');
xticks([cutSettings.Stretch.nbef+0.5, cutSettings.Stretch.nbef+Settings.Global.nbet+0.5]);
yticks([0 30 60]);
%                 yticklabels([]);
xticklabels([]);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength',[0.035 0.01]);

plot([npt-50 npt], [50 50], '-k', 'LineWidth', 4)