function plot_Raster(Settings, Data, doi, neui)

R = Data(doi);

pltcols = 0.9*jet(8);

nreach = length(R.Time.reachstarts);

start_vs_peak = R.Time.reachstarts-R.Time.reachpeaks;
stop_vs_peak = R.Time.reachstops-R.Time.reachpeaks;

st_peak = cell(nreach,1);
for i = 1:nreach
    st_peak{i} = R.N.SpkTimes{neui}(R.N.SpkTimes{neui} > R.Time.reachpeaks(i)-Settings.Raster.winaround & R.N.SpkTimes{neui} < R.Time.reachpeaks(i)+Settings.Raster.winaround) - R.Time.reachpeaks(i);
end

rdirs = R.Reach.reachdirs;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs >= 360) = rdirs(rdirs >= 360) - 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360;
rdirs(rdirs < 0) = rdirs(rdirs < 0) + 360; % For good measure
rdirs = round(rdirs/45)+1;

% Fix for the new targdirs reflection
rdirs = target_idx_swap(rdirs);

udir = unique(rdirs);
numdir = udir;
ndir = length(udir);

plottick = 1;
figure('Renderer', 'painters'); hold on;

for diri = 1:ndir
    thisdir = find(rdirs == udir(diri));
    nthisdir = length(thisdir);
    sorti = 1:nthisdir;

    for diri2 = 1:nthisdir
        for ticki = 1:length(st_peak{thisdir(sorti(diri2))})
            % Rasters
            plot([st_peak{thisdir(sorti(diri2))}(ticki), st_peak{thisdir(sorti(diri2))}(ticki)], [plottick-0.5, plottick+0.5], '-', 'Color', pltcols(numdir(diri),:), 'LineWidth', 2);
        end
        plot(start_vs_peak(thisdir(sorti(diri2))), plottick, '>k', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        plot(stop_vs_peak(thisdir(sorti(diri2))), plottick, '<k', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        plottick = plottick+1;
    end
    plottick = plottick+4;
    if diri < ndir
        plot([-Settings.Raster.winaround, Settings.Raster.winaround], [plottick-2.5, plottick-2.5], '-k');
    end
end

plot([0, 0], [0, plottick], '--k');

yticks([]);
ylim([0, plottick-2.5]);
xlim([-0.5 0.5]);
xticks(-1:0.5:1);
xticklabels([]);
set(gca,'TickDir', 'out');
set(gca,'TickLength',[0.035 0.01]);
set(gca, 'YColor', 'w');
set(gca,'YDir','reverse')
end