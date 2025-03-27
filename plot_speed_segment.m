function plot_speed_segment(Settings, D, dati, usetime)

gkernstr = ['SmoothG' num2str(1000*Settings.Plot.gkern,'%0.2i')];

figure('Position', [10 10 2400 200]);
clf;
hold on;
plot(D(dati).Time.ktime, D(dati).Kin.(gkernstr).Speed, 'LineWidth', 3);

arrow = D(dati).Time.dircue(find(D(dati).Time.dircue > usetime, 1, 'first'));
move = D(dati).Time.gocue(find(D(dati).Time.gocue > usetime, 1, 'first'));
back = D(dati).Time.gocue(find(D(dati).Time.gocue > usetime+2, 1, 'first'));

xticks([arrow move back]);
xticklabels({})

xlim([move-2 back+3]);

yticks([0 200 400]);

ax = gca;
ax.YAxis.FontSize = 16;

ylabel('Fingertip Speed (mm/s)', 'FontSize', 20);

patch([arrow arrow arrow+0.1 arrow+0.1], [200 400 400 200], 'k', 'LineStyle', 'none');