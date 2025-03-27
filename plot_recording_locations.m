function plot_recording_locations(Settings)

locs.center = 2;
locs.anterior = 1;
locs.posterior = 3;
locs.medial = 4;

lwid = 0.2;
lgap = 0.05;

figure('Renderer','Painters');
hold on;
for typei = 1:Settings.Global.Locs.ntracktypes
    plot(locs.(Settings.Global.Locs.possiblenames{typei})(1)*[1 1], [0 14], '-b');
end

ntype = length(Settings.Global.Locs.possiblenames);

nall = length(Settings.Global.Locs.depths);

% Plot the highlights first
hlh = 0.1;
hlw = 0.1;

allID = Settings.Global.Locs.tracks;

for neui = 1:Settings.Global.allnneu
    allID{neui} = [allID{neui} num2str(Settings.Global.Locs.depths(neui))];
end

alli = 1;
while alli <= nall    
    xloc = locs.(Settings.Global.Locs.tracks{alli});
    
    % check ahead 4 to see if it's a multi channel
    nmulti = 1;
    for checki = 1:4
        if (alli <= nall-checki) && nmulti == checki
            if (Settings.Global.Locs.depths(alli+checki) == Settings.Global.Locs.depths(alli)) && strcmp(allID{alli+checki}, allID{alli})
                nmulti = checki+1;
            end
        end
    end

    colstr = '-k';
    if nmulti == 3
        plot([xloc-lwid/2 xloc+lwid/2], [Settings.Global.Locs.depths(alli) Settings.Global.Locs.depths(alli)], colstr);
        plot([xloc-lwid/2 xloc+lwid/2]-(lwid+lgap), [Settings.Global.Locs.depths(alli) Settings.Global.Locs.depths(alli)], colstr);
        plot([xloc-lwid/2 xloc+lwid/2]+(lwid+lgap), [Settings.Global.Locs.depths(alli) Settings.Global.Locs.depths(alli)], colstr);
    elseif nmulti == 2
        plot([xloc-lwid/2 xloc+lwid/2]-(lwid+lgap)/2, [Settings.Global.Locs.depths(alli) Settings.Global.Locs.depths(alli)], colstr);
        plot([xloc-lwid/2 xloc+lwid/2]+(lwid+lgap)/2, [Settings.Global.Locs.depths(alli) Settings.Global.Locs.depths(alli)], colstr);
    else
        plot([xloc-lwid/2 xloc+lwid/2], [Settings.Global.Locs.depths(alli) Settings.Global.Locs.depths(alli)], colstr);
    end
    
    alli = alli+nmulti;
end


xlim([0.5, 4.5]);
ylim([0, 13]);


xticks(1:4);
yticks([-5 0 5 10]);

xticklabels({'Ant', 'Ctr', 'Post', 'Med'});

set(gca,'TickLength',[0.035 0.01]);
set(gca,'TickDir', 'out')
set(gca,'GridColor',[0 0 0]);
set(gca,'LineWidth',0.5);

plot(locs.center, 0, 'ok', 'MarkerSize', 20, 'LineWidth', 4);

grid on