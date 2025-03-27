function normdat = my_normalize(dat)

normdat = nan(size(dat));
mdat = nanmean(dat,1);
sdat = nanstd(dat,0,1);
for dimi = 1:size(dat,2)
    normdat(:,dimi) = (dat(:,dimi)-mdat(dimi))/sdat(dimi);
end