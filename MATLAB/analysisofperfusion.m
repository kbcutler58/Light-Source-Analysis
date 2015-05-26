


troughs = zeros(size(locs3),1);

for m = 1:size(locs3);
    [~,troughs(:,m)] = min(abs(mlocs- locs3(m)));
    m+1;
end

