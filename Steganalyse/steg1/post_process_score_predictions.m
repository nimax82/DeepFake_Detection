
mean_score_cfa = NaN(numel(score_color) / 10, 1);
mean_score_color = NaN(numel(score_color) / 10, 1);

local_sum_cfa = 0;
local_sum_color = 0;

for i=1 : numel(score_color)
    local_sum_cfa = local_sum_cfa + score_cfa(i);
    local_sum_color = local_sum_color + score_color(i);
    if mod(i,10) == 0
        mean_score_cfa(i/10) = local_sum_cfa / 10;
        mean_score_color(i/10) = local_sum_color / 10;
        local_sum_cfa = 0;
        local_sum_color = 0;
    end
end
 
        