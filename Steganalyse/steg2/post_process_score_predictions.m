
mean_score_cfa = NaN(numel(score_cfa) / 22, 1);
mean_score_color = NaN(numel(score_color) / 22, 1);

local_sum_cfa = 0;
local_sum_color = 0;

for i=1 : numel(score_color)
    local_sum_cfa = local_sum_cfa + score_cfa(i);
    local_sum_color = local_sum_color + score_color(i);
    if mod(i,22) == 0
        mean_score_cfa(i/22) = local_sum_cfa / 22;
        mean_score_color(i/22) = local_sum_color / 22;
        local_sum_cfa = 0;
        local_sum_color = 0;
    end
end
 
        