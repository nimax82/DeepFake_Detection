num_frame = 10;
mean_score_cfa = NaN(numel(score_cfa) / num_frame, 1);
mean_score_color = NaN(numel(score_color) / num_frame, 1);

local_sum_cfa = 0;
local_sum_color = 0;

for i=1 : numel(score_color)
    local_sum_cfa = local_sum_cfa + score_cfa(i);
    local_sum_color = local_sum_color + score_color(i);
    if mod(i,num_frame) == 0
        mean_score_cfa(i/num_frame) = local_sum_cfa / num_frame;
        mean_score_color(i/num_frame) = local_sum_color / num_frame;
        local_sum_cfa = 0;
        local_sum_color = 0;
    end
end
 
        