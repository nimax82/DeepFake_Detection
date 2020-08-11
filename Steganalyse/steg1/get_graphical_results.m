cfa_orig_seen = 0;
cfa_orig_new = 0;
cfa_fake_pnew = 0;
cfa_fake_tnew = 0;

color_orig_seen = 0;
color_orig_new = 0;
color_fake_pnew = 0;
color_fake_tnew = 0;

%x(end+1) = newval

for i=1:136
    if graphic_labels(i) == 1
        cfa_orig_seen(end+1) = mean_score_cfa(i);
        color_orig_seen(end+1) = mean_score_color(i);
    elseif graphic_labels(i) == 2
        cfa_orig_new(end+1) = mean_score_cfa(i);
        color_orig_new(end+1) = mean_score_color(i);
    elseif graphic_labels(i) == 3
        cfa_fake_pnew(end+1) = mean_score_cfa(i);
        color_fake_pnew(end+1) = mean_score_color(i);
    elseif graphic_labels(i) == 4
        cfa_fake_tnew(end+1) = mean_score_cfa(i);
        color_fake_tnew(end+1) = mean_score_color(i);
    end
end

cfa_orig_seen(1) = [];
cfa_orig_new(1) = [];
cfa_fake_pnew(1) = [];
cfa_fake_tnew(1) = [];

color_orig_seen(1) = [];
color_orig_new(1) = [];
color_fake_pnew(1) = [];
color_fake_tnew(1) = [];

cfa_orig_seen(length(cfa_orig_seen): 73) = 0;
cfa_orig_new(length(cfa_orig_new): 73) = 0;
cfa_fake_tnew(length(cfa_fake_tnew): 73) = 0;

color_orig_seen(length(color_orig_seen): 73) = 0;
color_orig_new(length(color_orig_new): 73) = 0;
color_fake_tnew(length(color_fake_tnew): 73) = 0;



figure
%plot(cfa_orig_seen, cfa_orig_new,...
%     cfa_fake_pnew, cfa_fake_tnew);
plot(cfa_orig_seen);
hold on
plot(cfa_orig_new);
hold on
plot(cfa_fake_pnew);
hold on
plot(cfa_fake_tnew);
legend('original & seen', 'orignal & new','fake & partially new', 'fake & totally new')
title('feature [19] - CFA')
xlabel('video')
ylabel('score value (mean)')

figure
plot(color_orig_seen);
hold on
plot(color_orig_new);
hold on
plot(color_fake_pnew);
hold on
plot(color_fake_tnew);
legend('original & seen', 'orignal & new','fake & partially new', 'fake & totally new')
title('feature [16] - Color')
xlabel('video')
ylabel('score value (mean)')
