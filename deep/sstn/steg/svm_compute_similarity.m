


%expected_labels = load('features/test/test_labels.mat');
%test_labels = expected_labels.expected_labels;

%test_labels = expected_labels;

nb_frame = 15;
predict_labels_color_video = NaN(numel(label_color) / nb_frame, 1);
predict_labels_cfa_video = NaN(numel(label_cfa) / nb_frame, 1);

sum_color = 0;
sum_cfa = 0;

j = 0;
for i=1 : numel(label_color)
    sum_color = sum_color + label_color(i);
    sum_cfa = sum_cfa + label_cfa(i);
    if mod(i,nb_frame) == 0
        j = j + 1;
       
        predict_labels_color_video(j) = int8(sum_color > (nb_frame / 2));
        predict_labels_cfa_video(j) = int8(sum_cfa > (nb_frame / 2));
       
        sum_color = 0;
        sum_cfa = 0;
        
    end
end

% compare

similarity_color = (predict_labels_color_video == expected_labels_video);
result_color = sum(double(similarity_color));

similarity_cfa = int8((predict_labels_cfa_video == expected_labels_video));
result_cfa = sum(double(similarity_cfa));




