


expected_labels = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/labels/gen_raw_labels_video.mat');
test_labels = expected_labels.expected_labels_video;

label_color = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/steg/prediction/raw_svm_prediction.mat').label_color;

nb_frame = 15;
predict_labels_color_video = NaN(numel(label_color) / nb_frame, 1);

sum_color = 0;

j = 0;
for i=1 : numel(label_color)
    sum_color = sum_color + label_color(i);
    if mod(i,nb_frame) == 0
        j = j + 1;
       
        predict_labels_color_video(j) = int8(sum_color > (nb_frame / 2));
       
        sum_color = 0;
    end
end

% compare

similarity_color = (predict_labels_color_video == test_labels);
result_color = sum(double(similarity_color));



