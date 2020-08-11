clear;
clc;
svm_prediction = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/steg/prediction/P_svm_prediction.mat').label_color;
deep_prediction = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/deep/predictions/P_deep_predictions.mat').classification_labels_frame;
expected_labels_frame = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/labels/test_labels_frame.mat').expected_labels_frame;
fusion_prediction = NaN(size(svm_prediction));

nb_frame = 15;
count_middle =0;
count_mid_svm =0;
count_mid_deep =0;
middle_labels = NaN(517,1);
count_svm_right = 0;
count_deep_right = 0;

j = 0;
for i=1 : numel(svm_prediction)
    current_fusion = svm_prediction(i) + deep_prediction(i);
    fusion_prediction(i) = double(current_fusion > 1);
    if current_fusion == 1
        j = j + 1;
        count_middle = count_middle + 1;
        count_mid_svm = count_mid_svm + svm_prediction(i);
        count_mid_deep = count_mid_deep + deep_prediction(i);
        middle_labels(j) =  expected_labels_frame(i);
        if expected_labels_frame(i) == svm_prediction(i)
            count_svm_right = count_svm_right + 1;
        else
            count_deep_right = count_deep_right + 1;
        end
    end
end

%frame prediction to video
expected_labels_video = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/labels/test_labels_video.mat').expected_labels_video;
result_labels_fusion = NaN(numel(expected_labels_video), 1);
sum_fusion = 0;
j = 0;
for i=1 : numel(fusion_prediction)
    sum_fusion = sum_fusion + fusion_prediction(i);
    if mod(i,nb_frame) == 0
        j = j + 1;
        result_labels_fusion(j) = int8(sum_fusion > (nb_frame / 2));
        sum_fusion = 0;
    end
end

% compare

similarity = (result_labels_fusion == expected_labels_video);
result_fusion = sum(double(similarity));