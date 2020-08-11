clear
clc;
svm_prediction = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/steg/prediction/All_svm_prediction.mat').label_color;
deep_prediction = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/deep/predictions/all_deep_predictions.mat').classification_labels_frame;
fusion_prediction = NaN(size(svm_prediction));

nb_frame = 15;

for i=1 : numel(svm_prediction)
    current_fusion = svm_prediction(i) + deep_prediction(i);
    fusion_prediction(i) = double(current_fusion >= 1);
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