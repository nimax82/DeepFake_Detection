clear;
clc;
svm_prediction = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/steg/prediction/P_svm_prediction.mat').label_color;
deep_prediction = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/deep/predictions/P_deep_predictions.mat').classification_labels_frame;
expected_labels_frame = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/labels/test_labels_frame.mat').expected_labels_frame;

nb_video = 100;
nb_frame = 15;

fusion_prediction = NaN(nb_video, 1);

for i=1 : nb_video
    proba_svm = 0;
    proba_deep = 0;
    for j=1 : nb_frame
        proba_svm = proba_svm + svm_prediction((i-1)*nb_frame + j);
        proba_deep = proba_deep + deep_prediction((i-1)*nb_frame + j);
    end
    confi_svm = abs(proba_svm - nb_frame/2);
    confi_deep = abs(proba_deep - nb_frame/2);
    
    if (confi_svm >= confi_deep)
        fusion_prediction(i) = double(proba_svm >= nb_frame/2);
    else
        fusion_prediction(i) = double(proba_deep >= nb_frame/2);
    end
end

expected_labels_video = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/labels/test_labels_video.mat').expected_labels_video;
similarity = (fusion_prediction == expected_labels_video);
result_fusion = sum(double(similarity));