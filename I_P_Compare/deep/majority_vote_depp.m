expected_labels_frame = load('labels/test_labels_frame.mat').expected_labels_frame;
classification_labels = load('deep/predictions/I_deep_predictions').classification_labels_frame;

nb_video = 100;
mj_classif_labels = NaN(nb_video, 1);
mj_expected_labels = NaN(nb_video, 1);
nb_frame = 15;


sum_e = 0;
sum_c = 0;

j = 0;
for i=1 : numel(classification_labels)
    sum_e = sum_e + expected_labels_frame(i);
    sum_c = sum_c + classification_labels(i);
    if mod(i,nb_frame) == 0
        j = j + 1;
        mj_classif_labels(j) = int8(sum_c > (nb_frame / 2));
        mj_expected_labels(j) = int8(sum_e > (nb_frame / 2));
        sum_e = 0;
        sum_c = 0;
    end
end

% compare

mj_similarity = (mj_classif_labels == mj_expected_labels);
mj_result = sum(double(mj_similarity));





