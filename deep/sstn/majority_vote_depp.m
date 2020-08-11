%expected_labels = load('features/test/test_labels.mat');
%test_labels = expected_labels.expected_labels;

classif_labels_video = NaN(nb_video, 1);
expected_labels_video = NaN(nb_video, 1);
nb_frame = 15;


sum_e = 0;
sum_c = 0;

j = 0;
for i=1 : numel(classification_labels_frame)
    sum_e = sum_e + expected_labels_frame(i);
    sum_c = sum_c + classification_labels_frame(i);
    if mod(i,nb_frame) == 0
        j = j + 1;
        classif_labels_video(j) = int8(sum_c > (nb_frame / 2));
        expected_labels_video(j) = int8(sum_e > (nb_frame / 2));
        sum_e = 0;
        sum_c = 0;
    end
end

% compare

mj_similarity = (classif_labels_video == expected_labels_video);
mj_result = sum(double(mj_similarity));





