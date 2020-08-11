%imported_cfa = load('labels/class_labels_cfa.mat');
%labels_cfa = imported_cfa.label_cfa;

imported_color = load('sstn/prediction/svm_predictions.mat');
svm_predictions = imported_color.label_color;

imported_deep = load('sstn/prediction/deep_predictions.mat');
deep_predictions = imported_deep.deep_predictions;

nb_frame_per_video = 15;
total_nb_frame = size(svm_predictions,1);
nb_blocs = total_nb_frame / 5;
nb_video = nb_blocs / 3;

label_frames = NaN(total_nb_frame, 1);

label_blocs = NaN(nb_blocs, 1);

label_videos = NaN(nb_video, 1);

sum_bloc = 0;
for i=1 : total_nb_frame
    value_frame = svm_predictions(i) + deep_predictions(i);
    label_frames(i) = double(value_frame > 1);
    sum_bloc = sum_bloc + label_frames(i);
    if mod(i,5) == 0
        label_blocs(i/5) = double(sum_bloc >= 3);
        sum_bloc = 0;
    end
end


sum_video = 0;
for i=1 : nb_blocs
    sum_video = sum_video + label_blocs(i);
    if mod(i,3) == 0
        label_videos(i/3) = double(sum_video >= 2);
        sum_video = 0;
    end
end

import_expected_labels = load('sstn/steg/features/test/test_labels_frame');
expected_labels_frame = import_expected_labels.expected_labels_frame;

expected_labels_video = NaN(nb_video, 1);
for i=1 : total_nb_frame
    if mod(i,nb_frame_per_video) == 0
        expected_labels_video(i/nb_frame_per_video) = expected_labels_frame(i);
    end
end

similarity = (label_videos == expected_labels_video);
result = sum(double(similarity));
