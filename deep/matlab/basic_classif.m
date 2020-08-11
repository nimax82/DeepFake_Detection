% path_image = "/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/svm/test/id29_id31_0009/face/0014.tiff";
% I = imread(path_image);
% I = imresize(I,inputSize(1:2));
% [label,scores] = classify(net,I);

folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/svm/test';
folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
nb_frame = 15;
%expected_labels_init = NaN(nb_video, 1);
expected_labels = NaN(nb_video * nb_frame, 1);

classification_labels = NaN(nb_video * nb_frame, 1);

for i = 1  : nb_video
    %init
    expected_video_label = NaN;
    %Grab folder that include sets of images
    folder = folder_container(i);
    images_root = strcat(folder.folder, "/" , folder.name);
    filenames=dir(fullfile(images_root, 'face', '0*.tiff'));
    local_nb_frame = numel(filenames);
    assert(nb_frame <= local_nb_frame);
    % save_expected label for current video
    if count(folder.name, 'id') == 1
        expected_video_label = 0;
    elseif count(folder.name, 'id') == 2
        expected_video_label = 1;
    end
    % classify all frames of the current video
    for n = 1:nb_frame
        fullname=fullfile(filenames(n).folder, filenames(n).name);
        X = imread(fullname);
        X = imresize(X,inputSize(1:2));
        [label,scores] = classify(net,X);
        classification_labels((i - 1) * nb_frame + n) = double(label) - 1;
        expected_labels((i - 1) * nb_frame + n) = expected_video_label;
    end
end