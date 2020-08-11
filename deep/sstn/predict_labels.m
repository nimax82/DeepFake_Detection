% path_image = "/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/svm/test/id29_id31_0009/face/0014.tiff";
% I = imread(path_image);
% I = imresize(I,inputSize(1:2));
% [label,scores] = classify(net,I);

folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/test';
folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
nb_frame = 15;
inputSize = net.Layers(1).InputSize;


labels_path = 'steg/features/test/';
%expected_labels_video = load(labels_path + 'test_labels_video.mat').expected_labels_video;
%load('steg/features/test/test_labels_frame.mat');
expected_labels_frame = load('steg/features/test/test_labels_frame.mat').expected_labels_frame;

classification_labels_frame = NaN(nb_video * nb_frame, 1);
classification_labels_video = NaN(nb_video, 1);

for i = 1  : nb_video
    %Grab folder that include sets of images
    folder = folder_container(i);
    images_root = strcat(folder.folder, "/" , folder.name);
    filenames=dir(fullfile(images_root, '0*.tiff'));
    % classify all frames of the current video
    for n = 1:nb_frame
        fullname=fullfile(filenames(n).folder, filenames(n).name);
        X = imread(fullname);
        X = imresize(X,inputSize(1:2));
        [label,scores] = classify(net,X);
        classification_labels_frame((i - 1) * nb_frame + n) = double(label) - 1;
    end
end

similarity_frame = (classification_labels_frame == expected_labels_frame);
result = sum(double(similarity_frame));
deep_predictions = classification_labels_frame;