net = load('deep/networks/all_frame_network').net;
inputSize = net.Layers(1).InputSize;

folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/I_P_Compare/all/test';
folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
nb_frame = 15;

%expected_labels_video = load(labels_path + 'test_labels_video.mat').expected_labels_video;
%load('steg/features/test/test_labels_frame.mat');
expected_labels_frame = load('labels/test_labels_frame.mat').expected_labels_frame;

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