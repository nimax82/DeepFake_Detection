%clc;
%clear;
%folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/training';
folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/test';
%save_labels_path = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/stag3/features/training/labels.mat';



folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
nb_frame = 15;
expected_labels_init = NaN(nb_video, 1);
expected_labels = NaN(nb_video * nb_frame, 1);

for i = 1  : nb_video
    %Grab folder that include sets of images
    %folder = folder_path(i);
    folder = folder_container(i);
%     if folder.name(1) == '.'
%         %avoid to compute './' and '../'
%         continue
%     end

    if count(folder.name, 'id') == 1
        expected_labels_init(i) = 0;
    elseif count(folder.name, 'id') == 2
        expected_labels_init(i) = 1;
    end
end

for j=1 : nb_video * nb_frame
    id = ceil(j/nb_frame)
    expected_labels(j) = expected_labels_init(id);
end