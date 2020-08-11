%clc;
%clear;
%'/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images/training';
folder_path = '/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-synthesis/images/test_base';
save_labels_path = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/staganalysis/features/expected_labels/labels.mat';



folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
expected_labels = NaN(nb_video, 1);

for i = 1  : nb_video
    %Grab folder that include sets of images
    %folder = folder_path(i);
    folder = folder_container(i);
%     if folder.name(1) == '.'
%         %avoid to compute './' and '../'
%         continue
%     end

    if count(folder.name, 'id') == 1
        expected_labels(i) = 0;
    elseif count(folder.name, 'id') == 2
        expected_labels(i) = 1;
    end
end