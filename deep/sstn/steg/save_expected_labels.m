%clc;
%clear;
%'/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images/training';
folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/test';
video_labels_path = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/deep/sstn/steg/features/test/test_labels_video.mat';
frame_labels_path = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/deep/sstn/steg/features/test/test_labels_frame.mat';


nb_frame = 15;
folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
expected_labels_video = NaN(nb_video, 1);
expected_labels_frame = NaN(nb_video * nb_frame, 1);

for i = 1  : nb_video
    %Grab folder that include sets of images
    %folder = folder_path(i);
    folder = folder_container(i);
%     if folder.name(1) == '.'
%         %avoid to compute './' and '../'
%         continue
%     end

    if count(folder.name, 'id') == 1
        expected_labels_video(i) = 0;
    elseif count(folder.name, 'id') == 2
        expected_labels_video(i) = 1;
    end
end

label = expected_labels_video(1);
for i=1:nb_video * nb_frame
    expected_labels_frame(i) = label;
    if mod(i, nb_frame) == 0
        if i == nb_video * nb_frame
            break;
        end
        label = expected_labels_video((i / nb_frame) + 1);
    end  
end

save(video_labels_path, 'expected_labels_video')
save(frame_labels_path, 'expected_labels_frame')