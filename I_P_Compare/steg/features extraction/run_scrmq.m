clc;
clear;
%folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/I_P_Compare/all/train_val';
%folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/I_P_Compare/all/test';
folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/faceForensics/dataSet_c23/gen_test/images/P';

%save_color = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/steg/features/training/all/color.mat';
save_color = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/steg/features/test/P/Gen_color.mat';

folder_container = dir(fullfile(folder_path, '[fo]_*'));
nb_video = numel(folder_container);
frame_by_video = 15;
features_color = NaN(nb_video * frame_by_video,18157);%2073
%features_color = load(save_color).features_color;

 
tic
for i = 1  : nb_video
    i
    %Grab folder that include sets of images
    folder = folder_container(i);
    images_root = strcat(folder.folder, "/" , folder.name);
    filenames=dir(fullfile(images_root, '0*.tiff'));
    local_nb_frame = numel(filenames);
    assert(frame_by_video <= local_nb_frame);
    for n = 1:frame_by_video
%         if (filenames(n).name(1) == '.')
%             continue
%         end
        fullname=fullfile(filenames(n).folder, filenames(n).name);
        X = imread(fullname);
        %imshow(X);
        %tic
        features_color((i-1) * frame_by_video + n, :) = struct2array(SCRMQ1(double(X)));
        %toc 
        save(save_color, 'features_color')
    end
end
toc