clc;
clear;
%folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/training';
folder_path = '/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/test';

%save_color = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/stag3/features/training/color/color.mat';
%save_cfa = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/stag3/features/training/cfa/cfa.mat';
save_color = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/stag3/features/test/color/color_testv2.mat';
save_cfa = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/stag3/features/test/cfa/cfa_testv2.mat';

%features;% = cell(numel(folder_path), 100);

folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
frame_by_video = 15;
%features_color = NaN(nb_video * frame_by_video,18157);%2073
%features_cfa = NaN(nb_video * frame_by_video,10323);
features_color = load(save_color).features_color;
features_cfa = load(save_cfa).features_cfa;
 

for i = 2  : nb_video
    %Grab folder that include sets of images
    %folder = folder_path(i);
    folder = folder_container(i);
%     if folder.name(1) == '.'
%         %avoid to compute './' and '../'
%         continue
%     end
    images_root = strcat(folder.folder, "/" , folder.name);
    filenames=dir(fullfile(images_root, 'face', '0*.tiff'));
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
        features_cfa((i-1) * frame_by_video + n, :) = struct2array(SRMQ1C2cfaRGGB(double(X), 2,'color'));
        %toc 
        save(save_color, 'features_color')
        save(save_cfa, 'features_cfa')
    end
end