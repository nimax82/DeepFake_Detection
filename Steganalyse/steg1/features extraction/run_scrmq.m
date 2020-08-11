clc;
clear;
%'/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images/training';
folder_path = '/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-synthesis/images/test_base';

save_color = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/staganalysis/features/test/color/color_test.mat';
save_cfa = '/Users/test/Documents/Stage_Maxime/DNGI_Detection/staganalysis/features/test/cfa/cfa_test.mat';
%folder_path = '/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-synthesis/images/test_base';

%features;% = cell(numel(folder_path), 100);

folder_container = dir(fullfile(folder_path, 'id*'));
nb_video = numel(folder_container);
frame_by_video = 10;
%features_color = NaN(nb_video * frame_by_video,2073);
%features_cfa = NaN(nb_video * frame_by_video,10323);
features_color = load(save_color).features_color;
features_cfa = load(save_cfa).features_cfa;

% folder_container = dir(fullfile(folder_path, 'id*'));
% %N = 10;
% N = numel(folder_container);
% features_color = repmat( struct( 'x', 1 ), N, 1 );
% features_cfa = repmat( struct( 'x', 1 ), N, 1 );
% for ii=1:N
%     features_color(ii).x(1)=1;  
%     features_cfa(ii).x(1)=1;
% end  

for i = 1  : nb_video
    %Grab folder that include sets of images
    %folder = folder_path(i);
    folder = folder_container(i);
%     if folder.name(1) == '.'
%         %avoid to compute './' and '../'
%         continue
%     end
    
    filenames=dir(fullfile(folder.folder, folder.name, '*.tiff'));
    local_nb_video = numel(filenames);
    assert(frame_by_video == local_nb_video);
    for n = 1:local_nb_video
        if (filenames(n).name(1) == '.')
            continue
        end
        fullname=fullfile(filenames(n).folder, filenames(n).name);
        X = imread(fullname);
        %imshow(X);
        %tic
        features_color((i-1) * frame_by_video + n, :) = struct2array(SCRMQ1(double(X),2,'color'));
        features_cfa((i-1) * frame_by_video + n, :) = struct2array(SRMQ1C2cfaRGGB(double(X), 2,'color'));
        %toc 
        save(save_color, 'features_color')
        save(save_cfa, 'features_cfa')
    end
end

% id_to_delete = ones(1,N);
% for i = 1 : N
%     id_to_delete(i) = isstruct(features(i).x);
% end
% structfun(@(x) x('x' == 1), features, 'UniformOutput', false)
%features = rmfield(features,  

% X = imread('/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images/training/id0_0002/0029.tiff');
% tic
% SCRMQ1Features = SCRMQ1(double(X));
% toc