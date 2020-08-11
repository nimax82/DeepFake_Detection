clc;
clear;
% Estimate and compare PRNU from two sets of images
%vor2
% --- Paremeters --- %

%Celeb-synthesis
data_path = '/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images';

path_mode = strcat(data_path, '/byGroup');

videos_folder = dir(path_mode);
nb_grp = 8;
nb_comp = floor((nb_grp * (nb_grp-1)) / 2);

% --- Main loop --- %
results = cell(numel(videos_folder), floor((nb_grp * (nb_grp-1)) / 2) + 1);
for i = 1 : numel(videos_folder)
    %Grab folder that include sets of images
    folder = videos_folder(i);
    if folder.name(1) == '.'
        %avoid to compute './' and '../'
        continue
    end
    %Initialization of the data sets
    sets = open_images(path_mode, folder.name, nb_grp);
    %Compute prnu for each set
    x_size = size(imread(sets{1}),1);
    y_size = size(imread(sets{1}),2);
    fg = zeros(nb_grp, x_size, y_size); 
    for j = 1:nb_grp
        fg(j,:,:) = getFingerPrintFromSet(sets(j));
    end
    %Statistic comparison
    detection_array = cell(nb_comp, 6);
    offset = 0;
    for c = 1:nb_grp-1
        fg_a = reshape(fg(c,:,:), x_size, y_size);
        for k = c+1:nb_grp
            fg_b = reshape(fg(k,:,:), x_size, y_size);
            crossC = crosscorr(fg_a, fg_b);
            index = offset + (k - c);
            %a = PCE(crossC)
            detection_array(index, :) = struct2cell(PCE(crossC));
        end
        offset = offset + nb_grp - c;
    end
  
    results(i,1) = {folder.name};
    test = reshape(detection_array(:,3), 1, nb_comp);
    for z = 1:nb_comp
        results(i, z+1) = test(z);
    end    
    %results(i)[2,29] = folder.name;
    %imshow(uint8(255*fg1));
    %Save result as struct element
    %results.(folder.name) = get_struct_result(detection);
end 

%id_to_remove = cellfun('isempty',results);
%id_to_remove = id_to_remove(:,1);
results(1,:) = [];
results(1,:) = [];
%results = results(cellfun('isempty',results) == 0)
%results = results(~cellfun('isempty',results));

% % Exportation of the results
% 
% %handle file openning
% file_name = strcat('results_', mode, '.json');
% edit results.json
% fileID = fopen(file_name,'w');
% %format results as json structure
% results_f = jsonencode(results);
% %write in the file
% fprintf(fileID,results_f);
% fclose(fileID);

% --- Internal functions --- %

%Read and store all image paths structure
function sets = open_images(data_path, video_id, nb_gpr)
    %Hold all images path in arrays
    sets = cell(nb_gpr,20);
    %Read all images on the two sub-folders
    for id = 1:nb_gpr
        %construct sub-folder name
        patch_id = strcat('patch', int2str(id-1));
        current_folder = strcat(data_path, '/', video_id, '/', patch_id);
        %read content
        filenames=dir(fullfile(current_folder,'*.tiff'));
        set = cell(1, (numel(filenames)));
        %store content in 'set' value
        for n = 1:numel(filenames)
            if (filenames(n).name(1) == '.')
                continue
            end
            fullname=fullfile(current_folder,filenames(n).name);
            set(n) = cellstr(fullname);
        end
        set = set(~cellfun('isempty',set));
        sets(id,:) = set;
    end
%set1 = sets(1,:);
%set2 = sets(2,:);
end

%Compute fingerPrint from a set of images
function fingerPrint = getFingerPrintFromSet(Images)
    RP = getFingerprint(Images);
    RP = rgb2gray1(RP);
    sigmaRP = std2(RP);
    fingerPrint = WienerInDFT(RP,sigmaRP);
end

%Build and return a struct storing video_id and prnu comparison result
function res_struct = get_struct_result(detection)
    field1 = 'peakheight';
    field2 = 'PeakLocation';
    field3 = 'PCE';
    field4 = 'pvalue';
    field5 = 'P_FA';
    field6 = 'log10P_FA';
    
    res_struct = struct(field1, detection.peakheight ,...
                        field2, detection.PeakLocation ,...
                        field3, detection.PCE ,...
                        field4, detection.pvalue ,...
                        field5, detection.P_FA ,...
                        field6, detection.log10P_FA ...
                        ); 
end
