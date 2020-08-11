clc;
clear;
% Estimate and compare PRNU from two sets of images

% --- Paremeters --- %

mode = 't';
data_path = '/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/Celeb-real/images';

if strcmp(mode, 's')
    path_mode = strcat(data_path, '/spacial');
else
    path_mode = strcat(data_path, '/temporal');
end

videos_folder = dir(path_mode);

% --- Main loop --- %
results = struct;
for i = 1 : numel(videos_folder)
    %Grab folder that include the two sets of images
    folder = videos_folder(i);
    if folder.name(1) == '.'
        %avoid to compute './' and '../'
        continue
    end
    %Initialization of the data sets
    [set1, set2] = open_images(path_mode, folder.name);
    %Compute prnu for each set
    fg1 = getFingerPrintFromSet(set1);
    fg2 = getFingerPrintFromSet(set2);
    %Statistic comparison
    C = crosscorr(fg1, fg2);
    imshow(uint8(255*fg1));
    detection = PCE(C);
    %Save result as struct element
    results.(folder.name) = get_struct_result(detection);
end 

% Exportation of the results

%handle file openning
file_name = strcat('results_', mode, '.json');
edit results.json
fileID = fopen(file_name,'w');
%format results as json structure
results_f = jsonencode(results);
%write in the file
fprintf(fileID,results_f);
fclose(fileID);

% --- Internal functions --- %

%Read and store all image paths in two distinct structures (one for each patchs)
function [set1, set2] = open_images(data_path, video_id)
%Hold all images path in two distinct arrays
sets = cell(2,100);
%Read all images on the two sub-folders
for id = 1:2
    %construct sub-folder name
    patch_id = strcat('patch', int2str(id-1));
    current_folder = strcat(data_path, '/', video_id, '/', patch_id);
    %read content
    filenames=dir(fullfile(current_folder,'*.png'));
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
set1 = sets(1,:);
set2 = sets(2,:);
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