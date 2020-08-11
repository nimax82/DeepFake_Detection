clc;
clear;

% --- Forensic DataBase ---

% data_path = '/Volumes/VERBATIM HD/Stage_Maxime/faceForensics/dataSet_c23';
% orig_folder = strcat(data_path, "/", "original_sequences/actors/c23/images/03__talking_angry_couch");
% deep_folder = strcat(data_path, "/", "manipulated_sequences/DeepFakeDetection/c23/images/03_04__talking_angry_couch__T04P6ELC");
% 
% current_folder = orig_folder;
% images = dir(fullfile(current_folder, '0*.png'));
% orig_peaks = get_peaks(images);
% 
% current_folder = deep_folder;
% images = dir(fullfile(current_folder, '0*.png'));
% fake_peaks = get_peaks(images);

% --- Celeb DataBase ---

data_path = '/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2';
deep_folder = strcat(data_path, '/Celeb-synthesis/images/fft/id6_id16_0004');
orig_folder = strcat(data_path, '/Celeb-real/images/fft/id6_0004');

current_folder = orig_folder;
images = dir(fullfile(strcat(current_folder, '/patch0'), '0*.png'));
images = [images; dir(fullfile(strcat(current_folder, '/patch1'), '0*.png'))];
orig_peaks = get_peaks(images);

current_folder = deep_folder;
images = dir(fullfile(strcat(current_folder, '/patch0'), '0*.png'));
images = [images; dir(fullfile(strcat(current_folder, '/patch1'), '0*.png'))];
fake_peaks = get_peaks(images);


function peaks = get_peaks(images)

    high_pass_Kernel = [1,0,-1;2,0,-2;1,0,-1]; 

    peaks = zeros(numel(images) - 1, 1);
    
    img1 = imread(strcat(images(1).folder, "/" ,images(1).name));
    img1 = rgb2gray(img1);
    img1 = uint8(conv2(img1, high_pass_Kernel, 'same'));
    last_fft = fftshift(fft2(double(img1)));
    last_fft = 255 * mat2gray(abs(last_fft));
    
    img2 = imread(strcat(images(2).folder, "/", images(2).name));
    img2 = rgb2gray(img2);
    current_fft = fftshift(fft2(double(img2)));
    current_fft = 255 * mat2gray(abs(current_fft));
    
    %diff = abs(current_fft - last_fft);
    %peak = mean(mean(diff));
    %peaks(1) = mean(mean(diff));
    
    peaks(1) = sum(sum((round(current_fft, 3) == round(last_fft,3))));
    
    
    last_fft = current_fft;
    for img_id = 3 : numel(images)
        img = imread(strcat(images(img_id).folder, "/", images(img_id).name));
        img = rgb2gray(img);
        img = uint8(conv2(img, high_pass_Kernel, 'same'));
        current_fft = fftshift(fft2(double(img)));
        current_fft = 255 * mat2gray(abs(current_fft));
        
        %diff = abs(current_fft - last_fft);
        
        %sum
        %peaks(img_id - 1) = sum(sum(diff));
        
        %mean
        %peaks(img_id - 1) = mean(mean(diff));
        
        %peak = max(peak, mean(mean(diff)));
        
        %diff max
        %max_diff = max(diff, max_diff);
        
        %nb differences
        peaks(img_id - 1) = sum(sum((round(current_fft, 3) == round(last_fft,3))));
        
        last_fft = current_fft;
    end
    %peak = mean(mean(max_diff));
end