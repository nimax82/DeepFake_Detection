clc;
clear;

data_path = '/Volumes/VERBATIM HD/Stage_Maxime/Celeb-DF-v2/image_couples';
videos_folder = dir(data_path);

nb_frame = 2;

%create frame clusters
clusters = struct;
for i = 1 : numel(videos_folder)
    %Grab folder that include the two sets of images
    folder = videos_folder(i);
    if folder.name(1) == '.'
        %avoid to compute './' and '../'
        continue
    end
    clusters.(folder.name) = open_images(data_path, folder.name, nb_frame);
end


field_names = fieldnames(clusters);
images_result = cell(nb_frame * length(field_names));
for i=1 : length(field_names)
    current_clusters_name = char(field_names(i));
    %for frame_num=1 : nb_frame
        frame_paths = getfield(clusters, current_clusters_name);
        ffts = get_ffts(frame_paths, nb_frame);%return list of ffts (cell)
        result = format_fft(frame_paths, ffts, nb_frame);
    %end
    
    
end


% --- Internal functions --- %

function clusters = open_images(data_path, video_id, nb_frame)
    current_folder = strcat(data_path, '/', video_id);
    sub_real =  strcat(current_folder, '/', "real");
    sub_fake =  strcat(current_folder, '/', "synthesis");
    nb_fake = 0;
    
    elements_in_dir = dir(fullfile(current_folder,'synthesis'));
    
    for i = 1 : numel(elements_in_dir)
        folder_fake = elements_in_dir(i);
        if folder_fake.name(1) ~= '.'
            %avoid to compute './' and '../'
            nb_fake = nb_fake + 1;
        end
    end
    clusters = cell(nb_frame, nb_fake + 1);
    
    %open real frame
    frame_dir_real = dir(fullfile(sub_real,'0*.tiff'));
    frame_dir_fake = dir(fullfile(sub_fake,'id*'));
    for i=1 : nb_frame
        name = frame_dir_real(i).name;
        clusters(i,1) = cellstr(strcat(sub_real, "/", name));
        for j=1: nb_fake
            concat_path = strcat(sub_fake, "/" ,frame_dir_fake(j).name);
            sub_frame_dir_fake =  dir(fullfile(concat_path,'0*.tiff'));
            name = sub_frame_dir_fake(j).name;
            clusters(i,j+1) = cellstr(strcat(concat_path, "/", name));
        end
        
    end
end

function ffts = get_ffts(frame_paths, nb_frame)
    ffts = cell(nb_frame , length(frame_paths));
    for num_frame=1 : nb_frame
        for path_id = 1:length(frame_paths)
            img = imread(char(frame_paths(path_id)));
%             ////
%  FFT
            F = rgb2gray(img);
            F = fft2(double(F));
            F = fftshift(F);
            
            %F = abs(F); % Get the magnitude
            %F = fftshift(log(F+1)); % Use log, for perceptual scaling, and +1 since log(0) is undefined
            %imagesc(F);
            %F = mat2gray(F); % Use mat2gray to scale the image between 0 and 1
            ffts(num_frame, path_id) = {F};
%            
%             /////
%  WAVELET                
           %[a2,h2,v2,d2] = haart2(img);
           
%            [cA,cH,cV,cD] = dwt2(img,'haar');
%            
%            [size_x, size_y, channel] = size(img);
%            [size_x_wave, size_y_wave, channel_wave] = size(cD);
%            F = ones(size_x, size_y, channel);
%            indice_x_start = ceil(1/2 * size_x) - ceil(1/2 * size_x_wave);
%            indice_x_end = indice_x_start + size_x_wave - 1;
%            
%            indice_y_start = ceil(1/2 * size_y) - ceil(1/2 * size_y_wave);
%            indice_y_end = indice_y_start + size_y_wave - 1;
%            
%            F(indice_x_start : indice_x_end,...
%                 indice_y_start : indice_y_end, :)...
%                 = cD;
%             
%             ffts(num_frame, path_id) = {F};




            
            
            %imshow(F,[]);
            %imwrite(F,'myGray.png');
            
            %img_fft = fft2(img) / 715;
            %imshow(uint8(fftshift(abs(img_fft)))) ;
            %ffts(num_frame, path_id) = fftshift(img_fft);
        end
    
    end
    
    
end

function result = format_fft(frame_paths, ffts, nb_frame)
    %[size_x, size_y, channel] = size(cell2mat(ffts(1,1)));
    %result = zeros(size_x * 2,size_y * length(frame_paths));
   
    
%     result_cell = cell(1, 3);
%     for i=1: length(frame_paths)
%         r = rgb2gray(imread(char(frame_paths(1,i))));
%         F = cell2mat(ffts(1, i));    
%         %result_cell(i) = {[r;im2uint8(F)]};
%         result_cell(i) = {[r;uint8(abs(F))]};
%         imshow([r;uint8(abs(F))]);
%     end
%     
%     result = im2uint8(cell2mat(result_cell));
%     %imshow(im2uint8(F));
%     %
%     
%     out_name = strsplit(char(frame_paths(1,1)),'/');
%     out_name = strcat(out_name, ".tiff");
%     out_path = strcat('spectre_compare/results/', char(out_name(7)));
%     %mapTrace=colormap(winter);
%     imwrite(result, winter, out_path);


%   ////////////////////////

    orig_fft = real(cell2mat(ffts(1, 1)));
    fake_fft1 = real(cell2mat(ffts(1, 2)));
    fake_fft2 = real(cell2mat(ffts(1, 3)));
    fake_fft3 = real(cell2mat(ffts(1, 4)));
    
    %imagesc(abs(orig_fft));
    
    test1 = orig_fft == fake_fft1;
    equal_test1 = min(min(test1));
    diff_test1 = orig_fft - fake_fft1;
    mean_diff_test1 = mean(mean(diff_test1));
    
    
    test2 = orig_fft == fake_fft2;
    equal_test2 = min(min(test2));
    diff_test2 = orig_fft - fake_fft2;
    mean_diff_test2 = mean(mean(diff_test2));
    
  
    test3 = orig_fft == fake_fft3;
    equal_test3 = min(min(test3));
    diff_test3 = orig_fft - fake_fft3;
    mean_diff_test3 = mean(mean(diff_test3));
    
%   ////////////////////////

%     tiledlayout(3,2);
%     for i=1: length(frame_paths)
%         
%         ax1 = nexttile;
%         r = imread(char(frame_paths(1,i)));
%         imagesc(r);
%         ax2 = nexttile;
%         F = cell2mat(ffts(1, i)); 
%         imagesc(abs(F));
%         
%     end
    
%   ////////////////////////

    subplot(2, 3, 1);
    fft_orig = 255 * mat2gray(abs(orig_fft));
    imagesc(fft_orig);
    title('orignal', 'FontSize', 20);
    subplot(2, 3, 2);
    fft_fake1 = 255 * mat2gray(abs(fake_fft1));
    imagesc(fft_fake1);
    title('fake 1', 'FontSize', 20);
    
    subplot(2, 3, 3);
    fft_fake2 = 255 * mat2gray(abs(fake_fft2));
    %imshow(fft_fake2);
    %imshow(log(fft_fake2), []);
    %imshow(uint8(abs(fft_fake2)));
    imagesc(log(fft_fake2));
    title('fake 2', 'FontSize', 20)
    
    subplot(2, 3, 4);
    fft_fake3 = 255 * mat2gray(abs(fake_fft3));
    imagesc(fft_fake3);
    title('fake 3', 'FontSize', 20);
    axis on;
    
    
    
    




%     subplot(2, 3, 3);
%     scaledFFTr = 255 * mat2gray(real(orig_fft));
%     imshow(log(scaledFFTr), []);
%     title('Log of Real Part of Spectrum', 'FontSize', 20)
%     subplot(2, 3, 3);
%     scaledFFTi = mat2gray(imag(orig_fft));
%     imshow(log(scaledFFTi), []);
%     axis on;
%     title('Log of Imaginary Part of Spectrum', 'FontSize', 20)

result = 0;
end
