clc;
clear;

data_path = '/Volumes/VERBATIM HD/Stage_Maxime/faceForensics/picked_for_fft/channel';
orig_folder = strcat(data_path, "/", "orig");
fake_folder = strcat(data_path, "/", "fake");


images_orig = dir(fullfile(orig_folder, 'a.tiff'));
images_fake = dir(fullfile(fake_folder, 'a.tiff'));

mean_orig = get_result(orig_folder, images_orig);
mean_fake = get_result(fake_folder, images_fake);


function mean_s = get_result(images_folder, images)
    
    image_path = images_folder + "/" + images(1).name;
    img1 = imread(image_path);
    YCBCR = rgb2ycbcr(img1);
    
    process_img(img1, 'r', 'g', 'b');
    %process_img(YCBCR, 'Y', 'CB', 'CR');
    mean_s = 0;
end    

function process_img(img, c1, c2, c3)

    fontSize = 20;
    low_pass_Kernel = ones(3,3) ./ 9;
    high_pass_Kernel = [1,0,-1;2,0,-2;1,0,-1];
   
    %red
    img_r = img(:,:,1);
    
    %high_pass
    img_r_low = uint8(conv2(img_r, low_pass_Kernel, 'same'));
    img_r_high = uint8(conv2(img_r, high_pass_Kernel, 'same'));
    %low_pass
    fft_r_low = log(abs(fftshift(fft2(img_r_low))));
    fft_r_high = log(abs(fftshift(fft2(img_r_high))));
    
   
    %green
    img_g = img(:,:,2);
    %high_pass
    img_g_low = uint8(conv2(img_g, low_pass_Kernel, 'same'));
    img_g_high = uint8(conv2(img_g, high_pass_Kernel, 'same'));
    %low_pass
    fft_g_low = log(abs(fftshift(fft2(img_g_low))));
    fft_g_high = log(abs(fftshift(fft2(img_g_high))));
    %
    %fft_g = log(abs(fftshift(fft2(img_g))));
    
    
    %blue
    img_b = img(:,:,3);
    %high_pass
    img_b_low = uint8(conv2(img_r, low_pass_Kernel, 'same'));
    img_b_high = uint8(conv2(img_r, high_pass_Kernel, 'same'));
    %low_pass
    fft_b_low = log(abs(fftshift(fft2(img_b_low))));
    fft_b_high = log(abs(fftshift(fft2(img_b_high))));
    %
    %fft_b = log(abs(fftshift(fft2(img_b))));
    
    
    
    
%    Display
     figure;
     subplot(3, 3, 1);
     imshow(img_r_low);
     %imagesc(fft_r);
     axis on;
     title(strcat(c1 , ' low img'), 'FontSize', fontSize)
     
     subplot(3, 3, 2);
     imshow(img_g_low);
     %imagesc(fft_r);
     axis on;
     title(strcat(c2 , ' low img'), 'FontSize', fontSize)
     
     subplot(3, 3, 3);
     imshow(img_b_low);
     %imagesc(fft_r);
     axis on;
     title(strcat(c3 , ' low img'), 'FontSize', fontSize)
     
     subplot(3, 3, 4);
     imshow(fft_r_low, [0 10]);
     axis on;
     title(strcat(c1 , ' low fft'), 'FontSize', fontSize)
     
     subplot(3, 3, 5);
     imshow(fft_g_low, [0 10]);
     axis on;
     title(strcat(c2 , ' low fft'), 'FontSize', fontSize)
     
     subplot(3, 3, 6);
     imshow(fft_b_low, [0 10]);
     axis on;
     title(strcat(c3 , ' low fft'), 'FontSize', fontSize)
     
     %display profile
    subplot(3, 3, 7);
    radialProfile = get_profile(img_r_low);
    plot(radialProfile, 'b-', 'LineWidth', 2);
    grid on;
    title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
    
    subplot(3, 3, 8);
    radialProfile = get_profile(img_g_low);
    plot(radialProfile, 'b-', 'LineWidth', 2);
    grid on;
    title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
    
    subplot(3, 3, 9);
    radialProfile = get_profile(img_b_low);
    plot(radialProfile, 'b-', 'LineWidth', 2);
    grid on;
    title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
     
   %   Display high
     figure;
     subplot(3, 3, 1);
     imshow(img_r_high);
     %imagesc(fft_r);
     axis on;
     title(strcat(c1 , ' high img'), 'FontSize', fontSize)
     
     subplot(3, 3, 2);
     imshow(img_g_high);
     %imagesc(fft_r);
     axis on;
     title(strcat(c2 , ' low img'), 'FontSize', fontSize)
     
     subplot(3, 3, 3);
     imshow(img_b_high);
     %imagesc(fft_r);
     axis on;
     title(strcat(c3 , ' high img'), 'FontSize', fontSize)
     
     subplot(3, 3, 4);
     imshow(fft_r_high, [0 10]);
     axis on;
     title(strcat(c1 , ' high fft'), 'FontSize', fontSize)
     
     subplot(3, 3, 5);
     imshow(fft_g_high, [0 10]);
     axis on;
     title(strcat(c2 , ' high fft'), 'FontSize', fontSize)
     
     subplot(3, 3, 6);
     imshow(fft_b_high, [0 10]);
     axis on;
     title(strcat(c3 , ' high fft'), 'FontSize', fontSize)
    
     
     %%%
    %display profile
    subplot(3, 3, 7);
    radialProfile = get_profile(img_r_high);
    plot(radialProfile, 'b-', 'LineWidth', 2);
    grid on;
    title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
    
    subplot(3, 3, 8);
    radialProfile = get_profile(img_g_high);
    plot(radialProfile, 'b-', 'LineWidth', 2);
    grid on;
    title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
    
    subplot(3, 3, 9);
    radialProfile = get_profile(img_b_high);
    plot(radialProfile, 'b-', 'LineWidth', 2);
    grid on;
    title('Average Radial Profile of Spectrum', 'FontSize', fontSize)
     
end 

function profile = get_profile(img)
 % Get the average radial profile
    grayImage = img; 
    [rows, columns, numberOfColorChannels] = size(grayImage);
    fftOriginal = fft2(double(grayImage));
    shiftedFFT = fftshift(fftOriginal);
    shiftedFFTMagnitude = abs(shiftedFFT);
    midRow = rows/2+1;
    midCol = columns/2+1;
    maxRadius = ceil(sqrt(rows^2 + columns^2));
    profile = zeros(maxRadius, 1);
    count = zeros(maxRadius, 1);
    for col = 1 : columns
      for row = 1 : rows
        radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2);
        thisIndex = ceil(radius) + 1;
        profile(thisIndex) = profile(thisIndex) + shiftedFFTMagnitude(row, col);
        count(thisIndex) = count(thisIndex) + 1;
      end
    end
    % Get average
    profile = profile ./ count;

end
     
     
     
% % % % % % % %      
%  
%      subplot(2, 3, 1);
%      imshow(fft_r, [-1 1]);
%      %imagesc(fft_r);
%      axis on;
%      title('red', 'FontSize', fontSize)
%      
%      subplot(1, 3, 2);
%      imshow(fft_g, [0 10]);
%      %imagesc(fft_g);
%      axis on;
%      title('green', 'FontSize', fontSize)
%      
%      subplot(1, 3, 3);
%      imshow(fft_b, [0 10]);
%      %imagesc(fft_b);
%      axis on;
%      title('blue', 'FontSize', fontSize)
% % % % % % %     
      
%     image_path = images_folder + "/" + images(2).name;
%     img2 = imread(image_path);
%     fft2s = fftshift(fft2(double(img2)));
%     fft2s = 255 * mat2gray(abs(fft2s));
% 
%     image_path = images_folder + "/" + images(3).name;
%     img3 = imread(image_path);
%     fft3s = fftshift(fft2(double(img3)));
%     fft3s = 255 * mat2gray(abs(fft3s));
% 
%     image_path = images_folder + "/" + images(4).name;
%     img4 = imread(image_path);
%     fft4s = fftshift(fft2(double(img4)));
%     fft4s = 255 * mat2gray(abs(fft4s));

%      %Display images
%     figure;
%     subplot(2, 4, 1);
%     imshow(img1)
%     axis on;
%     title('Image 1', 'FontSize', fontSize)
%     subplot(2, 4, 2);
%     imshow(img2);
%     axis on;
%     title('Image 2', 'FontSize', fontSize)
%     subplot(2, 4, 3);
%     imshow(img3);
%     axis on;
%     title('Image 3', 'FontSize', fontSize)
% 
%     subplot(2, 4, 4);
%     imshow(img4);
%     axis on;
%     title('Image 4', 'FontSize', fontSize)
% 
%     %Display fft
% 
%     subplot(2, 4, 5);
%     imagesc(fft1s);
%     axis on;
%     title('fft 1', 'FontSize', fontSize)
% 
%     subplot(2, 4, 6);
%     imagesc(fft2s);
%     axis on;
%     title('fft 2', 'FontSize', fontSize)
% 
%     subplot(2, 4, 7);
%     imagesc(fft3s);
%     axis on;
%     title('fft 3', 'FontSize', fontSize)
% 
%     subplot(2, 4, 8);
%     imagesc(fft4s);
%     axis on;
%     title('fft 4', 'FontSize', fontSize)
% 
%     %diff mean
%     diff1 = abs(fft2s - fft1s);
%     diff2 = abs(fft3s - fft2s); 
%     diff3 = abs(fft4s - fft3s);
%     sum = diff1 + diff2 + diff3;
%     sum = sum ./ 3;
%     
%     figure;
%     imagesc(sum);
%     
%     mean_s = mean(mean(mean(sum)));
%%%%%
% 
% image_path = orig_folder + "/" + images_orig(1).name;
% img1 = imread(image_path);
% fft1s = fftshift(fft2(double(img1)));
% fft1s = 255 * mat2gray(abs(fft1s));
% 
% image_path = orig_folder + "/" + images_orig(2).name;
% img2 = imread(image_path);
% fft2s = fftshift(fft2(double(img2)));
% fft2s = 255 * mat2gray(abs(fft2s));
% 
% image_path = orig_folder + "/" + images_orig(3).name;
% img3 = imread(image_path);
% fft3s = fftshift(fft2(double(img3)));
% fft3s = 255 * mat2gray(abs(fft3s));
% 
% image_path = orig_folder + "/" + images_orig(4).name;
% img4 = imread(image_path);
% fft4s = fftshift(fft2(double(img4)));
% fft4s = 255 * mat2gray(abs(fft4s));
% 
%  %Display images
% figure;
% subplot(2, 4, 1);
% imshow(img1)
% axis on;
% title('Image 1', 'FontSize', fontSize)
% subplot(2, 4, 2);
% imshow(img2);
% axis on;
% title('Image 2', 'FontSize', fontSize)
% subplot(2, 4, 3);
% imshow(img3);
% axis on;
% title('Image 3', 'FontSize', fontSize)
% 
% subplot(2, 4, 4);
% imshow(img4);
% axis on;
% title('Image 4', 'FontSize', fontSize)
% 
% %Display fft
% 
% subplot(2, 4, 5);
% imagesc(fft1s);
% axis on;
% title('fft 1', 'FontSize', fontSize)
% 
% subplot(2, 4, 6);
% imagesc(fft2s);
% axis on;
% title('fft 2', 'FontSize', fontSize)
% 
% subplot(2, 4, 7);
% imagesc(fft3s);
% axis on;
% title('fft 3', 'FontSize', fontSize)
% 
% subplot(2, 4, 8);
% imagesc(fft4s);
% axis on;
% title('fft 4', 'FontSize', fontSize)
% 
