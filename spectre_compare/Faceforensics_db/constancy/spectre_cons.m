clc;
clear;

data_path = '/Volumes/VERBATIM HD/Stage_Maxime/faceForensics/picked_for_fft/constancy';
orig_folder = strcat(data_path, "/", "orig");
fake_folder = strcat(data_path, "/", "fake");


images_orig = dir(fullfile(orig_folder, '0*.png'));
images_fake = dir(fullfile(fake_folder, '0*.tiff'));

mean_orig = get_result(orig_folder, images_orig);
mean_fake = get_result(fake_folder, images_fake);


function mean_s = get_result(images_folder, images)
    fontSize = 20;
    
    image_path = images_folder + "/" + images(1).name;
    img1 = imread(image_path);
    fft1s = fftshift(fft2(double(img1)));
    fft1s = 255 * mat2gray(abs(fft1s));

    image_path = images_folder + "/" + images(2).name;
    img2 = imread(image_path);
    fft2s = fftshift(fft2(double(img2)));
    fft2s = 255 * mat2gray(abs(fft2s));

    image_path = images_folder + "/" + images(3).name;
    img3 = imread(image_path);
    fft3s = fftshift(fft2(double(img3)));
    fft3s = 255 * mat2gray(abs(fft3s));

    image_path = images_folder + "/" + images(4).name;
    img4 = imread(image_path);
    fft4s = fftshift(fft2(double(img4)));
    fft4s = 255 * mat2gray(abs(fft4s));

     %Display images
    figure;
    subplot(2, 4, 1);
    imshow(img1)
    axis on;
    title('Image 1', 'FontSize', fontSize)
    subplot(2, 4, 2);
    imshow(img2);
    axis on;
    title('Image 2', 'FontSize', fontSize)
    subplot(2, 4, 3);
    imshow(img3);
    axis on;
    title('Image 3', 'FontSize', fontSize)

    subplot(2, 4, 4);
    imshow(img4);
    axis on;
    title('Image 4', 'FontSize', fontSize)

    %Display fft

    subplot(2, 4, 5);
    imagesc(fft1s);
    axis on;
    title('fft 1', 'FontSize', fontSize)

    subplot(2, 4, 6);
    imagesc(fft2s);
    axis on;
    title('fft 2', 'FontSize', fontSize)

    subplot(2, 4, 7);
    imagesc(fft3s);
    axis on;
    title('fft 3', 'FontSize', fontSize)

    subplot(2, 4, 8);
    imagesc(fft4s);
    axis on;
    title('fft 4', 'FontSize', fontSize)

    %diff mean
    diff1 = abs(fft2s - fft1s);
    diff2 = abs(fft3s - fft2s); 
    diff3 = abs(fft4s - fft3s);
    sum = diff1 + diff2 + diff3;
    sum = sum ./ 3;
    
    figure;
    imagesc(sum);
    
    mean_s = mean(mean(mean(sum)));

end
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
