clc;
clear;

fontSize = 20;
data_path = '/Volumes/VERBATIM HD/Stage_Maxime/faceForensics/picked_for_fft/fake_youtube';
swap_folder = strcat(data_path, "/", "Faceswap");
deep_folder = strcat(data_path, "/", "Deepfake");
f2f_folder = strcat(data_path, "/", "Face2Face");
nText_folder = strcat(data_path, "/", "NeuralTexture");

current_folder = nText_folder;

images = dir(fullfile(current_folder, '0*.png'));



image_path = current_folder + "/" + images(1).name;
img1 = imread(image_path);
fft1s = fftshift(fft2(double(img1)));
fft1s = 255 * mat2gray(abs(fft1s));

image_path = current_folder + "/" + images(2).name;
img2 = imread(image_path);
fft2s = fftshift(fft2(double(img2)));
fft2s = 255 * mat2gray(abs(fft2s));

image_path = current_folder + "/" + images(3).name;
img3 = imread(image_path);
fft3s = fftshift(fft2(double(img3)));
fft3s = 255 * mat2gray(abs(fft3s));

 %Display images
figure;
subplot(2, 3, 1);
imshow(img1)
axis on;
title('Image 1', 'FontSize', fontSize)
subplot(2, 3, 2);
imshow(img2);
axis on;
title('Image 2', 'FontSize', fontSize)
subplot(2, 3, 3);
imshow(img3);
axis on;
title('Image 3', 'FontSize', fontSize)

%Display fft
subplot(2, 3, 4);
imagesc(fft1s);
axis on;
title('fft 1', 'FontSize', fontSize)
subplot(2, 3, 5);
imagesc(fft2s);
axis on;
title('fft 2', 'FontSize', fontSize)
subplot(2, 3, 6);
imagesc(fft3s);
axis on;
title('fft 3', 'FontSize', fontSize)

