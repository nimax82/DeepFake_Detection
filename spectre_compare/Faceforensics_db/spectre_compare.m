clc;
clear;

fontSize = 20;
data_path = '/Volumes/VERBATIM HD/Stage_Maxime/faceForensics/picked_for_fft';
orig_folder = strcat(data_path, "/", "orig");
fake_folder = strcat(data_path, "/", "fake");

image_name = cell(5);

id_image = ['a', 'b', 'c', 'd', 'e'];


for i=1 : 5
    image_name(i) = cellstr(id_image(i) + ".tiff");  
end

for i=1 : 5
    path_orig = orig_folder + "/" + char(image_name(i));
    img_orig = imread(path_orig);
    fft_orig = fftshift(fft2(double(img_orig)));
    path_fake = fake_folder + "/" + char(image_name(i));
    img_fake = imread(path_fake);
    fft_fake = fftshift(fft2(double(img_fake)));
    
    gray_orig = rgb2gray(img_orig);
    gray_fake = rgb2gray(img_fake);
    
    fft_orig = 255 * mat2gray(abs(fft_orig));
    fft_fake = 255 * mat2gray(abs(fft_fake));
    
    %fft_orig = abs(fft_orig);
    %fft_fake = abs(fft_fake);
    
    %Display images
    figure;
    subplot(2, 2, 1);
    imshow(img_orig)
    axis on;
    title('Original Image', 'FontSize', fontSize)
    subplot(2, 2, 2);
    imshow(img_fake);
    axis on;
    title('Fake Image', 'FontSize', fontSize)
    
    %Display fft
    subplot(2, 2, 3);
    %imshow(log(fft_orig), []);
    imagesc(fft_orig);
    axis on;
    title('Original fft', 'FontSize', fontSize)
    subplot(2, 2, 4);
    %imshow(log(fft_fake), []);
    %fft_fake(:,:,3) = fft_fake(:,:,3) * 10;
    %fft_fake(:,:,2) = fft_fake(:,:,2) * 10;
    %fft_fake(:,:,1) = fft_fake(:,:,1) * 10;
    imagesc(fft_fake); %colormap(winter);
    %imshow(log(abs(fft_fake) + 1),[]);
    axis on;
    title('Fake fft', 'FontSize', fontSize) 
    
    
    figure;
    w = size(fft_fake,1);
    h = size(fft_fake,2);
    fft_orig = imresize(fft_orig, [w, h]);
    imagesc(2*abs(fft_orig - fft_fake));
end


