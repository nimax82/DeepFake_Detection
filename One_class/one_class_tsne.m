imageFileName = '/Volumes/VERBATIM_HD/Stage_Maxime/one_class/img';
labelFileName = '/Volumes/VERBATIM_HD/Stage_Maxime/one_class/labels/labels.mat';
labels = load(labelFileName).labels;

%convert images to 1d vector (Nx299x299x3 -> Nx1x268203);
image_celeb = get_img_1d(fullfile(imageFileName, 'celeb'), 'id*', 5); %celeb orig
image_ff_orig = get_img_1d(fullfile(imageFileName, 'ff'), 'gen_*', 5);%ff orig
image_ff_fake = get_img_1d(fullfile(imageFileName, 'fake'), 'gen_*', 5);%ff fake

%concatenate vectors
images = [image_celeb; image_ff_orig; image_ff_fake];
Y = tsne(images,'Algorithm','barneshut', 'Distance', 'cosine','NumPCAComponents',50);
gscatter(Y(:,1),Y(:,2),labels)


function images_1d = get_img_1d(imgFolder, sub_name, nb_frame)
    folder_container = dir(fullfile(imgFolder, sub_name));
    nb_video = numel(folder_container);
    images_1d = NaN(nb_video * nb_frame, 16384);
    for i = 1  : nb_video
        folder = folder_container(i);
        images_root = strcat(folder.folder, "/" , folder.name);
        filenames=dir(fullfile(images_root, '0*.tiff'));
        local_nb_frame = numel(filenames);
        assert(nb_frame == local_nb_frame);
        for n = 1:nb_frame
            fullname=fullfile(filenames(n).folder, filenames(n).name);
            X = imread(fullname);
            X = rgb2gray(X);
            X = imresize(X,[128 128]);
            images_1d(((i-1)*nb_frame) + n, :) = reshape(X, 1, []);
        end
    end
end