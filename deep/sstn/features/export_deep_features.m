
inputSize = net.Layers(1).InputSize;

path_train_img = "/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/limited_train_val";
imdsTrain = imageDatastore(path_train_img, 'IncludeSubfolders', true,'FileExtensions','.tiff');
imdsTrain = augmentedImageDatastore(inputSize(1:2),imdsTrain);

path_test_img = "/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/limited_test";
imdsTest = imageDatastore(path_test_img, 'IncludeSubfolders', true,'FileExtensions','.tiff');
imdsTest = augmentedImageDatastore(inputSize(1:2),imdsTest);


layer = 'avg_pool';
deep_featuresTrain = activations(net,imdsTrain,layer,'OutputAs','rows');
deep_featuresTest = activations(net,imdsTest,layer,'OutputAs','rows');