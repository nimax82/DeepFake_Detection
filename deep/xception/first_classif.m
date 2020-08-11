clc;
clear;

%Network
net = xception;
lgraph = layerGraph(net);
[learnableLayer, classLayer] = findLayersToReplace(lgraph);
inputSize = net.Layers(1).InputSize;

% learnable layer
newLearnableLayer = fullyConnectedLayer(2, ...
        'Name','new_fc', ...
        'WeightLearnRateFactor',10, ...
        'BiasLearnRateFactor',10);
    
lgraph = replaceLayer(lgraph,learnableLayer.Name,newLearnableLayer);

%classification layer
newClassLayer = classificationLayer('Name','new_classoutput');
lgraph = replaceLayer(lgraph,classLayer.Name,newClassLayer);

%set up data for training
path_img = "/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/deep/train_val";
imds = imageDatastore(path_img, 'IncludeSubfolders', true,'FileExtensions','.tiff');
    %labels
labels = zeros(numel(imds.Files), 1);
for i=1:numel(imds.Files)
    name = imds.Files(i);
    if count(name, 'id') == 1
        labels(i) = 0;
    elseif count(name, 'id') == 2
        labels(i) = 1;
    end
end

cat_labels = categorical(labels);
imds.Labels = cat_labels;
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.75);

augimdsTrain = augmentedImageDatastore(inputSize(1:2), imdsTrain);
augimdsValidation = augmentedImageDatastore(inputSize(1:2), imdsValidation);


%training
miniBatchSize = 32;
valFrequency = floor(numel(augimdsTrain.Files)/miniBatchSize);
options = trainingOptions('adam', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',8, ...
    'InitialLearnRate',0.0002, ...
    'ValidationData',augimdsValidation, ...
    'Shuffle','every-epoch', ...
    'ValidationFrequency',valFrequency, ...
    'Verbose',false, ...
    'Plots','training-progress');


net = trainNetwork(augimdsTrain,lgraph,options);

% 
% image = imread("/Volumes/VERBATIM_HD/Stage_Maxime/Celeb-DF-v2/svm/training/id0_0001/face/0000.tiff");
% image = imresize(image,inputSize(1:2));
% classNames = net.Layers(end).ClassNames;
% 
% [label,scores] = classify(net,image);
% sc = 100*scores(classNames == label);