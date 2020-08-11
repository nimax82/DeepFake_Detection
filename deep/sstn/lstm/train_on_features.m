clc;
clear;

%load features
deep_featuresTrain = load('../features/deep_features_train.mat').deep_featuresTrain;
deep_featuresTest = load('../features/deep_features_test.mat').deep_featuresTest;

svm_featuresTrain = load('../steg/features/training/color.mat').features_color;
svm_featuresTest = load('../steg/features/test/color_test.mat').features_color;

%concat features
featuresTrain = [deep_featuresTrain, svm_featuresTrain];
featuresTest = [deep_featuresTest, svm_featuresTest];
%save('./featuresTrain', 'featuresTrain');
%save('./featuresTest', 'featuresTest');

%load labels
training_labels = load('../steg/features/training/train_labels_frame.mat').expected_labels_frame;
test_labels = load('../steg/features/test/test_labels_frame.mat').expected_labels_frame;

%create network
concat_net = xception;
lgraph = layerGraph(concat_net);
[learnableLayer, classLayer] = findLayersToReplace(lgraph);
inputSize = concat_net.Layers(1).InputSize;

% learnable layer
newLearnableLayer = fullyConnectedLayer(2, ...
        'Name','new_fc', ...
        'WeightLearnRateFactor',10, ...
        'BiasLearnRateFactor',10);
    
lgraph = replaceLayer(lgraph,learnableLayer.Name,newLearnableLayer);

%classification layer
newClassLayer = classificationLayer('Name','new_classoutput');
lgraph = replaceLayer(lgraph,classLayer.Name,newClassLayer);

%format data
nb_train_frame = 134 * 25;
nb_frame = 179 * 25 ;
nb_valid_frame = nb_frame - nb_train_frame;

featuresValidation = featuresTrain(nb_train_frame + 1:nb_frame,:);
featuresTrain = featuresTrain(1:nb_train_frame,:);

labels_validation = training_labels(nb_train_frame + 1:nb_frame,:);
labels_train = training_labels(1:nb_train_frame,:);
cat_labels_train = categorical(labels_train);
cat_labels_validation = categorical(labels_validation);
%train_val_ds = datastore('./lstm/featuresTrain.mat', 'Type','tabulartext');


miniBatchSize = 32;
valFrequency = floor(nb_train_frame/miniBatchSize);
options = trainingOptions('adam', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',8, ...
    'InitialLearnRate',0.0002, ...
    'ValidationData',{featuresValidation,cat_labels_validation}, ...
    'Shuffle','every-epoch', ...
    'ValidationFrequency',valFrequency, ...
    'Verbose',false, ...
    'Plots','training-progress');

concat_net = trainNetwork(featuresTrain, cat_labels_train, lgraph, options);
