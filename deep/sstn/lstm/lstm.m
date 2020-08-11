clc;
clear; 

%load features
featuresTrain = load('./featuresTrain.mat').featuresTrain;
featuresTest = load('./featuresTest.mat').featuresTest;

%load labels
training_labels = load('../steg/features/training/train_labels_frame.mat').expected_labels_frame;
test_labels = load('../steg/features/test/test_labels_frame.mat').expected_labels_frame;


inputSize = 5;
sequenceNum = size(training_labels,1) / inputSize;
numHiddenUnits = 256;
numClasses = 2;
miniBatchSize = sequenceNum / 10;


%prepare data
%featuresTrain = uint8(featuresTrain);
%sequences = cell(sequenceNum,1);

sequences = mat2cell(featuresTrain, [size(training_labels,1), 1]);

% for i=1 : size(training_labels,1)
%     currentSeq = cell(inputSize,1);
%     for j=1 : inputSize
%         currentSeq(j) = featuresTrain((i-1) * inputSize + j, :);
%     end
%     sequences(i) = currentSeq;
% end

%define lstm
layers = [ ...
    sequenceInputLayer(inputSize)
    %averagePooling2dLayer(2, 'Name', 'reduceFusion')
    bilstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'GradientThreshold',1, ...
    'MaxEpochs',8, ...
    'MiniBatchSize',miniBatchSize, ...
    'SequenceLength','longest', ...
    'Shuffle','never', ...
    'Verbose',0, ...
    'Plots','training-progress');

net = trainNetwork(featuresTrain,training_labels,layers,options);