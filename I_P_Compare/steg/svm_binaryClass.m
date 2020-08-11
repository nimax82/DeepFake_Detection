clc;
clear;

% --- training --- %
imported_color = load('features/raw/train/color.mat');
train_f_color = imported_color.features_color;

imported_labels = load('/Users/test/Documents/Stage_Maxime/DNGI_Detection/I_P_Compare/labels/gen_raw_train_labels_frame.mat');
training_labels = imported_labels.expected_labels_frame;
%labels = zeros(1, size(train_f_cfa, 1));

SVMModel_Color = fitcsvm(train_f_color, training_labels);

% --- Classification --- %
import_test_color = load('features/raw/test/color.mat');
test_f_color = import_test_color.features_color;

[label_color, score_color] = predict(SVMModel_Color, test_f_color);
