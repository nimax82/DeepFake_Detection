% --- training --- %

%imported_cfa = load('features/training/cfa.mat');
%train_f_cfa = imported_cfa.features_cfa;

imported_color = load('features/training/color.mat');
train_f_color = imported_color.features_color;

imported_labels = load('features/training/train_labels_frame.mat');
training_labels = imported_labels.expected_labels_frame;
%labels = zeros(1, size(train_f_cfa, 1));

%SVMModel_CFA = fitcsvm(train_f_cfa, training_labels);
SVMModel_Color = fitcsvm(train_f_color, training_labels);

% --- Classification --- %
%import_test_cfa = load('features/test/cfa_test.mat');
%test_f_cfa = import_test_cfa.features_cfa;

import_test_color = load('features/test/color_test.mat');
test_f_color = import_test_color.features_color;

%[label_cfa, score_cfa] = predict(SVMModel_CFA, test_f_cfa);
[label_color, score_color] = predict(SVMModel_Color, test_f_color);
