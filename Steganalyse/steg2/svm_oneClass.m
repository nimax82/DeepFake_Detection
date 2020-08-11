% --- training --- %

 imported_cfa = load('features/training/cfa/cfa.mat');
 train_f_cfa = imported_cfa.features_cfa;

imported_color = load('features/training/color/color.mat');
train_f_color = imported_color.features_color;

labels = zeros(1, size(train_f_cfa, 1));

SVMModel_CFA = fitcsvm(train_f_cfa,labels);
SVMModel_Color = fitcsvm(train_f_color, labels);

% --- Classification --- %
import_test_cfa = load('features/test/cfa/cfa_test.mat');
test_f_cfa = import_test_cfa.features_cfa;

import_test_color = load('features/test/color/color_test.mat');
test_f_color = import_test_color.features_color;

[label_cfa, score_cfa] = predict(SVMModel_CFA, test_f_cfa);
[label_color, score_color] = predict(SVMModel_Color, test_f_color);
