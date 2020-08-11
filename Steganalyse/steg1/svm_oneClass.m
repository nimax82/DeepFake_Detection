% --- training --- %

% imported_cfa = load('features/f_train_cfa.mat');
% t_features_cfa = imported_cfa.train_f_cfa;
% 
% imported_color = load('features/f_train_color.mat');
% t_features_color = imported_color.train_f_color;
% 
% labels = zeros(1, size(t_features_cfa, 1));
% 
% SVMModel_CFA = fitcsvm(t_features_cfa,labels);
% SVMModel_Color = fitcsvm(t_features_color, labels);

% --- Classification --- %
import_test_cfa = load('features/f_test_cfa.mat');
test_features_cfa = import_test_cfa.f_test_cfa;

import_test_color = load('features/f_test_color.mat');
test_features_color = import_test_color.f_test_color;

[label_cfa, score_cfa] = predict(SVMModel_CFA, test_features_cfa);
[label_color, score_color] = predict(SVMModel_Color, test_features_color);
