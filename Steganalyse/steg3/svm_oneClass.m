% --- training --- %

imported_cfa = load('features/training/cfa/cfa.mat');
train_f_cfa = imported_cfa.features_cfa;

imported_color = load('features/training/color/color.mat');
train_f_color = imported_color.features_color;

imported_labels = load('features/training/training_labels.mat');
training_labels = imported_labels.expected_labels;

labels = zeros( 50 * 20, 1);

formatted_f_color = NaN(50 * 20,18157);
formatted_f_cfa =  NaN(50 * 20,10323);

j=0;
for i=1 : numel(training_labels)
    if training_labels(i) == 0
        j = j+1;
        formatted_f_color(j, :) = train_f_color(i,:);
        formatted_f_cfa(j, :) = train_f_cfa(i,:);
    end
end

SVMModel_CFA = fitcsvm(formatted_f_cfa, labels);
SVMModel_Color = fitcsvm(formatted_f_color, labels);

% --- Classification --- %
import_test_cfa = load('features/test/cfa/cfa_test.mat');
test_f_cfa = import_test_cfa.features_cfa;

import_test_color = load('features/test/color/color_test.mat');
test_f_color = import_test_color.features_color;

[label_cfa, score_cfa] = predict(SVMModel_CFA, test_f_cfa);
[label_color, score_color] = predict(SVMModel_Color, test_f_color);
