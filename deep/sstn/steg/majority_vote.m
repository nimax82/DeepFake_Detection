expected_labels = load('features/test/test_labels_frame.mat');
test_labels = expected_labels.expected_labels_frame;

%test_labels = expected_labels;

nb_frame = 15;
formated_labels_color = NaN(numel(label_color) / nb_frame, 1);
formated_test_labels = NaN(numel(label_color) / nb_frame, 1);
sum_color = 0;
sum_test = 0;
count_middle_color = 0;
j = 0;
for i=1 : numel(label_color)
    sum_color = sum_color + label_color(i);
    sum_test = sum_test + test_labels(i);
    if mod(i,nb_frame) == 0
        j = j + 1;
        count_middle_color = count_middle_color + double(sum_color == (nb_frame / 2));
        
        formated_labels_color(j) = int8(sum_color > (nb_frame / 2));
        
        formated_test_labels(j) = int8(sum_test > (nb_frame / 2));
        
        sum_color = 0;
        sum_test = 0;
    end
end

% compare

similarity_color = (formated_labels_color == formated_test_labels);
result_color = sum(double(similarity_color));



