%expected_labels = load('features/test/test_labels.mat');
%test_labels = expected_labels.expected_labels;

test_labels = expected_labels;

nb_frame = 15;
formated_labels_color = NaN(numel(label_color) / nb_frame, 1);
formated_labels_cfa = NaN(numel(label_color) / nb_frame, 1);
formated_test_labels = NaN(numel(label_color) / nb_frame, 1);
formated_labels_binary = NaN(numel(label_color) / nb_frame, 1);
sum_color = 0;
sum_test = 0;
sum_cfa = 0;
sum_bin = 0;
count_middle_color = 0;
count_middle_cfa = 0;
j = 0;
for i=1 : numel(label_color)
    sum_color = sum_color + label_color(i);
    sum_cfa = sum_cfa + label_cfa(i);
    sum_test = sum_test + test_labels(i);
    sum_bin = sum_bin + label_color(i) + label_cfa(i);
    if mod(i,nb_frame) == 0
        j = j + 1;
        count_middle_color = count_middle_color + double(sum_color == (nb_frame / 2));
        count_middle_cfa = count_middle_cfa + double(sum_cfa == (nb_frame / 2));
        formated_labels_color(j) = int8(sum_color > (nb_frame / 2));
        formated_labels_cfa(j) = int8(sum_cfa > (nb_frame / 2));
        formated_test_labels(j) = int8(sum_test > (nb_frame / 2));
        formated_labels_binary = int8(sum_bin > nb_frame);
        sum_color = 0;
        sum_cfa = 0;
        sum_test = 0;
        sum_bin = 0;
    end
end

% compare

similarity_color = (formated_labels_color == formated_test_labels);
result_color = sum(double(similarity_color));

similarity_cfa = int8((formated_labels_cfa == formated_test_labels));
result_cfa = sum(double(similarity_cfa));


similarity_bin = (formated_labels_binary == formated_test_labels);
result_bin = sum(double(similarity_bin));




