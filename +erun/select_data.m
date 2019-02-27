datasets = {'pinky40', 'pinky100'};
fprintf('\n**********choose the data to use**********\n');
for m=1:length(datasets)
    fprintf('%d: %s\n', m, datasets{m});
end
fprintf('********************************************\n');

data_id = input('data ID: ');
while true
    if any(data_id==[1, 2])
        data_name = datasets{data_id};
        fprintf('you selected data %s\n', data_name);
        break;
    else
        data_id = input('please type a valid data ID: ');
    end
end
