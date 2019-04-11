function [Aem, segment_ids] = load_Aem(obj, mscan)
%% load Aem given the scan ID 

if ~exist('mscan','var')
    mscan = obj.scan_id;
end

%% determine number of batches 
dir_em = fullfile(obj.output_folder, sprintf('segmentation_%d_zblur_%d', ...
    obj.em_segmentation, obj.em_zblur));

if ~exist(dir_em, 'dir')
    % construct Aem for each scan
    obj.construct_Aem();
end

fprintf('loading EM footprints for scan %d.\n', obj.scan_id); 
% load segmentation IDs
segment_ids = cell(20, 1);
Aem = cell(1, 20); 
mbatch = 0;
while true
    tmp_file = fullfile(dir_em, sprintf('batch_%d_ids.mat', mbatch+1));
    if exist(tmp_file, 'file')
        % load segment_ids 
        temp = load(tmp_file);
        segment_ids{mbatch+1} = temp.segment_id;
        
        % load Aem 
        tmp_file = fullfile(dir_em, sprintf('batch_%d_scan_%d_Aem.mat', mbatch+1, mscan));
        temp = load(tmp_file); 
        Aem{mbatch+1} = temp.Aem; 
        
        mbatch = mbatch + 1; 
    else
        break; 
    end
end
segment_ids = segment_ids(1:mbatch); 
Aem = Aem(1:mbatch); 

assignin('base', 'current_scan_id_for_em', mscan); 
assignin('base', 'Aem', Aem); 
fprintf('\nDone!\n'); 