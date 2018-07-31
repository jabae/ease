function Aem = get_Aem_scan(obj, mscan)
%% get EM masks in a specified scan 
%{
%}

%% inputs: 
%{
    scan_id: scan ID 
%}

%% outputs: 
%{
    Aem: num_slice X 1 cell array
%}

%% author: 
%{
    Pengcheng Zhou 
    Columbia University, 2018 
    zhoupc1988@gmail.com
%}

if ~exist('mslice','var')
    mscan = obj.scan_id;
end
Aem = cell(obj.num_slices, 1);

fprintf('Projecting EM segments to all slices in scan %d\n', mscan);
for mslice=1:obj.num_slices
    Aem{mslice} = obj.extract_em_segments(mscan , mslice);
end

assignin('base', 'current_scan_id_for_em', mscan); 
assignin('base', 'Aem', Aem); 
fprintf('\nDone!\n'); 
end