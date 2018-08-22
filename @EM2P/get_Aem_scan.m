function Aem = get_Aem_scan(obj, mscan, ssub)
%% get EM masks in a specified scan 
%{
%}

%% inputs: 
%{
    scan_id: scan ID 
    ssub: int, downsampling factor 
%}

%% outputs: 
%{
    Aem: d * K matrix, matrix form of Aem. here d can be downsampled by a factor of ssub 
%}

%% author: 
%{
    Pengcheng Zhou 
    Columbia University, 2018 
    zhoupc1988@gmail.com
%}

if ~exist('mscan','var')
    mscan = obj.scan_id;
end
if ~exist('ssub', 'var') || isempty(ssub)
    ssub = (diff(obj.FOV_stack(1:2))+1) / (diff(obj.FOV(1:2))+1) ; 
end
Aem = cell(obj.num_slices, 1);

fprintf('Projecting EM segments to all slices in the selected scan %d\n', mscan);
for mslice=1:obj.num_slices
    Aem{mslice} = obj.extract_em_segments(mscan , mslice, ssub);
end

% convert cell to mat
Aem = sparse(cell2mat(Aem));

assignin('base', 'current_scan_id_for_em', mscan); 
assignin('base', 'Aem', Aem); 
fprintf('\nDone!\n'); 
end