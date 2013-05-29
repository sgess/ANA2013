function data = load_E200_data(path,head,doYAG)

% locate data
data_path = [head path];

% comb files to detemine if scan or simple DAQ
isscan = 0;
all_file = dir(fullfile(data_path,'*.mat'));
scan = strfind([all_file(:).name],'scan_info');
daq  = strfind([all_file(:).name],'filenames');
if ~isempty(scan)
    isscan = 1;
elseif isempty(scan) && isempty(daq)
    error(['No DAQ information found in ' path]);    
elseif length(daq) ~= 1
    error('Single DAQ should only have one "filenames" file');
end

j = 0;
% organize .mat files
for i = 1:length(all_file)
    a = strfind(all_file(i).name,'scan_info');
    b = strfind(all_file(i).name,'filenames');
    if ~isempty(a) && isscan
        DAQ_info = load([data_path all_file(i).name]);
    elseif ~isempty(b) && ~isscan
        DAQ_info = load([data_path all_file(i).name]);
    elseif isempty(a) && isempty(b)
        j=j+1;
        d(j) = load([data_path all_file(i).name]); 
    end
end

% check to make sure we have expected files
if isscan && length(d) ~= length(DAQ_info.scan_info) 
    error('Number of image files and .mat files not equal');
elseif ~isscan && length(d) ~= 1
    error('Number of image files and .mat files not equal');
end

% get information of scan length, number of shots, and camera backgrounds
n_step = length(d);
n_shot = d(1).param.n_shot;
backgrounds = d(1).cam_back;

% get camera paths and scan info
scan_val = zeros(n_step,1);
scan_pv  = cell(n_step,1);
cams = d(1).param.cams(:,1);
for i = 1:length(cams)
    for j = 1:n_step
        if isscan
            cam_files.(cams{i}){j} = [head DAQ_info.scan_info(j).(cams{i})];
            scan_val(j) = DAQ_info.scan_info(j).Control_PV;
            scan_pv(j) = {DAQ_info.scan_info(j).Control_PV_name};
        else
            cam_files.(cams{i}){j} = [head DAQ_info.filenames.(cams{i})];
            scan_val(j) = 0;
            scan_pv = '';
        end
    end
end

% Get dataset ID
[~,b] = strtok(d(1).param.save_name,'_');
[a,~] = strtok(b,'_');
dataset_ID = str2num(a);

% extract EPICS data
data.EPICS = EXTRACT_EPICS(d,n_step,scan_val,scan_pv,dataset_ID);

% extract AIDA data
if d(1).param.aida_daq
    data.AIDA = EXTRACT_AIDA(d,n_step,n_shot,scan_val,scan_pv,dataset_ID);
end

% extract YAG data
if doYAG
    data.YAG = EXTRACT_YAG(cam_files,backgrounds,n_step,n_shot,scan_val,scan_pv,dataset_ID);
end

% match pulseid
i_lo = data.EPICS.PATT_SYS1_1_PULSEID(1:(end-1));
i_hi = data.EPICS.PATT_SYS1_1_PULSEID(2:end);
k = find(i_hi < i_lo,1,'first');
pid_lo = data.EPICS.PATT_SYS1_1_PULSEID(1:k);
pid_hi = data.EPICS.PATT_SYS1_1_PULSEID((k+1):end);
y_lo = data.YAG.pulse_id(1:(end-1));
y_hi = data.YAG.pulse_id(2:end);
j = find(y_hi < y_lo,1,'first');
yid_lo = data.YAG.pulse_id(1:j);
yid_hi = data.YAG.pulse_id((j+1):end);
[~,~,ib_lo] = intersect(yid_lo,pid_lo);
[~,~,ib_hi] = intersect(yid_hi,pid_hi);
%data.YAG.epics_index = [ib_lo; ib_hi+k];