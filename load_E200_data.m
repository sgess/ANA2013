function data = load_E200_data(path,head,doYAG,doCELOSS,doCEGAIN)

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

% get information of scan length and number of shots
n_step = length(d);
n_shot = d(1).param.n_shot;

% extract EPICS data
data.EPICS = EXTRACT_EPICS(d,n_step,DAQ_info,isscan);

% extract AIDA data
if d(1).param.aida_daq
    data.AIDA = EXTRACT_AIDA(d,n_step,n_shot,DAQ_info,isscan);
end

    % find matched pulse id indices
%     for i=1:n_step
%         [~,IA,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(:,i),data.aida.pulse_id(:,i),'rows','stable');
%         data.aida.EPID_ind(:,i) = IA;
%         data.epics.APID_ind(IA,i) = IB;
%         
%         data.aida.py_sort(:,i) = data.epics.py_sort(IA,i);
%         [~,~,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(data.epics.py_sort(1:data.epics_shots(i),i),i),data.aida.pulse_id(:,i),'rows','stable');        
%         data.aida.py_ind(:,i) = IB;
% 
%     end

if doYAG
    
    data.YAG = EXTRACT_YAG(path_file,head,isscan);
    
    data.YAG.EPID_ind = zeros(n_shot,n_step);
    data.YAG.py_sort  = zeros(n_shot,n_step);
    data.YAG.py_ind   = zeros(n_shot,n_step);
    data.epics.IPID_ind = zeros(max_epics_shots,n_step);
    if aida_daq
        data.YAG.APID_ind = zeros(n_shot,n_step);
        data.aida.IPID_ind = zeros(n_shot,n_step);
    end
    
    for i=1:n_step
        
        [~,IA,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(:,i),data.YAG.pulse_id(:,i),'rows','stable');
        data.YAG.EPID_ind(IB,i) = IA;
        data.epics.IPID_ind(IA,i) = IB;
        data.YAG.py_sort(:,i) = data.epics.py_sort(IA,i);
        [~,~,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(data.epics.py_sort(1:data.epics_shots(i),i),i),data.YAG.pulse_id(:,i),'rows','stable');        
        data.YAG.py_ind(:,i) = IB;
        if aida_daq
            [~,IA,IB] = intersect(data.aida.pulse_id(:,i),data.YAG.pulse_id(:,i),'rows','stable');
            data.YAG.APID_ind(IB,i) = IA;
            data.aida.IPID_ind(IA,i) = IB;
        end
    end
end
    
if doCEGAIN
    
    data.CEGAIN = EXTRACT_CEGAIN(path_file,head,isscan);
    data.CEGAIN.EPID_ind = zeros(n_shot,n_step);
    data.CEGAIN.py_sort  = zeros(n_shot,n_step);
    data.CEGAIN.py_ind   = zeros(n_shot,n_step);
    data.epics.IPID_ind = zeros(max_epics_shots,n_step);
    if aida_daq
        data.CEGAIN.APID_ind = zeros(n_shot,n_step);
        data.aida.IPID_ind = zeros(n_shot,n_step);
    end
    
    for i=1:n_step
        
        [~,IA,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(:,i),data.CEGAIN.pulse_id(:,i),'rows','stable');
        data.CEGAIN.EPID_ind(IB,i) = IA;
        data.epics.IPID_ind(IA,i) = IB;
        data.CEGAIN.py_sort(:,i) = data.epics.py_sort(IA,i);
        [~,~,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(data.epics.py_sort(1:data.epics_shots(i),i),i),data.CEGAIN.pulse_id(:,i),'rows','stable');        
        data.CEGAIN.py_ind(:,i) = IB;
        if aida_daq
            [~,IA,IB] = intersect(data.aida.pulse_id(:,i),data.CEGAIN.pulse_id(:,i),'rows','stable');
            data.CEGAIN.APID_ind(IB,i) = IA;
            data.aida.IPID_ind(IA,i) = IB;
        end
    end
    
end
    
if doCELOSS
    
    data.CELOSS = EXTRACT_CELOSS(path_file,head,isscan);
    data.CELOSS.EPID_ind = zeros(n_shot,n_step);
    data.CELOSS.py_sort  = zeros(n_shot,n_step);
    data.CELOSS.py_ind   = zeros(n_shot,n_step);
    data.epics.IPID_ind = zeros(max_epics_shots,n_step);
    if aida_daq
        data.CELOSS.APID_ind = zeros(n_shot,n_step);
        data.aida.IPID_ind = zeros(n_shot,n_step);
    end
    
    for i=1:n_step
        
        [~,IA,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(:,i),data.CELOSS.pulse_id(:,i),'rows','stable');
        data.CELOSS.EPID_ind(IB,i) = IA;
        data.epics.IPID_ind(IA,i) = IB;
        data.CELOSS.py_sort(:,i) = data.epics.py_sort(IA,i);
        [~,~,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(data.epics.py_sort(1:data.epics_shots(i),i),i),data.CELOSS.pulse_id(:,i),'rows','stable');        
        data.CELOSS.py_ind(:,i) = IB;
        if aida_daq
            [~,IA,IB] = intersect(data.aida.pulse_id(:,i),data.CELOSS.pulse_id(:,i),'rows','stable');
            data.CELOSS.APID_ind(IB,i) = IA;
            data.aida.IPID_ind(IA,i) = IB;
        end
    end
    
end