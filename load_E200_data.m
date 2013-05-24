function data = load_E200_data(path_file,head,isscan,doyag,doceloss,docegain)

% load data containing filenames
load(path_file);

% if data is from a scan . . .
if isscan
    
    save_path  = [head char(scan_info(1).save_path(1))];
    all_file = dir(fullfile(save_path,'*.mat'));
    j = 1;
    for i=1:length(all_file)
        a = strfind(all_file(i).name,'filenames');
        b = strfind(all_file(i).name,'scan_info');
        % identify .mat files containg BSA data and concatanate
        if isempty(a) && isempty(b)
            d(j) = load([save_path '/' all_file(i).name]);
            j=j+1;
        end
    end
    
    % check to make sure we have expected files
    if length(d) ~= length(scan_info); error('Number of image files and .mat files not equal'); end;
 
end

% get information of scan length and number of shots
n_step = length(d);
n_shot = d(1).param.n_shot;
for i=1:n_step 
    data.epics_shots(i) = length(d(i).epics_data);
    data.energy(i) = scan_info(i).Control_PV;
end
max_epics_shots = max(data.epics_shots);

% get epics fields
E_var_name  = fieldnames(d(1).epics_data);
N_epics_var = length(E_var_name);
for i = 1:N_epics_var
    if strcmp(E_var_name{i}(end),'_')
        epics_vars{i} = E_var_name{i}(1:(end-1));
    else
        epics_vars{i} = E_var_name{i};
    end
    data.epics.(epics_vars{i}) = zeros(max_epics_shots,n_step);
end
data.epics.APID_ind = zeros(max_epics_shots,n_step);
data.epics.PY_ind = zeros(max_epics_shots,n_step);

% assign epics data
for i=1:n_step
    for j=1:data.epics_shots(i)
        for k=1:N_epics_var       
            data.epics.(epics_vars{k})(j,i) = d(i).epics_data(j).(E_var_name{k});
        end
    end
end

% sort pyro data
for i=1:n_step
    [~,data.epics.py_sort(1:data.epics_shots(i),i)] = sort(data.epics.BLEN_LI20_3014_BRAW(1:data.epics_shots(i),i));
end


% get aida info
aida_daq = d(1).param.aida_daq;
if aida_daq
    data.aida.pulse_id = zeros(n_shot,n_step);
    data.aida.EPID_ind = zeros(n_shot,n_step);
    data.aida.py_sort  = zeros(n_shot,n_step);
    data.aida.py_ind   = zeros(n_shot,n_step);
    
    N_bpms = length(d(1).aida_data(1).bpms);
    N_toro = length(d(1).aida_data(1).toro);
    N_klys = length(d(1).aida_data(1).klys);
    
    % loop over BPMs
    for i = 1:N_bpms
        name = d(1).aida_data(1).bpms(i).name;
        bpms_name{i} = regexprep(name, ':', '_');
        data.aida.(bpms_name{i}).x = zeros(n_shot,n_step);
        data.aida.(bpms_name{i}).y = zeros(n_shot,n_step);
        data.aida.(bpms_name{i}).tmit = zeros(n_shot,n_step);
        data.aida.(bpms_name{i}).stat = zeros(n_shot,n_step);
        data.aida.(bpms_name{i}).good = zeros(n_shot,n_step);
    end
    
    % loop over Toros
    for i = 1:N_toro
        name = d(1).aida_data(1).toro(i).name;
        toro_name{i} = regexprep(name, ':', '_');
        data.aida.(toro_name{i}).tmit = zeros(n_shot,n_step);
        data.aida.(toro_name{i}).stat = zeros(n_shot,n_step);
        data.aida.(toro_name{i}).good = zeros(n_shot,n_step);
    end
        
    % loop over Klys
    for i = 1:N_klys
        name = d(1).aida_data(1).toro(i).name;
        klys_name{i} = regexprep(name, ':', '_');
        data.aida.(klys_name{i}).phas = zeros(n_shot,n_step);
        data.aida.(klys_name{i}).stat = zeros(n_shot,n_step);
    end
    
    % assign aida data
    for i=1:n_step
        for j=1:n_shot
            data.aida.pulse_id(j,i) = d(i).aida_data(j).pulse_id;
            for k=1:N_bpms              
                data.aida.(bpms_name{k}).x(j,i) = d(i).aida_data(j).bpms(k).x;
                data.aida.(bpms_name{k}).y(j,i) = d(i).aida_data(j).bpms(k).y;
                data.aida.(bpms_name{k}).tmit(j,i) = d(i).aida_data(j).bpms(k).tmit;
                data.aida.(bpms_name{k}).stat(j,i) = d(i).aida_data(j).bpms(k).stat;
                data.aida.(bpms_name{k}).good(j,i) = d(i).aida_data(j).bpms(k).goodmeas;
            end
            for l=1:N_toro
                data.aida.(toro_name{l}).tmit(j,i) = d(i).aida_data(j).toro(l).tmit;
                data.aida.(toro_name{l}).stat(j,i) = d(i).aida_data(j).toro(l).stat;
                data.aida.(toro_name{l}).good(j,i) = d(i).aida_data(j).toro(l).goodmeas;
            end
            for m=1:N_klys
                data.aida.(klys_name{m}).phas(j,i) = d(i).aida_data(j).klys(m).phase;
                data.aida.(klys_name{m}).stat(j,i) = d(i).aida_data(j).klys(m).stat;
            end
        end
    end
    
    % find matched pulse id indices
    for i=1:n_step
        [~,IA,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(:,i),data.aida.pulse_id(:,i),'rows','stable');
        data.aida.EPID_ind(:,i) = IA;
        data.epics.APID_ind(IA,i) = IB;
        
        data.aida.py_sort(:,i) = data.epics.py_sort(IA,i);
        [~,~,IB] = intersect(data.epics.PATT_SYS1_1_PULSEID(data.epics.py_sort(1:data.epics_shots(i),i),i),data.aida.pulse_id(:,i),'rows','stable');        
        data.aida.py_ind(:,i) = IB;

    end
end

if doyag
    
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
    
if docegain
    
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
    
if doceloss
    
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