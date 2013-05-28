function EPICS = EXTRACT_EPICS(d,n_step,scan_val,scan_pv)

% loop over data structure and determine number of epics shots
epics_shots = zeros(n_step,1);
for i=1:n_step 
    epics_shots(i) = length(d(i).epics_data);
end
tot_epics_shots = sum(epics_shots);

% get epics fields
E_var_name  = fieldnames(d(1).epics_data);
N_epics_var = length(E_var_name);
for i = 1:N_epics_var
    if strcmp(E_var_name{i}(end),'_')
        epics_vars{i} = E_var_name{i}(1:(end-1));
    else
        epics_vars{i} = E_var_name{i};
    end
    EPICS.(epics_vars{i}) = zeros(tot_epics_shots,1);
end

EPICS.scan_val = zeros(tot_epics_shots,1);
EPICS.scan_pv  = cell(tot_epics_shots,1);
shots = [0; cumsum(epics_shots)];
% assign epics data
for i=1:n_step
    
    start_ind = shots(i);
    end_ind = shots(i+1);
    EPICS.scan_val((start_ind+1):end_ind) = scan_val(i);
    EPICS.scan_pv((start_ind+1):end_ind) = scan_pv(i);
    
    for j=1:epics_shots(i)
        for k=1:N_epics_var       
            EPICS.(epics_vars{k})(start_ind+j) = d(i).epics_data(j).(E_var_name{k});
        end
    end
end