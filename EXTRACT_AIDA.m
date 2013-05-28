function AIDA = EXTRACT_AIDA(d,n_step,n_shot,DAQ_info,isscan)

AIDA.pulse_id = zeros(n_shot,n_step);
AIDA.EPID_ind = zeros(n_shot,n_step);
AIDA.py_sort  = zeros(n_shot,n_step);
AIDA.py_ind   = zeros(n_shot,n_step);

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

