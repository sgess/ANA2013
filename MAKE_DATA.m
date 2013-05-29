clear all;
%ar_qs = '/nas/nas-li20-pm01/E200/2013/20130428/E200_10794/';
%save_name = 'slim_10794.mat';

%rb_qs = '/nas/nas-li20-pm01/E200/2013/20130429/E200_10911/';
%save_name = 'slim_10911.mat';

rb_qs2 = '/nas/nas-li20-pm01/E200/2013/20130429/E200_10915/';
save_name = 'slim_10915.mat';

%save_dir = '/Users/sgess/Desktop/plots/2013/April28/';
save_dir = '/Users/sgess/Desktop/plots/2013/April29/';

save_data = '/Users/sgess/Desktop/data/2013/slims/';

head = '/Volumes/PWFA 4big';
%head = '/Users/sgess/Desktop/FACET/2013/DATA';

doyag = 1;
doceloss = 0;
docegain = 0;

data = load_E200_data(rb_qs2,head,doyag);
 
EPICS_PID = data.EPICS.PATT_SYS1_1_PULSEID;
EPICS_SCANSTEP = data.EPICS.scan_step;
EPICS_DATASET = data.EPICS.dataset;
IMAGE_PID = data.YAG.pulse_id;
IMAGE_SCANSTEP = data.YAG.scan_step;

[data.EPICS.UID, data.YAG.UID] = assign_UID(EPICS_PID,EPICS_SCANSTEP,EPICS_DATASET,IMAGE_PID,IMAGE_SCANSTEP);
[~,ii,ie] = intersect(data.YAG.UID,data.EPICS.UID);
isequal(data.EPICS.PATT_SYS1_1_PULSEID(ie),data.YAG.pulse_id(ii))
save([save_data save_name],'data');
%%
clear all;
%disp = '/nas/nas-li20-pm01/E200/2013/20130428/E200_10783/';
%save_name = 'slim_10783.mat';
disp = '/nas/nas-li20-pm01/E200/2013/20130429/E200_10888/';
save_name = 'slim_10888.mat';

head = '/Volumes/PWFA 4big';
doyag = 1;
data = load_E200_data(disp,head,doyag);
savE = 0;
%save_dir = '/Users/sgess/Desktop/plots/2013/April28/';
save_dir = '/Users/sgess/Desktop/plots/2013/April29/';

EPICS_PID = data.EPICS.PATT_SYS1_1_PULSEID;
EPICS_SCANSTEP = data.EPICS.scan_step;
EPICS_DATASET = data.EPICS.dataset;
IMAGE_PID = data.YAG.pulse_id;
IMAGE_SCANSTEP = data.YAG.scan_step;

[data.EPICS.UID, data.YAG.UID] = assign_UID(EPICS_PID,EPICS_SCANSTEP,EPICS_DATASET,IMAGE_PID,IMAGE_SCANSTEP);

[eta_max, eta_cent, eta_fmin, eta_fmax] = DISPANA(data,save_dir,savE);

save_data = '/Users/sgess/Desktop/data/2013/slims/';
save([save_data save_name],'data');


% EPICS_PID = data.EPICS.PATT_SYS1_1_PULSEID;
% EPICS_SCANSTEP = data.EPICS.scan_step;
% EPICS_DATASET = data.EPICS.dataset;
% IMAGE_PID = data.YAG.pulse_id;
% IMAGE_SCANSTEP = data.YAG.scan_step;
% AIDA_PID = data.AIDA.pulse_id;
% AIDA_SCANSTEP = data.AIDA.scan_step;
% 
% [data.EPICS.UID, data.YAG.UID, data.AIDA.UID] = assign_UID(EPICS_PID,EPICS_SCANSTEP,EPICS_DATASET,IMAGE_PID,IMAGE_SCANSTEP,AIDA_PID,AIDA_SCANSTEP);
% [~,ii,ie] = intersect(data.YAG.UID,data.EPICS.UID);
% isequal(data.EPICS.PATT_SYS1_1_PULSEID(ie),data.YAG.pulse_id(ii))
% 
% [~,ia,ie] = intersect(data.AIDA.UID,data.EPICS.UID);
% isequal(data.EPICS.PATT_SYS1_1_PULSEID(ie),data.AIDA.pulse_id(ii))


%%
pyros = [];
pyro2 = [];
specs = [];
spec2 = [];
for i = 1:length(data.epics_shots);
pyros = [pyros; data.epics.BLEN_LI20_3014_BRAW(data.YAG.EPID_ind(data.YAG.py_ind(:,i),i),i)];
pyro2 = [pyro2; data.epics.BLEN_LI20_3014_BRAW(data.YAG.EPID_ind(:,i),i)];
specs = [specs data.YAG.spectra(:,data.YAG.py_ind(:,i),i)];
spec2 = [spec2 data.YAG.spectra(:,:,i)];
end

[a,b] = sort(pyros);
[c,d] = sort(pyro2);
figure(1);
imagesc(specs(:,b));
figure(2);
imagesc(spec2);

%%
pyros = [];
yag_specs = [];
cel_specs = [];
ceg_specs = [];
for i = 1:length(data.epics_shots);
%pyros = [pyros; data.epics.BLEN_LI20_3014_BRAW(data.YAG.EPID_ind(data.YAG.py_ind(:,i),i),i)];
%yag_specs = [yag_specs data.YAG.spectra(:,data.YAG.py_ind(:,i),i)];
%cel_specs = [cel_specs data.CELOSS.spec(:,data.YAG.py_ind(:,i),i)];
%ceg_specs = [ceg_specs data.CEGAIN.spec(:,data.YAG.py_ind(:,i),i)];

end
