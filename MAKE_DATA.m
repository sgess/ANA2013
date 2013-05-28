clear all;
disp = '/nas/nas-li20-pm01/E200/2013/20130428/E200_10783/';
ar_qs = '/nas/nas-li20-pm01/E200/2013/20130428/E200_10794/';
pr = '/nas/nas-li20-pm01/E200/2013/20130514/E200_11159/';

save_dir = '/Users/sgess/Desktop/plots/2013/April28/';

%head = '/Volumes/PWFA 4big';
head = '/Users/sgess/Desktop/FACET/2013/DATA';

doyag = 1;
doceloss = 0;
docegain = 0;

data = load_E200_data(ar_qs,head,doyag);
 
savE = 1;
%[eta_max, eta_cent, eta_fmin, eta_fmax] = DISPANA(data,save_dir,savE);
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
