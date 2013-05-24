disp = '/Volumes/PWFA 4big/nas/nas-li20-pm01/E200/2013/20130428/E200_10783/E200_10783_scan_info.mat';
ar_qs = '/Volumes/PWFA 4big/nas/nas-li20-pm01/E200/2013/20130428/E200_10794/E200_10794_scan_info.mat';
pr = '/Volumes/PWFA 4big/nas/nas-li20-pm01/E200/2013/20130514/E200_11159/E200_11159_scan_info.mat';

head = '/Volumes/PWFA 4big';
isscan = 1;
doyag = 1;
doceloss = 1;
docegain = 1;

data = load_E200_data(ar_qs,head,isscan,doyag,doceloss,docegain);
%[eta_max, eta_cent, eta_fmin, eta_fmax] = DISPANA(data);
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
pyros = [pyros; data.epics.BLEN_LI20_3014_BRAW(data.YAG.EPID_ind(data.YAG.py_ind(:,i),i),i)];
yag_specs = [yag_specs data.YAG.spectra(:,data.YAG.py_ind(:,i),i)];
cel_specs = [cel_specs data.CELOSS.spec(:,data.YAG.py_ind(:,i),i)];
ceg_specs = [ceg_specs data.CEGAIN.spec(:,data.YAG.py_ind(:,i),i)];

end
