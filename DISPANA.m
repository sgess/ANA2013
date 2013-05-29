function [eta_max, eta_cent, eta_fmin, eta_fmax] = DISPANA(data,save_dir,savE)

steps = length(unique(data.EPICS.scan_step));
lineouts = data.YAG.spectra;
line_out = zeros(length(data.YAG.axis),steps);
for i = 1:steps
    line_out(:,i) = mean(data.YAG.spectra(:,data.YAG.scan_step == i),2);
end


x_line = data.YAG.axis;

energy = unique(data.EPICS.scan_val);
delta = energy/20.35e3;
dE = data.EPICS.scan_val/20.35e3;
dY = data.YAG.scan_val/20.35e3;

% Dispersion for spectrum maximum
[max_spec,max_pos] = max(line_out);
x_max = x_line(max_pos);
d_max = polyfit(delta,x_max,2);
eta_max = d_max(2);
max_fit = d_max(1)*delta.^2 + d_max(2)*delta + d_max(3);

% Dispersion for spectrum centroid
x_mat = repmat(x_line,1,steps);
cent_spec = sum(x_mat.*line_out)./sum(line_out);
cent_ind = superfind(cent_spec,x_line');
cent_val = diag(line_out(cent_ind,:));
d_cent = polyfit(delta',cent_spec,2);
eta_cent = d_cent(2);
cent_fit = d_cent(1)*delta.^2 + d_cent(2)*delta + d_cent(3);

% Dispersion for FWHM edge
fwhms = zeros(1,steps);
low_i = zeros(1,steps);
high_i = zeros(1,steps);
for i = 1:steps; [fwhms(i),low_i(i),high_i(i)] = FWHM(x_line,line_out(:,i)); end;

fx_min = x_line(low_i);
fmin_val = diag(line_out(low_i,:));
d_fmin = polyfit(delta,fx_min,2);
eta_fmin = d_fmin(2);
fmin_fit = d_fmin(1)*delta.^2 + d_fmin(2)*delta + d_fmin(3);

fx_max = x_line(high_i);
fmax_val = diag(line_out(high_i,:));
d_fmax = polyfit(delta,fx_max,2);
eta_fmax = d_fmax(2);
fmax_fit = d_fmax(1)*delta.^2 + d_fmax(2)*delta + d_fmax(3);

figure(1);
subplot(2,1,1);
plot(x_line,line_out);
legend(strcat(num2str(energy,'%0.1f'),' MeV'),'location','northwest');
hold on;
plot(x_max,max_spec,'k*',cent_spec,cent_val,'r*',fx_min,fmin_val,'g*',fx_max,fmax_val,'c*');
axis([x_line(1) x_line(end) 0 max(max(line_out))+10]);
hold off;
xlabel('X (mm)','fontsize',14);
title('Average of SYAG Spectra','fontsize',16);
subplot(2,1,2);
plot(delta,x_max,'k*',delta,cent_spec,'r*',delta,fx_min,'g*',delta,fx_max,'c*');
legend('Spectrum Maximum','Spectrum Centroid','FWHM Low','FWHM High','location','northwest');
hold on;
plot(delta,max_fit,'k:',delta,cent_fit,'r:',delta,fmin_fit,'g:',delta,fmax_fit,'c:');
hold off;
xlabel('\delta','fontsize',14);
ylabel('X (mm)','fontsize',14);
title(['\eta_{max} = ' num2str(eta_max,'%0.2f') ', \eta_{cent} = ' num2str(eta_cent,'%0.2f') ...
    ', \eta_{low} = ' num2str(eta_fmin,'%0.2f') ', \eta_{high} = ' num2str(eta_fmax,'%0.2f')],'fontsize',16);
if savE; saveas(gca,[save_dir 'yag_disp.pdf']); end;

%shots = length(data.aida.EPID_ind);
%shots = length(data.YAG.EPID_ind);
shots = length(data.YAG.pulse_id);

%dd = repmat(delta,1,shots);
%de = dd';

% Dispersion for spectrum maximum
[~,Smax_pos] = max(lineouts,[],1);
Smax_pos = squeeze(Smax_pos);
Sx_max = x_line(Smax_pos);
%Sd_max = polyfit(de,Sx_max,2);
Sd_max = polyfit(dY,Sx_max,2);
Seta_max = Sd_max(2);
Smax_fit = Sd_max(1)*delta.^2 + Sd_max(2)*delta + Sd_max(3);

% Dispersion for spectrum centroid
Sx_mat = repmat(x_line,[1,shots]);
Scent_spec = sum(Sx_mat.*lineouts)./sum(lineouts);
Scent_spec = squeeze(Scent_spec);
Sd_cent = polyfit(dY,Scent_spec',2);
Seta_cent = Sd_cent(2);
Scent_fit = Sd_cent(1)*delta.^2 + Sd_cent(2)*delta + Sd_cent(3);

% Dispersion for FWHM
fwhms = zeros(shots,1);
low_i = zeros(shots,1);
high_i = zeros(shots,1);
for j = 1:shots
    [fwhms(j),low_i(j),high_i(j)] = FWHM(x_line,lineouts(:,j));
end
Sx_fmin = x_line(low_i);
Sx_fmax = x_line(high_i);
Sd_fmin = polyfit(dY,Sx_fmin,2);
Sd_fmax = polyfit(dY,Sx_fmax,2);
Seta_fmin = Sd_fmin(2);
Seta_fmax = Sd_fmax(2);
Sfmin_fit = Sd_fmin(1)*delta.^2 + Sd_fmin(2)*delta + Sd_fmin(3);
Sfmax_fit = Sd_fmax(1)*delta.^2 + Sd_fmax(2)*delta + Sd_fmax(3);

%ax_2050 = data.aida.BPMS_LI20_2050.x;
%ax_2445 = data.aida.BPMS_LI20_2445.x;
% ex_2445 = [];
% for i=1:steps
%     ex_2445 = [ex_2445, data.epics.BPMS_LI20_2445_X(data.YAG.EPID_ind(:,i),i)];
% end
[~,ii,ie] = intersect(data.YAG.UID,data.EPICS.UID);
ex_2445 = data.EPICS.BPMS_LI20_2445_X(ie);

% Dispersion from BPMs
%pax_2050=polyfit(de,ax_2050,2);
%pax_2445=polyfit(de,ax_2445,2);
pex_2445=polyfit(dY,ex_2445,2);
%ax_2050_fit = pax_2050(1)*delta.^2 + pax_2050(2)*delta + pax_2050(3);
%ax_2445_fit = pax_2445(1)*delta.^2 + pax_2445(2)*delta + pax_2445(3);
ex_2445_fit = pex_2445(1)*delta.^2 + pex_2445(2)*delta + pex_2445(3);



figure(2);
subplot(2,1,1);
%plot(de(:),Sx_max(:),'k*',de(:),Scent_spec(:),'r*',de(:),Sx_fmin(:),'g*',de(:),Sx_fmax(:),'c*');
plot(dY,Sx_max,'k*',dY,Scent_spec,'r*',dY,Sx_fmin,'g*',dY,Sx_fmax,'c*');
l=legend('SYAG Max','SYAG Cent','SYAG FWHM Lo','SYAG FWHM Hi','location','northwest');
hold on;
plot(delta,Smax_fit,'k',delta,Scent_fit,'r',delta,Sfmin_fit,'g',delta,Sfmax_fit,'c');
hold off;
xlabel('\delta','fontsize',14);
ylabel('X (mm)','fontsize',14);
title(['Fit to YAG Data: \eta_{max} = ' num2str(Seta_max,'%0.2f') ', \eta_{cent} = ' num2str(Seta_cent,'%0.2f')...
    ', \eta_{low} = ' num2str(Seta_fmin,'%0.2f') ', \eta_{high} = ' num2str(Seta_fmax,'%0.2f')],'fontsize',16);



subplot(2,1,2);
plot(dY,ex_2445,'m*');
l=legend('BPM 2445 EPICS','location','northwest');
hold on;
plot(delta,ex_2445_fit,'c');
hold off
xlabel('\delta','fontsize',14);
ylabel('X (mm)','fontsize',14);
title(['Fit to BPM Data: \eta_{2445e} = ' num2str(pex_2445(2),'%0.2f')],'fontsize',16);
if savE; saveas(gca,[save_dir 'bpm_disp.pdf']); end;

% subplot(2,1,2);
% plot(de(:),ax_2050(:),'g*',de(:),ax_2445(:),'b*',de(:),ex_2445(:),'m*');
% l=legend('BPM 2050','BPM 2445 AIDA','BPM 2445 EPICS','location','northwest');
% hold on;
% plot(delta,ax_2050_fit,'k',delta,ax_2445_fit,'r',delta,ex_2445_fit,'c');
% hold off
% xlabel('\delta','fontsize',14);
% ylabel('X (mm)','fontsize',14);
% title(['Fit to BPM Data: \eta_{2050} = ' num2str(pax_2050(2),'%0.2f') ', \eta_{2445a} = ' num2str(pax_2445(2),'%0.2f')...
%     ', \eta_{2445e} = ' num2str(pex_2445(2),'%0.2f')],'fontsize',16);
% if savE; saveas(gca,[save_dir 'bpm_disp.pdf']); end;