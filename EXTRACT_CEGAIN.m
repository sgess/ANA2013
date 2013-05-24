function CEGAIN = EXTRACT_CEGAIN(path_file,head,isscan)

pixels = 1:1392;

% Get Axis
E = E200_cher_get_E_axis('20130423', 'CEGAIN', 0, pixels);

% load data containing filenames
load(path_file);

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
    n_step = length(d);
    n_shot = d(1).param.n_shot;
    
    for i = 1:n_step
        scan_info(i).CEGAIN = [head scan_info(i).CEGAIN];
    end
    
    N_pix = length(pixels);
    
    % allocate image matrices
    spec = zeros(N_pix,n_shot,n_step);
    filt_spec = zeros(N_pix,n_shot,n_step);
    filt_spec2 = zeros(N_pix,n_shot,n_step);
    
    % allocate pulse_id
    pulse_id = zeros(n_shot,n_step);
    
    % background
    back = uint16(d(1).cam_back.CEGAIN.img);
    backs = repmat(back,[1,1,n_shot]);
    
    % allocate scalar data
    E_ACC = zeros(n_shot,n_step);
    E_UNAFFECTED = zeros(n_shot,n_step);
    E_EMAX = zeros(n_shot,n_step);
    E_EMAX_ind = zeros(n_shot,n_step);
    E_EMAX2 = zeros(n_shot,n_step);
    E_EMAX_ind2 = zeros(n_shot,n_step);
    E_EMAX3 = zeros(n_shot,n_step);
    E_EMAX_ind3 = zeros(n_shot,n_step);
    
    % get image data
    for i = 1:n_step
        
        % read images
        [im_dat,~,pulse_id(:,i)] = E200_readImages(scan_info(i).CEGAIN);
        im_dat = im_dat - backs;
        
        for j = 1:n_shot
            
            img = double(im_dat(:,:,j));
            band = mean(img(:,E>50),2);
            bands = repmat(band,[1,N_pix]);
            img = img - bands;
            spec(:,j,i) = sum(img,1);
            
            r = 3;     % Adjust for desired window size
            k = 5;     % Select the kth largest element
            
            A = zeros([size(img), r^2]);
            for s=1:r^2
                w = zeros(r);
                w(s) = 1;
                A(:,:,s) = filter2(w, img);
            end
            B = sort(A,3);
            img2 = squeeze(B(:,:,k));
            
            E_ACC(j,i) = sum(sum(img2(:,E>22)));
            E_UNAFFECTED(j,i) = sum(sum(img2(:,E<22 & E>19)));
            
            filt_img = filter2(ones(25,10)/250, img2);
            filt_spec(:,j,i) = sum(filt_img,1);
            filt_spec2(:,j,i) = max(filt_img);
            
            a = cumsum(sum(filt_img(:,51:1392),1));
            ind = find(a>0.995*max(a), 1);
            if isempty(ind)
                E_EMAX(j,i) = E(1);
                E_EMAX_ind(j,i) = 1;
            else
                E_EMAX(j,i) = E(ind);
                E_EMAX_ind(j,i) = ind;
            end
            
            b = cumsum(filt_spec2(:,j,i));
            ind = find(b>0.995*max(b), 1);
            if isempty(ind)
                E_EMAX2(j,i) = E(1);
                E_EMAX_ind2(j,i) = 1;
            else
                E_EMAX2(j,i) = E(ind);
                E_EMAX_ind2(j,i) = ind;
            end
            
            ind = find(filt_spec2(:,j,i)>10., 1, 'last');
            if isempty(ind)
                E_EMAX3(j,i) = E(1);
                E_EMAX_ind3(j,i) = 1;
            else
                E_EMAX3(j,i) = E(ind);
                E_EMAX_ind3(j,i) = ind;
            end
        end
        clear('im_dat');
        
    end
    
end

CEGAIN.pulse_id = pulse_id;
CEGAIN.axis = E';
CEGAIN.spec = spec;
CEGAIN.filt_spec = filt_spec;
CEGAIN.filt_spec2 = filt_spec2;
CEGAIN.E_ACC = E_ACC;
CEGAIN.E_UNAFFECTED = E_UNAFFECTED;
CEGAIN.E_EMAX = E_EMAX;
CEGAIN.E_EMAX_ind = E_EMAX_ind;
CEGAIN.E_EMAX2 = E_EMAX2;
CEGAIN.E_EMAX_ind2 = E_EMAX_ind2;
CEGAIN.E_EMAX3 = E_EMAX3;
CEGAIN.E_EMAX_ind3 = E_EMAX_ind3;