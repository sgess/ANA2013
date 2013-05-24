function CELOSS = EXTRACT_CELOSS(path_file,head,isscan)

pixels = 1:1392;

% Get Axis
E = E200_cher_get_E_axis('20130423', 'CELOSS', 0, pixels);

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
        scan_info(i).CELOSS = [head scan_info(i).CELOSS];
    end
    
    N_pix = length(pixels);

    % allocate image matrices   
    spec = zeros(N_pix,n_shot,n_step);

    % allocate pulse_id
    pulse_id = zeros(n_shot,n_step);
    
    % background
    back = uint16(d(1).cam_back.CELOSS.img);
    backs = repmat(back,[1,1,n_shot]);
    
    % allocate scalar data
    E_DECC = zeros(n_shot,n_step);
    E_UNAFFECTED2 = zeros(n_shot,n_step);
    E_EMIN = zeros(n_shot,n_step);
    
    % get image data
    for i = 1:n_step
        
        % read images
        [im_dat,~,pulse_id(:,i)] = E200_readImages(scan_info(i).CELOSS);
        im_dat = im_dat - backs;
        
        for j = 1:n_shot
            
            img = double(im_dat(:,:,j));
            spec(:,j,i) = sum(img,1);
            
            E_DECC(j,i) = sum(sum(img(:,E<19)));
            E_UNAFFECTED2(j,i) = sum(sum(img(:,E<22 & E>19)));

            a = cumsum(sum(img(:,51:1392),1));
            ind = find(a<0.01*max(a), 1, 'last');
            if isempty(ind)
                E_EMIN(j,i) = E(1392);
            else
                E_EMIN(j,i) = E(ind);
            end
        end
    end
end

CELOSS.pulse_id = pulse_id;
CELOSS.axis = E';
CELOSS.spec = spec;
CELOSS.E_DECC = E_DECC;
CELOSS.E_UNAFFECTED2 = E_UNAFFECTED2;
CELOSS.E_EMIN = E_EMIN;