function YAG = EXTRACT_YAG(path_file,head,isscan)

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
        scan_info(i).YAG = [head scan_info(i).YAG];
    end
    
    % Get YAG information
    res = 9.62;
    x_ctr = 671;
    y_ctr = 416;
    LineLim = [0.6; 1.0; -3.5; 3.8;];
    
    % create image axes
    xx = res/1000*((1:1392)-x_ctr);
    yy = -res/(sqrt(2)*1000)*((1:1040)-y_ctr);
    
    % Match lineout limits to axes
    if LineLim(1) < yy(end); LineLim(1) = yy(end)+0.1; end;
    if LineLim(2) > yy(1); LineLim(2) = yy(1)-0.1; end;
    if LineLim(3) < xx(1); LineLim(3) = xx(1)+0.1; end;
    if LineLim(4) > xx(end); LineLim(4) = xx(end)-0.1; end;
    
    % Map lineout to pixels
    PixLim(1) = round(-LineLim(1)*1000*sqrt(2)/res)+y_ctr;
    PixLim(2) = round(-LineLim(2)*1000*sqrt(2)/res)+y_ctr;
    PixLim(3) = round(LineLim(3)*1000/res)+x_ctr;
    PixLim(4) = round(LineLim(4)*1000/res)+x_ctr;
    
    % create lineout axis
    x_line = xx(PixLim(3):PixLim(4))';
    N_pix = length(x_line);
    
    % allocate image matrices
    lineouts = zeros(N_pix,n_shot*n_step);
    
    % allocate pulse_id
    pulse_id = zeros(n_shot*n_step,1);
    pix_sum = zeros(n_shot*n_step,1);
    fwhms = zeros(n_shot*n_step,1);
    x_lo = zeros(n_shot*n_step,1);
    x_hi = zeros(n_shot*n_step,1);
    x_max = zeros(n_shot*n_step,1);
    x_cent = zeros(n_shot*n_step,1);
    x_rms = zeros(n_shot*n_step,1);
    scan_val = zeros(n_shot*n_step,1);
    
    % get bg
    back = uint16(rot90(d(1).cam_back.YAG.img,2));
    backs = repmat(back,[1,1,n_shot]);
    
    % get image data
    for i = 1:n_step
        
        start_ind = (i-1)*n_shot;
        end_ind = i*n_shot;
        
        % fill in scan value
        scan_val((start_ind+1):end_ind) = scan_info(i).Control_PV;
        
        % read images
        [im_dat,~,pulse_id((start_ind+1):end_ind)] = E200_readImages(scan_info(i).YAG);
        
        % orient images and subtract bg
        im_dat = flipdim(flipdim(im_dat,2),1) - backs;
        
        % sum pixels
        pix_sum((start_ind+1):end_ind) = sum(sum(im_dat));
        
        % all line outs
        lines = squeeze(mean(im_dat(PixLim(2):PixLim(1),PixLim(3):PixLim(4),:),1));
        lineouts(:,(start_ind+1):end_ind) = lines;
        
        for j = 1:n_shot
            [fwhms(start_ind+j),low_i,high_i] = FWHM(x_line,lineouts(:,start_ind+j));
            x_lo(start_ind+j) = x_line(low_i);            
            x_hi(start_ind+j) = x_line(high_i);            
            [~,max_pos] = max(lineouts(:,start_ind+j));
            x_max(start_ind+j) = x_line(max_pos);
            x_cent(start_ind+j) = sum(x_line.*lineouts(:,start_ind+j))./sum(lineouts(:,start_ind+j));
            x_rms(start_ind+j) = sqrt(sum(lineouts(:,start_ind+j).*(x_line - x_cent(start_ind+j)).^2.)/sum(lineouts(:,start_ind+j)));
        end
        % all image data
        clear('im_dat');
        
    end
    
end

YAG.spectra  = lineouts;
YAG.pulse_id = pulse_id;
YAG.axis     = x_line;
YAG.pix_sum  = pix_sum;
YAG.fwhms    = fwhms;
YAG.x_lo     = x_lo;
YAG.x_hi     = x_hi;
YAG.x_max    = x_max;
YAG.x_cent   = x_cent;
YAG.x_rms    = x_rms;
YAG.scan_val = scan_val;