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
    x_line = xx(PixLim(3):PixLim(4));
    N_pix = length(x_line);
    
    % allocate image matrices
    lineouts = zeros(N_pix,n_shot,n_step);
    
    % allocate pulse_id
    pulse_id = zeros(n_shot,n_step);
    
    % get image data
    for i = 1:n_step
        
        % read images
        [im_dat,~,pulse_id(:,i)] = E200_readImages(scan_info(i).YAG);
        % get bg
        back = uint16(rot90(d(1).cam_back.YAG.img,2));
        backs = repmat(back,[1,1,size(im_dat,3)]);
        
        % orient images and subtract bg
        im_dat = flipdim(flipdim(im_dat,2),1) - backs;
        
        % all line outs
        lines = squeeze(mean(im_dat(PixLim(2):PixLim(1),PixLim(3):PixLim(4),:),1));
        lineouts(:,:,i) = lines;
        
        % all image data
        clear('im_dat');
        
    end
    
end

YAG.spectra = lineouts;
YAG.pulse_id = pulse_id;
YAG.axis = x_line';