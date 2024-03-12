% SEP
sep_flag       = 1; % input('Interdigitations (0) or SEP (1)? ');
if sep_flag == 1
    
    disp('---> SEP selected')
    while isempty(sep_type)
        disp('Error: input type of exit pathways')
        sep_type = input('Type of SEP: canal (0) or opening (1)? ');
    end
    
    if sep_type == 0
        %%% Canal version
        if nSEP == 1
            ellipse( y_c-sep_half_width:y_c+sep_half_width, x_c:sep_right_c-1) =1; % central SE
        elseif nSEP == 2
            sep_high_2 = y_c + 30;
            sep_low_2 = y_c - 30;
            
            ellipse( sep_low_2-sep_half_width:sep_low_2+sep_half_width, x_c:sep_right_m-1) = 1; %high mid SEP
            ellipse( sep_high_2-sep_half_width:sep_high_2+sep_half_width, x_c:sep_right_m-1) = 1; %low mid SEP
            
        elseif nSEP == 3
            
            sep_high_m = y_c + 40; %floor(y_r * 2/5);
            sep_low_m = y_c - 40; %floor(y_r * 2/5);
            
            ellipse( y_c-sep_half_width:y_c+sep_half_width, x_c:sep_right_c-1) =1; % central SEP
            ellipse( sep_low_m-sep_half_width:sep_low_m+sep_half_width, x_c:sep_right_m-1) = 1; %high mid SEP
            ellipse( sep_high_m-sep_half_width:sep_high_m+sep_half_width, x_c:sep_right_m-1) = 1; %low mid SEP
            
        elseif nSEP == 4
            sep_high = y_c + 45;
            sep_low = y_c - 45;
            sep_high_m = y_c + 15;
            sep_low_m = y_c - 15;
            
            ellipse( sep_low_m-sep_half_width:sep_low_m+sep_half_width, x_c:sep_right_m-1) = 1; %high mid SEP
            ellipse( sep_high_m-sep_half_width:sep_high_m+sep_half_width, x_c:sep_right_m-1) = 1; %low mid SEP
            ellipse( sep_low-sep_half_width:sep_low+sep_half_width, x_c:sep_right_h-1) = 1; % high SEP
            ellipse( sep_high-sep_half_width:sep_high+sep_half_width, x_c:sep_right_h-1) = 1; %low SEP
            
        elseif nSEP == 5
            sep_high = y_c + 60; %floor(y_r * 2/3);
            sep_high_m = y_c + 30; %floor(y_r * 1/3);
            sep_low_m = y_c - 30; %floor(y_r * 1/3);
            sep_low = y_c - 60; %floor(y_r * 2/3);
            
            ellipse( y_c-sep_half_width:y_c+sep_half_width, x_c:sep_right_c-1)               = 1; % central SEP
            ellipse( sep_low_m-sep_half_width:sep_low_m+sep_half_width, x_c:sep_right_m-1)   = 1; %high mid SEP
            ellipse( sep_high_m-sep_half_width:sep_high_m+sep_half_width, x_c:sep_right_m-1) = 1; %low mid SEP
            ellipse( sep_low-sep_half_width:sep_low+sep_half_width, x_c:sep_right_h-1)       = 1; % high SEP
            ellipse( sep_high-sep_half_width:sep_high+sep_half_width, x_c:sep_right_h-1)     = 1; %low SEP
            
            
        elseif nSEP == 6
            sep_high_top_6 = y_c + 60; %floor(y_r * 4/5);
            sep_high_6 = y_c + 35; %floor(y_r * 2/5);
            sep_high_m_6 = y_c + 12; %floor(y_r * 2/5);
            sep_low_m_6 = y_c - 12; %floor(y_r * 2/5);
            sep_low_6 = y_c - 35; %floor(y_r * 2/5);
            sep_low_top_6 = y_c - 60; %floor(y_r * 4/5); % coordinates for the SEPs to be used below
            
            ellipse( sep_low_m_6-sep_half_width:sep_low_m_6+sep_half_width, x_c:sep_right_m-1) = 1; %high mid SEP
            ellipse( sep_high_m_6-sep_half_width:sep_high_m_6+sep_half_width, x_c:sep_right_m-1) = 1; %low mid SEP
            ellipse( sep_low_6-sep_half_width:sep_low_6+sep_half_width, x_c:sep_right_h-1) = 1; % high SEP
            ellipse( sep_high_6-sep_half_width:sep_high_6+sep_half_width, x_c:sep_right_h-1) = 1; %low SEP
            ellipse( sep_high_top_6-sep_half_width:sep_high_top_6+sep_half_width, x_c:sep_right_top-1) = 1; %low mid SEP
            ellipse( sep_low_top_6-sep_half_width:sep_low_top_6+sep_half_width, x_c:sep_right_top-1) = 1; %low mid SEP
            
        elseif nSEP == 7
            sep_high_top = y_c + 60; %floor(y_r * 4/5);
            sep_high = y_c + 40; %floor(y_r * 2/5);
            sep_high_m = y_c + 20; %floor(y_r * 2/5);
            sep_low_m = y_c - 20; %floor(y_r * 2/5);
            sep_low = y_c - 40; %floor(y_r * 2/5);
            sep_low_top = y_c - 60; %floor(y_r * 4/5); % coordinates for the SEPs to be used below
            
            ellipse( y_c-sep_half_width:y_c+sep_half_width, x_c:sep_right_c-1) =1; % central SEP
            ellipse( sep_low_m-sep_half_width:sep_low_m+sep_half_width, x_c:sep_right_m-1) = 1; %high mid SEP
            ellipse( sep_high_m-sep_half_width:sep_high_m+sep_half_width, x_c:sep_right_m-1) = 1; %low mid SEP
            ellipse( sep_low-sep_half_width:sep_low+sep_half_width, x_c:sep_right_h-1) = 1; % high SEP
            ellipse( sep_high-sep_half_width:sep_high+sep_half_width, x_c:sep_right_h-1) = 1; %low SEP
            ellipse( sep_high_top-sep_half_width:sep_high_top+sep_half_width, x_c:sep_right_top-1) = 1; %low mid SEP
            ellipse( sep_low_top-sep_half_width:sep_low_top+sep_half_width, x_c:sep_right_top-1) = 1; %low mid SEP
        end
        
        %%%Costruzione del bordo del SNA
        san_border = logical(imdilate(ellipse, strel('diamond', 1)) - ellipse); %% border
        if nSEP == 1
            san_border([y_c-sep_half_width:y_c+sep_half_width, y_c-sep_half_width:y_c+sep_half_width],sep_right_c) = 0;
        elseif nSEP == 2
            san_border([sep_low_2-sep_half_width:sep_low_2+sep_half_width, sep_high_2-sep_half_width:sep_high_2+sep_half_width], sep_right_m) = 0;
        elseif nSEP == 3
            san_border([y_c-sep_half_width:y_c+sep_half_width, y_c-sep_half_width:y_c+sep_half_width-1],sep_right_c) = 0;
            san_border([sep_low_m-sep_half_width:sep_low_m+sep_half_width, sep_high_m-sep_half_width:sep_high_m+sep_half_width], sep_right_m) = 0;
        elseif nSEP == 4
            san_border([sep_low_m-sep_half_width:sep_low_m+sep_half_width, sep_high_m-sep_half_width:sep_high_m+sep_half_width], sep_right_m) = 0;
            san_border([sep_low-sep_half_width:sep_low+sep_half_width, sep_high-sep_half_width:sep_high+sep_half_width], sep_right_h) = 0;
        elseif nSEP == 5
            san_border([y_c-sep_half_width:y_c+sep_half_width, y_c-sep_half_width:y_c+sep_half_width],sep_right_c) = 0;
            san_border([sep_low_m-sep_half_width:sep_low_m+sep_half_width, sep_high_m-sep_half_width:sep_high_m+sep_half_width], sep_right_m) = 0;
            san_border([sep_low-sep_half_width:sep_low+sep_half_width, sep_high-sep_half_width:sep_high+sep_half_width], sep_right_h) = 0; % open the right SEPs
        elseif nSEP == 6
            san_border([sep_low_m_6-sep_half_width:sep_low_m_6+sep_half_width, sep_high_m_6-sep_half_width:sep_high_m_6+sep_half_width], sep_right_m) = 0;
            san_border([sep_low_6-sep_half_width:sep_low_6+sep_half_width, sep_high_6-sep_half_width:sep_high_6+sep_half_width], sep_right_h) = 0; % open the right SEPs
            san_border([sep_low_top_6-sep_half_width:sep_low_top_6+sep_half_width, sep_high_top_6-sep_half_width:sep_high_top_6+sep_half_width], sep_right_top) = 0;
        elseif nSEP == 7
            san_border([y_c-sep_half_width:y_c+sep_half_width, y_c-sep_half_width:y_c+sep_half_width],sep_right_c) = 0;
            san_border([sep_low_m-sep_half_width:sep_low_m+sep_half_width, sep_high_m-sep_half_width:sep_high_m+sep_half_width], sep_right_m) = 0;
            san_border([sep_low-sep_half_width:sep_low+sep_half_width, sep_high-sep_half_width:sep_high+sep_half_width], sep_right_h) = 0; % open the right SEPs
            san_border([sep_low_top-sep_half_width:sep_low_top+sep_half_width, sep_high_top-sep_half_width:sep_high_top+sep_half_width], sep_right_top) = 0;
        end
        
        %%% ENCIRCLE LPMs
        %         san_border(35:55, x_c) = 1;
        %         san_border(35, x_c:115) = 1;
        %         san_border(55, x_c:111) = 1;
        %
        %         san_border(65:115, x_c) = 1;
        %         san_border(65, x_c:115) = 1;
        %         san_border(115, x_c:112) = 1;
        %
        %         san_border(125:166, x_c) = 1;
        %         san_border(125, x_c:115) = 1;
        %         san_border(166, x_c:115) = 1;
        %%%
        
        disp('---> Canal SEP selected')
        sep_str = 'SEPc';
        
    elseif sep_type == 1
        %%% Ellipse opening version
        san_border = logical(imdilate(ellipse, strel('diamond', 1)) - ellipse); %% border
        %%% NO SAN OUTPUT
        %         san_border([sep_low-sep_half_width:sep_low+sep_half_width, sep_high-sep_half_width:sep_high+sep_half_width], :) = 0; % open the left SEPs
        %         san_border([sep_low-sep_half_width:sep_low+sep_half_width, sep_high-sep_half_width:sep_high+sep_half_width], :) = 0; % open the right SEPs
        
        disp('---> open SEP selected')
        sep_str = 'SEPo';
    end
    
elseif sep_flag == 0
    sep_str = 'INT';
    disp('---> Interdigitations selected')
end
% SAN
atrial_tissue(ellipse) =  idx_san; % SAN (inside the border)
atrial_tissue(san_border) = idx_fat; % SAN border

atrial_tissue(100, 67) = idx_atr;
atrial_tissue(100, 68) = idx_fat;
if sep_type == 1
    atrial_tissue(100, 93) = idx_atr;
    atrial_tissue(100, 92) = idx_fat;
end
