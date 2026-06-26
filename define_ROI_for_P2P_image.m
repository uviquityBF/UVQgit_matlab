
function A = define_ROI_for_P2P_image( )
    %clear & get files
    clear; close all;
    program_name = 'define ROI for P2P'
%     dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL01\id3.7_vary_align_v_Z\P2P';
    dirs.d = 'G:\Shared drives\Corp Main\Engineering\LAB Studies\on Samples\GaugeLots\GL00\am123\sd01\Standardized Testing';
    
    
    cd(dirs.d); [filenames,PathName] = uigetfile({'*.tiff'},program_name,'Multiselect','on');
    Nf = length(filenames);
    if isstr(filenames)
        filename = filenames; clear filenames;
        filenames{1} = filename;
        Nf=1;
    end
    for kf=1:Nf
        FileName=filenames{kf};
        [output_str,Ps_x,Ps_x] = select_ROI_from_Image(PathName,FileName);
        if kf==1    fID = fopen([PathName,'\tmp_ROIs.csv'],'w'); end
        fprintf(fID,[output_str,'\r\n']);
    end
    fclose(fID);
    
end

function [output_str,Ps_x,Ps_y] = select_ROI_from_Image(PathName,FileName)
     program_name = 'define ROI for P2P'

     %open TIFF  file 
    if strcmp(FileName(end-3:end),'tiff')  %for TIFF file input
        I0 = imread([PathName,FileName]);
        [Nr,Nc,Nlayers] = size(I0);
        if (Nlayers>=3)
            I1 = I0(:,:,3);                  %take Blue layer
        else
            I1 = I0;
        end
    end
    
    % show image and select it
    hf = figure;    himg  = imagesc(log10(double(I1))); 
    h = msgbox('Click opposite Corners of ROI',program_name);
    uiwait(h);
    [Ps_x(1),Ps_y(1)] = ginput(1);
    hold all; plot(Ps_x(1),Ps_y(1),'+');
    [Ps_x(2),Ps_y(2)] = ginput(1);
    hold all; plot(Ps_x(2),Ps_y(2),'+');
  
    Ps_x = round(sort(Ps_x),0);   
    Ps_y = round(sort(Ps_y),0);
 
    output_str = [FileName,', ',num2str(Ps_x(1)),', ',num2str(Ps_y(1)),', ',num2str(Ps_x(2)),', ',num2str(Ps_y(2))];
    disp(output_str)  
    close(hf);
end
