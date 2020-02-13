function [bvec,TR,TE]=readb(folder)
    if nargout == 1
        docs=dir(fullfile(folder,'*.bval'));
        if ~isempty(docs)
            data = fileread(fullfile(folder,docs.name(1)));
            bvec=str2num(data);
            return
        end
    end
    
    files=dir(fullfile(folder,'*.dcm'));
    files = {files.name}';
    sfiles=size(files,1);
    
    infolast=dicominfo(fullfile(folder,files{end}));
    numberofslices=infolast.Private_2001_100a;
    TE=infolast.EchoTime;
    TR=infolast.RepetitionTime;

    bvec=zeros(sfiles/numberofslices,1);
        
    for ii=1:sfiles/numberofslices
        info=dicominfo(fullfile(folder,files{ii}));
        bvalsqrt=sqrt(info.DiffusionBValue);
        bvec(ii)=(bvalsqrt*info.DiffusionGradientOrientation(1))^2+(bvalsqrt*info.DiffusionGradientOrientation(2))^2+(bvalsqrt*info.DiffusionGradientOrientation(3))^2;
    end        
    bvec=round(bvec);
end


