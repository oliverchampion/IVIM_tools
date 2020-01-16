function [S,f,D,SSE,rsquare] = IVIMfixed(selection,bvec,TR,TE,allb)

%%%%%%%%
% This program calculates S, f, D, SSE and rsquare for each voxel in the
% input matrix "selection". The program uses the IVIM-fixed aproach from
% the article, meaning D* is fixed to a set value (0.07 mm2/s in this case)
%
% Input:
% Selection = n*m matrix, where the first dimension includes the data from
% the different voxels and the second dimension is the data as function of
% bvalues.
%
% bvec = a vector containing the b-values of each data point. Must be m
% long
%
% TR = repetition time of the acquisition: used to calculate tissue fraction
% instead of signal fraction
%
% TE = echo time of the acquisition: used to calculate tissue fraction
% instead of signal fraction
%
% allb = 1 or 0; 1 means all b-values are taken along seperately in the
% fit, whereas 0 means b-value data is grouped and weighted according to
% the number of acquisitions per b-value. This is only important for data
% where multiple acquisitions occured per b-value (i.e. DTI data)
%
%
% Output (m long vectors, repressenting the voxel values):
% S = signal at b=0 s/mm2 (fitted)
%
% f = perfusion fraction (tissue fraction)
%
% D = diffusion coefficient
%
% SSE = sum square due to error
%
% rsquare = adjusted R-squared
%
% We advice not to take along all data (i.e. let data be a o*p*q*m 
% matrix which repressents a 3D o*p*q voxel space with a 4th b-value 
% dimension (m); 
% alldata=reshape(data,size(data,1)*size(data,2)*size(data,3),size(data,4)), 
% but to mask your data before
% using IVIM fixed (i.e. selection=alldata(mask)).
%%%%%%%%%

%% This are tissue specific parameters. They are set for panctreas. May be changed for other organs.
if nargin==2
    TR=100000;
    TE=0;
    allb=0;
end

Dstergeuss=0.07; % this will be the fixed value of D*
% T1=725; %--> http://doi.org/10.1148/radiol.2303021331
% T2=43; %--> http://doi.org/10.1148/radiol.2303021331
% T1b=1932; %--> http://doi.org/10.1002/mrm.20605
% T2b=275; %--> http://doi.org/10.1002/mrm.20605
T1=10000; 
T2=0.001; 
T1b=10000; 
T2b=0.001; 

% fit constraints [initial guess, min, max]
fcon=[0.1, 0.001, 1];
Dcon=[0.0017, 0.0005, 0.006];

%% here the data is sorted in blocks of b-values. If allb=0, the average signal intensity per b-value is taken for repeated measures in that b-value in order to speed fitting.
fail=0;
if allb==1
    [bvec, order]=sort(bvec);
    selection=selection(:,order);
    selection=transpose(selection);
else
    a=unique(bvec);
    
    weights=zeros(size(a,2),1);
    
    for ii=1:size(weights,1)
        weights(ii)=sum(bvec==a(ii));
    end
    weights=round(weights/min(weights));
    
    bvec=a;
    selection=transpose(selection);
    
    
    qq=1;
    for ii=1:size(bvec,2)
        for kk=1:weights(ii)
            bvec2(qq)=bvec(ii);
            selection2(qq,:)=selection(ii,:);
            qq=qq+1;
        end
    end
    %% updating data en b-vector in case data is averaged
    bvec=bvec2;
    selection=selection2;
    clear selection2 bvec2
end

%% initiating parameters
ssel=size(selection);
S=zeros(ssel(2),1);
f=zeros(ssel(2),1);
D=zeros(ssel(2),1);
SSE=zeros(ssel(2),1);
rsquare=zeros(ssel(2),1);

bvecbu=transpose(bvec);
%% looping over voxels. Can be parfor loop in case of 1 patient. I use parfor over the patients to minimize overhead.
for k=1:ssel(2)
    try
        bvec=bvecbu;
        data1=selection(:,k);
        % when data is missing (0), through away in fit. Data can be missing due to registration of data at edge of FOV for specific b-values. Furthermore, data can be masked for bad slices effected by heartbeats. A masking method was described in http://doi.org/10.1097/RLI.0000000000000225
        bvec(data1==0)=[];
        data1(data1==0)=[];
        
        % fitoptions, including constraints
        sbiexpf=fitoptions('Method','NonlinearLeastSquares','robust','on','Startpoint',[Dcon(1) data1(1) fcon(1)],'Lower',[Dcon(2) 0.1*data1(1) fcon(2)],'Upper',[Dcon(3) 10*data1(1) fcon(3)]);
        
        % fit model from http://doi.org/10.1002/mrm.22565, equation 2
        modelbiexpffree=fittype('a*((1-f)*(1-exp(-TR/T1))*exp(-x*D-TE/T2)+f*(1-exp(-TR/T1b))*exp(-x*(D+Dster)-TE/T2b))/((1-f)*exp(-TE/T2)*(1-exp(-TR/T1))+f*exp(-TE/T2b)*(1-exp(-TR/T1b)))','options',sbiexpf,'problem',{'Dster','TR', 'TE', 'T1', 'T2', 'T1b','T2b'});
        
        % fit
        [c21, gof]=fit(bvec,data1,modelbiexpffree,'problem',{Dstergeuss, TR, TE, T1, T2, T1b, T2b});
        
        D(k)=c21.D;
        f(k)=c21.f;
        S(k)=c21.a;
        SSE(k)=gof.sse;
        rsquare(k)=gof.adjrsquare;
    catch
        % in case an error occured, give negative number and add 1 to fail
        D(k)=-0.00001;
        f(k)=-0.00001;
        S(k)=-0.00001;
        fail=fail+1;
        SSE(k)=-0.00001;
        rsquare(k)=-0.00001;
    end
end

% when voxels are rejected (i.e. 0), then the fit will fail due to too many
% variables compared to data. This will tell you how often that occured.
sprintf('%d pixels failed due to too much rejection',fail)

end
