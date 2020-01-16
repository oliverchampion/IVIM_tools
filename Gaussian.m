function [S,Madc,sigma,SSE,rsquare] = Gaussian(selection,bvec,allb)

%%%%%%%%
% This program calculates S, Madc, sigma, SSE and rsquare for each voxel in the
% input matrix "selection". The program uses the Gaussian aproach from
% the article. Fit model from http://doi.org/10.1002/mrm.10578
%
% Input:
% Selection = n*m matrix, where the first dimension includes the data from
% the different voxels and the second dimension is the data as function of
% bvalues.
%
% bvec = a vector containing the b-values of each data point. Must be m
% long
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
% Madc = apparent diffusion coefficient distribution maxima
%
% sigma = ADC distribution width
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
%%
%
% Code is written by Oliver Gurney-Champion
% o.j.gurney-chapion@amsterdamumc.nl
% 
%%
%%%%%%%%%

% fit constraints [initial guess, min, max]
Madccon=[0.0001, 0.00001, 0.01];
sigmacon=[0.0016, 0.0001, 0.1];


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
Madc=zeros(ssel(2),1);
sigma=zeros(ssel(2),1);
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
        % fitoptions
        gaus=fitoptions('Method','NonlinearLeastSquares','robust','on','Lower',[Madccon(2) 0.1*data1(1) sigmacon(2)],'Upper',[Madccon(3) 10*data1(1) sigmacon(3)],'Startpoint',[Madccon(1) data1(1) sigmacon(1)],'MaxIter',1000);
        % fit model from http://doi.org/10.1002/mrm.10578
        modelgaus=fittype('a*((1+erf(D/sigma/sqrt(2)-x*sigma/sqrt(2)))/(1+erf(D/sigma/sqrt(2)))*exp(-x*D+1/2*x^2*sigma^2))','options',gaus);
        % fit
        [c21, gof]=fit(bvec,data1,modelgaus);
        
        Madc(k)=c21.D;
        S(k)=c21.a;
        sigma(k)=c21.sigma;
        SSE(k)=gof.sse;
        rsquare(k)=gof.adjrsquare;
    catch
        % in case an error occured, give negative number and add 1 to fail
        Madc(k)=-0.00001;
        sigma(k)=-0.00001;
        S(k)=-0.00001;
        SSE(k)=-0.00001;
        rsquare(k)=-0.00001;
    end
end

% when voxels are rejected (i.e. 0), then the fit will fail due to too many
% variables compared to data. This will tell you how often that occured.

sprintf('%d pixels failed due to too much rejection',fail)

end
