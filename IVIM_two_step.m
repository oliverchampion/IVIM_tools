function [f,D,Dster,SSE,rsquare] = IVIM_two_step(selection,bvec,cutoff,TR,TE)

%%%%%%%%
% This program calculates S, f, D, Dstr, SSE and rsquare for each voxel in the
% input matrix "selection". The program uses the IVIM-adaptive aproach from
% the article. The original work on which this implementation was built can
% be found here: http://doi.org/10.1002/jmri.22003
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
%
% Output (m long vectors, repressenting the voxel values):
% S = signal at b=0 s/mm2 (fitted)
%
% f = perfusion fraction (tissue fraction)
%
% D = diffusion coefficient
%
% Dstr = pseudo-diffusion
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

%% important parameter
%% The following parameter indicates the minimum of datapoints to be used. As in my data, low bvalues were acquiered 9-16 times, and the higher b-values only 4, the data sorting step will downscale to 1 mean acquisition per high b-value and 2-4 mean acquisitions for the lower bvalues. having minimal of 10 b-values effectively takes along the lowes 4 b-values (b=0 (4x), b=10 (2x), b=20 (2x) and b=3 (2x). This is the bare minimum to get a fit. However as you can see, this is very acquisition dependent and it might well be that for another acquisition this parameter can be decreased to 4.

if nargin==3
    TR=100000;
    TE=0;
end

%% This are tissue specific parameters. They are set for panctreas. May be changed for other organs.
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
Dstercon=[0.05, 0.006, 0.2];

%% here the data is sorted in blocks of b-values to speed the fitting. For the adaptice approach it is a must, as otherwise it steps over each acquisition.
fail=0;
a=unique(bvec);
weights=zeros(size(a,2),1);
for ii=1:size(weights,1)
    weights(ii)=sum(bvec==a(ii));
end
weights=round(weights/min(weights));
weights=ones(size(weights));

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
bvec=bvec2;
selection=selection2;
clear selection2 bvec2

%% initializing parameters
ssel=size(selection);
S=zeros(ssel(2),1);
f=zeros(ssel(2),1);
D=zeros(ssel(2),1);
Dster=zeros(ssel(2),1);
SSE=zeros(ssel(2),1);
rsquare=zeros(ssel(2),1);
bvecbu=transpose(bvec);


for k=1:ssel(2)
    try
        bvec=bvecbu;
        data1=selection(:,k);
        bvec(data1==0)=[];
        data1(data1==0)=[];
        lts=log(data1);
        %% here we loop over different combinations of b-values used as cutoff to determine the D and f before determining D*.
            %% first get D
            a=polyfit(bvec(bvec>cutoff),lts(bvec>cutoff),1);
            D_temp=-a(1);
            %% correcting for T1 and T2 effects of blood and tissue
            aa=(1-exp(-TR/T1))*exp(-TE/T2);
            bb=((1-exp(-TR/T1b))*exp(-TE/T2b));
            %% then use the b=0 measurement together with the intersept from the polynominal fit to find the perfusion fraction
            s0=data1(1);
            fff=(s0-exp(a(2)))/s0;
            %% applying correction for T1 and T2 effects
            f_temp=(fff*aa)/(bb+aa*fff-bb*fff);
            %% restrict the perfusion fraction to possitive values
            if f_temp<0
                f_temp=0;
            end
            
            
            %% then do the IVIM fit with D and f fixed to find D* and get a SSE parameter
            sbiexpD=fitoptions('Method','NonlinearLeastSquares','Lower',0.01,'Upper',4,'Startpoint',0.20,'MaxIter',1000);
            modelbiexpD=fittype('a*((1-f)*(1-exp(-TR/T1))*exp(-x*D-TE/T2)+f*(1-exp(-TR/T1b))*exp(-x*(Dster)-TE/T2b))/((1-f)*exp(-TE/T2)*(1-exp(-TR/T1))+f*exp(-TE/T2b)*(1-exp(-TR/T1b)))','options',sbiexpD,'problem',{'TR', 'TE', 'T1', 'T2', 'T1b','T2b','a','D','f'});
            
            [c21, geof]=fit(double(bvec),double(data1),modelbiexpD,'problem',{double(TR),double(TE), T1, T2, T1b, T2b,double(s0),D_temp,f_temp});
            %                     figure
            %                     semilogy(bvec,data1)
            %                     hold on
            %                     semilogy(bvec,s0*((1-f_temp(bs-1))*(1-exp(-TR/T1))*exp(-bvec.*D_temp(bs-1)-TE/T2)+f_temp(bs-1)*(1-exp(-TR/T1b))*exp(-bvec.*(Dstertemp(bs-1))-TE/T2b))/((1-f_temp(bs-1))*exp(-TE/T2)*(1-exp(-TR/T1))+f_temp(bs-1)*exp(-TE/T2b)*(1-exp(-TR/T1b))))
            Dstertemp=c21.Dster;
            sse=geof.sse;
            %           rsqr(bs-1)=geof.rsquare;
        %% after looping over all possible combinations of splitting the data, select the results from the cutoff value that gave the smallest SSE
        SSE(k)=geof.sse;
        rsquare(k)=geof.rsquare;
        f(k)=f_temp;
        D(k)=D_temp;
        Dster(k)=Dstertemp;
    catch
        %% if fails
        D(k)=-0.00001;
        f(k)=-0.00001;
        Dster(k)=-0.00001;
        fail=fail+1;
        SSE(k)=-0.00001;
        rsquare(k)=-0.00001;
    end
end

% when voxels are rejected (i.e. 0), then the fit will fail due to too many
% variables compared to data. This will tell you how often that occured.

sprintf('%d pixels failed due to too much rejection',fail)

end
