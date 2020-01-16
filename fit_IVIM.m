function [S,f,D,Dstr,SSE,rsquare] = fit_IVIM(selection,bvec,TR,TE,T1,T2,T1b,T2b)
%%
% Slow, voxel-wise model fitting of IVIM using the T1/T2 componsated IVIM
% fit by Lemke et al: https://doi.org/10.1002/mrm.22565 
%
% Code is written by Oliver Gurney-Champion
% o.j.gurney-chapion@amsterdamumc.nl
% 
%%


ssel=size(selection);

S=zeros(ssel(2),1);
f=zeros(ssel(2),1);
D=zeros(ssel(2),1);
SSE=zeros(ssel(2),1);
rsquare=zeros(ssel(2),1);
Dstr=zeros(ssel(2),1);

for k=1:ssel(2)
    
        data1=selection(:,k);

                       
        sbiexpf=fitoptions('Method','NonlinearLeastSquares','robust','on','Startpoint',[0.001 0.08 1 0.05],'Lower',[0.0005 0.005 0 0.001],'Upper',[0.006 0.2 2 0.999]);
        
        modelbiexpffree=fittype('a*((1-f)*(1-exp(-TR/T1))*exp(-x*D-TE/T2)+f*(1-exp(-TR/T1b))*exp(-x*(D+Dster)-TE/T2b))/((1-f)*exp(-TE/T2)*(1-exp(-TR/T1))+f*exp(-TE/T2b)*(1-exp(-TR/T1b)))','options',sbiexpf,'problem',{'TR', 'TE', 'T1', 'T2', 'T1b','T2b'});
        
        [c21, gof]=fit(bvec,data1,modelbiexpffree,'problem',{TR, TE, T1, T2, T1b, T2b});
        
        D(k)=c21.D;
        f(k)=c21.f;
        S(k)=c21.a;
        Dstr(k)=c21.Dster;
        SSE(k)=gof.sse;
        rsquare(k)=gof.adjrsquare;

end


end
