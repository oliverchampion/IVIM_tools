% Bayesian stretched exponential fitting by Oliver Gurney-Champion based on Sebastiano Barbieri's BP IVIM fit
% If you found this software useful, please cite 
% Barbieri et al., Impact of the calculation algorithm on biexponential fitting 
% of diffusion-weighted MRI in upper abdominal organs,
% Magn Reson Med. 2016 May;75(5):2175-84. doi: 10.1002/mrm.25765. 

function [fitBi,falpha,xialpha] = fitStrExponentialBP(S,b,varargin)
%fitStrExponentialBP Perfom a stretched-exponential fit to the DWI data of one
% voxel according to the model 
% S(b) = S_0*exp(-(x*D)^alpha) . 
% Also computes the root-mean-square error for each fit.
%
% Bayesian Probability with Markov Chain Monte Carlo Integration is used
% for the fit
%
% Usage:
%   fitBi = fitBiExponentialSegmented(S,b,varargin)
%       S: vector of DWI values corresponding to the b values specified by
%              the b vector for one image voxel.
%       b: the vector of employed b values.
%       varargin: optional vector of initial parameter estimates [D,alpha]
%       fit.D: the fitted D value according to the stretched-exponential model.
%       fit.alpha: the fitted alpha value according to the stretched-exponential model.
%       fit.rmse: the root-mean-square error of the fit.
%       fit.time: the duration of the fit.
%       fit.failed: wether the fit produced an acceptable result or not.

S0 = S(1);

if S0 == 0
  fitBi.D = 0;
  fitBi.alpha = 0;
  fitBi.S0 = 0;
  fitBi.rmse = 0;
  fitBi.time = 0;
  fitBi.failed = 1;
else
  %compute fit
  tic;
  %magnitude estimates
  magD = 0.001;
  magalpha = 0.1;
  magP = [magD,magalpha];
  %priors
  priorD = @(D) (D > 0 & D < 0.006);
  prioralpha = @(alpha) (alpha > 0 & alpha < 2);
  %likelihood
  % p = [Dp,Dt,Fp]
  SNorm = S/S0;
  N = length(S);
  Q = @(p) sum( ((exp(-(b*p(1)).^p(2))) -SNorm).^2 );
  likelihood = @(p) (Q(p)/2).^(-N/2);
  %posterior
  post = @(p) likelihood(p)*priorD(p(1))*prioralpha(p(2));
  
 	try
    %MCMC
    nSamples = 3000;
    if ~isempty(varargin)
      initialP = [varargin{1}(1),varargin{1}(2)];
      samples = slicesample(initialP,nSamples,'pdf',post,'width',magP,'burnin',200);
    else
      samples = slicesample(magP,nSamples,'pdf',post,'width',magP,'burnin',200);
    end
    
    [fD,xiD] = ksdensity(samples(:,1),'npoints',1000);
    [falpha,xialpha] = ksdensity(samples(:,2),'npoints',1000);
    [~,idMax] = max(fD);
    tmp.D = xiD(idMax);
    [~,idMax] = max(falpha);
    tmp.alpha = xialpha(idMax);

    fitBi.D = tmp.D;
    fitBi.alpha = tmp.alpha;
    fitBi.S0 = S0;
    fitBi.time = toc;
    diff = S-fitBi.S0*(exp(-(b*fitBi.D).^fitBi.alpha));
    fitBi.rmse = sqrt(mean(diff.^2));
    if fitBi.D <= 0 || fitBi.D >= 0.006 || fitBi.alpha < 0 || fitBi.alpha > 2
      fitBi.failed = 1;
    else
      fitBi.failed = 0;
    end
  catch
    fitBi.D = 0;
    fitBi.alpha = 0;
    fitBi.S0 = 0;
    fitBi.rmse = 0;
    fitBi.time = 0;
    fitBi.failed = 1;
  end  
end

end
