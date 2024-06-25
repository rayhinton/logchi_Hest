close all; clear all; clc
%addpath('/Users/hvimalajeewa2/Documents/Tamu_Documents/TAMU/DemosNew')
%addpath('/Users/hvimalajeewa2/Documents/Tamu_Documents/TAMU/DistVar/supCodes/function/')
addpath('./MatlabFunctions/')

lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;
seed = 10;
rng(seed, "twister")

J = 13; n = 2^J;         % sample size
L = 1;

filt=[ -0.075765714789341  -0.029635527645954 ...
       0.497618667632458   0.803738751805216 ...
       0.297857795605542  -0.099219543576935 ...
      -0.012603967262261   0.032223100604071]; %Symmlet 4  

pairs = nchoosek(1 :J-1, 2); 
a = 3; b = 9;
pairs = pairs(find( pairs(:,1) >= a & pairs(:,2 ) <=b ),:);

% ismean >> 0 - median, 1 = mean
 ismean = 0;

H     = .20:.1:.8;%linspace(0.1,.7, 5);
nrep = 1;

H_est = zeros(length(H), nrep); H_est_old= zeros(length(H), nrep);

for j = 1%: length(H)
    for i = 1: nrep
        data = MakeFBMNew(2^J, H(j));
        
        wddata = dwtr(data, J - L, filt);
        
        [h_hat E] = MomentMatchHurst_new(wddata, pairs,L, ismean)
        %H_est(j, i) = h_hat;
    
        [slope, levels, log2spec] = waveletspectra_new(data, L, filt, a, b, 0,0)
        H_est_old(j, i) = log2spec
    
    end
end 