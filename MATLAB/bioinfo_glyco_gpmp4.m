%% Description 
% Data set: glycomics
% Method: multi-profile alignment with Gaussian process prior
% 
% ----- data_matrix_glycomics.mat -----
% glyData: binned matrix 23 x 1000 x 3000 (sample x RT points x mz bins)
% timeRes: registered time index
% nbSamp: # of samples (23)
% nbRT: # of RT points (1000)
% nbMZ: # of mz bins (3000)
% mzLow: mz lower bound for each bin
% mzHi: mz upper bound for each bin


%% Required utilities
% ----- attached functions -----
% coda() to calculate MCQ values
% 
% ----- functions from elsewhere -----
% inv_posdef(), randnorm(), scale_rows(), ndsum() from Tom Minka's Lightspeed toolbox, 
%    downloaded at http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
% randraw() from File Exchange at MATLAB Central, 
%    downloaded at http://www.mathworks.com/matlabcentral/fileexchange/7309
% bsplinebasis() from Scott Gaffney's CCToolbox, 
%    downloaded at http://www.ics.uci.edu/~sgaffney/software/CCT/ 
% apcluster() by Frey Lab,
%    downloaded at http://www.psi.toronto.edu/index.php?q=affinity%20propagation
% GPML toolbox (v3.1) by Carl Edward Rasmussen and Hannes Nickisch
%    downloaded at http://www.gaussianprocess.org/gpml/code/matlab/doc/index.html


%% Load the glycomic data set 
clear all; close all; clc;

load data_matrix_glycomics.mat

[nbSamp,nbRT,nbMZ] = size(glyData);


%% 1st-phase screening based on MCQ 
binICs = glyData;
winQ = 3;
% considering the worst case 
mcq = coda(squeeze(binICs(1,:,:)),winQ);
for i = 2:nbSamp
    mcq = min([mcq; coda(squeeze(binICs(i,:,:)),winQ)],[],1);
end

idxQ = find(mcq>=0.9); 
binQICs = binICs(:,:,idxQ);
mcqQ = mcq(idxQ);
nbQBin = numel(idxQ);

% scale the remaining chromatograms (profile is more important than abosolute amount)
for b = 1:nbQBin
    for i = 1:nbSamp
        binQICs(i,:,b) = binQICs(i,:,b)./ndsum(binQICs(i,:,b),2);
    end
end


%% 2nd-phase screening based on reproducibility among samples (with xcorr) 
xcMat = zeros(nbSamp,nbSamp,nbQBin); 
xcAve = zeros(1,nbQBin);
rootE = zeros(1,nbSamp);   % root of energy used for normalization
winX = 100;

for b = 1:nbQBin
    for i = 1:nbSamp
        rootE(i) = sqrt(sum((binQICs(i,:,b).*binQICs(i,:,b))));
        xcMat(i,i,b) = 1;
    end
    for i = 1:nbSamp 
        for j = i+1:nbSamp
            tmpXC = xcorr(binQICs(i,:,b),binQICs(j,:,b),winX); 
            idxXC = maxind(tmpXC);
            if isempty(idxXC)
                xcMat(i,j,b) = 0;
            else
                idxXC(find(tmpXC(idxXC)<0.5*max(tmpXC(idxXC)))) = [];
                [~,idxx] = min(abs(idxXC-(winX+1)));
                xcMat(i,j,b) = tmpXC(idxXC(idxx))/(rootE(i)*rootE(j)); 
            end
            xcMat(j,i,b) = xcMat(i,j,b);
        end
    end
    xcAve(b) = (ndsum(xcMat(:,:,b),1:2)-nbSamp)/(nbSamp*nbSamp-nbSamp);
end

idxC = find(xcAve>=0.85); 
binQCICs = binQICs(:,:,idxC);
mcqQC = mcqQ(idxC);
nbQCBin = numel(idxC);
xcAveQC = xcAve(idxC);


%% Identify exemplars using affinity propagation (correlation coeff. as similarity) 
ccQC = zeros(nbQCBin);
for b = 1:nbSamp
    ccQC = ccQC + corr(squeeze(binQCICs(b,:,:))); 
end
ccQC = ccQC/nbSamp;

sim = zeros(nbQCBin*nbQCBin-nbQCBin,3);
cnt = 1;
for b = 1:nbQCBin
    for p = [1:b-1,b+1:nbQCBin]
        sim(cnt,1)=b; sim(cnt,2)=p; sim(cnt,3)=ccQC(b,p);
        cnt=cnt+1;
    end
end
prefSim = mean(sim(:,3)); 
[idxExmp,~,~,~] = apcluster(sim,prefSim);
nbExmp = numel(unique(idxExmp)); 
exmpICs = binQCICs(:,:,unique(idxExmp)); 


%% Agglomerative clustering of the exemplars (based on overlapping level) 
tmpNum = nbExmp;
tmpICs = exmpICs;
tmpIdx = cell(1,tmpNum);
for b = 1:tmpNum
    tmpIdx{b} = b;
end

for l = 1:nbExmp-1
    tmpDis = zeros(tmpNum*(tmpNum-1)/2,3); 
    cnt = 1;
    for b = 1:tmpNum-1
        for p = b+1:tmpNum
            tmpDis(cnt,1) = b;
            tmpDis(cnt,2) = p;
            for i = 1:nbSamp
                tmpDis(cnt,3) = tmpDis(cnt,3) + sum(min(squeeze(tmpICs(i,:,[b,p])),[],2));
            end
            cnt = cnt+1;
        end
    end
    [layer(l).ovp, idxPair] = min(tmpDis(:,3));
    tmpICs(:,:,tmpDis(idxPair,1)) = tmpICs(:,:,tmpDis(idxPair,1))+tmpICs(:,:,tmpDis(idxPair,2));
    tmpICs(:,:,tmpDis(idxPair,2)) = [];
    tmpIdx{tmpDis(idxPair,1)} = [tmpIdx{tmpDis(idxPair,1)},tmpIdx{tmpDis(idxPair,2)}];
    tmpIdx(tmpDis(idxPair,2)) = [];
    tmpNum = tmpNum-1;
    layer(l).num = tmpNum;
    layer(l).idx = tmpIdx;
end

nbEIC = 4;  % pre-defined value
EICs = zeros(nbSamp,nbRT,nbEIC);
for b=1:nbEIC
    EICs(:,:,b) = ndsum(exmpICs(:,:,layer(numel(layer)-nbEIC+1).idx{b}),3);
end
EICs = EICs*10;     % scale to a range of [0,10]

EICs = permute(EICs,[2 1 3]);
timeGrid = (1:nbRT)';

clear chrom glyData agg aggPair bb binEdge binICs binQCICs binQICs ccQC cnt ... 
    idxAgg idxC idxExmp idxQ idxXC idxx mcq mcqQ mcqQC ... 
    mzHi mzLow nbExmp nbMZ nbBin nbQBin nbQCBin p prefAgg prefSim rootE ... 
    sim tmpXC winQ winX xcAve xcAveQC xcMat xcMin xcMinQC i j 


%% Compile internal standard reference
timeIS = textread('time_standard_glycomics.txt');

timeGrid = (1:nbRT)';
timeGridScaled = timeGrid./nbRT;

timeScaled = (timeIS-10)./50;   % scale RT [10,60] to [0,1]
timeRef = nanmean(timeScaled,2);   % reference time for the internal standard

meanMapTime = zeros(nbRT,nbSamp);
stdMapTime = zeros(nbRT,nbSamp);


%% GP hyperparameters set 
covfunc = @covSEiso; hyp.cov = log([0.05; 0.1]); 
likfunc = @likGauss; hyp.lik = log(0.05); 
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [1; 0]; 

for i = 1:nbSamp
    idxTime = find(~isnan(timeScaled(:,i)));    % ignore absent peaks
    [meanMapTime(:,i), stdMapTime(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, ...
        timeScaled(idxTime,i), timeRef(idxTime), timeGridScaled);
end

meanMapTime = meanMapTime.*nbRT;
stdMapTime = stdMapTime.*nbRT;
varMapTime = stdMapTime.^2;


%% B-spline
order = 3;
timeExt = timeGrid;

denKnotReg = 0.5;  % density of knots for prototype function (0.25--0.75)
denKnotMap = 0.025;  % density of knots for mapping function (<= 0.1)
ptStart = timeExt(1);
ptEnd = timeExt(end);
lenBS = length(timeExt);
nbKnotReg = ceil(lenBS*denKnotReg);
knotsReg = unique(linspace(ptStart,ptEnd,nbKnotReg));
knotsReg = [knotsReg(1)*ones(1,order) knotsReg(2:(end-1)) ...
    knotsReg(end)*ones(1,order)];
nbReg = length(knotsReg)-order;
BSReg = bsplinebasis(knotsReg,order,timeExt);
muReg = zeros(nbReg,1);

nbMap = ceil((timeGrid(end)-timeGrid(1))*denKnotMap);
varKnot = unique(round(linspace(timeGrid(1),timeGrid(end),nbMap)))';
meanMap = meanMapTime(varKnot,:);  % mean of mapping function coeff by GP
meanMap(1,:) = varKnot(1);          % fixed initial point
meanMap(end,:) = varKnot(end);      % fixed ending point

for i = 1:nbSamp
    map(i).coeff = meanMap(:,i);
    map(i).acpt = zeros(4,nbRT-1);
end


%% Hyperparameters
hyMuScale = 1;
hyMuShift = 0;
hyTauScale = 1/0.5;
hyTauShift = 1/0.5;
hyShapeScale = 0.1;
hyRateScale = 1;
hyShapeShift = 0.1;
hyRateShift = 1;
hyShapePsi = 0.1;
hyRatePsi = 1;
hyShapeEpsilon = 0.1;
hyRateEpsilon = 0.2;
SigmaReg = diag([2*ones(1,nbReg-1) 1],0) + diag(-1*ones(1,nbReg-1),1) ...
    +diag(-1*ones(1,nbReg-1),-1);


%% MCMC setting/initialization
nbMCMC = 15000;

% Space allocation MCMC runs
spMuScale = zeros(1,nbMCMC);    % a0
spMuShift = zeros(1,nbMCMC);    % c0
spScale = zeros(nbSamp,nbMCMC);  % ai
spShift = zeros(nbSamp,nbMCMC);  % ci
spTauScale = zeros(1,nbMCMC);   % 1/var(ai)
spTauShift = zeros(1,nbMCMC);   % 1/var(ci)
spTauEpsilon = zeros(1,nbMCMC); % 1/var(ei)
spTauPsi = zeros(1,nbMCMC);     % 1/var for regression coeff
spReg = zeros(nbReg,nbEIC,nbMCMC);    % regression coeff (prototype function)
spMap = zeros(nbMap,nbSamp,nbMCMC);% mapping function coeff

% Initial value assignment
spMuScale(1) = hyMuScale;
spMuShift(1) = hyMuShift;
spScale(:,1) = hyMuScale*ones(nbSamp,1);
spShift(:,1) = hyMuShift*ones(nbSamp,1);
spTauScale(1) = hyShapeScale/hyRateScale;
spTauShift(1) = hyShapeShift/hyRateShift;
spTauEpsilon(1) = hyShapeEpsilon/hyRateEpsilon;
spTauPsi(1) = hyShapePsi/hyRatePsi;
spMap(:,:,1) = meanMap;     % mean of mapping function coeff by GP

% Metropolis step
stepMH1 = 3;
stepMH2 = 10;


%% Run MCMC
BSTilt = repmat(BSReg,nbSamp,1); % BS_i (space allocation)
vecScale = ones(nbSamp*nbRT,1);
vecShift = ones(nbSamp*nbRT,1);
vecEICs = zeros(nbSamp*nbRT,nbEIC);
for b = 1:nbEIC
    vecEICs(:,b) = reshape(EICs(:,:,b),nbSamp*nbRT,1);
end

rng default     % reset the random seed
tic             % initialize timer

for mc = 2:nbMCMC
    % Matrix manipulation
    idxMat = interp1(varKnot, spMap(:,:,mc-1), timeGrid);
    BSTilt = interp1(timeExt,BSReg,idxMat(:));
    repScale = repmat(spScale(:,mc-1)',nbRT,1);
    vecScale = repScale(:);
    SBSTilt = scale_rows(BSTilt,vecScale);  % a_i*BS_i from lightspeed
    repShift = repmat(spShift(:,mc-1)',nbRT,1);
    vecShift = repShift(:);
    %% Gibbs sampling goes below
    % regression coefficients of prototype function
    invCovReg = SigmaReg*spTauPsi(mc-1);
    tmpCov = inv_posdef(invCovReg + (SBSTilt'*SBSTilt)*spTauEpsilon(mc-1)); % from lightspeed
    for b = 1:nbEIC
        tmpMuVec = tmpCov*SBSTilt'*(vecEICs(:,b)-vecShift)*spTauEpsilon(mc-1);
        spReg(:,b,mc) = randnorm(1,tmpMuVec,[],tmpCov);   % from lightspeed
    end
    % a0
    tmpVar = (hyTauScale + nbSamp*spTauScale(mc-1))^(-1);
    tmpMu = tmpVar*(hyMuScale*hyTauScale + sum(spScale(:,mc-1))*spTauScale(mc-1));
    spMuScale(mc) = tmpMu + sqrt(tmpVar)*randn(1);
    % c0
    tmpVar = (hyTauShift + nbSamp*spTauShift(mc-1))^(-1);
    tmpMu = tmpVar*(hyMuShift*hyTauShift + sum(spShift(:,mc-1))*spTauShift(mc-1));
    spMuShift(mc) = tmpMu + sqrt(tmpVar)*randn(1);
    % (ai, ci)
    for i = 1:nbSamp
        tmpMat = [reshape(BSTilt((i-1)*nbRT+1:i*nbRT,:)*spReg(:,:,mc),nbRT*nbEIC,1), ...
            ones(nbRT*nbEIC,1)];
        tmpCov = inv_posdef(diag([spTauScale(mc-1) spTauShift(mc-1)]) + spTauEpsilon(mc-1)*tmpMat'*tmpMat); % from lightspeed
        tmpMuVec = tmpCov*(diag([spTauScale(mc-1) spTauShift(mc-1)])*[spMuScale(mc);spMuShift(mc)]...
            + spTauEpsilon(mc-1)*tmpMat'*reshape(EICs(:,i,:),nbRT*nbEIC,1));
        tmpSp = randnorm(1,tmpMuVec,[],tmpCov);     % from lightspeed
        spScale(i,mc) = tmpSp(1);
        spShift(i,mc) = tmpSp(2);
    end
    % 1/var(ei)
    tmpShape = hyShapeEpsilon + 0.5*nbRT*nbSamp*nbEIC;
    repScale = repmat(spScale(:,mc)',nbRT,1);
    vecScale = repScale(:);
    SBSTilt = scale_rows(BSTilt,vecScale);  % from lightspeed
    repShift = repmat(spShift(:,mc)',nbRT,1);
    vecShift = repShift(:);
    vecTICHat = reshape(SBSTilt*spReg(:,:,mc),nbRT*nbSamp*nbEIC,1) + repmat(vecShift,nbEIC,1);
    tmpRate = hyRateEpsilon + 0.5*(vecEICs(:)-vecTICHat)'*(vecEICs(:)-vecTICHat);
    spTauEpsilon(mc) = gamrnd(tmpShape,1/tmpRate);
    % 1/var(ai)
    tmpShape = hyShapeScale + 0.5*nbSamp;
    tmpRate = hyRateScale + 0.5*sum((spScale(:,mc)-spMuScale(mc)).^2);
    spTauScale(mc) = gamrnd(tmpShape,1/tmpRate);
    % 1/var(ci)
    tmpShape = hyShapeShift + 0.5*nbSamp;
    tmpRate = hyRateShift + 0.5*sum((spShift(:,mc)-spMuShift(mc)).^2);
    spTauShift(mc) = gamrnd(tmpShape,1/tmpRate);
    % 1/var for the prototype function
    tmpShape = hyShapePsi + 0.5*nbReg*nbEIC;
    tmpMat1 = spReg(:,1,mc)*spReg(:,1,mc)';
    for b=2:nbEIC
        tmpMat1 = tmpMat1 + spReg(:,b,mc)*spReg(:,b,mc)';
    end
    tmpRate = hyRatePsi + 0.5*trace(tmpMat1*SigmaReg);
    spTauPsi(mc) = gamrnd(tmpShape,1/tmpRate);
    %% Metropolis-Hastings algo
    for i = 1:nbSamp
        tmpMap = map(i).coeff;
        tmpIdx = interp1(varKnot,tmpMap,timeGrid);
        tmpBSReg = interp1(timeExt,BSReg,tmpIdx);
        for b = 1:nbEIC
            tmpEICs(:,b) = spScale(i,mc)*tmpBSReg*spReg(:,b,mc) + spShift(i,mc);
        end
        tmpEvaln = -0.5*spTauEpsilon(mc)* ndsum((squeeze(EICs(:,i,:))-tmpEICs).^2, [1 2]);
        % generating blocks 
        rInd = randi(3);
        switch rInd
            case 1
                rBound = 1;
            case 2
                rBound = 0.5;
            case 3
                rBound = 0.25;
        end
        idxBound = find(rand(1,nbMap-3)<rBound)+2;  % block ends at nbMap-1
        blkMH = [2,idxBound; idxBound-1,nbMap-1];
        nbBlock = size(blkMH,2);
        rStep = randi(2);
        switch rStep
            case 1
                stepMH = stepMH1;
            case 2
                stepMH = stepMH2;
        end
        if mc <= 200
            stepMH = 30;    % propose big move in early MCMC iterations
        end
        for m = 1:nbBlock
            tmpMapProp = tmpMap;
            % identify moveable range
            tmpLB = tmpMap(blkMH(1,m)-1)-tmpMap(blkMH(1,m));
            tmpUB = tmpMap(blkMH(2,m)+1)-tmpMap(blkMH(2,m));
            % uniform proposal reflective on the boundary
            unBound = true;
            stepProp = stepMH*(2*rand(1)-1);
            while unBound
                if stepProp > tmpUB
                    stepProp = 2*tmpUB-stepProp;
                elseif stepProp < tmpLB
                    stepProp = 2*tmpLB-stepProp;
                else
                    unBound = false;
                end
            end
            tmpMapProp(blkMH(1,m):blkMH(2,m)) = tmpMapProp(blkMH(1,m):blkMH(2,m)) + stepProp;
            tmpIdxProp = interp1(varKnot,tmpMapProp,timeGrid);
            idxEval = varKnot(blkMH(1,m)-1):varKnot(blkMH(2,m)+1);  % range of interest
            pRatioln = -0.5*sum( ((tmpIdxProp(idxEval)-meanMapTime(idxEval,i)).^2 ...
               - (tmpIdx(idxEval)-meanMapTime(idxEval,i)).^2) ...
               ./ varMapTime(idxEval,i) );    % log of prior odds
            tmpBSReg = interp1(timeExt,BSReg,tmpIdxProp);
            for b = 1:nbEIC
                tmpEICs(:,b) = spScale(i,mc)*tmpBSReg*spReg(:,b,mc) + spShift(i,mc);
            end
            tmpEvaPropln = -0.5*spTauEpsilon(mc)* ndsum((squeeze(EICs(:,i,:))-tmpEICs).^2, [1 2]);
            lRatioln = tmpEvaPropln-tmpEvaln;
            switch rStep
                case 1
                    map(i).acpt(2,varKnot(blkMH(1,m):blkMH(2,m))-1) = map(i).acpt(2,varKnot(blkMH(1,m):blkMH(2,m))-1) + 1;
                case 2
                    map(i).acpt(4,varKnot(blkMH(1,m):blkMH(2,m))-1) = map(i).acpt(4,varKnot(blkMH(1,m):blkMH(2,m))-1) + 1;
            end
            if rand(1) < min([1, exp(lRatioln+pRatioln)])
                tmpMap = tmpMapProp;
                tmpIdx = tmpIdxProp;
                tmpEvaln = tmpEvaPropln;
                switch rStep
                    case 1
                        map(i).acpt(1,varKnot(blkMH(1,m):blkMH(2,m))-1) = map(i).acpt(1,varKnot(blkMH(1,m):blkMH(2,m))-1) + 1;
                    case 2
                        map(i).acpt(3,varKnot(blkMH(1,m):blkMH(2,m))-1) = map(i).acpt(3,varKnot(blkMH(1,m):blkMH(2,m))-1) + 1;
                end
            end
        end
        map(i).coeff = tmpMap;
        spMap(:,i,mc) = tmpMap;
    end
    if rem(mc,500)==0
        fprintf('Iteration %d, time %d   \r', mc,toc/60);
    end
end


%% Retention time correction
postKnot = ndsum(spMap(:,:,5001:15000),3)/10000;    % initial 5000 samples as burn-in
postMap = interp1(varKnot, postKnot, timeGrid);

crtRT = zeros(nbRT,nbSamp);
for i = 1:nbSamp
    crtRT(:,i) = interp1(timeGrid,timeRes,postMap(:,i),'pchip');
end

for i = 1:nbSamp
    in_text = ['sima\G1_' num2str(i) '.txt'];
    out_text = ['gpmp4\G1_' num2str(i) '.txt'];
    tmp = textread(in_text);
    newRT = interp1(timeRes,crtRT(:,i),tmp(:,4),'pchip');
    tmp(:,4) = round(newRT*100)/100;
    dlmwrite(out_text, tmp, 'delimiter', '\t', 'precision', 10, 'newline', 'pc');
end
