%% Description 
% Data set: proteomics 
% Method: single-profile alignment (base peak)
% 
% ----- bpc_proteomics.mat -----
% TICs: base peak chromatogram 20 x 1000 (sample x RT points)
% timeRes: registered time index


%% Required utilities
% ----- attached functions -----
% priorRatioLn() to calculate the prior ratio (in log) when GP prior is not available
% 
% ----- functions from elsewhere -----
% inv_posdef(), randnorm(), scale_rows(), ndsum() from Tom Minka's Lightspeed toolbox, 
%    downloaded at http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
% randraw() from File Exchange at MATLAB Central, 
%    downloaded at http://www.mathworks.com/matlabcentral/fileexchange/7309
% bsplinebasis() from Scott Gaffney's CCToolbox, 
%    downloaded at http://www.ics.uci.edu/~sgaffney/software/CCT/ 


%% Load the proteomic data set 
clear all; close all; clc;

load bpc_proteomics.mat


%% Generate the chromatograms to be aligned
[nbRT,nbSamp] = size(TICs);
TICs = TICs./(2*10^7);  % scale to a range of [0,10]
timeGrid = (1:nbRT)';


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

for i = 1:nbSamp
    map(i).coeff = varKnot;
    map(i).acpt = zeros(4,nbRT-1);
end


%% hyperparameters
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
hyStd = 0.1;
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
spReg = zeros(nbReg,nbMCMC);    % regression coeff (prototype function)
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
spReg(:,1) = muReg;
spMap(:,:,1) = repmat(varKnot,1,nbSamp);

% Metropolis step
stepMH1 = 3;
stepMH2 = 10;

%% Run MCMC
BSTilt = repmat(BSReg,nbSamp,1); % BS_i (space allocation)
vecScale = ones(nbSamp*nbRT,1);
vecShift = ones(nbSamp*nbRT,1);
vecTIC = TICs(:);

rng default     % reset the random seed
tic             % initialize timer

for mc = 2:nbMCMC
    % Matrix manipulation
    idxMat = interp1(varKnot, spMap(:,:,mc-1), timeGrid);
    BSTilt = interp1((1:lenBS)',BSReg,idxMat(:));
    repScale = repmat(spScale(:,mc-1)',nbRT,1);
    vecScale = repScale(:);
    SBSTilt = scale_rows(BSTilt,vecScale);  % a_i*BS_i from lightspeed
    repShift = repmat(spShift(:,mc-1)',nbRT,1);
    vecShift = repShift(:);
    %% Gibbs sampling goes below
    % regression coefficients of prototype function
    invCovReg = SigmaReg*spTauPsi(mc-1);
    tmpCov = inv_posdef(invCovReg + (SBSTilt'*SBSTilt)*spTauEpsilon(mc-1)); % from lightspeed
    tmpMuVec = tmpCov*SBSTilt'*(vecTIC-vecShift)*spTauEpsilon(mc-1);
    spReg(:,mc) = randnorm(1,tmpMuVec,[],tmpCov);   % from lightspeed
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
        tmpMat = [BSTilt((i-1)*nbRT+1:i*nbRT,:)*spReg(:,mc) ones(nbRT,1)];
        tmpCov = inv_posdef(diag([spTauScale(mc-1) spTauShift(mc-1)]) + spTauEpsilon(mc-1)*tmpMat'*tmpMat); % from lightspeed
        tmpMuVec = tmpCov*(diag([spTauScale(mc-1) spTauShift(mc-1)])*[spMuScale(mc);spMuShift(mc)]...
            + spTauEpsilon(mc-1)*tmpMat'*TICs(:,i));
        tmpSp = randnorm(1,tmpMuVec,[],tmpCov);     % from lightspeed
        spScale(i,mc) = tmpSp(1);
        spShift(i,mc) = tmpSp(2);
    end
    % 1/var(ei)
    tmpShape = hyShapeEpsilon + 0.5*nbRT*nbSamp;
    repScale = repmat(spScale(:,mc)',nbRT,1);
    vecScale = repScale(:);
    SBSTilt = scale_rows(BSTilt,vecScale);  % from lightspeed
    repShift = repmat(spShift(:,mc)',nbRT,1);
    vecShift = repShift(:);
    vecTICHat = SBSTilt*spReg(:,mc) + vecShift;
    tmpRate = hyRateEpsilon + 0.5*(vecTIC-vecTICHat)'*(vecTIC-vecTICHat);
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
    tmpShape = hyShapePsi + 0.5*nbReg;
    tmpRate = hyRatePsi + 0.5*spReg(:,mc)'*SigmaReg*spReg(:,mc);
    spTauPsi(mc) = gamrnd(tmpShape,1/tmpRate);
    %% Metropolis-Hastings algo
    for i = 1:nbSamp
        tmpMap = map(i).coeff;
        tmpIdx = interp1(varKnot,tmpMap,timeGrid);
        tmpTIC = spScale(i,mc)*interp1(timeExt,BSReg,tmpIdx)*spReg(:,mc) + spShift(i,mc);
        tmpEvaln = -0.5*spTauEpsilon(mc)*sum((TICs(:,i)-tmpTIC).^2);
        % generate blocks 
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
            pRatioln = priorRatioLn(tmpMap, varKnot, blkMH(1,m):blkMH(2,m), stepProp, hyStd);
            tmpIdx = interp1(varKnot,tmpMapProp,timeGrid);
            tmpTIC = spScale(i,mc)*interp1(timeExt,BSReg,tmpIdx)*spReg(:,mc) + spShift(i,mc);
            tmpEvaPropln = -0.5*spTauEpsilon(mc)*sum((TICs(:,i)-tmpTIC).^2);
            lRatioln = tmpEvaPropln-tmpEvaln;
            switch rStep
                case 1
                    map(i).acpt(2,varKnot(blkMH(1,m):blkMH(2,m))-1) = map(i).acpt(2,varKnot(blkMH(1,m):blkMH(2,m))-1) + 1;
                case 2
                    map(i).acpt(4,varKnot(blkMH(1,m):blkMH(2,m))-1) = map(i).acpt(4,varKnot(blkMH(1,m):blkMH(2,m))-1) + 1;
            end
            if rand(1) < min([1, exp(lRatioln+pRatioln)])
                tmpMap = tmpMapProp;
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
    in_text = ['sima\P1_' num2str(i) '.txt'];
    out_text = ['sp\P1_' num2str(i) '.txt'];
    tmp = textread(in_text);
    newRT = interp1(timeRes,crtRT(:,i),tmp(:,4),'pchip');
    tmp(:,4) = round(newRT*100)/100;
    dlmwrite(out_text, tmp, 'delimiter', '\t', 'precision', 10, 'newline', 'pc');
end
