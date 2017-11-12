function mcq = coda(data,winSize)

[nbRT,nbBin] = size(data);
hWinSize = (winSize-1)/2;

dataSmoothed = zeros(nbRT-winSize+1,nbBin);

for i=1:winSize
    dataSmoothed = dataSmoothed + data(i:(nbRT-winSize+i),:);
end

data([1:hWinSize, nbRT-hWinSize+1:nbRT],:) = [];
nbRT = nbRT-winSize+1;

normChrom = sqrt(sum(data.^2));
dataScaled = data./repmat(normChrom,nbRT,1);

meanSmoothed = mean(dataSmoothed);
stdSmoothed = std(dataSmoothed);

dataSmoothedScaled = (dataSmoothed-repmat(meanSmoothed,nbRT,1))./repmat(stdSmoothed,nbRT,1);

mcq = (1/sqrt(nbRT-1))*sum(dataScaled.*dataSmoothedScaled);
