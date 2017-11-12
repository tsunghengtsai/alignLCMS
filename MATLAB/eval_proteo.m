%% Peak matching result by SIMA
% in_file = 'sima_mz01_rt75.txt';   % SIMA matching output

[idxPeak, names, idxIdv, massAll, rtAll] = textread(in_file,'%d %s %d %f %f');

nbIdvPeak = numel(names);
runs = zeros(nbIdvPeak,1);
for i = 1:numel(names)
    runs(i) = sscanf(names{i}, ['P1_' '%d' '.txt']);
end
idxPeak = idxPeak+1;
idxIdv = idxIdv+1;

nbPeak = idxPeak(end);
uniRuns = unique(runs);
nbRun = numel(uniRuns);

sizePeak = zeros(nbPeak,1);

idxDrop = [];
idxSinSize = [];
for i = 1:nbPeak
    tmpIdx = find(idxPeak==i);
    sizePeak(i) = numel(tmpIdx);
    if numel(tmpIdx)==1
        idxDrop = [idxDrop tmpIdx];
        idxSinSize = [idxSinSize i];
    end
end

idxPeak(idxDrop) = [];
runs(idxDrop) = [];
idxIdv(idxDrop) = [];
massAll(idxDrop) = [];
rtAll(idxDrop) = [];


%% load ground-truth data and evaluate the matching result
load ground_truth_proteomics.mat

for g = 1:numel(ground)
    tmpNb = numel(ground(g).idxSamp);
    match{g,1} = nan(1,tmpNb);
    for p = 1:tmpNb
        tmpIdx = find(runs==ground(g).idxSamp(p) & idxIdv==ground(g).idxPeak(p));
        if numel(tmpIdx)>1
            disp('error');
        elseif numel(tmpIdx)==1
            match{g,1}(p) = idxPeak(tmpIdx);
        end
    end
end

prec = 0;
rec = 0;
for g = 1:numel(ground)
    nonNan = match{g,1}(~isnan(match{g,1}));
    numPR = numel(find(~isnan(match{g,1})));
    pk = unique(nonNan);
    if numel(pk) ~= 0
        denP = sum(sizePeak(pk));
        denR = numel(pk)*numel(match{g,1});
        prec = prec + numPR/denP;
        rec = rec + numPR/denR;
    end
end

prec = prec/numel(ground);
rec = rec/numel(ground);
