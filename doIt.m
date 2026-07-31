clear all
% close all

%%%%%%%%%%%%%%%%%%%%%
%% Set up environment
%%%%%%%%%%%%%%%%%%%%%
workDir = fileparts(mfilename('fullpath'));
%% %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ================ CONTROL PANEL ================
% Age-structured (Leslie-matrix) population model. Each yearly age class
% ages up one bin per year, subject to age-specific mortality; births
% enter age 0 in proportion to age-specific fertility. Edit and hit Run.
ageMax     = 150;      % y, oldest age class (plus-group: survivors accumulate here)
steepAgeFrac        = 0.95;   % fraction of ageMax past which mortality accelerates sharply (senescence cliff)
steepMortalityScale = 0.02;   % extra-hazard scale added past steepAgeFrac*ageMax
steepMortalityRate  = 1;      % extra-hazard growth rate (per year) past steepAgeFrac*ageMax; higher = sharper cliff
nYears     = 150;     % y, simulation horizon
popInit    = 8e9;   % total starting population, distributed across ages per the real-world
                       % 2021-2023 world age structure below (itself NOT the model's own stable age
                       % distribution, so the model's convergence to its own shape over time is visible)

% Mortality: per-capita annual death rate, AGE-DEPENDENT, rising with age (Gompertz-like).
deathRateBase  = 0.002;   % baseline hazard at age 0
deathRateSlope = 0.01;    % additional hazard at age ageMax (hazard = base + slope*(age/ageMax)^2)

% Fertility: per-capita annual fertility rate, AGE-DEPENDENT, triangular over a fertile age window
% (zero outside it, so births only come from ages in [fertileMin fertileMax]).
fertileMin       = 15;    % y, youngest fertile age
fertileMax       = 45;    % y, oldest fertile age
fertilityPeakAge = 28;    % y, age of peak fertility

% Net reproduction rate (R0) target -- expected offspring per person over
% a lifetime, given the mortality schedule above. Fertility is auto-scaled
% to hit this target, which is what makes the population's LONG-RUN total
% steady rather than exponentially growing/shrinking. R0=1 -> steady
% state; R0>1 -> long-run growth; R0<1 -> long-run decline.
targetR0 = 1.0;
%% =================================================
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%% Run simulation
%%%%%%%%%%%%%%%%%%%
ageVec = (0:ageMax)';
nAge   = numel(ageVec);

%%% age-specific rates (both birth and death vary with age, not constant across the population)
mortality = deathRateBase + deathRateSlope * (ageVec/ageMax).^2;   % per-capita annual hazard, increases with age
steepAgeThreshold = steepAgeFrac * ageMax;
pastSteepAge = ageVec > steepAgeThreshold;
mortality(pastSteepAge) = mortality(pastSteepAge) + steepMortalityScale * (exp(steepMortalityRate*(ageVec(pastSteepAge)-steepAgeThreshold)) - 1);   % sharp senescence cliff past steepAgeFrac*ageMax
survival  = exp(-mortality);                                       % per-capita annual survival probability

fertileAges = ageVec >= fertileMin & ageVec <= fertileMax;
fertilityShape = zeros(nAge,1);
halfWidth = max(fertilityPeakAge-fertileMin, fertileMax-fertilityPeakAge);
fertilityShape(fertileAges) = max(0, 1 - abs(ageVec(fertileAges)-fertilityPeakAge)/halfWidth);   % triangular, zero outside fertile window

% scale fertility so the net reproduction rate matches targetR0: R0 = sum
% over ages of (probability of surviving birth-to-age) x (fertility at age).
survivorship = cumprod([1; survival(1:end-1)]);   % l(age) = P(alive at age | born)
R0_unscaled  = sum(survivorship .* fertilityShape);
fertility    = fertilityShape * (targetR0 / R0_unscaled);   % age-dependent fertility rate

%%% initial age distribution: current real-world world age structure (CIA World Factbook,
%%% 2021-2023 estimates -- 0-14: 25.2%, 15-24: 15.3%, 25-54: 40.6%, 55-64: 9.2%, 65+: 9.7%),
%%% assigned uniformly across the years within each published bin. The open-ended 65+ bin has no
%%% published within-bin breakdown, so its share is tapered across ages 65:ageMax by the model's OWN
%%% survivorship curve (renormalized to that bin's real total) rather than spread flat all the way to
%%% ageMax, which would be unrealistic.
worldAgeBins  = [0 14; 15 24; 25 54; 55 64; 65 ageMax];   % y, [lo hi] inclusive
worldAgeShare = [0.252 0.153 0.406 0.092 0.097];          % fraction of world population in each bin
N0 = zeros(nAge,1);
for b = 1:size(worldAgeBins,1)
    inBin = ageVec >= worldAgeBins(b,1) & ageVec <= worldAgeBins(b,2);
    if b < size(worldAgeBins,1)
        N0(inBin) = worldAgeShare(b) / nnz(inBin);          % uniform within bin
    else
        w = survivorship(inBin); w = w / sum(w);            % taper the open-ended 65+ bin by the model's own survivorship shape
        N0(inBin) = worldAgeShare(b) * w;
    end
end

%%% cohort simulation (Leslie matrix, applied step by step)
N = zeros(nAge, nYears+1);
N(:,1) = popInit * N0;   % initial age distribution: real-world world age structure (see above), not uniform

for t = 1:nYears
    survivors = N(1:end-1,t) .* survival(1:end-1);   % age up one year, age-dependent survival
    N(2:end,t+1) = survivors;
    N(end,t+1)   = N(end,t+1) + N(end,t)*survival(end);   % plus-group: oldest bin retains its own survivors
    N(1,t+1)     = sum(N(:,t) .* fertility);              % births enter age 0, age-dependent fertility
end

totalPop = sum(N,1);
%% %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%% Plot results
%%%%%%%%%%%%%%%%%%%
years = 0:nYears;
fig = figure('Position',[100 100 1000 800],'Color','k'); ht = tiledlayout(2,2); ht.Padding = 'compact'; ht.TileSpacing = 'compact';
title(ht, sprintf('age-structured population model (R0 target = %.2f)', targetR0), 'Color','w')

nexttile
plot(ageVec, mortality, 'r-', ageVec, fertility, 'g-', 'LineWidth', 1.5)
xlabel('age'); ylabel('per-capita annual rate')
title('age-dependent rates')
legend({'death rate','fertility rate'}, 'Location','best')
grid on
mortalityCapNoCliff = deathRateBase + deathRateSlope * (ageMax/ageMax)^2;   % baseline mortality at ageMax with steepMortalityScale=0 (no cliff); caps the axis so an active cliff spike doesn't squish the rest of the plot
yl = ylim; yl(2) = max(max(fertility), mortalityCapNoCliff); ylim(yl);

nexttile
plot(years, totalPop, 'w-', 'LineWidth', 1.5)
xlabel('year'); ylabel('total population')
title('total population over time')
grid on
yl = ylim; ylim([0 yl(2)])

nexttile
hold on
plot(ageVec, N(:,1)/sum(N(:,1)), 'c-', 'LineWidth', 1.5)
plot(ageVec, N(:,end)/sum(N(:,end)), 'y-', 'LineWidth', 1.5)
xlabel('age'); ylabel('fraction of population')
title('age distribution: initial vs. final')
legend({'initial (world)','final (stable)'}, 'Location','best')
grid on

ax = nexttile;
Nnorm = N ./ sum(N,1);   % normalize each year (column) to sum to 1 -- per-year age-distribution shape
imagesc(years, ageVec, Nnorm); set(gca,'YDir','normal')
xlabel('year'); ylabel('age')
title('age-stratified population (normalized per year)')
colorbar
colormap gray
%% %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%% Export figure
%%%%%%%%%%%%%%%%%
figDir = fullfile(workDir,'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end
arrayfun(@(a) axtoolbar(a,'Visible','off'), findall(gcf,'Type','axes'));   % suppress the axes toolbar icon from the export
exportgraphics(gcf, fullfile(figDir,'populationModel.png'), 'Resolution',300)
%% %%%%%%%%%%%%%%%
