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
rateSource = 'empirical';   % 'empirical' (current real-world age-specific rates -- World, UN
                             % World Population Prospects via Our World in Data; NOT rescaled to
                             % targetR0, so its own R0 emerges from the real data as-is) |
                             % 'parameterized' (hand-tuned Gompertz+cliff mortality, triangular
                             % fertility, rescaled to targetR0). See README "Rate source" for details/sources.
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
switch rateSource
case 'parameterized'
    mortality = deathRateBase + deathRateSlope * (ageVec/ageMax).^2;   % per-capita annual hazard, increases with age
    steepAgeThreshold = steepAgeFrac * ageMax;
    pastSteepAge = ageVec > steepAgeThreshold;
    mortality(pastSteepAge) = mortality(pastSteepAge) + steepMortalityScale * (exp(steepMortalityRate*(ageVec(pastSteepAge)-steepAgeThreshold)) - 1);   % sharp senescence cliff past steepAgeFrac*ageMax

    fertileAges = ageVec >= fertileMin & ageVec <= fertileMax;
    fertilityShape = zeros(nAge,1);
    halfWidth = max(fertilityPeakAge-fertileMin, fertileMax-fertilityPeakAge);
    fertilityShape(fertileAges) = max(0, 1 - abs(ageVec(fertileAges)-fertilityPeakAge)/halfWidth);   % triangular, zero outside fertile window

    mortalityCapRef = deathRateBase + deathRateSlope * (ageMax/ageMax)^2;   % baseline mortality at ageMax with steepMortalityScale=0 (no cliff); caps the rates-plot axis

case 'empirical'
    % Real-world age-specific rates (World; Our World in Data, sourced from UN World Population
    % Prospects). Fertility: 2023 age-specific fertility rate by 5-year age band (births per 1000
    % women per year). Mortality: 2021 life-table probability of dying within each age band. Both
    % converted below to this model's per-capita ANNUAL rate convention (mortality as a hazard,
    % with survival = exp(-mortality); fertility as a per-capita, not per-woman, rate).
    fertAgeBins = [10 14; 15 19; 20 24; 25 29; 30 34; 35 39; 40 44; 45 49; 50 54];   % y, [lo hi] inclusive; zero outside
    fertAsfr    = [1.053 39.005 115.997 127.830 93.972 51.179 17.612 3.187 0.194];   % births per 1000 women per year
    femaleFrac  = 0.5;                                                              % share of each age cohort assumed female (unisex model; see populationModel.md assumption 2)

    mortAgeBins = [0 0; 1 4; 5 9; 10 14; 15 19; 20 24; 25 29; 30 34; 35 39; 40 44; 45 49; 50 54; 55 59; 60 64; 65 69; 70 74; 75 79; 80 84];
    mortQx      = [0.028159427 0.009127571 0.003440582 0.002693729 0.004317795 0.005813476 0.006770404 0.008635458 0.01193883 0.017490493 0.024098558 0.03552827 0.05312812 0.07979178 0.11539516 0.17025621 0.2427777 0.3592315];   % P(die within band | alive at band start)

    fertilityShape = zeros(nAge,1);
    for b = 1:size(fertAgeBins,1)
        inBin = ageVec >= fertAgeBins(b,1) & ageVec <= fertAgeBins(b,2);
        fertilityShape(inBin) = (fertAsfr(b)/1000) * femaleFrac;   % per-capita annual rate
    end

    mortality = nan(nAge,1);
    for b = 1:size(mortAgeBins,1)
        inBin = ageVec >= mortAgeBins(b,1) & ageVec <= mortAgeBins(b,2);
        bandWidth = mortAgeBins(b,2) - mortAgeBins(b,1) + 1;
        mortality(inBin) = -log(1 - mortQx(b)) / bandWidth;   % per-capita annual hazard from the band's death probability
    end
    % extrapolate past the data's oldest covered age (84) with a Gompertz (log-linear) fit to the
    % last 5 empirical bands, rather than leaving mortality undefined out to ageMax
    lastCovered = mortAgeBins(end,2);
    fitAges = mean(mortAgeBins(end-4:end,:),2);
    fitHaz  = -log(1 - mortQx(end-4:end)') ./ (mortAgeBins(end-4:end,2)-mortAgeBins(end-4:end,1)+1);
    p = polyfit(fitAges, log(fitHaz), 1);
    beyondData = ageVec > lastCovered;
    mortality(beyondData) = exp(polyval(p, ageVec(beyondData)));

    mortalityCapRef = mortality(ageVec==lastCovered);   % cap the rates-plot axis at the last directly-observed (non-extrapolated) hazard
end
survival = exp(-mortality);   % per-capita annual survival probability

% net reproduction rate (R0) this run actually has: sum over ages of (probability of surviving
% birth-to-age) x (fertility at age). In 'parameterized' mode fertility is rescaled so this hits
% targetR0 exactly; in 'empirical' mode it's left as the real data implies (see rateSource above).
survivorship = cumprod([1; survival(1:end-1)]);   % l(age) = P(alive at age | born)
if strcmp(rateSource,'parameterized')
    R0_unscaled = sum(survivorship .* fertilityShape);
    fertility   = fertilityShape * (targetR0 / R0_unscaled);   % age-dependent fertility rate, rescaled to targetR0
else
    fertility = fertilityShape;   % real-world fertility rate, used as-is (not rescaled)
end
R0 = sum(survivorship .* fertility);   % actual net reproduction rate this run ends up with

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
title(ht, sprintf('age-structured population model (%s rates, R0 = %.2f)', rateSource, R0), 'Color','w')

nexttile
plot(ageVec, mortality, 'r-', ageVec, fertility, 'g-', 'LineWidth', 1.5)
xlabel('age'); ylabel('per-capita annual rate')
title('age-dependent rates')
legend({'death rate','fertility rate'}, 'Location','best')
grid on
yl = ylim; yl(2) = max(max(fertility), mortalityCapRef); ylim(yl);

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
