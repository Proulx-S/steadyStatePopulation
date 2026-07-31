clear all
close all

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
ageMax     = 90;      % y, oldest age class (plus-group: survivors accumulate here)
steepAgeFrac        = 0.95;   % fraction of ageMax past which mortality accelerates sharply (senescence cliff)
steepMortalityScale = 0.05;   % extra-hazard scale added past steepAgeFrac*ageMax
steepMortalityRate  = 1;      % extra-hazard growth rate (per year) past steepAgeFrac*ageMax; higher = sharper cliff
nYears     = 1500;     % y, simulation horizon
popInit    = 100000;   % total starting population, distributed uniformly across ages
                       % (deliberately NOT the stable age distribution, so its emergence over time is visible)

% Mortality: per-capita annual death rate, AGE-DEPENDENT, rising with age (Gompertz-like).
deathRateBase  = 0.004;   % baseline hazard at age 0
deathRateSlope = 0.08;    % additional hazard at age ageMax (hazard = base + slope*(age/ageMax)^2)

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
targetR0 = 1.1;
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

%%% cohort simulation (Leslie matrix, applied step by step)
N = zeros(nAge, nYears+1);
N(:,1) = popInit/nAge;   % uniform initial age distribution

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
figure; ht = tiledlayout(2,2); ht.Padding = 'compact'; ht.TileSpacing = 'compact';
title(ht, sprintf('age-structured population model (R0 target = %.2f)', targetR0))

nexttile
plot(ageVec, mortality, 'r-', ageVec, fertility, 'g-', 'LineWidth', 1.5)
xlabel('age'); ylabel('per-capita annual rate')
title('age-dependent rates')
legend({'death rate','fertility rate'}, 'Location','best')
grid on; axis square

nexttile
plot(years, totalPop, 'w-', 'LineWidth', 1.5)
xlabel('year'); ylabel('total population')
title('total population over time')
grid on; axis square

nexttile
hold on
plot(ageVec, N(:,1)/sum(N(:,1)), 'c-', 'LineWidth', 1.5)
plot(ageVec, N(:,end)/sum(N(:,end)), 'y-', 'LineWidth', 1.5)
xlabel('age'); ylabel('fraction of population')
title('age distribution: initial vs. final')
legend({'initial (uniform)','final (stable)'}, 'Location','best')
grid on; axis square

ax = nexttile;
Nnorm = N ./ sum(N,1);   % normalize each year (column) to sum to 1 -- per-year age-distribution shape
imagesc(ageVec, years, Nnorm.'); set(gca,'YDir','normal')
xlabel('age'); ylabel('year')
title('age-stratified population (normalized per year)')
colorbar
axis square
%% %%%%%%%%%%%%%%%%%
