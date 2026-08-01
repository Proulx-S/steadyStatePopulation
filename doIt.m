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
rateSource = 'parameterized';   % 'empirical' (current real-world age-specific rates -- World;
                             % fertility from UN World Population Prospects, mortality from WHO
                             % Global Health Observatory, both via Our World in Data) |
                             % 'parameterized' (hand-tuned Gompertz+cliff mortality, triangular
                             % fertility). Either way, whether fertility gets rescaled is controlled
                             % by targetR0 below, not by this choice. See README "Rate source" and
                             % populationModel.md Sec. 9 for details/sources/citations.
% ageMax is the SIMULATION's oldest age class (plus-group); rateFit's ageRef/steepAge below are
% fixed absolute ages, independent of it (an earlier version normalized them BY ageMax, so changing
% ageMax to reduce plus-group truncation silently changed the fitted rates too -- fixed). ageMax can
% now be set arbitrarily large without changing mortality(age)/fertility(age) at all -- only
% ageMaxDisplay controls what the plots show.
ageMax        = 1000;   % y, oldest age class (plus-group: survivors accumulate here)
ageMaxDisplay = 100;      % y, x/y-axis limit for age in the plots ONLY; does not affect the simulation at all
nYears        = 500;     % y, simulation horizon
popInit       = 8e9;     % total starting population, distributed across ages per the real-world
                          % 2021-2023 world age structure below (itself NOT the model's own stable age
                          % distribution, so the model's convergence to its own shape over time is visible)

% 'parameterized' mode's mortality/fertility rate constants, bundled into one struct (self-populating
% defaults, like a no-arg opts function) rather than a dozen loose CONTROL PANEL variables: every
% field here was fit jointly to the empirical curves as a cohesive set, not independently tuned, so
% keeping them together avoids silently mixing values from different fits. See getFittedRateParams
% below (or populationModel.md Sec. 10) for what each field means, its exact optimizer value, and how
% it was derived. To explore off the fitted defaults, override individual fields after the call, e.g.
% rateFit = getFittedRateParams(); rateFit.steepAge = 20;
rateFit = getFittedRateParams();

% Net reproduction rate (R0) target -- expected offspring per person over a lifetime, given the
% mortality schedule above. If non-empty, fertility is rescaled to hit this target exactly, which
% is what makes the population's LONG-RUN total steady (targetR0=1) rather than growing/shrinking
% (>1/<1) -- see populationModel.md for why. If EMPTY, fertility is left as fertilityShape gives it
% (unscaled) and R0 is simply estimated/reported (in the figure title) instead of targeted --
% pairs naturally with rateSource='empirical', to see what the real data's own rates imply.
targetR0 = [];
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
    mortality = rateFit.infantMortalityScale * exp(-ageVec/rateFit.infantMortalityDecay) + ...   % elevated hazard at birth, decaying over early childhood
                rateFit.deathRateBase + rateFit.deathRateSlope * (ageVec/rateFit.ageRef).^2;      % per-capita annual hazard, increases with age (ageRef fixed, NOT ageMax -- see CONTROL PANEL note)
    pastSteepAge = ageVec > rateFit.steepAge;
    mortality(pastSteepAge) = mortality(pastSteepAge) + rateFit.steepMortalityScale * (exp(rateFit.steepMortalityRate*(ageVec(pastSteepAge)-rateFit.steepAge)) - 1);   % sharp senescence cliff past steepAge

    fertileAges = ageVec >= rateFit.fertileMin & ageVec <= rateFit.fertileMax;
    fertilityShape = zeros(nAge,1);
    halfWidth = max(rateFit.fertilityPeakAge-rateFit.fertileMin, rateFit.fertileMax-rateFit.fertilityPeakAge);
    fertilityShape(fertileAges) = rateFit.fertilityPeakRate * max(0, 1 - abs(ageVec(fertileAges)-rateFit.fertilityPeakAge)/halfWidth);   % triangular, zero outside fertile window

case 'empirical'
    % Real-world age-specific rates (World; via Our World in Data). Fertility: 2023 age-specific
    % fertility rate by 5-year age band (births per 1000 women per year), source: UN World
    % Population Prospects (2024). Mortality: 2021 life-table probability of dying within each age
    % band, source: WHO Global Health Observatory -- NOT UN WPP (corrected; see populationModel.md
    % Sec. 9 for full citations). Both converted below to this model's per-capita ANNUAL rate
    % convention (mortality as a hazard, with survival = exp(-mortality); fertility as a per-capita,
    % not per-woman, rate).
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

end
survival = exp(-mortality);   % per-capita annual survival probability

% net reproduction rate (R0) this run actually has: sum over ages of (probability of surviving
% birth-to-age) x (fertility at age). If targetR0 is non-empty, fertility is rescaled so this hits
% targetR0 exactly (works for either rateSource); if targetR0 is empty, fertility is used as-is and
% R0 is simply estimated/reported, not targeted.
survivorship = cumprod([1; survival(1:end-1)]);   % l(age) = P(alive at age | born)
if isempty(targetR0)
    fertility = fertilityShape;   % don't rescale -- estimate/report R0 as the shape (or real data) implies
else
    R0_unscaled = sum(survivorship .* fertilityShape);
    fertility   = fertilityShape * (targetR0 / R0_unscaled);   % age-dependent fertility rate, rescaled to targetR0
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
fig = figure('Color','k'); ht = tiledlayout(2,2); ht.Padding = 'compact'; ht.TileSpacing = 'compact';
title(ht, sprintf('age-structured population model (%s rates, R0 = %.2f)', rateSource, R0), 'Color','w')

nexttile
plot(ageVec, mortality, 'r-', ageVec, fertility, 'g-', 'LineWidth', 1.5)
set(gca,'YScale','log')   % log axis: mortality spans several orders of magnitude even just up
                           % to ageMaxDisplay; a linear axis can only show a small slice of that range
                           % at once. Ages/values where fertility is exactly zero simply don't plot on
                           % a log axis (log(0) is undefined) -- this, and the xlim below, only affect
                           % the DISPLAY, not the simulation, which always uses the full
                           % mortality(age)/fertility(age) arrays (out to the real, large ageMax).
xlabel('age'); ylabel('per-capita annual rate (log scale)')
title('age-dependent rates')
legend({'death rate','fertility rate'}, 'Location','best')
grid on
xlim([0 ageMaxDisplay])

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
xlim([0 ageMaxDisplay])

ax = nexttile;
Nnorm = N ./ sum(N,1);   % normalize each year (column) to sum to 1 -- per-year age-distribution shape
imagesc(years, ageVec, Nnorm); set(gca,'YDir','normal')
xlabel('year'); ylabel('age')
title('age-stratified population (normalized per year)')
colorbar
colormap gray
ylim([0 ageMaxDisplay])
%% %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%% Export figure
%%%%%%%%%%%%%%%%%
figDir = fullfile(workDir,'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end
arrayfun(@(a) axtoolbar(a,'Visible','off'), findall(gcf,'Type','axes'));   % suppress the axes toolbar icon from the export
exportgraphics(gcf, fullfile(figDir,'populationModel.png'), 'Resolution',300)
%% %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%% Function definitions
%%%%%%%%%%%%%%%%%%%%%%%
function opts = getFittedRateParams()
% getFittedRateParams  Default rateSource='parameterized' mortality/fertility rate constants.
%   opts = getFittedRateParams() returns every 'parameterized'-mode mortality/fertility parameter,
%   fully populated with its least-squares fit to the rateSource='empirical' curves (see
%   populationModel.md Sec. 10 for the full derivation). Override individual fields on the
%   returned struct to explore off the fitted defaults, e.g.
%       rateFit = getFittedRateParams(); rateFit.steepAge = 20;
%
%   Mortality (Eq. 1): a 3-term hazard -- decaying infant-mortality excess, a slowly-rising
%   baseline, and a senescence cliff past an absolute age -- fit jointly in log-hazard space via
%   fminsearch. ageRef/steepAge are fixed absolute ages (NOT fractions of the simulation's own
%   ageMax), so this fitted curve doesn't change if ageMax does.
%   Fertility (Eq. 3): a triangular pulse, fit via pointwise squared error PLUS a penalty on the
%   resulting net reproduction rate R0's mismatch from the empirical mode's own R0 -- a
%   pointwise-only fit matches shape well but still undershoots R0 by ~5%.
    opts.infantMortalityScale = 0.0129;   % excess hazard at age 0 from the infant-mortality term (optimizer: 0.0128711)
    opts.infantMortalityDecay = 1.21;     % y, its decay time constant (hazard ~ exp(-age/this), so ~gone within a few years) (optimizer: 1.20517)
    opts.deathRateBase        = 0.00066;  % baseline hazard at age 0, excluding the infant term (optimizer: 0.000659275)
    opts.deathRateSlope       = 0.00057;  % additional hazard at age ageRef (hazard = base + slope*(age/ageRef)^2) (optimizer: 0.000572547)
    opts.ageRef               = 150;      % y, fixed reference age deathRateSlope's term was fit against -- NOT ageMax
    opts.steepAge             = 11.01;    % y, absolute age past which mortality accelerates sharply (senescence cliff) -- NOT a fraction of ageMax (optimizer: 0.0733907*150=11.0086)
    opts.steepMortalityScale  = 0.000227; % extra-hazard scale added past steepAge (optimizer: 0.000226789)
    opts.steepMortalityRate   = 0.0836;   % extra-hazard growth rate (per year) past steepAge; higher = sharper cliff (optimizer: 0.083606)
    opts.fertileMin           = 12.1;     % y, youngest fertile age (optimizer: 12.0851)
    opts.fertileMax           = 42.6;     % y, oldest fertile age (optimizer: 42.6161)
    opts.fertilityPeakAge     = 27.0;     % y, age of peak fertility (optimizer: 27.0495)
    opts.fertilityPeakRate    = 0.0713;   % per-capita annual rate at fertilityPeakAge, i.e. the unscaled shape's own peak height (optimizer: 0.07127)
end
