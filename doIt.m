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
rateSource = 'conventional';   % 'conventional' (literal Siler 1979 mortality + Hadwiger 1940
                                % fertility -- the standard demographic functional forms; see
                                % populationModel.md Sec. 11) | 'rough' (this project's own
                                % hand-built 3-term mortality + triangular fertility -- simpler
                                % closed form, but a nonstandard shape; Secs. 1/3) | 'empirical'
                                % (current real-world age-specific rates -- World; fertility from
                                % UN World Population Prospects, mortality from WHO Global Health
                                % Observatory, both via Our World in Data). Whichever is picked,
                                % whether fertility gets rescaled is controlled by targetR0 below,
                                % not by this choice. See README "Rate source" and
                                % populationModel.md Secs. 9/11 for details/sources/citations.

% ageMax is the SIMULATION's oldest age class (plus-group); the rough/conventional forward models
% below use fixed absolute reference ages internally, independent of it (an earlier version
% normalized them BY ageMax, so changing ageMax to reduce plus-group truncation silently changed
% the fitted rates too -- fixed). ageMax can be set arbitrarily large without changing
% mortality(age)/fertility(age) at all -- only ageMaxDisplay controls what the plots show.
ageMax        = 1000;   % y, oldest age class (plus-group: survivors accumulate here)
ageMaxDisplay = 100;     % y, x/y-axis limit for age in the plots ONLY; does not affect the simulation at all
nYears        = 500;     % y, simulation horizon
popInit       = 8e9;     % total starting population, distributed across ages per the real-world
                          % 2021-2023 world age structure below (itself NOT the model's own stable age
                          % distribution, so the model's convergence to its own shape over time is visible)

% Net reproduction rate (R0) target -- expected offspring per person over a lifetime, given the
% mortality schedule above. If non-empty, fertility is rescaled to hit this target exactly, which
% is what makes the population's LONG-RUN total steady (targetR0=1) rather than growing/shrinking
% (>1/<1) -- see populationModel.md for why. If EMPTY, fertility is left as-is and R0 is simply
% estimated/reported (in the figure title) instead of targeted -- pairs naturally with
% rateSource='empirical' or 'conventional', to see what the real data's own rates imply.
targetR0 = [];
%% =================================================
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%% Run simulation
%%%%%%%%%%%%%%%%%%%
ageVec = (0:ageMax)';
nAge   = numel(ageVec);

%%% Fit the rough and conventional models to the CURRENT empirical data. This always runs,
%%% regardless of rateSource, so switching modes is instant and every fitted parameter struct is
%%% always available to inspect or manually tweak, e.g.:
%%%     roughMortParams.steepAge = 20;   mortality = mortalityRough(ageVec, roughMortParams);
%%% If the empirical data itself changes (edit empiricalRates below), these calls automatically
%%% refit both models to it -- no separate manual re-fitting step required.
[mortalityEmp, fertilityEmp] = empiricalRates(ageVec);
survivorshipEmp = cumprod([1; exp(-mortalityEmp(1:end-1))]);   % for weighting the fertility fits' R0-matching term
R0emp = sum(survivorshipEmp .* fertilityEmp);

roughMortParams = fitMortalityRough(ageVec, mortalityEmp);
convMortParams  = fitMortalityConventional(ageVec, mortalityEmp);

survivalRough      = exp(-mortalityRough(ageVec, roughMortParams));
survivorshipRough  = cumprod([1; survivalRough(1:end-1)]);
roughFertParams    = fitFertilityRough(ageVec, fertilityEmp, survivorshipRough, R0emp);

survivalConv        = exp(-mortalityConventional(ageVec, convMortParams));
survivorshipConv    = cumprod([1; survivalConv(1:end-1)]);
convFertParams       = fitFertilityConventional(ageVec, fertilityEmp, survivorshipConv, R0emp);

%%% age-specific rates (both birth and death vary with age, not constant across the population)
switch rateSource
case 'rough'
    mortality      = mortalityRough(ageVec, roughMortParams);
    fertilityShape = fertilityRough(ageVec, roughFertParams);
case 'conventional'
    mortality      = mortalityConventional(ageVec, convMortParams);
    fertilityShape = fertilityConventional(ageVec, convFertParams);
case 'empirical'
    mortality      = mortalityEmp;
    fertilityShape = fertilityEmp;
end
survival = exp(-mortality);   % per-capita annual survival probability

% net reproduction rate (R0) this run actually has: sum over ages of (probability of surviving
% birth-to-age) x (fertility at age). If targetR0 is non-empty, fertility is rescaled so this hits
% targetR0 exactly (works for any rateSource); if targetR0 is empty, fertility is used as-is and
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
                           % a log axis (log(0) is undefined) -- this, and the xlim/ylim below, only
                           % affect the DISPLAY, not the simulation, which always uses the full
                           % mortality(age)/fertility(age) arrays (out to the real, large ageMax).
xlabel('age'); ylabel('per-capita annual rate (log scale)')
title('age-dependent rates')
legend({'death rate','fertility rate'}, 'Location','best')
grid on
xlim([0 ageMaxDisplay])
ylim([1e-6 1])   % floor the log axis explicitly: fertilityConventional (Hadwiger) is smooth and
                  % technically nonzero at every age > 1, unlike fertilityRough/empirical's hard
                  % cutoff -- near age 0 it decays to astronomically small (e.g. 1e-39) but genuine
                  % values, which would otherwise blow the auto-scaled axis out to 40 orders of
                  % magnitude and squash every rate that actually matters into an unreadable sliver

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
function mortality = mortalityRough(ageVec, p)
% mortalityRough  This project's own 3-term mortality hazard: decaying infant excess + slowly
%   rising baseline + thresholded senescence cliff (populationModel.md Eq. 1).
%   mortality = mortalityRough(ageVec) uses hardcoded default parameters, fit to the current
%   empirical mortality curve (see fitMortalityRough).
%   mortality = mortalityRough(ageVec, p) uses the given parameter struct instead (same fields as
%   the default: infantMortalityScale, infantMortalityDecay, deathRateBase, deathRateSlope,
%   ageRef, steepAge, steepMortalityScale, steepMortalityRate).
    if nargin < 2 || isempty(p)
        p.infantMortalityScale = 0.012589;  % excess hazard at age 0 (fitMortalityRough: 0.012589407)
        p.infantMortalityDecay = 1.2306;    % y, decay time constant (fitMortalityRough: 1.2306087)
        p.deathRateBase        = 0.00064457;% baseline hazard at age 0 (fitMortalityRough: 0.00064457167)
        p.deathRateSlope       = 0.00057295;% added baseline hazard at age ageRef (fitMortalityRough: 0.00057294626)
        p.ageRef               = 150;       % y, fixed reference age -- NOT ageMax
        p.steepAge             = 11.0;      % y, absolute cliff-onset age (fitMortalityRough: 11.000)
        p.steepMortalityScale  = 0.00023021;% extra-hazard scale past steepAge (fitMortalityRough: 0.00023020666)
        p.steepMortalityRate   = 0.083487;  % extra-hazard growth rate past steepAge (fitMortalityRough: 0.083486942)
    end
    mortality = p.infantMortalityScale * exp(-ageVec/p.infantMortalityDecay) + ...
                p.deathRateBase + p.deathRateSlope * (ageVec/p.ageRef).^2;
    pastSteepAge = ageVec > p.steepAge;
    mortality(pastSteepAge) = mortality(pastSteepAge) + p.steepMortalityScale * ...
        (exp(p.steepMortalityRate*(ageVec(pastSteepAge)-p.steepAge)) - 1);
end

function fertility = fertilityRough(ageVec, p)
% fertilityRough  This project's own triangular fertility pulse (populationModel.md Eq. 3).
%   fertility = fertilityRough(ageVec) uses hardcoded default parameters, fit to the current
%   empirical fertility curve (see fitFertilityRough).
%   fertility = fertilityRough(ageVec, p) uses the given parameter struct instead (fields:
%   fertileMin, fertileMax, fertilityPeakAge, fertilityPeakRate).
    if nargin < 2 || isempty(p)
        p.fertileMin        = 12.992;  % y, youngest fertile age (fitFertilityRough: 12.992296)
        p.fertileMax        = 42.652;  % y, oldest fertile age (fitFertilityRough: 42.652)
        p.fertilityPeakAge  = 27.037;  % y, age of peak fertility (fitFertilityRough: 27.037019)
        p.fertilityPeakRate = 0.071173;% per-capita rate at fertilityPeakAge (fitFertilityRough: 0.071173414)
    end
    fertileAges = ageVec >= p.fertileMin & ageVec <= p.fertileMax;
    fertility = zeros(size(ageVec));
    halfWidth = max(p.fertilityPeakAge-p.fertileMin, p.fertileMax-p.fertilityPeakAge);
    fertility(fertileAges) = p.fertilityPeakRate * max(0, 1 - abs(ageVec(fertileAges)-p.fertilityPeakAge)/halfWidth);
end

function mortality = mortalityConventional(ageVec, p)
% mortalityConventional  Literal Siler (1979) competing-hazards mortality law:
%   h(a) = a1*exp(-a2*a) + a3 + a4*exp(a5*a) -- a decaying infant/juvenile hazard, a true constant
%   background hazard, and an UNTHRESHOLDED exponentially-rising senescent hazard (unlike
%   mortalityRough's ad hoc cliff and rising-quadratic background; see populationModel.md Sec. 11).
%   mortality = mortalityConventional(ageVec) uses hardcoded default parameters, fit to the
%   current empirical mortality curve (see fitMortalityConventional).
%   mortality = mortalityConventional(ageVec, p) uses the given parameter struct instead (fields:
%   a1, a2, a3, a4, a5).
    if nargin < 2 || isempty(p)
        p.a1 = 0.012163;     % infant/juvenile hazard scale at age 0 (fitMortalityConventional: 0.012162645)
        p.a2 = 0.75018;      % infant/juvenile hazard decay rate, 1/y (fitMortalityConventional: 0.75017859)
        p.a3 = 0.00043079;   % constant background hazard (fitMortalityConventional: 0.0004307947)
        p.a4 = 9.1977e-05;   % senescent hazard scale (fitMortalityConventional: 9.1977067e-05)
        p.a5 = 0.083486;     % senescent hazard growth rate, 1/y (fitMortalityConventional: 0.083485611)
    end
    mortality = p.a1*exp(-p.a2*ageVec) + p.a3 + p.a4*exp(p.a5*ageVec);
end

function fertility = fertilityConventional(ageVec, p)
% fertilityConventional  Literal Hadwiger (1940) fertility function:
%   ASFR(a) = (a*b/c)*(c/a)^1.5 * exp(-b^2*(c/a + a/c - 2)) -- a smooth, right-skewed unimodal
%   curve (unlike fertilityRough's piecewise-linear triangle; see populationModel.md Sec. 11).
%   Singular at age 0 (division by age); ages at or below 1 are set to exactly zero rather than
%   evaluating the formula there (real fertility is negligible at those ages regardless, and the
%   formula's age->0 limit is 0 but numerically ill-conditioned right at it -- Inf*0).
%   fertility = fertilityConventional(ageVec) uses hardcoded default parameters, fit to the
%   current empirical fertility curve (see fitFertilityConventional).
%   fertility = fertilityConventional(ageVec, p) uses the given parameter struct instead (fields:
%   a, b, c).
    if nargin < 2 || isempty(p)
        p.a = 0.62758;    % overall level (fitFertilityConventional: 0.62758139)
        p.b = 2.7152;     % spread/shape (fitFertilityConventional: 2.7151723)
        p.c = 28.199;     % y, modal (peak) age (fitFertilityConventional: 28.198696)
    end
    fertility = zeros(size(ageVec));
    valid = ageVec > 1;
    x = ageVec(valid);
    fertility(valid) = (p.a*p.b/p.c) .* (p.c./x).^1.5 .* exp(-p.b^2*(p.c./x + x/p.c - 2));
end

function [mortality, fertility] = empiricalRates(ageVec)
% empiricalRates  Current real-world age-specific mortality/fertility rates for the World (see
%   populationModel.md Sec. 9 for full citations/derivation). Fertility: 2023 age-specific
%   fertility rate by 5-year age band (UN World Population Prospects, via Our World in Data).
%   Mortality: 2021 life-table probability of dying within each age band (WHO Global Health
%   Observatory, via Our World in Data) -- NOT UN WPP. Both converted to this model's per-capita
%   ANNUAL rate convention (mortality as a hazard, with survival = exp(-mortality); fertility as a
%   per-capita, not per-woman, rate). Mortality beyond the data's oldest covered age (84) is
%   extrapolated with a Gompertz (log-linear) fit to the last 5 empirical bands, out to ageVec's
%   own range. This is the fitting TARGET for fitMortalityRough/Conventional and
%   fitFertilityRough/Conventional below, and is used directly when rateSource='empirical'.
    nAge = numel(ageVec);
    fertAgeBins = [10 14; 15 19; 20 24; 25 29; 30 34; 35 39; 40 44; 45 49; 50 54];   % y, [lo hi] inclusive; zero outside
    fertAsfr    = [1.053 39.005 115.997 127.830 93.972 51.179 17.612 3.187 0.194];   % births per 1000 women per year
    femaleFrac  = 0.5;                                                              % share of each age cohort assumed female (unisex model; see populationModel.md assumption 2)

    mortAgeBins = [0 0; 1 4; 5 9; 10 14; 15 19; 20 24; 25 29; 30 34; 35 39; 40 44; 45 49; 50 54; 55 59; 60 64; 65 69; 70 74; 75 79; 80 84];
    mortQx      = [0.028159427 0.009127571 0.003440582 0.002693729 0.004317795 0.005813476 0.006770404 0.008635458 0.01193883 0.017490493 0.024098558 0.03552827 0.05312812 0.07979178 0.11539516 0.17025621 0.2427777 0.3592315];   % P(die within band | alive at band start)

    fertility = zeros(nAge,1);
    for b = 1:size(fertAgeBins,1)
        inBin = ageVec >= fertAgeBins(b,1) & ageVec <= fertAgeBins(b,2);
        fertility(inBin) = (fertAsfr(b)/1000) * femaleFrac;   % per-capita annual rate
    end

    mortality = nan(nAge,1);
    for b = 1:size(mortAgeBins,1)
        inBin = ageVec >= mortAgeBins(b,1) & ageVec <= mortAgeBins(b,2);
        bandWidth = mortAgeBins(b,2) - mortAgeBins(b,1) + 1;
        mortality(inBin) = -log(1 - mortQx(b)) / bandWidth;   % per-capita annual hazard from the band's death probability
    end
    lastCovered = mortAgeBins(end,2);
    fitAges = mean(mortAgeBins(end-4:end,:),2);
    fitHaz  = -log(1 - mortQx(end-4:end)') ./ (mortAgeBins(end-4:end,2)-mortAgeBins(end-4:end,1)+1);
    pTail = polyfit(fitAges, log(fitHaz), 1);
    beyondData = ageVec > lastCovered;
    mortality(beyondData) = exp(polyval(pTail, ageVec(beyondData)));
end

function xBest = multiStartFminsearch(obj, x0list, opts, nPolish)
% multiStartFminsearch  Runs fminsearch from each row of x0list (each chained nPolish times, since
%   Nelder-Mead can stall after a single pass), and returns the lowest-objective result. Guards
%   against Nelder-Mead's sensitivity to its starting simplex landing in a locally-good-but-
%   globally-worse optimum -- seen concretely fitting fertilityRough's triangular shape, where a
%   single fixed starting point sometimes converged with ~50% higher error (and a visibly
%   R0-mismatched result) than a different, equally-generic starting point found from the same
%   empirical data.
    bestObj = inf; xBest = x0list(1,:);
    for i = 1:size(x0list,1)
        x = x0list(i,:);
        for k = 1:nPolish
            x = fminsearch(obj, x, opts);
        end
        o = obj(x);
        if o < bestObj
            bestObj = o; xBest = x;
        end
    end
end

function p = fitMortalityRough(ageVec, mortEmp)
% fitMortalityRough  Least-squares fit of mortalityRough to a target mortality curve, in
%   log-hazard space (mortality spans several orders of magnitude with age). The objective calls
%   mortalityRough itself (not a re-derived copy of its formula), so the fit can never drift out
%   of sync with the forward model it's fitting. Uses multi-start fminsearch (Nelder-Mead; see
%   multiStartFminsearch), reparametrized through exp() to keep every rate/scale parameter
%   positive; ageRef is fixed at 150 (not fit -- see populationModel.md Sec. 1a's note on why
%   mortalityRough's reference age is a fixed constant, not tied to ageVec's own range).
%   p = fitMortalityRough(ageVec, mortEmp) returns a struct with the same fields as
%   mortalityRough's own default parameters, fit so mortalityRough(ageVec,p) best matches mortEmp.
    ageRef = 150;
    toParams = @(x) struct('infantMortalityScale',exp(x(1)), 'infantMortalityDecay',exp(x(2)), ...
        'deathRateBase',exp(x(3)), 'deathRateSlope',exp(x(4)), 'ageRef',ageRef, ...
        'steepAge',exp(x(5)), 'steepMortalityScale',exp(x(6)), 'steepMortalityRate',exp(x(7)));
    obj = @(x) sum((log(mortalityRough(ageVec,toParams(x))) - log(mortEmp)).^2);
    x0list = [log(0.02) log(2)   log(0.001) log(0.001) log(10) log(0.0005) log(0.08);   % generic starting points, not tuned to any one dataset -- steepAge varied since its threshold kink is the likeliest source of a bad local optimum
              log(0.02) log(2)   log(0.001) log(0.001) log(15) log(0.0005) log(0.08);
              log(0.02) log(2)   log(0.001) log(0.001) log(20) log(0.0005) log(0.08)];
    opts = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-12,'TolX',1e-12);
    x = multiStartFminsearch(obj, x0list, opts, 2);
    p = toParams(x);
end

function p = fitMortalityConventional(ageVec, mortEmp)
% fitMortalityConventional  Least-squares fit of the literal Siler (1979) mortality law to a
%   target mortality curve, in log-hazard space. The objective calls mortalityConventional
%   itself, so the fit can never drift out of sync with the forward model it's fitting. Uses
%   multi-start fminsearch, reparametrized through exp() to keep all 5 parameters positive
%   (Siler's own requirement).
%   p = fitMortalityConventional(ageVec, mortEmp) returns a struct with fields a1..a5 (see
%   mortalityConventional), fit so mortalityConventional(ageVec,p) best matches mortEmp.
    toParams = @(x) struct('a1',exp(x(1)), 'a2',exp(x(2)), 'a3',exp(x(3)), 'a4',exp(x(4)), 'a5',exp(x(5)));
    obj = @(x) sum((log(mortalityConventional(ageVec,toParams(x))) - log(mortEmp)).^2);
    x0list = [log(0.02) log(0.5) log(0.001) log(0.0005) log(0.08);   % generic starting points
              log(0.05) log(1.0) log(0.001) log(0.0002) log(0.10);
              log(0.01) log(0.3) log(0.002) log(0.0010) log(0.06)];
    opts = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-12,'TolX',1e-12);
    x = multiStartFminsearch(obj, x0list, opts, 2);
    p = toParams(x);
end

function p = fitFertilityRough(ageVec, fertEmp, survivorship, R0target)
% fitFertilityRough  Least-squares fit of fertilityRough to a target fertility curve, PLUS a
%   penalty on the resulting net reproduction rate R0's mismatch from R0target, weighted by the
%   given survivorship curve (normally the paired mortality model's own, since that's what R0
%   actually gets computed against once plugged into the simulation -- see populationModel.md
%   Sec. 10). A pointwise-only fit (no R0 penalty) matches shape well but tends to undershoot R0,
%   since it doesn't preserve the curve's area; the penalty corrects that. The objective calls
%   fertilityRough itself, so the fit can never drift out of sync with the forward model it's
%   fitting -- an earlier version re-derived the triangle formula inline and got its asymmetric
%   truncation subtly wrong (a symmetric half-width with no independent left-edge cutoff), which
%   silently fit the wrong shape while still looking like it converged.
%   Uses multi-start fminsearch (see multiStartFminsearch): the triangle's left edge
%   (fertileMin) turns out to be only weakly identified once the peak/right-edge/height are set
%   (the empirical data is already near-zero there), so a single fixed starting point can converge
%   to a visibly worse fit depending on exactly where it lands relative to that flat direction.
%   p = fitFertilityRough(ageVec, fertEmp, survivorship, R0target) returns a struct with fields
%   fertileMin, fertileMax, fertilityPeakAge, fertilityPeakRate (see fertilityRough).
    toParams = @(x) struct('fertileMin',x(1), 'fertilityPeakAge',x(1)+exp(x(2)), ...
        'fertileMax',x(1)+exp(x(2))+exp(x(3)), 'fertilityPeakRate',exp(x(4)));
    lambda = 1e5;
    obj = @(x) sum((fertilityRough(ageVec,toParams(x)) - fertEmp).^2) + ...
        lambda*(sum(survivorship.*fertilityRough(ageVec,toParams(x))) - R0target)^2;
    x0list = [12, log(15), log(15), log(0.06);   % generic fertile-window starting points, varying the left/right half-widths
              15, log(13), log(17), log(0.07);
              10, log(18), log(13), log(0.05);
              18, log(9),  log(15), log(0.07);
              20, log(7),  log(20), log(0.06)];
    opts = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-14,'TolX',1e-14);
    x = multiStartFminsearch(obj, x0list, opts, 3);
    p = toParams(x);
end

function p = fitFertilityConventional(ageVec, fertEmp, survivorship, R0target)
% fitFertilityConventional  Least-squares fit of the literal Hadwiger (1940) fertility function to
%   a target fertility curve, PLUS the same R0-matching penalty as fitFertilityRough (see there).
%   The objective calls fertilityConventional itself, so the fit can never drift out of sync with
%   the forward model it's fitting. Uses multi-start fminsearch (see multiStartFminsearch) for the
%   same robustness reason as fitFertilityRough.
%   p = fitFertilityConventional(ageVec, fertEmp, survivorship, R0target) returns a struct with
%   fields a, b, c (see fertilityConventional).
    toParams = @(x) struct('a',exp(x(1)), 'b',exp(x(2)), 'c',exp(x(3)));
    lambda = 1e5;
    obj = @(x) sum((fertilityConventional(ageVec,toParams(x)) - fertEmp).^2) + ...
        lambda*(sum(survivorship.*fertilityConventional(ageVec,toParams(x))) - R0target)^2;
    x0list = [log(3.4) log(2.5) log(22.2);   % Hadwiger's (1940) own commonly-cited starting point, plus generic perturbations
              log(2.0) log(2.0) log(25.0);
              log(5.0) log(3.0) log(20.0)];
    opts = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',1e-14,'TolX',1e-14);
    x = multiStartFminsearch(obj, x0list, opts, 2);
    p = toParams(x);
end
