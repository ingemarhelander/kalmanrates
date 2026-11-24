% plots.m
% Figure 1: objective value and gradient norm after each outer iteration of
%           kalmanOptimizeMain, on a shared x-axis. Iteration 0 is the start point.
% Figure 2: market quotes and the model's one-day-ahead prediction for one day
%           at a time (slider), and histograms of the daily error (RMSE) in the
%           optimization period and the out-of-sample period.
% Run kalmanOptimizeMain first; it saves data/kalmanOptimizeResult.mat.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
burnIn = 20; % First filter days left out of the histograms (start-up transient)

R = load(dataPath('kalmanOptimizeResult.mat'));
h = R.history;
it = 0:numel(h.fvals)-1;

tl = tiledlayout(figure('Name', 'Konvergens'), 2, 1, 'TileSpacing', 'compact');
ax1 = nexttile(tl);
plot(ax1, it, h.fvals, 'o-');
ylabel(ax1, 'målfunktion f');
grid(ax1, 'on');

ax2 = nexttile(tl);
semilogy(ax2, it, h.normOfgrads, 'o-');
ylabel(ax2, '||g||');
xlabel(ax2, 'yttre iteration');
grid(ax2, 'on');
linkaxes([ax1 ax2], 'x');

%% Market vs model per day
[K, times, ricAll] = loadKalmanModel(R.modelFile);
optIdx = selectDays(times, R.periods.optimizeStart, R.periods.optimizeEnd);
oosIdx = selectDays(times, R.periods.outOfSampleStart, R.periods.outOfSampleEnd);
idx = optIdx(1):oosIdx(end);
r = filterDays(subsetKalmanModel(K, idx), R.xOpt);

d = struct(); % kalmanOptimizeMain leaves a vector d in the workspace
[d.tenor, order] = sort(cellfun(@ricTenorYears, ricAll(1:K.nu)));
d.z = 100*r.z(:, order);                                 % %
d.pred = 100*r.pred(:, order);
d.rmse = 1e4*sqrt(mean((r.z - r.pred).^2, 2, 'omitnan')); % bp
d.days = floor(times(idx));
d.isOos = ismember(idx, oosIdx)';
isOpt = ismember(idx, optIdx)';
isOpt(1:burnIn) = false;
d.ylim = [min(d.z(burnIn+1:end,:), [], 'all') max(d.z(burnIn+1:end,:), [], 'all')] + [-0.1 0.1]; % same for all days

fig = figure('Name', 'Marknad mot modell');
tl = tiledlayout(fig, 2, 1);
tl.OuterPosition = [0 0.06 1 0.94]; % room for the slider
axC = nexttile(tl);
axH = nexttile(tl); hold(axH, 'on');
edges = linspace(0, prctile(d.rmse(isOpt | d.isOos), 99), 40); % top 1 % left out
histogram(axH, d.rmse(isOpt), edges, 'Normalization', 'probability', 'DisplayName', 'optimering (in-sample)');
histogram(axH, d.rmse(d.isOos), edges, 'Normalization', 'probability', 'DisplayName', 'out-of-sample');
xlabel(axH, 'RMSE per dag (bp)'); ylabel(axH, 'andel dagar');
legend(axH); grid(axH, 'on');

n = numel(d.days);
s = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.1 0.01 0.8 0.03], ...
              'Min', 1, 'Max', n, 'Value', burnIn+1, 'SliderStep', [1 20]/(n-1));
addlistener(s, 'ContinuousValueChange', @(~, ~) drawDay(axC, d, round(s.Value))); % updates while dragging
drawDay(axC, d, burnIn+1);

%% ------------------------------ Functions -------------------------------

function drawDay(ax, d, t)
  have = ~isnan(d.z(t,:));
  plot(ax, d.tenor(have), d.z(t,have), 'o', d.tenor(have), d.pred(t,have), '-');
  if d.isOos(t), tag = 'out-of-sample'; else, tag = 'in-sample'; end
  title(ax, sprintf('%s (%s), RMSE %.1f bp', datestr(d.days(t), 'yyyy-mm-dd'), tag, d.rmse(t)));
  ylim(ax, d.ylim);
  xlabel(ax, 'löptid (år)'); ylabel(ax, 'ränta (%)');
  legend(ax, 'marknad', 'modell, prognos en dag framåt', 'Location', 'southeast');
  grid(ax, 'on');
end

function y = ricTenorYears(ric)
% 'EUREONSW=' -> 1/52, 'EUREON3M=' -> 0.25, 'EUREON2Y=' -> 2, overnight -> 1/365.
  tok = regexp(ric, 'EUREON(\d*|S)([WMY])=', 'tokens', 'once');
  if isempty(tok)
    y = 1/365;
    return;
  end
  n = str2double(tok{1});
  if isnan(n), n = 1; end
  y = n * struct('W', 7/365, 'M', 1/12, 'Y', 1).(tok{2});
end
