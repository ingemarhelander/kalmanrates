classdef testLineSearch < matlab.unittest.TestCase
% testLineSearch - Testbänk för line_search utan riktig data och mex-filer.
%
% Kör från projektroten med:  setupPaths; runtests('tests')
%
% line_search.m är i dag en prototyp med backtracking och ska ersättas av
% polynommetoden i designspecifikationens avsnitt 5.2.4. Testerna kontrollerar
% därför vad metoden ska uppfylla, inte hur prototypen råkar arbeta. Antal
% halveringar och villkoret fNew > f testas alltså inte.
%
% Testerna är uppdelade i två grupper:
%   Grupp A  Gäller både prototypen och den färdiga metoden. Ska alltid gå igenom.
%   Grupp B  Krav som prototypen inte uppfyller än. De hoppas över och visas som
%            Incomplete tills PolynomialLineSearchReady sätts till true.
%
% Bruset gör att värdena från stubben inte kan jämföras direkt. När testerna
% avgör om en punkt verkligen är bättre används därför trueValue, som räknar
% ut det brusfria värdet.

    properties (Constant)
        % Sätt till true när polynommetoden är på plats, så körs grupp B.
        PolynomialLineSearchReady = false

        % Övre gräns för antalet kolumner i evals (startpunkten inräknad).
        % Maxantalet utvärderingar är ett öppet beslut i designspecifikationen,
        % så gränsen är generös. Justera den när beslutet är taget.
        MaxEvals = 60

        % Frön för brustesterna. Flera frön gör att ett test inte kan gå
        % igenom bara för att ett enskilt frö råkade ge gynnsamt brus.
        Seeds = 1:20

        % Brusnivå i brustesterna, liten jämfört med vinsten i ett bra steg.
        NoiseStd = 0.05

        % Lägre brusnivå för de icke-konkava funktionerna, vars värden längs
        % linjen är av storleksordningen 0.1 till 1.
        NonConcaveNoiseStd = 0.005
    end

    properties (TestParameter)
        % Testerna med denna parameter körs en gång per målfunktion.
        quarticWeight = struct('kvadratisk', 0, 'kvartisk', 0.1)

        % Icke-konkava former längs sökriktningen, se nonConcaveProblem.
        lineShape = {'monoton', 'tvaToppar', 'inflexion'}
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            % Sökvägarna gäller bara under testkörningen och återställs sedan.
            % setupPaths lägger inte till scripts/functions, och stubben måste
            % ligga före model/ så att den används i stället för den riktiga
            % kalmanObj. PathFixture lägger mappar först på sökvägen.
            testDir = fileparts(mfilename('fullpath'));
            projectRoot = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'scripts', 'functions')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testDir, 'stubs')));
        end
    end

    methods (TestMethodSetup)
        function restoreRandomState(testCase)
            % Testerna sätter egna frön. Återställ slumptalsgeneratorn efteråt
            % så att de inte påverkar varandra eller användarens session.
            state = rng;
            testCase.addTeardown(@() rng(state));
        end
    end

    % ------------------------------ Grupp A ---------------------------------
    methods (Test)
        function stubIsUsedInsteadOfRealObjective(testCase)
            % Om detta fallerar anropar övriga tester den riktiga kalmanObj.
            expected = fullfile('tests', 'stubs', 'kalmanObj.m');
            testCase.assertTrue(endsWith(which('kalmanObj'), expected), ...
                'Stubben ligger inte först på sökvägen.');
        end

        function outputsFollowInterface(testCase)
            [K, x] = defaultProblem();
            result = runLineSearch(K, x);

            % Den nya punkten ska ligga på linjen genom x i riktningen d.
            testCase.verifyEqual(result.xNew, x + result.alpha*result.d, 'AbsTol', 1e-12);

            % Startpunkten ska vara första kolumnen i evals.
            testCase.verifyEqual(result.evals.X(:, 1), x);
            testCase.verifyEqual(result.evals.F(1), result.f);
            testCase.verifyEqual(result.evals.G(:, 1), result.g);

            % X, F och G ska ha en kolumn per utvärderad punkt.
            nEvals = size(result.evals.X, 2);
            testCase.verifySize(result.evals.F, [1, nEvals]);
            testCase.verifySize(result.evals.G, [numel(x), nEvals]);
        end

        function returnedPointIsAnEvaluatedPoint(testCase)
            % Med brus skulle en ny utvärdering ge andra värden. Därför visar
            % testet att fNew och gNew verkligen kommer från evals.
            [K, x] = defaultProblem();
            K.noiseStdF = testCase.NoiseStd;
            K.noiseStdG = testCase.NoiseStd;
            rng(1);
            result = runLineSearch(K, x);

            distances = vecnorm(result.evals.X - result.xNew);
            [smallestDistance, column] = min(distances);
            testCase.verifyLessThan(smallestDistance, 1e-12);
            testCase.verifyEqual(result.fNew, result.evals.F(column));
            testCase.verifyEqual(result.gNew, result.evals.G(:, column));
        end

        function firstStepIsLimitedByPredictedIncrease(testCase)
            % g'*d > 100, så steget kortas till predikterad ökning 100.
            [K, x] = defaultProblem();
            result = runLineSearch(K, x);
            steps = stepLengths(result);
            testCase.assertGreaterThan(result.g'*result.d, 100);
            testCase.verifyEqual(steps(2), 100/(result.g'*result.d), 'RelTol', 1e-10);
        end

        function firstStepIsAtMostOne(testCase)
            % Liten lutning ger 100/(g'*d) > 1. Då används hela steget 1.
            [K, x] = defaultProblem();
            K.xStar = 0.1*K.xStar;
            result = runLineSearch(K, x);
            steps = stepLengths(result);
            testCase.assertLessThan(result.g'*result.d, 100);
            testCase.verifyEqual(steps(2), 1, 'RelTol', 1e-10);
        end

        function firstStepIsAtLeastLowerBound(testCase)
            % Mycket stor lutning ger 100/(g'*d) < 1e-5. Då används undre gränsen.
            [K, x] = defaultProblem();
            K.xStar = 1e4*K.xStar;
            result = runLineSearch(K, x);
            steps = stepLengths(result);
            testCase.assertGreaterThan(result.g'*result.d, 1e7);
            testCase.verifyEqual(steps(2), 1e-5, 'RelTol', 1e-10);
        end

        function allStepsStayWithinUnitInterval(testCase)
            % find_search_direction skalar d så att alla steg i [0, 1] håller
            % sig inom de undre gränserna. Längre steg är därför inte tillåtna.
            [K, x] = defaultProblem();
            result = runLineSearch(K, x);
            steps = stepLengths(result);
            testCase.verifyGreaterThanOrEqual(min(steps), -1e-12);
            testCase.verifyLessThanOrEqual(max(steps), 1 + 1e-12);
            testCase.verifyLessThanOrEqual(result.alpha, 1);
        end

        function objectiveIncreasesWithoutNoise(testCase, quarticWeight)
            [K, x] = defaultProblem();
            K.quarticWeight = quarticWeight;
            result = runLineSearch(K, x);
            testCase.verifyGreaterThan(trueValue(K, result.xNew), trueValue(K, x));
        end

        function objectiveIncreasesWithNoise(testCase, quarticWeight)
            % Bruset är litet jämfört med avståndet till maximum, så det
            % brusfria värdet ska öka för varje frö.
            [K, x] = defaultProblem();
            K.quarticWeight = quarticWeight;
            K.noiseStdF = testCase.NoiseStd;
            K.noiseStdG = testCase.NoiseStd;
            for seed = testCase.Seeds
                rng(seed);
                result = runLineSearch(K, x);
                testCase.verifyGreaterThan(trueValue(K, result.xNew), trueValue(K, x), ...
                    sprintf('Det brusfria värdet ökade inte med frö %d.', seed));
            end
        end

        function noiseNearMaximumDoesNotRuinSolution(testCase)
            % Nära maximum är den möjliga vinsten mindre än bruset. Då kan det
            % brusfria värdet sjunka något, men bara i storleksordningen brus.
            K = defaultProblem();
            K.noiseStdF = testCase.NoiseStd;
            K.noiseStdG = testCase.NoiseStd;
            x = K.xStar + 1e-3;
            tolerance = 5*testCase.NoiseStd;
            for seed = testCase.Seeds
                rng(seed);
                result = runLineSearch(K, x);
                testCase.verifyGreaterThan(trueValue(K, result.xNew), ...
                    trueValue(K, x) - tolerance, ...
                    sprintf('Lösningen försämrades för mycket med frö %d.', seed));
            end
        end

        function descentDirectionDoesNotWorsenObjective(testCase)
            % g'*d < 0: alla steg framåt är sämre, så metoden ska avslutas utan
            % att i praktiken lämna startpunkten.
            [K, x] = defaultProblem();
            [~, g] = kalmanObj(x, K, 0, 1:numel(x));
            result = runLineSearch(K, x, d = -g);
            testCase.assertLessThan(result.g'*result.d, 0);
            verifyNotWorseAndFinished(testCase, K, x, result);
        end

        function orthogonalDirectionDoesNotWorsenObjective(testCase)
            % g'*d = 0: 100/(g'*d) blir Inf, vilket inte får ge problem. På en
            % konkav funktion ger d vinkelrät mot g ingen förbättring.
            [K, x] = defaultProblem();
            [~, g] = kalmanObj(x, K, 0, 1:numel(x));
            d = zeros(size(g));
            d(1:2) = [g(2); -g(1)];
            result = runLineSearch(K, x, d = d);
            testCase.assertEqual(result.g'*result.d, 0);
            verifyNotWorseAndFinished(testCase, K, x, result);
        end

        function stopsWhenNoPointImproves(testCase)
            % Ett startvärde över funktionens maximum kan aldrig förbättras.
            % Metoden måste då avbrytas av maxgränsen.
            [K, x] = defaultProblem();
            result = runLineSearch(K, x, f = 1e6);
            testCase.verifyLessThanOrEqual(size(result.evals.X, 2), testCase.MaxEvals);
        end

        function nanOnLongStepsIsNotAccepted(testCase)
            % Startpunkten ligger innanför radien men första steget utanför.
            % Metoden ska då hitta en giltig och bättre punkt på ett kortare steg.
            [K, x] = defaultProblem();
            K.nanRadius = 4;
            result = runLineSearch(K, x);

            % Säkerställer att testet faktiskt träffade NaN-området.
            testCase.assertTrue(any(isnan(result.evals.F)), ...
                'Ingen utvärdering gav NaN, så testet kontrollerar ingenting.');
            testCase.verifyTrue(isfinite(result.fNew));
            testCase.verifyTrue(all(isfinite(result.gNew)));
            testCase.verifyGreaterThan(trueValue(K, result.xNew), trueValue(K, x));
        end

        function nonConcaveObjectiveIncreasesWithoutNoise(testCase, lineShape)
            % Utan konkavitet kan vi inte kräva maximum, men steget ska ändå
            % ge en förbättring och hålla sig inom [0, 1].
            [K, x, d] = nonConcaveProblem(lineShape);
            result = runLineSearch(K, x, d = d);
            steps = stepLengths(result);
            testCase.verifyGreaterThan(trueValue(K, result.xNew), trueValue(K, x));
            testCase.verifyGreaterThanOrEqual(min(steps), -1e-12);
            testCase.verifyLessThanOrEqual(max(steps), 1 + 1e-12);
            testCase.verifyLessThanOrEqual(size(result.evals.X, 2), testCase.MaxEvals);
        end

        function nonConcaveObjectiveWithNoiseIsNotWorsened(testCase, lineShape)
            % Med brus kan en punkt se bättre ut än den är. Det brusfria värdet
            % får därför bara sjunka i storleksordningen brus.
            [K, x, d] = nonConcaveProblem(lineShape);
            K.noiseStdF = testCase.NonConcaveNoiseStd;
            K.noiseStdG = testCase.NonConcaveNoiseStd;
            tolerance = 5*testCase.NonConcaveNoiseStd;
            for seed = testCase.Seeds
                rng(seed);
                result = runLineSearch(K, x, d = d);
                testCase.verifyGreaterThan(trueValue(K, result.xNew), ...
                    trueValue(K, x) - tolerance, ...
                    sprintf('Det brusfria värdet sjönk för mycket med frö %d.', seed));
            end
        end

        function sameSeedGivesSameResult(testCase)
            [K, x] = defaultProblem();
            K.noiseStdF = testCase.NoiseStd;
            K.noiseStdG = testCase.NoiseStd;
            rng(7);
            first = runLineSearch(K, x);
            rng(7);
            second = runLineSearch(K, x);
            testCase.verifyEqual(second.alpha, first.alpha);
            testCase.verifyEqual(second.xNew, first.xNew);
            testCase.verifyEqual(second.evals, first.evals);
        end
    end

    % ------------------------------ Grupp B ---------------------------------
    methods (Test)
        function plusInfinityIsNotAccepted(testCase)
            % Icke-ändliga värden får aldrig räknas som en förbättring.
            % Prototypen accepterar +Inf eftersom Inf > f är sant.
            assumePolynomialLineSearch(testCase);
            [K, x] = defaultProblem();
            K.infRadius = 4;
            result = runLineSearch(K, x);
            testCase.assertTrue(any(isinf(result.evals.F)), ...
                'Ingen utvärdering gav Inf, så testet kontrollerar ingenting.');
            testCase.verifyTrue(isfinite(result.fNew));
        end

        function nanEverywhereKeepsStartPoint(testCase)
            % Startpunkten är då den enda giltiga punkten och ska returneras.
            % Prototypen returnerar i stället den sista punkten, med f = NaN.
            assumePolynomialLineSearch(testCase);
            [K, x] = defaultProblem();
            [f, g] = kalmanObj(x, K, 0, 1:numel(x));
            K.nanRadius = -1;  % norm(r) > -1 gäller alltid, alltså NaN överallt
            result = runLineSearch(K, x, f = f, g = g);
            testCase.verifyEqual(result.xNew, x);
            testCase.verifyEqual(result.fNew, f);
            testCase.verifyEqual(result.gNew, g);
        end

        function alphaIsCloseToMaximumOnLine(testCase, quarticWeight)
            % Polynommetoden ska hitta maximum längs linjen. Backtracking
            % stannar vid första förbättringen och uppfyller inte detta.
            assumePolynomialLineSearch(testCase);
            [K, x] = defaultProblem();
            K.quarticWeight = quarticWeight;
            result = runLineSearch(K, x);
            alphaStar = maximumOnLine(K, x, result.d);
            testCase.verifyEqual(result.alpha, alphaStar, 'RelTol', 0.1);
        end

        function alphaIsCloseToMaximumOnLineWithNoise(testCase, quarticWeight)
            % Enskilda körningar kan avvika, så testet använder medianen.
            assumePolynomialLineSearch(testCase);
            [K, x] = defaultProblem();
            K.quarticWeight = quarticWeight;
            [~, d] = kalmanObj(x, K, 0, 1:numel(x));  % brusfri riktning, samma för alla frön
            alphaStar = maximumOnLine(K, x, d);
            K.noiseStdF = testCase.NoiseStd;
            K.noiseStdG = testCase.NoiseStd;
            alphas = zeros(size(testCase.Seeds));
            for i = 1:numel(testCase.Seeds)
                rng(testCase.Seeds(i));
                result = runLineSearch(K, x, d = d);
                alphas(i) = result.alpha;
            end
            testCase.verifyEqual(median(alphas), alphaStar, 'RelTol', 0.2);
        end

        function nonConcaveStepEndsAtLocalMaximumOrEndpoint(testCase, lineShape)
            % Det globala maximumet krävs inte. Valt steg ska vara en förbättring
            % och antingen ha riktningsderivata nära noll, eller ligga i
            % ändpunkten alpha = 1 med en lutning som fortfarande pekar uppåt.
            assumePolynomialLineSearch(testCase);
            [K, x, d] = nonConcaveProblem(lineShape);
            result = runLineSearch(K, x, d = d);
            slope = trueSlope(K, result.xNew, d);
            tolerance = 0.2;
            testCase.verifyGreaterThan(trueValue(K, result.xNew), trueValue(K, x));
            if result.alpha >= 1 - 1e-12
                testCase.verifyGreaterThan(slope, -tolerance);
            else
                testCase.verifyLessThan(abs(slope), tolerance);
            end
        end
    end
end

% ---------------------------- Hjälpfunktioner -------------------------------

function [K, x] = defaultProblem()
% Standardproblem: fyra dimensioner, utan brus och utan ogiltiga områden.
% Med d = g och startpunkten x = 0 är g'*d cirka 163. Första steget kortas
% därmed till cirka 0.6, vilket är längre än maximum längs linjen (cirka 0.27).
    K.xStar = [1; -2; 0.5; 3];
    K.H = diag([1, 2, 3, 4]);
    K.quarticWeight = 0;
    K.noiseStdF = 0;
    K.noiseStdG = 0;
    K.nanRadius = Inf;
    K.infRadius = Inf;
    K.lineProfile = [];
    K.lineDirection = [];
    x = zeros(4, 1);
end

function [K, x, d] = nonConcaveProblem(shape)
% Problem där funktionen längs sökriktningen är icke-konkav. Med x = 0 och
% enhetsvektorn d gäller f(x + alpha*d) = phi(alpha), så phi är formen på linjen.
%   monoton    ökar hela vägen men växlar mellan konvex och konkav, bäst i alpha = 1
%   tvaToppar  två lokala maxima nära alpha = 1/6 och 5/6, med en dal mellan dem
%   inflexion  konvex början och sedan konkav, med maximum nära alpha = 0.69
    [K, x] = defaultProblem();
    d = ones(4, 1)/2;
    K.lineDirection = d;
    switch shape
        case 'monoton'
            K.lineProfile = @(s) deal(s + sin(4*pi*s)/(8*pi), 1 + 0.5*cos(4*pi*s));
        case 'tvaToppar'
            K.lineProfile = @(s) deal(sin(3*pi*s)/(3*pi) + 0.1*s, cos(3*pi*s) + 0.1);
        case 'inflexion'
            K.lineProfile = @(s) deal(0.1*s + 2*s^2 - 2*s^3, 0.1 + 4*s - 6*s^2);
    end
end

function result = runLineSearch(K, x, options)
% Enda stället där testerna anropar line_search. Om gränssnittet ändras, till
% exempel om K ersätts av ett function handle, räcker det att ändra här.
%
% Som standard utvärderas f och g i x, och riktningen blir d = g. Testerna kan
% ersätta dessa med namngivna argument, till exempel runLineSearch(K, x, d = -g).
    arguments
        K
        x
        options.d = []
        options.f = []
        options.g = []
    end

    mu = 0;
    varInd = 1:numel(x);
    [f, g] = kalmanObj(x, K, mu, varInd);
    if ~isempty(options.f), f = options.f; end
    if ~isempty(options.g), g = options.g; end
    d = g;
    if ~isempty(options.d), d = options.d; end

    % evalc fångar utskrifterna från line_search så att testloggen blir läsbar.
    result.printed = evalc(['[result.alpha, result.xNew, result.fNew, ' ...
        'result.gNew, result.evals] = line_search(x, f, g, d, K, mu, varInd);']);

    % Indata sparas för att testerna ska kunna jämföra med dem.
    result.x = x;
    result.f = f;
    result.g = g;
    result.d = d;
end

function value = trueValue(K, x)
% Målfunktionens brusfria värde, som testerna använder för att avgöra om en
% punkt verkligen är bättre. Ogiltiga områden stängs av, så att värdet alltid
% finns. Slumptalsströmmen sparas så att anropet inte påverkar bruset i testet.
    K.noiseStdF = 0;
    K.noiseStdG = 0;
    K.nanRadius = Inf;
    K.infRadius = Inf;
    state = rng;
    value = kalmanObj(x, K, 0, 1:numel(x));
    rng(state);
end

function slope = trueSlope(K, x, d)
% Brusfri riktningsderivata g'*d i punkten x.
    K.noiseStdF = 0;
    K.noiseStdG = 0;
    state = rng;
    [~, g] = kalmanObj(x, K, 0, 1:numel(x));
    rng(state);
    slope = g'*d;
end

function steps = stepLengths(result)
% Räknar ut steglängden för varje kolumn i evals.X enligt designspecifikationen:
% alpha_j = d'*(x_j - x)/(d'*d). evals lagrar inte steglängderna själv.
    d = result.d;
    steps = (d'*(result.evals.X - result.x))/(d'*d);
end

function alphaStar = maximumOnLine(K, x, d)
% Hittar det brusfria maximumet längs linjen x + alpha*d för alpha i [0, 1].
% fminbnd minimerar, så funktionen byter tecken.
    alphaStar = fminbnd(@(alpha) -trueValue(K, x + alpha*d), 0, 1, ...
        optimset('TolX', 1e-8));
end

function verifyNotWorseAndFinished(testCase, K, x, result)
% Gemensam kontroll för riktningar som inte går uppåt.
    testCase.verifyLessThanOrEqual(size(result.evals.X, 2), testCase.MaxEvals);
    testCase.verifyGreaterThanOrEqual(trueValue(K, result.xNew), trueValue(K, x) - 1e-6);
end

function assumePolynomialLineSearch(testCase)
% Hoppar över testet tills polynommetoden finns. Testet visas då som Incomplete
% och räknas inte som fallerat.
    testCase.assumeTrue(testCase.PolynomialLineSearchReady, ...
        'Väntar på polynommetoden. Sätt PolynomialLineSearchReady till true.');
end
