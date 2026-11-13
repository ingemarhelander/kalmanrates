function [stop, reason] = evaluate(fit, alphaHat, alphas, ok, F, fj, f0, opts)
% evaluate - Stoppkriterier för linjesökningen (komponent 10): ska line search
% stanna innan alphaHat utvärderas?
%
% Anropas efter att den senaste punkten alpha_j = alphas(end) utvärderats och
% polynomet anpassats om, men innan nästa steg alphaHat utvärderas. Stoppar om
%   1. lutningen längs d har planat ut, |f'_j| <= c |f'_0|, och någon giltig
%      punkt har högre uppmätt F än startpunkten, max(F(ok)) > F(1). Vi söker
%      inte den bästa punkten, bara en bättre. Returpunkten (högst uppmätt F)
%      väljs i line_search_main, så villkoret betyder att den är bättre än
%      startpunkten,
%   2. alphaHat ligger för nära en redan utvärderad punkt (Req. 31), och
%      någon giltig punkt har högre uppmätt F än startpunkten. Utan
%      förbättring stoppar vi inte här, eftersom startpunkten då skulle
%      returneras med alpha = 0. line_search_main ersätter i stället
%      alphaHat med ett kortare steg (opts.shrinkFactor, opts.minStepRel),
%   3. den förutsagda ökningen p(alphaHat) - max_j p(alpha_j) är mindre än
%      opts.noiseFactor gånger bruset, och någon giltig punkt har högre uppmätt F än startpunkten.
%      Med få punkter kan fit.sigmaF vara modellfel snarare än brus och bli
%      stort; utan kravet på förbättring kunde sökningen då stoppa i
%      startpunkten fast f'_0 > 0.
% Kriterierna prövas i den ordningen; reason anger det första som uppfylls.
%
% fit:    anpassning till alla giltiga punkter (alphas(end) är alltid giltig här)
% alphas: steglängd för alla utvärderade punkter, även ogiltiga (NaN/Inf),
%         så att vi inte går tillbaka till en punkt som redan misslyckats.
% ok:     logisk vektor, giltiga punkter.
% F:      uppmätt F för alla utvärderade punkter, startpunkten först.
% fj, f0: uppmätt riktningsderivata g'*d i alpha_j respektive startpunkten.
% opts.c:       konstant i lutningsvillkoret. För kvadratisk f tappas högst
%               andelen c^2 av den möjliga ökningen längs d.
% opts.dMinRel: minsta avstånd relativt största steget. Relativt stegen och
%               inte alphaMax, som kan vara mycket större än stegen.
% opts.noiseFactor: faktor på fit.sigmaF i kriterium 3. Större värde ger
%               tidigare stopp.
% opts.minPoints: kriterium 1 och 3 får stoppa först när minst så många
%               punkter (startpunkten inräknad) är utvärderade.
% Standardvärden sätts i line_search_main (default_opts).
%
% Bruset i riktningsderivatan antas vara litet jämfört med c|f'_0|, så
% lutningsvillkoret har ingen brusspärr (jfr fit.sigmaG om det behövs).
% maxEvals hanteras av loopen i line_search_main, inte här.

    stop = true;
    improved = max(F(ok)) > F(1);    % returpunkten är bättre än startpunkten
    enough = numel(alphas) >= opts.minPoints;   % tillräckligt många punkter för kriterium 1 och 3
    if enough && abs(fj) <= opts.c*abs(f0) && improved
        reason = sprintf('lutningen har planat ut, |f''_j| <= %.2g |f''_0|, och F har ökat', opts.c);
        return
    end
    if improved && too_close(alphaHat, alphas, opts.dMinRel)
        reason = 'för nära en utvärderad punkt (Req. 31)';
        return
    end
    % fit.sigmaF är ~0 vid exakt anpassning. Då stoppar bara de andra kriterierna.
    if enough && improved && poly_value(fit, alphaHat) - max(poly_value(fit, alphas(ok))) < opts.noiseFactor*fit.sigmaF
        reason = sprintf('förutsagd ökning mindre än %.2g * bruset, och F har ökat', opts.noiseFactor);
        return
    end
    stop = false;
    reason = '';
end
