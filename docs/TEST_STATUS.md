# TATA85 – teststatus och verifieringsplan

**Testansvarig:** Isak Ekbom  
**Senast inventerat:** 2026-11-16  
**Källor:** `main`, GitHub Issues, `docs/utveckling.md`, krav- och designspecifikationer.

> **CI-status (2026-11-16):** Den utökade femfilerssviten har
> passerat GitHub Actions på commit `107368a6` och `79214cf9`.
> Detta verifierar att de automatiserade testerna passerar,
> men inte att hela Kalmanmodellen är verifierad.

## Statusdefinitioner

| Status | Betydelse |
|---|---|
| **Approved** | Funktionaliteten har verifierats mot definierat beteende med relevant passerat test, dokumenterad miljö och commit. |
| **In testing** | Implementerad men ännu inte fullständigt verifierad, eller verifiering pågår. |
| **Fix required** | Ett reproducerbart problem eller misslyckat relevant test är dokumenterat. |
| **Not implemented** | Den avsedda funktionaliteten saknas fortfarande. |

`Design pending` kan läggas som extra notering när förväntat beteende ännu inte beslutats; det ersätter inte teststatus.

## Funktionsöversikt (preliminär)

| Funktion | Status | Källor / tester | Nästa steg |
|---|---|---|---|
| Dataurval och PCA vid kort historik | **In testing** | [#1](https://github.com/ingemarhelander/kalmanrates/issues/1), `testDataLimit.m` | CI godkänd. Komplettera med verkliga data. |
| Linjesökning – syntetisk integration | **In testing** | [#10](https://github.com/ingemarhelander/kalmanrates/issues/10), `testLineSearchMain.m` | Tidigare CI godkänd; verifiera på verklig Kalman-målfunktion. |
| ML-baserat polynomgradval | **In testing** | [#17](https://github.com/ingemarhelander/kalmanrates/issues/17), `testMlPolyMain.m` | CI godkänd; verifiera ytterligare gränsfall. |
| Likelihood ratio-baserat gradval | **In testing** | [#14](https://github.com/ingemarhelander/kalmanrates/issues/14), `testSelectPolyDegreeLr.m` | CI godkänd. Kontrollera statistiska kriterier. |
| Brusreducering/kvadratisk approximation | **In testing** | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2), `test_denoise_quadratic_approximation.m` | Tidigare CI godkänd; komplettera med fler brusnivåer och integration. |
| Stoppkriterier i linjesökningen | **In testing** / Design pending | [#36](https://github.com/ingemarhelander/kalmanrates/issues/36), [#38](https://github.com/ingemarhelander/kalmanrates/issues/38) | Fastställ vilket beteende som ska gälla, kontrollera aktuella testfall. |
| Modifierad CG/sökriktning | **In testing** | [#5](https://github.com/ingemarhelander/kalmanrates/issues/5), `find_search_direction.m` | Testa syntetisk kvadratisk funktion, ortogonalitet och gränser. |
| Föreslagen Gram–Schmidt-ändring | **Not implemented** (ändringen ej verifierad) | [#35](https://github.com/ingemarhelander/kalmanrates/issues/35) | Implementera och skapa tester för (nära) linjärt beroende. |
| Kalmanfilter och RTS-smoothing | **In testing** | `clsKalmanIr.m`, `kalmanSub.m` | Oberoende facit för liten modell; kontrollera kovarianser. |
| AD/MEX-gradient | **In testing** | `testAdCounts.m`, `checkMexRegression.m`, `docs/utveckling.md` | Kör på rätt Windows/MATLAB/MEX och kontrollera finit differens. |
| Full optimeringsloop och verkliga EUR-data | **In testing** | [#18](https://github.com/ingemarhelander/kalmanrates/issues/18), `kalmanOptimizeMain.m` | Kör fler iterationer; jämför backtracking och polynomlinjesökning. |

**Observera:** `Fix required` tilldelas först när ett faktiskt problem bekräftats. Omfattande krav-för-krav-kartläggning är nästa inventeringssteg.

## Automatisk testning

GitHub Actions-workflow: `.github/workflows/matlab-tests.yml`, körs vid varje push till `main` och manuellt via Actions. Fem MEX-oberoende testfiler ingår i den aktiva CI-testsviten och har passerat:

```matlab
setupPaths
files = {'tests/testDataLimit.m', ...
         'tests/testLineSearchMain.m', ...
         'tests/testMlPolyMain.m', ...
         'tests/testSelectPolyDegreeLr.m', ...
         'tests/test_denoise_quadratic_approximation.m'};
results = runtests(files);
assertSuccess(results)
```

**Ännu inte i CI:** `testAdCounts.m` (kräver kompatibel MEX-fil), `test_poly_fit_ml.m`, `test_mult_q_poly_fits.m`, `test_select_poly_degree.m` (äldre fristående skript som bör ses över och konverteras till `matlab.unittest`). Detta innebär inte att skripten är felaktiga, bara att de inte är del av den nya ordinarie sviten.

## Godkännandekriterier

En rad får byta till **Approved** när:

1. Beteendet är fastställt mot gällande krav/design (eller tydligt beskrivet avvikelsebeslut).
2. Relevanta enhets-/integrationstester har passerat på aktuell kod.
3. Väsentliga numeriska egenskaper och toleranser har verifierats.
4. Bevis finns med commit, testkommando, miljö, datum och resultatlänk.

### Verifieringslogg

| Datum | Funktion | Commit | Miljö / data | Testkörning | Slutsats |
|---|---|---|---|---|---|
| 2026-11-16 | Syntetisk linjesökning och denoising (begränsat testurval) | `61f05edf` | GitHub Actions, Ubuntu, MATLAB R2026b | [MATLAB Tests](https://github.com/ingemarhelander/kalmanrates/actions/runs/43884499822) | CI success; **inte** fullständig funktionell verifiering |
| 2026-11-16 | Utökad MATLAB-testsvit (5 filer) | `79214cf9` | Ubuntu, MATLAB R2026b | [Godkänd CI-körning](https://github.com/ingemarhelander/kalmanrates/actions/runs/43886139997) | Passed – ytterligare systemtestning återstår |

## Testplan framåt, i prioritetsordning

1. **Klart:** Utökad femfilerssvit körd och godkänd i CI.
   Resultatet dokumenterat i verifieringsloggen.
2. Identifiera testluckor i Req. 22–31 och följ upp [#38](https://github.com/ingemarhelander/kalmanrates/issues/38) vid ändrat stoppkriterium.
3. Komplettera enhetstester för CG/Gram–Schmidt och Kalman/RTS med kända facit.
4. Genomför integrationskontroll med riktiga EUR-data och jämför två linjesökningsmetoder från samma startpunkt.
5. Kör MEX/AD-regression på rätt Windows-version och så småningom hela EUR-historiken.

När en design ändras: markera `Design pending`, identifiera berörda tester och återgå till **In testing** tills ny verifiering är genomförd.
