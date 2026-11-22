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

Testavgränsning: Testningen fokuserar på gruppens egna implementationer, förändringar av Jörgens ursprungskod och integrationen mellan dessa. Oförändrad ursprungskod kvalitetssäkras inte från grunden.

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
| Modifierad CG/sökriktning | **In testing** | [#5](https://github.com/ingemarhelander/kalmanrates/issues/5), `find_search_direction.m`, `testFindSearchDirection.m` (4 tester godkända i CI) | Följ upp Julianas kodstädning och komplettera vid behov med ytterligare specialfall. |
| Föreslagen Gram–Schmidt-ändring | **Not implemented** (ändringen ej verifierad) | [#35](https://github.com/ingemarhelander/kalmanrates/issues/35) | Implementera och skapa tester för (nära) linjärt beroende. |
| AD/MEX-gradient | **In testing** | `testAdCounts.m`, `checkMexRegression.m`, `docs/utveckling.md` | Kör på rätt Windows/MATLAB/MEX och kontrollera finit differens. |
| Full optimeringsloop och verkliga EUR-data | **In testing** | [#18](https://github.com/ingemarhelander/kalmanrates/issues/18), `kalmanOptimizeMain.m` | Kör fler iterationer; jämför backtracking och polynomlinjesökning. |

**Observera:** `Fix required` tilldelas först när ett faktiskt problem bekräftats. Omfattande krav-för-krav-kartläggning är nästa inventeringssteg.

## Automatisk testning


GitHub Actions-workflow: `.github/workflows/matlab-tests.yml`, körs vid varje push till `main` och manuellt via Actions.

Den aktuella CI-testsviten innehåller fem MEX-oberoende testfiler:

```matlab
setupPaths;
testFiles = {
    'tests/testDataLimit.m'
    'tests/testLineSearchMain.m'
    'tests/testSelectPolyDegree.m'
    'tests/test_denoise_quadratic_approximation.m'
    'tests/testFindSearchDirection.m'
};
results = runtests(testFiles);
disp(table(results));
assertSuccess(results);
```

**Senaste verifierade CI-resultat (2026-11-22):**
- Fyra tester för sökriktning/Gram–Schmidt (#5) passerade.
- Övriga tester passerade utom `testSelectPolyDegree/tooLowDegreeIsRaised`.
- Hela CI-körningen har därför status **Failed**.
- Felet i polynomgradvalet är reproducerat och utreds separat.

[Körningslogg](https://github.com/ingemarhelander/kalmanrates/actions/runs/45042917181)

**Inte i CI:** MEX-beroende tester som `testAdCounts.m`, vilka kräver en kompatibel testmiljö.


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
| 2026-11-16 | Linjesökningens stoppkriterier (`testEvaluate.m`), utökad sexfilerssvit | `84a475d` | GitHub Actions, Ubuntu, MATLAB R2026b | [Godkänd CI-körning](https://github.com/ingemarhelander/kalmanrates/actions/runs/43897123079) | Passed – direkta tester av stoppkriterierna godkända; designändring #38 återstår |
| 2026-11-22 | Sökriktning/Gram–Schmidt (#5), fyra tester | `794a49ef` | GitHub Actions, Ubuntu, MATLAB R2026b | [CI-körning](https://github.com/ingemarhelander/kalmanrates/actions/runs/45042917181) | 4/4 tester passerade. Hela CI är röd på grund av ett separat fel i polynomgradvalet. |

## Testplan framåt, i prioritetsordning

1. **Klart:** Utökad femfilerssvit körd och godkänd i CI.
   Resultatet dokumenterat i verifieringsloggen.
2. Identifiera testluckor i Req. 22–31 och följ upp [#38](https://github.com/ingemarhelander/kalmanrates/issues/38) vid ändrat stoppkriterium.
3. 3. Komplettera vid behov testerna för gruppens egna implementationer av CG/Gram–Schmidt och följ upp Issue #5.
4. Genomför integrationskontroll med riktiga EUR-data och jämför två linjesökningsmetoder från samma startpunkt.
5. Kör MEX/AD-regression på rätt Windows-version och så småningom hela EUR-historiken.

När en design ändras: markera `Design pending`, identifiera berörda tester och återgå till **In testing** tills ny verifiering är genomförd.
