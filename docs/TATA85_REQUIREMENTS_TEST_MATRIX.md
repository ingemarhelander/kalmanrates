# TATA85 – systematisk krav- och testtäckningsanalys

> **Förslag att lägga i `docs/REQUIREMENTS_TEST_MATRIX.md`**. Första inventering den 16 november 2026. Underlag: `Kravspecifikation/requirements1_1.tex` (version **1.1, 2026-11-13**), `Designspecifikation/designspec.tex`, GitHub `main`, öppna och stängda GitHub Issues, `tests/` och senaste femfilers-CI. Ingen ny MATLAB-körning eller MEX-test har genomförts inom denna analys.

## Resultat och tolkning

- **54 numrerade krav** i v1.1, varav **4 är uttryckligen överstrukna** (#4, #5, #36 och #42). Dessa fyra är inte aktiva acceptanskriterier utan versionshistorik.
- **50 aktiva krav**. Preliminär klassning: **Approved: 1**, **In testing: 36**, **Not implemented: 13**. Inga krav märks Fix required utan reproducerat fel.
- Prioritet för de aktiva kraven: **P1: 34**, **P2: 12**, **P3: 4**. Statusarna nedan beskriver *belagd testtäckning*, inte en fullständig kvalitetssäkring.
- **Approved är medvetet restriktivt.** Ett passerat CI-test bekräftar ett begränsat testscenario; det räcker inte ensamt till full acceptans för hela kravområden.
- Kravspecifikation **v1.0 och v1.1 skiljer sig materiellt**: #52–54 är tillägg; #19–20 och #27 har ändrats; #4, #5, #36, #42 är överstrukna. Följ v1.1 tills en ännu nyare godkänd version finns.

## Viktiga iakttagelser före nya tester

1. **Täckningslucka trots stängt Issue:** [#36](https://github.com/ingemarhelander/kalmanrates/issues/36) anger ett färdigt `tests/testEvaluate.m`, men filen finns **inte i aktuell `main`** enligt GitHub-trädet. Följ upp filens merge/läge innan ny implementation.
2. **Design i rörelse:** [#38](https://github.com/ingemarhelander/kalmanrates/issues/38) föreslår att stoppkriterier ändras. Slutlig acceptans för #31/#54 bör vänta på beslut om förväntat beteende.
3. **Viktiga P1-krav utan dokumenterad end-to-end-verifiering:** #10, #17–28, #34, #43–44, #46–49, #52–54. Några har passerade syntetiska tester, men deras fulla syften är ännu inte verifierade.
4. **Ingen belagd multistart:** #52 kräver flera startpunkter och #35 minst 9/10 giltiga körningar; `kalmanOptimizeMain.m` använder ett sparat initialvärde och saknar belagd multistart-harness.
5. **Denoising-testet mäter förbättrad riktning**, inte direkt reduktion av brus i både f och g. Därmed täcker det bara en del av #43.
6. **Marknadsdata och MEX-filer ligger utanför Git**. GitHub CI (Ubuntu) kör fem MEX-oberoende filer, så filter/MEX, helmodell, full historik och out-of-sample kräver en särskild testmiljö.
7. **#24:** `find_step_length.m` beskriver sig som tillfällig grid-search (1000 punkter), så det är implementerad prototyp men inte ännu verifierad numerisk maximering.

## Fullständig kravmatris

Kolumnen **Status** är en konservativ första bedömning. **In testing** omfattar både belagd implementation utan full verifiering och funktioner där kvarvarande täckning måste kontrolleras. **Not implemented** betyder att den efterfrågade funktionen eller leveransen *inte har kunnat identifieras i granskade källor*, inte ett definitivt påstående om all gruppens lokala kod.

| Req. | Prio | Funktion / krav | Status | Implementation / belägg | Testbelägg | GitHub Issue | Nästa verifiering |
|---|:---:|---|---|---|---|---|---|
| **1** | 1 | Dokument följer anvisad standard | **In testing** | Projektets TeX-dokument | — | — | Dokumentgranskning mot kursens formalia |
| **2** | 1 | Användarvänligt gränssnitt | **In testing** | README och körskript | — | — | Definiera användartest: installera, konfigurera, starta, tolka fel |
| **3** | 1 | Stoppa vid konvergens | **In testing** | kalmanOptimizeMain.m: gradient, maxIter, stalls | — | [#19](https://github.com/ingemarhelander/kalmanrates/issues/19) | Testa alla stopporsaker på styrda funktioner och riktig modell |
| **4** | 1 | Gränssnitt likelihood + gradient [struket] | **Struket** | kalmanObj.m | MEX-kontroller utanför CI | — | Ingen separat acceptans enligt v1.1 |
| **5** | 2 | Generiskt objektivgränssnitt [struket] | **Struket** | Stubb i tests/stubs/kalmanObj.m | testLineSearchMain.m | — | Ingen separat acceptans enligt v1.1 |
| **6** | 1 | Valbar tidsperiod för marknadsdata | **In testing** | modelRunConfig.m, selectObservations.m | testDataLimit.m – CI | [#1](https://github.com/ingemarhelander/kalmanrates/issues/1) | Verifiera på riktiga data och tidsgränser |
| **7** | 1 | OIS-data från handledare | **In testing** | measurement/, data/raw/ | — | — | Spåra datakälla, konventioner, saknade värden och instrumentpriser |
| **8** | 3 | Jämföra antal PCA-komponenter | **Not implemented** | PCA-kod finns, studie ej belagd | PCA-del i testDataLimit.m | — | Parameterstudie: komponentantal kontra likelihood, runtime, stabilitet |
| **9** | 3 | UKF eller alternativt Kalmanfilter (prio 3) | **In testing** | clsKalmanIr.m, MEX-kod | MEX-regression (manuell) | — | Bekräfta vald filtervariant och testa mot liten referensmodell |
| **10** | 1 | Modifierad CG med Gram–Schmidt | **In testing** | find_search_direction.m | Inga riktade CI-tester | [#5](https://github.com/ingemarhelander/kalmanrates/issues/5), [#35](https://github.com/ingemarhelander/kalmanrates/issues/35) | Verifiera ascent, ortogonalitet, randskalning och degenerering |
| **11** | 2 | BFGS-benchmark med båda linjesökningar | **Not implemented** | Ingen BFGS-väg identifierad | — | — | Beslut om prioritet; implementera och jämför mot baslinje |
| **12** | 1 | Steepest descent utan denoising | **In testing** | kalmanOptimizeSD.m | Inga riktade tester | — | Jämför med båda linjesökningar; reproducera körning |
| **13** | 2 | Jämföra sökriktningsalgoritmer | **Not implemented** | Flera alternativ i kod; ingen studie belagd | — | [#5](https://github.com/ingemarhelander/kalmanrates/issues/5) | Jämför kvalitet, evalueringsantal, tid, brus på samma data |
| **14** | 3 | Alternativ BFGS (prio 3) | **Not implemented** | Ingen identifierad implementering | — | — | Besluta om det prioriteras eller skjuts upp |
| **15** | 3 | Filtrera f och g i denoising (prio 3) | **In testing** | denoise.m innehåller lokal anpassning | test_denoise_quadratic_approximation.m | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Klargör om denna filtrering avses och verifiera separat |
| **16** | 2 | Välja lämplig delrumsdimension | **In testing** | config.noisy_size i kalmanOptimizeMain.m | — | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Studie över dimension, runtime, robusthet, felreduktion |
| **17** | 1 | Välja lämpligt w (ändrat från sigma) | **In testing** | denoise.m och config | Test av riktning, ej parameterstudie | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Definiera w och acceptanskriterium enligt v1.1 |
| **18** | 1 | Projicera till/från delrum | **In testing** | B' * X och B * qsmall i denoise.m | Indirekt denoise-test | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Direkt round-trip-test för ortonormalt basrum |
| **19** | 1 | Lös eq. 96–97 med fixed-point | **In testing** | denoise.m använder parameter fixed_or_exact | Endast riktningstest | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Verifiera faktisk fixed-point-gren, residual och konvergens |
| **20** | 1 | Lös eq. 96–97 med invers/pseudoinvers | **In testing** | denoise.m: fixed_or_exact | Endast riktningstest | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Verifiera exakt lösning, rankbrist, residual |
| **21** | 1 | Sökriktning i delrum | **In testing** | denoise.m returnerar d_prime | test_denoise_quadratic_approximation.m – CI | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Kontrollera riktningsvinkel, ascent och randvillkor |
| **22** | 1 | Första extra punkt i linjesökning | **In testing** | initial_step_length.m | Indirekt testLineSearchMain.m | [#10](https://github.com/ingemarhelander/kalmanrates/issues/10), [#11](https://github.com/ingemarhelander/kalmanrates/issues/11) | Direkt test av slope=0, negativ slope, min/max alpha |
| **23** | 1 | ML-anpassa polynom till f och g | **In testing** | poly_fit_ml.m, ml_poly_main.m | testMlPolyMain.m – CI | [#10](https://github.com/ingemarhelander/kalmanrates/issues/10), [#17](https://github.com/ingemarhelander/kalmanrates/issues/17) | Oberoende facit, variansvikter, rankbrist, 2 punkter |
| **24** | 1 | Hitta polynommaximum i 1D | **In testing** | find_step_length.m (markerad tillfällig gridmetod) | Indirekt testLineSearchMain.m | [#10](https://github.com/ingemarhelander/kalmanrates/issues/10), [#16](https://github.com/ingemarhelander/kalmanrates/issues/16) | Känt analytiskt maximum, randmaximum, flera stationära punkter |
| **25** | 1 | Dynamiskt polynomgradval | **In testing** | select_poly_degree_lr.m, select_poly_degree_lag1.m | testMlPolyMain.m, testSelectPolyDegreeLr.m – CI | [#14](https://github.com/ingemarhelander/kalmanrates/issues/14), [#17](https://github.com/ingemarhelander/kalmanrates/issues/17) | Kontrollera metodval, precision, kalibrering, runtime |
| **26** | 1 | Residualer för f och riktningsderivata | **In testing** | poly_fit_ml.m: rF, rG | Indirekt testMlPolyMain.m | [#17](https://github.com/ingemarhelander/kalmanrates/issues/17) | Kontrollera residualer mot explicit beräkning |
| **27** | 1 | Statistiskt residualkriterium (reviderat) | **In testing** | LR och lag1-implementationer | testSelectPolyDegreeLr.m – CI | [#14](https://github.com/ingemarhelander/kalmanrates/issues/14), [#17](https://github.com/ingemarhelander/kalmanrates/issues/17) | Enas om primär metod och kontrollera falsklarm/förmåga |
| **28** | 1 | Fallback grad 2 vid få punkter | **In testing** | qMin=2 och fåpunktslogik | testMlPolyMain.m, testSelectPolyDegreeLr.m – CI | [#14](https://github.com/ingemarhelander/kalmanrates/issues/14), [#17](https://github.com/ingemarhelander/kalmanrates/issues/17) | Granska 2–3 punkter, singularitet, No-finite inputs |
| **29** | 2 | Bestäm tröskel för gradval | **In testing** | whiteness_threshold.m, LR kritisk gräns | Indirekt testMlPolyMain.m | [#14](https://github.com/ingemarhelander/kalmanrates/issues/14) | Verifiera statistisk nivå över frön och antal punkter |
| **30** | 2 | Hysteresis vid minskande grad | **In testing** | select_poly_degree_lag1.m, ej LR | Inga dedikerade CI-tester | [#15](https://github.com/ingemarhelander/kalmanrates/issues/15) | Sekvenstester med varierande data; avgör omfattning |
| **31** | 2 | Minsta avstånd mellan punkter | **In testing** | evaluate.m, line_search_main.m | #36 stängt men testEvaluate.m saknas i main | [#12](https://github.com/ingemarhelander/kalmanrates/issues/12), [#36](https://github.com/ingemarhelander/kalmanrates/issues/36), [#38](https://github.com/ingemarhelander/kalmanrates/issues/38) | Kontrollera #36-merge, uppdatera tester vid designändring |
| **32** | 2 | Minneskrav för fjärrdator | **In testing** | 250-dagarsläge; remote-instruktioner | Ingen mätning i CI | [#1](https://github.com/ingemarhelander/kalmanrates/issues/1) | Protokollera peak RAM och tid med fulla data på simhall |
| **33** | 2 | Dokumentation möjliggör vidareutveckling | **In testing** | README, docs/utveckling.md, kodoversikt | — | — | Praktiskt återuppbyggnadstest på ny miljö |
| **34** | 1 | Full optimering utan manuella steg | **In testing** | kalmanOptimizeMain.m | Endast separata testdelar | [#18](https://github.com/ingemarhelander/kalmanrates/issues/18), [#19](https://github.com/ingemarhelander/kalmanrates/issues/19) | Kör end-to-end till dokumenterad stoppreason |
| **35** | 2 | Minst 9 av 10 startpunkter fungerar | **Not implemented** | En startpunkt hårdkodad i setup_problem | — | [#19](https://github.com/ingemarhelander/kalmanrates/issues/19) | Definiera 10 startpunkter; 10 körningar och logg |
| **36** | 1 | Giltiga startpunkter [struket] | **Struket** | — | — | — | Ersatt/sammanfört med Req. 35 |
| **37** | 2 | Rapportera ogiltig indata | **In testing** | En del validering i selectObservations.m | testDataLimit.m täcker begränsningar | — | Negativa varianser, dimensioner, tomma dataset och NaN |
| **38** | 1 | Total tidsinsats 2160 h ±10% | **In testing** | Tidsrapport i Drive | Administrativ kontroll | — | Jämför rapporterad tid med krav |
| **39** | 1 | Leverans genom Git | **In testing** | Repo finns och används | — | — | Verifiera slutleverans/tagg och installationsanvisning |
| **40** | 1 | Slutpresentation | **Not implemented** | Ingen slutpresentation granskad | — | — | Administrativ leveranskontroll |
| **41** | 1 | Slutrapport | **Not implemented** | Ingen slutrapport granskad | — | — | Administrativ leveranskontroll |
| **42** | 1 | Ändliga värden för giltig indata [struket] | **Struket** | — | — | — | Följs praktiskt av Req. 35 |
| **43** | 1 | Denoising minskar numeriskt brus | **In testing** | denoise.m | Syntetiskt riktningstest i CI | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Mät faktiskt brus i f/g, inte enbart riktningsvinkel |
| **44** | 1 | Reproducerbarhet samma start och data | **In testing** | config, rng i enskilda tester | CI kör syntetiska prov | — | Repetera hela optimeringen och jämför parameter-/målskillnad |
| **45** | 1 | Git används för versionshantering | **Approved** | Aktivt Git-repository och commits | Git-historik | — | Bibehåll repo-spårning |
| **46** | 1 | Tester och definition of done per delsystem | **In testing** | tests/ och TEST_STATUS.md | 5 CI-testfiler, inte alla delsystem | — | Definiera DoD för prissättning, Kalman, optimering |
| **47** | 1 | Denoising på kända brusiga testfunktioner | **In testing** | test_denoise_quadratic_approximation.m | CI-passerad för 3 funktioner | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Fler frön, brusnivåer och mått för kvalitet |
| **48** | 1 | Jämförelse med AI-kodad referens | **Not implemented** | Inget motsvarande jämförelsetest identifierat | — | — | Avgränsa delsystem, oberoende referens, matchande indata |
| **49** | 1 | Utanför träningsperioden (out-of-sample) | **Not implemented** | PCA-urval finns, ej OOS-utvärdering | testDataLimit.m täcker endast PCA-kantfall | — | Definiera holdout, mått och läckagekontroll |
| **50** | 2 | Verifiera ekvation 86 | **Not implemented** | Inget riktat testsuite-test identifierat | — | — | Syntetiskt jämviktsresidualtest för y, Sbar, alpha |
| **51** | 2 | Verifiera ekvationer 96–97 | **Not implemented** | Ingen explicit residualassertion identifierad | — | [#2](https://github.com/ingemarhelander/kalmanrates/issues/2) | Testa stationaritet i båda lösningsmetoderna |
| **52** | 1 | Starta optimering från flera startpunkter | **Not implemented** | setup_problem använder första lagrade optimum | — | [#19](https://github.com/ingemarhelander/kalmanrates/issues/19) | Inför konfigurerad multistart före Req. 35 |
| **53** | 1 | Pre-filtrering med Kalman | **Not implemented** | Design anger planerad pre-filtrering | Inga identifierade tester | — | Avgränsa interface/algoritm och planera syntetiskt referenstest |
| **54** | 1 | Stoppkriterium för ny linjesökning | **In testing** | evaluate.m, line_search_main.m | Indirekt test; #36 stängt men testfil saknas | [#12](https://github.com/ingemarhelander/kalmanrates/issues/12), [#23](https://github.com/ingemarhelander/kalmanrates/issues/23), [#36](https://github.com/ingemarhelander/kalmanrates/issues/36), [#38](https://github.com/ingemarhelander/kalmanrates/issues/38) | Besluta aktuell stoppregel; återställ direkt testtäckning |

## CI-täckning och ännu ej inkopplade testfiler

| Kategori | Filer | Verifiering |
|---|---|---|
| I CI, passerade på commit `79214cf9` | `testDataLimit.m`, `testLineSearchMain.m`, `testMlPolyMain.m`, `testSelectPolyDegreeLr.m`, `test_denoise_quadratic_approximation.m` | [Godkänd GitHub Actions-körning](https://github.com/ingemarhelander/kalmanrates/actions/runs/43886139997) |
| MEX-beroende | `testAdCounts.m` | Kräver kompatibel Windows-MEX / separat runner |
| Äldre fristående skript (ej i CI) | `test_mult_q_poly_fits.m`, `test_poly_fit_ml.m`, `test_select_poly_degree.m` | Behöver kontrolleras/ev. konverteras till `matlab.unittest` |
| Manuella data-/MEX-kontroller | `checkLaptopSample.m`, `checkMexRegression.m` | Kräver faktiska `.mat`-filer och MEX |
| Gap i `main` | `testEvaluate.m` | Issue #36 anger klar, men fil saknas i aktuell trädlista |

## Föreslagna åtgärder – utan att börja koda tester

| Prio | Föreslagen åtgärd | Krav | Issue-hantering |
|---|---|---|---|
| **P0** | Bestäm ny officiell stoppregel, kontrollera var `testEvaluate.m` tog vägen och vilka gamla tester som fortfarande gäller | #31, #54 | Uppdatera **#36 och #38**, undvik nytt Issue innan läget klarnat |
| **P0** | Besluta hur multistart och kravet 9/10 ska operationaliseras, inklusive startpunkter och mål | #35, #52, #34 | Utöka befintligt #19 om relevant eller skapa *ett* tydligt systemtest-Issue |
| **P0** | Dokumentera DoD och testfacit för CG/Gram–Schmidt: ascent, ortogonalitet, bound feasibility, rankbrist | #10 | Skriv konkret acceptans i **#5 / #35** |
| **P0** | Kontrollera denoising:s två lösningsmetoder och residualer mot eq. 86/96/97; särskilj approximation från verklig brusreduktion | #17–21, #43, #50–51 | Komplettera **#2** med verifieringsdelar eller separat test-Issue |
| **P1** | Testplan för hela optimeringsloopen med verklig EUR-data och två linjesökningar | #3, #34, #44 | Länka **#18**; skriv mätbara acceptanskrav |
| **P1** | Granska statistiska gradvalsmekanismer, fastställ metod/kravnivå samt verifiera polynommaximum | #23–30 | Konsolidera **#14–17**; undvik dubblering |
| **P1** | Kontrollera vilken pre-filtreringsalgoritm som är tänkt, och om Req. 53 ligger inom årets scope | #53 | Skapa Issue *först om inget befintligt äger funktionen* |
| **P1** | Förbered OOS-upplägg: tidsuppdelning, läckage, utfallsmått och baslinje | #49 | Nytt Issue om ingen motsvarighet finns |
| **P2** | Dokumentera baseline mot AI-kodad lösning och referensdata | #48 | Eget validerings-Issue om det står kvar i kravspec |
| **P2** | Separera leverans-/formaliakrav från kodtestplan och ge dem ägare | #1–2, #32–33, #38–41, #45 | Projektledare/dokumentationsansvarig snarare än testansvarig |

## Föreslagen metod när designen ändras

1. Ändra **inte** godkännandestatus direkt baserat på en stängd Issue. Koppla först krav → kod → test → faktisk testlogg/commit.
2. Markera berörda krav som **Design pending** (notering) och behåll **In testing** tills beteendet är fastställt.
3. Dokumentera i befintligt Issue vilka testfall som blir ogiltiga och vilka som ska läggas till.
4. Vid ny implementering: definiera testfall, tolerans och reproducerbar indata innan funktionen får **Approved**.
5. För varje godkännande: logga commit, MATLAB/MEX-version, data, testkommando, resultatlänk och kvarvarande undantag i `docs/TEST_STATUS.md`.

## Avgränsningar / öppna frågor

- Detta är en **statisk källkods- och ärendegranskning**, inte en körning av MATLAB. En separat filinventering och riktade tester behövs för säkert godkännande.
- Version 1.1 i den uppladdade ZIP-filen används som norm. Bekräfta gärna med dokumentationsansvarig att den också är den senaste *godkända* kravbaslinjen. Google Drive innehåller även kommenterade kopior och utkast.
- För Req. 19–20 och 50–51 måste ekvationsreferenserna och avsedd numerisk tolerans bekräftas innan de markeras som godkända.
- Den här filen är en separat *analysmatris*. Behåll `docs/TEST_STATUS.md` som kortfattad operativ statusöversikt; länka hit därifrån för detaljer.
