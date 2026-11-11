# Utvecklingsdokumentation

Tillbaka till [startguiden](../README.md).
Den här sidan är för Isak och den som behöver förstå eller bygga om koden.

De ursprungliga MEX-filerna kommer från Jörgen. I det delade paketet är
`mexPortfolio` oförändrad. Kalman-filen har byggts om för AD-fixen och
`testAdCountsMex` är en ny testhjälpare.

## För Isak: bygginformation vid behov

Nya byggen ingår normalt inte i gruppens arbetsflöde: kollegorna använder de
färdiga filerna från Isak. Informationen nedan sparas om en framtida kodändring
eller annan MATLAB-version skulle kräva omkompilering. Isak ansvarar då för
bygget, som kräver en kompatibel C++-kompilator och Boost.
Det lokala bygget har verifierats med MSVC 2022 och Boost 1.67.0.
`compileMex` använder `BOOST_ROOT` om den miljövariabeln är satt, annars
`%USERPROFILE%/Documents/boost_1_67_0`. Annan sökväg kan anges som argument.
Spara referensen med den gamla MEX-filen **innan** den kompileras om:

```matlab
setupPaths
checkMexRegression('baseline') % Sparar gamla värden och gradienter
compileMex                    % Bygger den rättade MEX-filen och testhjälparen
which mexKalmanLogLikelihood -all
results = runtests('tests');
assertSuccess(results)
checkMexRegression('verify')   % Jämför gamla och nya resultat
```

Verifieringen följer `maxObservations`, alltså 250 dagar som standard.
Referensen sparas då i `data/adCountReference_sample250.mat`. Om den redan har
sparats med den gamla MEX-filen, hoppa över `baseline`. Den gamla AD-varningen kan visas
under referenskörningen, men ska försvinna efter omkompileringen.
Flera storlekar kan väljas med ett tredje argument, exempelvis
`checkMexRegression('baseline', dataPath('fullReference.mat'), [100 250 500 Inf])`.
Använd samma fil och storlekar med `'verify'`. Välj hela historiken på den starkare datorn.

Efter godkända kontroller delar Isak de tre MEX-filerna i [startguiden](../README.md), tillsammans
med MATLAB-release och plattform, via gruppens Drive. Övriga behöver då inte
installera Boost eller en C++-kompilator och behöver inte köra `compileMex`
eller jämförelsen mellan gamla och nya binärer.

Lokalt verifierat 2026-11-11 med MATLAB R2026b, MSVC 2022 och Boost 1.67.0:
alla sju MATLAB-tester passerar. För 250 dagar är målfunktion och gradient
identiska med den gamla binärens resultat, gradientkontrollen passerar och
AD-varningen är borta. Ett söksteg i den egna optimeringen höjer målfunktionen.
Hela optimeringsloopen och hela historiken är ännu inte verifierade.

## Ändringar i Jörgens ursprungliga kod

Den ursprungliga koden kommer från Jörgen. På `issue-1-komprimera-data` har följande
ändrats för att kunna välja antal dagar och göra mindre körningar på en laptop:

- **Gemensam inställning:** `modelRunConfig.m` styr antalet dagar för både
  modellinitiering och optimering. Ett positivt heltal väljer de första dagarna;
  `Inf` använder hela historiken. Standard är 250 dagar.
- **Tidigare dataurval:** dagarna väljs före PCA och instrumentberäkningarna.
  Kurvor och datum hålls ihop, och instrumentindex anpassas till instrumenten i
  urvalet. De gamla avstängda blocken för sent dataurval har tagits bort.
- **Sparade modeller:** `loadKalmanSample.m` och `limitKalmanObservations.m`
  kortar alla tidsberoende fält tillsammans. Parameterdimensionerna behålls så
  att befintliga sparade parametrar fortfarande passar modellen.
- **Kortare PCA-perioder:** `pcaIRS.m` hanterar urval utan data efter
  träningsperioden. Minst åtta dagar före träningsgränsen krävs fortfarande
  (2020-01-01 för EUR).
- **Rättad AD-räknare:** i `AD_As_bdx::count()` har `n * Av.size()` ersatts
  med `Av.size()`, eftersom varje lagrat matriselement ger en gradientpost.
  De fasta korrigeringarna `12323500` och `3799737` har tagits bort eftersom de
  bara passade en viss datamängd. Felutskriften visar nu faktiskt och beräknat
  antal separat, utan missvisande osignerade differenser. Minnesreserven behålls.
- **Uppdaterat MEX-bygge:** projektets lokala AD-header används så att fixen
  kommer med vid kompilering. Även `sparseMatrix.h` hämtas från projektet.
  `compileMex.m` bygger modellen och testhjälparen med en C++-kompilator och
  Boost-headerfiler. De externa XiMath-/XiFin-sökvägarna och länkningen till
  det gamla kompilerade Boost-biblioteket behövs inte för detta bygge.
  Pi-konstanten beräknas utan det kompilatorberoende makrot `M_PI`.
- **Separata resultat:** korta körningar får tillägget `_sample250` (eller
  valt antal dagar), så att resultaten från fulla körningar inte skrivs över.
- **Nya kontroller:** tester har lagts till för dataurval, PCA och AD-räknaren.
  `checkMexRegression.m` kan spara referensvärden med den gamla binären före
  omkompilering och jämföra målfunktion och gradient efteråt för valt antal
  dagar. Hela historiken kräver ett uttryckligt val. Utvalda gradienter
  kontrolleras även med finita differenser för 250 dagar.
