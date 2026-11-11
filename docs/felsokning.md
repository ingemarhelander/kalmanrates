# Felsökning

Tillbaka till [startguiden](../README.md).

## MATLAB-version eller operativsystem skiljer sig

Det delade MEX-paketet är verifierat på 64-bitars Windows med MATLAB R2026b.
Kontrollera din miljö med `version('-release')` och `computer('arch')`.
Stäm av med Isak om du använder en annan MATLAB-version. MEX-filer byggda i en
nyare MATLAB-release är inte garanterade att fungera i äldre releaser.
Se [MathWorks om versionskompatibilitet](https://www.mathworks.com/help/matlab/matlab_external/version-compatibility.html).

Windows-filer med ändelsen `.mexw64` fungerar inte på Mac/Linux. Där behövs
plattformsspecifika MEX-filer, även för Jörgens `mexPortfolio`, vars källkod
inte finns i projektet.

## MATLAB hittar inte en fil

- Öppna projektmappen i MATLAB och kör `setupPaths`.
- Lägg `EUR.mat`, `EUR100.mat`, `kalmanDebug.mat` och `kalmanXoptAll.mat` i `data/`.
- Packa upp det uppdaterade ZIP-paketet från `TATA85/kod` i projektmappen.
  Behåll mapparna `mex/` och `tests/mex/` från paketet.
- `which mexKalmanLogLikelihood -all` ska visa filen i projektets egen `mex/`-mapp.

Git innehåller inte MAT-data eller kompilerade MEX-filer. Enbart `git pull`
hämtar därför inte dessa filer. Paketets `BUILDINFO.txt` innehåller byggmiljö
och filhashar. En lokal kopia av paketet kan finnas i `dist/` på Isaks dator;
paketet delas med gruppen via Drive.

## AD-varningen visas fortfarande

Den gamla modellbinären från Jörgen innehåller inte AD-fixen.
Hämta det uppdaterade paketet från Isak på Drive och ersätt motsvarande MEX-filer.
Kör `clear mexKalmanLogLikelihood testAdCountsMex` eller starta om MATLAB.
Kontrollera sedan sökvägen med `which mexKalmanLogLikelihood -all`.

## Saknade DLL-filer

Om MATLAB rapporterar saknade `MSVCP`-/`VCRUNTIME`-DLL-filer kan Microsoft
Visual C++ Redistributable för x64 behövas. Kollegorna behöver inte installera
en C++-kompilator eller Boost för att köra färdiga MEX-filer.
Se [MathWorks om saknade runtime-bibliotek](https://www.mathworks.com/help/matlab/matlab_external/invalid-mex-file-error.html).

## Saknade toolboxar

| Körning | Extra toolbox |
| --- | --- |
| De sju testerna, `kalmanOptimizeMain`, `kalmanOptimizeSD` | Ingen |
| `runReducedIrCurvesV2`, `checkLaptopSample` (RTS, använder `pdf`) | Statistics and Machine Learning Toolbox |
| `kalmanOptimizeMatlab` (använder `fmincon`) | Optimization Toolbox |

## Testerna stoppar

Alla tester ska passera innan du börjar optimera. Spara felmeddelandet och
kontakta Isak om något test misslyckas. Vanliga orsaker är saknad testhjälpare
eller MEX-filer för fel MATLAB-version.

## Omkompilering behövs

Isak ansvarar för bygget och det delade MEX-paketet. Nya byggen ingår normalt
inte i kollegornas arbetsflöde. Se [bygginstruktionerna](utveckling.md).

## Kort historik och minne

Ny modellinitiering med sex PCA-komponenter kräver minst åtta dagar före
träningsgränsen (2020-01-01 för EUR). Ett sparat Kalman-objekt kan däremot
begränsas utan att PCA görs om.

Även en kort körning kan behöva läsa in hela MAT-filerna innan urvalet görs.
MEX-kodens extra minnesreserv är kvar. Om laptopen får slut på minne, prova
färre dagar eller kör på simhall.
