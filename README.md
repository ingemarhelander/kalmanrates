# Kalmanrates

Skattning av EUR-räntekurvan med Kalmanfilter. Koden kommer ursprungligen från Jörgen.
Gruppen kan köra 250 dagar på egna laptops och hela historiken på simhall.

**Isak ansvarar för det uppdaterade MEX-paketet som krävs för att köra koden på sin egen dator.** De andra använder
färdiga filer och behöver normalt inte installera C++ eller Boost.

## Förutsättningar

- MATLAB och Git samt åtkomst till projektets GitHub-repo och Google Drive.
- Det delade paketet är verifierat på **64-bitars Windows med MATLAB R2026b**.
  Har du annan version eller Mac/Linux? Isåfall är det oklart om det kommer funka, men jag hoppas det /Isak.
- Standardoptimeringen nedan kräver inga extra toolboxar.
  För RTS och `fmincon`, se [toolboxkrav](docs/felsokning.md#saknade-toolboxar).

## Engångssetup på din laptop

### 1. Hämta koden

Klona projektet från github.

Har du redan projektet, kör följande från din projektmapp i terminalen:

```powershell
git fetch origin
git switch issue-1-komprimera-data
git pull --ff-only origin issue-1-komprimera-data
```

### 2. Hämta filer från Google Drive

Allt finns i **`TATA85/kod`**:

- Lägg `EUR.mat`, `EUR100.mat`, `kalmanDebug.mat` och `kalmanXoptAll.mat` i `data/`.
- Packa upp **`kalmanrates-mex-win64-R2026b-20261111.zip`** direkt i projektmappen
  (mappen med denna README). Behåll mapparna i ZIP-filen och ersätt äldre MEX-filer.

Kontrollera att filerna ligger så här:

```text
mex/mexKalmanLogLikelihood.mexw64
mex/mexPortfolio.mexw64
tests/mex/testAdCountsMex.mexw64
```

**MAT-data och MEX-filer följer inte med Git.** Använd det uppdaterade ZIP-paketet;
de äldre fristående modellfilerna från Jörgen saknar AD-fixen.

### 3. Kontrollera installationen

Öppna projektmappen i MATLAB och kör:

```matlab
setupPaths
clear mexKalmanLogLikelihood testAdCountsMex
which mexKalmanLogLikelihood -all
results = runtests('tests');
assertSuccess(results)
```

`which` ska visa filen i projektets egen `mex/`-mapp och alla sju tester ska passera.
Om något stoppar: se [felsökning](docs/felsokning.md).

## Kör 250 dagar på din laptop

Behåll `config.maxObservations = 250` i `modelRunConfig.m` (standard).
När setup är klar räcker följande i MATLAB från projektmappen vid varje körning:

```matlab
setupPaths
kalmanOptimizeMain
```

## Kör hela historiken på simhall

1. [Anslut till simhall](docs/fjarrskrivbord.md).
2. Hämta samma branch, data och uppdaterade MEX-filer där. Stäm av MATLAB-versionen. Laptop-paketet är byggt för R2026b.
3. Ändra **lokalt på simhall** till `config.maxObservations = Inf` i `modelRunConfig.m`.
   Behåll 250 som gemensam standard i Git.
4. Kör `setupPaths` och `kalmanOptimizeMain` i MATLAB från projektmappen.

Inställningen gäller även de andra optimeringsskripten och `runReducedIrCurvesV2`.
Korta körningars resultatfiler får tillägget `_sample250` (eller valt antal dagar); fulla körningar använder
de ursprungliga resultatnamnen. Hela optimeringsloopen och hela historiken återstår att verifiera.

## Mer dokumentation

- [Felsökning och toolboxkrav](docs/felsokning.md)
- [Bygginformation för Isak och ändringar i Jörgens kod](docs/utveckling.md)
- [Projektstruktur, flödesschema och fillista](docs/kodoversikt.md)
- [Fjärrskrivbord och filöverföring till simhall](docs/fjarrskrivbord.md)
