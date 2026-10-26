# Kalmanrates

Skattning av EUR-räntekurvan (OIS) som en tillståndsmodell med Kalmanfilter. Forwardkurvan = PCA-faktorer + trappsteg vid ECB:s räntemöten. Parametrar skattas med maximum likelihood, gradienten beräknas med automatisk differentiering i C++.

## Kom igång

1. Klona ned projektet med `git clone`.
2. Hämta datafilerna (markerade med \*) separat från Google Drive (Projekt TATA85/Kod):
   - `.mat`-filerna läggs i `data/`
   - `.mexw64`-filerna läggs i `mex/`

   De ignoreras av git och ska inte laddas upp på GitHub pga storlek.
3. Kör skripten i `scripts/`, t.ex. `runReducedIrCurvesV2`. Varje skript anropar `setupPaths` automatiskt, som lägger alla projektmappar på MATLAB-sökvägen. Om du kör delar av koden manuellt, kör `setupPaths` först.

## Struktur

```
kalmanrates/
├── setupPaths.m     lägger till projektmapparna på MATLAB-sökvägen
├── dataPath.m       sökväg till fil i data/, t.ex. dataPath('EUR100.mat')
├── scripts/         skript som körs
├── model/           Kalmanmodellen
├── mex/             C++/MEX-kod och kompilerade filer
├── measurement/     instrumentklasser och dataladdning
├── data/            .mat-filer (*), raw/ innehåller Excel-rådata
└── docs/            dokument och figurer
```

## Filer

**scripts/**
- `runTenorMultiIrCurves.m` – skattar forwardkurvor `f` och `π` per dag, sparar `data/EUR100.mat`
- `runReducedIrCurvesV2.m` – huvudskript: kör Init + Cls, RTS-iterationer, plottar
- `runReducedIrCurvesInit.m` – PCA, ECB-datum och instrumentdata per dag
- `runReducedIrCurvesCls.m` – startparametrar, skapar Kalman-objektet `K`
- `kalmanOptimizeMatlab.m` – ML-skattning med `fmincon`
- `kalmanOptimizeSD.m` – ML-skattning med brantaste lutning

**model/**
- `clsKalmanIr.m` – Kalman-klassen (systemmatriser, filter, målfunktion, plottar)
- `kalmanSub.m` – ett RTS-steg (Riccati-rekursion)
- `taylorInstrPC.m` – instrumentpris och Jacobian för €STR/OIS
- `pcaIRS.m` – PCA av forwardkurvor
- `marketInfo.m` – marknadskonventioner per valuta
- `kalmanObj.m` / `kalmanObjMin.m` – målfunktion + gradient för optimering

**mex/**
- `mexKalmanLogLikelihood.cpp` – MEX-gränssnitt för log-likelihood + gradient
- `kalmanLogLikelihoodAd.h` – Kalmanfilter med automatisk differentiering
- `sparseMatrix.h` – gles matrisklass
- `compileMex.m` – kompilering (kräver handledarens bibliotek, körs från `mex/`)
- `mexKalmanLogLikelihood.mexw64` \* – kompilerad version
- `mexPortfolio.mexw64` \* – handledarens bibliotek för marknadsdata och instrument
- `old/kalmanLogLikelihoodAd2.h` – äldre version, används inte

**measurement/**
- `createMatFromExcel.m` – Excel (`data/raw/`) → `.mat` (`data/`)
- `populatePortfolioWithData.m` – laddar data i `mexPortfolio`
- `createInstrumentsClasses.m` – skapar instrumentobjekt för en dag
- `cls*.m` – instrumentklasser (OIS, IBOR, IRS, FRA, m.fl.)

**data/**
- `raw/EUR.xlsx` – rådata (BID, ASK, IBOR)
- `raw/datesEURECB.xlsx` – ECB:s mötesdatum
- `EUR.mat` \* – rådata i MATLAB-format
- `EUR100.mat` \* – forwardkurvor 2005–2025 (287 MB)
- `kalmanDebug.mat` \* – sparat Kalman-objekt
- `kalmanXoptAll.mat` \* – optimerade parametrar
- `axRtsOpt.mat` \* – starttillstånd från RTS

**docs/**
- `projektbeskrivning.pdf` – teori och ekvationer
- `ecbCutsHikes.pdf` – figur: marknadsimplicerad vs. ECB:s styrränta (skapas av `runReducedIrCurvesV2.m`)

\* Ignoreras av git, hämtas separat från Google Drive.
