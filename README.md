# Kalmanrates

Skattning av EUR-räntekurvan (OIS) som en tillståndsmodell med Kalmanfilter. Forwardkurvan = PCA-faktorer + trappsteg vid ECB:s räntemöten. Parametrar skattas med maximum likelihood, gradienten beräknas med automatisk differentiering i C++.

## Kom igång

1. Klona ned projektet med `git clone`.
2. Hämta datafilerna (markerade med \*) separat från Google Drive (Projekt TATA85/Kod) och lägg dem i rotkatalogen – de ignoreras av git och ska inte laddas upp på GitHub pga storlek.

## Filer

**Skript**
- `runTenorMultiIrCurves.m` – skattar forwardkurvor `f` och `π` per dag, sparar `EUR100.mat`
- `runReducedIrCurvesV2.m` – huvudskript: kör Init + Cls, RTS-iterationer, plottar
- `runReducedIrCurvesInit.m` – PCA, ECB-datum och instrumentdata per dag
- `runReducedIrCurvesCls.m` – startparametrar, skapar Kalman-objektet `K`
- `kalmanOptimizeMatlab.m` – ML-skattning med `fmincon`
- `kalmanOptimizeSD.m` – ML-skattning med brantaste lutning

**Modell**
- `clsKalmanIr.m` – Kalman-klassen (systemmatriser, filter, målfunktion, plottar)
- `kalmanSub.m` – ett RTS-steg (Riccati-rekursion)
- `taylorInstrPC.m` – instrumentpris och Jacobian för €STR/OIS
- `pcaIRS.m` – PCA av forwardkurvor
- `marketInfo.m` – marknadskonventioner per valuta
- `kalmanObj.m` / `kalmanObjMin.m` – målfunktion + gradient för optimering

**C++ / MEX**
- `mexKalmanLogLikelihood.cpp` – MEX-gränssnitt för log-likelihood + gradient
- `kalmanLogLikelihoodAd.h` – Kalmanfilter med automatisk differentiering
- `kalmanLogLikelihoodAd2.h` – äldre version, används inte
- `sparseMatrix.h` – gles matrisklass
- `compileMex.m` – kompilering (kräver handledarens bibliotek)
- `mexKalmanLogLikelihood.mexw64` \* – kompilerad version
- `mexPortfolio.mexw64` \* – handledarens bibliotek för marknadsdata och instrument

**measurement/**
- `EUR.xlsx` – rådata (BID, ASK, IBOR)
- `datesEURECB.xlsx` – ECB:s mötesdatum
- `createMatFromExcel.m` – Excel → `.mat`
- `populatePortfolioWithData.m` – laddar data i `mexPortfolio`
- `createInstrumentsClasses.m` – skapar instrumentobjekt för en dag
- `cls*.m` – instrumentklasser (OIS, IBOR, IRS, FRA, m.fl.)

**Data** \*
- `EUR.mat` – rådata i MATLAB-format
- `EUR100.mat` – forwardkurvor 2005–2025 (287 MB)
- `kalmanDebug.mat` – sparat Kalman-objekt
- `kalmanXoptAll.mat` – optimerade parametrar
- `axRtsOpt.mat` – starttillstånd från RTS

**Dokument**
- `projektbeskrivning.pdf` – teori och ekvationer
- `ecbCutsHikes.pdf` – figur: marknadsimplicerad vs. ECB:s styrränta

\* Ignoreras av git, hämtas separat.
