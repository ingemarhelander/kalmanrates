# Kodöversikt

Tillbaka till [startguiden](../README.md).
Filnamn och sökvägar nedan utgår från projektmappen.

## Struktur

```
kalmanrates/
├── setupPaths.m      lägger till projektmapparna på MATLAB-sökvägen
├── dataPath.m        sökväg till fil i data/, t.ex. dataPath('EUR100.mat')
├── modelRunConfig.m  antal observationer: 250 eller Inf för all data
├── scripts/          körskript, kontrollskript och functions/ för optimering
├── model/            Kalmanmodellen och funktioner för dataurval
├── mex/              C++/MEX-kod och kompilerade filer
├── measurement/      instrumentklasser och dataladdning
├── data/             .mat-filer (*), raw/ innehåller Excel-rådata
├── tests/            tester för dataurval och PCA
├── dist/             lokalt MEX-paket för delning, ignoreras av Git
└── docs/             dokument och figurer
```

## Hur filerna hänger ihop

Heldragna pilar = anrop, streckade pilar = läser/sparar data.
Vid korta körningar får resultatfilerna tillägget `_sample250` (eller valt antal dagar).

```mermaid
flowchart TD
    subgraph DATA["data/"]
        XLSX["raw/EUR.xlsx"]
        ECB["raw/datesEURECB.xlsx"]
        EURMAT[("EUR.mat")]
        EUR100[("EUR100.mat")]
        KDEBUG[("kalmanDebug.mat")]
        XOPT[("kalmanXoptAll.mat")]
        AXRTS[("axRtsOpt.mat<br/>eller axRtsOpt_sample250.mat")]
    end

    subgraph MEAS["measurement/"]
        CMFE["createMatFromExcel.m"]
        PPWD["populatePortfolioWithData.m"]
        CIC["createInstrumentsClasses.m"]
        CLS["cls*.m<br/>OIS, IBOR, IRS, FRA,<br/>FxSwap, FutureIR, ..."]
    end

    subgraph MEX["mex/"]
        MEXP[["mexPortfolio.mexw64"]]
        MEXK[["mexKalmanLogLikelihood.mexw64"]]
        CPP["mexKalmanLogLikelihood.cpp<br/>kalmanLogLikelihoodAd.h<br/>sparseMatrix.h"]
        COMP["compileMex.m"]
    end

    subgraph MODEL["model/"]
        MI["marketInfo.m"]
        PCA["pcaIRS.m"]
        KIR["clsKalmanIr.m"]
        TAY["taylorInstrPC.m"]
        KSUB["kalmanSub.m"]
        KOBJ["kalmanObj.m"]
        KOBJMIN["kalmanObjMin.m"]
        SELECT["selectObservations.m"]
        LIMIT["limitKalmanObservations.m"]
        LOAD["loadKalmanSample.m"]
    end

    subgraph SCRIPTS["scripts/"]
        TENOR["runTenorMultiIrCurves.m"]
        V2["runReducedIrCurvesV2.m<br/>huvudskript"]
        INIT["runReducedIrCurvesInit.m"]
        CLSS["runReducedIrCurvesCls.m"]
        OPTM["kalmanOptimizeMatlab.m<br/>fmincon"]
        OPTSD["kalmanOptimizeSD.m<br/>brantaste lutning"]
        OPTMAIN["kalmanOptimizeMain.m"]
        OPTFUN["functions/<br/>sökriktning och steglängd"]
        CHECK["checkLaptopSample.m"]
    end

    subgraph TESTS["tests/"]
        TEST["testDataLimit.m"]
    end

    CONFIG{{"modelRunConfig.m<br/>maxObservations"}}
    SETUP{{"setupPaths.m + dataPath.m"}}
    PDF["docs/ecbCutsHikes.pdf<br/>eller ecbCutsHikes_sample250.pdf"]

    XLSX -.-> CMFE -.->|sparar| EURMAT

    TENOR --> MI
    TENOR --> PPWD
    TENOR --> CIC
    TENOR -.->|sparar| EUR100

    EURMAT -.-> PPWD
    PPWD --> MEXP
    CIC --> CLS
    CIC --> MEXP
    CLS --> MEXP
    TENOR --> MEXP
    INIT --> MEXP

    V2 --> INIT
    V2 --> CLSS
    V2 --> MEXK
    V2 --> KSUB
    V2 -.->|sparar| AXRTS
    V2 -.->|sparar| PDF

    INIT --> MI
    INIT --> CONFIG
    INIT --> SELECT
    INIT --> PCA
    INIT --> PPWD
    INIT --> CIC
    EUR100 -.-> INIT
    ECB -.-> INIT

    CLSS --> KIR
    KIR --> TAY
    KSUB --> KIR

    KDEBUG -.-> LOAD
    LOAD --> CONFIG
    LOAD --> LIMIT --> SELECT
    OPTM --> LOAD
    OPTSD --> LOAD
    OPTMAIN --> LOAD
    XOPT -.-> OPTM
    XOPT -.-> OPTSD
    XOPT -.-> OPTMAIN
    OPTM --> KOBJMIN --> KOBJ
    OPTSD --> KOBJ
    OPTMAIN --> KOBJ
    OPTMAIN --> OPTFUN --> KOBJ
    KOBJ --> MEXK

    CHECK --> INIT
    CHECK --> CLSS
    CHECK --> KSUB
    CHECK --> LOAD
    CHECK --> KOBJ
    CHECK --> KOBJMIN
    TEST --> SELECT
    TEST --> LIMIT
    TEST --> PCA
    TEST --> KIR

    COMP -.->|kompilerar| CPP
    CPP -.-> MEXK

    SETUP -.-> SCRIPTS
    SETUP -.-> CMFE
```

## Filer

**Rotmappen/**

- `setupPaths.m` / `dataPath.m` – sökvägar till projektets mappar och data
- `modelRunConfig.m` – gemensam inställning för antalet observationer

**scripts/**

- `runTenorMultiIrCurves.m` – skattar forwardkurvor `f` och `π` per dag, sparar `data/EUR100.mat`
- `runReducedIrCurvesV2.m` – huvudskript: kör Init + Cls, RTS-iterationer, plottar
- `runReducedIrCurvesInit.m` – väljer dagar före PCA och hämtar instrumentdata för dessa dagar
- `runReducedIrCurvesCls.m` – startparametrar, skapar Kalman-objektet `K`
- `kalmanOptimizeMain.m` – ML-skattning med egen gradientmetod och linjesökning
- `kalmanOptimizeMatlab.m` – ML-skattning med `fmincon`
- `kalmanOptimizeSD.m` – ML-skattning med brantaste lutning
- `functions/` – hjälpfunktioner för sökriktning och linjesökning; `noisy_optimization.m` är en platshållare
- `checkLaptopSample.m` – kontrollerar modellen, ett RTS-steg och gradienter utan att spara resultat
- `checkMexRegression.m` – sparar referens före omkompilering och jämför värden och gradienter efteråt

De tre optimeringsskripten laddar det ursprungliga sparade Kalman-objektet via `loadKalmanSample.m`.

**model/**

- `clsKalmanIr.m` – Kalman-klassen (systemmatriser, filter, målfunktion, plottar)
- `kalmanSub.m` – ett RTS-steg (Riccati-rekursion)
- `taylorInstrPC.m` – instrumentpris och Jacobian för €STR/OIS
- `pcaIRS.m` – PCA av forwardkurvor, även när urvalet saknar en period utanför träningsdata
- `marketInfo.m` – marknadskonventioner per valuta
- `kalmanObj.m` / `kalmanObjMin.m` – målfunktion + gradient för optimering
- `selectObservations.m` – väljer de första observationerna och bestämmer resultatfilernas tillägg
- `limitKalmanObservations.m` – kortar ett sparat Kalman-objekt utan att ändra parameterdimensionerna
- `loadKalmanSample.m` – laddar sparad modell och tillämpar gränsen från `modelRunConfig.m`

**tests/**

- `testDataLimit.m` – testar dataurval, oförändrade parametrar och PCA med kort historik
- `testAdCounts.m` – kör räknartesterna i `mex/testAdCountsMex.cpp`, som byggs av `compileMex`

**mex/**

- `mexKalmanLogLikelihood.cpp` – MEX-gränssnitt för log-likelihood + gradient
- `kalmanLogLikelihoodAd.h` – Kalmanfilter med automatisk differentiering
- `sparseMatrix.h` – gles matrisklass
- `compileMex.m` – bygger modellen och testhjälparen med lokal AD-/matriskod och Boost-headerfiler
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
- `axRtsOpt_sample250.mat` – motsvarande resultat från en kort körning (antalet dagar styr namnet)
- `adCountReference_sample250.mat` – genererad referens för MEX-kontroll med 250 dagar (namnet följer urvalet)

**docs/**

- `utveckling.md` – bygginstruktioner, verifiering och ändringar i Jörgens kod
- `kodoversikt.md` – struktur, anropsschema och fillista
- `fjarrskrivbord.md` – anslutning och filöverföring till simhall
- `felsokning.md` – kompatibilitet, saknade filer och toolboxar

- `projektbeskrivning.pdf` – teori och ekvationer
- `ecbCutsHikes.pdf` – figur: marknadsimplicerad vs. ECB:s styrränta (skapas av `runReducedIrCurvesV2.m`)
- `ecbCutsHikes_sample250.pdf` – motsvarande figur från en kort körning

\* Ignoreras av Git. Hämta indata och delade MEX-filer från Drive enligt [startguiden](../README.md). Resultatfiler skapas av skripten.
