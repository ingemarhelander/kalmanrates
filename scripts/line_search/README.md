# Updates från Oskar

## Workflow:
1. Generera mockdata
2. Generera massa kurvor i mult_q_poly_fits
    2.1. mult_q_poly_fits anropar i en loop poly_fit_ml som i sin tur faktiskt genererar ml-uppskattningarna
    2.2. mult_q_poly_fits anropar lag1_autocorr och residual_trend_test för att få mer data till ml-uppskattningarna
3. Skicka dessa kurvor till new_select_poly_degree, som väljer det bäst lämpade gradtalet givet kurvorna och 
   den tillhörande datan som mult_q_poly_fits skickar med sig.

## Ändringar som har gjorts:

### lr_statistic.m:
Ändrat lite - tror inte man behövde köra (m-q-2) som faktor, utan bara m.

### mult_q_poly_fits.m:
Genererar en ml-uppskattad kurva givet data F, S, alpha, qmax.
Gör en kurva per gradtal. Sparar alla kurvor i fits. Sparar även följande info:
    q (gradtal för just denna ml-uppskattning)
    fits{q}.sigmaF (standardavvikelsen på kurvan jämfört med den observerade mätdatans funktionsvärden. Beräknas i poly_fit_ml.m)
    fits{q}.sigmaG (standardavvikelsen på kurvan jämfört med den observerade mätdatans derivator. Beräknas i poly_fit_ml.m)
    fits{q}.rhoF (autokorrelationen mellan ml-kurvan och den 'verkliga' kurvan. Beräknas i lag1_autocorr.m)
    fits{q}.rhoG (autokorrelationen mellan ml-kurvans derivata och den 'verkliga' kurvans derivata. Beräknas i lag1_autocorr.m)
    fits{q}.pTrendF (ett p-värde som baserat på statistik säger om f:s residualer är korrelerade, alltså inte bara brus. Beräknas i residual_trend_test.m)
    fits{q}.pTrendG (ett p-värde som baserat på statistik säger om g:s residualer är korrelerade, alltså inte bara brus. Beräknas i residual_trend_test.m)
    fits{q}.iterations;
Denna funktion genererar massa data om varje gradtal-ml-uppskattning men lämnar sen över till andra skript

### new_select_poly_degree.m:
Denna rackare tar in samtliga ml-kurvor från mult_q_polt_fits och gör ETT SMART VAL av gradtal man ska gå vidare med.
Den har just nu två stycken metoder som den kan göra ETT SMART VAL utifrån.

#### Metod 1: Threshold
Välj det lägsta gradtal som uppfyller vissa villkor. Just nu ser de ut som följande:
%   Selects the lowest degree satisfying all diagnostic thresholds.
%   Right now: abs(rhoF) <=0.2
%              abs(rhoG) <=0.2
%              ptrendF >= 0.05
%              ptrendG >= 0.05
%              converged == true

#### Metod 2: Score
Ger varje gradtal ett score, baserat på följande fyra faktorer:
%   Scores the fits using residual variance, residual autocorrelation,
%   trend-test p-values and model complexity. The lowest score wins.
%              1. Residual variance (sigmaF, sigmaG)
%              2. Residual autocorrelation (rhoF, rhoG)
%              3. Trend-test p-values (pTrendF, pTrendG)
%              4. Model complexity (degree)

Varje faktor får dessutom en vikt, med hur stort dess bidrag är till totala scoret. Just nu är vikterna:
%   options.weights         = [0.35 0.30 0.25 0.10]
%       Weights correspond to variance, autocorrelation, trend and
%       complexity, respectively. They are normalized to sum to one.

### poly_fit_ml.m:
Genererar en ml-kurva för ett visst gradtal givet data F, S, alpha, q. 
Sparar massa info i fit.() som exempelvis varianser, antal iterationer, koefficienter osv

### tests/test_select_poly_degree.m
Testskript som gör att man kan testa select_poly_degree som standalone unit.
Genererar mockdata o skriver output så man kan se förbättringarna hos våra parametrar
Skriver också vilket gradtal den hade valt

### tests/test_mult_q_poly_fits.m
Testskript som gör att man kan testa mult_q_poly_fits som standalone unit.
Genererar mockdata o skriver output så man kan se förbättringarna hos våra parametrar

### tests/test_poly_fit_ml.m
Hmm mest en relik just nu. Den testar inte poly_fit_ml ens, utan den anropar mult_q_poly_fits.m.
Jag behåller dock filen för att man i framtiden kan vilja ha ett testskript som testar poly_fit_ml ensamt.