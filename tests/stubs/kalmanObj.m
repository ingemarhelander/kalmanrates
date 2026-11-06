function [f, grad] = kalmanObj(x, K, mu, varInd)
% kalmanObj (STUBB) - Syntetisk målfunktion som ersätter model/kalmanObj.m i testerna.
%
% Filen har samma namn och signatur som den riktiga kalmanObj. line_search
% anropar därför stubben utan att själv behöva ändras, och testerna kräver
% varken EUR100.mat eller mex-filer.
%
% VIKTIGT: Lägg inte tests/stubs i setupPaths. Då skulle stubben kunna ersätta
% den riktiga kalmanObj i vanliga körningar. Testklassen lägger själv till
% mappen, och tar bort den igen när testerna är klara.
%
% Målfunktionen är konkav med känt maximum i K.xStar. Med r = x - K.xStar är
%
%   f(x) = -0.5*r'*H*r - c*sum(r.^4)
%
% där H är positivt definit och c >= 0. Med c = 0 är funktionen kvadratisk.
% Med c > 0 kan längs en linje inget andragradspolynom beskriva den exakt,
% vilket behövs för att testa att polynomgraden anpassas efter datan.
%
% Eftersom K bara skickas vidare av line_search används K här för inställningar:
%   K.xStar         punkten där funktionen har sitt maximum (kolumnvektor)
%   K.H             positivt definit matris i den kvadratiska termen
%   K.quarticWeight vikten c i fjärdegradstermen
%   K.noiseStdF     standardavvikelse för bruset på funktionsvärdet
%   K.noiseStdG     standardavvikelse för bruset på varje gradientelement
%   K.nanRadius     f och grad blir NaN om norm(r) > nanRadius (Inf = aldrig)
%   K.infRadius     f blir +Inf om norm(r) > infRadius (Inf = aldrig)
%   K.lineProfile   tom för den konkava funktionen ovan, annars ett handtag
%                   [phi, dphi] = profile(s) för en icke-konkav funktion
%   K.lineDirection enhetsvektor u som hör till K.lineProfile
%
% Med K.lineProfile ersätts funktionen ovan av en icke-konkav funktion med
% s = u'*x och p = x - s*u (den del av x som är vinkelrät mot u):
%
%   f(x) = phi(s) - 0.5*p'*p
%
% Längs linjen x = alpha*u är p = 0, så där gäller f = phi(alpha). Testerna kan
% alltså välja exakt vilken form funktionen har längs sökriktningen.
%
% mu och varInd används inte. De finns bara för att signaturen ska stämma.

    r = x - K.xStar;

    % Brusfritt funktionsvärde och exakt gradient.
    if isempty(K.lineProfile)
        f = -0.5*(r'*K.H*r) - K.quarticWeight*sum(r.^4);
        grad = -K.H*r - 4*K.quarticWeight*r.^3;
    else
        u = K.lineDirection;
        s = u'*x;
        p = x - s*u;
        [phi, dphi] = K.lineProfile(s);
        f = phi - 0.5*(p'*p);
        grad = dphi*u - p;
    end

    % Bruset dras alltid, även när standardavvikelsen är noll. Då förbrukas
    % slumptalen på samma sätt oavsett inställning.
    f = f + K.noiseStdF*randn();
    grad = grad + K.noiseStdG*randn(size(x));

    % Efterliknar att Kalmanfiltret kan ge ogiltiga värden i vissa punkter.
    distance = norm(r);
    if distance > K.infRadius
        f = Inf;
    end
    if distance > K.nanRadius
        f = NaN;
        grad = NaN(size(x));
    end
end
