function alphas = initial_points_main(alpha1)
% initial_points_main - De fyra första steglängderna i line_search_main.
%
% alpha1 kommer från initial_step_length_main. Utöver den mäts 3*alpha1/4,
% alpha1/2 och alpha1/4. Med startpunkten (alpha = 0) blir det 5 punkter
% innan första polynomanpassningen.
%
% Varför 5 punkter:
%   - Statistiken: polynomet av grad 2 har 3 koefficienter och dess
%     derivata 2. Med färre än 4 punkter kan polynomet gå exakt genom alla
%     F (m = 3) eller alla lutningar (m = 2). ML-skattningen i poly_fit_ml
%     ger då den sortens sigma värdet 0 och oändlig vikt, och den andra
%     sortens mätningar ignoreras.
%   - Modellen: select_poly_degree tillåter högst grad m - 2. Med 5 punkter
%     kan LR-testet välja grad 3 redan i första anpassningen. Det behövs
%     när F inte är kvadratisk mellan 0 och alpha1, t.ex. när lutningen
%     ökar allt snabbare; ett andragradspolynom får då fel lutning i 0.
%
% Ordningen är fallande. line_search_main sänker alphaMax till en punkt
% som är sämre än F0; med fallande ordning ligger nästa förvalda punkt
% ändå alltid i [0, alphaMax].
%
% In:   alpha1  första steglängden, i (0, 1]
% Ut:   alphas  [alpha1, 3*alpha1/4, alpha1/2, alpha1/4], mäts i den ordningen

    alphas = alpha1*[1, 3/4, 1/2, 1/4];
end
