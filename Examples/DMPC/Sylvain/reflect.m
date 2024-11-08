function P_sym = reflect(P, Q)
    % Créer des variables pour le point symétrique
    %P_sym = sdpvar(size(P, 1), size(P, 2));
    
    % Calcul du point symétrique de P par rapport à Q
    P_sym = 2*P - Q;
    % d= 2;
    % PQ = Q-P;
    % P_sym = P - PQ;%*d/norm(PQ);

end