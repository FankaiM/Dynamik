%Aufraeumen
close all;
clear;

%Lade Sollbahn im Arbeitsraum in Workspace
load( 'Gelenkwinkelbahn.mat' );

%Erzeuge Roboterstruktur
rob = erstelle_roboter();

%% Berechnung der Lagerkraefte und Lagermomente

lagerkraft = zeros(3,rob.N_Q);
lagermoment = zeros(3,rob.N_Q);

lagerkraft_norm = zeros(rob.N_Q,length(T));
lagermoment_norm = zeros(rob.N_Q,length(T));

for z = 1:length(T)
    
    rob.q = Q(:,z);
    rob.dot_q = dot_Q(:,z);
    rob.ddot_q = ddot_Q(:,z);
    
    rob=berechne_dk_positionen_vektorkette(rob);
    rob=berechne_dk_geschwindigkeiten(rob);
    rob=berechne_dk_beschleunigungen(rob);
    
    for i = 1:rob.N_Q
        g = 7-i;
        schwerkraft = rob.kl(g).A_i0 * rob.kl(g).m * rob.B0_g;
        %Absolutbeschleunigung des Schwerpunkts:
        rob.kl(i).Bi_ddot_r_s = rob.kl(i).Bi_ddot_r_i... % Beschleunigung des Ursprunges
                               + tilde(rob.kl(i).Bi_omega)*tilde(rob.kl(i).Bi_omega)*rob.kl(i).Bi_r_s... % zentripetale beschleunigung
                               + tilde(rob.kl(i).Bi_dot_omega)*rob.kl(i).Bi_r_s; % tangentiale Beschleunigung
        if g == 6
            lagermoment(:,g) = rob.kl(g).I_o * rob.kl(g).Bi_dot_omega + tilde(rob.kl(g).Bi_omega)*rob.kl(g).I_o*rob.kl(g).Bi_omega...% Drallaenderung
                                + rob.kl(g).m*tilde(rob.kl(g).Bi_r_s)*rob.kl(g).Bi_ddot_r_i...
                                - tilde(rob.kl(g).Bi_r_s)*schwerkraft; % Moment der Schwerkraft
            lagerkraft(:,g) = rob.kl(g).m * rob.kl(g).Bi_ddot_r_s - schwerkraft;
        else
            lagermoment(:,g) = rob.kl(g).I_o * rob.kl(g).Bi_dot_omega + tilde(rob.kl(g).Bi_omega)*rob.kl(g).I_o*rob.kl(g).Bi_omega...% Drallaenderung
                                + rob.kl(g).m*tilde(rob.kl(g).Bi_r_s)*rob.kl(g).Bi_ddot_r_i...
                                - tilde(rob.kl(g).Bi_r_s)*schwerkraft...% Moment von Schwerkraft
                                - tilde(rob.kl(g+1).Bv_r_vi)*rob.kl(g+1).A_iv'*(-lagerkraft(:,g+1))...% Moment von Lagerkraft von naechsten Gelenk
                                - rob.kl(g+1).A_iv'*(-lagermoment(:,g+1));% Lagermoment von naechsten Gelenk
            lagerkraft(:,g) = rob.kl(g).m * rob.kl(g).Bi_ddot_r_s - schwerkraft-rob.kl(g+1).A_iv'*(-lagerkraft(:,g+1)); % - schwerkraft - Lagerkraft von naechsten Gelenk
        end
        lagermoment_norm(g,z) = sqrt(lagermoment(1,g)^2+lagermoment(2,g)^2+lagermoment(3,g)^2);
        lagerkraft_norm(g,z) = sqrt(lagerkraft(1,g)^2+lagerkraft(2,g)^2+lagerkraft(3,g)^2);
    end 
end

figure('Name','Kraftverlauf')
plot(T,lagerkraft_norm(1,:),...
     T,lagerkraft_norm(2,:),...
     T,lagerkraft_norm(3,:),...
     T,lagerkraft_norm(4,:),...
     T,lagerkraft_norm(5,:),...
     T,lagerkraft_norm(6,:));
h=legend( '$F_{q_1}$','$F_{q_2}$','$F_{q_3}$','$F_{q_4}$','$F_{q_5}$','$F_{q_6}$','Location','northwest');
h.Interpreter='latex';
xlabel( 't in [s]','Interpreter','latex');
ylabel( '$F_q$ in [N]','Interpreter','latex');
 
figure('Name','Momentverlauf')
plot(T,lagermoment_norm(1,:),...
     T,lagermoment_norm(2,:),...
     T,lagermoment_norm(3,:),...
     T,lagermoment_norm(4,:),...
     T,lagermoment_norm(5,:),...
     T,lagermoment_norm(6,:));
h=legend( '$M_{q_1}$','$M_{q_2}$','$M_{q_3}$','$M_{q_4}$','$M_{q_5}$','$M_{q_6}$','Location','northwest');
h.Interpreter='latex';
xlabel( 't in [s]','Interpreter','latex');
ylabel( '$M_q$ in [Nm]','Interpreter','latex');
