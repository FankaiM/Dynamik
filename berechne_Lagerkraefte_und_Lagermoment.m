%Aufraeumen
close all;
clear all;

%Lade Sollbahn im Arbeitsraum in Workspace
load( 'Gelenkwinkelbahn.mat' );

%Erzeuge Roboterstruktur
rob = erstelle_roboter();

%% Berechnung der Lagerkraefte und Lagermomente

lagerkraft = zeros(rob.N_Q,length(T),3);
lagermoment = zeros(rob.N_Q,length(T),3);

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
        if g == 6
            lagermoment(g,z,:) = rob.kl(g).I_o * rob.kl(g).Bi_dot_omega + tilde(rob.kl(g).Bi_omega)*rob.kl(g).I_o*rob.kl(g).Bi_omega...% Drallaenderung
                                - tilde(rob.kl(g).Bi_r_s)*schwerkraft; % Moment der Schwerkraft
            lagerkraft(g,z,:) = -schwerkraft;
        else
            lagermoment(g,z,:) = rob.kl(g).I_o * rob.kl(g).Bi_dot_omega + tilde(rob.kl(g).Bi_omega)*rob.kl(g).I_o*rob.kl(g).Bi_omega...% Drallaenderung
                                - tilde(rob.kl(g).Bi_r_s)*schwerkraft...% Moment von Schwerkraft
                                - tilde(rob.kl(g+1).Bv_r_vi)*rob.kl(g+1).A_iv'*[lagerkraft(g+1,z,1);lagerkraft(g+1,z,2);lagerkraft(g+1,z,3)]... % Moment von Lagerkraft von naechsten Gelenk
                                - rob.kl(g+1).A_iv'*[lagermoment(g+1,z,1);lagermoment(g+1,z,2);lagermoment(g+1,z,3)]; % Lagermoment von naechsten Gelenk
            lagerkraft(g,z,:) = -schwerkraft-rob.kl(g+1).A_iv'*[lagerkraft(g+1,z,1);lagerkraft(g+1,z,2);lagerkraft(g+1,z,3)]; % - schwerkraft - Lagerkraft von naechsten Gelenk
        end
    end 
end

figure('Name','Lagerkraftverlauf');
subplot(2,3,1);
plot(T, lagerkraft(1,:,1),...
     T, lagerkraft(1,:,2),...
     T, lagerkraft(1,:,3));
subplot(2,3,2);
plot(T, lagerkraft(2,:,1),...
     T, lagerkraft(2,:,2),...
     T, lagerkraft(2,:,3));
subplot(2,3,3);
plot(T, lagerkraft(3,:,1),...
     T, lagerkraft(3,:,2),...
     T, lagerkraft(3,:,3));
subplot(2,3,4);
plot(T, lagerkraft(4,:,1),...
     T, lagerkraft(4,:,2),...
     T, lagerkraft(4,:,3));
subplot(2,3,5);
plot(T, lagerkraft(5,:,1),...
     T, lagerkraft(5,:,2),...
     T, lagerkraft(5,:,3));
subplot(2,3,6)
plot(T, lagerkraft(6,:,1),...
     T, lagerkraft(6,:,2),...
     T, lagerkraft(6,:,3));
legend( 'F_x','F_y','F_z','Location','northwest');
xlabel( 't in [s]');
ylabel( 'F in [kg]');

  
figure();
plot( T, lagermoment(1,:,3), ...
      T, lagermoment(2,:,3), ...
      T, lagermoment(3,:,3), ...
      T, lagermoment(4,:,3), ...
      T, lagermoment(5,:,3), ...
      T, lagermoment(6,:,3) );
h=legend( '$M_{q_1}$','$M_{q_2}$','$M_{q_3}$','$M_{q_4}$','$M_{q_5}$','$M_{q_6}$','Location','southeast');
h.Interpreter='latex';
xlabel( 't in [s]','Interpreter','latex');
ylabel( '$M_q$ in [Nm]','Interpreter','latex');
