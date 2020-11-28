function[results, SOGLIA_LLRs] = Simula_identifica_Byz_mod_def(sim_param)

Nprove = sim_param.Nprove;
PH1 = sim_param.PH1;
Pd_Hp = sim_param.Pd_Hp;
Pfa_Hp = sim_param.Pfa_Hp;
Pd_Bp = sim_param.Pd_Bp;
Pfa_Bp = sim_param.Pfa_Bp;
Pmal = sim_param.Pmal;
N = sim_param.N;
K = sim_param.K;
L = sim_param.L;
L2 = sim_param.L2;
T = sim_param.T;
gammas = sim_param.gammas;

%PH1 = 0.5;
PH0 = 1-PH1;
%Pd_Hp = 0.9;
%Pfa_Hp = 0.1;
%Pd_Bp = 0.9;
%Pfa_Bp = 0.1;
%Pmal = 0.5;

Pd_H = Pd_Hp;
Pfa_H = Pfa_Hp;
Pd_B = Pmal*(1-Pd_Bp)+(1-Pmal)*Pd_Bp;
Pfa_B = Pmal*(1-Pfa_Bp)+(1-Pmal)*Pfa_Bp;

%N = 7; %Numero di nodi totali
%K = 5; %Numero di nodi H
M = N-K; %Numero di nodi Bizantini
%L = 2; %Soglia per decidere ipotesi H1
%T = 20; %Finestra di osservazione per identificare i Bizantini
%gamma = 5; %Soglia sul numero di volte che una decisione e' diversa da quella globale per indentificare i Bizantini 

%Nprove = 1000; %Sarebbe il numero di simulazioni
Nsoglie_LLR = sim_param.Nsoglie_LLR;
Nerr_h = 0;
Nerr_b = 0;
N0=0;
N1=0;
Nerr_hr = zeros(length(gammas),1);
Nerr_br = zeros(length(gammas),1);
Nerr_hr_LLR = zeros(Nsoglie_LLR,1);
Nerr_br_LLR = zeros(Nsoglie_LLR,1);
Nerr_H = zeros(length(gammas),1);
Nerr_B = zeros(length(gammas),1);
Nerr_H_LLR = zeros(Nsoglie_LLR,1);
Nerr_B_LLR = zeros(Nsoglie_LLR,1);
for np = 1:Nprove
    if rem(np,500) == 0
        fprintf('Simulazione %d su %d\n',np,Nprove);
    end;
    rd = rand(1,T);
    P = zeros(1,T);
    P(rd < PH1) = 1;
    UH = zeros(K,T);
    UB = zeros(M,T);
    D = zeros(1,T);
    LLRs_OUT = zeros(N,T);
    for t = 1:T
        if P(t) == 1
            UH(:,t) = 1;
            GH = rand(K,1);
            UH(GH < 1-Pd_H,t) = 0;
            
            UB(:,t) = 1;
            GB = rand(M,1);
            UB(GB < 1-Pd_B,t) = 0;
        else
            GH = rand(K,1);
            UH(GH < Pfa_H,t) = 1;
            
            GB = rand(M,1);
            UB(GB < Pfa_B,t) = 1;
        end;
        %Costruzione delle log-likehood 
        Prob_err = Pmal*M/N;
        %Primo passo: calcolo le affidabilita' delle stime
        U_ALL = [UH(:,t);UB(:,t)];
        Num_ones = length(find(U_ALL == 1));
        Num_zeros = length(find(U_ALL == 0));
        alpha = M/N;
        P1 = (1-alpha)*Pfa_H+alpha*Pfa_B;
        P2 = (1-alpha)*Pd_H+alpha*Pd_B;
        PUH0 = ((1-Prob_err)*P1+Prob_err*(1-P1))^Num_ones * ((1-Prob_err)*(1-P1)+Prob_err*P1)^Num_zeros;
        PUH1 = ((1-Prob_err)*P2+Prob_err*(1-P2))^Num_ones * ((1-Prob_err)*(1-P2)+Prob_err*P2)^Num_zeros;
        for dec = 1:N
            if U_ALL(dec) == 0
                PUH0d = PUH0/(((1-Prob_err)*(1-P1)+Prob_err*P1));
                PUH1d = PUH1/((1-Prob_err)*(1-P2)+Prob_err*P2);
                Px0U = (1-Prob_err)*((1-P1)*(1-PH1)*PUH0d+(1-P2)*PH1*PUH1d);
                Px1U = Prob_err*(P1*(1-PH1)*PUH0d+P2*PH1*PUH1d);
            else
                PUH0d = PUH0/((1-Prob_err)*P1+Prob_err*(1-P1));
                PUH1d = PUH1/((1-Prob_err)*P2+Prob_err*(1-P2));
                Px0U = Prob_err*((1-P1)*(1-PH1)*PUH0d+(1-P2)*PH1*PUH1d);
                Px1U = (1-Prob_err)*(P1*(1-PH1)*PUH0d+P2*PH1*PUH1d);
            end;
            LLRs_OUT(dec,t) = abs(log(Px0U/Px1U));
        end;
        %Decisione
        if sum([UH(:,t);UB(:,t)]) >= L
            D(t) = 1;
        else
            D(t) = 0;
        end;
    end
    REL = sum(LLRs_OUT,2);
    if  (max(REL)- min(REL)) > 0 %Nsoglie_LLR > 0
        %SOGLIA_LLRs =
        %[min(REL)+(max(REL)-min(REL))/Nsoglie_LLR:(max(REL)-min(REL))/Nsoglie_LLR:max(REL)]; %da uno strano errore (sembra che non generi sempre un vettore di Nsoglie_LLR elementi)
        for i=1:Nsoglie_LLR
         SOGLIA_LLRs(i) = min(REL) + ((max(REL)-min(REL))/Nsoglie_LLR)*i;
        end
    else
        SOGLIA_LLRs = mean(REL)*ones(1,Nsoglie_LLR);
    end;
    %disp(SOGLIA_LLRs);
    for is = 1:Nsoglie_LLR
        SOGLIA_LLR = SOGLIA_LLRs(is);
        Nerr_H_LLR(is) = Nerr_H_LLR(is) + length(find(REL(1:K) < SOGLIA_LLR));
        Nerr_B_LLR(is) = Nerr_B_LLR(is) + length(find(REL(K+1:N) >= SOGLIA_LLR));
    end;
    Dall_H = repmat(D,K,1);
    Errs_H = xor(Dall_H,UH);
    Dall_B = repmat(D,M,1);
    Errs_B = xor(Dall_B,UB);
    eta_H = sum(Errs_H,2);
    eta_B = sum(Errs_B,2);
    for ig = 1:length(gammas)
        gamma = gammas(ig);
        Nerr_H(ig) = Nerr_H(ig) + length(find(eta_H > gamma));
        Nerr_B(ig) = Nerr_B(ig) + length(find(eta_B <= gamma));
    end;
    indx = find(P == 1);
    Nerr_h = Nerr_h + length(find(D(indx) == 0));
    indx = find(P == 0);
    Nerr_b = Nerr_b + length(find(D(indx) == 1));
    % variabili che servono nel caso non simmetrico
    N0 = N0 + numel(find(P==0));
    N1 = N1 + numel(find(P==1));
    
    %Valutazione delle prestazioni dopo rimozione
    for ig = 1:length(gammas)
        Dr = 0*D;
        gamma = gammas(ig);
        indxH = find(eta_H <= gamma);
        indxB = find(eta_B <= gamma);
        if length(indxH)+length(indxB) == 0
            indxH = 1:length(eta_H);
            indxB = 1:length(eta_B);
        end;        
        for t = 1:T
            %Decisione
            if sum([UH(indxH,t);UB(indxB,t)]) >= (length(indxH)+length(indxB))*(L2/100)
                Dr(t) = 1;
            else
                Dr(t) = 0;
            end;
        end
        indx = find(P == 1);
        Nerr_hr(ig) = Nerr_hr(ig) + length(find(Dr(indx) == 0));
        indx = find(P == 0);
        Nerr_br(ig) = Nerr_br(ig) + length(find(Dr(indx) == 1));
    end;
    %Valutazione delle prestazioni dopo rimozione nel caso LLR
    for is = 1:Nsoglie_LLR
        Dr = 0*D;
        SOGLIA_LLR = SOGLIA_LLRs(is);
        indxH = find(REL(1:K) >= SOGLIA_LLR);
        indxB = find(REL(K+1:N) >= SOGLIA_LLR);
        if length(indxH)+length(indxB) == 0
            indxH = 1:length(eta_H);
            indxB = 1:length(eta_B);
        end;
        for t = 1:T
            %%%%Decisione
            if sum([UH(indxH,t);UB(indxB,t)]) >= (length(indxH)+length(indxB))*(L2/100)
                Dr(t) = 1;
            else
                Dr(t) = 0;
            end;
%            %%%%%% Decisione diversa
%            if D(t) == 1
%                if sum([UH(indxH,t);UB(indxB,t)]) >= (length(indxH)+length(indxB))*(L2/100)
%                 Dr(t) = 1;
%             else
%                 Dr(t) = 0;
%             end;
%            else 
%                  if sum([UH(indxH,t);UB(indxB,t)]) <= (length(indxH)+length(indxB))*(L2/100)
%                 Dr(t) = 0;
%             else
%                 Dr(t) = 1;
%             end;
               
        end
        indx = find(P == 1);
        Nerr_hr_LLR(is) = Nerr_hr_LLR(is) + length(find(Dr(indx) == 0));
        indx = find(P == 0);
        Nerr_br_LLR(is) = Nerr_br_LLR(is) + length(find(Dr(indx) == 1));
    end;
end;

% if(Pd_H == 1 - Pfa_H)
results.PD = 1-Nerr_h/Nprove/T;
results.PFA = Nerr_b/Nprove/T;
for ig = 1:length(gammas)
    results.PDr(ig) = 1-Nerr_hr(ig)/Nprove/T;
    results.PFAr(ig) = Nerr_br(ig)/Nprove/T;
    results.PD_IDB(ig) = 1-Nerr_H(ig)/K/Nprove;
    results.PFA_IDB(ig) = Nerr_B(ig)/M/Nprove;
    results.P_ISO_H(ig) = Nerr_H(ig)/K/Nprove;
    results.P_ISO_B(ig) = 1 - Nerr_B(ig)/M/Nprove;
end;
for is = 1:Nsoglie_LLR
    results.PD_IDB_LLR(is) = 1-Nerr_H_LLR(is)/K/Nprove; % 1 - P_ISO^H 
    results.PFA_IDB_LLR(is) = Nerr_B_LLR(is)/M/Nprove;  % 1 - P_ISO^B = P_NONISO^B 
    results.P_ISO_H_LLR(is) = Nerr_H_LLR(is)/K/Nprove; % P_ISO^H 
    results.P_ISO_B_LLR(is) = 1 - Nerr_B_LLR(is)/M/Nprove;  % P_ISO^B 
    results.PDr_LLR(is) = 1-Nerr_hr_LLR(is)/Nprove/T;
    results.PFAr_LLR(is) = Nerr_br_LLR(is)/Nprove/T;
end;
% else
% results.PD = 1-Nerr_h/N1;
% results.PFA = Nerr_b/N0;
% for ig = 1:length(gammas)
%     results.PDr(ig) = 1-Nerr_hr(ig)/N1;
%     results.PFAr(ig) = Nerr_br(ig)/N0;
%     results.PD_IDB(ig) = 1-Nerr_H(ig)/K/Nprove;
%     results.PFA_IDB(ig) = Nerr_B(ig)/M/Nprove;
%     results.P_ISO_H(ig) = Nerr_H(ig)/K/Nprove;
%     results.P_ISO_B(ig) = 1 - Nerr_B(ig)/M/Nprove;
% end;
% for is = 1:Nsoglie_LLR
%     results.PD_IDB_LLR(is) = 1-Nerr_H_LLR(is)/K/Nprove; % 1 - P_ISO^H 
%     results.PFA_IDB_LLR(is) = Nerr_B_LLR(is)/M/Nprove;  % 1 - P_ISO^B = P_NONISO^B 
%     results.P_ISO_H_LLR(is) = Nerr_H_LLR(is)/K/Nprove; % P_ISO^H 
%     results.P_ISO_B_LLR(is) = 1 - Nerr_B_LLR(is)/M/Nprove;  % P_ISO^B 
%     results.PDr_LLR(is) = 1-Nerr_hr_LLR(is)/N1;
%     results.PFAr_LLR(is) = Nerr_br_LLR(is)/N0;
% end;
end
%end

    
    
    
   



