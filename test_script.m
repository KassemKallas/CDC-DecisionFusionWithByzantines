


sim_param.PH1 = 0.5;  % leave as it is



sim_param.N = 100;
sim_param.K = 60;
sim_param.T = 10;
sim_param.gammas = 0:sim_param.T;
sim_param.Nprove = 50000;

%folder = 'Nprove_50000\N=100\'; %'div_TH_AR\N=100\';
folder = 'Trials_for_SPL_50000_N=100_K=60_T=10\';


for TH= 30:10:70
    sim_param.L = TH;

    for TH_AR = 30:10:70 % questa � sempre in percentuale (TH_AR/100)
    sim_param.L2 = TH_AR;

    for Pd = 80:10:90

        sim_param.Pd_Hp = Pd/100;%detection at honest
        sim_param.Pfa_Hp = (100 - Pd)/100;%fa at honest
        sim_param.Pd_Bp = sim_param.Pd_Hp;%detection at Byzantine
        sim_param.Pfa_Bp = sim_param.Pfa_Hp;%false alarm at byzantine

        Pfa = 100 - Pd;

        mkdir([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR)]);
        save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\parametersTEST.mat'], 'sim_param');

        %-------------------parameters for LLR scheme-------------------------

        sim_param.Nsoglie_LLR =  6; % gam_var; % different LLR TH
        N_pmal= 6; %cz now we will have only 5 values 0.6 to 1.0
        Pmal=zeros(1,N_pmal);


        i_pmal=1;
        for p_mal_dec=5:1:6 + N_pmal - 1 %10 %0.6 to 1.0

            sim_param.Pmal = p_mal_dec/10;
            Pmal(i_pmal) = p_mal_dec/10;

            fprintf('Pmal LLR = %f\n',sim_param.Pmal);

            %--------------------------------------------------------------------------

            [results, SOGLIA_LLRs] = Simula_identifica_Byz_mod_def(sim_param);
            PD_IDB = results.PD_IDB;
            PFA_IDB = results.PFA_IDB;
            P_ISO_H= results.P_ISO_H; %ISO VALUE OH HONEST AT VARSHNEY SCHEME
            P_ISO_B= results.P_ISO_B;%ISO VALUE OF BYZANTINES AT VARSHNEY
            PD_IDB_LLR = results.PD_IDB_LLR;
            PFA_IDB_LLR = results.PFA_IDB_LLR;
            P_ISO_H_LLR = results.P_ISO_H_LLR;
            P_ISO_B_LLR = results.P_ISO_B_LLR;
            PFA = results.PFA;
            PD = results.PD;
            PFAr = results.PFAr; % Varshney
            PDr = results.PDr;
            PFAr_LLR = results.PFAr_LLR; % LLR
            PDr_LLR = results.PDr_LLR;
            PH1 = sim_param.PH1;


            % WITHOUT CHOOSING THE THRESHOLD WHISH GIVES THE MINIMUM

            PERR = PFAr + 1-PDr;
            PERR_LLR = PFAr_LLR+1-PDr_LLR;
%
%             % per il caso Pd \neq 1 - Pfa � necessario??
%             PERR = (1 - PH1)*PFAr + PH1*(1-PDr);
%             PERR_LLR = (1 - PH1)*(PFAr_LLR)+ PH1*(1-PDr_LLR);



            save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\resultsTEST_Pmal=',num2str(Pmal(i_pmal)),'.mat'], 'results', 'PERR', 'SOGLIA_LLRs', 'PERR_LLR');



            %plot(1-PD_IDB,1-PFA_IDB);title('ROC : PB^{iso} vs PH^{iso}');grid;hold;plot(1-PD_IDB_LLR,1-PFA_IDB_LLR,'*');
            %Prestazioni dopo rimozione
            %PERR = min(PFAr+1-PDr);

            %---------------------NOTE:=========================================
            %  hnest iso = 1-PD_IDB_LLR and Byzantine iso=1-PFA_IDB_LLR

            % honest_iso_llr(i_pmal)= 1-PD_IDB_LLR; %honest_iso_llr(i_pmal,i_gam)= 1-PD_IDB_LLR(gam_var);
            %
            % byzan_iso_llr(i_pmal)=1-PFA_IDB_LLR;

            honest_iso_llr{1,i_pmal}= P_ISO_H_LLR; %ISO value of honest at LLR scheme
            %
            byzan_iso_llr{1,i_pmal}=P_ISO_B_LLR;% ISO VALUE OF BYZANTINES AT LLR




            %
            % PFA_before{1,i_pmal}=PFA;
            % PD_before{1,i_pmal}=PD;


            payoff_LLR{1,i_pmal}= PERR_LLR;
            payoff{1,i_pmal}= PERR;
            soglia{1,i_pmal} = SOGLIA_LLRs;% the point of retrieving LLR threshold.
            %fprintf('Perr minima per lo schema Varshney = %f\n',PERR);

            varsh_iso_ho{1,i_pmal}=P_ISO_H;
            varsh_iso_by{1,i_pmal}=P_ISO_B;

            llr_iso_ho{1,i_pmal}=P_ISO_H_LLR;
            llr_iso_by{1,i_pmal}=P_ISO_B_LLR;


            fprintf('Perr per lo schema LLR = %f\n',PERR_LLR);

            i_pmal=i_pmal+1;
        end
        %i_gam= i_gam+1;
        %end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GAME matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%% Varshney scheme %%%%%

        GI= [];
        GI_PISO = [];
        GI_ISO_ratioBH = [];

        for i=1:1:N_pmal
            GI(:,i) = payoff{1,i}';   %M(:) mette tutti gli elementi in fila

            GI_PISO(:,i) = (varsh_iso_ho{1,i} - varsh_iso_by{1,i})';
%             if (varsh_iso_ho{1,i} == zeros(1,sim_param.T +1) && varsh_iso_by{1,i} == zeros(1,sim_param.T +1)) % NON MI ASPETTO DI DOVER GESTIRE IL CASO DI SOLO DENOMINATORE = 0 PERCHE? I BIZANTINI RIMOSSI NON SARANNO MAI MENO DEGLI ONESTI RiMOSSI
%                GI_ISO_ratioBH(:,i) = 1;
%              else
               GI_ISO_ratioBH(:,i) = (varsh_iso_by{1,i}./varsh_iso_ho{1,i})';
        %end
           % i=i+1;
        end

        save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\Vgame.mat'],'GI');

         %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\V_Honest_Pr_iso.mat'],'GI_HONEST_ISO');

         %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\V_Byz_Pr_iso.mat'],'GI_BYZ_ISO');

        save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_Pr_iso.mat'],'GI_PISO');
        save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_iso_ratio.mat'],'GI_ISO_ratioBH');
%

        %%%%% LLR scheme %%%%%

        GII= [];
%         GII_HONEST_ISO=[];
%         GII_BYZ_ISO=[];
        GII_PISO = [];
        GII_ISO_ratioBH = [];

        for i=1:N_pmal
            GII(:,i) = payoff_LLR{1,i}';   %M(:) mette tutti gli elementi in fila
%             GII_HONEST_ISO(:,i)=llr_iso_ho{1,i}';
%             GII_BYZ_ISO(:,i)=llr_iso_by{1,i}';
            GII_PISO(:,i) = (llr_iso_ho{1,i} - llr_iso_by{1,i})';
%             if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0)
%               GII_ISO_ratioBH(:,i) = 1;
%             else
              GII_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
            %end
            i=i+1;
        end

        save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLRgame.mat'],'GII');

         %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');

         %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');

        save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR), '\LLR_Pr_iso.mat'],'GII_PISO');
        save([folder,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR), '\LLR_iso_ratio.mat'], 'GII_ISO_ratioBH');


        %%%%%%%%%%%%section to save ISOLATION PROBAS%%%%%%%%%%%%%%%%%%%%%%



    end
    end
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% fine parte principale %%%%%%%%%%%%%%%%%%%
%
%
% %
% %
% clear all
%
% %%%%%%%%%%%%%%%%%%%%% altre prove in cui vario M e K %%%%%%%%%%%%%%%%%% %%%%%%%%%%%%
% % (da inserire in un nuovo for)
%
% sim_param.K = 54;
%
% folder1 = 'Nprove_50000\N=100\MoreOnest_K=54\'; % 'div_TH_AR\N=100\MoreOnest_K=54\';
%
% % li riscrivo dopo il clear all
% sim_param.PH1 = 0.5;  % leave as it is
% sim_param.N = 100;
% sim_param.T = 4;
% sim_param.gammas = 0:sim_param.T;
%
%
%
% sim_param.Nprove = 50000;
%
%
% for TH=  50 %30:10:70
%     sim_param.L = TH;
%
%     for TH_AR = 30:10:70
%     sim_param.L2 = TH_AR;
%
%     for Pd = 80:10:90
%
%         sim_param.Pd_Hp = Pd/100;%detection at honest
%         sim_param.Pfa_Hp = (100 - Pd)/100;%fa at honest
%         sim_param.Pd_Bp = sim_param.Pd_Hp;%detection at Byzantine
%         sim_param.Pfa_Bp = sim_param.Pfa_Hp;%false alarm at byzantine
%
%         Pfa = 100 - Pd;
%
%         mkdir([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR)]);
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\parametersTEST.mat'], 'sim_param');
%
%
%         sim_param.Nsoglie_LLR =  6; % gam_var; % different LLR TH
%         N_pmal= 5;
%         Pmal=zeros(1,N_pmal);
%
%
%         i_pmal=1;
%         for p_mal_dec=6:1:6 + N_pmal - 1 %10
%
%             sim_param.Pmal = p_mal_dec/10;
%             Pmal(i_pmal) = p_mal_dec/10;
%
%             %fprintf('Threshold LLR Gamma = %f\n',gam_var);
%             fprintf('Pmal LLR = %f\n',sim_param.Pmal);
%
%             %--------------------------------------------------------------------------
%
%             [results, SOGLIA_LLRs] = Simula_identifica_Byz_mod_def(sim_param);
%             % PD_IDB = results.PD_IDB;
%             % PFA_IDB = results.PFA_IDB;
%             P_ISO_H= results.P_ISO_H;
%             P_ISO_B= results.P_ISO_B;
%             % PD_IDB_LLR = results.PD_IDB_LLR;
%             % PFA_IDB_LLR = results.PFA_IDB_LLR;
%             P_ISO_H_LLR = results.P_ISO_H_LLR;
%             P_ISO_B_LLR = results.P_ISO_B_LLR;
%             PFA = results.PFA;
%             PD = results.PD;
%             PFAr = results.PFAr; % Varshney
%             PDr = results.PDr;
%             PFAr_LLR = results.PFAr_LLR; % LLR
%             PDr_LLR = results.PDr_LLR;
%
%
%             % WITHOUT CHOOSING THE THRESHOLD WHISH GIVES THE MINIMUM
%
%             PERR = PFAr + 1-PDr;
%             PERR_LLR = PFAr_LLR+1-PDr_LLR;
%
%
%
%             save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\resultsTEST_Pmal=',num2str(Pmal(i_pmal)),'.mat'], 'results', 'PERR', 'SOGLIA_LLRs', 'PERR_LLR');
%
%             honest_iso_llr{1,i_pmal}= P_ISO_H_LLR; %honest_iso_llr(i_pmal,i_gam)= 1-PD_IDB_LLR(gam_var);
%             %
%             byzan_iso_llr{1,i_pmal}=P_ISO_B_LLR;
%
%             payoff_LLR{1,i_pmal}= PERR_LLR;
%             payoff{1,i_pmal}= PERR;
%             soglia{1,i_pmal} = SOGLIA_LLRs;
%
%             varsh_iso_ho{1,i_pmal}=P_ISO_H;
%             varsh_iso_by{1,i_pmal}=P_ISO_B;
%
%             llr_iso_ho{1,i_pmal}=P_ISO_H_LLR;
%             llr_iso_by{1,i_pmal}=P_ISO_B_LLR;
%
%             %fprintf('Perr minima per lo schema Varshney = %f\n',PERR);
%             fprintf('Perr per lo schema LLR = %f\n',PERR_LLR);
%
%             i_pmal=i_pmal+1;
%         end
%
%
%
%         %%%%% Varshney scheme %%%%%
%
%         GI= [];
%         GI_PISO = [];
%         GI_ISO_ratioBH = [];
%
%
%         for i=1:N_pmal
%             GI(:,i) = payoff{1,i}';   %M(:) mette tutti gli elementi in fila
% %             GI_HONEST_ISO(:,i)=varsh_iso_ho{1,i}';
% %             GI_BYZ_ISO(:,i)= varsh_iso_by{1,i}';
%             GI_PISO(:,i) = (varsh_iso_ho{1,i} - varsh_iso_by{1,i})';
% %               if (varsh_iso_ho{1,i} == 0 && varsh_iso_by{1,i} == 0) % NON MI ASPETTO DI DOVER GESTIRE IL CASO DI SOLO DENOMINATORE = 0 PERCHE? I BIZANTINI RIMOSSI NON SARANNO MAI MENO DEGLI ONESTI RiMOSSI
% %                GI_ISO_ratioBH(:,i) = 1;
% %              else
%                GI_ISO_ratioBH(:,i) = (varsh_iso_by{1,i}./varsh_iso_ho{1,i})';
%             %end
%             i=i+1;
%         end
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\Vgame.mat'],'GI');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_Pr_iso.mat'],'GI_PISO');
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_iso_ratio.mat'],'GI_ISO_ratioBH');
%
%
%
%         %%%%% LLR scheme %%%%%
%
%         GII= [];
% %         GII_HONEST_ISO=[];
% %         GII_BYZ_ISO=[];
%         GII_PISO = [];
%         GII_ISO_ratioBH = [];
%
%         for i=1:N_pmal
%             GII(:,i) = payoff_LLR{1,i}';   %M(:) mette tutti gli elementi in fila
% %             GII_HONEST_ISO(:,i)=llr_iso_ho{1,i}';
% %             GII_BYZ_ISO(:,i)=llr_iso_by{1,i}';
%             GII_PISO(:,i) = (llr_iso_ho{1,i} - llr_iso_by{1,i})';
% %             if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0)
% %               GII_ISO_ratioBH(:,i) = 1;
% %             else
%               GII_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
%             %end
%             i=i+1;
%         end
%
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLRgame.mat'],'GII');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLR_Pr_iso.mat'],'GII_PISO');
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLR_iso_ratio.mat'],'GII_ISO_ratioBH');
%
%
%
%
%
%     end
% end
% end
%
%
%
%
% %
% %
% clear all
%
% %%%%%%%%%%%%%%%%%%%%% altre prove in cui vario M e K %%%%%%%%%%%%%%%%%% %%%%%%%%%%%%
% % (da inserire in un nuovo for)
%
% sim_param.K = 52;
%
% folder1 = 'Nprove_50000\N=100\MoreOnest_K=52\'; % 'div_TH_AR\N=100\MoreOnest_K=54\';
%
% % li riscrivo dopo il clear all
% sim_param.PH1 = 0.5;  % leave as it is
% sim_param.N = 100;
% sim_param.T = 4;
% sim_param.gammas = 0:sim_param.T;
%
%
%
% sim_param.Nprove = 50000;
%
%
% for TH=  50 %30:10:70
%     sim_param.L = TH;
%
%     for TH_AR = 30:10:70
%     sim_param.L2 = TH_AR;
%
%     for Pd = 80:10:90
%
%         sim_param.Pd_Hp = Pd/100;%detection at honest
%         sim_param.Pfa_Hp = (100 - Pd)/100;%fa at honest
%         sim_param.Pd_Bp = sim_param.Pd_Hp;%detection at Byzantine
%         sim_param.Pfa_Bp = sim_param.Pfa_Hp;%false alarm at byzantine
%
%         Pfa = 100 - Pd;
%
%         mkdir([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR)]);
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\parametersTEST.mat'], 'sim_param');
%
%
%         sim_param.Nsoglie_LLR =  6; % gam_var; % different LLR TH
%         N_pmal= 5;
%         Pmal=zeros(1,N_pmal);
%
%
%         i_pmal=1;
%         for p_mal_dec=6:1:6 + N_pmal - 1 %10
%
%             sim_param.Pmal = p_mal_dec/10;
%             Pmal(i_pmal) = p_mal_dec/10;
%
%             %fprintf('Threshold LLR Gamma = %f\n',gam_var);
%             fprintf('Pmal LLR = %f\n',sim_param.Pmal);
%
%             %--------------------------------------------------------------------------
%
%             [results, SOGLIA_LLRs] = Simula_identifica_Byz_mod_def(sim_param);
%             % PD_IDB = results.PD_IDB;
%             % PFA_IDB = results.PFA_IDB;
%             P_ISO_H= results.P_ISO_H;
%             P_ISO_B= results.P_ISO_B;
%             % PD_IDB_LLR = results.PD_IDB_LLR;
%             % PFA_IDB_LLR = results.PFA_IDB_LLR;
%             P_ISO_H_LLR = results.P_ISO_H_LLR;
%             P_ISO_B_LLR = results.P_ISO_B_LLR;
%             PFA = results.PFA;
%             PD = results.PD;
%             PFAr = results.PFAr; % Varshney
%             PDr = results.PDr;
%             PFAr_LLR = results.PFAr_LLR; % LLR
%             PDr_LLR = results.PDr_LLR;
%
%
%             % WITHOUT CHOOSING THE THRESHOLD WHISH GIVES THE MINIMUM
%
%             PERR = PFAr + 1-PDr;
%             PERR_LLR = PFAr_LLR+1-PDr_LLR;
%
%
%
%             save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\resultsTEST_Pmal=',num2str(Pmal(i_pmal)),'.mat'], 'results', 'PERR', 'SOGLIA_LLRs', 'PERR_LLR');
%
%             honest_iso_llr{1,i_pmal}= P_ISO_H_LLR; %honest_iso_llr(i_pmal,i_gam)= 1-PD_IDB_LLR(gam_var);
%             %
%             byzan_iso_llr{1,i_pmal}=P_ISO_B_LLR;
%
%             payoff_LLR{1,i_pmal}= PERR_LLR;
%             payoff{1,i_pmal}= PERR;
%             soglia{1,i_pmal} = SOGLIA_LLRs;
%
%             varsh_iso_ho{1,i_pmal}=P_ISO_H;
%             varsh_iso_by{1,i_pmal}=P_ISO_B;
%
%             llr_iso_ho{1,i_pmal}=P_ISO_H_LLR;
%             llr_iso_by{1,i_pmal}=P_ISO_B_LLR;
%
%             %fprintf('Perr minima per lo schema Varshney = %f\n',PERR);
%             fprintf('Perr per lo schema LLR = %f\n',PERR_LLR);
%
%             i_pmal=i_pmal+1;
%         end
%
%
%
%         %%%%% Varshney scheme %%%%%
%
%         GI= [];
%         GI_PISO = [];
%         GI_ISO_ratioBH = [];
%
%
%         for i=1:N_pmal
%             GI(:,i) = payoff{1,i}';   %M(:) mette tutti gli elementi in fila
% %             GI_HONEST_ISO(:,i)=varsh_iso_ho{1,i}';
% %             GI_BYZ_ISO(:,i)= varsh_iso_by{1,i}';
%             GI_PISO(:,i) = (varsh_iso_ho{1,i} - varsh_iso_by{1,i})';
% %               if (varsh_iso_ho{1,i} == 0 && varsh_iso_by{1,i} == 0) % NON MI ASPETTO DI DOVER GESTIRE IL CASO DI SOLO DENOMINATORE = 0 PERCHE? I BIZANTINI RIMOSSI NON SARANNO MAI MENO DEGLI ONESTI RiMOSSI
% %                GI_ISO_ratioBH(:,i) = 1;
% %              else
%                GI_ISO_ratioBH(:,i) = (varsh_iso_by{1,i}./varsh_iso_ho{1,i})';
%             %end
%             i=i+1;
%         end
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\Vgame.mat'],'GI');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_Pr_iso.mat'],'GI_PISO');
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_iso_ratio.mat'],'GI_ISO_ratioBH');
%
%
%
%         %%%%% LLR scheme %%%%%
%
%         GII= [];
% %         GII_HONEST_ISO=[];
% %         GII_BYZ_ISO=[];
%         GII_PISO = [];
%         GII_ISO_ratioBH = [];
%
%         for i=1:N_pmal
%             GII(:,i) = payoff_LLR{1,i}';   %M(:) mette tutti gli elementi in fila
% %             GII_HONEST_ISO(:,i)=llr_iso_ho{1,i}';
% %             GII_BYZ_ISO(:,i)=llr_iso_by{1,i}';
%             GII_PISO(:,i) = (llr_iso_ho{1,i} - llr_iso_by{1,i})';
% %             if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0)
% %               GII_ISO_ratioBH(:,i) = 1;
% %             else
%               GII_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
%             %end
%             i=i+1;
%         end
%
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLRgame.mat'],'GII');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLR_Pr_iso.mat'],'GII_PISO');
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLR_iso_ratio.mat'],'GII_ISO_ratioBH');
%
%
%
%
%
%     end
% end
% end
%
%
%
%
% %
% %
% clear all
%
% %%%%%%%%%%%%%%%%%%%%% altre prove in cui vario M e K %%%%%%%%%%%%%%%%%% %%%%%%%%%%%%
% % (da inserire in un nuovo for)
%
% sim_param.K = 53;
%
% folder1 = 'Nprove_50000\N=100\MoreOnest_K=53\'; % 'div_TH_AR\N=100\MoreOnest_K=54\';
%
% % li riscrivo dopo il clear all
% sim_param.PH1 = 0.5;  % leave as it is
% sim_param.N = 100;
% sim_param.T = 4;
% sim_param.gammas = 0:sim_param.T;
%
%
%
% sim_param.Nprove = 50000;
%
%
% for TH=  50 %30:10:70
%     sim_param.L = TH;
%
%     for TH_AR = 30:10:70
%     sim_param.L2 = TH_AR;
%
%     for Pd = 80:10:90
%
%         sim_param.Pd_Hp = Pd/100;%detection at honest
%         sim_param.Pfa_Hp = (100 - Pd)/100;%fa at honest
%         sim_param.Pd_Bp = sim_param.Pd_Hp;%detection at Byzantine
%         sim_param.Pfa_Bp = sim_param.Pfa_Hp;%false alarm at byzantine
%
%         Pfa = 100 - Pd;
%
%         mkdir([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR)]);
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\parametersTEST.mat'], 'sim_param');
%
%
%         sim_param.Nsoglie_LLR =  6; % gam_var; % different LLR TH
%         N_pmal= 5;
%         Pmal=zeros(1,N_pmal);
%
%
%         i_pmal=1;
%         for p_mal_dec=6:1:6 + N_pmal - 1 %10
%
%             sim_param.Pmal = p_mal_dec/10;
%             Pmal(i_pmal) = p_mal_dec/10;
%
%             %fprintf('Threshold LLR Gamma = %f\n',gam_var);
%             fprintf('Pmal LLR = %f\n',sim_param.Pmal);
%
%             %--------------------------------------------------------------------------
%
%             [results, SOGLIA_LLRs] = Simula_identifica_Byz_mod_def(sim_param);
%             % PD_IDB = results.PD_IDB;
%             % PFA_IDB = results.PFA_IDB;
%             P_ISO_H= results.P_ISO_H;
%             P_ISO_B= results.P_ISO_B;
%             % PD_IDB_LLR = results.PD_IDB_LLR;
%             % PFA_IDB_LLR = results.PFA_IDB_LLR;
%             P_ISO_H_LLR = results.P_ISO_H_LLR;
%             P_ISO_B_LLR = results.P_ISO_B_LLR;
%             PFA = results.PFA;
%             PD = results.PD;
%             PFAr = results.PFAr; % Varshney
%             PDr = results.PDr;
%             PFAr_LLR = results.PFAr_LLR; % LLR
%             PDr_LLR = results.PDr_LLR;
%
%
%             % WITHOUT CHOOSING THE THRESHOLD WHISH GIVES THE MINIMUM
%
%             PERR = PFAr + 1-PDr;
%             PERR_LLR = PFAr_LLR+1-PDr_LLR;
%
%
%
%             save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\resultsTEST_Pmal=',num2str(Pmal(i_pmal)),'.mat'], 'results', 'PERR', 'SOGLIA_LLRs', 'PERR_LLR');
%
%             honest_iso_llr{1,i_pmal}= P_ISO_H_LLR; %honest_iso_llr(i_pmal,i_gam)= 1-PD_IDB_LLR(gam_var);
%             %
%             byzan_iso_llr{1,i_pmal}=P_ISO_B_LLR;
%
%             payoff_LLR{1,i_pmal}= PERR_LLR;
%             payoff{1,i_pmal}= PERR;
%             soglia{1,i_pmal} = SOGLIA_LLRs;
%
%             varsh_iso_ho{1,i_pmal}=P_ISO_H;
%             varsh_iso_by{1,i_pmal}=P_ISO_B;
%
%             llr_iso_ho{1,i_pmal}=P_ISO_H_LLR;
%             llr_iso_by{1,i_pmal}=P_ISO_B_LLR;
%
%             %fprintf('Perr minima per lo schema Varshney = %f\n',PERR);
%             fprintf('Perr per lo schema LLR = %f\n',PERR_LLR);
%
%             i_pmal=i_pmal+1;
%         end
%
%
%
%         %%%%% Varshney scheme %%%%%
%
%         GI= [];
%         GI_PISO = [];
%         GI_ISO_ratioBH = [];
%
%
%         for i=1:N_pmal
%             GI(:,i) = payoff{1,i}';   %M(:) mette tutti gli elementi in fila
% %             GI_HONEST_ISO(:,i)=varsh_iso_ho{1,i}';
% %             GI_BYZ_ISO(:,i)= varsh_iso_by{1,i}';
%             GI_PISO(:,i) = (varsh_iso_ho{1,i} - varsh_iso_by{1,i})';
% %               if (varsh_iso_ho{1,i} == 0 && varsh_iso_by{1,i} == 0) % NON MI ASPETTO DI DOVER GESTIRE IL CASO DI SOLO DENOMINATORE = 0 PERCHE? I BIZANTINI RIMOSSI NON SARANNO MAI MENO DEGLI ONESTI RiMOSSI
% %                GI_ISO_ratioBH(:,i) = 1;
% %              else
%                GI_ISO_ratioBH(:,i) = (varsh_iso_by{1,i}./varsh_iso_ho{1,i})';
%             %end
%             i=i+1;
%         end
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\Vgame.mat'],'GI');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_Pr_iso.mat'],'GI_PISO');
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\V_iso_ratio.mat'],'GI_ISO_ratioBH');
%
%
%
%         %%%%% LLR scheme %%%%%
%
%         GII= [];
% %         GII_HONEST_ISO=[];
% %         GII_BYZ_ISO=[];
%         GII_PISO = [];
%         GII_ISO_ratioBH = [];
%
%         for i=1:N_pmal
%             GII(:,i) = payoff_LLR{1,i}';   %M(:) mette tutti gli elementi in fila
% %             GII_HONEST_ISO(:,i)=llr_iso_ho{1,i}';
% %             GII_BYZ_ISO(:,i)=llr_iso_by{1,i}';
%             GII_PISO(:,i) = (llr_iso_ho{1,i} - llr_iso_by{1,i})';
% %             if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0)
% %               GII_ISO_ratioBH(:,i) = 1;
% %             else
%               GII_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
%             %end
%             i=i+1;
%         end
%
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLRgame.mat'],'GII');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%         save([folder1, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLR_Pr_iso.mat'],'GII_PISO');
%         save([folder1,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\LLR_iso_ratio.mat'],'GII_ISO_ratioBH');
%
%
%
%
%
%     end
% end
% end
%
%
%
%
% clear all
%
% %...e per K=60.....
%
% sim_param.K = 60;
%
% folder2 = 'Nprove_50000\N=100\MoreOnest_K=60\'; %'div_TH_AR\N=100\MoreOnest_K=60\';
%
% % li riscrivo dopo il clear all
% sim_param.PH1 = 0.5;  % leave as it is
% sim_param.N = 100;
% sim_param.T = 4;
% sim_param.gammas = 0:sim_param.T;
%
% sim_param.Nprove = 50000;
%
%
% for TH= 50 %30:10:70
%    sim_param.L = TH;
%
%      for TH_AR = 30:10:70
%     sim_param.L2 = TH_AR;
%
%        for Pd = 80:10:90
%
%         sim_param.Pd_Hp = Pd/100;%detection at honest
%         sim_param.Pfa_Hp = (100 - Pd)/100;%fa at honest
%         sim_param.Pd_Bp = sim_param.Pd_Hp;%detection at Byzantine
%         sim_param.Pfa_Bp = sim_param.Pfa_Hp;%false alarm at byzantine
%
%         Pfa = 100 - Pd;
%
%         mkdir([folder2,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR)]);
%         save([folder2,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR), '\parametersTEST.mat'], 'sim_param');
%
%
%         sim_param.Nsoglie_LLR =  6; % gam_var; % different LLR TH
%         N_pmal= 5;
%         Pmal=zeros(1,N_pmal);
%
%
%         i_pmal=1;
%         for p_mal_dec=6:1:6 + N_pmal - 1 %10
%
%             sim_param.Pmal = p_mal_dec/10;
%             Pmal(i_pmal) = p_mal_dec/10;
%
%             %fprintf('Threshold LLR Gamma = %f\n',gam_var);
%             fprintf('Pmal LLR = %f\n',sim_param.Pmal);
%
%             %--------------------------------------------------------------------------
%
%             [results, SOGLIA_LLRs] = Simula_identifica_Byz_mod_def(sim_param);
%             % PD_IDB = results.PD_IDB;
%             % PFA_IDB = results.PFA_IDB;
%             P_ISO_H= results.P_ISO_H;
%             P_ISO_B= results.P_ISO_B;
%             % PD_IDB_LLR = results.PD_IDB_LLR;
%             % PFA_IDB_LLR = results.PFA_IDB_LLR;
%             P_ISO_H_LLR = results.P_ISO_H_LLR;
%             P_ISO_B_LLR = results.P_ISO_B_LLR;
%             PFA = results.PFA;
%             PD = results.PD;
%             PFAr = results.PFAr; % Varshney
%             PDr = results.PDr;
%             PFAr_LLR = results.PFAr_LLR; % LLR
%             PDr_LLR = results.PDr_LLR;
%
%
%             % WITHOUT CHOOSING THE THRESHOLD WHISH GIVES THE MINIMUM
%
%             PERR = PFAr + 1-PDr;
%             PERR_LLR = PFAr_LLR+1-PDr_LLR;
%
%
%
%
%             save([folder2, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR),'\resultsTEST_Pmal=',num2str(Pmal(i_pmal)),'.mat'], 'results', 'PERR', 'SOGLIA_LLRs', 'PERR_LLR');
%
%             honest_iso_llr{1,i_pmal}= P_ISO_H_LLR; %honest_iso_llr(i_pmal,i_gam)= 1-PD_IDB_LLR(gam_var);
%             %
%             byzan_iso_llr{1,i_pmal}=P_ISO_B_LLR;
%
%
%             %
%             % PFA_before{1,i_pmal}=PFA;
%             % PD_before{1,i_pmal}=PD;
%
%
%             payoff_LLR{1,i_pmal}= PERR_LLR;
%             payoff{1,i_pmal}= PERR;
%             soglia{1,i_pmal} = SOGLIA_LLRs;
%
%
%             varsh_iso_ho{1,i_pmal}=P_ISO_H;
%             varsh_iso_by{1,i_pmal}=P_ISO_B;
%
%             llr_iso_ho{1,i_pmal}=P_ISO_H_LLR;
%             llr_iso_by{1,i_pmal}=P_ISO_B_LLR;
%
%
%             %fprintf('Perr minima per lo schema Varshney = %f\n',PERR);
%             fprintf('Perr per lo schema LLR = %f\n',PERR_LLR);
%
%             i_pmal=i_pmal+1;
%         end
%
%
%         %%%%% Varshney scheme %%%%%
%
%         GI= [];
%
%         GI_PISO = [];
%         GI_ISO_rationBH = [];
%
%
%         for i=1:N_pmal
%             GI(:,i) = payoff{1,i}';   %M(:) mette tutti gli elementi in fila
%              GI_PISO(:,i) = (varsh_iso_ho{1,i} - varsh_iso_by{1,i})';
% %              if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0) % NON MI ASPETTO DI DOVER GESTIRE IL CASO DI SOLO DENOMINATORE = 0 PERCHE? I BIZANTINI RIMOSSI NON SARANNO MAI MENO DEGLI ONESTI RiMOSSI
% %                GI_ISO_ratioBH(:,i) = 1;
% %              else
%                GI_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
%             %end
%
%             i=i+1;
%         end
%
%          save([folder2,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa), '\prove_L2=',num2str(TH_AR),'\Vgame.mat'],'GI');
%
%           %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%          save([folder2,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR), '\V_Pr_iso.mat'],'GI_PISO');
%          save([folder2,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa), '\prove_L2=',num2str(TH_AR), '\V_iso_ratio.mat'],'GI_ISO_ratioBH');
%
%         %%%%% LLR scheme %%%%%
%
%         GII= [];
%         GII_PISO = [];
%         GII_ISO_ratioBH = [];
%
%         for i=1:N_pmal
%             GII(:,i) = payoff_LLR{1,i}';   %M(:) mette tutti gli elementi in fila
% %
%             GII_PISO(:,i) = (llr_iso_ho{1,i} - llr_iso_by{1,i})';
%
% %             if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0)
% %               GII_ISO_ratioBH(:,i) = 1;
% %             else
%               GII_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
%             %end
%             i=i+1;
%
%         end
%
%
%
%         save([folder2, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR), '\LLRgame.mat'],'GII');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Honest_Pr_iso.mat'],'GII_HONEST_ISO');
%
%          %save(['prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Byz_Pr_iso.mat'],'GII_BYZ_ISO');
%
%         save([folder2,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR), '\LLR_Pr_iso.mat'],'GII_PISO');
%         save([folder2,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\prove_L2=',num2str(TH_AR), '\LLR_iso_ratio.mat'],'GII_ISO_ratioBH');
%
%
%     end
%      end
% end
%
%
%
%
%
%
%
% %
% % clear all
% %
% % %...e per K=64.....
% %
% % sim_param.K = 64;
% %
% % folder3 = 'N=100\MoreOnest_K=64\';
% %
% % % li riscrivo dopo il clear all
% % sim_param.PH1 = 0.5;  % leave as it is
% % sim_param.N = 100;
% % sim_param.T = 4;
% % sim_param.gammas = 0:sim_param.T;
% %
% % sim_param.Nprove = 5000;
% %
% %
% % for TH= 30:10:70
% %    sim_param.L = TH;
% %
% %     %Pd=80;
% %     for Pd = 60:10:80
% %
% %         sim_param.Pd_Hp = Pd/100;
% %         sim_param.Pfa_Hp = 1 - Pd/100;
% %         sim_param.Pd_Bp = sim_param.Pd_Hp;
% %         sim_param.Pfa_Bp = sim_param.Pfa_Hp;
% %
% %
% %         save([folder3,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\parametersTEST.mat'], 'sim_param');
% %
% %
% %         sim_param.Nsoglie_LLR =  6; % gam_var; % different LLR TH
% %         N_pmal= 5;
% %         Pmal=zeros(1,N_pmal);
% %
% %
% %         i_pmal=1;
% %         for p_mal_dec=6:1:6 + N_pmal - 1 %10
% %
% %             sim_param.Pmal = p_mal_dec/10;
% %             Pmal(i_pmal) = p_mal_dec/10;
% %
% %             %fprintf('Threshold LLR Gamma = %f\n',gam_var);
% %             fprintf('Pmal LLR = %f\n',sim_param.Pmal);
% %
% %             %--------------------------------------------------------------------------
% %
% %             [results, SOGLIA_LLRs] = Simula_identifica_Byz_mod_def(sim_param);
% %             % PD_IDB = results.PD_IDB;
% %             % PFA_IDB = results.PFA_IDB;
% %             P_ISO_H= results.P_ISO_H;
% %             P_ISO_B= results.P_ISO_B;
% %             % PD_IDB_LLR = results.PD_IDB_LLR;
% %             % PFA_IDB_LLR = results.PFA_IDB_LLR;
% %             P_ISO_H_LLR = results.P_ISO_H_LLR;
% %             P_ISO_B_LLR = results.P_ISO_B_LLR;
% %             PFA = results.PFA;
% %             PD = results.PD;
% %             PFAr = results.PFAr; % Varshney
% %             PDr = results.PDr;
% %             PFAr_LLR = results.PFAr_LLR; % LLR
% %             PDr_LLR = results.PDr_LLR;
% %
% %
% %             % WITHOUT CHOOSING THE THRESHOLD WHISH GIVES THE MINIMUM
% %
% %             PERR = PFAr + 1-PDr;
% %             PERR_LLR = PFAr_LLR+1-PDr_LLR;
% %
% %
% %
% %
% %             save([folder3, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\resultsTEST_Pmal=',num2str(Pmal(i_pmal)),'.mat'], 'results', 'PERR', 'SOGLIA_LLRs', 'PERR_LLR');
% %
% %             honest_iso_llr{1,i_pmal}= P_ISO_H_LLR; %honest_iso_llr(i_pmal,i_gam)= 1-PD_IDB_LLR(gam_var);
% %             %
% %             byzan_iso_llr{1,i_pmal}=P_ISO_B_LLR;
% %
% %
% %             %
% %             % PFA_before{1,i_pmal}=PFA;
% %             % PD_before{1,i_pmal}=PD;
% %
% %
% %             payoff_LLR{1,i_pmal}= PERR_LLR;
% %             payoff{1,i_pmal}= PERR;
% %             soglia{1,i_pmal} = SOGLIA_LLRs;
% %
% %
% %             varsh_iso_ho{1,i_pmal}=P_ISO_H;
% %             varsh_iso_by{1,i_pmal}=P_ISO_B;
% %
% %             llr_iso_ho{1,i_pmal}=P_ISO_H_LLR;
% %             llr_iso_by{1,i_pmal}=P_ISO_B_LLR;
% %
% %
% %             %fprintf('Perr minima per lo schema Varshney = %f\n',PERR);
% %             fprintf('Perr per lo schema LLR = %f\n',PERR_LLR);
% %
% %             i_pmal=i_pmal+1;
% %         end
% %
% %
% %         %%%%% Varshney scheme %%%%%
% %
% %         GI= [];
% %
% %         GI_PISO = [];
% %         GI_ISO_rationBH = [];
% %
% %
% %         for i=1:N_pmal
% %             GI(:,i) = payoff{1,i}';   %M(:) mette tutti gli elementi in fila
% %              GI_PISO(:,i) = (varsh_iso_ho{1,i} - varsh_iso_by{1,i})';
% % %              if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0) % NON MI ASPETTO DI DOVER GESTIRE IL CASO DI SOLO DENOMINATORE = 0 PERCHE? I BIZANTINI RIMOSSI NON SARANNO MAI MENO DEGLI ONESTI RiMOSSI
% % %                GI_ISO_ratioBH(:,i) = 1;
% % %              else
% %                GI_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
% %             %end
% %
% %             i=i+1;
% %         end
% %
% %          save([folder3,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\Vgame.mat'],'GI');
% %
% %          save([folder3,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\V_Pr_iso.mat'],'GI_PISO');
% %          save([folder3,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\V_iso_ratio.mat'],'GI_ISO_ratioBH');
% %
% %         %%%%% LLR scheme %%%%%
% %
% %         GII= [];
% %         GII_PISO = [];
% %         GII_ISO_ratioBH = [];
% %
% %         for i=1:N_pmal
% %             GII(:,i) = payoff_LLR{1,i}';   %M(:) mette tutti gli elementi in fila
% % %
% %             GII_PISO(:,i) = (llr_iso_ho{1,i} - llr_iso_by{1,i})';
% %
% % %             if (llr_iso_ho{1,i} == 0 && llr_iso_by{1,i} == 0)
% % %               GII_ISO_ratioBH(:,i) = 1;
% % %             else
% %               GII_ISO_ratioBH(:,i) = (llr_iso_by{1,i}./llr_iso_ho{1,i})';
% %             %end
% %             i=i+1;
% %
% %         end
% %
% %
% %
% %         save([folder3, 'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLRgame.mat'],'GII');
% %         save([folder3,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_Pr_iso.mat'],'GII_PISO');
% %         save([folder3,'prove_L=',num2str(TH),'\',num2str(Pd),'_',num2str(Pfa),'\LLR_iso_ratio.mat'],'GII_ISO_ratioBH');
% %
% %
% %     end
% % end
% %
%
%
