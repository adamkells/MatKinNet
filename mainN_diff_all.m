close all
clear all
addpath('Data') %folder containing all the data
addpath('support') %folder containing all the functions called in the main script
potential_type = 0;
%% CLUSTERING %%%%%%%%%%%%%%%
load('testsystem_nonoise_deeper_2.mat')
dd=0;
% important parameters to present
T=N*10;
n_sim=80;

one_vec=ones(1,length(K));
INV_K=(inv(eq*one_vec-K));

color_scheme=[0,0,1;1,0,0;0,1,0;1,1,0;0,1,1;1,0,1;0,0,0;0.5,0.5,0.5;0,0,0.5;0.5,0,0;0,0.5,0;0.5,0.5,0;0,0.5,0.5;0.5,0,0.5];

% committor and end points
[end_points]=ep_choice(K,K_eig_R);
[committor]=compute_commit(K',end_points(2,:)); % find committor for each state
[~,tmp2]=sort(committor); % order all nodes from the committor

for NCLUS=4 %number of clusters to look for
    counter=0;
    px=0;
    NCLUS
    red_method=0;
    
    for param=0
        close all
        clear cc
        dd = dd+1;
        counter=counter+1;
        
        % initial boundary guesses, this could be improved
        boundary(1,:)=randi([2,floor(N/(NCLUS-1))-1],[1,n_sim]);
        for i=1:n_sim
            for j=2:(NCLUS-1)
                boundary(j,i)=randi([boundary((j-1),i)+2,(floor(j*N/(NCLUS-1))-1)],1);
            end
        end
        
        % This code creates a clustering matrix A from the boundaries which
        % specifies which cluster each node belongs to
        A=zeros([N,NCLUS,n_sim]);
        for i=1:n_sim
            A(tmp2(1:boundary(1,i)),1,i)=1;
            for j=2:(NCLUS-1)
                A(tmp2(boundary((j-1),i)+1:boundary(j,i)),j,i)=1;
            end
            A(tmp2(boundary(end,i)+1:end),NCLUS,i)=1;
        end
        
        
        kemeny_latest=0;
        % calculate the kemeny value of each starting configuration
        for i=1:n_sim
            [kemeny_latest(i)]=kemeny_boundary(K,INV_K,eq,A(:,:,i),red_method,param);
        end
        
        % initializing the temperature to be very small
        temp=linspace(1,10,n_sim);
        temp=temp*0.001;
        
        % count of how many proposed moves are accepted, used to calculate the
        % neccesary change to the temperature
        switchcount=0;
        attswitch=0;
        switchcount2(1:n_sim)=0;
        switchcount3(1:n_sim)=0;
        
        % preset the modularity value
        Q=0;
        t=0;
        
        [a,b]=max(kemeny_latest);
        % if the optimum is better than the previous best, update
        best_split=A(:,:,b);
        best_kem_yet=a;
        
        % start the actual parallel tempering
        while t<T
            %counter to display progress
                            if mod(t,1000)==0
            
                                disp(['Completion is:', num2str(100*t/T), '%']);
                                disp(['Best value is:', num2str(best_kem_yet)]);
                                disp(['Modularity is:', num2str(Q)]);
                                %disp(['Acceptance is:', num2str(px)]);
                            end
            
            t=t+1;
            
            %% new configs for each parallel sim
            for i=1:n_sim
                
                [A_new]=make_new_config(i,A,Adj);
                [kemenyR_new]=kemeny_boundary(K,INV_K,eq,A_new(:,:,i),red_method,param);
                
                if kemenyR_new>kemeny_latest(i)
                    switchcount=switchcount+1;
                    switchcount2(i)=switchcount2(i)+1;
                    kemeny_latest(i)=kemenyR_new;
                    A(:,:,i)=A_new(:,:,i);
                elseif kemenyR_new<=kemeny_latest(i)
                    attswitch=attswitch+1;
                    val=rand(1);
                    condition = exp((kemenyR_new-kemeny_latest(i))*(1/temp(i)));
                    if val<condition
                        switchcount=switchcount+1;
                        switchcount3(i)=switchcount3(i)+1;
                        kemeny_latest(i)=kemenyR_new;
                        A(:,:,i)=A_new(:,:,i);
                    end
                end
                
            end
            
            %% tempering
            if mod(t,10)==0
                for j=1:(n_sim-1)
                    ii=j+1;
                    jj=j;
                    %keyboard
                    if kemeny_latest(ii)>kemeny_latest(jj)
                        kemeny_latest([ii jj])=kemeny_latest([jj ii]);
                        A(:,:,[ii jj])=A(:,:,[jj, ii]);
                    else
                        val=rand(1);
                        condition = exp((kemeny_latest(jj)-kemeny_latest(ii))*((1/temp(jj))-(1/temp(ii))));
                        if condition<val
                            kemeny_latest([ii jj])=kemeny_latest([jj ii]);
                            A(:,:,[ii jj])=A(:,:,[jj, ii]);
                        end
                    end
                end
            end
            
            %% temperature updating
            % next I update temperature during the simulation to keep 50% acceptance
            if mod(t,10)==0
                if sum(switchcount3)==0
                    switchcount3(1)=1;
                end
                px=sum(switchcount3)/(attswitch); % currect accepted fraction since last update
                
                attswitch=0;
                switchcount3(1:n_sim)=0; % reset the counter of accepted configurations
                correction=log(px)/log(0.5);
                temp=temp*correction;
            end
            
            %% state sweep
            % after each timestep, sweep all the states and look for the optimum
            [a,b]=max(kemeny_latest);
            % if the optimum is better than the previous best, update
            if a>best_kem_yet
                t
                best_split=A(:,:,b);
                best_kem_yet=a
                % compute x from A
                for iii=1:size(A,1)
                    x(iii)=find(best_split(iii,:));
                end
                
                % Computation of the modularity of the best clustering
                m = sum(sum(Adj));
                Q = 0;
                COMu = unique(x);
                for jj=1:length(COMu)
                    Cj = find(x==COMu(jj));
                    Ec = sum(sum(Adj(Cj,Cj)));
                    Et = sum(sum(Adj(Cj,:)));
                    if Et>0
                        Q = Q + Ec/m-(Et/m)^2;
                    end
                end
            end
        end
%         
%         if mod(t,500)==0
%             for i=50:-1:0
%                 A(:,:,end-i)=best_split;
%             end
%         end
        
        [~,~,~,taus]=kemeny_boundary(K,INV_K,eq,best_split,red_method,param);
        %taus
        %% plotting
        % for network type stuff
        if potential_type==0
            figure(counter+1)
            h = plot(G,'Layout','force');
            for i=1:NCLUS
                highlight(h,find(best_split(:,i)),'NodeColor',color_scheme(i,:))
            end
            title('Best splitting','FontSize', 18)
            saveas(gcf,[num2str(NCLUS) '_pol_' num2str(param) '_figure.eps'],'epsc')
            
            cc(:,:,counter)=best_split;
            % for 2D potential
        elseif potential_type==1
            cluster_V=zeros(N2);
            for jjj=1:NCLUS
                [a,b]=find(best_split(:,jjj));
                new_ind1=mod(a,N2);
                new_ind1(new_ind1==0)=N2;
                new_ind2=1+(a-new_ind1)/N2;
                for iii=1:length(a)
                    cluster_V(new_ind1(iii),new_ind2(iii))=jjj;
                end
            end
            cc(:,:,counter)=best_split;
            dd(:,:,counter)=cluster_V;
            % for 1D potential
        elseif potential_type==2
            cc(:,:,counter)=best_split;
            figure(counter+1)
            plot(v,'linewidth',2)
            
            hold on
            state_c = sum(cc(:,:,counter));
            for stemp = 1:(NCLUS-1)
                stem(sum(state_c(1:stemp)),v(sum(state_c(1:stemp))),'linewidth',2)
            end
            ylim([0,7])
            xlabel('X','fontsize',20)
            ylabel('V(X) [kcal/mol]','linewidth',20)
            ax = gca;
            ax.FontSize = 16;
            box off
            saveas(gcf,[num2str(NCLUS) num2str(counter) 'hs' num2str(param) num2str(potential_type) '_figure_nn.eps'],'epsc')
            %saveas(gcf,[num2str(counter) 'lin' num2str(param) num2str(potential_type) '_figure.fig'],'fig')
        end
        
        % check the peq of states for the new clustering
        [kemenyR_new,R,PEQ]=kemeny_boundary(K,INV_K,eq,best_split,red_method,param);
        disp(PEQ)
        T2 = expm(R*slow_rels(1))
        save([num2str(NCLUS) '_pol_' num2str(param) '_data.mat'])
    end
    
    P=expm(K');
        [chi]=PCCA_plus(P,NCLUS);
        for i=1:size(chi,1)
            [~,ind_tmp]=max(chi(i,:));
            chi(i,:)=0;
            chi(i,ind_tmp)=1;
        end
        if potential_type==0
            counter = counter + 1;
            figure(counter+1)
            h = plot(G,'Layout','force');
            for i=1:NCLUS
                highlight(h,find(chi(:,i)),'NodeColor',color_scheme(i,:))
            end
            title('Best splitting','FontSize', 18)
            saveas(gcf,['PCCA_' num2str(NCLUS) '_pol_' num2str(param) '_figure.png'],'png')
            %cc(:,:,counter)=best_split;
            % for 2D potential
        elseif potential_type==2
            counter = counter + 1;
            cc(:,:,counter)=chi;
    
            figure(counter+1)
            plot(v,'linewidth',2)
    
            hold on
            state_c = sum(cc(:,:,counter));
            for stemp = 1:(NCLUS-1)
                stem(sum(state_c(1:stemp)),v(sum(state_c(1:stemp))),'linewidth',2)
            end
            ylim([0,7])
            xlabel('X','fontsize',20)
            ylabel('V(X) [kcal/mol]','linewidth',20)
            ax = gca;
            ax.FontSize = 16;
            box off
            saveas(gcf,['PCCA_' num2str(NCLUS) 'lin' num2str(potential_type) '_figure_nn.eps'],'epsc')
            %saveas(gcf,[num2
       end
end

save(['results_' num2str(potential_type) '.mat'])
