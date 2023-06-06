function delta  = MDEA(Data, Stripesize, Rule, ST, EN, PLOT)

% Korosh Mahmoodi, 02/28/2023.

% Data: is the time series of interest. The minimum lentgh of data for MDEA
% is (2*10^4, 1).

% Stripesize: is the size of the strips for the Data normalized to [0, 1].
% suggested Stripesize is 0.01 ;

% Rules: can be 1,-1, or 0. It represents the way that diffusion trajectory
% (Diff) will be built using the extracted events:
% - Rule = 1: assigns +1 to the events and accumulates them. 
% - Rule = -1: assigns randomly +1 or -1 to the events and accumulates
% them.
% - Rule = 0: fills the region between the events either by +1 or -1
% (randomly) and accumulates them.
% The relationship between the scaling evaluated with each of these rules
% and the temporal complexity index mu is described in details in:
% [1] Grigolini, P., Palatella, L., & Raffelli, G. (2001). Asymmetric anomalous diffusion: An efficient way to detect memory in time series. Fractals,
% 09(04),439–449. https://doi.org/10.1142/s0218348x01000865
% [2] Scafetta, N., & Grigolini, P. (2002). Scaling detection in time series: Diffusion
% entropy analysis. Phys. Rev. E, 66, 036130. https://doi.org/10.1103/PhysRevE.66.036130

% For Rule = 1 (the most common one) mu=1+delta for mu<2 and delta=1/(mu-1)for
% mu >=2. 


% PLOT = 1 creates the Shannon Entropy of the diffusion trajectory vs. log of the window size: 
% The slope of the linear part of the graph is the scaling delta.
% ST is the beggining of the linear estimate; suggested ST = 0.4. 
% EN is the ending of the linear estimate;  suggested EN = 0.8.
% Note that ST and EN should measure the linear part on the S(l) vs. log(l)
% graph (asymptotic regime.)
                           
% Using this code please cite:
%[3] Mahmoodi, K., Kerick, S.E., Grigolini, P., Franaszczuk, P.J., & West, B.J. Complexity Synchronization, 
% https://arxiv.org/ftp/arxiv/papers/2210/2210.16675.pdf


%%% Normalizing the Data to [0 1]
                                   Data = Data - min(Data) ;
                                   Data = Data ./ max(Data) ;
                                   
                                   Lenghtdata = length(Data) ;
                                   
%%% Extracting events using stripes; Each time that the Data passes from
%%% one stripe to a different one gets recorded as an event.
Ddata = Data./(Stripesize); %% This projects Data to the interval [0 1/Stripesize]; consequently, whole numbers would limit stripes. So, using floor and ceil commands, we can determine when the time series passes from one stripe to a different one (line 57).
Event = zeros(Lenghtdata, 1);

k = 1 ;
Event(1) = 1 ;
StartEvent = zeros() ;

for i = 2 : Lenghtdata
    if ( Ddata(i) < floor(Ddata(i-1)) ) || (Ddata(i) > ceil(Ddata(i-1))  )
                        Event(i) = 1;
                        StartEvent(k) = i ;
k = k + 1;

    else

    end
end

%%% Creating diffusion trajectory (Diff) from extracted events (Event.)                        

                        StartEvent = StartEvent(StartEvent~=0) ;
                        if Rule == 1
                             Diff = cumsum(Event) ;   
                        end
                        
                        if Rule == -1
                            State0 = zeros( Lenghtdata, 1) ;
                                 for yy =  1 : Lenghtdata
                                     
                                     if Event(yy) == 1
                                         r = rand ; 
                                         if r <0.5
                                             State0(yy) = 1 ;
                                         else
                                         State0(yy) = -1 ;

                                         end
                                         
                                         
                                     end
                                     
                                     
                                 end
                                     Diff = cumsum(State0) ;   
                            
                        end
                        
                        if Rule == 0 
                            State0 = zeros( Lenghtdata, 1) ;
                                 for yy =  1 : Lenghtdata
                                     
                                     if Event(yy) == 1
                                         r = rand ; 
                                         if r <0.5
                                             State0(yy) = 1 ;
                                         else
                                         State0(yy) = -1 ;

                                         end
                                         
                                     end
                                     
                                 end
                                 State00 = zeros(Lenghtdata, 1) ;
                                 State0(1) = 1 ;
                                 for ee = 2 : Lenghtdata
                                     if State0(ee) ==  0
                                                State00(ee) = State00(ee - 1) ;
                                     else
                                               State00(ee) = State0(ee ) ;

                                     end
                                     
                                 end
                                          Diff = cumsum(State00) ;   

                        end
                      
%%% Evaluating the Shannon Entropy of the diffusion trajectory Diff.
% For each given window size ll, we divide Diff to the appropriate slices
% and treat each slice as an independent trajectory.  

ll =   floor(log(length(Diff))/log(1.2))   - 5 ;
Delh = zeros(1, ll);

de = zeros(1, ll);
DE = zeros(1, ll);


% makes an array of l of powers of 2
for i = 1 : ll
    Delh(i)=  floor(1.2^i);
end

    for q = 1:length(Delh)

        SliceNum = 1* length(StartEvent) ; %%% TimeStep / 100   ;
        
       del =  Delh(q) ;
       
  HH = zeros(SliceNum , 1) ;
  enn = length(Diff) ;
  enn2 = length(StartEvent) ;
  i = 1 ;
while StartEvent(i) + del <  enn 
    
        iiii = StartEvent(i) ;
                 HH(i) = Diff(iiii + del) - Diff(iiii) ;
        
   i = i + 1 ;
   if i == enn2
       break ;
   end
   
end

XF = HH(HH~=0);  
nbins = floor( max( abs(XF) ) /1 ) ; 
nbins(nbins==0) = [];
[counts] = hist(XF, nbins) ; % hist(XF,100); % hist(XF,nbins);

% discount any entry that is zero
counts = counts(counts ~= 0);
Pc = counts./sum(counts);

% this is the integral for Entropy (S)
DE0 = -sum((Pc).*log(Pc)/log(10));
        
        DE(:,q) = DE0;
    end


for t = 1 : length(Delh)
    de(t) = log(Delh(t))/log(10);
end

 Starr = round(ST *length(de)) ; 
 endd =  round(EN *length(de)); 
 DE0 = DE(   Starr    :  endd  );
de0 = de(    Starr   : endd  );

% performs a polyfit of degree 1 to find the slope, delta
FitLine = polyfit(de0,DE0, 1);
delta = FitLine(1) ;    % this parameter is the scaling

% % % plot
if PLOT == 1  
figure
subplot(1,2,1)
  plot(Ddata)   
xlabel('t'), ylabel('X(t)');
legend('X(t)','Location','northwest');
title('Signal');
subplot(1,2,2)
plot(de(3:length(de)-3),DE(3:length(DE)-3),'+',de0,FitLine(1)*de0+FitLine(2),'r--','LineWidth',1.5);
xlabel('log(l)'), ylabel('S(l)');
%legend(['N =' num2str(N) ', ' 'kT =' num2str(sprintf('%.3f',kT)),', ' '\delta = ' num2str( sprintf('%.3f',delta ))],'Location','northeast');

legend(['\delta = ' num2str( sprintf('%.3f',delta ))],'Location','northwest');

end  

end
