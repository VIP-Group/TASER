% =========================================================================
% -- Triangular Approximate SEmidefinite Relaxation (TASER) 
% -- JED in large SIMO Simulator
% -------------------------------------------------------------------------
% -- (c) 2016 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@cornell.edu and oc66@cornell.edu
% -------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our 
% -- paper: 
% -- Oscar Castañeda, Tom Goldstein, and Christoph Studer,
% -- "Data Detection in Large Multi-Antenna Wireless Systems via
% -- Approximate Semidefinite Relaxation," 
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- vol. 63, no. 12, pp. 2334-2346, Dec. 2016.
% =========================================================================

function TASER_SIMO_JED_sim(varargin)

  % -- set up default/custom parameters
  
  if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
        
    % set default simulation parameters     
    par.runId = 0;         % simulation ID (used to reproduce results)
    par.MR = 16;           % receive antennas 
    par.MT = 1;            % transmit antennas (set not larger than MR!) 
    par.Time = 16;         % time slots
    par.mod = 'QPSK';      % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 1e4;      % number of Monte-Carlo trials (transmissions)
    par.simName = ...      % simulation name (used for saving results)
      ['NCERR_', num2str(par.MR), 'x', num2str(par.MT), '_', ...
       num2str(par.Time), 'TS_', par.mod, '_', ...
       num2str(par.trials), 'Trials'];    
    par.SNRdB_list = ...   % list of SNR [dB] values to be simulated
      -10:2:10; 
    par.detector = ...     % define detector(s) to be simulated
      {'SIMO', ...         % You can also use 'SDR' and 'ML-JED'. We recommend
       'SIMO-CHEST',...    % using them for small systems (the execution
       'TASER'};           % time for large systems is of several hours!) 
                           % To use the exact SDR detector, you need to install CVX.
                           % You can download CVX here: http://cvxr.com/cvx/download/
    par.CHEST = 'yes';     % enable ('yes') or disable ('no') channel estimation 
    par.tmax = 20;         % Number of TASER iterations
    par.alphaScale = 0.99; % Alpha scale for TASER's step size.
    
  else
      
    disp('use custom simulation settings and parameters...')    
    par = varargin{1};     % only argument is par structure
    
  end

  % -- initialization
  
  % use runId random seed (enables reproducibility)
  rng(par.runId); 

  % set up Gray-mapped constellation alphabet (according to IEEE 802.11)
  switch (par.mod)
    case 'BPSK',
      par.symbols = [ -1 1 ];
    case 'QPSK', 
      par.symbols = [ -1-1i,-1+1i, ...
                      +1-1i,+1+1i ];
    case '16QAM',
      par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
                      -1-3i,-1-1i,-1+3i,-1+1i, ...
                      +3-3i,+3-1i,+3+3i,+3+1i, ...
                      +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
      par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                      +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];                         
    case '8PSK',
      par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
                      exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
                      exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
                      exp(1i*2*pi/8*4), exp(1i*2*pi/8*5)];
  end

  % extract average symbol energy
  par.Es = mean(abs(par.symbols).^2); 
  
  % precompute bit labels
  par.Q = log2(length(par.symbols)); % number of bits per symbol
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

  % track simulation time
  time_elapsed = 0;
  
  % -- start simulation 
  
  % initialize result arrays (detector x SNR)
  % packet error rate:
  res.PER = zeros(length(par.detector),length(par.SNRdB_list));
  % symbol error rate:
  res.SER = zeros(length(par.detector),length(par.SNRdB_list));
  % bit error rate:
  res.BER = zeros(length(par.detector),length(par.SNRdB_list));

  % trials loop
  tic
  for t=1:par.trials
      
    % generate random bit stream (antenna x bit x time slots)
    bits = randi([0 1],par.MT,par.Q,par.Time);
    bits(:,:,1)=0; % nail down first signal 
  
    % generate transmit symbol
    for time=1:par.Time
      idx(:,time) = bi2de(bits(:,:,time),'left-msb')+1;
      S(:,time) = par.symbols(idx(:,time));
    end
  
    % generate iid Gaussian channel matrix & noise matrices
    % receive noise:
    N = sqrt(0.5)*(randn(par.MR,par.Time)+1i*randn(par.MR,par.Time));
    % channel:
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    % CHEST noise:
    NH = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    
    % transmit over noiseless channel (will be used later)
    X = H*S;
  
    % SNR loop
    for k=1:length(par.SNRdB_list)
      
      % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
      N0 = par.MT*par.Es*10^(-par.SNRdB_list(k)/10);
      
      % transmit data over noisy channel
      Y = X+sqrt(N0)*N;    
      Hest = H + sqrt(N0)*NH; % channel under estimation errors
            
      % algorithm loop      
      for d=1:length(par.detector)
          
        switch (par.detector{d})     % select algorithms
          case 'SIMO',               % SIMO LB with perfect CSI (no CHEST errors)              
            [idxhat,bithat] = SIMO(par,H,Y,S);   
          case 'SIMO-CHEST',         % SIMO LB with imperfect CSI
            [idxhat,bithat] = SIMO(par,Hest,Y,S);               
          case 'SDR',                % noncoherent SDR detection
            [idxhat,bithat] = SDR(par,Y);            
          case 'TASER',              % Triangular Approximate 
                                     % SEmidefinite Relaxation (TASER)
            [idxhat,bithat] = TASER(par,Y);            
          case 'ML-JED',             % ML detection using sphere decoding
            [idxhat,bithat] = ML(par,Y);            
          otherwise,
            error('par.detector type not defined.')      
        end

        % -- compute error metrics
        err = (idx~=idxhat);
        res.PER(d,k) = res.PER(d,k) + any(err(:));
        res.SER(d,k) = res.SER(d,k) + sum(err(:))/par.MT/par.Time;    
        res.BER(d,k) = res.BER(d,k) + ...
                         sum(bits(:)~=bithat(:))/(par.MT*par.Time*par.Q);      
      
      end % algorithm loop
                 
    end % SNR loop    
    
    % keep track of simulation time    
    if toc>10
      time=toc;
      time_elapsed = time_elapsed + time;
      fprintf('estimated remaining simulation time: %3.0f min.\n', ...
              time_elapsed*(par.trials/t-1)/60);
      tic
    end      
  
  end % trials loop

  % normalize results
  res.PER = res.PER/par.trials;
  res.SER = res.SER/par.trials;
  res.BER = res.BER/par.trials;
  res.time_elapsed = time_elapsed;
  
  % -- save final results (par and res structures)
    
  save([ par.simName '_' num2str(par.runId) ],'par','res');    
    
  % -- show results (generates fairly nice Matlab plot) 
  
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  figure(1)
  for d=1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.PER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.PER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('packet error rate (PER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-3 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)
  
end

% -- set of detector functions 

%% SIMO bound
function [idxhat,bithat] = SIMO(par,H,Y,S)
  Z = Y-H*S;
  for time=1:par.Time        
    for m=1:par.MT
      hm = H(:,m);
      yhat = Z(:,time)+hm*S(m,time);
      xhat(m,1) = hm'*yhat/norm(hm,2)^2;    
    end 
    [~,idxhat(:,time)] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
    bithat(:,:,time) = par.bits(idxhat(:,time),:);  
  end
end

%% detection via exact semidefinite relaxation (SDR)
%  You need to install CVX to use this.
function [idxhat,bithat] = SDR(par,Y)
  
  % -- Re-express the problem as in the Massive MU-MIMO case.
  y = Y(:,1)*conj(par.symbols(1));
  H = Y(:,2:end);        
  switch par.mod
    case 'BPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];
      % -- preprocessing for SDR  
      T = -[HR'*HR HR'*yR ; yR'*HR yR'*yR];        
      N = par.Time;        
    case 'QPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];
      % -- preprocessing for SDR  
      T = -[HR'*HR HR'*yR ; yR'*HR yR'*yR];        
      N = par.Time*2-1;     
    otherwise
      error('modulation type not supported')
  end
  
  % -- solve SDP via CVX
  cvx_begin quiet
    variable S(N,N) symmetric;
    S == semidefinite(N);       
    minimize( trace( T*S ) );
    diag(S) == 1;              
  cvx_end

  % -- post processing
  [V,U] = eig(S);
  root = V*sqrt(U);  
  sest = sign(root(:,end));  
  shat(1) = conj(par.symbols(1));
  
  switch par.mod
    case 'BPSK'
      shat(2:par.Time,1) = sest(1:par.Time-1); 
    case 'QPSK'        
      shat(2:par.Time,1) = sest(1:par.Time-1) + 1i*sest(par.Time:end-1);  
      shat = conj(shat);
    otherwise
      error('modulation type not supported')
  end
  
  % -- compute outputs  
  mt=1;
  [~,idxhat(mt,:)] = min(abs(shat*ones(1,length(par.symbols))-ones(par.Time,1)*par.symbols).^2,[],2);
  bithat(mt,:,:) = par.bits(idxhat(mt,:),:)';  
  
end

%% detection via Triangular Approximate SEmidefinite Relaxation (TASER)
function [idxhat,bithat] = TASER(par,Y)

  % -- Re-express the problem as in the Massive MU-MIMO case.
  y = Y(:,1)*conj(par.symbols(1));
  H = Y(:,2:end);  
      
  switch par.mod
    case 'BPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];
      % -- preprocessing for SDR  
      T = -[HR'*HR HR'*yR ; yR'*HR yR'*yR];           
    case 'QPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];
      % -- preprocessing for SDR 
      T = -[HR'*HR HR'*yR ; yR'*HR yR'*yR];        
    otherwise
      error('modulation type not supported')
  end
  
  % -- preconditioning for SDR
  DInv = diag(1./sqrt(abs(diag(T))));
  Ttilde = DInv*T*DInv;
  stepsize = par.alphaScale/norm(Ttilde,2);
  
  % -- use standard gradient on non-convex problem  
  gradf = @(L) 2*tril(L*Ttilde);
  proxg = @(L,t) prox_normalizer(L,diag(DInv).^-1);
  
  % Initialize Ltilde  
  Ltilde = diag(diag(DInv).^-1);  
  
  % -- Fast Iterative Soft Thresholding [Beck & Tebouille, 2009]  
  for k = 1:par.tmax
    Ltilde = proxg(Ltilde-stepsize*gradf(Ltilde)); % compute proxy    
  end  
  
  % -- post processing  
  sest = sign(Ltilde(end,:))';
  shat(1) = conj(par.symbols(1));
  
  switch par.mod
    case 'BPSK'   
      shat(2:par.Time,1) = sest(1:par.Time-1); 
    case 'QPSK'        
      shat(2:par.Time,1) = sest(1:par.Time-1) + 1i*sest(par.Time:end-1);  
      shat = conj(shat);
    otherwise
      error('modulation type not supported')
  end
  
  % -- compute outputs  
  mt=1;
  [~,idxhat(mt,:)] = min(abs(shat*ones(1,length(par.symbols))-ones(par.Time,1)*par.symbols).^2,[],2);
  bithat(mt,:,:) = par.bits(idxhat(mt,:),:)';  
  
end

% normalize columns of Z to have norm equal to its corresponding scale
function Q = prox_normalizer(Z,scale)
  [N,~] = size(Z);
  Q = Z.*(ones(N,1)*(1./sqrt(sum(abs(Z).^2,1)).*scale'));
end


%% ML detection using sphere decoding
function [idxhat,bithat] = ML(par,Y)

  % -- Re-express the problem as in the Massive MU-MIMO case.
  y = Y(:,1)*conj(par.symbols(1));
  H = Y(:,2:end);  
      
  switch par.mod
    case 'BPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];
      % -- preprocessing
      T = [HR'*HR HR'*yR ; yR'*HR yR'*yR];
      N = par.Time;
    case 'QPSK'
      % -- convert to real domain                  
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];
      % -- preprocessing
      T = [HR'*HR HR'*yR ; yR'*HR yR'*yR];
      N = par.Time*2-1;        
    otherwise
      error('modulation type not supported')
  end
  
  symbols = [ -1 1 ];
  
  % -- initialization  
  Radius = inf;
  PA = zeros(N,1); % path
  ST = zeros(N,length(symbols)); % stack  

  % -- preprocessing
  R = chol(max(eig(T)+1)*eye(N)-T);
  
  % -- add root node to stack
  Level = N; 
  ST(Level,:) = abs(R(Level,Level)*symbols.').^2;
  
  % -- begin sphere decoder
  while ( Level<=N )          
    % -- find smallest PED in boundary    
    [minPED,idx] = min( ST(Level,:) );
    
    % -- only proceed if list is not empty
    if minPED<inf
      ST(Level,idx) = inf; % mark child as tested        
      NewPath = [ idx ; PA(Level+1:end,1) ]; % new best path
      
      % -- search child
      if ( minPED<Radius )
        % -- valid candidate found
        if ( Level>1 )                  
          % -- expand this best node
          PA(Level:end,1) = NewPath;
          Level = Level-1; % downstep
          DF = R(Level,Level+1:end) * symbols(PA(Level+1:end,1)).';         
          ST(Level,:) = minPED + abs(R(Level,Level)*symbols.'+DF).^2;
        else
          % -- valid leaf found     
          idxML = NewPath;
          bitML = par.bits(idxML',:);
          % -- update radius (radius reduction)
          Radius = minPED;    
        end
      end      
    else
      % -- no more childs to be checked
      Level=Level+1;      
    end    
  end  

 
  % -- post processing
  sest = symbols(idxML).';
  shat(1) = conj(par.symbols(1));
  
  switch par.mod
    case 'BPSK'  
      shat(2:par.Time,1) = -sest(1:par.Time-1); 
    case 'QPSK'
      shat(1) = par.symbols(1);
      shat(2:par.Time,1) = -sest(1:par.Time-1) + 1i*sest(par.Time:end-1);  
    otherwise
      error('modulation type not supported')
  end
  
  % -- compute outputs  
  mt=1;
  [~,idxhat(mt,:)] = min(abs(shat*ones(1,length(par.symbols))-ones(par.Time,1)*par.symbols).^2,[],2);
  bithat(mt,:,:) = par.bits(idxhat(mt,:),:)';   
  
end