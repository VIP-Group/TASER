% =========================================================================
% -- Triangular Approximate SEmidefinite Relaxation (TASER) 
% -- Massive MU-MIMO Simulator
% -------------------------------------------------------------------------
% -- (c) 2016 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@cornell.edu and oc66@cornell.edu
% -------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our 
% -- paper: Oscar Castañeda, Tom Goldstein, Christoph Studer,
% -- "Data Detection in Large Multi-Antenna Wireless Systems via
% -- Approximate Semidefinite Relaxation," IEEE Transactions on Circuits
% -- and Systems I, Nov. 2016
% =========================================================================

function TASER_MIMO_sim(varargin)

  % -- set up default/custom parameters
  
  if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
        
    % set default simulation parameters     
    par.runId = 0;         % simulation ID (used to reproduce results)
    par.MR = 32;           % receive antennas 
    par.MT = 32;           % transmit antennas (set not larger than MR!)         
    par.mod = 'QPSK';      % modulation type: 'BPSK','QPSK','16QAM','64QAM'             
    par.trials = 1e4;      % number of Monte-Carlo trials (transmissions)
    par.simName = ...      % simulation name (used for saving results)
      ['ERR_', num2str(par.MR), 'x', num2str(par.MT), '_', ...
        par.mod, '_', num2str(par.trials),'Trials'];
    par.SNRdB_list = ...   % list of SNR [dB] values to be simulated        
      6:2:16;              
    par.detector = ...     % define detector(s) to be simulated      
      {'SIMO','TASER', ... % You can also use 'SDR' and 'ML'. We recommend
       'MMSE'};            % using them for small systems (the execution
                           % time for large systems is of several hours!)
                           % To use the exact SDR detector, you need to install CVX.
                           % You can download CVX here: http://cvxr.com/cvx/download/
    par.tmax = 100;        % Number of TASER iterations
    par.alphaScale = 0.99; % Alpha scale for TASER's step size.
    %Step size used for different systems:
    %-------------------------------------------------------
    % MR / MT | Example system | Modulation | par.alphaScale
    % ratio   | (MRxMT)        | scheme     |
    %-------------------------------------------------------
    % 1       | 32x32          | BPSK       | 0.99
    % 2       | 64x32          | BPSK       | 0.95
    % 4       | 128x32         | BPSK       | 0.8
    % 8       | 256x32         | BPSK       | 0.75
    % 1       | 32x32          | QPSK       | 0.99
    % 2       | 64x32          | QPSK       | 0.99
    % 4       | 128x32         | QPSK       | 0.99
    % 8       | 256x32         | QPSK       | 0.85
    %-------------------------------------------------------
  else
      
    disp('use custom simulation settings and parameters...')    
    par = varargin{1};     % only argument is par structure
    
  end

  % -- initialization
  
  % use runId random seed (enables reproducibility)
  rng(par.runId,'twister'); 

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
  % vector error rate:
  res.VER = zeros(length(par.detector),length(par.SNRdB_list)); 
  % symbol error rate:
  res.SER = zeros(length(par.detector),length(par.SNRdB_list));
  % bit error rate:
  res.BER = zeros(length(par.detector),length(par.SNRdB_list));

  % generate random bit stream (antenna x bit x trial)
  bits = randi([0 1],par.MT,par.Q,par.trials);

  % trials loop
  tic
  for t=1:par.trials
  
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
  
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
  
    % SNR loop
    for k=1:length(par.SNRdB_list)
      
      % compute noise variance 
      % (average SNR per receive antenna is: SNR=MT*Es/N0)
      N0 = par.MT*par.Es*10^(-par.SNRdB_list(k)/10);
      
      % transmit data over noisy channel
      y = x+sqrt(N0)*n;
    
      % algorithm loop      
      for d=1:length(par.detector)

        switch (par.detector{d})     % select algorithms
          case 'SIMO',               % SIMO lower bound detector
            [idxhat,bithat] = SIMO(par,H,y,s);          
          case 'SDR',                % Detection via exact SDR
            [idxhat,bithat] = SDR(par,H,y);         
          case 'TASER',              % Triangular Approximate 
                                     % SEmidefinite Relaxation (TASER)
            [idxhat,bithat] = TASER(par,H,y); 
          case 'MMSE',              % unbiased MMSE detector
            [idxhat,bithat] = MMSE(par,H,y,N0);
          case 'ML',                 % ML detection using sphere decoding
            [idxhat,bithat] = ML(par,H,y);                                
          otherwise,
            error('par.detector type not defined.')      
        end

        % -- compute error metrics
        err = (idx~=idxhat);
        res.VER(d,k) = res.VER(d,k) + any(err);
        res.SER(d,k) = res.SER(d,k) + sum(err)/par.MT;    
        res.BER(d,k) = res.BER(d,k) + ...
                         sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);                   
        
      end % algorithm loop
                 
    end % SNR loop    
    
    % keep track of simulation time    
    if toc>10
      time = toc;
      time_elapsed = time_elapsed + time;
      fprintf('estimated remaining simulation time: %3.0f min.\n', ...
                time_elapsed*(par.trials/t-1)/60);
      tic
    end      
  
  end % trials loop
  
  % normalize results
  res.VER = res.VER/par.trials;
  res.SER = res.SER/par.trials;
  res.BER = res.BER/par.trials;
  res.time_elapsed = time_elapsed;
  
  % -- save final results (par and res structures)

  save([ par.simName '_' num2str(par.runId) ],'par','res');
  
  % -- show results (generates fairly nice Matlab plot) 
      
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  figure(1)
  for d = 1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.VER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.VER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('vector error rate (VER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-3 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)
    
end

% -- set of detector functions 

%% SIMO lower bound
function [idxhat,bithat] = SIMO(par,H,y,s)
  z = y-H*s;
  xhat = zeros(par.MT,1);
  for m=1:par.MT
    hm = H(:,m);
    yhat = z+hm*s(m,1);
    xhat(m,1) = hm'*yhat/norm(hm,2)^2;    
  end 
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);  
end

%% detection via exact semidefinite relaxation (SDR)
%  You need to install CVX to use this.
function [idxSDP,bitSDP] = SDR(par,H,y)

  switch par.mod
    case 'QPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];   
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = 2*par.MT+1; 
    case 'BPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];  
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = par.MT+1; 
    otherwise
      error('modulation type not supported')
  end  
  
  % -- solve SDP via CVX
  cvx_begin quiet
    variable A(N,N) symmetric;
    A == semidefinite(N);       
    minimize( trace( T*A ) );
    diag(A) == 1;              
  cvx_end
  
  % -- post processing
  [V,S] = eig(A);
  root = V*sqrt(S);
  
  sRhat = sign(root(:,end));  
  switch par.mod
    case 'QPSK'
      shat = sRhat(1:par.MT,1)+1i*sRhat(par.MT+1:end-1,1);
    case 'BPSK'  
      shat = sRhat(1:par.MT,1);
    otherwise
      error('modulation type not supported')
  end
  
  % -- compute outputs
  [~,idxSDP] = min(abs(shat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bitSDP = par.bits(idxSDP,:);
  
end

%% detection via Triangular Approximate SEmidefinite Relaxation (TASER)
function [idxSDP,bitSDP] = TASER(par,H,y)

  switch par.mod
    case 'QPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];   
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = 2*par.MT+1; 
    case 'BPSK'
      % -- convert to real domain
      yR = [ real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];  
      % -- preprocessing for SDR  
      T = [HR'*HR , -HR'*yR ; -yR'*HR yR'*yR ];
      N = par.MT+1; 
    otherwise
      error('not supported')
  end
  
  D = diag(diag(T).^-.5);
  That = D*T*D;
  stepsize = par.alphaScale/norm(That,2);

  % -- use standard gradient on non-convex problem  
  gradf = @(R) 2*triu(That*R);
  proxg = @(R,t) prox_normalizer(R,diag(D).^-1);
  
  x0 = eye(N);
  x0 = prox_normalizer(x0,diag(D).^-1);
      
  % -- Fast Iterative Soft Thresholding [Beck & Tebouille, 2009] 
  xk = x0;
  for k = 1:par.tmax
    xk = proxg(xk-stepsize*gradf(xk)); % compute proxy    
  end
  A = xk;
  
  % -- post processing
  rightvec = D*A(:,end);
  sRhat = sign(rightvec)/sign(rightvec(end));
  switch par.mod
    case 'QPSK'
      shat = sRhat(1:par.MT,1)+1i*sRhat(par.MT+1:end-1,1);
    case 'BPSK'  
      shat = sRhat(1:par.MT,1);
    otherwise
      error('not supported')
  end
  
  % -- compute outputs
  [~,idxSDP] = min(abs(shat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
  bitSDP = par.bits(idxSDP,:);
  
end

% normalize rows of Z to 1 to satisfy diag(Z*Z')=1
function Q = prox_normalizer(Z,scale)
 [N,~] = size(Z); 
  Q = Z.*((1./sqrt(sum(abs(Z).^2,2)).*scale)*ones(1,N));  
end

%% MMSE detector
function [idxhat,bithat] = MMSE(par,H,y,N0)
  W = (H'*H+(N0/par.Es)*eye(par.MT))\(H');
  xhat = W*y;
  G = real(diag(W*H));
  [~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
  bithat = par.bits(idxhat,:);
end

%% ML detection using sphere decoding
function [idxML,bitML] = ML(par,H,y)

  % -- initialization  
  Radius = inf;
  PA = zeros(par.MT,1); % path
  ST = zeros(par.MT,length(par.symbols)); % stack  

  % -- preprocessing
  [Q,R] = qr(H,0);  
  y_hat = Q'*y;    
  
  % -- add root node to stack
  Level = par.MT; 
  ST(Level,:) = abs(y_hat(Level)-R(Level,Level)*par.symbols.').^2;
  
  % -- begin sphere decoder
  while ( Level<=par.MT )          
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
          DF = R(Level,Level+1:end) * par.symbols(PA(Level+1:end,1)).';
          ST(Level,:) = minPED + abs(y_hat(Level)-R(Level,Level)*par.symbols.'-DF).^2;
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
  
end
