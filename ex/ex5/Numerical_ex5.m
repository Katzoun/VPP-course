% Viterbi decoding via forward DP
clear; clc

K  = 4;                 % constraint length
Gs = [1 1 0 1;          % g0 = 1101
      1 0 1 0;          % g1 = 1010
      0 1 1 0];         % g2 = 0110
Ngs = size(Gs,1);        % number of generators = 3  (rate 1/3)
S  = 2^(K-1);           % number of states = 8 (3-bit state)
x0 = 0;                 % initial state bits = 000 (as integer 0)

% Received bit string (105 bits = 35*3)
rx_str = '000010101010010101010100110101101001011010100010100110101101011110011001011100110101011010010010001010110';
rx_vec = rx_str - '0';
if mod(length(rx_vec), Ngs) ~= 0
    error('Length of received string must be divisible by number of generators');
end


N = length(rx_vec)/Ngs;  % number of trellis diagram stages (input bits) = 35

rx = reshape(rx_vec, Ngs, N)' ;


%% precompute next state and output parity bits. 
% State is a (K-1)-bit register: [x(n-1) x(n-2) x(n-3)] as integer 0..7

next_state = zeros(S,2);     % next_state(s+1,u+1) in 1..S
out_bits   = zeros(S,2,Ngs);  % out_bits(s+1,u+1,gen)

for x = 0:S-1
    % expand s -> vector of K-1 bits [b1 b2 b3] (Most significant bit first)
    sbits = bitget(x, K-1:-1:1);    % for example s=5 -> [1 0 1]
    for u = 0:1  %input
        reg = [u sbits];            % full reg length K
        % outputs: parity = sum(G .* reg) mod 2 for each generator row
        for gen_idx = 1:Ngs
            out_bits(x+1,u+1,gen_idx) = mod(sum(Gs(gen_idx,:).*reg), 2);
        end
        % next state = top K-1 bits after shift (drop last)
        ns_bits = reg(1:K-1);
        ns = ns_bits*[4;2;1];
        next_state(x+1,u+1) = ns;
    end
end

%% ---------- Forward DP (Viterbi) ----------
INF = 1e9;
D  = INF*ones(S, N+1);      % D(s, k): best cost to ARRIVE in state s after k stages
BP = zeros(S, N+1);          % BP(s, k): input bit that led to state s at stage k
PR = zeros(S, N+1);          % PR(s, k): predecessor state index (1..S)

% init
D(:,1) = INF;
D(x0+1,1) = 0;

% iterate stages k=1..N
for k = 1:N
    r = rx(k,:);           % received Ng-bit vector at stage k
    for x = 1:S            % predecessor state index at stage k-1
        if D(x,k) >= INF/2, continue; end
        for u = 0:1
            ns = next_state(x, u+1);             % next state (1..S)
            expected = squeeze(out_bits(x, u+1, :))'; % expected parity
            hamming = sum(xor(expected, r));               % Hamming distance
            cand = D(x,k) + hamming;
            if cand < D(ns+1, k+1)   %ns+1 because of indexing. ns = 0 is state [0 0 0], ns = 5 [0 1 1] ...
                D(ns+1, k+1) = cand;
                BP(ns+1, k+1) = u;         % input used to arrive to ns at k
                PR(ns+1, k+1) = x;        % predecessor state
            end
        end
    end
end

%% ---------- Termination & backtrace ----------
% choose best terminal state (or set sN=1 to force 000)
[minval, minpos] = min(D(:, N+1));

u_hat = zeros(N,1);
x = minpos;
for k = N:-1:1
    u_hat(k) = BP(x, k+1);    % input that arrived to state s at stage k
    x = PR(x, k+1);    % go to predecessor state at stage k
end

decoded_bits = char(u_hat.' + '0');  % 35-bit string

%% symbol mapping 
pad = mod(5 - mod(N,5), 5);
u5  = reshape([u_hat; zeros(pad,1)], 5, []).';
u5_str = join(string(u5), "");  % Each row as a string
vals = bin2dec(u5_str);      % indices 0..31

symbols = 'A':'Z';
symbols = symbols(:);  % makes it a column vector

message = symbols(flipud(vals))'  % map indices to letters

