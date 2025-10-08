% Viterbi decoding (hard decision) via forward DP
clear; clc

%% ---------- Problem data (from Ex5) ----------
K  = 4;                 % constraint length
Gs = [1 1 0 1;          % g0 = 1101
      1 0 1 0;          % g1 = 1010
      0 1 1 0];         % g2 = 0110
Ngs = size(Gs,1);        % number of generators = 3  (rate 1/3)
S  = 2^(K-1);           % number of states = 8 (3-bit state)
s0 = 0;                 % initial state bits = 000 (as integer 0)

% Received bit string (105 bits = 35*3)
rx_str = '000010101010010101010100110101101001011010100010100110101101011110011001011100110101011010010010001010110';
rx_vec = rx_str - '0';
if mod(length(rx_vec), Ngs) ~= 0
    error('Length of received string must be divisible by number of generators');
end


N = length(rx_vec)/Ngs;  % number of trellis stages (input bits) = 35

rx = reshape(rx_vec, Ngs, N).';  % each row = [r0 r1 r2] 

%% precompute next state and output parity bits. 
% State is a (K-1)-bit register: [x(n-1) x(n-2) x(n-3)] as integer 0..7

next_state = zeros(S,2);     % next_state(s+1,u+1) in 1..S
out_bits   = zeros(S,2,Ngs);  % out_bits(s+1,u+1,gen)

for s = 0:S-1
    % expand s -> vector of K-1 bits [b1 b2 b3] (MSB first)
    sbits = bitget(s, K-1:-1:1);    % e.g., s=5 -> [1 0 1]
    for u = 0:1  %input
        reg = [u sbits];            % full reg length K
        % outputs: parity = sum(G .* reg) mod 2 for each generator row
        for gen_idx = 1:Ngs
            out_bits(s+1,u+1,gen_idx) = mod(sum(Gs(gen_idx,:).*reg), 2);
        end
        % next state = top K-1 bits after shift (drop last)
        ns_bits = reg(1:K-1);
        ns = ns_bits*[4;2;1];
        next_state(s+1,u+1) = ns + 1;
    end
end

%% ---------- Forward DP (Viterbi) ----------
INF = 1e9;
D  = INF*ones(S, N+1);      % D(s, k): best cost to ARRIVE in state s after k stages
BP = zeros(S, N+1, 'uint8');% BP(s, k): input bit that led to state s at stage k
PR = zeros(S, N+1, 'uint8');% PR(s, k): predecessor state index (1..S)

% init
D(:,1) = INF;
D(s0+1,1) = 0;

% iterate stages k=1..N
for k = 1:N
    r = rx(k,:);           % received Ng-bit vector at stage k
    for s = 1:S            % predecessor state index at stage k-1
        if D(s,k) >= INF/2, continue; end
        for u = 0:1
            ns = next_state(s, u+1);             % next state (1..S)
            tx = squeeze(out_bits(s, u+1, :)).'; % expected parity
            ham = sum(xor(tx, r));               % Hamming distance
            cand = D(s,k) + ham;
            if cand < D(ns, k+1)
                D(ns, k+1) = cand;
                BP(ns, k+1) = uint8(u);          % input used to arrive to ns at k
                PR(ns, k+1) = uint8(s);          % predecessor state
            end
        end
    end
end

%% ---------- Termination & backtrace ----------
% choose best terminal state (or set sN=1 to force 000)
[~, sN] = min(D(:, N+1));

u_hat = zeros(N,1,'uint8');
s = sN;
for k = N:-1:1
    u_hat(k) = BP(s, k+1);    % input that arrived to state s at stage k
    s        = PR(s, k+1);    % go to predecessor state at stage k
    if s==0, error('Unreachable state during backtrace'); end
end

decoded_bits = char(u_hat.' + '0');  % 35-bit string

%u_hat = flipud(u_hat);
pad = mod(5 - mod(N,5), 5);
u5  = reshape([u_hat; zeros(pad,1,'uint8')], 5, []).';

u5_str = join(string(u5), "")  % Each row as a string

vals = bin2dec(u5_str)       % indices 0..31

% fprintf('Decoded input bits (%d): %s\n', N, decoded_bits);
% fprintf('Grouped indices (0..31):\n'); disp(vals.')


symbols = ['A':'Z'];
symbols = symbols(:);  % make it a column vector

message = symbols(flipud(vals))  % map indices to letters

