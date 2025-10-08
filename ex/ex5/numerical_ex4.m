% Viterbi decoding (hard decision) via forward DP
% Ex5: K=4, g0=1101, g1=1010, g2=0110, initial state 000
clear; clc

%% ---------- Problem data (from Ex5) ----------
K  = 4;                 % constraint length
Gs = [1 1 0 1;          % g0 = 1101
      1 0 1 0;          % g1 = 1010
      0 1 1 0];         % g2 = 0110
Ng = size(Gs,1);        % number of generators = 3  (rate 1/3)
S  = 2^(K-1);           % number of states = 8 (3-bit state)
s0 = 0;                 % initial state bits = 000 (as integer 0)

% Received hard-decision bit string (105 bits = 35*3)
rx_str = '000010101010010101010100110101101001011010100010100110101101011110011001011100110101011010010010001010110';
if mod(length(rx_str), Ng) ~= 0
    error('Length of received string must be divisible by number of generators');
end
N = length(rx_str)/Ng;  % number of trellis stages (input bits) = 35

% Convert to numeric (N x Ng) array
rx = reshape(rx_str-'0', Ng, N).';  % each row = [r0 r1 r2] at time k

%% ---------- Precompute trellis: next state & output for (state,input) ----------
% State is a (K-1)-bit register: [x(n-1) x(n-2) x(n-3)] as integer 0..7
% On input u in {0,1}, new full register is [u, state_bits], outputs are dot over GF(2)
next_state = zeros(S,2);     % next_state(s+1,u+1) in 1..S
out_bits   = zeros(S,2,Ng);  % out_bits(s+1,u+1,gen)

for s = 0:S-1
    % expand s -> vector of K-1 bits [b1 b2 b3] (MSB first)
    sbits = bitget(s, K-1:-1:1);    % e.g., s=5 -> [1 0 1]
    for u = 0:1
        reg = [u sbits];            % full reg length K
        % outputs: parity = sum(G .* reg) mod 2 for each generator row
        for g = 1:Ng
            out_bits(s+1,u+1,g) = mod(sum(Gs(g,:).*reg), 2);
        end
        % next state = top K-1 bits after shift (drop last)
        ns_bits = reg(1:K-1);
        ns = 0;
        for j = 1:K-1
            ns = bitor(ns, bitshift(ns_bits(j), K-1-j));
        end
        next_state(s+1,u+1) = ns + 1; % 1-based indexing
    end
end

%% ---------- Forward DP (Viterbi) ----------
INF = 1e9;
D = INF*ones(N+1, S);        % D(k,s) = best path metric (cost-to-arrive) at time k in state s
BP = zeros(N, S, 'uint8');   % backpointer: which predecessor input (0/1) led to best cost

D(1, :) = INF;               % before any input
D(1, s0+1) = 0;              % cost at initial state is 0

% Iterate k = 1..N (each consumes one input bit, produces Ng parity bits)
for k = 1:N
    r = squeeze(rx(k, :));   % received Ng-bit vector at time k
    for s = 1:S
        if D(k, s) >= INF/2, continue; end  % unreachable
        for u = 0:1
            ns  = next_state(s, u+1);
            tx  = squeeze(out_bits(s, u+1, :)).';
            ham = sum( xor(tx, r) );       % Hamming distance for this branch
            cand = D(k, s) + ham;
            if cand < D(k+1, ns)
                D(k+1, ns) = cand;
                BP(k, ns)  = uint8(u);     % remember input bit that led here
            end
        end
    end
end

%% ---------- Termination & backtrace ----------
% If final state is not specified, choose the best one
[~, sN] = min(D(N+1, :));

u_hat = zeros(N,1,'uint8');
s = sN;
for k = N:-1:1
    u_hat(k) = BP(k, s);
    % find predecessor state: invert next_state mapping
    % brute-force (S=8 -> tiny): find s_prev,u such that next_state(s_prev,u)=s
    found = false;
    for sp = 1:S
        for u = 0:1
            if next_state(sp, u+1) == s && BP(k, s) == u
                s = sp;
                found = true;
                break
            end
        end
        if found, break; end
    end
end

decoded_bits = char(u_hat.' + '0');  % 35-bit string

%% ---------- Optional: group by 5 bits for your letter table ----------
L = length(u_hat);
pad = mod(5 - mod(L,5), 5);
u_pad = [u_hat; zeros(pad,1,'uint8')];
u5 = reshape(u_pad, 5, []).';         % rows of 5 bits
vals = bin2dec(u5, 'left-msb');         % 0..31 (use your table from Ex5.pdf)

% Print results
fprintf('Decoded input bits (%d): %s\n', N, decoded_bits);
fprintf('Grouped by 5 bits (index 0..31 per row):\n');
disp(vals.')