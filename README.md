# VPP Course - Dynamic Programming Exercises

MATLAB implementations of Dynamic Programming problems covering stochastic optimal control, shortest path algorithms, and multi-objective optimization.

## Course Content

### Exercise 2: Resource Allocation Problem
**Files:** `ex2/ex2.m`, `ex2/numerical_ex2.m`, `ex2/numerical_ex2_drunk.m`, `ex2/analytic_ex2.m`

**Problem:** Optimal resource consumption over a finite horizon
- **Objective:** Maximize utility from consuming a limited resource over N=3 time steps
- **Cost function:** `-2*sqrt(u)` (utility increases with square root of consumption)
- **Dynamics:** `x_{k+1} = x_k - u_k` (resource depletion)
- **Discount factor:** α = 0.8

**Variations:**
- **Deterministic version** (`numerical_ex2.m`): Standard finite-horizon DP with backward recursion
- **Stochastic "drunk" version** (`numerical_ex2_drunk.m`): Adds uncertainty where consumed resource may be further reduced by a random factor with probability `p=0.7` the reduction is 0, otherwise `λ=0.3` times the consumption

**Key concepts:** Bellman equation, cost-to-go function, optimal policy computation

### Exercise 3: Machine Maintenance Problem (Finite Horizon)
**Files:** `ex3/numerical_ex3.m`

**Problem:** Optimal maintenance policy for a degrading machine
- **States:** 5 levels of machine condition (1=new to 5=broken)
- **Actions:** Do nothing (0), Repair (1), Buy new machine (2)
- **Horizon:** N = 10 steps

**State transitions:**
- **Do nothing:** Machine gradually degrades (probabilistic transitions)
- **Repair:** Improves machine state with high probability
- **Buy new:** Returns to perfect condition

**Cost structure:**
- Production value decreases as machine degrades: `[1, 0.9, 0.6, 0.4, 0]`
- Action costs: Do nothing (0), Repair (-0.4), Buy new (-1)

**Key concepts:** Markov Decision Process (MDP), transition probability matrices, state-dependent policies

### Exercise 4: Inventory Management with Delays
**Files:** `ex4/numerical_ex4.m`

**Problem:** Multi-period inventory control with order delays
- **Capacity:** M = 4 units
- **Horizon:** N = 15 days
- **State space:** Current inventory (x) and yesterday's order arriving today (s)

**Dynamics:**
- Orders arrive one day after being placed
- Random daily demand: {0, 1, 2, 3} units with probabilities `[0.1, 0.3, 0.5, 0.1]`
- `x_{k+1} = min(max(x_k - w_k, 0) + s_k, M)`

**Costs:**
- Inventory holding cost: `sqrt(x)`
- Ordering cost: Fixed `0.2` + quantity ordered
- Revenue from sales: `3 * min(x, w)` (negative cost)

**Key concepts:** State augmentation, delayed controls, inventory policies, stockout penalties

### Exercise 5: Viterbi Decoding (Convolutional Codes)
**Files:** `ex5/Numerical_ex5.m`

**Problem:** Maximum likelihood decoding of convolutional codes using the Viterbi algorithm
- **Constraint length:** K = 4
- **Code rate:** 1/3 (one input bit produces three output bits)
- **Generators:** Three polynomials defining the encoding process
  - g₀ = 1101, g₁ = 1010, g₂ = 0110
- **States:** 8 possible states (2³ from 3-bit shift register)

**Task:** Decode a 105-bit received sequence (35 information bits × 3 code bits)

**Algorithm:**
1. **Precompute trellis:** For each state and input, compute next state and output bits
2. **Forward DP (Viterbi):** Find minimum Hamming distance path through trellis
3. **Backtracking:** Recover transmitted bit sequence

**Key concepts:** Dynamic programming on graphs, shortest path in trellis, error correction, forward-backward algorithm

### Exercise 6: A* Pathfinding Algorithm
**Files:** `ex6/cv_06_s1.m`, `ex6/cv_06a.mat`, `ex6/cv_06b.mat`

**Problem:** Shortest path in a spatial graph using A* algorithm
- Graph nodes with coordinates and neighborhood information
- Heuristic-guided search (A* = Dijkstra + heuristic)

**Implementation:**
- Open set for nodes to be evaluated
- Closed set for already evaluated nodes
- Cost functions: g(n) = actual cost from start, h(n) = heuristic to goal
- f(n) = g(n) + h(n) guides node expansion

**Visualization:** Plots found path on map loaded from `.mat` files

**Key concepts:** Informed search, heuristic functions, optimal pathfinding, graph traversal

### Exercise 7: Multi-Objective Knapsack Problem
**Files:** `ex7/ex7.m`

**Problem:** Knapsack problem with dual objectives (cost and reliability)
- **Capacity:** W = 12 units
- **Items:** N = 10 items with weight, value, and success probability
- **Objectives:**
  1. Maximize total value (cost)
  2. Maximize probability of success

**Approach:** Pareto-optimal frontier computation
- Transform probability to log-space for additive computation
- Store non-dominated solutions at each DP stage
- Final result: Pareto front of cost vs. reliability trade-offs

**Terminal cost:** Penalty based on unused capacity

**Key concepts:** Multi-objective optimization, Pareto optimality, scalarization, dominated solutions pruning

### Exercise 8: Linear Quadratic Regulator (LQR)
**Files:** `ex8/ex8.m`

**Problem:** Optimal control of a 2D point mass with gravity
- **State:** `[x, vx, y, vy, g]` (position, velocity, gravity estimate)
- **Control:** `[ux, uy]` (accelerations in x and y directions)
- **Horizon:** N = 50 steps with dt = 0.1

**Dynamics:**
```matlab
x_{k+1} = A*x_k + B*u_k + w_k
```
where w ~ N(0, Σ) is Gaussian noise

**Cost:**
- Quadratic terminal cost: `x_N' * Q_N * x_N`
- Control effort: `u_k' * R * u_k` with λ = 0.005

**Solution:** Riccati recursion to compute optimal feedback gains
- Backward recursion: `K_k = A'*(K_{k+1} - K_{k+1}*B*(B'*K_{k+1}*B + R)^{-1}*B'*K_{k+1})*A + Q`
- Optimal control: `u_k = L_k * x_k` where `L_k = -(B'*K_{k+1}*B + R)^{-1}*B'*K_{k+1}*A`

**Key concepts:** LQR, Riccati equation, optimal feedback control, continuous state spaces

### Exercise 9: Infinite Horizon Value Iteration
**Files:** `ex9/numerical_ex9.m`

**Problem:** Machine maintenance with infinite horizon (stationary policy)
- Same setup as Exercise 3, but with infinite time horizon
- **Discount factor:** α = 0.9
- Goal: Find stationary optimal policy μ*(x)

**Algorithm:** Value Iteration
1. Initialize J(x) = 0 for all states
2. Iterate until convergence:
   ```
   J_{new}(x) = min_u { g(x,u) + α * Σ P(x'|x,u) * J_old(x') }
   ```
3. Extract policy: `μ*(x) = argmin_u { ... }`

**Convergence:** Stops when `||J_{new} - J_old|| < ε` with ε = 1e-6

**Key concepts:** Stationary policies, value iteration, contraction mapping, infinite horizon discounted cost

## Theoretical Background

This course covers fundamental concepts in **Dynamic Programming** and **Optimal Control**:

1. **Finite Horizon Problems** (Ex 2-5, 7-8)
   - Deterministic and stochastic dynamics
   - Backward recursion (Bellman equation)
   - Cost-to-go functions

2. **Infinite Horizon Problems** (Ex 9)
   - Discounted cost criterion
   - Value iteration and policy iteration
   - Stationary optimal policies

3. **Specialized Applications**
   - **MDPs:** Discrete state/action spaces (Ex 3, 9)
   - **LQR:** Continuous states with quadratic costs (Ex 8)
   - **Graph problems:** Shortest paths and decoding (Ex 5, 6)
   - **Multi-objective:** Pareto optimality (Ex 7)

## Implementation

**Technology:** MATLAB/Octave  
**Methods:** Grid-based discretization, backward recursion, value iteration

Each exercise is self-contained:
```matlab
cd ex/ex2
ex2  % or numerical_ex2, etc.
```

Output includes optimal cost-to-go values, optimal policies, and simulation results.
