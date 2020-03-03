# Automata Minimization

Given a Moore/Mealy automaton, the code finds an equivalent
compressed version. This is done by compressing similar states
into large superstates, which behave the same as the initial
states.
Note that the approach in this solution is neither optimal from
the correctness standpoint (it's a heuristic), nor is it very
efficient. The worst case runtime of the algorithm is 
O(n^2*m), where n is the number of states and m is the number of
outward transitions per state (number of outgoing edges).