# NAND minimization

The problem in question is to represent a disjunctive form, i.e. a logic expression of the form:\
y = (A_1 A_2 ... A_n) v (B_1 B_2 ... B_m) v ... \
using an equivalent logic expression, which uses only 2-port NANDs and the number of NANDs is minimized.

The types of operations normally used for this conversion are: 
1. De Morgan's law: (AvB)' = A'B' 
2. Negation using NAND: A' = (A nand A) 
3. Double negation: A = A'' 

Utilizing these 3 operations one can convert a disjunctive form into a form using NANDs, and the sequence of operations determines the number of NANDs in the end. We can illustrate an example of how this is done on a simple form: 

y = A v B'C (double negation) \
y = (A v B'C)'' (De Morgan) \
y = (A'(B'C)')' (definition of NAND) \
y = A' nand (B'C)' (negation using NAND) (definition of NAND) \
y = (A nand A) nand (B' nand C) (negation using NAND) \
y = (A nand A) nand ((B nand B) nand C) 

In general, the way we eliminate disjunctions or conjuctions is double negation: 

AB = (AB)'' = (A nand B)' = (A nand B) nand (A nand B) \
AvB = (AvB)'' = (A'B')' = A' nand B' 

In case of disjunction it is not always optimal to write A' = A nand A. In the case where A is a conjuction, A' gives us the NAND form immediatly, where as A needs a double negation. Therefore, we would use more than 4 times the operations if we didn't utilize the negation. \
In the case where A is a disjunction negating it as A' = A nand A is optimal. Either approach(negating first or negating last) will yield the same form: 

(AvB)' = ((AvB)')'' = (A'B')'' = (A' nand B')' \
(AvB)' = ((AvB)'')' = ((A'B')')' = (A' nand B')' 

In the case of conjuction, it's always optimal to write it as was done above, because there is no negation to take advantage of(the negation is over an expression containing NAND). 

# Algorithm #

Let's now describe the algorithm for generating the optimal solution. The input is given in a disjunctive form:



