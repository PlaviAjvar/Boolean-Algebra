# NAND minimization

The problem in question is to represent a disjunctive form, i.e. a logic expression of the form:

y = (A_1 A_2 ... A_n) v (B_1 B_2 ... B_m) v ... 

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

y = (A_1 A_2 ... A_n) v (B_1 B_2 ... B_m) v ... v (Z_1 Z_2 ... Z_p)

We will split the problem into two separate problems:
1. minimizing the number of NANDs in a conjuction
2. minimizing the number of NANDs in the disjunctive form

The reason we can do this is that negation always covers an entire expression, and can never cover only part of it. In other words, we can only negate a segment of the expression containing entire conjuctions separated by disjunctions.\
This means that whenever we apply double negation to a segment containing more than one conjuction, we have to transform the disjunction into NAND, since we don't have access to any of the conjuctions(they are in parentheses). We will only be able to negate conjuctions when there is only one conjuction in the segment we are double negating.

The first part of the problem is to solve the problem for conjuctions only. Afterwards, we know the optimal NAND expression for each of the conjuctions.\
In the disjunctive form, we can look at the conjuctions as indivisible blocks where each of them has a certain cost, represented by the number of NANDs used to represent it minimally. Taking these costs into consideration, we will try to optimize the sequence of negations over disjunctions to minimize the number of NANDs.

# First part #

Let's first describe how to minimize the number of NANDs in conjuctions. Let's say our conjuction is given by:

y = A_1 A_2 ... A_n

Let's define a function F_kon(i,j) which gives us the minimal number of NANDs needed to represent a segment [i,j] of the expression. In order to transform A[i...j] into a NAND form, we need to apply double negation over it, and split the expression around one of the conjuctions.

We can choose(because of the associativity of conjuction) which conjuction we will split by. Let's say we split by the conjuction between A_k and A_(k+1). We group the first k arguments and last n-k arguments, and then we double negate the k-th conjuction. We can now write a recursive relation for the optimal solution using this particular split:

F_kon(i,j)|k = 2\*F_kon(i,k) + 2\*F_kon(k+1,j) + 3

We are using the following property:

AB = (AB)'' = (A nand B)' = (A nand B) nand (A nand B)

F_kon(i,k) represents the optimal solution for A, F_kon(k+1,j) represents the optimal solution for B, and we need 3 more NANDs to calculate the above expression.\
Now it's easy to modify this in order to find the optimal solution. We simply try every possible k, split over it, and calculate the optimal solution using that split. Then we minimize over all these. Ergo:

F_kon(i,j) = min{F_kon(i,j)|k} = min{2\*F_kon(i,k) + 2\*F_kon(k+1,j) + 3}, for k = i,i+1,...,j-1

The base case is given by:

F_kon(i,i) = 0, if A[i] isn't negated\
F_kon(i,i) = 1, otherwise

The second case exists in the situation where a negated argument is given in the input, and A' = A nand A.\
All we need to do now is write the recursion and memoize it, and we have a dynamic programing solution for our problem. Since we are generally interested in both the number of NANDs in the optimal sequence, and the sequence itself, we will add another table into our program which will be used to backtrack the solution.

Let's analyze the complexity of the above algorithm. Since there are O(n^2) states, and each of them takes O(n) time to execute, the algorithm runs in O(n^3) time and O(n^2) space. 

# Second part #

Now that we've described how to minimie the conjuctions, let's describe the method of minimizing the disjunctive form. It can be written as:

y = A_1 v A_2 v ... A_m

where each of the A_i is a conjuction of arbitrarily many inputs. The approach for minimizing this is similar to the one before. First, note that:

AvB = (AvB)'' = (A'B')' = A' nand B'

Based on this define a function F_dis(i,j) which is the minimal number of NANDs which form an equivalent expression to the one on segment [i,j]. We will again apply double negation to the expression and let's again assume this is done around some index k. We can write:

F_dis(i,j)|k = 2\*F_dis(i,k) + 2\*F_dis(k+1,j) + 3, if k != i and k != j-1\
F_dis(i,j)|k = F_dis(i,k) + 2\*F_dis(k+1,j) + 2, if k = i\
F_dis(i,j)|k = 2\*F_dis(i,k) + F_dis(k+1,j) + 2, if k = j-1\
F_dis(i,j)|k = F_dis(i,k) + F_dis(k+1, j) + 1, j = i+1 = k+1

The absence of the factor 2 in the two special cases is due to the fact that we have a conjuction as one of the parameters. It is better in this case not to negate it twice and then again negate it within the expression, but rather use the fact that it is negated to our advantage, gaining a 4 times reduction in length.

Like in the previous part, the optimal solution is given by:

F_dis(i,j) = min{F_dis(i,j)|k}, for k = i,i+1,...,j-1

The base case reduces to the previous problem and is given by:

F_dis(i,i) = F_kon(1, len(A_i)) / 2

The division rounds down, which is intended as we need to remove the second operand and the NAND in the middle(since we are utilizing the inherent negation created by De Morgan to our advantage).\
Note that this doesn't hold only in one edge case, which is when there are zero disjunctions in the expression. In this case the one conjuction in the expression is not inherently negated, and we therefore need to handle this case separatelly. It reduces exactly to the problem in the previous section.

Again, we memoize the recursion and add a backtracking table for reconstructing the solution, and the problem is solved. The second part uses the first one only in the base case, and there is no additional computation since the lookup of the first table is O(1). Since the algorithm is nearly the same, the asymptotic complexity of the second part is the same as the first one O(n^3) time and O(n^2) space.
