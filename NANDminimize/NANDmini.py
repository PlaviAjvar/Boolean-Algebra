# The task is to minimize the number of NAND-s in the expression obtained
# by transforming the initial one, utilizing De-Morgan's laws and properties of NAND.

inf = 10 ** 18  # artificial infinity, larger than any feasible length of expression


# utility function which optimizes nand usage for both cases(conjuctive and disjunctive)
# just need to define different base cases
def minimize_util(low_idx, high_idx, min_nand, backtrack):
    # base case
    if low_idx == high_idx:
        return min_nand[low_idx][high_idx]
    # if we've already calculated the value
    if min_nand[low_idx][high_idx] != inf:
        return min_nand[low_idx][high_idx]

    for split_point in range(low_idx, high_idx):
        # calculate minimal cost for split at split point
        split_cost = 2 * minimize_util(low_idx, split_point, min_nand, backtrack) +\
                     2 * minimize_util(split_point + 1, high_idx, min_nand, backtrack) + 3
        # if found better solution, update
        if split_cost < min_nand[low_idx][high_idx]:
            min_nand[low_idx][high_idx] = split_cost
            backtrack[low_idx][high_idx] = split_point

    return min_nand[low_idx][high_idx]


# generates optimal expression of conjuction from calculated backtracking table
def expression_conj(backtrack, conj):
    # utility function for recursivelly generating expression
    def expression_util(low_idx, high_idx):
        if low_idx == high_idx:
            if conj[low_idx] > 0:
                single_var = chr(conj[low_idx])
            else:
                single_var = "(" + chr(-conj[low_idx]) + "\u22BC" + chr(-conj[low_idx]) + ")"
            return single_var
        # AB = (AB)'' = (A nand B)' = (A nand B) nand (A nand B)
        left_expr = expression_util(low_idx, backtrack[low_idx][high_idx])
        right_expr = expression_util(backtrack[low_idx][high_idx] + 1, high_idx)
        if len(left_expr) > 1:  # unless the expression is a single letter, put parentheses around it
            left_expr = "(" + left_expr + ")"
        if len(right_expr) > 1:
            right_expr = "(" + right_expr + ")"

        half_expr = "(" + left_expr + "\u22BC" + right_expr + ")"  # (A nand B)
        full_expr = half_expr + "\u22BC" + half_expr
        return full_expr

    return expression_util(0, len(conj) - 1)


# finds minimal equivalent expression of conjuction using only NANDs
# argument is list of numbers corresponding to the indices of the input variables
def conjunction_nand(conj):
    num_var = len(conj)
    # minimal number of NANDs to represent segment in conjuction
    min_nand = [[inf] * num_var for i in range(num_var)]
    backtrack = [[-1] * num_var for i in range(num_var)]  # store optimal splitting point of each segment
    # base cases
    for i in range(num_var):
        if conj[i] > 0:
            min_nand[i][i] = 0   # zero NANDs are needed to write expression with a single non-negated variable
        else:
            min_nand[i][i] = 1   # one NAND needed when there is a single negation
        backtrack[i][i] = i

    minimize_util(0, num_var-1, min_nand, backtrack)
    return expression_conj(backtrack, conj)


# generates optimal expression of disjunctive form from calculated backtracking table
def expression_disj(backtrack, conj_nand):
    # utility function for recursivelly generating expression
    def expression_util(low_idx, high_idx):
        if low_idx == high_idx:
            return conj_nand[low_idx]
        # AvB = (AvB)'' = (A'B')' = (A nand A) nand (B nand B)
        left_expr = expression_util(low_idx, backtrack[low_idx][high_idx])
        right_expr = expression_util(backtrack[low_idx][high_idx] + 1, high_idx)
        if len(left_expr) > 1:  # unless the expression is a single letter, put parentheses around it
            left_expr = "(" + left_expr + ")"
        if len(right_expr) > 1:
            right_expr = "(" + right_expr + ")"

        first_half = "(" + left_expr + "\u22BC" + left_expr + ")"  # (A nand A)
        second_half = "(" + right_expr + "\u22BC" + right_expr + ")"  # (B nand B)
        full_expr = first_half + "\u22BC" + second_half
        return full_expr

    return expression_util(0, len(conj_nand) - 1)


# finds equivalent of disjunctive form with least number of NAND's
def disjunction_nand(conj_cost, conj_nand):
    num_var = len(conj_cost)
    # minimal number of NANDs to represent segment in disjunctive form
    min_nand = [[inf] * num_var for i in range(num_var)]
    backtrack = [[-1] * num_var for i in range(num_var)]  # store optimal splitting point of each segment
    # base cases
    for i in range(num_var):
        min_nand[i][i] = conj_cost[i]  # minimal cost of conjuction in disjunctive form
        backtrack[i][i] = i

    minimize_util(0, num_var - 1, min_nand, backtrack)
    return expression_disj(backtrack, conj_nand)


# the first argument is the list of lists(list of conjuctions)
# the inner lists contain the indexes of the logical input variables
# it is already assumed there are disjunctions between the conjuctive forms
# the second argument represents the number of variables in the initial logic function
def min_nand_form(dnf):
    # calculate min cost for conjuctions individually
    conj_nand = [conjunction_nand(conj) for conj in dnf]
    conj_cost = [len(nand_ex) for nand_ex in conj_nand]
    # find minimum equivalent disjunctive form taking into consideration costs of conjuctions
    nand_form = disjunction_nand(conj_cost, conj_nand)
    return nand_form


# convert DNF in string form to DNF in list form(suitable for the algorithm above)
# input should be in form "ADEvB'CvG'H" with letter v representing disjunction
# any character not 'v' or space will be considered a logical variable
# negations are specified with '
def parse_dnf(dnf_str):
    # we first strip all whitespace
    dnf_str = "".join(dnf_str)
    conj_list = dnf_str.split('v')  # splits expression by disjunction
    conj_num = len(conj_list)

    dnf = [[] for i in range(conj_num)]
    for conj_idx in range(conj_num):
        conj_len = len(conj_list[conj_idx])
        for i in range(conj_len):
            if conj_list[conj_idx][i].isalnum():
                # pass the ascii value of the character into the dnf list
                # if the argument is negated, pass it's negative ascii value into the dnf
                if i < conj_len-1 and conj_list[conj_idx][i+1] == "'":
                    dnf[conj_idx].append(-ord(conj_list[conj_idx][i]))
                else:
                    dnf[conj_idx].append(ord(conj_list[conj_idx][i]))

    return dnf

# TODO
# 1. modify for conjuctive form
# 2. modify for RPN
# 3. implement MDNF with this
# 4. add tester for equivalency of logic functions
# 5. draw tree of operations


def main():
    print("Input your logical expression in disjunctive form:\n\ny = ")
    expr = input()
    dnf = parse_dnf(expr)
    min_form = min_nand_form(dnf)
    print("The minimal form using only NANDs is:\ny =", min_form)
    nand_count = min_form.count('\u22BC');
    print("\n\nThe number of NANDs in the optimal solution is", nand_count)


if __name__ == "__main__":
    main()
