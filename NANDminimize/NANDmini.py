# The task is to minimize the number of NANDs/NORs in the expression obtained
# by transforming the initial one, utilizing De-Morgan's laws and properties of NAND/NOR.

inf = 10 ** 18  # artificial infinity, larger than any feasible length of expression


# utility function which optimizes nand/nor usage for both cases(conjuctive and disjunctive)
# just need to define different base cases
# adv_flag tells us if we can use DeMorgan negation to our advantage
# either we are using nand on dnf, or nor on cnf, and in both cases in the second layer
def minimize_util(low_idx, high_idx, min_op, backtrack, adv_flag):
    # base case
    if low_idx == high_idx:
        return min_op[low_idx][high_idx]
    # if we've already calculated the value
    if min_op[low_idx][high_idx] != inf:
        return min_op[low_idx][high_idx]

    for split_point in range(low_idx, high_idx):
        left_cost = 2 * minimize_util(low_idx, split_point, min_op, backtrack, adv_flag) + 1
        # if one of the operands is a first layer problem, and we can exploit negation
        # then half the nands/nors dissapear, and we use it only once since there is no need to set A' = A op A
        # special case when left operand is a single variable and we cant exploit this feature
        if left_cost != 1 and low_idx == split_point and adv_flag:
            left_cost = left_cost // 4

        right_cost = 2 * minimize_util(split_point + 1, high_idx, min_op, backtrack, adv_flag) + 1
        if right_cost != 1 and split_point + 1 == high_idx and adv_flag:
            right_cost = right_cost // 4

        # calculate minimal cost for split at split point
        split_cost = left_cost + right_cost + 1
        # if found better solution, update
        if split_cost < min_op[low_idx][high_idx]:
            min_op[low_idx][high_idx] = split_cost
            backtrack[low_idx][high_idx] = split_point

    return min_op[low_idx][high_idx]

# generates optimal expression of arbitrary normal form from calculated backtracking table
def nf_expression(backtrack, expr, is_nand, is_dnf):
    # define operation character(nand/nor)
    op = "\u22BC"
    if not is_nand:
        op = "\u22BD"

    # utility function for generating negation of first layer problem
    # let's assume we have conjuctive form without loss of generality
    # we can use the negation to our advantage
    # we have the following cases:
    # 1. conjuction is one letter - we have to form it's negation A' = A nand A
    # 2. conjuction is three letters(negation), and therefore (A nand A)' = A
    # 3. conjuction is of the form A nand A(any case where both halves of A are nonempty)
    # last case is same as second, except we need to add parentheses in the end
    def factor_negate(factor):
        if len(factor) == 1:
            return "(" + factor + op + factor + ")"
        elif len(factor) == 3:
            return factor[1]
        return factor[0:len(factor) // 2]

    # returns (A op A) as string
    def negation_form(expr, op):
        if len(expr) == 1:
            return "(" + expr + op + expr + ")"
        return "((" + expr + ")" + op + "(" + expr + "))"

    # returns (A op B) as string
    def op_form(left_expr, right_expr, op):
        if len(left_expr) != 1:
            left_expr = "(" + left_expr + ")"
        if len(right_expr) != 1:
            right_expr = "(" + right_expr + ")"
        return "(" + left_expr + op + right_expr + ")"

        # utility function for recursivelly generating expression
    def expression_util(low_idx, high_idx):
        if low_idx == high_idx: # reduces to first layer problem
            return expr[low_idx]
        # AvB = (AvB)'' = (A'B')' = (A nand A) nand (B nand B)
        # AvB = (AvB)'' = (A nor B)' = (A nor B) nor (A nor B)
        # AB = (AB)'' = (A nand B)' = (A nand B) nand (A nand B)
        # AB = (AB)'' = (A' v B')' = A' nor B' = (A nor A) nor (B nor B)
        # in cases 1 and 4 we can exploit the negation on the left and right side
        left_expr = expression_util(low_idx, backtrack[low_idx][high_idx])
        right_expr = expression_util(backtrack[low_idx][high_idx] + 1, high_idx)
        # cases 1 and 4
        if is_nand == is_dnf:
            first_half = negation_form(left_expr, op) # (A op A)
            second_half = negation_form(right_expr, op)
        # cases 2 and 3
        else:
            first_half = op_form(left_expr, right_expr, op)  # (A op B)
            second_half = first_half

        if is_nand == is_dnf and low_idx == backtrack[low_idx][high_idx]:
            first_half = factor_negate(left_expr)

        if is_nand == is_dnf and high_idx == backtrack[low_idx][high_idx] + 1:
            second_half = factor_negate(right_expr)

        full_expr = first_half + op + second_half
        return full_expr

    return expression_util(0, len(expr) - 1)


# finds minimal equivalent expression of conjuction/disjunction using only NANDs/NORs
# argument is list of numbers corresponding to the indices of the input variables
def first_layer(factor, is_nand, is_dnf):
    num_var = len(factor)
    # minimal number of NANDs/NORs to represent segment in conjuctive/disjunctive factor
    min_op = [[inf] * num_var for i in range(num_var)]
    backtrack = [[-1] * num_var for i in range(num_var)]  # store optimal splitting point of each segment
    # base cases
    for i in range(num_var):
        if factor[i] > 0:
            # zero NANDs/NORs are needed to write expression with a single non-negated variable
            min_op[i][i] = 0
        else:
            # one NAND/NOR needed when there is a single negation
            min_op[i][i] = 1
        backtrack[i][i] = i

    op = "\u22BC"
    if not is_nand:
        op = "\u22BD"
    expr = [chr(x) if x > 0 else chr(-x) + op + chr(-x) for x in factor]

    minimize_util(0, num_var-1, min_op, backtrack, False)
    # the first layer problem can be considered the oposite of the second layer one
    return nf_expression(backtrack, expr, is_nand, not is_dnf)


# finds equivalent of disjunctive/conjuctive form with least number of NANDs/NORs
def second_layer(expr, cost, is_nand, is_dnf):
    num_var = len(cost)
    # minimal number of NANDs/NORs to represent segment in disjunctive/conjuctive form
    min_op = [[inf] * num_var for i in range(num_var)]
    backtrack = [[-1] * num_var for i in range(num_var)]  # store optimal splitting point of each segment
    # base cases
    for i in range(num_var):
        min_op[i][i] = cost[i]  # minimal cost of conjuction/disjunction in disjunctive/conjuctive form
        backtrack[i][i] = i

    minimize_util(0, num_var - 1, min_op, backtrack, (is_nand == is_dnf))
    return nf_expression(backtrack, expr, is_nand, is_dnf)


# the first argument is the list of lists(list of conjuctions/disjunctions)
# the inner lists contain the indexes of the logical input variables
# it is already assumed there are disjunctions between the conjuctive forms
# the second argument represents the number of variables in the initial logic function
def min_log_form(nf, is_nand, is_dnf):
    # depending on if form is dnf calculate minimum conjuction or disjunction cost individually
    expr = [first_layer(factor, is_nand, is_dnf) for factor in nf]
    cost = [len(nand_ex) for nand_ex in expr]
    # calculate optimal form given the expressions and costs for first layer problem
    form = ""
    if len(nf) > 1:
        form = second_layer(expr, cost, is_nand, is_dnf)
    elif len(nf) == 1:
        form = expr[0]
    return form



# convert DNF in string form to DNF in list form(suitable for the algorithm above)
# input should be in form "ADEvB'CvG'H" with letter v representing disjunction
# any character not 'v' or space will be considered a logical variable
# negations are specified with ' after the variable
def parse_dnf(dnf_str):
    # we first strip all whitespace
    dnf_str = "".join(dnf_str.split())
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

# convert CNF in string form to CNF in list form(suitable for the algorithm above)
# input should be in form "(AvBvC)(D'vEvF')G" with letter v representing disjunction
# any character not 'v', parentheses or space will be considered a logical variable
# negations are specified with ' after the variable
def parse_cnf(cnf_str):
    # tests if there is a conjuction between the 2 characters
    def is_lim(left, right):
        left_lim = False
        right_lim = False
        # it should be a right parentheses or a variable
        if left != "(" and left != "v" and left != " ":
            left_lim = True
        if right != ")" and right != "v" and right != "'" and right != " ":
            right_lim = True
        return (left_lim and right_lim)

    # we will turn the problem into equivalent dnf and call parse_dnf
    # remove all "v"
    cnf_str = cnf_str.replace("v", " ")
    # for all delimitters ")(", "AB", "A(" or ")A" add "v" in the middle
    new_cnf = ""
    for i in range(len(cnf_str)):
        if i > 0 and is_lim(cnf_str[i-1], cnf_str[i]):  # two delimitters next to eachother
            new_cnf += "v"
        new_cnf += cnf_str[i]
    cnf_str = new_cnf

    # remove all parentheses
    cnf_str = cnf_str.replace("(", "")
    cnf_str = cnf_str.replace(")", "")
    # remove all whitespace
    cnf_str = cnf_str.replace(" ", "")

    return parse_dnf(cnf_str)

# TODO
# 2. modify for RPN
# 4. add tester for equivalency of logic functions
# 5. draw tree of operations
# 6. improve scalability(integers, not chars, modify input output)
# 9. make input nicer, maybe file support
# 10. DNF verify


def main():
    print("Input type of operation(nand, nor): ")
    op = input().lower()
    is_nand = (op == "nand")
    print("\nInput your logical expression in conjuctive or disjunctive form:\n\ny = ")
    expr = input()
    if expr.count('(') > 0:
        nf = parse_cnf(expr)
        min_form = min_log_form(nf, is_nand, False)
    else:
        nf = parse_dnf(expr)
        min_form = min_log_form(nf, is_nand, True)
    print("The minimal form is:\n\ny =", min_form)
    if is_nand:
        print("\nThe number of NANDs in the optimal solution is", min_form.count('\u22BC'))
    else:
        print("\nThe number of NORs in the optimal solution is", min_form.count('\u22BD'))


if __name__ == "__main__":
    main()