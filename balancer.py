from sympy import Matrix, lcm
from fractions import Fraction
from utils import parse_term, split_equation, gcd_list
from molar_mass import check_molmass
import molmass

def count_elements(formula):
    check_molmass()
    return {s: i.count for s, i in molmass.Formula(formula).composition().items()}

def parse_chemical_equation(equation):
    lhs, rhs = split_equation(equation)
    reactants = [parse_term(r.strip()) for r in lhs.split('+')]
    products  = [parse_term(p.strip()) for p in rhs.split('+')]
    return reactants, products

def balance_equation(reactants, products):
    all_elements = sorted(set(e for _, f in reactants + products for e in count_elements(f)))
    n_r, n_p = len(reactants), len(products)

    matrix = []
    for element in all_elements:
        row = [count_elements(f).get(element, 0) for _, f in reactants]
        row += [-count_elements(f).get(element, 0) for _, f in products]
        matrix.append(row)

    A = Matrix(matrix)
    null = A.nullspace()
    if not null:
        raise ValueError("Could not balance the equation.")

    solution = null[0]
    fractions = [Fraction(float(x)).limit_denominator() for x in solution]
    lcm_val = lcm([f.denominator for f in fractions])
    coeffs = [int(f * lcm_val) for f in fractions]
    if coeffs[0] < 0:
        coeffs = [-c for c in coeffs]
    divisor = gcd_list(coeffs)
    coeffs = [c // divisor for c in coeffs]

    balanced_reactants = [(coeffs[i], reactants[i][1]) for i in range(n_r)]
    balanced_products  = [(coeffs[i + n_r], products[i][1]) for i in range(n_p)]
    return balanced_reactants, balanced_products

def format_balanced_equation(reactants, products):
    def fmt(side):
        return ' + '.join([f"{'' if c == 1 else c}{f}" for c, f in side])
    return f"{fmt(reactants)} → {fmt(products)}"
