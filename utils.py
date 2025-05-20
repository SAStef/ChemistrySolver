import re

def split_equation(equation):
    for sep in ['→', '->', '-->']:
        if sep in equation:
            return [s.strip() for s in equation.split(sep)]
    raise ValueError("Invalid equation format. Use '->' or '→'.")

def parse_term(term):
    term = term.strip()
    match = re.match(r'^(\d+)\s*(.+)$', term)
    return (int(match.group(1)), match.group(2)) if match else (1, term)

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def gcd_list(numbers):
    result = numbers[0]
    for num in numbers[1:]:
        result = gcd(result, num)
    return result
