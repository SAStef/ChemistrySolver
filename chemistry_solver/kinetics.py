"""
Chemical Kinetics Problem Solver Module

This module provides functions to solve various chemical kinetics problems.
"""
import math
from fractions import Fraction

def calculate_rate_constant_from_half_life(half_life):
    """
    Calculate rate constant (k) from half-life for first-order reactions.
    
    Args:
        half_life (float): Half-life value
    
    Returns:
        float: Rate constant k
    """
    if half_life <= 0:
        raise ValueError("Half-life must be positive")
    
    # k = ln(2) / t1/2
    return math.log(2) / half_life

def calculate_half_life_from_rate_constant(rate_constant):
    """
    Calculate half-life from rate constant (k) for first-order reactions.
    
    Args:
        rate_constant (float): Rate constant value
    
    Returns:
        float: Half-life
    """
    if rate_constant <= 0:
        raise ValueError("Rate constant must be positive")
    
    # t1/2 = ln(2) / k
    return math.log(2) / rate_constant

def calculate_concentration_after_time(initial_concentration, rate_constant, time):
    """
    Calculate concentration after a specified time for first-order reactions.
    
    Args:
        initial_concentration (float): Initial concentration
        rate_constant (float): Rate constant k
        time (float): Time elapsed
    
    Returns:
        float: Final concentration
    """
    if initial_concentration < 0:
        raise ValueError("Initial concentration cannot be negative")
    if rate_constant <= 0:
        raise ValueError("Rate constant must be positive")
    if time < 0:
        raise ValueError("Time cannot be negative")
    
    # [A]t = [A]0 * e^(-kt)
    return initial_concentration * math.exp(-rate_constant * time)

def calculate_fraction_remaining(rate_constant, time):
    """
    Calculate the fraction of initial concentration remaining after time.
    
    Args:
        rate_constant (float): Rate constant k
        time (float): Time elapsed
    
    Returns:
        float: Fraction remaining (between 0 and 1)
    """
    if rate_constant <= 0:
        raise ValueError("Rate constant must be positive")
    if time < 0:
        raise ValueError("Time cannot be negative")
    
    # Fraction = e^(-kt)
    return math.exp(-rate_constant * time)

def as_simplified_fraction(decimal, max_denominator=100):
    """
    Convert a decimal to a simplified fraction if possible.
    
    Args:
        decimal (float): Decimal value to convert
        max_denominator (int): Maximum denominator to consider
    
    Returns:
        tuple or None: (numerator, denominator) if a "nice" fraction exists, None otherwise
    """
    if decimal <= 0 or decimal >= 1:
        return None
    
    try:
        frac = Fraction(decimal).limit_denominator(max_denominator)
        return (frac.numerator, frac.denominator)
    except (ValueError, ZeroDivisionError):
        return None

def solve_first_order_kinetics(half_life=None, rate_constant=None, time=None, 
                              fraction_remaining=None, initial_concentration=None, 
                              final_concentration=None):
    """
    Solve first-order kinetics problems given any two parameters.
    
    Args:
        half_life (float, optional): Half-life
        rate_constant (float, optional): Rate constant k
        time (float, optional): Time elapsed
        fraction_remaining (float, optional): Fraction of initial concentration remaining
        initial_concentration (float, optional): Initial concentration
        final_concentration (float, optional): Final concentration
    
    Returns:
        dict: Solution with all calculated parameters
    """
    # Count how many parameters are provided
    params = [half_life, rate_constant, time, fraction_remaining, 
              initial_concentration, final_concentration]
    provided = sum(p is not None for p in params)
    
    if provided < 2:
        raise ValueError("At least two parameters must be provided to solve the problem")
    
    # Solution steps
    steps = ["First-order kinetics equations:"]
    steps.append("1. k = ln(2) / t₁/₂")
    steps.append("2. t₁/₂ = ln(2) / k")
    steps.append("3. [A]t = [A]₀ × e^(-kt)")
    steps.append("4. Fraction remaining = [A]t / [A]₀ = e^(-kt)")
    steps.append("5. t = ln([A]₀/[A]t) / k = -ln(fraction) / k")
    
    # Create a result dictionary
    result = {
        "half_life": half_life,
        "rate_constant": rate_constant,
        "time": time,
        "fraction_remaining": fraction_remaining,
        "initial_concentration": initial_concentration,
        "final_concentration": final_concentration,
        "fraction_as_tuple": None,
        "steps": steps
    }
    
    # Step 1: If we have half-life but not rate constant, calculate rate constant
    if half_life is not None and rate_constant is None:
        rate_constant = calculate_rate_constant_from_half_life(half_life)
        result["rate_constant"] = rate_constant
        steps.append(f"Calculating rate constant from half-life:")
        steps.append(f"k = ln(2) / t₁/₂ = ln(2) / {half_life} = {rate_constant:.6f}")
    
    # Step 2: If we have rate constant but not half-life, calculate half-life
    elif rate_constant is not None and half_life is None:
        half_life = calculate_half_life_from_rate_constant(rate_constant)
        result["half_life"] = half_life
        steps.append(f"Calculating half-life from rate constant:")
        steps.append(f"t₁/₂ = ln(2) / k = ln(2) / {rate_constant} = {half_life:.6f}")
    
    # Step 3: Calculate fraction remaining if we have rate constant and time
    if rate_constant is not None and time is not None and fraction_remaining is None:
        fraction_remaining = calculate_fraction_remaining(rate_constant, time)
        result["fraction_remaining"] = fraction_remaining
        steps.append(f"Calculating fraction remaining from rate constant and time:")
        steps.append(f"Fraction = e^(-kt) = e^(-{rate_constant} × {time}) = {fraction_remaining:.6f}")
        
        # Try to represent as simplified fraction
        fraction_as_tuple = as_simplified_fraction(fraction_remaining)
        result["fraction_as_tuple"] = fraction_as_tuple
    
    # Step 4: Calculate time if we have rate constant and fraction remaining
    elif rate_constant is not None and fraction_remaining is not None and time is None:
        if fraction_remaining <= 0 or fraction_remaining >= 1:
            raise ValueError("Fraction remaining must be between 0 and 1 exclusively")
        
        time = -math.log(fraction_remaining) / rate_constant
        result["time"] = time
        steps.append(f"Calculating time from rate constant and fraction remaining:")
        steps.append(f"t = -ln(fraction) / k = -ln({fraction_remaining}) / {rate_constant} = {time:.6f}")
    
    # Step 5: If we have initial and final concentrations but not fraction remaining
    if initial_concentration is not None and final_concentration is not None and fraction_remaining is None:
        if initial_concentration <= 0:
            raise ValueError("Initial concentration must be positive")
        if final_concentration <= 0:
            raise ValueError("Final concentration must be positive")
        if final_concentration > initial_concentration:
            raise ValueError("Final concentration cannot be greater than initial concentration")
        
        fraction_remaining = final_concentration / initial_concentration
        result["fraction_remaining"] = fraction_remaining
        steps.append(f"Calculating fraction remaining from concentrations:")
        steps.append(f"Fraction = [A]t / [A]₀ = {final_concentration} / {initial_concentration} = {fraction_remaining:.6f}")
        
        # Try to represent as simplified fraction
        fraction_as_tuple = as_simplified_fraction(fraction_remaining)
        result["fraction_as_tuple"] = fraction_as_tuple
    
    # Step 6: Calculate final concentration if we have initial concentration and fraction remaining
    elif initial_concentration is not None and fraction_remaining is not None and final_concentration is None:
        final_concentration = initial_concentration * fraction_remaining
        result["final_concentration"] = final_concentration
        steps.append(f"Calculating final concentration:")
        steps.append(f"[A]t = [A]₀ × fraction = {initial_concentration} × {fraction_remaining} = {final_concentration:.6f}")
    
    # Step 7: Calculate initial concentration if we have final concentration and fraction remaining
    elif final_concentration is not None and fraction_remaining is not None and initial_concentration is None:
        if fraction_remaining <= 0:
            raise ValueError("Fraction remaining must be positive")
        
        initial_concentration = final_concentration / fraction_remaining
        result["initial_concentration"] = initial_concentration
        steps.append(f"Calculating initial concentration:")
        steps.append(f"[A]₀ = [A]t / fraction = {final_concentration} / {fraction_remaining} = {initial_concentration:.6f}")
    
    # Calculate any missing parameters if possible
    
    # If we now have rate constant and time but not fraction remaining
    if rate_constant is not None and time is not None and result["fraction_remaining"] is None:
        fraction_remaining = calculate_fraction_remaining(rate_constant, time)
        result["fraction_remaining"] = fraction_remaining
        steps.append(f"Calculating fraction remaining from rate constant and time:")
        steps.append(f"Fraction = e^(-kt) = e^(-{rate_constant} × {time}) = {fraction_remaining:.6f}")
        
        # Try to represent as simplified fraction
        fraction_as_tuple = as_simplified_fraction(fraction_remaining)
        result["fraction_as_tuple"] = fraction_as_tuple
    
    # If we now have initial concentration, fraction, but not final concentration
    if result["initial_concentration"] is not None and result["fraction_remaining"] is not None and result["final_concentration"] is None:
        final_concentration = result["initial_concentration"] * result["fraction_remaining"]
        result["final_concentration"] = final_concentration
        steps.append(f"Calculating final concentration:")
        steps.append(f"[A]t = [A]₀ × fraction = {result['initial_concentration']} × {result['fraction_remaining']} = {final_concentration:.6f}")
    
    return result