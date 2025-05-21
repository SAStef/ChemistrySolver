import math
from typing import Dict, List, Tuple, Any, Union
from fractions import Fraction

def calculate_rate_constant_from_half_life(half_life: float) -> float:
    """
    Calculate the rate constant from half-life for a first-order reaction.
    
    Parameters:
        half_life (float): The half-life of the reaction
        
    Returns:
        float: The rate constant k
    """
    # For first-order reactions: k = ln(2) / t_half
    return math.log(2) / half_life

def calculate_half_life_from_rate_constant(rate_constant: float) -> float:
    """
    Calculate the half-life from rate constant for a first-order reaction.
    
    Parameters:
        rate_constant (float): The rate constant k
        
    Returns:
        float: The half-life
    """
    # For first-order reactions: t_half = ln(2) / k
    return math.log(2) / rate_constant

def calculate_concentration_after_time(
    initial_concentration: float, 
    rate_constant: float, 
    time: float
) -> float:
    """
    Calculate the concentration after a given time using the first-order rate law.
    
    Parameters:
        initial_concentration (float): The initial concentration [A]₀
        rate_constant (float): The rate constant k
        time (float): The time elapsed
        
    Returns:
        float: The concentration at time t
    """
    # For first-order reactions: [A]t = [A]₀ * e^(-kt)
    return initial_concentration * math.exp(-rate_constant * time)

def calculate_fraction_remaining(
    rate_constant: float, 
    time: float
) -> float:
    """
    Calculate the fraction of the initial concentration remaining after a given time.
    
    Parameters:
        rate_constant (float): The rate constant k
        time (float): The time elapsed
        
    Returns:
        float: The fraction remaining at time t
    """
    # Fraction remaining = [A]t/[A]₀ = e^(-kt)
    return math.exp(-rate_constant * time)

def as_simplified_fraction(decimal: float, tolerance: float = 1e-10) -> Tuple[int, int]:
    """
    Convert a decimal to a simplified fraction with the format (numerator, denominator).
    
    Parameters:
        decimal (float): The decimal number to convert
        tolerance (float): The tolerance for approximation
        
    Returns:
        Tuple[int, int]: (numerator, denominator)
    """
    # Use Python's Fraction class with a limit
    f = Fraction(decimal).limit_denominator(1000)
    return f.numerator, f.denominator

def solve_first_order_kinetics(
    half_life: float = None,
    rate_constant: float = None,
    initial_concentration: float = None,
    final_concentration: float = None,
    time: float = None,
    fraction_remaining: float = None
) -> Dict[str, Any]:
    """
    Solve a first-order kinetics problem with any of the given parameters.
    
    Parameters:
        half_life (float, optional): The half-life of the reaction
        rate_constant (float, optional): The rate constant k
        initial_concentration (float, optional): The initial concentration [A]₀
        final_concentration (float, optional): The final concentration [A]₍
        time (float, optional): The time elapsed
        fraction_remaining (float, optional): The fraction remaining
        
    Returns:
        Dict[str, Any]: Dictionary with solution and steps
    """
    steps = []
    
    # Step 1: Calculate rate constant if not provided
    if rate_constant is None and half_life is not None:
        rate_constant = calculate_rate_constant_from_half_life(half_life)
        steps.append(f"Step 1: Calculate the rate constant (k) from the half-life:")
        steps.append(f"k = ln(2) / t₁/₂")
        steps.append(f"k = ln(2) / {half_life}")
        steps.append(f"k = {math.log(2):.6f} / {half_life}")
        steps.append(f"k = {rate_constant:.6f}")
    elif half_life is None and rate_constant is not None:
        half_life = calculate_half_life_from_rate_constant(rate_constant)
        steps.append(f"Step 1: Calculate the half-life from the rate constant:")
        steps.append(f"t₁/₂ = ln(2) / k")
        steps.append(f"t₁/₂ = ln(2) / {rate_constant}")
        steps.append(f"t₁/₂ = {math.log(2):.6f} / {rate_constant}")
        steps.append(f"t₁/₂ = {half_life:.6f}")
    
    # Step 2: Calculate remaining values based on what's provided
    if time is not None and fraction_remaining is None:
        fraction_remaining = calculate_fraction_remaining(rate_constant, time)
        steps.append(f"\nStep 2: Calculate the fraction remaining after {time} time units:")
        steps.append(f"fraction remaining = e^(-kt)")
        steps.append(f"fraction remaining = e^(-{rate_constant:.6f} × {time})")
        steps.append(f"fraction remaining = e^(-{rate_constant * time:.6f})")
        steps.append(f"fraction remaining = {fraction_remaining:.6f}")
        
        # Convert to simplified fraction
        num, den = as_simplified_fraction(fraction_remaining)
        if den <= 100:  # Only present as fraction if denominator is reasonable
            steps.append(f"As a simplified fraction: {num}/{den}")
    
    elif fraction_remaining is not None and time is None:
        time = -math.log(fraction_remaining) / rate_constant
        steps.append(f"\nStep 2: Calculate the time required to reach {fraction_remaining:.6f} of the initial amount:")
        steps.append(f"fraction remaining = e^(-kt)")
        steps.append(f"{fraction_remaining} = e^(-{rate_constant:.6f} × t)")
        steps.append(f"ln({fraction_remaining}) = -{rate_constant:.6f} × t")
        steps.append(f"t = -ln({fraction_remaining}) / {rate_constant:.6f}")
        steps.append(f"t = {-math.log(fraction_remaining):.6f} / {rate_constant:.6f}")
        steps.append(f"t = {time:.6f}")
    
    # Handle concentration calculations if provided
    if initial_concentration is not None:
        if final_concentration is None and fraction_remaining is not None:
            final_concentration = initial_concentration * fraction_remaining
            steps.append(f"\nStep 3: Calculate the final concentration:")
            steps.append(f"[A]t = [A]₀ × fraction remaining")
            steps.append(f"[A]t = {initial_concentration} × {fraction_remaining:.6f}")
            steps.append(f"[A]t = {final_concentration:.6f}")
        elif final_concentration is not None and fraction_remaining is None:
            fraction_remaining = final_concentration / initial_concentration
            steps.append(f"\nStep 3: Calculate the fraction remaining:")
            steps.append(f"fraction remaining = [A]t / [A]₀")
            steps.append(f"fraction remaining = {final_concentration} / {initial_concentration}")
            steps.append(f"fraction remaining = {fraction_remaining:.6f}")
            
            # Convert to simplified fraction
            num, den = as_simplified_fraction(fraction_remaining)
            if den <= 100:  # Only present as fraction if denominator is reasonable
                steps.append(f"As a simplified fraction: {num}/{den}")
    
    # Half-life calculation for the given time
    if time is not None and half_life is not None:
        half_lives_elapsed = time / half_life
        steps.append(f"\nAdditional Information:")
        steps.append(f"Number of half-lives elapsed: {time} / {half_life} = {half_lives_elapsed:.4f}")
        steps.append(f"Expected fraction remaining after {half_lives_elapsed:.4f} half-lives: (1/2)^{half_lives_elapsed:.4f} = {0.5 ** half_lives_elapsed:.6f}")
    
    return {
        "half_life": half_life,
        "rate_constant": rate_constant,
        "initial_concentration": initial_concentration,
        "final_concentration": final_concentration,
        "time": time,
        "fraction_remaining": fraction_remaining,
        "fraction_as_tuple": as_simplified_fraction(fraction_remaining) if fraction_remaining is not None else None,
        "steps": steps
    }