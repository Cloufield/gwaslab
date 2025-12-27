"""
Status code manipulation functions for GWASLab.

Status codes are 7-digit integers where:
- Digit 1-2: Build number (e.g., 13 for hg13, 19 for hg19, 38 for hg38)
- Digit 3-7: Status flags for variant processing

All functions use integer arithmetic for optimal performance.
"""
from typing import Union, Tuple, Optional, List
import pandas as pd
import re

# Constants
STATUS_CODE_LENGTH = 7
STATUS_CATEGORIES = [
    str(j + i) 
    for j in [1300000, 1800000, 1900000, 3800000, 3900000, 9700000, 9800000, 9900000] 
    for i in range(0, 100000)
]

# Integer dtype constants for type checking
INTEGER_DTYPES = ['int64', 'Int64', 'int32', 'Int32']


# ============================================================================
# Helper Functions
# ============================================================================

def _normalize_to_integer_series(status: Union[int, pd.Series, List[int]]) -> pd.Series:
    """
    Convert status to integer Series, handling various input types.
    
    Parameters:
    -----------
    status : int, Series, or array-like
        Status code(s) to normalize
        
    Returns:
    --------
    pd.Series
        Integer Series with status codes
    """
    if isinstance(status, pd.Series):
        if status.dtype in INTEGER_DTYPES:
            return status
        elif status.dtype.name == 'category':
            # Convert categorical to string first, then to integer
            return status.astype(str).astype(int)
        else:
            return status.astype(int)
    else:
        return pd.Series(status).astype(int)


def ensure_status_int(sumstats: pd.DataFrame, status_col: str = "STATUS") -> pd.DataFrame:
    """
    Ensure STATUS column is integer type (Int64), converting from categorical or other types if needed.
    
    Parameters:
    -----------
    sumstats : pd.DataFrame
        DataFrame with STATUS column
    status_col : str, default="STATUS"
        Name of the status column
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with STATUS column as Int64
    """
    if status_col in sumstats.columns:
        if sumstats[status_col].dtype.name == 'category':
            sumstats[status_col] = sumstats[status_col].astype(str).astype('Int64')
        elif sumstats[status_col].dtype not in ['int64', 'Int64', 'int32', 'Int32']:
            sumstats[status_col] = sumstats[status_col].astype('Int64')
        else:
            # Ensure it's Int64 (nullable integer)
            sumstats[status_col] = sumstats[status_col].astype('Int64')
    return sumstats


def _calculate_powers(digit: int) -> Tuple[int, int]:
    """
    Calculate powers of 10 for extracting/replacing digits at position.
    
    Parameters:
    -----------
    digit : int
        Digit position (1-indexed, 1=leftmost, 7=rightmost)
        
    Returns:
    --------
    tuple
        (power_left, power_right) for digit extraction
    """
    power_right = 10**(STATUS_CODE_LENGTH - digit)
    power_left = 10**(STATUS_CODE_LENGTH - digit + 1)
    return power_left, power_right


def _extract_digit(status_int: Union[int, pd.Series], digit: int) -> Union[int, pd.Series]:
    """
    Extract a specific digit from integer status code(s).
    
    Parameters:
    -----------
    status_int : int or Series
        Status code(s) as integer
    digit : int
        Digit position (1-indexed)
        
    Returns:
    --------
    int or Series
        Extracted digit(s) at position
    """
    power_left, power_right = _calculate_powers(digit)
    return (status_int // power_right) % 10


def _parse_pattern_for_digit_constraints(pattern: str) -> Optional[List[Tuple[int, List[int]]]]:
    """
    Parse regex pattern to extract digit position constraints.
    
    Parameters:
    -----------
    pattern : str
        Regex pattern like r'\\w\\w\\w\\w\\w[35]\\w'
        
    Returns:
    --------
    tuple or None
        (digit_positions, allowed_digits_list) if parseable, None otherwise
    """
    pattern_clean = pattern.strip('^$')
    
    # Check for literal digits (not in character classes) - if found, fall back to string matching
    # This must be done BEFORE parsing, as literal digits can't be handled by digit-based matching
    # Split by character classes and check the last part for literal digits
    parts = re.split(r'\[[0-9]+\]', pattern_clean)
    if len(parts) > 1:
        # Check the last part (after the last character class)
        last_part = parts[-1]
        # Remove all wildcards from the last part
        remaining = re.sub(r'\\\\w', '', last_part)
        # If there are any digits remaining, we have literal digits
        if remaining and re.search(r'[0-9]', remaining):
            # Has literal digits that can't be parsed - return None to force string matching
            return None
    
    digit_matches = list(re.finditer(r'\[([0-9]+)\]', pattern_clean))
    wildcard_parts = re.split(r'\[[0-9]+\]', pattern_clean)
    
    if not digit_matches or len(wildcard_parts) != len(digit_matches) + 1:
        return None
    
    constraints = []
    pos = 0
    
    for wildcard, match in zip(wildcard_parts, digit_matches):
        # Count wildcard patterns (\w)
        wildcard_count = len(re.findall(r'\\w', wildcard))
        pos += wildcard_count
        digit_pos = pos + 1
        
        if digit_pos < 1 or digit_pos > STATUS_CODE_LENGTH:
            return None
        
        allowed_digits = [int(d) for d in match.group(1)]
        constraints.append((digit_pos, allowed_digits))
        pos = digit_pos
    
    return constraints


# ============================================================================
# Core Status Manipulation Functions
# ============================================================================

def status_match(status: Union[int, pd.Series], digit: int, to_match: Union[int, List[int]]) -> Union[bool, pd.Series]:
    """
    Check if a specific digit in status code(s) matches given value(s).
    
    Parameters:
    -----------
    status : int or Series
        Status code(s) to check
    digit : int
        Digit position (1-indexed, 1=leftmost, 7=rightmost)
    to_match : int or list of int
        Digit value(s) to match against
        
    Returns:
    --------
    bool or Series
        True where digit matches, False otherwise
        
    Examples:
    --------
    >>> status_match(1234567, 6, [3, 5])
    False
    >>> status_match(1234537, 6, [3, 5])
    True
    """
    status_int = _normalize_to_integer_series(status)
    power_left, power_right = _calculate_powers(digit)
    middle = (status_int // power_right) % 10
    
    if len(to_match) == 1:
        return middle == to_match[0]
    else:
        to_match_set = set(to_match) if not isinstance(to_match, set) else to_match
        return middle.isin(to_match_set)


def vchange_status(status: Union[int, pd.Series], digit: int, before: Union[str, List[str]], after: Union[str, List[str]]) -> Union[int, pd.Series]:
    """
    Change specific digits in status code(s) based on mapping.
    
    Replaces digits at position 'digit' that match 'before' with corresponding
    values from 'after'. Uses integer arithmetic for optimal performance.
    
    Parameters:
    -----------
    status : int or Series
        Status code(s) to modify
    digit : int
        Digit position to change (1-indexed, 1=leftmost, 7=rightmost)
    before : str or list
        Digit values to match (e.g., "45" or ["4", "5"])
    after : str or list
        Replacement digit values (e.g., "12" or ["1", "2"])
        
    Returns:
    --------
    int or Series
        Modified status code(s) as integers
        
    Examples:
    --------
    >>> vchange_status(1234547, 6, "45", "12")
    1234517  # digit 6 changed from 4 to 1
    """
    # Normalize inputs
    if isinstance(before, list):
        before = ''.join(str(b) for b in before)
    if isinstance(after, list):
        after = ''.join(str(a) for a in after)
    
    before_str = str(before)
    after_str = str(after)
    
    # Early exit if no changes needed
    if before_str == after_str:
        return status if isinstance(status, pd.Series) else pd.Series(status)
    
    # Build sparse mapping dictionary (only non-identity mappings)
    mapping = {}
    for b_char, a_char in zip(before_str, after_str):
        b_int, a_int = int(b_char), int(a_char)
        if b_int != a_int:
            mapping[b_int] = a_int
    
    # Early exit if no actual mappings
    if not mapping:
        return status if isinstance(status, pd.Series) else pd.Series(status)
    
    # Convert status to integer Series
    status_int = _normalize_to_integer_series(status)
    
    # Extract components
    power_left, power_right = _calculate_powers(digit)
    prefix = status_int // power_left
    middle = _extract_digit(status_int, digit)
    suffix = status_int % power_right
    
    # Apply mapping
    middle_replaced = middle.map(mapping).fillna(middle).astype('int64')
    
    # Reconstruct status code
    result = prefix * power_left + middle_replaced * power_right + suffix
    return result


def set_status_digit(status: Union[int, pd.Series], digit: int, value: Union[int, str]) -> Union[int, pd.Series]:
    """
    Set a specific digit in status code(s) to a given value.
    
    Parameters:
    -----------
    status : int or Series
        Status code(s) to modify
    digit : int
        Digit position (1-indexed, 1=leftmost, 7=rightmost)
    value : int or str
        New digit value (0-9)
        
    Returns:
    --------
    int or Series
        Status code(s) with digit replaced (returns same type as input)
        
    Examples:
    --------
    >>> set_status_digit(1234567, 6, 0)
    1234507
    """
    # Preserve input type (scalar vs Series)
    is_scalar = not isinstance(status, pd.Series)
    
    if is_scalar:
        status_int = int(status)
    else:
        status_int = _normalize_to_integer_series(status)
    
    value = int(value)
    
    power_left, power_right = _calculate_powers(digit)
    
    if is_scalar:
        prefix = status_int // power_left
        suffix = status_int % power_right
        result = prefix * power_left + value * power_right + suffix
        return result
    else:
        prefix = status_int // power_left
        suffix = status_int % power_right
        result = prefix * power_left + value * power_right + suffix
        return result


def copy_status(from_status: Union[int, pd.Series], to_status: Union[int, pd.Series], digit: int) -> Union[int, pd.Series]:
    """
    Copy a specific digit from one status code to another.
    
    Parameters:
    -----------
    from_status : int or Series
        Source status code(s)
    to_status : int or Series
        Target status code(s)
    digit : int
        Digit position to copy (1-indexed)
        
    Returns:
    --------
    int or Series
        to_status with digit replaced from from_status
        
    Examples:
    --------
    >>> copy_status(1234567, 9999999, 6)
    9999969  # digit 6 (value 6) copied from first to second
    """
    from_status_int = _normalize_to_integer_series(from_status)
    to_status_int = _normalize_to_integer_series(to_status)
    
    power_left, power_right = _calculate_powers(digit)
    
    # Extract digit from source
    from_digit = _extract_digit(from_status_int, digit)
    
    # Extract prefix and suffix from target
    to_prefix = to_status_int // power_left
    to_suffix = to_status_int % power_right
    
    # Reconstruct with copied digit
    result = to_prefix * power_left + from_digit * power_right + to_suffix
    return result


def get_status_prefix(status: Union[int, pd.Series], digit: int) -> Union[int, pd.Series]:
    """
    Get the prefix (left part) of status code before a specific digit.
    
    Parameters:
    -----------
    status : int or Series
        Status code(s)
    digit : int
        Digit position (1-indexed)
        
    Returns:
    --------
    int or Series
        Prefix part (as integer, needs to be multiplied by appropriate power)
        Returns same type as input
        
    Examples:
    --------
    >>> get_status_prefix(1234567, 4)
    123  # left part before digit 4
    """
    is_scalar = not isinstance(status, pd.Series)
    if is_scalar:
        status_int = int(status)
        power_left, _ = _calculate_powers(digit)
        return status_int // power_left
    else:
        status_int = _normalize_to_integer_series(status)
        power_left, _ = _calculate_powers(digit)
        return status_int // power_left


def get_status_suffix(status: Union[int, pd.Series], digit: int) -> Union[int, pd.Series]:
    """
    Get the suffix (right part) of status code after a specific digit.
    
    Parameters:
    -----------
    status : int or Series
        Status code(s)
    digit : int
        Digit position (1-indexed)
        
    Returns:
    --------
    int or Series
        Suffix part (as integer)
        Returns same type as input
        
    Examples:
    --------
    >>> get_status_suffix(1234567, 4)
    567  # right part after digit 4
    """
    is_scalar = not isinstance(status, pd.Series)
    if is_scalar:
        status_int = int(status)
        _, power_right = _calculate_powers(digit)
        return status_int % power_right
    else:
        status_int = _normalize_to_integer_series(status)
        _, power_right = _calculate_powers(digit)
        return status_int % power_right


# ============================================================================
# Pattern Matching Functions
# ============================================================================

def match_status(status: Union[int, pd.Series, List[int]], pattern: str, na: bool = False) -> Union[bool, pd.Series]:
    """
    Match status codes against a regex pattern using integer arithmetic.
    
    For simple digit-based patterns, uses fast integer arithmetic.
    For complex patterns, falls back to string matching.
    
    Parameters:
    -----------
    status : int, Series, or array-like
        Status code(s) to match
    pattern : str
        Regex pattern (e.g., r'\\w\\w\\w\\w\\w[35]\\w' for digit 6 in [3,5])
    na : bool
        How to handle NA values (default: False = treat as non-matching)
        
    Returns:
    --------
    bool or Series
        True where pattern matches, False otherwise
        
    Examples:
    --------
    >>> match_status(1234537, r'\\w\\w\\w\\w\\w[35]\\w')
    True
    >>> match_status(1234567, r'\\w\\w\\w\\w\\w[35]\\w')
    False
        
    Note:
    -----
    For best performance with integer status codes, use status_match() directly
    for simple digit checks, or combine multiple status_match() calls with & or |.
    """
    # Normalize input to Series
    if not isinstance(status, pd.Series):
        status = pd.Series(status)
    
    # Fast path: if already integer, use directly
    if status.dtype in INTEGER_DTYPES:
        status_int = status
    else:
        # Handle Categorical or other dtypes
        status_str = status.astype(str)
        
        # Try to convert to integer, fall back to string matching if fails
        try:
            status_int = status_str.astype(int)
        except (ValueError, TypeError):
            # Fall back to string matching for non-numeric values
            status_str = status_str.str.zfill(STATUS_CODE_LENGTH)
            pat = f"^{pattern}$"
            return status_str.str.match(pat, na=na)
    
    # Try to parse pattern for digit-based matching (faster)
    constraints = _parse_pattern_for_digit_constraints(pattern)
    
    if constraints:
        # Use integer arithmetic for digit-based matching
        result = pd.Series(True, index=status_int.index, dtype=bool)
        
        for digit_pos, allowed_digits in constraints:
            digit_match = status_match(status_int, digit_pos, allowed_digits)
            result = result & digit_match
        
        if na:
            result = result.fillna(False)
        return result
    
    # Fall back to string matching for complex patterns
    status_str = status_int.astype(str).str.zfill(STATUS_CODE_LENGTH)
    pat = f"^{pattern}$"
    return status_str.str.match(pat, na=na)


# ============================================================================
# Legacy Functions (kept for backward compatibility)
# ============================================================================

def change_status(status: int, digit: int, after: int) -> int:
    """
    Legacy function: Set a specific digit in status code.
    
    Note: Prefer set_status_digit() for new code.
    """
    power_left, power_right = _calculate_powers(digit)
    prefix = status // power_left
    suffix = status % power_right
    return prefix * power_left + after * power_right + suffix


def schange_status(status: pd.Series, digit: int, after: int) -> pd.Series:
    """
    Legacy function: Set digit using pandas eval (slower).
    
    Note: Prefer set_status_digit() for new code.
    """
    prefix = status.floordiv(10**(STATUS_CODE_LENGTH - digit + 1))
    suffix = status.mod(10**(STATUS_CODE_LENGTH - digit))
    return pd.eval("prefix*10**(7-digit+1) + after*10**(7-digit) + suffix")
