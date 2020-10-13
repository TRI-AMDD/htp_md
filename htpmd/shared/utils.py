"""
Utility methods for shared package.
"""


def check_params(required_parameters, passed_params):
    """
    Checks for required parameters in passed dict and raises ValueError if not found.

    Args:
        required_parameters:
        passed_params:

    """
    for param_name in required_parameters:
        if param_name not in passed_params:
            raise ValueError(f'Required parameter {param_name} missing.')
