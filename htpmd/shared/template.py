"""
Template Module

Use this template to add a new shared computation module.

Guidelines:
    1. Name module after property/metric that is being computed.
    2. Name methods intuitively after what they compute and how they operate.
    3. Fill out docstring in method to specify usage.

"""


def my_new_method(trajectory, **params):
    """
    Description:
        <Enter description of the method here>

        Example:
            <Provide one or more examples on how to call the method>

    Version: <use semantic versioning for your method e.g. 1.0.0>

    Author:
        Name:                                           <name>
        Affiliation:                                    <optional, affiliation>
        Email:                                          <optional, email>

    Args:
        trajectory (trajectory.base.Trajectory):        trajectory to compute metric on
        **params:                                       Methodology specific parameters.
                                                        Required fields:
                                                            <name> (<type>)

    Returns:
        <type>:                                         <description>
    """
