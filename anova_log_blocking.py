"""anova_log_blocking.py

Blocking or zoning of well log curve data using the analyis of variance technique.

reference:
- Al-Adani, Nabil, 2012, Data blockign or zoning: Well-log-data application: Journal of Canadian Petroleum Technology.
- Gill, D., 1970, Application of a statistical zonation method to reservoir evaluation adn digitized log analysis: AAPG Bulletin, 54 (5), 719-729
"""

# standard libs
import copy
import statistics
import random

# 3rd party
import numpy
import binarytree

# local

def _mean_variance_within_zone(zone1:numpy.ndarray, zone2:numpy.ndarray):
    m1 = numpy.mean(zone1)
    m2 = numpy.mean(zone2)
    n1 = len(zone1)
    n2 = len(zone2)
    a = numpy.sum(numpy.square(zone1 - m1))
    b = numpy.sum(numpy.square(zone2 - m2))
    return (a + b) / (n1 + n2 - 2)

def _mean_variance_among_zones(zone1:numpy.ndarray, zone2:numpy.ndarray):
    m1 = numpy.mean(zone1)
    m2 = numpy.mean(zone2)
    n1 = len(zone1)
    n2 = len(zone2)
    m0 = (numpy.sum(zone1) + numpy.sum(zone2)) / (n1 + n2)  # overall average
    return n1 * (numpy.square(m1-m0)) + n2 * (numpy.square(m2-m0))

def _anova_breakpoint(arr:numpy.ndarray, min_samples_in_zone:int):
    """determine the optimal breakpoint, i.e. the index with the largest ratio of variances.
    """
    if len(arr) < 2*min_samples_in_zone:
        kbest = None
    else:
        kbest = min_samples_in_zone  
        rbest = 0
        for k in range(min_samples_in_zone,len(arr)-min_samples_in_zone):
            z1 = arr[:k]
            z2 = arr[k:]
            if _mean_variance_among_zones(zone1=z1, zone2=z2) != 0.0:
                ratio_of_variances = 1 - (_mean_variance_within_zone(zone1=z1, zone2=z2) / _mean_variance_among_zones(zone1=z1, zone2=z2))
                if ratio_of_variances > rbest:
                    rbest = ratio_of_variances
                    kbest = k
    return kbest

def _anova_recursive_tree_build(node:binarytree.Node, a:numpy.ndarray, min_samples_in_zone:int):
    """anova
    recursive tree building
    """
    knot = node.value  # parent node value
    k = _anova_breakpoint(arr=a, min_samples_in_zone=min_samples_in_zone)
    if k:
        node.left = _anova_recursive_tree_build(node=binarytree.Node(value=knot), a=a[:k], min_samples_in_zone=min_samples_in_zone)
        node.right = _anova_recursive_tree_build(node=binarytree.Node(value=k+knot), a=a[k:], min_samples_in_zone=min_samples_in_zone)
    return node
        

def anova_zoning(input_array:numpy.ndarray, min_samples_in_zone=2):
    """anova zoning
    Analyis of Variance (ANOVA) ...

    references
    - Al-Adani, Nabil, 2012, Data Blocking or Zoning: Well-Log-Data Application: Journal of Canadian Petroleum Technology.
    """
    output_array = numpy.zeros_like(input_array)

    root = _anova_recursive_tree_build(node=binarytree.Node(value=0), a=input_array, min_samples_in_zone=min_samples_in_zone)
    # print("binary tree for anova zoning")
    # print(root)
    # NOTE: leaf node values correspond to breakpoints
    # print(root.leaves)
    breakpoints = []
    for leaf in root.postorder:
        if not leaf.left and not leaf.right:
            breakpoints.append(leaf.value)
    breakpoints.append(len(input_array))
    for n in range(len(breakpoints)-1):
        wstart = breakpoints[n]
        wend = breakpoints[n+1]
        output_array[wstart:wend] = numpy.median(a=input_array[wstart:wend])
    return output_array

# utilities
def number_of_zones(arr:numpy.ndarray) -> int:
    d = numpy.diff(arr)
    n = numpy.count_nonzero(d)
    return int(n) + 1
