=======================
ANOVA WELL LOG BLOCKING
=======================

Blocking well log curve data is a method of upscaling high resolution petrophysical data (relative to conventional seismic acquisition data) to a lower resolution by intelligently blocking zones of geologic packages. Each geologic package is expected have similar properties and the objective is to remove small variations within the package and replace with a single value producing a "blocked" curve.

For this example, a well from the `Teapot`_ open source dataset has been used. Information from the LAS file is given below.

.. _Teapot: https://wiki.seg.org/wiki/Teapot_dome_3D_survey

.. code-block::

    Mnemonic  Unit  Value                       Description             
    --------  ----  -----                       -----------             
    STRT      F     4332.0                      START DEPTH             
    STOP      F     449.0                       STOP DEPTH              
    STEP      F     -1.0                        STEP                    
    NULL            -999.25                     NULL VALUE              
    COMP            DEPARTMENT OF ENERGY        COMPANY                 
    WELL            48X-28                      WELL                    
    FLD             NAVAL PETROLEUM RESERVE #3  FIELD                   
    LOC             490' FSL, 2449' FWL         LOCATION                
    CNTY            NATRONA                     COUNTY                  
    STAT            WYOMING                     STATE                   
    CTRY                                        COUNTRY                 
    API             49-025-23195                API NUMBER              
    UWI                                         UNIQUE WELL ID          
    DATE            15-Mar-2004                 LOG DATE {DD-MMM-YYYY}  
    SRVC            Schlumberger                SERVICE COMPANY        
    LATI      DEG                               LATITUDE                
    LONG      DEG                               LONGITUDE               
    GDAT                                        GeoDetic Datum          
    SECT            28                          Section                 
    RANG            78W                         Range                   
    TOWN            39N                         Township                

    Mnemonic  Unit  Value          Description                                                  
    --------  ----  -----          -----------                                                  
    RUN             1              RUN NUMBER                                                   
    PDAT            GROUND LEVEL   Permanent Datum                                              
    EPD:1     F     5105.0         Elevation of Permanent Datum above Mean Sea Level            
    EPD:2     F     5105.0         Elevation of tool zero above Mean Sea Level                  
    LMF             KELLY BUSHING  Logging Measured From (Name of Logging Elevation Reference)  
    APD                            Elevation of Depth Reference (LMF) above Permanent Datum     

These examples use the Gamma Ray curve, but any curve could be used.

.. code-block:: python

    df = pandas.read_csv('data/490252319500_GammaRay.csv')
    df = df.set_index('DEPT')

    depth = df.index.to_numpy()
    curve = df['GR'].to_numpy()

Constant Thickness
------------------

The simplest method of log blocking is to use a constant thickness. This method replaces each window of constant thickness, or constant number of samples in this case, with the median value in that window.

.. code-block:: python

    def constant_thickness_zoning(input_array:numpy.ndarray, nsamples:int):
        output_array = numpy.zeros_like(input_array)

        for n in range(0,len(input_array),nsamples):
            output_array[n:n+nsamples] = numpy.median(input_array[n:n+nsamples])
        # handle the last samples
        if n+nsamples < len(input_array):
            output_array[n+nsamples:] = numpy.median(input_array[n+nsamples:])
    return output_array

For example, window lengths of 3, 6, 12, and 24 are shown below:

.. image:: images/constant_thickness_blocking.png
    :alt: constant thickness log blocking example

Analysis of Variance
--------------------
A more robust approach is to use the analysis of variance statistical method to determine the index at which two zones within a region are statistically different. Al-Adani (2012) provides a basic summary: 

    1. Select a zone break point to divide into two new zones. Each zone should include at least two sample data.
    2. Calculate the *mean variance within zones (MVWZ)* and *mean variance among zones (MVAZ)*
    3. Compute the *ratio of variances (R)*

The mean variance within zones is defined as:

:math: `MVWZ = \frac{\sum_{i}^{n_1}\left ( X_i-\overline{X_1} \right )^{2}+\sum_{i}^{n_2}\left ( X_i-\overline{X_2} \right )^{2}}{n_1+n_2-2}`
    

and written in python as a function

.. code-block:: python

    def _mean_variance_within_zone(zone1:numpy.ndarray, zone2:numpy.ndarray):
        m1 = numpy.mean(zone1)
        m2 = numpy.mean(zone2)
        n1 = len(zone1)
        n2 = len(zone2)
        a = numpy.sum(numpy.square(zone1 - m1))
        b = numpy.sum(numpy.square(zone2 - m2))
        return (a + b) / (n1 + n2 - 2)

The mean variance among zones is defined as:

:math: `MVAZ = n_1\left ( \overline{X_1}-\overline{X} \right )^{2}+n_2\left ( \overline{X_2}-\overline{X} \right )^{2}`


and written in python as a function

.. code-block::python

    def _mean_variance_among_zones(zone1:numpy.ndarray, zone2:numpy.ndarray):
        m1 = numpy.mean(zone1)
        m2 = numpy.mean(zone2)
        n1 = len(zone1)
        n2 = len(zone2)
        m0 = (numpy.sum(zone1) + numpy.sum(zone2)) / (n1 + n2)  # overall average
        return n1 * (numpy.square(m1-m0)) + n2 * (numpy.square(m2-m0))

To determine the breakpoint, all possible "splits" or division into two zones are tested. The breakpoint is the index with the largest ratio of variances, defined as:

:math: `R = 1 - \frac{MVWZ}{MVAZ}`

The python function below allows for an additional paramter to be set which defines the minimum number of samples in window or zone, i.e. no zones should be smaller than this parameter.

.. code-block:: python

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

**HERE IS SOMETHING DIFFERENT**
In reading through the procedure, particularly the first step *Select a zone break point to divide into two new zones*, one may postulate the best data structure for this is a binary tree. The implementation here recursively builds a binary tree (using the third party library and open source project `binarytree`_) where the leaf nodes are the breakpoints in order from left to right. 

.. _binarytree: https://binarytree.readthedocs.io/en/main/index.html
    
implemented in python

.. code-block:: python

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

For example:

.. code-block:: python

    df = pandas.read_csv('data/490252319500_GammaRay.csv')
    df = df.set_index('DEPT')

    # if missing values, use forward fill and backward fill
    if df.isnull().values.any():
        df = df.fillna(method="ffill")
        df = df.fillna(method="bfill")

    depth = df.index.to_numpy()
    curve = df['GR'].to_numpy()

    # select a smaller depth range
    z_lower_indx = numpy.argwhere(depth > 2800)
    idx0 = z_lower_indx[0,0]
    z_upper_indx = numpy.argwhere(depth < 2900)
    idx1 = z_upper_indx[-1,0]

    depth = depth[idx0:idx1]
    curve = curve[idx0:idx1]

    root = _anova_recursive_tree_build(node=binarytree.Node(value=0), a=curve, min_samples_in_zone=12)
    print(root)

.. code-block::
    
        _________0_________
       /                   \
      0___              ____51
     /    \            /      \
    0     _19        _51       77
         /   \      /   \
        19    31   51    64

.. code-block:: python
    
    root = _anova_recursive_tree_build(node=binarytree.Node(value=0), a=curve, min_samples_in_zone=6)
    print(root)

.. code-block::

        _______________________________0_______________
       /                                               \
      0__                                     __________51_________
     /   \                                   /                     \
    0     7___                             _51___               ____77
         /    \                           /      \             /      \
        7     _14___                     51      _57         _77       90
             /      \                           /   \       /   \
            14      _20___                     57    66    77    83
                   /      \
                  20      _26___
                         /      \
                        26      _37
                               /   \
                              37    43


Finally, applying the anova zonation in a function

.. code-block:: python

    def anova_zoning(input_array:numpy.ndarray, min_samples_in_zone=2):
        """anova zoning
        """
        output_array = numpy.zeros_like(input_array)

        root = _anova_recursive_tree_build(node=binarytree.Node(value=0), a=input_array, min_samples_in_zone=min_samples_in_zone)
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

Using the same curve as the constant thickness log blocking example above, one can see that the blocked zones using the analysis of variance statistic are more closely representative of geologic packages.

.. image:: images/anova_min_samples_in_zone.png
    :alt: anova log blocking example

References
----------
- Al-Adani, Nabil, 2012, Data Blockign or Zoning: Well-Log-Data Application: Journal of Canadian Petroleum Technology.


