Installation
------------

The installation requires a few simple steps and works on most Linux or Mac systems.

Download and extract the source file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can download the code directly from
`github <https://github.com/felixplasser/matilda>`_.

Use the latest release or the current development version.

Setup the path specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To use MATILDA the ``PATH`` and ``PYTHONPATH`` variables have to be adjusted.
For this purpose you can use the provided bash script `setpaths.bash`

.. code:: bash

    #!/bin/bash
    MATILDA=/PhotoChem/programs/MATILDA/GIT

    export PATH=$MATILDA/scripts:$PATH
    export PYTHONPATH=$MATILDA:$PYTHONPATH

Simply change the ``MATILDA`` variable to point to your actual installation.

Copy the above lines into your .bashrc file or run:

::

    source setpaths.bash

If you are not using bash, please set the environment variables using your preferred shell.

External packages
~~~~~~~~~~~~~~~~~

The following external packages are used by MATILDA and require a separate installation:

- `python3-numpy <http://numpy.scipy.org/>`_ - for basic numerical manipulations
- `python3-openbabel <http://openbabel.org/wiki/Python>`_ - for extended file-parsing capabilities of molecular structure files
- `python3-matplotlib <http://matplotlib.sourceforge.net/>`_ *(optional)* - for plotting of graphs

The first three are usually readily available with the standard installation tools, e.g. ``apt-get``, ``yum`` etc.
Alternatively, they may be downloaded from the URLs specified.
If no integrated installation is performed, then it is necessary to add these libraries to the ``PYTHONPATH`` (see above).

Using anaconda
~~~~~~~~~~~~~~

A straightforward and universal way of installing most of the required packages is through the use of Anaconda.

First download the `anaconda distribution <https://www.anaconda.com/distribution/>`_ and do the installation. Then run the commands:

::

    conda install numpy matplotlib
    conda install -c openbabel openbabel

Testing
~~~~~~~

MATILDA does not yet contain automated test routines.
Sample jobs with reference results are given in the ``EXAMPLES`` directory.
