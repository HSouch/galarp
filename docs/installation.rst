Installation
============


Stable Version
--------------

``GalaRP`` exists in a stable version on pip:

.. code-block:: python

    python -m pip install galarp


Development
-----------


Development is active on `Github <https://github.com/HSouch/galarp/>`_, and the Authors encourage
updates and suggestions to improve and enhance its functionality.

To install the current development version:

.. code-block:: python

    python -m pip install git+https://github.com/HSouch/galarp


Gala and GSL
------------

``gala`` requires the Gnu Scientific Library (GSL) to have full functionality. 
It is recommended to install gala first, and ensure GSL is running properly before trying to
install ``galarp``.

Please refer to  `Gala's Installation Page <https://gala.adrian.pw/en/latest/install.html>`_.

To check that gsl is working, run the following in your virtual environment:

.. code-block:: python

    import gala
    import gala.potential as gp

The second import statement directly calls GSL, so if it runs without any errors, everything
has been installed correctly on the ``gala`` side.


Running Tests
-------------

If contributing code via GitHub, please ensure the code is testing compliant by adding tests to the proper locations,
following the same testing paradigm. Furthermore, check that everything is running properly with ``nox``.

>>> pip install nox
>>> nox

This will call a standard build test, a runthrough of ``pytest``, and then linting with ``ruff``. Make sure things 
complete before contributing code.