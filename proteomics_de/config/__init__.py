"""proteomics_de.config — the experimental design and the shared constants.

Two modules live here:

* :mod:`proteomics_de.config.design` reads ``config/sample_sheet.tsv`` and is the
  single source of truth for *which* samples exist, which group each belongs to,
  and in which order the intensity columns must appear. Everything that today
  hardcodes ``Intensity 31578`` / ``ctrl_31578`` / ``["control","control",
  "treated","treated"]`` should ask this module instead.
* :mod:`proteomics_de.config.constants` re-exports the numeric thresholds, the
  RNG seed, and the figure caveat string that are currently duplicated across
  ``foldchange.py``, ``limma_test.py``, ``limma_test.R`` and ``viz/style.py``.

``config.yaml`` (also in this directory) stays what it has always been:
declarative documentation of the design. It is not parsed by these modules --
``sample_sheet.tsv`` is the machine-readable contract.

Note that ``proteomics_de/`` itself deliberately has no ``__init__.py``; it is an
implicit namespace package. Only the subpackages are regular packages.
"""
