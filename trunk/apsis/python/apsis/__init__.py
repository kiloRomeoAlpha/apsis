#!/usr/bin/env python

# $Id: __init__.py,v 1.7 2002/08/06 20:26:29 anderson Exp $
# ---------------------------------------------------------------------
# init file for Pipeline package of the ACS Pipeline

__version__      = '$Revision: 1.7 $ '[11:-3]
__version_date__ = '$Date: 2002/08/06 20:26:29 $ '[7:-3]


__all__ = ['ingest','align','combDither','pyblot','combFilter','detectionCatalog','filterCatalog','colorCatalog']

try:
    print "Package Apsis, release "+__release__ + ", loading..."
except NameError:
    print "Loading a non-release of the Apsis package."

