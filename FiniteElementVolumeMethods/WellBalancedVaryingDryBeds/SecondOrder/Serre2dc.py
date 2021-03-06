# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Serre2dc', [dirname(__file__)])
        except ImportError:
            import _Serre2dc
            return _Serre2dc
        if fp is not None:
            try:
                _mod = imp.load_module('_Serre2dc', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Serre2dc = swig_import_helper()
    del swig_import_helper
else:
    import _Serre2dc
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0



def mallocPy(n):
    return _Serre2dc.mallocPy(n)
mallocPy = _Serre2dc.mallocPy

def conc(a, b, c, n, m, k, d):
    return _Serre2dc.conc(a, b, c, n, m, k, d)
conc = _Serre2dc.conc

def writetomem(x, i, f):
    return _Serre2dc.writetomem(x, i, f)
writetomem = _Serre2dc.writetomem

def readfrommem(x, i):
    return _Serre2dc.readfrommem(x, i)
readfrommem = _Serre2dc.readfrommem

def deallocPy(x):
    return _Serre2dc.deallocPy(x)
deallocPy = _Serre2dc.deallocPy

def getufromG(h, G, bed, hMbeg, hMend, GMbeg, GMend, uMbeg, uMend, wMbeg, wMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbbc, u, hhbc, Ghbc, whbc, bedhbc):
    return _Serre2dc.getufromG(h, G, bed, hMbeg, hMend, GMbeg, GMend, uMbeg, uMend, wMbeg, wMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbbc, u, hhbc, Ghbc, whbc, bedhbc)
getufromG = _Serre2dc.getufromG

def RKstep(a, b, n):
    return _Serre2dc.RKstep(a, b, n)
RKstep = _Serre2dc.RKstep

def uhall(x, h, u, n, nBC, dx):
    return _Serre2dc.uhall(x, h, u, n, nBC, dx)
uhall = _Serre2dc.uhall

def hall(x, h, n, nBC, dx):
    return _Serre2dc.hall(x, h, n, nBC, dx)
hall = _Serre2dc.hall

def HankEnergyall(x, h, u, g, n, nBC, dx):
    return _Serre2dc.HankEnergyall(x, h, u, g, n, nBC, dx)
HankEnergyall = _Serre2dc.HankEnergyall

def hALLW(hbc, n, dx):
    return _Serre2dc.hALLW(hbc, n, dx)
hALLW = _Serre2dc.hALLW

def uhALLW(hbc, ubc, n, dx):
    return _Serre2dc.uhALLW(hbc, ubc, n, dx)
uhALLW = _Serre2dc.uhALLW

def HamilW(hbc, ubc, n, dx):
    return _Serre2dc.HamilW(hbc, ubc, n, dx)
HamilW = _Serre2dc.HamilW

def evolvewrapForcing(G, h, bed, hMbeg, hMend, wMbeg, wMend, GMbeg, GMend, uMbeg, uMend, bMbeg, bMend, hMbeg1, hMend1, wMbeg1, wMend1, GMbeg1, GMend1, uMbeg1, uMend1, g, dx, dt, n, nGhBC, unBC, bnBC, nGhhbc, nubc, nbhc, theta, hhbc, whbc, Ghbc, bedhbc, ubc, x, t, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9):
    return _Serre2dc.evolvewrapForcing(G, h, bed, hMbeg, hMend, wMbeg, wMend, GMbeg, GMend, uMbeg, uMend, bMbeg, bMend, hMbeg1, hMend1, wMbeg1, wMend1, GMbeg1, GMend1, uMbeg1, uMend1, g, dx, dt, n, nGhBC, unBC, bnBC, nGhhbc, nubc, nbhc, theta, hhbc, whbc, Ghbc, bedhbc, ubc, x, t, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9)
evolvewrapForcing = _Serre2dc.evolvewrapForcing

def getufromGsplit(h, G, bed, hMbeg, hMend, GMbeg, GMend, uMbeg, uMend, wMbeg, wMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbbc, u, hhbc, Ghbc, whbc, bedhbc):
    return _Serre2dc.getufromGsplit(h, G, bed, hMbeg, hMend, GMbeg, GMend, uMbeg, uMend, wMbeg, wMend, bMbeg, bMend, theta, dx, n, m, nGhBC, unBC, nbBC, nGhbc, nubc, nbbc, u, hhbc, Ghbc, whbc, bedhbc)
getufromGsplit = _Serre2dc.getufromGsplit
# This file is compatible with both classic and new-style classes.


