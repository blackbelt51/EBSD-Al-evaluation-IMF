# -*- mode: python ; coding: utf-8 -*-
# PyInstaller-Spezifikation: Windows-Programmordner "EBSD-Auswertung" (onedir, ohne Konsole)
#   python -m PyInstaller --noconfirm --clean EBSD-Auswertung.spec
import os
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, collect_submodules

hidden = collect_submodules("mtex", filter=lambda n: ".docpages" not in n)
hidden += ["matplotlib.backends.backend_qtagg", "matplotlib.backends.backend_agg", "matplotlib.backends.backend_pdf",
           "matplotlib.backends.backend_svg", "matplotlib.backends.backend_ps"]
datas = collect_data_files("mtex", excludes=["docpages/**", "**/__pycache__/**"])
datas += [("icon.png", "."), ("icon.ico", ".")]
binaries = collect_dynamic_libs("finufft")
import finufft as _f
_libs = os.path.join(os.path.dirname(os.path.dirname(_f.__file__)), "finufft.libs")
if os.path.isdir(_libs):
    binaries += [(os.path.join(_libs, f), "finufft.libs") for f in os.listdir(_libs)]

excludes = ["tkinter", "_tkinter", "IPython", "jupyter_client", "ipykernel", "notebook", "pytest", "sphinx",
            "pyvista", "vtk", "vtkmodules", "h5py", "PyQt5", "PyQt6", "PySide2", "pandas",
            "PySide6.QtNetwork", "PySide6.QtQml", "PySide6.QtQuick", "PySide6.QtQuickWidgets", "PySide6.QtWebEngineCore",
            "PySide6.QtWebEngineWidgets", "PySide6.QtMultimedia", "PySide6.QtPdf", "PySide6.QtSql", "PySide6.QtTest",
            "PySide6.QtOpenGL", "PySide6.QtOpenGLWidgets", "PySide6.QtDesigner", "PySide6.QtHelp", "PySide6.Qt3DCore"]

a = Analysis(
    ["ebsd_gui.py"],
    pathex=["."],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden,
    excludes=excludes,
    noarchive=False,
    module_collection_mode={"mtex": "pyz+py"},    # Quelltexte mitliefern: numba-Cache (cache=True) braucht sie
)
# nicht benötigte Qt-Teile entfernen: Software-OpenGL (≈ 20 MB) und Übersetzungen außer Deutsch
a.binaries = [b for b in a.binaries if not b[0].lower().endswith("opengl32sw.dll")]
a.datas = [d for d in a.datas if not ("translations" in d[0].replace("\\", "/") and d[0].endswith(".qm")
                                      and not d[0].endswith("_de.qm"))]
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="EBSD-Auswertung",
    icon="icon.ico",
    console=False,
    debug=False,
    strip=False,
    upx=False,
)
coll = COLLECT(exe, a.binaries, a.datas, name="EBSD-Auswertung", strip=False, upx=False)
