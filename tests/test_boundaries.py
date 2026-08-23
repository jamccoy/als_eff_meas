"""Architectural invariants.

Cheap, and they catch the drift review does not: a dependency leaking into a
part that promised not to have one.
"""

from endstation.testing.boundary import assert_qt_free, modules_in

PACKAGE = modules_in("als_eff_meas")


def test_the_package_needs_no_gui_toolkit():
    """This is a library. It has no window, and importing it must not act as
    though it might -- a reduction run on a headless beamline machine or in CI
    should never touch a toolkit.

    A subprocess import rather than a source scan, so a *transitive* import is
    caught too.
    """
    assert_qt_free(PACKAGE, reason="als_eff_meas is a library with no window.")


def test_the_walk_finds_something():
    """Guard against discovery silently emptying and passing vacuously."""
    assert "als_eff_meas.geometry" in PACKAGE, PACKAGE


def test_importing_the_package_does_not_read_the_disk():
    """The bug this package was ported to fix.

    In the legacy scripts, importing the configuration *ran* ``dark_val()`` and
    ``gather_file()`` -- so reading the configuration touched the beamline tree,
    and the module could not be imported at all without the private data
    mounted. `config.py` is now inert on import, and this keeps it that way.

    Checked by importing with ``ALS_DATA_ROOT`` pointed somewhere that does not
    exist: if anything reached for the tree at import time it would raise.
    """
    import subprocess
    import sys

    probe = subprocess.run(
        [sys.executable, "-c", "import als_eff_meas.config"],
        capture_output=True,
        text=True,
        check=False,
        env={"PATH": "/usr/bin:/bin", "ALS_DATA_ROOT": "/nonexistent", "HOME": "/nonexistent"},
    )
    assert probe.returncode == 0, probe.stderr
