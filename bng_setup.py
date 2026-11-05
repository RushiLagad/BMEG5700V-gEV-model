"""
Point BioNetGen at its own Perl2 library directory before any model is built.

Why this file exists: on Windows the perl that ends up first on PATH is often the
MSYS / Git-for-Windows one. BNG2.pl locates its modules with FindBin, which cannot
resolve a C:\\ script path under that perl, so it looks for BNGUtils.pm in the
current working directory and dies with "Can't locate BNGUtils.pm in @INC".
BNG2.pl checks the BNGPATH environment variable before falling back to FindBin,
so setting it fixes the lookup.

Import this module before building a model. It does nothing if BNGPATH is already
set, and nothing harmful on Linux or macOS.
"""
import os
import pathlib
import sys


def configure_bngpath(verbose=False):
    """Set BNGPATH to the bundled BioNetGen directory. Returns the path, or None."""
    if os.environ.get("BNGPATH"):
        return os.environ["BNGPATH"]
    try:
        import bionetgen
    except ImportError:
        return None
    root = pathlib.Path(bionetgen.__file__).parent
    order = {"win32": ["bng-win"], "darwin": ["bng-mac"]}.get(sys.platform, ["bng-linux"])
    order += [d for d in ("bng-win", "bng-mac", "bng-linux") if d not in order]
    for name in order:
        d = root / name
        if (d / "Perl2" / "BNGUtils.pm").exists():
            # forward slashes on purpose: the MSYS perl does not take backslashes in @INC
            os.environ["BNGPATH"] = d.as_posix()
            if verbose:
                print(f"BNGPATH = {os.environ['BNGPATH']}")
            return os.environ["BNGPATH"]
    return None


configure_bngpath()
