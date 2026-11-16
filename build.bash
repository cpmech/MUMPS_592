#!/bin/bash
#
# build.bash — build this MUMPS checkout and install it
#
# The sources are copied to ${BUILD_DIR}/MUMPS_<version> and compiled there, so the checkout
# stays clean and no object file, library or generated Makefile.inc lands next to the sources.
# The copy is removed after a successful install (KEEP_BUILD=1 keeps it for inspection).
#
# The build files in zscripts/ are the ones this workstation uses for russell_sparse: same
# ordering (METIS + PORD), same PLAT, so the libraries are lib*_cpmech.{a,so} and the four
# arithmetic variants (s, d, c, z) are built, static and shared. These are the libraries
# RUSSELL_LOCAL_PREFIX expects under $PREFIX/{lib,include}/mumps.
#
# Usage:
#
#     ./build.bash            # Intel oneAPI (ifx/icx) with MKL; the workstation default
#     ./build.bash 0          # gfortran/clang with the system BLAS and LAPACK
#
# Environment:
#
#     PREFIX       installation prefix (default /usr/local)
#     BUILD_DIR    directory that holds the build copy (default /tmp)
#     TARGETS      what to build, from d, dshared, z, zshared (default: all four); a probe that
#                  only calls dmumps_c needs "d dshared" (about 40 s, the shared libraries a C
#                  probe links) or "d" alone if the probe is linked with ifx against the archives
#     JOBS         make jobs (default 8)
#     KEEP_BUILD   1 keeps the build copy after installation
#
# Examples:
#
#     ./build.bash                                  # /usr/local, MKL (asks for sudo)
#     PREFIX=$HOME/mumps ./build.bash               # a private install, no sudo
#     BUILD_DIR=$HOME/tmp ./build.bash 0            # gfortran/OpenBLAS, build in $HOME/tmp
#     TARGETS="d dshared" ./build.bash               # the real arithmetics, enough for a probe
#
set -euo pipefail

USE_INTEL_MKL="${1:-1}"
case "$USE_INTEL_MKL" in
    0|1) ;;
    *)
        echo "usage: $0 [0|1]   (1 = Intel MKL, the default; 0 = system BLAS/LAPACK)" >&2
        exit 1
        ;;
esac

SRC="$(cd "$(dirname "$0")" && pwd)"
PREFIX="${PREFIX:-/usr/local}"
BUILD_DIR="${BUILD_DIR:-/tmp}"
TARGETS="${TARGETS:-d dshared z zshared}"
JOBS="${JOBS:-8}"
INCDIR="$PREFIX/include/mumps"
LIBDIR="$PREFIX/lib/mumps"

case "$JOBS" in
    ''|*[!0-9]*|0)
        echo "error: JOBS must be a positive integer (got '$JOBS')" >&2
        exit 1
        ;;
esac
for target in $TARGETS; do
    case "$target" in
        d|z|dshared|zshared) ;;
        *)
            echo "error: unknown target '$target' in TARGETS (use d, dshared, z or zshared)" >&2
            exit 1
            ;;
    esac
done

# Version from the file that the release ships; the same string is compared at run time by
# libmumps against the headers the client compiled with.
if [ ! -f "$SRC/VERSION" ]; then
    echo "error: $SRC/VERSION not found; run this script from the MUMPS checkout" >&2
    exit 1
fi
VERSION="$(awk 'NR == 1 {print $2}' "$SRC/VERSION")"
if [ -z "$VERSION" ]; then
    echo "error: cannot read the version from $SRC/VERSION" >&2
    exit 1
fi
WORK="$BUILD_DIR/MUMPS_$VERSION"

for entry in Makefile src include libseq PORD; do
    if [ ! -e "$SRC/$entry" ]; then
        echo "error: $SRC/$entry not found; run this script from the MUMPS checkout" >&2
        exit 1
    fi
done

# Root containers and writable prefixes do not need sudo; /usr/local does.
SUDO=""
if [ "$EUID" != 0 ] && [ ! -w "$(dirname "$PREFIX")" ] && [ ! -w "$PREFIX" ]; then
    SUDO="sudo"
fi

# Orderings: the build files link METIS by name
# (grep -c, not grep -q: with pipefail, -q closes the pipe and ldconfig's SIGPIPE fails the test)
METIS_FOUND="$(ldconfig -p 2>/dev/null | grep -c 'libmetis\.so' || true)"
if [ "$METIS_FOUND" -eq 0 ]; then
    echo "error: libmetis not found (the build links -lmetis)" >&2
    echo "       install it first: sudo pacman -S metis" >&2
    exit 1
fi

# Compiler and build file for the chosen BLAS/LAPACK backend
if [ "$USE_INTEL_MKL" = 1 ]; then
    SETVARS="/opt/intel/oneapi/setvars.sh"
    if [ ! -f "$SETVARS" ]; then
        echo "error: $SETVARS not found; install Intel oneAPI or run '$0 0'" >&2
        exit 1
    fi
    set +u  # setvars.sh reads variables before setting them; -u would abort it
    # shellcheck disable=SC1090
    source "$SETVARS" >/dev/null
    set -u
    for tool in ifx icx; do
        if ! command -v "$tool" >/dev/null; then
            echo "error: $tool not found after sourcing $SETVARS" >&2
            exit 1
        fi
    done
    BUILD_FILE="$SRC/zscripts/MakefileMKL.inc"
    BACKEND="Intel oneAPI (ifx/icx) with MKL"
else
    for tool in gfortran clang; do
        if ! command -v "$tool" >/dev/null; then
            echo "error: $tool not found; install the GNU Fortran and Clang toolchains" >&2
            exit 1
        fi
    done
    BUILD_FILE="$SRC/zscripts/Makefile.inc"
    BACKEND="gfortran/clang with the system BLAS and LAPACK"
fi

echo "... MUMPS $VERSION"
echo "... source:  $SRC"
echo "... build:   $WORK"
echo "... install: $PREFIX"
echo "... backend: $BACKEND"

# Fresh copy of the working tree: the checkout may be modified (this is a fork), so the build
# must see the working tree, not a commit. Build products are not copied.
if [ -d "$WORK" ]; then
    echo "... removing previous $WORK"
    rm -rf "$WORK"
fi
mkdir -p "$WORK"
tar -C "$SRC" --exclude=.git --exclude='*.o' --exclude='*.a' --exclude='*.so' --exclude='*.mod' -cf - . \
    | tar -C "$WORK" -xf -
cp "$BUILD_FILE" "$WORK/Makefile.inc"

$SUDO mkdir -p "$LIBDIR" "$INCDIR"

# One target per arithmetic and linkage: the top-level Makefile is not reentrant across
# targets (object files are not recompiled when only the output target changes), so each
# build starts from a clean tree. The targets are independent of each other, so TARGETS can
# take a subset; "d" builds the two real arithmetics, "z" the two complex ones.
cd "$WORK"
build_variant () {
    local target="$1" ext="$2"
    echo "... make -j$JOBS $target"
    make -j"$JOBS" "$target"
    $SUDO cp -av lib/lib*"$ext" "$LIBDIR/"
    # The headers go in while the tree is built: mumps_int_def.h is generated into include/
    # by every build and removed by make clean (src/Makefile), and every client header
    # includes it in turn.
    $SUDO cp -av include/*.h "$INCDIR/"
}

for target in $TARGETS; do
    case "$target" in
        d|z) ext=.a ;;
        *) ext=.so ;;
    esac
    build_variant "$target" "$ext"
    make clean
done

# Register the library directory with the dynamic linker when installing system-wide.
if [ "$PREFIX" = "/usr/local" ] && [ "${LDCONFIG:-1}" != 0 ]; then
    echo "$LIBDIR" | $SUDO tee /etc/ld.so.conf.d/mumps.conf >/dev/null
    $SUDO ldconfig 2> >(grep -v 'is not an ELF file\|is not a symbolic link' >&2 || true)
fi

if [ "${KEEP_BUILD:-0}" != 1 ]; then
    rm -rf "$WORK"
else
    echo "... build copy kept at $WORK"
fi

echo "... installed in $LIBDIR:"
$SUDO ls -1 "$LIBDIR" | sed 's/^/    /'
echo "... headers in $INCDIR"
