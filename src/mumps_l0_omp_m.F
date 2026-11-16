C
C  This file is part of MUMPS 5.9.2, released
C  on Fri Nov 13 10:00:55 UTC 2026
C
C
C  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
C  =============================================================================
C  Modification for thread safety (Phase 0), 2026-11-17
C  Author: Dorival Pedroso (Kleo completed the edit)
C  Nature: removed the !$OMP THREADPRIVATE directive and the dead module
C          globals; the module body is now empty and thread ids come from
C          OpenMP runtime queries into PRIVATE locals. No public interface
C          changed.
C  =============================================================================
C
      MODULE MUMPS_L0_OMP_M
      END MODULE MUMPS_L0_OMP_M
