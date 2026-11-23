// Build stamp: ALWAYS recompiled, never cached.
//
// Why this file exists: the stamp used to live in robot_safety.cpp as
// __DATE__ " " __TIME__. With the incremental build (mkbuild.py) that file
// is skipped when untouched, so the DLL kept reporting the timestamp of the
// last robot_safety.cpp compile while spray_traj.cpp had in fact changed --
// i.e. the UI showed "16:12" for a 16:40 DLL, which is exactly the confusion
// this stamp was added to prevent.
//
// Keeping it in a 3-line TU costs ~0.3s to compile and makes the value always
// equal to the current link.
extern "C" const char* spraystamp_build_time(void) {
    return __DATE__ " " __TIME__;
}
