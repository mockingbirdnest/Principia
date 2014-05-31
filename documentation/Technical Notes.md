```cpp
template<typename AutonomousRightHandSideComputation,
         typename RightHandSideComputation>
class SymplecticIntegrator {
```

This change requires to have templates all the way down, and doesn't result in a significant performance improvement:

Before:
```
.\Release\benchmarks_tests.exe --benchmark_repetitions=5 --benchmark_min_time=300          .
Benchmarking on 1 X 3310 MHz CPU
2014/05/30-20:51:41
Benchmark                           Time(ns)    CPU(ns) Iterations
------------------------------------------------------------------
BM_SolveHarmonicOscillator        1388241978 1227819635         51                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1220045434 1215559792         50                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1214497281 1212439772         50                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1226465770 1223047840         50                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1231751867 1225231854         50                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator_mean   1256726528 1220847667        251                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator_stddev   66665752    5858502        251                                 1.37019e-013, 1.37057e-013

```
After:
```
.\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30
Benchmarking on 1 X 3310 MHz CPU
2014/05/31-23:15:17
Benchmark                           Time(ns)    CPU(ns) Iterations
------------------------------------------------------------------
BM_SolveHarmonicOscillator        1198510902 1193407650          6                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1200832249 1201207700          5                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1203193594 1198087680          5                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1198585925 1196007667          6                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator        1192404764 1188207617          6                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator_mean   1198469241 1195079089         28                                 1.37019e-013, 1.37057e-013
BM_SolveHarmonicOscillator_stddev    3587060    4384616         28                                 1.37019e-013, 1.37057e-013

```
This is for a harmonic oscillator where the force/velocity computation is trivial.  For a real problem (e.g., gravitation) the impact would be even smaller.
