
// .\Release\benchmarks.exe --benchmark_repetitions=5 --benchmark_min_time=30 --benchmark_filter=HarmonicOscillator                                                                                                                 // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/03/09-23:14:45
// Benchmark                                                                                           Time(ns)    CPU(ns) Iterations
// ----------------------------------------------------------------------------------------------------------------------------------
// BM_SolveHarmonicOscillator<&Integrator::Leapfrog>                                                 2485034999 2485615933          3                               +4.15606773735669070e-07 m, +4.16264394754398830e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Leapfrog>                                                 2451504525 2444015667          3                               +4.15606773735669070e-07 m, +4.16264394754398830e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Leapfrog>                                                 2446403957 2438815633          3                               +4.15606773735669070e-07 m, +4.16264394754398830e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Leapfrog>                                                 2449195521 2438815633          3                               +4.15606773735669070e-07 m, +4.16264394754398830e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Leapfrog>                                                 2451225496 2444015667          3                               +4.15606773735669070e-07 m, +4.16264394754398830e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Leapfrog>_mean                                            2456672900 2450255707          3                               +4.15606773735669070e-07 m, +4.16264394754398830e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Leapfrog>_stddev                                            14297953   17832400          0                               +4.15606773735669070e-07 m, +4.16264394754398830e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::PseudoLeapfrog>                                           2383059002 2381615267          3                               +4.15606773750305800e-07 m, +4.16261893789685730e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::PseudoLeapfrog>                                           2387063480 2376415233          3                               +4.15606773750305800e-07 m, +4.16261893789685730e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::PseudoLeapfrog>                                           2388716039 2392015333          3                               +4.15606773750305800e-07 m, +4.16261893789685730e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::PseudoLeapfrog>                                           2382651862 2386815300          3                               +4.15606773750305800e-07 m, +4.16261893789685730e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::PseudoLeapfrog>                                           2381467637 2371215200          3                               +4.15606773750305800e-07 m, +4.16261893789685730e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::PseudoLeapfrog>_mean                                      2384591604 2381615267          3                               +4.15606773750305800e-07 m, +4.16261893789685730e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::PseudoLeapfrog>_stddev                                       2792553    7353958          0                               +4.15606773750305800e-07 m, +4.16261893789685730e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order2Optimal>                          2450820820 2433615600          3                               +2.01686230361618200e-07 m, +2.02003871729904490e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order2Optimal>                          2450978812 2444015667          3                               +2.01686230361618200e-07 m, +2.02003871729904490e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order2Optimal>                          2454569237 2449215700          3                               +2.01686230361618200e-07 m, +2.02003871729904490e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order2Optimal>                          2459043718 2444015667          3                               +2.01686230361618200e-07 m, +2.02003871729904490e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order2Optimal>                          2454426265 2449215700          3                               +2.01686230361618200e-07 m, +2.02003871729904490e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order2Optimal>_mean                     2453967771 2444015667          3                               +2.01686230361618200e-07 m, +2.02003871729904490e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order2Optimal>_stddev                      3005807    5696351          0                               +2.01686230361618200e-07 m, +2.02003871729904490e-07 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Ruth1983>                                                 2600918270 2600016667          3                               +7.13595849077819370e-14 m, +8.22154844204447950e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Ruth1983>                                                 2589792308 2589616600          3                               +7.13595849077819370e-14 m, +8.22154844204447950e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Ruth1983>                                                 2591876245 2579216533          3                               +7.13595849077819370e-14 m, +8.22154844204447950e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Ruth1983>                                                 2597641683 2589616600          3                               +7.13595849077819370e-14 m, +8.22154844204447950e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Ruth1983>                                                 2592631794 2594816633          3                               +7.13595849077819370e-14 m, +8.22154844204447950e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Ruth1983>_mean                                            2594572060 2590656607          3                               +7.13595849077819370e-14 m, +8.22154844204447950e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Ruth1983>_stddev                                             4088399    6898624          0                               +7.13595849077819370e-14 m, +8.22154844204447950e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order3Optimal>                          2593483345 2589616600          3                               +6.86607888600310190e-14 m, +5.68867869477074350e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order3Optimal>                          2588724272 2584416567          3                               +6.86607888600310190e-14 m, +5.68867869477074350e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order3Optimal>                          2595668849 2589616600          3                               +6.86607888600310190e-14 m, +5.68867869477074350e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order3Optimal>                          2613592284 2600016667          3                               +6.86607888600310190e-14 m, +5.68867869477074350e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order3Optimal>                          2593292769 2579216533          3                               +6.86607888600310190e-14 m, +5.68867869477074350e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order3Optimal>_mean                     2596952304 2588576593          3                               +6.86607888600310190e-14 m, +5.68867869477074350e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order3Optimal>_stddev                      8621630    6898624          0                               +6.86607888600310190e-14 m, +5.68867869477074350e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta>          2777868168 2776817800          3                               +1.64770974642181050e-13 m, +1.64847302475124020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta>          2775879916 2761217700          3                               +1.64770974642181050e-13 m, +1.64847302475124020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta>          2776412384 2771617767          3                               +1.64770974642181050e-13 m, +1.64847302475124020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta>          2772056398 2776817800          3                               +1.64770974642181050e-13 m, +1.64847302475124020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta>          2779580649 2756017667          3                               +1.64770974642181050e-13 m, +1.64847302475124020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta>_mean     2776359503 2768497747          3                               +1.64770974642181050e-13 m, +1.64847302475124020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta>_stddev      2506048    8449054          0                               +1.64770974642181050e-13 m, +1.64847302475124020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions>        2673820711 2667617100          3                               +1.64717198214425760e-13 m, +1.64795260770844720e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions>        2672528800 2667617100          3                               +1.64717198214425760e-13 m, +1.64795260770844720e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions>        2672010001 2667617100          3                               +1.64717198214425760e-13 m, +1.64795260770844720e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions>        2670338449 2657217033          3                               +1.64717198214425760e-13 m, +1.64795260770844720e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions>        2673207498 2667617100          3                               +1.64717198214425760e-13 m, +1.64795260770844720e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions>_mean   2672381092 2665537087          3                               +1.64717198214425760e-13 m, +1.64795260770844720e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions>_stddev    1190577    4160027          0                               +1.64717198214425760e-13 m, +1.64795260770844720e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order4Optimal>                          2737638927 2750817633          3                               +8.21200746292660710e-14 m, +8.20767065423666510e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order4Optimal>                          2741549946 2719617433          3                               +8.21200746292660710e-14 m, +8.20767065423666510e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order4Optimal>                          2748253481 2745617600          3                               +8.21200746292660710e-14 m, +8.20767065423666510e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order4Optimal>                          2739052749 2735217533          3                               +8.21200746292660710e-14 m, +8.20767065423666510e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order4Optimal>                          2742142734 2709217367          3                               +8.21200746292660710e-14 m, +8.20767065423666510e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order4Optimal>_mean                     2741727567 2732097513          3                               +8.21200746292660710e-14 m, +8.20767065423666510e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order4Optimal>_stddev                      3651386   15634728          0                               +8.21200746292660710e-14 m, +8.20767065423666510e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order5Optimal>                          3027555681 3026419400          2                               +1.37018868473504090e-13 m, +1.37057032389975570e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order5Optimal>                          3070461250 3026419400          2                               +1.37018868473504090e-13 m, +1.37057032389975570e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order5Optimal>                          3030260063 3026419400          2                               +1.37018868473504090e-13 m, +1.37057032389975570e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order5Optimal>                          3032249189 3026419400          2                               +1.37018868473504090e-13 m, +1.37057032389975570e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order5Optimal>                          3026327190 3000419233          3                               +1.37018868473504090e-13 m, +1.37057032389975570e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order5Optimal>_mean                     3036366721 3019328445          2                               +1.37018868473504090e-13 m, +1.37057032389975570e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::McLachlanAtela1992Order5Optimal>_stddev                     16211696   11579480          0                               +1.37018868473504090e-13 m, +1.37057032389975570e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6A>                                       3227502800 3221420650          2                               +2.08746214758193100e-13 m, +2.09027239961301350e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6A>                                       3217163182 3205820550          2                               +2.08746214758193100e-13 m, +2.09027239961301350e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6A>                                       3217649675 3182420400          2                               +2.08746214758193100e-13 m, +2.09027239961301350e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6A>                                       3219068623 3205820550          2                               +2.08746214758193100e-13 m, +2.09027239961301350e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6A>                                       3222009516 3229220700          2                               +2.08746214758193100e-13 m, +2.09027239961301350e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6A>_mean                                  3220678759 3208940570          2                               +2.08746214758193100e-13 m, +2.09027239961301350e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6A>_stddev                                   3806983   16061286          0                               +2.08746214758193100e-13 m, +2.09027239961301350e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6B>                                       3231235243 3221420650          2                               +3.59743485001118300e-13 m, +3.60239615915247670e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6B>                                       3228966451 3190220450          2                               +3.59743485001118300e-13 m, +3.60239615915247670e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6B>                                       3221972442 3221420650          2                               +3.59743485001118300e-13 m, +3.60239615915247670e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6B>                                       3225649572 3221420650          2                               +3.59743485001118300e-13 m, +3.60239615915247670e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6B>                                       3236787654 3244820800          2                               +3.59743485001118300e-13 m, +3.60239615915247670e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6B>_mean                                  3228922272 3219860640          2                               +3.59743485001118300e-13 m, +3.60239615915247670e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6B>_stddev                                   5024474   17371536          0                               +3.59743485001118300e-13 m, +3.60239615915247670e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6C>                                       3222663975 3182420400          2                               +1.66231611808953520e-13 m, +1.66453656413878550e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6C>                                       3222845650 3237020750          2                               +1.66231611808953520e-13 m, +1.66453656413878550e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6C>                                       3221522666 3221420650          2                               +1.66231611808953520e-13 m, +1.66453656413878550e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6C>                                       3219599581 3213620600          2                               +1.66231611808953520e-13 m, +1.66453656413878550e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6C>                                       3222406722 3205820550          2                               +1.66231611808953520e-13 m, +1.66453656413878550e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6C>_mean                                  3221807719 3212060590          2                               +1.66231611808953520e-13 m, +1.66453656413878550e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order6C>_stddev                                   1193852   18058421          0                               +1.66231611808953520e-13 m, +1.66453656413878550e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8A>                                       4352442480 4344627850          2                               +5.11861386609524520e-13 m, +5.12544867659059380e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8A>                                       4332950092 4321227700          2                               +5.11861386609524520e-13 m, +5.12544867659059380e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8A>                                       4339432932 4321227700          2                               +5.11861386609524520e-13 m, +5.12544867659059380e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8A>                                       4336501933 4313427650          2                               +5.11861386609524520e-13 m, +5.12544867659059380e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8A>                                       4349410511 4336827800          2                               +5.11861386609524520e-13 m, +5.12544867659059380e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8A>_mean                                  4342147589 4327467740          2                               +5.11861386609524520e-13 m, +5.12544867659059380e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8A>_stddev                                   7517600   11463685          0                               +5.11861386609524520e-13 m, +5.12544867659059380e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8B>                                       4340371943 4329027750          2                               +1.74530528918026560e-13 m, +1.74762981863807450e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8B>                                       4344166136 4336827800          2                               +1.74530528918026560e-13 m, +1.74762981863807450e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8B>                                       4337970234 4329027750          2                               +1.74530528918026560e-13 m, +1.74762981863807450e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8B>                                       4336148716 4336827800          2                               +1.74530528918026560e-13 m, +1.74762981863807450e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8B>                                       4350216127 4336827800          2                               +1.74530528918026560e-13 m, +1.74762981863807450e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8B>_mean                                  4341774631 4333707780          2                               +1.74530528918026560e-13 m, +1.74762981863807450e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8B>_stddev                                   5001401    3821228          0                               +1.74530528918026560e-13 m, +1.74762981863807450e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8C>                                       4342942572 4336827800          2                               +4.31519403543134670e-13 m, +4.32123087312774600e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8C>                                       4347059108 4313427650          2                               +4.31519403543134670e-13 m, +4.32123087312774600e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8C>                                       4339219785 4297827550          2                               +4.31519403543134670e-13 m, +4.32123087312774600e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8C>                                       4330824471 4321227700          2                               +4.31519403543134670e-13 m, +4.32123087312774600e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8C>                                       4335027076 4321227700          2                               +4.31519403543134670e-13 m, +4.32123087312774600e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8C>_mean                                  4339014602 4318107680          2                               +4.31519403543134670e-13 m, +4.32123087312774600e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8C>_stddev                                   5712508   12673581          0                               +4.31519403543134670e-13 m, +4.32123087312774600e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8D>                                       4321775056 4305627600          2                               +2.44630704582249340e-13 m, +2.44953363148781020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8D>                                       4345233417 4336827800          2                               +2.44630704582249340e-13 m, +2.44953363148781020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8D>                                       4326605654 4313427650          2                               +2.44630704582249340e-13 m, +2.44953363148781020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8D>                                       4327863074 4290027500          2                               +2.44630704582249340e-13 m, +2.44953363148781020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8D>                                       4327442981 4313427650          2                               +2.44630704582249340e-13 m, +2.44953363148781020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8D>_mean                                  4329784036 4311867640          2                               +2.44630704582249340e-13 m, +2.44953363148781020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8D>_stddev                                   8026206   15124858          0                               +2.44630704582249340e-13 m, +2.44953363148781020e-13 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8E>                                       4313737608 4313427650          2                               +8.20819107127945810e-14 m, +8.20281342850393000e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8E>                                       4340174771 4305627600          2                               +8.20819107127945810e-14 m, +8.20281342850393000e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8E>                                       4343757368 4329027750          2                               +8.20819107127945810e-14 m, +8.20281342850393000e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8E>                                       4334179378 4329027750          2                               +8.20819107127945810e-14 m, +8.20281342850393000e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8E>                                       4322726227 4290027500          2                               +8.20819107127945810e-14 m, +8.20281342850393000e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8E>_mean                                  4330915070 4313427650          2                               +8.20819107127945810e-14 m, +8.20281342850393000e-14 m kg s^-1  // NOLINT(whitespace/line_length)
// BM_SolveHarmonicOscillator<&Integrator::Yoshida1990Order8E>_stddev                                  11166421   14799554          0                               +8.20819107127945810e-14 m, +8.20281342850393000e-14 m kg s^-1  // NOLINT(whitespace/line_length)

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include <algorithm>
#include <vector>

#include "base/not_null.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "testing_utilities/numerical_analysis.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using integrators::SPRKIntegrator;
using quantities::Abs;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Momentum;
using quantities::SIUnit;
using testing_utilities::ComputeHarmonicOscillatorForce;
using testing_utilities::ComputeHarmonicOscillatorVelocity;

namespace benchmarks {

using Integrator = SPRKIntegrator<Length, Momentum>;

void SolveHarmonicOscillatorAndComputeError(
    not_null<benchmark::State*> const state,
    not_null<Length*> const q_error,
    not_null<Momentum*> const p_error,
    Integrator::Scheme const& (Integrator::*scheme)() const) {
  std::vector<Integrator::SystemState> solution;
  Integrator integrator;
  Integrator::Parameters parameters;

  integrator.Initialize((integrator.*scheme)());

  parameters.initial.positions.emplace_back(SIUnit<Length>());
  parameters.initial.momenta.emplace_back(Momentum());
  parameters.initial.time = Time();
#ifdef _DEBUG
  parameters.tmax = 100.0 * SIUnit<Time>();
#else
  parameters.tmax = 1000.0 * SIUnit<Time>();
#endif
  parameters.Δt = 1.0E-4 * SIUnit<Time>();
  parameters.sampling_period = 1;
  integrator.Solve(&ComputeHarmonicOscillatorForce,
                   &ComputeHarmonicOscillatorVelocity,
                   parameters,
                   &solution);

  state->PauseTiming();
  *q_error = Length();
  *p_error = Momentum();
  for (std::size_t i = 0; i < solution.size(); ++i) {
    *q_error = std::max(*q_error,
                        Abs(solution[i].positions[0].value -
                            SIUnit<Length>() *
                            Cos(solution[i].time.value *
                                SIUnit<AngularFrequency>())));
    *p_error = std::max(*p_error,
                        Abs(solution[i].momenta[0].value +
                            SIUnit<Momentum>() *
                            Sin(solution[i].time.value *
                                SIUnit<AngularFrequency>())));
  }
  state->ResumeTiming();
}

template<Integrator::Scheme const& (Integrator::*scheme)() const>
void BM_SolveHarmonicOscillator(
    benchmark::State& state) {  // NOLINT(runtime/references)
  Length   q_error;
  Momentum p_error;
  while (state.KeepRunning()) {
    SolveHarmonicOscillatorAndComputeError(&state, &q_error, &p_error, scheme);
  }
  std::stringstream ss;
  ss << q_error << ", " << p_error;
  state.SetLabel(ss.str());
}

BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Leapfrog);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::PseudoLeapfrog);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order2Optimal);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Ruth1983);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order3Optimal);
BENCHMARK_TEMPLATE(
    BM_SolveHarmonicOscillator,
    &Integrator::CandyRozmus1991ForestRuth1990SynchronousMomenta);
BENCHMARK_TEMPLATE(
    BM_SolveHarmonicOscillator,
    &Integrator::CandyRozmus1991ForestRuth1990SynchronousPositions);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order4Optimal);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator,
                   &Integrator::McLachlanAtela1992Order5Optimal);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order6A);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order6B);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order6C);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8A);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8B);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8C);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8D);
BENCHMARK_TEMPLATE(BM_SolveHarmonicOscillator, &Integrator::Yoshida1990Order8E);

}  // namespace benchmarks
}  // namespace principia
