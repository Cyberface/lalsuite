/**
 *
 * \file
 * \brief Tabulated Quasi-Normal Mode Information for Ringdown
 *
 * This file contains listed final spin values, and corresponding Quasi-Normal-Mode(QNM) central frequencies and decay rates. The spins are on [-1,1] where values are the dimensionless Kerr parameters S/M^2. The frequencies and decay rates are in units of 1/(s*pi*M), where M is 1. To compare values with tables given at, for example, www.phy.olemiss.edu/~berti/ringdown/, multiply by 2*pi, and note that the decay rate corresponds to the imaginary part of the QNM frequency.
 *
 * Values for spins between -0.994 and 0.994 were sourced from the above website (2014), while qualitatively accurate values for the remaining spins where sourced from the implementation used in arxiv:1404:3197. Both references use the analytic representation of QNMs developed by Leaver in 1986, and for spin values within +-0.994, have identical values within numerical error.
 *
 * */

#ifndef _LALSIM_IMR_PHENOMD_H
#define _LALSIM_IMR_PHENOMD_H

/* CONSTANTS */

/**
 * Dimensionless frequency (Mf) at which define the end of the waveform
 */
#define f_CUT 0.2

/**
 * Dimensionless frequency (Mf) at which the inspiral amplitude
 * switches to the intermediate amplitude
 */
#define AMP_fJoin_INS 0.014

/**
 * Dimensionless frequency (Mf) at which the inspiral phase
 * switches to the intermediate phase
 */
#define PHI_fJoin_INS 0.018

/**
  * Minimal final spin value below which the waveform might behave pathological
  * because the ISCO frequency is too low. For more details, see the review wiki
  * page https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/PhenD_LargeNegativeSpins
  */
#define MIN_FINAL_SPIN -0.717

/**
  * A large mass ratio causes memory over-runs.
  * We test and put the limit an order of magnitude above that of previous waveform models (which were around q=100).
  */
#define MAX_ALLOWED_MASS_RATIO 5000


/* QNM data */

static const int QNMData_length = 1003;
static const double QNMData_a[] = { -1.0, -0.999, -0.998,-0.996,         -0.994, -0.992, -0.99, -0.988, -0.986, -0.984, -0.982, \
-0.98, -0.978, -0.976, -0.974, -0.972, -0.97, -0.968, -0.966, -0.964, \
-0.962, -0.96, -0.958, -0.956, -0.954, -0.952, -0.95, -0.948, -0.946, \
-0.944, -0.942, -0.94, -0.938, -0.936, -0.934, -0.932, -0.93, -0.928, \
-0.926, -0.924, -0.922, -0.92, -0.918, -0.916, -0.914, -0.912, -0.91, \
-0.908, -0.906, -0.904, -0.902, -0.9, -0.898, -0.896, -0.894, -0.892, \
-0.89, -0.888, -0.886, -0.884, -0.882, -0.88, -0.878, -0.876, -0.874, \
-0.872, -0.87, -0.868, -0.866, -0.864, -0.862, -0.86, -0.858, -0.856, \
-0.854, -0.852, -0.85, -0.848, -0.846, -0.844, -0.842, -0.84, -0.838, \
-0.836, -0.834, -0.832, -0.83, -0.828, -0.826, -0.824, -0.822, -0.82, \
-0.818, -0.816, -0.814, -0.812, -0.81, -0.808, -0.806, -0.804, \
-0.802, -0.8, -0.798, -0.796, -0.794, -0.792, -0.79, -0.788, -0.786, \
-0.784, -0.782, -0.78, -0.778, -0.776, -0.774, -0.772, -0.77, -0.768, \
-0.766, -0.764, -0.762, -0.76, -0.758, -0.756, -0.754, -0.752, -0.75, \
-0.748, -0.746, -0.744, -0.742, -0.74, -0.738, -0.736, -0.734, \
-0.732, -0.73, -0.728, -0.726, -0.724, -0.722, -0.72, -0.718, -0.716, \
-0.714, -0.712, -0.71, -0.708, -0.706, -0.704, -0.702, -0.7, -0.698, \
-0.696, -0.694, -0.692, -0.69, -0.688, -0.686, -0.684, -0.682, -0.68, \
-0.678, -0.676, -0.674, -0.672, -0.67, -0.668, -0.666, -0.664, \
-0.662, -0.66, -0.658, -0.656, -0.654, -0.652, -0.65, -0.648, -0.646, \
-0.644, -0.642, -0.64, -0.638, -0.636, -0.634, -0.632, -0.63, -0.628, \
-0.626, -0.624, -0.622, -0.62, -0.618, -0.616, -0.614, -0.612, -0.61, \
-0.608, -0.606, -0.604, -0.602, -0.6, -0.598, -0.596, -0.594, -0.592, \
-0.59, -0.588, -0.586, -0.584, -0.582, -0.58, -0.578, -0.576, -0.574, \
-0.572, -0.57, -0.568, -0.566, -0.564, -0.562, -0.56, -0.558, -0.556, \
-0.554, -0.552, -0.55, -0.548, -0.546, -0.544, -0.542, -0.54, -0.538, \
-0.536, -0.534, -0.532, -0.53, -0.528, -0.526, -0.524, -0.522, -0.52, \
-0.518, -0.516, -0.514, -0.512, -0.51, -0.508, -0.506, -0.504, \
-0.502, -0.5, -0.498, -0.496, -0.494, -0.492, -0.49, -0.488, -0.486, \
-0.484, -0.482, -0.48, -0.478, -0.476, -0.474, -0.472, -0.47, -0.468, \
-0.466, -0.464, -0.462, -0.46, -0.458, -0.456, -0.454, -0.452, -0.45, \
-0.448, -0.446, -0.444, -0.442, -0.44, -0.438, -0.436, -0.434, \
-0.432, -0.43, -0.428, -0.426, -0.424, -0.422, -0.42, -0.418, -0.416, \
-0.414, -0.412, -0.41, -0.408, -0.406, -0.404, -0.402, -0.4, -0.398, \
-0.396, -0.394, -0.392, -0.39, -0.388, -0.386, -0.384, -0.382, -0.38, \
-0.378, -0.376, -0.374, -0.372, -0.37, -0.368, -0.366, -0.364, \
-0.362, -0.36, -0.358, -0.356, -0.354, -0.352, -0.35, -0.348, -0.346, \
-0.344, -0.342, -0.34, -0.338, -0.336, -0.334, -0.332, -0.33, -0.328, \
-0.326, -0.324, -0.322, -0.32, -0.318, -0.316, -0.314, -0.312, -0.31, \
-0.308, -0.306, -0.304, -0.302, -0.3, -0.298, -0.296, -0.294, -0.292, \
-0.29, -0.288, -0.286, -0.284, -0.282, -0.28, -0.278, -0.276, -0.274, \
-0.272, -0.27, -0.268, -0.266, -0.264, -0.262, -0.26, -0.258, -0.256, \
-0.254, -0.252, -0.25, -0.248, -0.246, -0.244, -0.242, -0.24, -0.238, \
-0.236, -0.234, -0.232, -0.23, -0.228, -0.226, -0.224, -0.222, -0.22, \
-0.218, -0.216, -0.214, -0.212, -0.21, -0.208, -0.206, -0.204, \
-0.202, -0.2, -0.198, -0.196, -0.194, -0.192, -0.19, -0.188, -0.186, \
-0.184, -0.182, -0.18, -0.178, -0.176, -0.174, -0.172, -0.17, -0.168, \
-0.166, -0.164, -0.162, -0.16, -0.158, -0.156, -0.154, -0.152, -0.15, \
-0.148, -0.146, -0.144, -0.142, -0.14, -0.138, -0.136, -0.134, \
-0.132, -0.13, -0.128, -0.126, -0.124, -0.122, -0.12, -0.118, -0.116, \
-0.114, -0.112, -0.11, -0.108, -0.106, -0.104, -0.102, -0.1, -0.098, \
-0.096, -0.094, -0.092, -0.09, -0.088, -0.086, -0.084, -0.082, -0.08, \
-0.078, -0.076, -0.074, -0.072, -0.07, -0.068, -0.066, -0.064, \
-0.062, -0.06, -0.058, -0.056, -0.054, -0.052, -0.05, -0.048, -0.046, \
-0.044, -0.042, -0.04, -0.038, -0.036, -0.034, -0.032, -0.03, -0.028, \
-0.026, -0.024, -0.022, -0.02, -0.018, -0.016, -0.014, -0.012, -0.01, \
-0.008, -0.006, -0.004, -0.002, 0., 0.002, 0.004, 0.006, 0.008, 0.01, \
0.012, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 0.028, 0.03, \
0.032, 0.034, 0.036, 0.038, 0.04, 0.042, 0.044, 0.046, 0.048, 0.05, \
0.052, 0.054, 0.056, 0.058, 0.06, 0.062, 0.064, 0.066, 0.068, 0.07, \
0.072, 0.074, 0.076, 0.078, 0.08, 0.082, 0.084, 0.086, 0.088, 0.09, \
0.092, 0.094, 0.096, 0.098, 0.1, 0.102, 0.104, 0.106, 0.108, 0.11, \
0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, \
0.132, 0.134, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, \
0.152, 0.154, 0.156, 0.158, 0.16, 0.162, 0.164, 0.166, 0.168, 0.17, \
0.172, 0.174, 0.176, 0.178, 0.18, 0.182, 0.184, 0.186, 0.188, 0.19, \
0.192, 0.194, 0.196, 0.198, 0.2, 0.202, 0.204, 0.206, 0.208, 0.21, \
0.212, 0.214, 0.216, 0.218, 0.22, 0.222, 0.224, 0.226, 0.228, 0.23, \
0.232, 0.234, 0.236, 0.238, 0.24, 0.242, 0.244, 0.246, 0.248, 0.25, \
0.252, 0.254, 0.256, 0.258, 0.26, 0.262, 0.264, 0.266, 0.268, 0.27, \
0.272, 0.274, 0.276, 0.278, 0.28, 0.282, 0.284, 0.286, 0.288, 0.29, \
0.292, 0.294, 0.296, 0.298, 0.3, 0.302, 0.304, 0.306, 0.308, 0.31, \
0.312, 0.314, 0.316, 0.318, 0.32, 0.322, 0.324, 0.326, 0.328, 0.33, \
0.332, 0.334, 0.336, 0.338, 0.34, 0.342, 0.344, 0.346, 0.348, 0.35, \
0.352, 0.354, 0.356, 0.358, 0.36, 0.362, 0.364, 0.366, 0.368, 0.37, \
0.372, 0.374, 0.376, 0.378, 0.38, 0.382, 0.384, 0.386, 0.388, 0.39, \
0.392, 0.394, 0.396, 0.398, 0.4, 0.402, 0.404, 0.406, 0.408, 0.41, \
0.412, 0.414, 0.416, 0.418, 0.42, 0.422, 0.424, 0.426, 0.428, 0.43, \
0.432, 0.434, 0.436, 0.438, 0.44, 0.442, 0.444, 0.446, 0.448, 0.45, \
0.452, 0.454, 0.456, 0.458, 0.46, 0.462, 0.464, 0.466, 0.468, 0.47, \
0.472, 0.474, 0.476, 0.478, 0.48, 0.482, 0.484, 0.486, 0.488, 0.49, \
0.492, 0.494, 0.496, 0.498, 0.5, 0.502, 0.504, 0.506, 0.508, 0.51, \
0.512, 0.514, 0.516, 0.518, 0.52, 0.522, 0.524, 0.526, 0.528, 0.53, \
0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 0.544, 0.546, 0.548, 0.55, \
0.552, 0.554, 0.556, 0.558, 0.56, 0.562, 0.564, 0.566, 0.568, 0.57, \
0.572, 0.574, 0.576, 0.578, 0.58, 0.582, 0.584, 0.586, 0.588, 0.59, \
0.592, 0.594, 0.596, 0.598, 0.6, 0.602, 0.604, 0.606, 0.608, 0.61, \
0.612, 0.614, 0.616, 0.618, 0.62, 0.622, 0.624, 0.626, 0.628, 0.63, \
0.632, 0.634, 0.636, 0.638, 0.64, 0.642, 0.644, 0.646, 0.648, 0.65, \
0.652, 0.654, 0.656, 0.658, 0.66, 0.662, 0.664, 0.666, 0.668, 0.67, \
0.672, 0.674, 0.676, 0.678, 0.68, 0.682, 0.684, 0.686, 0.688, 0.69, \
0.692, 0.694, 0.696, 0.698, 0.7, 0.702, 0.704, 0.706, 0.708, 0.71, \
0.712, 0.714, 0.716, 0.718, 0.72, 0.722, 0.724, 0.726, 0.728, 0.73, \
0.732, 0.734, 0.736, 0.738, 0.74, 0.742, 0.744, 0.746, 0.748, 0.75, \
0.752, 0.754, 0.756, 0.758, 0.76, 0.762, 0.764, 0.766, 0.768, 0.77, \
0.772, 0.774, 0.776, 0.778, 0.78, 0.782, 0.784, 0.786, 0.788, 0.79, \
0.792, 0.794, 0.796, 0.798, 0.8, 0.802, 0.804, 0.806, 0.808, 0.81, \
0.812, 0.814, 0.816, 0.818, 0.82, 0.822, 0.824, 0.826, 0.828, 0.83, \
0.832, 0.834, 0.836, 0.838, 0.84, 0.842, 0.844, 0.846, 0.848, 0.85, \
0.852, 0.854, 0.856, 0.858, 0.86, 0.862, 0.864, 0.866, 0.868, 0.87, \
0.872, 0.874, 0.876, 0.878, 0.88, 0.882, 0.884, 0.886, 0.888, 0.89, \
0.892, 0.894, 0.896, 0.898, 0.9, 0.902, 0.904, 0.906, 0.908, 0.91, \
0.912, 0.914, 0.916, 0.918, 0.92, 0.922, 0.924, 0.926, 0.928, 0.93, \
0.932, 0.934, 0.936, 0.938, 0.94, 0.942, 0.944, 0.946, 0.948, 0.95, \
0.952, 0.954, 0.956, 0.958, 0.96, 0.962, 0.964, 0.966, 0.968, 0.97, \
0.972, 0.974, 0.976, 0.978, 0.98, 0.982, 0.984, 0.986, 0.988, 0.99, \
0.992, 0.994,       0.996, 0.998, 0.999, 1.0};

static const double QNMData_fring[] = {0.0464014,0.0464110,0.0464197,0.0464373,       0.0464526, 0.046473, 0.0464909, 0.0465084, 0.0465259, \
0.0465435, 0.0465611, 0.0465789, 0.0465966, 0.0466144, 0.0466322, \
0.0466501, 0.0466679, 0.0466858, 0.0467038, 0.0467217, 0.0467397, \
0.0467577, 0.0467757, 0.0467937, 0.0468117, 0.0468298, 0.0468479, \
0.046866, 0.0468842, 0.0469024, 0.0469205, 0.0469388, 0.046957, \
0.0469752, 0.0469935, 0.0470118, 0.0470302, 0.0470485, 0.0470669, \
0.0470853, 0.0471037, 0.0471221, 0.0471406, 0.0471591, 0.0471776, \
0.0471962, 0.0472147, 0.0472333, 0.0472519, 0.0472705, 0.0472892, \
0.0473079, 0.0473266, 0.0473453, 0.0473641, 0.0473829, 0.0474017, \
0.0474205, 0.0474393, 0.0474582, 0.0474771, 0.047496, 0.047515, \
0.047534, 0.047553, 0.047572, 0.047591, 0.0476101, 0.0476292, \
0.0476483, 0.0476675, 0.0476866, 0.0477058, 0.0477251, 0.0477443, \
0.0477636, 0.0477829, 0.0478022, 0.0478216, 0.0478409, 0.0478603, \
0.0478798, 0.0478992, 0.0479187, 0.0479382, 0.0479577, 0.0479773, \
0.0479969, 0.0480165, 0.0480361, 0.0480558, 0.0480755, 0.0480952, \
0.0481149, 0.0481347, 0.0481545, 0.0481743, 0.0481942, 0.048214, \
0.0482339, 0.0482539, 0.0482738, 0.0482938, 0.0483138, 0.0483339, \
0.0483539, 0.048374, 0.0483941, 0.0484143, 0.0484345, 0.0484547, \
0.0484749, 0.0484951, 0.0485154, 0.0485357, 0.0485561, 0.0485764, \
0.0485968, 0.0486173, 0.0486377, 0.0486582, 0.0486787, 0.0486992, \
0.0487198, 0.0487404, 0.048761, 0.0487817, 0.0488024, 0.0488231, \
0.0488438, 0.0488646, 0.0488854, 0.0489062, 0.048927, 0.0489479, \
0.0489688, 0.0489898, 0.0490108, 0.0490318, 0.0490528, 0.0490739, \
0.0490949, 0.0491161, 0.0491372, 0.0491584, 0.0491796, 0.0492008, \
0.0492221, 0.0492434, 0.0492647, 0.0492861, 0.0493075, 0.0493289, \
0.0493504, 0.0493719, 0.0493934, 0.0494149, 0.0494365, 0.0494581, \
0.0494797, 0.0495014, 0.0495231, 0.0495448, 0.0495666, 0.0495884, \
0.0496102, 0.0496321, 0.049654, 0.0496759, 0.0496979, 0.0497198, \
0.0497419, 0.0497639, 0.049786, 0.0498081, 0.0498302, 0.0498524, \
0.0498746, 0.0498969, 0.0499192, 0.0499415, 0.0499638, 0.0499862, \
0.0500086, 0.050031, 0.0500535, 0.050076, 0.0500986, 0.0501212, \
0.0501438, 0.0501664, 0.0501891, 0.0502118, 0.0502345, 0.0502573, \
0.0502801, 0.050303, 0.0503259, 0.0503488, 0.0503717, 0.0503947, \
0.0504177, 0.0504408, 0.0504639, 0.050487, 0.0505102, 0.0505334, \
0.0505566, 0.0505799, 0.0506032, 0.0506265, 0.0506499, 0.0506733, \
0.0506967, 0.0507202, 0.0507437, 0.0507673, 0.0507909, 0.0508145, \
0.0508381, 0.0508618, 0.0508856, 0.0509093, 0.0509332, 0.050957, \
0.0509809, 0.0510048, 0.0510288, 0.0510527, 0.0510768, 0.0511008, \
0.0511249, 0.0511491, 0.0511733, 0.0511975, 0.0512217, 0.051246, \
0.0512704, 0.0512947, 0.0513192, 0.0513436, 0.0513681, 0.0513926, \
0.0514172, 0.0514418, 0.0514664, 0.0514911, 0.0515158, 0.0515406, \
0.0515654, 0.0515902, 0.0516151, 0.05164, 0.051665, 0.05169, \
0.051715, 0.0517401, 0.0517652, 0.0517904, 0.0518156, 0.0518408, \
0.0518661, 0.0518914, 0.0519168, 0.0519422, 0.0519677, 0.0519931, \
0.0520187, 0.0520443, 0.0520699, 0.0520955, 0.0521212, 0.052147, \
0.0521727, 0.0521986, 0.0522244, 0.0522503, 0.0522763, 0.0523023, \
0.0523283, 0.0523544, 0.0523805, 0.0524067, 0.0524329, 0.0524592, \
0.0524855, 0.0525118, 0.0525382, 0.0525646, 0.0525911, 0.0526176, \
0.0526442, 0.0526708, 0.0526975, 0.0527242, 0.0527509, 0.0527777, \
0.0528045, 0.0528314, 0.0528583, 0.0528853, 0.0529123, 0.0529394, \
0.0529665, 0.0529936, 0.0530208, 0.0530481, 0.0530754, 0.0531027, \
0.0531301, 0.0531576, 0.053185, 0.0532126, 0.0532402, 0.0532678, \
0.0532954, 0.0533232, 0.0533509, 0.0533788, 0.0534066, 0.0534345, \
0.0534625, 0.0534905, 0.0535186, 0.0535467, 0.0535749, 0.0536031, \
0.0536313, 0.0536596, 0.053688, 0.0537164, 0.0537449, 0.0537734, \
0.0538019, 0.0538305, 0.0538592, 0.0538879, 0.0539167, 0.0539455, \
0.0539744, 0.0540033, 0.0540322, 0.0540613, 0.0540903, 0.0541195, \
0.0541486, 0.0541779, 0.0542072, 0.0542365, 0.0542659, 0.0542953, \
0.0543248, 0.0543544, 0.054384, 0.0544136, 0.0544433, 0.0544731, \
0.0545029, 0.0545328, 0.0545627, 0.0545927, 0.0546228, 0.0546529, \
0.054683, 0.0547132, 0.0547435, 0.0547738, 0.0548042, 0.0548346, \
0.0548651, 0.0548956, 0.0549262, 0.0549569, 0.0549876, 0.0550184, \
0.0550492, 0.0550801, 0.0551111, 0.0551421, 0.0551731, 0.0552043, \
0.0552354, 0.0552667, 0.055298, 0.0553293, 0.0553608, 0.0553922, \
0.0554238, 0.0554554, 0.055487, 0.0555188, 0.0555505, 0.0555824, \
0.0556143, 0.0556463, 0.0556783, 0.0557104, 0.0557425, 0.0557747, \
0.055807, 0.0558393, 0.0558717, 0.0559042, 0.0559367, 0.0559693, \
0.056002, 0.0560347, 0.0560675, 0.0561003, 0.0561333, 0.0561662, \
0.0561993, 0.0562324, 0.0562656, 0.0562988, 0.0563321, 0.0563655, \
0.0563989, 0.0564324, 0.056466, 0.0564996, 0.0565333, 0.0565671, \
0.056601, 0.0566349, 0.0566688, 0.0567029, 0.056737, 0.0567712, \
0.0568054, 0.0568398, 0.0568742, 0.0569086, 0.0569432, 0.0569778, \
0.0570124, 0.0570472, 0.057082, 0.0571169, 0.0571519, 0.0571869, \
0.057222, 0.0572572, 0.0572924, 0.0573278, 0.0573632, 0.0573986, \
0.0574342, 0.0574698, 0.0575055, 0.0575413, 0.0575771, 0.057613, \
0.057649, 0.0576851, 0.0577213, 0.0577575, 0.0577938, 0.0578302, \
0.0578666, 0.0579032, 0.0579398, 0.0579765, 0.0580132, 0.0580501, \
0.058087, 0.058124, 0.0581611, 0.0581983, 0.0582355, 0.0582728, \
0.0583102, 0.0583477, 0.0583853, 0.058423, 0.0584607, 0.0584985, \
0.0585364, 0.0585744, 0.0586125, 0.0586506, 0.0586888, 0.0587272, \
0.0587656, 0.058804, 0.0588426, 0.0588813, 0.05892, 0.0589589, \
0.0589978, 0.0590368, 0.0590759, 0.0591151, 0.0591543, 0.0591937, \
0.0592331, 0.0592727, 0.0593123, 0.059352, 0.0593918, 0.0594317, \
0.0594717, 0.0595118, 0.0595519, 0.0595922, 0.0596326, 0.059673, \
0.0597135, 0.0597542, 0.0597949, 0.0598357, 0.0598767, 0.0599177, \
0.0599588, 0.06, 0.0600413, 0.0600827, 0.0601242, 0.0601658, \
0.0602074, 0.0602492, 0.0602911, 0.0603331, 0.0603752, 0.0604174, \
0.0604597, 0.060502, 0.0605445, 0.0605871, 0.0606298, 0.0606726, \
0.0607155, 0.0607585, 0.0608016, 0.0608448, 0.0608881, 0.0609315, \
0.0609751, 0.0610187, 0.0610624, 0.0611063, 0.0611502, 0.0611943, \
0.0612385, 0.0612827, 0.0613271, 0.0613716, 0.0614162, 0.0614609, \
0.0615058, 0.0615507, 0.0615958, 0.0616409, 0.0616862, 0.0617316, \
0.0617771, 0.0618227, 0.0618684, 0.0619143, 0.0619603, 0.0620064, \
0.0620526, 0.0620989, 0.0621453, 0.0621919, 0.0622385, 0.0622853, \
0.0623323, 0.0623793, 0.0624265, 0.0624737, 0.0625212, 0.0625687, \
0.0626163, 0.0626641, 0.062712, 0.06276, 0.0628082, 0.0628565, \
0.0629049, 0.0629534, 0.0630021, 0.0630509, 0.0630998, 0.0631489, \
0.0631981, 0.0632474, 0.0632968, 0.0633464, 0.0633961, 0.063446, \
0.063496, 0.0635461, 0.0635964, 0.0636468, 0.0636973, 0.063748, \
0.0637988, 0.0638497, 0.0639008, 0.063952, 0.0640034, 0.0640549, \
0.0641066, 0.0641584, 0.0642103, 0.0642624, 0.0643147, 0.0643671, \
0.0644196, 0.0644723, 0.0645251, 0.0645781, 0.0646312, 0.0646845, \
0.0647379, 0.0647915, 0.0648453, 0.0648992, 0.0649532, 0.0650074, \
0.0650618, 0.0651163, 0.065171, 0.0652258, 0.0652808, 0.065336, \
0.0653913, 0.0654468, 0.0655024, 0.0655582, 0.0656142, 0.0656703, \
0.0657266, 0.0657831, 0.0658398, 0.0658966, 0.0659535, 0.0660107, \
0.066068, 0.0661255, 0.0661832, 0.0662411, 0.0662991, 0.0663573, \
0.0664157, 0.0664742, 0.066533, 0.0665919, 0.066651, 0.0667103, \
0.0667697, 0.0668294, 0.0668892, 0.0669493, 0.0670095, 0.0670699, \
0.0671305, 0.0671913, 0.0672523, 0.0673134, 0.0673748, 0.0674364, \
0.0674981, 0.0675601, 0.0676222, 0.0676846, 0.0677471, 0.0678099, \
0.0678729, 0.067936, 0.0679994, 0.068063, 0.0681268, 0.0681908, \
0.068255, 0.0683194, 0.0683841, 0.0684489, 0.068514, 0.0685793, \
0.0686448, 0.0687105, 0.0687765, 0.0688426, 0.068909, 0.0689757, \
0.0690425, 0.0691096, 0.0691769, 0.0692444, 0.0693122, 0.0693802, \
0.0694484, 0.0695169, 0.0695856, 0.0696546, 0.0697238, 0.0697932, \
0.0698629, 0.0699328, 0.070003, 0.0700734, 0.0701441, 0.0702151, \
0.0702862, 0.0703577, 0.0704294, 0.0705013, 0.0705735, 0.070646, \
0.0707188, 0.0707918, 0.0708651, 0.0709386, 0.0710124, 0.0710865, \
0.0711609, 0.0712355, 0.0713105, 0.0713857, 0.0714612, 0.0715369, \
0.071613, 0.0716893, 0.071766, 0.0718429, 0.0719201, 0.0719976, \
0.0720755, 0.0721536, 0.072232, 0.0723107, 0.0723898, 0.0724691, \
0.0725488, 0.0726287, 0.072709, 0.0727896, 0.0728705, 0.0729518, \
0.0730333, 0.0731152, 0.0731975, 0.07328, 0.0733629, 0.0734462, \
0.0735297, 0.0736136, 0.0736979, 0.0737825, 0.0738675, 0.0739528, \
0.0740384, 0.0741245, 0.0742109, 0.0742976, 0.0743847, 0.0744722, \
0.0745601, 0.0746483, 0.0747369, 0.0748259, 0.0749153, 0.075005, \
0.0750952, 0.0751857, 0.0752767, 0.075368, 0.0754597, 0.0755519, \
0.0756445, 0.0757374, 0.0758308, 0.0759246, 0.0760188, 0.0761135, \
0.0762086, 0.0763041, 0.0764001, 0.0764965, 0.0765933, 0.0766906, \
0.0767883, 0.0768865, 0.0769852, 0.0770843, 0.0771839, 0.077284, \
0.0773845, 0.0774855, 0.077587, 0.077689, 0.0777915, 0.0778945, \
0.077998, 0.078102, 0.0782065, 0.0783115, 0.078417, 0.0785231, \
0.0786297, 0.0787368, 0.0788445, 0.0789527, 0.0790614, 0.0791708, \
0.0792806, 0.0793911, 0.0795021, 0.0796137, 0.0797259, 0.0798386, \
0.079952, 0.080066, 0.0801805, 0.0802957, 0.0804115, 0.0805279, \
0.0806449, 0.0807626, 0.0808809, 0.0809999, 0.0811195, 0.0812398, \
0.0813608, 0.0814824, 0.0816048, 0.0817278, 0.0818515, 0.0819759, \
0.082101, 0.0822269, 0.0823534, 0.0824808, 0.0826088, 0.0827376, \
0.0828672, 0.0829975, 0.0831286, 0.0832605, 0.0833932, 0.0835267, \
0.083661, 0.0837962, 0.0839321, 0.0840689, 0.0842066, 0.0843451, \
0.0844845, 0.0846248, 0.084766, 0.084908, 0.085051, 0.0851949, \
0.0853398, 0.0854856, 0.0856323, 0.08578, 0.0859287, 0.0860784, \
0.0862292, 0.0863809, 0.0865336, 0.0866875, 0.0868423, 0.0869983, \
0.0871553, 0.0873135, 0.0874727, 0.0876331, 0.0877947, 0.0879574, \
0.0881213, 0.0882864, 0.0884527, 0.0886202, 0.088789, 0.0889591, \
0.0891304, 0.0893031, 0.089477, 0.0896523, 0.089829, 0.090007, \
0.0901865, 0.0903674, 0.0905497, 0.0907335, 0.0909188, 0.0911056, \
0.0912939, 0.0914838, 0.0916753, 0.0918684, 0.0920631, 0.0922595, \
0.0924576, 0.0926574, 0.092859, 0.0930623, 0.0932675, 0.0934745, \
0.0936834, 0.0938942, 0.0941069, 0.0943216, 0.0945384, 0.0947571, \
0.094978, 0.095201, 0.0954262, 0.0956536, 0.0958832, 0.0961151, \
0.0963494, 0.096586, 0.0968251, 0.0970667, 0.0973109, 0.0975576, \
0.0978069, 0.098059, 0.0983138, 0.0985715, 0.098832, 0.0990955, \
0.099362, 0.0996316, 0.0999043, 0.10018, 0.10046, 0.100742, 0.101028, \
0.101318, 0.101611, 0.101908, 0.102209, 0.102514, 0.102823, 0.103136, \
0.103454, 0.103775, 0.104102, 0.104432, 0.104768, 0.105109, 0.105454, \
0.105805, 0.106161, 0.106523, 0.106891, 0.107264, 0.107644, 0.10803, \
0.108422, 0.108822, 0.109228, 0.109642, 0.110063, 0.110493, 0.11093, \
0.111376, 0.111831, 0.112295, 0.112769, 0.113254, 0.113748, 0.114254, \
0.114772, 0.115302, 0.115845, 0.116401, 0.116972, 0.117559, 0.118161, \
0.118781, 0.119418, 0.120076, 0.120754, 0.121454, 0.122179, 0.12293, \
0.123709, 0.124519, 0.125362, 0.126243, 0.127165, 0.128132, 0.129151, \
0.130228, 0.131371, 0.132592, 0.133904, 0.135325, 0.136881, 0.138607, \
0.14056, 0.142833,      0.1456111,0.1493707,0.1521282,0.1579619};

static const double QNMData_fdamp[] = {0.0140098,0.0140102,0.0140106,0.0140114,       0.0140177, 0.0140154, 0.0140148, 0.014015, 0.0140156, \
0.0140164, 0.0140172, 0.0140181, 0.0140189, 0.0140198, 0.0140206, \
0.0140214, 0.0140223, 0.0140231, 0.0140239, 0.0140247, 0.0140256, \
0.0140264, 0.0140272, 0.014028, 0.0140288, 0.0140296, 0.0140305, \
0.0140313, 0.0140321, 0.0140329, 0.0140337, 0.0140345, 0.0140353, \
0.0140361, 0.0140369, 0.0140377, 0.0140385, 0.0140393, 0.0140401, \
0.0140409, 0.0140417, 0.0140425, 0.0140433, 0.014044, 0.0140448, \
0.0140456, 0.0140464, 0.0140472, 0.014048, 0.0140487, 0.0140495, \
0.0140503, 0.0140511, 0.0140519, 0.0140526, 0.0140534, 0.0140542, \
0.0140549, 0.0140557, 0.0140565, 0.0140572, 0.014058, 0.0140587, \
0.0140595, 0.0140603, 0.014061, 0.0140618, 0.0140625, 0.0140633, \
0.014064, 0.0140648, 0.0140655, 0.0140663, 0.014067, 0.0140677, \
0.0140685, 0.0140692, 0.01407, 0.0140707, 0.0140714, 0.0140722, \
0.0140729, 0.0140736, 0.0140744, 0.0140751, 0.0140758, 0.0140765, \
0.0140772, 0.014078, 0.0140787, 0.0140794, 0.0140801, 0.0140808, \
0.0140815, 0.0140822, 0.0140829, 0.0140837, 0.0140844, 0.0140851, \
0.0140858, 0.0140865, 0.0140872, 0.0140879, 0.0140885, 0.0140892, \
0.0140899, 0.0140906, 0.0140913, 0.014092, 0.0140927, 0.0140934, \
0.014094, 0.0140947, 0.0140954, 0.0140961, 0.0140967, 0.0140974, \
0.0140981, 0.0140988, 0.0140994, 0.0141001, 0.0141007, 0.0141014, \
0.0141021, 0.0141027, 0.0141034, 0.014104, 0.0141047, 0.0141053, \
0.014106, 0.0141066, 0.0141073, 0.0141079, 0.0141086, 0.0141092, \
0.0141098, 0.0141105, 0.0141111, 0.0141117, 0.0141124, 0.014113, \
0.0141136, 0.0141142, 0.0141149, 0.0141155, 0.0141161, 0.0141167, \
0.0141173, 0.014118, 0.0141186, 0.0141192, 0.0141198, 0.0141204, \
0.014121, 0.0141216, 0.0141222, 0.0141228, 0.0141234, 0.014124, \
0.0141246, 0.0141252, 0.0141257, 0.0141263, 0.0141269, 0.0141275, \
0.0141281, 0.0141287, 0.0141292, 0.0141298, 0.0141304, 0.0141309, \
0.0141315, 0.0141321, 0.0141326, 0.0141332, 0.0141338, 0.0141343, \
0.0141349, 0.0141354, 0.014136, 0.0141365, 0.0141371, 0.0141376, \
0.0141382, 0.0141387, 0.0141392, 0.0141398, 0.0141403, 0.0141408, \
0.0141414, 0.0141419, 0.0141424, 0.0141429, 0.0141435, 0.014144, \
0.0141445, 0.014145, 0.0141455, 0.014146, 0.0141466, 0.0141471, \
0.0141476, 0.0141481, 0.0141486, 0.0141491, 0.0141496, 0.01415, \
0.0141505, 0.014151, 0.0141515, 0.014152, 0.0141525, 0.014153, \
0.0141534, 0.0141539, 0.0141544, 0.0141549, 0.0141553, 0.0141558, \
0.0141562, 0.0141567, 0.0141572, 0.0141576, 0.0141581, 0.0141585, \
0.014159, 0.0141594, 0.0141599, 0.0141603, 0.0141608, 0.0141612, \
0.0141616, 0.0141621, 0.0141625, 0.0141629, 0.0141633, 0.0141638, \
0.0141642, 0.0141646, 0.014165, 0.0141654, 0.0141658, 0.0141662, \
0.0141667, 0.0141671, 0.0141675, 0.0141679, 0.0141683, 0.0141686, \
0.014169, 0.0141694, 0.0141698, 0.0141702, 0.0141706, 0.014171, \
0.0141713, 0.0141717, 0.0141721, 0.0141724, 0.0141728, 0.0141732, \
0.0141735, 0.0141739, 0.0141742, 0.0141746, 0.0141749, 0.0141753, \
0.0141756, 0.014176, 0.0141763, 0.0141766, 0.014177, 0.0141773, \
0.0141776, 0.014178, 0.0141783, 0.0141786, 0.0141789, 0.0141792, \
0.0141796, 0.0141799, 0.0141802, 0.0141805, 0.0141808, 0.0141811, \
0.0141814, 0.0141817, 0.0141819, 0.0141822, 0.0141825, 0.0141828, \
0.0141831, 0.0141833, 0.0141836, 0.0141839, 0.0141842, 0.0141844, \
0.0141847, 0.0141849, 0.0141852, 0.0141854, 0.0141857, 0.0141859, \
0.0141862, 0.0141864, 0.0141867, 0.0141869, 0.0141871, 0.0141874, \
0.0141876, 0.0141878, 0.014188, 0.0141882, 0.0141884, 0.0141887, \
0.0141889, 0.0141891, 0.0141893, 0.0141895, 0.0141897, 0.0141899, \
0.01419, 0.0141902, 0.0141904, 0.0141906, 0.0141908, 0.0141909, \
0.0141911, 0.0141913, 0.0141914, 0.0141916, 0.0141917, 0.0141919, \
0.014192, 0.0141922, 0.0141923, 0.0141925, 0.0141926, 0.0141927, \
0.0141929, 0.014193, 0.0141931, 0.0141932, 0.0141934, 0.0141935, \
0.0141936, 0.0141937, 0.0141938, 0.0141939, 0.014194, 0.0141941, \
0.0141942, 0.0141942, 0.0141943, 0.0141944, 0.0141945, 0.0141946, \
0.0141946, 0.0141947, 0.0141947, 0.0141948, 0.0141949, 0.0141949, \
0.0141949, 0.014195, 0.014195, 0.0141951, 0.0141951, 0.0141951, \
0.0141951, 0.0141952, 0.0141952, 0.0141952, 0.0141952, 0.0141952, \
0.0141952, 0.0141952, 0.0141952, 0.0141952, 0.0141952, 0.0141952, \
0.0141951, 0.0141951, 0.0141951, 0.014195, 0.014195, 0.014195, \
0.0141949, 0.0141949, 0.0141948, 0.0141948, 0.0141947, 0.0141946, \
0.0141946, 0.0141945, 0.0141944, 0.0141943, 0.0141942, 0.0141941, \
0.014194, 0.0141939, 0.0141938, 0.0141937, 0.0141936, 0.0141935, \
0.0141934, 0.0141933, 0.0141931, 0.014193, 0.0141929, 0.0141927, \
0.0141926, 0.0141924, 0.0141923, 0.0141921, 0.0141919, 0.0141918, \
0.0141916, 0.0141914, 0.0141912, 0.0141911, 0.0141909, 0.0141907, \
0.0141905, 0.0141903, 0.0141901, 0.0141898, 0.0141896, 0.0141894, \
0.0141892, 0.0141889, 0.0141887, 0.0141885, 0.0141882, 0.014188, \
0.0141877, 0.0141874, 0.0141872, 0.0141869, 0.0141866, 0.0141863, \
0.014186, 0.0141858, 0.0141855, 0.0141852, 0.0141848, 0.0141845, \
0.0141842, 0.0141839, 0.0141836, 0.0141832, 0.0141829, 0.0141826, \
0.0141822, 0.0141819, 0.0141815, 0.0141811, 0.0141808, 0.0141804, \
0.01418, 0.0141796, 0.0141792, 0.0141788, 0.0141784, 0.014178, \
0.0141776, 0.0141772, 0.0141768, 0.0141763, 0.0141759, 0.0141755, \
0.014175, 0.0141746, 0.0141741, 0.0141736, 0.0141732, 0.0141727, \
0.0141722, 0.0141717, 0.0141712, 0.0141707, 0.0141702, 0.0141697, \
0.0141692, 0.0141687, 0.0141681, 0.0141676, 0.0141671, 0.0141665, \
0.0141659, 0.0141654, 0.0141648, 0.0141642, 0.0141637, 0.0141631, \
0.0141625, 0.0141619, 0.0141613, 0.0141607, 0.01416, 0.0141594, \
0.0141588, 0.0141582, 0.0141575, 0.0141569, 0.0141562, 0.0141555, \
0.0141549, 0.0141542, 0.0141535, 0.0141528, 0.0141521, 0.0141514, \
0.0141507, 0.01415, 0.0141492, 0.0141485, 0.0141478, 0.014147, \
0.0141462, 0.0141455, 0.0141447, 0.0141439, 0.0141431, 0.0141424, \
0.0141416, 0.0141407, 0.0141399, 0.0141391, 0.0141383, 0.0141374, \
0.0141366, 0.0141357, 0.0141349, 0.014134, 0.0141331, 0.0141322, \
0.0141313, 0.0141304, 0.0141295, 0.0141286, 0.0141277, 0.0141267, \
0.0141258, 0.0141248, 0.0141239, 0.0141229, 0.0141219, 0.014121, \
0.01412, 0.014119, 0.014118, 0.0141169, 0.0141159, 0.0141149, \
0.0141138, 0.0141128, 0.0141117, 0.0141106, 0.0141095, 0.0141085, \
0.0141074, 0.0141062, 0.0141051, 0.014104, 0.0141029, 0.0141017, \
0.0141006, 0.0140994, 0.0140982, 0.014097, 0.0140958, 0.0140946, \
0.0140934, 0.0140922, 0.014091, 0.0140897, 0.0140885, 0.0140872, \
0.0140859, 0.0140846, 0.0140833, 0.014082, 0.0140807, 0.0140794, \
0.0140781, 0.0140767, 0.0140753, 0.014074, 0.0140726, 0.0140712, \
0.0140698, 0.0140684, 0.014067, 0.0140655, 0.0140641, 0.0140626, \
0.0140612, 0.0140597, 0.0140582, 0.0140567, 0.0140552, 0.0140536, \
0.0140521, 0.0140505, 0.014049, 0.0140474, 0.0140458, 0.0140442, \
0.0140426, 0.014041, 0.0140393, 0.0140377, 0.014036, 0.0140343, \
0.0140327, 0.014031, 0.0140292, 0.0140275, 0.0140258, 0.014024, \
0.0140223, 0.0140205, 0.0140187, 0.0140169, 0.0140151, 0.0140132, \
0.0140114, 0.0140095, 0.0140076, 0.0140057, 0.0140038, 0.0140019, \
0.014, 0.013998, 0.0139961, 0.0139941, 0.0139921, 0.0139901, \
0.0139881, 0.013986, 0.013984, 0.0139819, 0.0139798, 0.0139777, \
0.0139756, 0.0139735, 0.0139713, 0.0139691, 0.013967, 0.0139648, \
0.0139625, 0.0139603, 0.0139581, 0.0139558, 0.0139535, 0.0139512, \
0.0139489, 0.0139466, 0.0139442, 0.0139419, 0.0139395, 0.0139371, \
0.0139346, 0.0139322, 0.0139297, 0.0139273, 0.0139248, 0.0139223, \
0.0139197, 0.0139172, 0.0139146, 0.013912, 0.0139094, 0.0139068, \
0.0139041, 0.0139014, 0.0138987, 0.013896, 0.0138933, 0.0138906, \
0.0138878, 0.013885, 0.0138822, 0.0138793, 0.0138765, 0.0138736, \
0.0138707, 0.0138678, 0.0138648, 0.0138619, 0.0138589, 0.0138559, \
0.0138528, 0.0138498, 0.0138467, 0.0138436, 0.0138405, 0.0138373, \
0.0138341, 0.0138309, 0.0138277, 0.0138244, 0.0138212, 0.0138179, \
0.0138145, 0.0138112, 0.0138078, 0.0138044, 0.013801, 0.0137975, \
0.013794, 0.0137905, 0.013787, 0.0137834, 0.0137798, 0.0137762, \
0.0137726, 0.0137689, 0.0137652, 0.0137615, 0.0137577, 0.0137539, \
0.0137501, 0.0137462, 0.0137424, 0.0137384, 0.0137345, 0.0137305, \
0.0137265, 0.0137225, 0.0137184, 0.0137143, 0.0137102, 0.013706, \
0.0137018, 0.0136976, 0.0136933, 0.013689, 0.0136847, 0.0136803, \
0.0136759, 0.0136715, 0.013667, 0.0136625, 0.0136579, 0.0136533, \
0.0136487, 0.0136441, 0.0136394, 0.0136346, 0.0136298, 0.013625, \
0.0136202, 0.0136153, 0.0136103, 0.0136054, 0.0136004, 0.0135953, \
0.0135902, 0.0135851, 0.0135799, 0.0135747, 0.0135694, 0.0135641, \
0.0135587, 0.0135533, 0.0135479, 0.0135424, 0.0135368, 0.0135312, \
0.0135256, 0.0135199, 0.0135142, 0.0135084, 0.0135026, 0.0134967, \
0.0134907, 0.0134848, 0.0134787, 0.0134726, 0.0134665, 0.0134603, \
0.0134541, 0.0134478, 0.0134414, 0.013435, 0.0134285, 0.013422, \
0.0134154, 0.0134088, 0.0134021, 0.0133953, 0.0133885, 0.0133816, \
0.0133747, 0.0133676, 0.0133606, 0.0133534, 0.0133462, 0.013339, \
0.0133316, 0.0133242, 0.0133168, 0.0133092, 0.0133016, 0.013294, \
0.0132862, 0.0132784, 0.0132705, 0.0132625, 0.0132545, 0.0132464, \
0.0132382, 0.0132299, 0.0132216, 0.0132131, 0.0132046, 0.013196, \
0.0131874, 0.0131786, 0.0131697, 0.0131608, 0.0131518, 0.0131427, \
0.0131335, 0.0131242, 0.0131148, 0.0131054, 0.0130958, 0.0130861, \
0.0130764, 0.0130665, 0.0130566, 0.0130465, 0.0130363, 0.0130261, \
0.0130157, 0.0130052, 0.0129946, 0.0129839, 0.0129731, 0.0129622, \
0.0129512, 0.01294, 0.0129288, 0.0129174, 0.0129059, 0.0128942, \
0.0128825, 0.0128706, 0.0128586, 0.0128464, 0.0128342, 0.0128218, \
0.0128092, 0.0127965, 0.0127837, 0.0127707, 0.0127576, 0.0127443, \
0.0127309, 0.0127174, 0.0127036, 0.0126898, 0.0126757, 0.0126615, \
0.0126472, 0.0126326, 0.0126179, 0.0126031, 0.012588, 0.0125728, \
0.0125574, 0.0125418, 0.012526, 0.01251, 0.0124938, 0.0124774, \
0.0124609, 0.0124441, 0.0124271, 0.0124099, 0.0123925, 0.0123748, \
0.012357, 0.0123389, 0.0123205, 0.012302, 0.0122832, 0.0122641, \
0.0122448, 0.0122252, 0.0122054, 0.0121853, 0.012165, 0.0121443, \
0.0121234, 0.0121022, 0.0120807, 0.0120589, 0.0120368, 0.0120144, \
0.0119917, 0.0119686, 0.0119452, 0.0119215, 0.0118974, 0.011873, \
0.0118482, 0.011823, 0.0117975, 0.0117716, 0.0117452, 0.0117185, \
0.0116914, 0.0116638, 0.0116358, 0.0116074, 0.0115785, 0.0115492, \
0.0115193, 0.011489, 0.0114582, 0.0114269, 0.011395, 0.0113626, \
0.0113296, 0.0112961, 0.011262, 0.0112273, 0.011192, 0.0111561, \
0.0111195, 0.0110822, 0.0110442, 0.0110056, 0.0109662, 0.0109261, \
0.0108852, 0.0108435, 0.0108009, 0.0107576, 0.0107133, 0.0106682, \
0.0106221, 0.0105751, 0.010527, 0.010478, 0.0104278, 0.0103766, \
0.0103243, 0.0102707, 0.0102159, 0.0101599, 0.0101025, 0.0100438, \
0.00998358, 0.00992191, 0.00985869, 0.00979385, 0.00972732, \
0.00965903, 0.00958889, 0.00951682, 0.00944272, 0.0093665, \
0.00928805, 0.00920725, 0.00912399, 0.00903811, 0.00894949, \
0.00885795, 0.00876332, 0.00866542, 0.00856403, 0.00845893, \
0.00834985, 0.00823651, 0.0081186, 0.00799576, 0.0078676, 0.00773366, \
0.00759343, 0.00744631, 0.00729164, 0.00712861, 0.00695629, \
0.00677358, 0.00657914, 0.00637134, 0.00614819, 0.00590712, \
0.00564485, 0.00535699, 0.0050375, 0.00467763, 0.00426389, 0.00377349,    0.0031618,0.0023131,0.0016762,0.0002908};

/* Non-reviewd functions */

int IMRPhenomDPhaseFrequencySequence(
    REAL8Sequence *phases,
    REAL8Sequence *freqs,
    size_t ind_min,
    size_t ind_max,
    REAL8 m1,
    REAL8 m2,
    REAL8 chi1z,
    REAL8 chi2z,
    REAL8 Rholm,
    REAL8 Taulm,
    LALDict *extraParams
);

int IMRPhenomDAmpFrequencySequence(
    REAL8Sequence *phases,
    REAL8Sequence *freqs,
    size_t ind_min,
    size_t ind_max,
    REAL8 m1,
    REAL8 m2,
    REAL8 chi1z,
    REAL8 chi2z
);

#endif /* _LALSIM_IMR_PHENOMD_H */
