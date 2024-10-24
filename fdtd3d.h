#include <string>
#include <cmath>
#include <eigen3/Eigen/Dense>

extern std::string global_dirName;
const std::string PATH { "/home/itabashi/OneDrive/Lab/3D_Spherical_FDTD_HJ/" };

constexpr double C0 { 3.0e8 };
constexpr double MU0 { 4.0 * M_PI * 1.0e-7 };
constexpr double EPS0 { 1.0 / MU0 / C0 / C0 };
constexpr double CHARGE_e { 1.602e-19 };
constexpr double MASS_e { 9.1e-31 };
// constexpr double B0 { 50000e-9 };
constexpr double B0 { 0.0 };
constexpr double R0 { 6370.0e3 }; /* Radius of the Earth */

constexpr int Year { 2015 };
constexpr int Month { 6 };
constexpr int Day { 18 };
/* えびの 緯度経度 */
constexpr double Tx_Latitude { 32.067650 };
constexpr double Tx_Longitude { 130.828978 };

/* 調布 緯度経度 */
constexpr double Rx_Latitude { 35.65966 };
constexpr double Rx_Longitude { 139.54328 };

/* 解析領域 */
constexpr double Rr { 100.0e3 };
constexpr double Rth { 100.0e3 };
constexpr double Rph { 100.0e3 };
// constexpr double Rph { 300.0e3 };

// constexpr double Rr { 90.0e3 };
// constexpr double Rth { 50.0e3 };
// constexpr double Rph { 1000.0e3 };
// constexpr double dr { 0.5e3 };
// constexpr double Rdth { 0.5e3 };
// constexpr double Rdph { 0.5e3 };
constexpr double dr { 1.0e3 };
// constexpr double dr { 0.5e3 };
constexpr double Rdth { 1.0e3 };
constexpr double Rdph { 1.0e3 };
constexpr double dth { Rdth / R0 };
constexpr double dph { Rdph / R0 };
constexpr double thR { Rth / R0 };
constexpr double phR { Rph / R0 };

constexpr int Nr { int(Rr / dr) };
constexpr int Nth { int(Rth / Rdth) };
constexpr int Nph { int(Rph / Rdph) };

/*xiの導入が必要*/
// constexpr double f0 { 40.0e3 };
constexpr double f0 { 22.2e3 };
constexpr double Tmax { 0.005 };
// constexpr double Tmax { 0.01 };
constexpr double Ne_max { 1.0e10 };
constexpr double Omega { std::sqrt( CHARGE_e * CHARGE_e * Ne_max / MASS_e / EPS0 ) };
constexpr double xi { 1.0 / std::sqrt( 1.0 + Omega * Omega / 4.0 / C0 / C0 / ( 1.0 / dr / dr + 1.0 / ( R0 * dth ) / ( R0 * dth ) 
                    + 1.0 / ( R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph ) / ( R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph ) ) ) };
// constexpr double dt { xi * 0.9 / C0 / sqrt( 1.0 / dr / dr + 1.0 / (R0 * dth) / (R0 * dth) 
                        // + 1.0 / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph) / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph)  ) };
// constexpr double dt { xi * 0.5 / C0 / sqrt( 1.0 / dr / dr + 1.0 / (R0 * dth) / (R0 * dth) 
//                         + 1.0 / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph) / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph)  ) };
constexpr double dt {  0.9  / C0 / sqrt( 1.0 / dr / dr + 1.0 / (R0 * dth) / (R0 * dth) 
                        + 1.0 / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph) / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph)  ) };
constexpr int Nt { int(Tmax / dt) };

constexpr double OMG { f0 * 2. * M_PI }; /* 角周波数 */

/* PMLパラメタ */
constexpr int PML_L { 12 };
constexpr double PML_M { 4.0 };
constexpr double PML_R { 1.0e-6 };

/* 電流源パラメタ */
// constexpr double sigma_J { 12 * dt };
// constexpr double current_dt { 7.29756e-07 };
constexpr double t0 { 1.0 / f0 };
// constexpr double sigma_J { 18.0 * current_dt };
constexpr double sigma_J { 1.0 / 2.0 / M_PI / f0 };
// constexpr double t0 { 6.0 * sigma_J };
// constexpr double source_r { (PML_L + 1) * dr };
// constexpr double source_th { Rth / 2.0 };
// constexpr double source_ph { 50.0e3 };
constexpr double source_r { Rr / 2.0 };
constexpr double source_th { Rth / 2.0 };
constexpr double source_ph { Rph / 2.0 };

/* プラズマ領域パラメタ */
constexpr double Rr_iono_lower { 55.0e3 };
// constexpr double Rr_iono_lower { Rr / 2.0 };
constexpr double Rr_iono_upper { Rr };
constexpr double Rth_iono_lower { 0.0 };
constexpr double Rth_iono_upper { Rth };
constexpr double Rph_iono_lower { 0.0 };
constexpr double Rph_iono_upper { Rph };

// constexpr double Rr_iono_lower { 25.0e3 };
// constexpr double Rr_iono_upper { 50.0e3 };
// constexpr double Rth_iono_lower { 0.0 };
// constexpr double Rth_iono_upper { Rth };
// constexpr double Rph_iono_lower { 0.0 };
// constexpr double Rph_iono_upper { Rph };

constexpr int Nr_iono_lower { int(Rr_iono_lower / dr) + PML_L};
constexpr int Nr_iono_upper { int(Rr_iono_upper / dr) - 1 };
constexpr int Nth_iono_lower { int(Rth_iono_lower / Rdth)};
constexpr int Nth_iono_upper { int(Rth_iono_upper / Rdth) - 1};
constexpr int Nph_iono_lower { int(Rph_iono_lower / Rdph)};
constexpr int Nph_iono_upper { int(Rph_iono_upper / Rdph) - 1};
constexpr int Nr_iono { Nr_iono_upper - Nr_iono_lower + 1 };
constexpr int Nth_iono { Nth_iono_upper - Nth_iono_lower + 1 };
constexpr int Nph_iono { Nph_iono_upper - Nph_iono_lower + 1 };

constexpr double exp_Ne { 10 };
// constexpr double exp_Ne { 10.875 };
constexpr double exp_nu { 7.0 };

// constexpr double THETA { M_PI / 4.0 };
constexpr double THETA { 0.0 };
constexpr double PHI { M_PI / 4.0 };

/*観測パラメタ*/
// constexpr double obs_r { Rr_iono_lower - 5.0e3 };
// constexpr double obs_th { source_th };
// constexpr double obs_ph { 140.0e3 };
constexpr double obs_ph { source_ph };
constexpr double obs_r { source_r };
constexpr double obs_th { source_th };
// constexpr double obs_ph { 900.0e3 };
constexpr double obs_t_step{ 1.0e-3 };
constexpr int obs_Nr { int(obs_r / dr) };
constexpr int obs_Nth { int(obs_th / Rdth) };
constexpr int obs_Nph { int(obs_ph / Rdph) };


void update_Er(double ***Er, double ****Hth, double ****Hph, double ****Jr, double ****check, int n);
void update_Eth(double ****Eth, double ***Hr, double ****Hph, double ****Jth, double ****check, int n);
void update_Eph(double ****Eph, double ***Hr, double ****Hth, double ****Jph, double ****check, int n);
void update_Er_PML(double ***Er, double ****Dr, double ****Hth, double ****Hph, double ****Jr, double ****check, int n);
void update_Eth_PML(double ****Eth, double ****Dth, double ***Hr, double ****Hph, double ****Jth, double ****check, int n);
void update_Eph_PML(double ****Eph, double ****Dph, double ***Hr, double ****Hth, double ****Jph, double ****check, int n);
void update_Eth_tilde(double ****Eth_tilde, double ****Eth, double *CETH_TILDE_00, double ****check, int n);
void update_Eph_tilde(double ****Eph_tilde, double ****Eph, double *CEPH_TILDE_00, double ****check, int n);

void update_Dr_PML(double ****Drth1, double ****Drth2, double ****Drph, double ****Dr, double ***Hr, double ****Hth, double ****Hph,
                     double *CDRTH1_00, double *CDRTH1_01, double *CDRPH_00, double *CDRPH_01, double ****check, int n);
void update_Dth_PML(double ****Dthph, double ****Dthr, double ****Dthr_tilde, double ****Dth, double ***Hr, double ***Hph_tilde, double *CDTHPH_00, double *CDTHPH_01,
                     double *CDTHR_10, double *CDTHR_11, double *CDTHR_TILDE_00, double *CDTHR_TILDE_01, double ****check, int n);
void update_Dph_PML(double ****Dphr, double ****Dphr_tilde, double ****Dphth, double ****Dph, double ***Hr, double ***Hth_tilde,
                     double *CDPHR_10, double *CDPHR_11, double *CDPHR_TILDE_00, double *CDPHR_TILDE_01, double *CDPHTH_00,
                     double *CDPHTH_01, double ****check, int n);

void update_Hr(double ***Hr, double ****Eth, double ****Eph, double ****check, int n);
void update_Hth(double ****Hth, double ***Er, double ****Eph, double ****check, int n);
void update_Hph(double ****Hph, double ***Er, double ****Eth, double ****check, int n);
void update_Hr_PML(double ***Hr, double ***Hrth1, double ***Hrth2, double ***Hrph,
                     double ****Eth, double ****Eph, double *CHRTH1_00, double *CHRTH1_01,
                      double *CHRPH_00, double *CHRPH_01, double ****check, int n);
void update_Hth_PML(double ****Hth, double ***Hthr, double ***Hthph, double ****Hthr_tilde, double ***Er, double ****Eph_tilde, double *CHTHPH_00,
                 double *CHTHPH_01, double *CHTHR_TILDE_00, double *CHTHR_TILDE_01, double *CHTHR_10, double *CHTHR_11, double ****check, int n);
void update_Hph_PML(double ****Hph, double ***Hphr, double ***Hphth, double ****Hphr_tilde, double ***Er, double ****Eth_tilde,
                     double *CHPHTH_00, double *CHPHTH_01, double *CHPHR_TILDE_00, double *CHPHR_TILDE_01, double *CHPHR_10, double *CHPHR_11, double ****check, int n);
void update_Hth_tilde(double ***Hth_tilde, double ****Hth, double *CHTH_TILDE, double ****check, int n);
void update_Hph_tilde(double ***Hph_tilde, double ****Hph, double *CHPH_TILDE, double ****chcek, int n);

void update_Jr(double ****Jr, double ****Jth, double ****Jph, double ***Er, double ****Eth, double ****Eph, Eigen::Matrix3d ***S, Eigen::Matrix3d ***B, int n);
void update_Jth(double ****Jr, double ****Jth, double ****Jph, double ***Er, double ****Eth, double ****Eph, Eigen::Matrix3d ***S, Eigen::Matrix3d ***B, int n);
void update_Jph(double ****Jr, double ****Jth, double ****Jph, double ***Er, double ****Eth, double ****Eph, Eigen::Matrix3d ***S, Eigen::Matrix3d ***B, int n);

double PML_sigma_r(double i);
// double PML_sigma_r_0(double i);
// double PML_sigma_r_1(double i);
double PML_sigma_th(double j);
double PML_sigma_ph(double k);
double PML_sigma_r_tilde(double i);
// double PML_sigma_r_tilde_0(double i);
// double PML_sigma_r_tilde_1(double i);
double PML_C00(double sigma);
double PML_C01(double sigma);
double PML_C10(double r, double sigma);
double PML_C11(double r, double sigma);
void initialize_PML(double *CERTH1_00, double *CERTH1_01, double *CERPH_00, double *CERPH_01, double *CETHPH_00, double *CETHPH_01, double *CETHR_10, double *CETHR_11,
                     double *CETHR_TILDE_00, double *CETHR_TILDE_01, double *CEPHR_10, double *CEPHR_11, double *CEPHR_TILDE_00, double *CEPHR_TILDE_01, double *CEPHTH_00,
                     double *CEPHTH_01, double *CETH_TILDE_00, double *CETH_TILDE_01,
                     double *CEPH_TILDE_00, double *CEPH_TILDE_01, double *CHRTH1_00, double *CHRTH1_01, double *CHRPH_00, double *CHRPH_01, double *CHTHPH_00,
                     double *CHTHPH_01, double *CHTHR_TILDE_00, double *CHTHR_TILDE_01, double *CHTHR_10, double *CHTHR_11,
                     double *CHPHTH_00, double *CHPHTH_01, double *CHPHR_TILDE_00, double *CHPHR_TILDE_01, double *CHPHR_10, double *CHPHR_11, double *CHTH_TILDE, double *CHPH_TILDE);

double source_J(double t);
double cal_Ne(double exp_Ne);
double cal_nu(double exp_nu);
void initialize_Plasma(Eigen::Matrix3d ***S, Eigen::Matrix3d ***B);
Eigen::Vector3d transform_geomag(double Inc, double Dec, double th, double ph);

double **** allocate_4d(int dim1, int dim2, int dim3, int dim4, double initial_Value);
double *** allocate_3d(int dim1, int dim2, int dim3, double initial_Value);
double ** allocate_2d(int dim1, int dim2, double initial_Value);
double *allocate_1d(int dim1, double initial_Value);
void free_memory4d(double ****array, int dim1, int dim2, int dim3);
void free_memory3d(double ***array, int dim1, int dim2);
void free_memory2d(double **array, int dim1);
int cal_obs_n0();
void output_E(double ***Er, double ****Eth, double ****Eph, double ***Hr, double ****Hth, double ****Hph, 
                double ****Jr, double ****Jth, double ****Jph, double ****Dr, double ****Dth, double ****Dph, int n, int n0);
void output_pal();
void make_dir();

inline double r(double i){
    return R0 + i * dr;
}
inline double theta(double j){
    // return M_PI / 2.0 - j * dth;
    return M_PI / 2.0 - thR / 2.0 + j * dth;
}

inline double cot(double theta){
    return 1.0 / tan(theta);
}
