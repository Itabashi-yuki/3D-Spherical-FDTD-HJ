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
constexpr double B0 { 50000e-9 };
constexpr double R0 { 6370.0e3 }; /* Radius of the Earth */

/* 解析領域 */
// constexpr double Rr { 50.0e3 };
// constexpr double Rth { 100.0e3 };
// constexpr double Rph { 150.0e3 };


constexpr double Rr { 50.0e3 };
constexpr double Rth { 50.0e3 };
constexpr double Rph { 100.0e3 };
constexpr double dr { 0.5e3 };
constexpr double Rdth { 0.5e3 };
constexpr double Rdph { 0.5e3 };
constexpr double dth { Rdth / R0 };
constexpr double dph { Rdph / R0 };
constexpr double thR { Rth / R0 };
constexpr double phR { Rph / R0 };

constexpr int Nr { int(Rr / dr) };
constexpr int Nth { int(Rth / Rdth) };
constexpr int Nph { int(Rph / Rdph) };

/*xiの導入が必要*/
constexpr double Tmax { 0.001 };
// constexpr double Tmax { 0.1 };
constexpr double Ne_max { 1.0e9 };
constexpr double Omega { std::sqrt( CHARGE_e * CHARGE_e * Ne_max / MASS_e / EPS0 ) };
constexpr double xi { 1.0 / std::sqrt( 1.0 + Omega * Omega / 4.0 / C0 / C0 / ( 1.0 / dr / dr + 1.0 / ( R0 * dth ) / ( R0 * dth ) 
                    + 1.0 / ( R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph ) / ( R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph ) ) ) };
// constexpr double dt { 0.9 / C0 / sqrt( 1.0 / dr / dr + 1.0 / (R0 * dth) / (R0 * dth) 
//                         + 1.0 / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph) / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph)  ) };
constexpr double dt { xi * 0.9999 / C0 / sqrt( 1.0 / dr / dr + 1.0 / (R0 * dth) / (R0 * dth) 
                        + 1.0 / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph) / (R0 * std::sin(M_PI / 2.0 - thR / 2.0) * dph)  ) };
constexpr int Nt { int(Tmax / dt) };

/* PMLパラメタ */
constexpr int PML_L { 12 };
// constexpr int PML_L { 0 };
constexpr double PML_M { 4.0 };
constexpr double PML_R { 1.0e-6 };

/* 電流源パラメタ */
constexpr double f0 { 40.0e3 };
// constexpr double sigma_J { 12 * dt };
constexpr double current_dt { 7.29756e-07 };
constexpr double sigma_J { 12 * current_dt };
constexpr double t0 { 5.0 * sigma_J };
constexpr double source_r { Rr / 2.0 };
constexpr double source_th { Rth / 2.0 };
constexpr double source_ph { Rph / 4.0 };
// constexpr double source_ph { Rph / 2.0  + 10.0e3 };
// constexpr double source_r { 25.0e3 };
// constexpr double source_th { 35.0e3 };
// constexpr double source_ph { 45.0e3 };

/* プラズマ領域パラメタ */
constexpr double Rr_iono_lower { 20.0e3 };
constexpr double Rr_iono_upper { 30.0e3 };
constexpr double Rth_iono_lower { 20.0e3 };
constexpr double Rth_iono_upper { 30.0e3 };
constexpr double Rph_iono_lower { 70.0e3 };
constexpr double Rph_iono_upper { 80.0e3 };
constexpr int Nr_iono_lower { int(Rr_iono_lower / dr) };
constexpr int Nr_iono_upper { int(Rr_iono_upper / dr) };
constexpr int Nth_iono_lower { int(Rth_iono_lower / Rdth)};
constexpr int Nth_iono_upper { int(Rth_iono_upper / Rdth)};
constexpr int Nph_iono_lower { int(Rph_iono_lower / Rdph)};
constexpr int Nph_iono_upper { int(Rph_iono_upper / Rdph)};

constexpr double exp_Ne { 0.0 };
constexpr double exp_nu { 7.0 };

constexpr double THETA { M_PI / 2.0 };
constexpr double PHI { 0.0 };

/*観測パラメタ*/
constexpr double obs_r { source_r };
constexpr double obs_th { source_th };
// constexpr double obs_ph { 140.0e3 };
constexpr double obs_ph { source_ph + 10.0e3 };
constexpr double obs_t_step{ 1.0e-3 };

void update_Er(double ***Er, double ****Hth, double ****Hph, double ****Jr, double ****check, int n);
void update_Eth(double ****Eth, double ***Hr, double ****Hph, double ****Jth, double ****check, int n);
void update_Eph(double ****Eph, double ***Hr, double ****Hth, double ****Jph, double ****check, int n);
void update_Er_PML(double ****Erth1, double ****Erth2, double ****Erph, double ***Er, double ***Hr, double ****Hth, double ****Hph,
                     double *CERTH1_00, double *CERTH1_01, double *CERPH_00, double *CERPH_01, double ****check, int n);
void update_Eth_PML(double ****Ethph, double ****Ethr, double ****Ethr_tilde, double ****Eth, double ***Hr, double ***Hph_tilde, double *CETHPH_00, double *CETHPH_01,
                     double *CETHR_10, double *CETHR_11, double *CETHR_TILDE_00, double *CETHR_TILDE_01, double ****check, int n);
void update_Eph_PML(double ****Ephr, double ****Ephr_tilde, double ****Ephth, double ****Eph, double ***Hr, double ***Hth_tilde,
                     double *CEPHR_10, double *CEPHR_11, double *CEPHR_TILDE_00, double *CEPHR_TILDE_01, double *CEPHTH_00, double *CEPHTH_01, double ****check, int n);
void update_Eth_tilde(double ****Eth_tilde, double ****Eth, double *CETH_TILDE_00, double ****check, int n);
void update_Eph_tilde(double ****Eph_tilde, double ****Eph, double *CEPH_TILDE_00, double ****check, int n);

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

double **** allocate_4d(int dim1, int dim2, int dim3, int dim4, double initial_Value);
double *** allocate_3d(int dim1, int dim2, int dim3, double initial_Value);
double ** allocate_2d(int dim1, int dim2, double initial_Value);
double *allocate_1d(int dim1, double initial_Value);
void free_memory4d(double ****array, int dim1, int dim2, int dim3);
void free_memory3d(double ***array, int dim1, int dim2);
void free_memory2d(double **array, int dim1);
int cal_obs_n0();
void output_E(double ***Er, double ****Eth, double ****Eph, double ***Hr, double ****Hth, double ****Hph, int n, int n0);
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
