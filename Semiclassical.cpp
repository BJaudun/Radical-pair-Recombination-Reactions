#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <tuple>
#include <numeric>
#include <fstream>
#include <string>
#include <chrono>


using namespace std;
using namespace std::chrono;

struct Vector3d {
    double x, y, z;

    constexpr Vector3d() : x(0), y(0), z(0) {}
    constexpr Vector3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    constexpr Vector3d(const Vector3d& v) : x{v.x}, y{v.y}, z{v.z} {}
    constexpr Vector3d(Vector3d&& v) : x{v.x}, y{v.y}, z{v.z} {}

    constexpr Vector3d& operator=(const Vector3d& v) 
    {
        x = v.x; y = v.y; z = v.z;
        return *this;
    } 

    constexpr Vector3d& operator=(Vector3d&& v) 
    {
        x = v.x; y = v.y; z = v.z;
        return *this;
    } 

    // Addition of two vectors
    constexpr Vector3d operator+(const Vector3d& v) const {
        return Vector3d(x + v.x, y + v.y, z + v.z);
    }

    constexpr Vector3d& operator+=(const Vector3d& v)  {
        x += v.x;  y += v.y; z += v.z;
        return *this;
    }

    // Subtraction of two vectors
    constexpr Vector3d operator-(const Vector3d& v) const {
        return Vector3d(x - v.x, y - v.y, z - v.z);
    }

    // Scalar multiplication
    constexpr Vector3d operator*(double scalar) const {
        return Vector3d(x * scalar, y * scalar, z * scalar);
    }

    // Scalar division
    constexpr Vector3d operator/(double scalar) const {
        return Vector3d(x / scalar, y / scalar, z / scalar);
    }

    // Dot product
    constexpr double dot(const Vector3d& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    // Cross product
    constexpr Vector3d cross(const Vector3d& v) const {
        return Vector3d(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    // Vector norm (magnitude)
    double norm() const {
        return sqrt(x * x + y * y + z * z);
    }

    // Normalize the vector
    Vector3d normalized() const {
        double n = norm();
        return Vector3d(x / n, y / n, z / n);
    }

    // Multiply two vectors element-wise
    constexpr Vector3d cwiseProduct(const Vector3d& v) const {
        return Vector3d(x * v.x, y * v.y, z * v.z);
    }

};

constexpr double pi = 3.1415927;
constexpr double gyromagnetic_ratio = -1.001;
constexpr double applied_magnetic_field_strength = 0;
constexpr double w_i = -gyromagnetic_ratio * applied_magnetic_field_strength;
constexpr Vector3d omega_bold_i(0, 0, w_i);
vector<double> a_ik = {-0.999985,-0.7369246,0.511210,-0.0826998,0.0655341,-0.562082,-0.905911,0.357729,0.358593,0.869386,-0.232996,0.0388327,0.661931,-0.930856,-0.893077,-0.0594001};

Vector3d random_initial_vector(){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    double theta = 2 * pi * uniform01(generator);
    double phi = acos(1 - 2 * uniform01(generator));
    double sqrt_0_75 = sqrt(0.75);
    return {sqrt_0_75 * sin(phi) * cos(theta), 
            sqrt_0_75 * sin(phi) * sin(theta), 
            sqrt_0_75 * cos(phi)};
}


Vector3d Si_t(const Vector3d& omega_bold, const Vector3d& Si_0, double t) {
    double omega_norm = omega_bold.norm();
    Vector3d omega_normalised{omega_bold / omega_norm};
    double cos_omega_t = cos(omega_norm * t);
    double sin_omega_t = sqrt(1 - cos_omega_t*cos_omega_t);

    Vector3d Si_parallel = omega_normalised * (omega_normalised.dot(Si_0));

    return (Si_parallel) + ((Si_0 - Si_parallel) * cos_omega_t) + ((omega_normalised.cross(Si_0)) * sin_omega_t);
}


Vector3d Iik_t(const Vector3d& omega_bold, const Vector3d& Iik_0, double t) {
    double omega_norm = omega_bold.norm();
    Vector3d omega_normalised{omega_bold / omega_norm};
    double cos_omega_t = cos(omega_norm * t);
    double sin_omega_t = sqrt(1 - cos_omega_t*cos_omega_t); 

    Vector3d Iik_parallel = omega_normalised * (omega_normalised.dot(Iik_0));

    return Iik_parallel + ((Iik_0 - Iik_parallel) * cos_omega_t) + (omega_normalised.cross(Iik_0)* sin_omega_t);
}


Vector3d Si_omega_bold(Vector3d omega_bold_i, const vector<Vector3d>& a_ik_local, const vector<Vector3d>& Iik, int nuclei_number) {
    Vector3d sum_aik_Iik{omega_bold_i};

    for (int i = 0; i < nuclei_number; ++i) {
        sum_aik_Iik += a_ik_local[i].cwiseProduct(Iik[i]);
    }
    return sum_aik_Iik;
}


inline Vector3d Iik_omega_bold(Vector3d a_ik_local, Vector3d Si) {
    return a_ik_local.cwiseProduct(Si);
}

void solve_coupled_eq(double total_time, int time_steps, int nuclei_number, Vector3d Si_0, const vector<Vector3d>& Iik_0, const vector<Vector3d>& a_ik_local, vector<Vector3d>& Si_t_array) {

    vector<Vector3d> Iik = Iik_0;  
    double current_time = 0;
    const double time_increment = total_time / time_steps;

    Si_t_array[0] = Si_0; 
    auto Si_t_current{Si_0}; 
        
    for (int step = 0; step < time_steps - 1; ++step) {
   
        auto S_omega_bold{Si_omega_bold(omega_bold_i, a_ik_local, Iik, nuclei_number)};
        auto Si_mid = Si_t(S_omega_bold, Si_t_current, time_increment / 2.0);

        for (int j = 0; j < nuclei_number; ++j) {
            auto I_omega_bold = Iik_omega_bold(a_ik_local[j], Si_mid);
            Iik[j] = Iik_t(I_omega_bold, Iik[j], time_increment);
        }

        S_omega_bold = Si_omega_bold(omega_bold_i, a_ik_local, Iik, nuclei_number);
        Si_t_current = Si_t(S_omega_bold, Si_mid, time_increment / 2.0);  

        Si_t_array[step + 1] = Si_t_current;
        
        current_time += time_increment;
    }
}


tuple<vector<double>, vector<double>, vector<double>> monte_carlo_integration(double total_time, int time_steps, int nuclei_number, int iterations_Iik, int iterations_Si, const vector<Vector3d>& a_ik_local) {
    vector<double> integral_S_xx(time_steps, 0.0);
    vector<double> integral_S_xy(time_steps, 0.0);
    vector<double> integral_S_zz(time_steps, 0.0);

    vector<double> integral_Iik_xx(time_steps, 0.0);
    vector<double> integral_Iik_xy(time_steps, 0.0);
    vector<double> integral_Iik_zz(time_steps, 0.0);

    vector<Vector3d> Iik_0(nuclei_number);

    vector<Vector3d> Si_t_array(time_steps);

    for (int j = 0; j < iterations_Si; ++j) {
        const Vector3d Si_0  = random_initial_vector();

        fill(integral_Iik_xx.begin(), integral_Iik_xx.end(), 0.0);
        fill(integral_Iik_xy.begin(), integral_Iik_xy.end(), 0.0);
        fill(integral_Iik_zz.begin(), integral_Iik_zz.end(), 0.0);

        for (int i = 0; i < iterations_Iik; ++i) {
            for (int f = 0; f < nuclei_number; ++f) {
                Iik_0[f] = random_initial_vector(); 
            }

            solve_coupled_eq(total_time, time_steps, nuclei_number, Si_0, Iik_0, a_ik_local, Si_t_array);

            for (int l = 0; l < time_steps; ++l) {
                integral_Iik_xx[l] += (Si_0.x * Si_t_array[l].x);
                integral_Iik_xy[l] += (Si_0.x * Si_t_array[l].y);
                integral_Iik_zz[l] += (Si_0.z * Si_t_array[l].z);
            }
        }

        for (int k = 0; k < time_steps; ++k) {
            integral_S_xx[k] += integral_Iik_xx[k] / iterations_Iik;
            integral_S_xy[k] += integral_Iik_xy[k] / iterations_Iik;
            integral_S_zz[k] += integral_Iik_zz[k] / iterations_Iik;
        }
    }

    vector<double> electron_tensor_xx(time_steps);
    vector<double> electron_tensor_xy(time_steps);
    vector<double> electron_tensor_zz(time_steps);

    double inv_iterations_Si = 2.0 / iterations_Si;
    for (int p = 0; p < time_steps; ++p) {
        electron_tensor_xx[p] = integral_S_xx[p] * inv_iterations_Si;
        electron_tensor_xy[p] = integral_S_xy[p] * inv_iterations_Si;
        electron_tensor_zz[p] = integral_S_zz[p] * inv_iterations_Si;
    }

    return make_tuple(electron_tensor_xx, electron_tensor_xy, electron_tensor_zz);
}


int main(int argc, char* argv[]) {
    // Define parameters
    double total_time = 30; 
    int time_steps = 1000;
    int nuclei_number = std::atoi(argv[2]);
    int iterations_Iik = 100000;
    int iterations_Si = 100000;

    // Convert a_ik array to vector of Vector3d
    vector<Vector3d> a_ik_local;
    for (int a = 0; a < nuclei_number; ++a) {
    a_ik_local.emplace_back(a_ik[a], a_ik[a], a_ik[a]);
    }
    
    auto start = high_resolution_clock::now();
    // Run Monte Carlo Integration
    auto [electron_tensor_xx, electron_tensor_xy, electron_tensor_zz] = monte_carlo_integration(total_time, time_steps, nuclei_number, iterations_Iik, iterations_Si, a_ik_local);


    std::ofstream ofs;
    std::string filename = argv[1] + std::to_string(nuclei_number) + "_" + argv[3] + ".csv";

    ofs.open(filename);

    for (int s = 0 ; s < time_steps; ++s)
    {
        ofs << s << ',' << electron_tensor_xx[s] << ',' << electron_tensor_xy[s] << ',' << electron_tensor_zz[s] << '\n';
    }

    ofs.close();


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    std::cout << duration.count() << endl;

    return 0;
}
