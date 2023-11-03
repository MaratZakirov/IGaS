#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
#include <set>
#include <list>
#include <fstream>
#include <iomanip>

#define SPACENUM 2
#define NPARTICLES 100
#define CELL 0.01
#define RADIUS 0.01
#define G -0.5
#define WIDTH 200
#define HEIGTH 200
#define STEP 0.05

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> uniform01(0, 1);
std::uniform_real_distribution<> uniform11(-1, 1);

class vec_in_space {
public:
    double item[SPACENUM] = {0.0, 0.0};
    double& operator[](int i) {
        return item[i];
    }
    double norm() {
        double s = 0;
        for (int i = 0; i < SPACENUM; i++)
            s += item[i]*item[i];
    return std::sqrt(s);
}
};

vec_in_space operator+(vec_in_space a, vec_in_space b) {
    vec_in_space c;
    for (int j = 0; j < SPACENUM; j++)
        c[j] = a[j] + b[j];
    return c;
}

vec_in_space operator-(vec_in_space a, vec_in_space b) {
    vec_in_space c;
    for (int j = 0; j < SPACENUM; j++)
        c[j] = a[j] - b[j];
    return c;
}

vec_in_space operator* (vec_in_space a, vec_in_space b) {
    vec_in_space c;
    for (int j = 0; j < SPACENUM; j++)
        c[j] = a[j] * b[j];
    return c;
} 

vec_in_space operator* (double a, vec_in_space b) {
    vec_in_space c;
    for (int j = 0; j < SPACENUM; j++)
        c[j] = a * b[j];
    return c;
}

double dot(vec_in_space a, vec_in_space b)
{
    double c = 0;
    for (int j = 0; j < SPACENUM; j++)
        c += a[j] * b[j];
    return c;
}

class particle {
public:
    vec_in_space X;
    vec_in_space V;
    double m;
    double r;
    bool iscol;
};

class scene {
public:
    double step = STEP;
    double H = HEIGTH;
    double W = WIDTH;
    double g = G;
    double cell_size = CELL;
    double t = 0.0;

    // map from address to recent particle
    std::unordered_map<long, std::set<int>> mask;

    long par2addr(double x, double y) {
        return long(x/cell_size) + long(y/cell_size) * (1 + int(W/cell_size));
    }

    std::vector<particle> particles;
    void init_particles() {
        for (int i = 0; i < NPARTICLES; i++) {
            // Initialize X, V, m, r for each particle
            particle p;
            p.X[0] = W * uniform01(gen); // Set X[0] to a value
            p.X[1] = (H/3) * uniform01(gen); // Set X[1] to a value
            p.V[0] = 6 * uniform11(gen); // Set V[0] to a value
            p.V[1] = 6 * uniform11(gen); // Set V[1] to a value
            p.m = 1.0;    // Set m to a value
            p.r = RADIUS;    // Set r to a value
            bool skip = false;
            for (double x = (p.X[0] - p.r - cell_size);
                 x <= (p.X[0] + p.r + cell_size);
                 x += cell_size) {
                for (double y = (p.X[1] - p.r - cell_size);
                     y <= (p.X[1] + p.r + cell_size);
                     y += cell_size) {
                    std::set<int>& idx = mask[par2addr(x, y)];
                    if (idx.size() == 0) {
                        idx.insert(i);
                    } else {
                        skip = true;
                        if (skip) break;
                    }
                }
                if (skip) break;
            }
            if (skip) continue;
            particles.push_back(p);
        }
        mask.clear();
    }

    scene() {
        init_particles();
    }

    void move() {
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            if (G == 0) {
                p.X = p.X + step * p.V;
            } else {
                vec_in_space g_accel;
                g_accel[0] = 0;
                g_accel[1] = G;
                p.X = p.X + step * p.V + 0.5 * p.m * (step*step) * g_accel;
                p.V = p.V + step * g_accel;
            }
        }
    }

    void assign_masks() {
        //delete previous mask
        mask.clear();
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            p.iscol = false;
            for (double x = (p.X[0] - p.r - cell_size);
                 x <= (p.X[0] + p.r + cell_size);
                 x += cell_size) {
                for (double y = (p.X[1] - p.r - cell_size);
                     y <= (p.X[1] + p.r + cell_size);
                     y += cell_size) {
                    std::set<int>& idx = mask[par2addr(x, y)];
                    idx.insert(i);
                }
            }
        }
    }

    // Check collsion is exist exchange velocities
    void collide(particle& p_1, particle& p_2) {
        vec_in_space di = p_2.X - p_1.X;
        auto m1 = p_1.m;
        auto m2 = p_2.m;
        auto v1 = p_1.V;
        auto v2 = p_2.V;
        vec_in_space d = (1.0/((p_2.X - p_1.X).norm())) * (p_2.X - p_1.X);
        vec_in_space v1_new = v1 - 2 * m2/(m1 + m2) * dot(v1 - v2, d) * d;
        vec_in_space v2_new = v2 - 2 * m1/(m1 + m2) * dot(v2 - v1, d) * d;
        p_1.V = v1_new;
        p_2.V = v2_new;
        std::cout << "Collided!" << std::endl;
    }

    // Use mask to iterate
    int collision(double t) {
        //collision counter
        int col_cnt = 0;
        // firsrt obtain set of pairs
        std::set<std::tuple<int, int>> pairs;
        for (auto it = mask.begin(); it != mask.end(); it++) {
            std::set<int>& setpart = it->second;
            for (auto idx1 = setpart.begin(); idx1 != setpart.end(); idx1++) {
                for (auto idx2 = idx1; idx2 != setpart.end(); idx2++) {
                    // check collision
                    if (*idx1 != *idx2)
                    {
                        particle& p_i = particles[*idx1];
                        particle& p_j = particles[*idx2];

                        auto v1 = p_i.V;
                        auto v2 = p_j.V;
                        auto di = p_j.X - p_i.X;
                        auto r1 = p_i.r;
                        auto r2 = p_j.r;

                        double nor = di.norm();
                        if ((nor - (r1+r2) * 1.1) < (step * abs(dot(v1-v2, di))/nor))
                        //if ((p_i.r + p_j.r) > (p_i.X - p_j.X).norm())
                        {
                            std::tuple<int, int> pair(std::min(*idx1, *idx1), std::max(*idx1, *idx2));
                            pairs.insert(pair);
                        }
                    }
                }
            }
        }
        for (auto it = pairs.begin(); it != pairs.end(); it++)
        {
            auto pair = *it;
            int i = std::get<0>(pair);
            int j = std::get<1>(pair);
            particle& p_i = particles[i];
            particle& p_j = particles[j];
            if (p_i.iscol || p_j.iscol)
                continue;
            collide(p_i, p_j);
            p_i.iscol = true;
            p_j.iscol = true;
            col_cnt += 1;
        }
        return col_cnt;
    }

    void reflection() {
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            vec_in_space X_n = p.X + step * p.V;
            if ((X_n[0] - p.r < 0) || (X_n[0] + p.r > W))
                p.V[0] *= -1;
            if ((X_n[1] - p.r < 0) || (X_n[1] + p.r > H))
                p.V[1] *= -1;
        }
    }

    void dump_state()
    {
        std::ofstream simFile;
        if (t == 0) {
            simFile.open("simulation.csv");
            simFile << "t\tpnum\tx\ty\tvx\tvy\tr\tm\tcol" << std::endl;
        } else {
            simFile.open("simulation.csv", std::ios::app);
        }
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            simFile << std::fixed << std::setprecision(2) << t << '\t' << i << '\t' << p.X[0] << '\t' << p.X[1] << '\t' << p.V[0] << '\t' << p.V[1] << '\t' << p.r << '\t' << p.m << '\t' << p.iscol << std::endl;
        }
        simFile.close();
    }

    void step_scene() {
        reflection();
        assign_masks();
        int colnum = collision(t);
        move();
        if ((colnum > 0) || (int(t / step) % 100 == 0))
            dump_state();
        t += step;
    }
};

int main() {
    std::cout << "Dimensions: " << SPACENUM << std::endl;
    std::cout << "N particles: " << NPARTICLES << std::endl;
    std::cout << "cell size: " << CELL << std::endl;
    std::cout << "Radius of particles: " << RADIUS << std::endl;
    std::cout << "Gravity: " << G << std::endl;
    std::cout << "WIDTH: " << WIDTH << std::endl;
    std::cout << "HEIGTH: " << HEIGTH << std::endl;
    std::cout << "STEP: " << STEP << std::endl;
    auto sc = scene();
    for (int epoch=0; epoch < 300000; epoch++)
    {
        if (epoch % 1000 == 0)
        {
            std::cout << "epoch: " << epoch << std::endl;
        }
        sc.step_scene();
    }
    std::cout << "Done";
    return 0;
}