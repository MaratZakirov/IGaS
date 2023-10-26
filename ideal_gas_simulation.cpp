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
#define NPARTICLES 10

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

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
};

class scene {
public:
    double step = 0.05;
    double H = 500;
    double W = 100;
    double g = -0.5;
    double cell_size = 1.0;
    double t = 0.0;

    // map from address to recent particle
    std::unordered_map<long, std::set<int>> mask;

    long par2addr(double x, double y) {
        return long(x/cell_size) + long(y/cell_size) * (1 + int(W/cell_size));
    }

    std::vector<particle> particles;
    std::vector<particle> init_particles() {
        for (int i = 0; i < NPARTICLES; i++) {
            // Initialize X, V, m, r for each particle
            particle p;
            p.X[0] = W * dis(gen); // Set X[0] to a value
            p.X[1] = H * dis(gen); // Set X[1] to a value
            p.V[0] = 5 * dis(gen); // Set V[0] to a value
            p.V[1] = 5 * dis(gen); // Set V[1] to a value
            p.m = 1.0;    // Set m to a value
            p.r = 2.0;    // Set r to a value
            particles.push_back(p);
        }
        return particles;
    }

    scene() {
        init_particles();
    }

    void move() {
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            p.X = p.X + step * p.V;
        }
    }

    void assign_masks() {
        //delete previous mask
        mask.clear();
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            for (double x = (p.X[0] - p.r); x <= (p.X[0] + p.r); x += cell_size) {
                for (double y = (p.X[1] - p.r);y <= (p.X[1] + p.r); y += cell_size) {
                    std::set<int>& idx = mask[par2addr(p.X[0], p.X[1])];
                    idx.insert(i);
                }
            }
        }
    }

    // Check collsion is exist exchange velocities
    void collide(particle& p_1, particle& p_2) {
        if (abs(p_1.r - p_2.r) < (p_1.X - p_2.X).norm()) {
            vec_in_space di = p_2.X - p_1.X;
            auto m1 = p_1.m;
            auto m2 = p_2.m;
            auto v1 = p_1.V;
            auto v2 = p_2.V;
            p_1.V = v1 - 2 * m2/(m1 + m2) * dot(v1 - v2, di) / (di.norm()*di.norm()) * di;
            p_2.V = v2 - 2 * m1/(m1 + m2) * dot(v1 - v2, di) / (di.norm()*di.norm()) * di;
        }
    }

    // Use mask to iterate
    void collision() {
        for (auto it = mask.begin(); it != mask.end(); it++) {
            std::set<int>& setpart = it->second;
            for (auto idx1 = setpart.begin(); idx1 != setpart.end(); idx1++) {
                for (auto idx2 = setpart.begin(); idx2 != setpart.end(); idx2++) {
                    // check collision
                    particle& p_i = particles[*idx1];
                    particle& p_j = particles[*idx2];
                    collide(p_i, p_j);
                }
            }
        }
    }

    void reflection() {
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            if ((p.X[0] - p.r < 0) || (p.X[0] + p.r > W))
                p.X[0] *= -1;
            if ((p.X[1] - p.r < 0) || (p.X[1] + p.r > H))
                p.X[1] *= -1;
        }
    }

    void dump_state()
    {
        std::ofstream outputFile;
        if (t == 0) {
            outputFile.open("simulation.csv");
            outputFile << "t\tpnum\tx\ty\tvx\tvy\tr\tm" << std::endl;
        } else {
            outputFile.open("simulation.csv", std::ios::app);
        }
        for (int i = 0; i < particles.size(); i++)
        {
            particle& p = particles[i];
            outputFile << std::fixed << std::setprecision(2) << t << '\t' << i << '\t' << p.X[0] << '\t' << p.X[1] << '\t' << p.V[0] << '\t' << p.V[1] << '\t' << p.r << '\t' << p.m << std::endl;
        }
        outputFile.close();
    }

    void step_scene() {
        reflection();
        assign_masks();
        collision();
        move();
        dump_state();
        t += step;
    }
};

int main() {
    auto sc = scene();
    for (int epoch=0; epoch < 10000; epoch++)
    {
        if (epoch % 100 == 0)
        {
            std::cout << "epoch: " << epoch << std::endl;
        }
        sc.step_scene();
    }
    std::cout << "Done";
    return 0;
}