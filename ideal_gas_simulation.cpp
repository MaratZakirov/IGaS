#include <iostream>
#include <vector>
#include <cmath>
#include <random>

#define SPACENUM 2
#define NPARTICLES 100

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

class vec_in_space {
public:
    double item[SPACENUM] = {0.0, 0.0};

    double norm(vec_in_space a, vec_in_space b) {
        double s = 0;
        for (int i = 0; i < SPACENUM; i++)
            s += item[i];
        return std::sqrt(s);
    }

    double& operator[](int i) {
        return item[i];
    }
};

vec_in_space operator+(vec_in_space a, vec_in_space b) {
    vec_in_space c;
    for (int j = 0; j < SPACENUM; j++)
        c[j] = a[j] + b[j];
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

class particle {
public:
    vec_in_space X;
    vec_in_space V;
    double m;
    double r;
};

class scene {
public:
    double step=0.05;
    double H = 500;
    double W = 100;
    double g = -0.5;
    std::hash_map<long, std::list<int>> mask;

    std::tuple<int, int> particle2maskaddr(particle& p) {
        long addr;
        for (int j = 0; j < SPACENUM; j++)
            addr += int(p.X[j]) * j * 1000;
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
        
    }

    void collision() {
        ;
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

    void step_scene() {
        move();
        reflection();
        collision();
    }
};

int main() {
    std::vector<int> a;
    std::cout << "Hellgggggo World!";
    return 0;
}