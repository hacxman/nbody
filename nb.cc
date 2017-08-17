#include <boost/timer/timer.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;

struct V3d {
  float x;
  float y;
  float z;

  V3d operator - (const V3d &rhs) {
    return V3d({.x = rhs.x - x,
                .y = rhs.y - y,
                .z = rhs.z - z});
  }

  V3d operator / (float rhs) {
    return V3d({.x = x / rhs,
                .y = y / rhs,
                .z = z / rhs});
  }

  V3d & operator += (const V3d &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }

  float norm() {
    return sqrt(x*x + y*y + z*z);
  }
};

V3d operator * (const V3d &rhs, float lhs) {
  return V3d({.x = rhs.x * lhs,
              .y = rhs.y * lhs,
              .z = rhs.z * lhs});
}

V3d operator * (float lhs, const V3d &rhs) {
  return V3d({.x = rhs.x * lhs,
              .y = rhs.y * lhs,
              .z = rhs.z * lhs});
}

std::ostream & operator << (std::ostream &os, const V3d &rhs) {
  os << "[" << rhs.x << " " << rhs.y << " " << rhs.z << "]";
  return os;
}


const float G = 0.12312;

struct Particle {
  V3d position;
  float mass;
  V3d velocity;
  int idx;
  V3d f;
  static int _maxidx;

//  Particle(const Particle &p) = default;

  void setidx() {
    idx = Particle::_maxidx++;
  }

  bool operator == (const Particle &rhs) {
    return idx == rhs.idx;
  }

  V3d force(const struct Particle &p2) {
    return G * mass * p2.mass * ( position - p2.position)/pow(( position - p2.position).norm(), 3);
  }
};

int main(void) {
  boost::timer::auto_cpu_timer _timer;
  vector<Particle> particles;
  const float dt = 0.01;

  const int CNT = 4000;

  for (int i = 0; i < CNT; i++) {
    Particle p = {.position = {i, i, i},
                  .mass = 1,
                  .velocity = {1, 0, 0},
                  .idx = i,
                  .f = {0, 0, 0}};
    particles.push_back(p);
  }

  int _i = 0;
  const int _itop = 0x10;
  cpu_timer timer;
  long long int C = 0;
  for (int t = 0; t < _itop * 10; t++, _i++) {
    int p_cnt = particles.size();
#pragma omp parallel for
#pragma acc kernels
    for (int i = 0; i < p_cnt; i++) {
      auto &p = particles[i];
//    for (auto &p : particles) {
      for (int j = 0; j < p_cnt; j++) {
//      for (auto p2 : particles) {
        //if (p == p2) continue;
        auto p2 = particles[j];
        if (i == j) continue;

        p.f += p.force(p2);
      }
    }
    if (_i > _itop) {
      cout << fixed << setprecision( 2 ) << float(C)/ 1e06 << "MC" << endl;
    }
    C+=particles.size() * (particles.size() - 1);
#pragma omp parallel for
#pragma acc kernels
    for (int i = 0; i < p_cnt; i++) {
      auto &p = particles[i];
    //for (auto p : particles) {
      //cout << p.f << "; ";
      p.velocity += (p.f / p.mass) * dt;
      p.position += p.velocity * dt;
      if (_i > _itop) {
      //  cout << p.position << " "; //<< p.f << "; ";
      //  cout << fixed << setprecision( 2 ) << float(C)/ 1e06 << "MC/s" << endl;
      }
      //cout << p.f << "; ";

      //p.f = {0, 0, 0};
    }
    if (_i > _itop) {
      cout << endl;
      _i = 0;
    }
  }
  C+=particles.size();
  cpu_times const elapsed_times(timer.elapsed());
  nanosecond_type const elapsed(elapsed_times.system
    + elapsed_times.user);
  cout << fixed << setprecision( 2 ) << float(C*1000)/float(elapsed_times.wall) << "MC/s" << " " << elapsed << " " << C << endl;
}

