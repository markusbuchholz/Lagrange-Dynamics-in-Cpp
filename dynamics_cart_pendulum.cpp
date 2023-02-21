#include <iostream>
#include <tuple>
#include <vector>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//----------- system dynamic parameters --------------------

float M = 2.0;
float m = 1.0;
float R = 2.5;
float b = 0.005; // cart friction
float d = 0.001; // pendulum friction
float k = 1.0;
float g = 9.81;
float dt = 0.001;

//-----------------------------------------------------------
// dx
float function1(float x, float theta, float dx, float dtheta)
{

    return dx;
}

//-----------------------------------------------------------
// dtheta
float function2(float x, float theta, float dx, float dtheta)
{

    return dtheta;
}

//-----------------------------------------------------------
// dv
float function3(float x, float theta, float dx, float dtheta)
{

    auto fa = m * R * dtheta * dtheta * std::sin(theta) + m * g * std::sin(theta) * std::cos(theta) - k * x - d * dx + (float)b / R * dtheta * std::cos(theta);
    auto fb = M + m * std::sin(theta) * std::sin(theta);

    return (float)fa / fb;
}

//-----------------------------------------------------------
// domega
float function4(float x, float theta, float dx, float dtheta)
{

    auto fa = -m * R * dtheta * dtheta * std::sin(theta) * std::cos(theta) - (m + M) * g * std::sin(theta) + k * x * std::cos(theta) + d * dx * std::cos(theta) - (1 + (float)M / m) * (float)b / R * dtheta;
    auto fb = R * (M + m * std::sin(theta) * std::sin(theta));

    return (float)fa / fb;
}

//-----------------------------------------------------------

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> methodRuneKuttaDynamics()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;
    std::vector<float> diffEq3;
    std::vector<float> diffEq4;

    std::vector<float> time;

    // init values
    float x1 = 2.0;
    float x2 = (M_PI / 180) * 30;
    float x3 = 0;
    float x4 = 0;
    float t = 0.0; // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    diffEq3.push_back(x3);
    diffEq4.push_back(x4);
    time.push_back(t);

    for (int ii = 0; ii < 30000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2, x3, x4);
        float k12 = function2(x1, x2, x3, x4);
        float k13 = function3(x1, x2, x3, x4);
        float k14 = function4(x1, x2, x3, x4);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k23 = function3(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);
        float k24 = function4(x1 + dt / 2 * k11, x2 + dt / 2 * k12, x3 + dt / 2 * k13, x4 + dt / 2 * k14);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k33 = function3(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);
        float k34 = function4(x1 + dt / 2 * k21, x2 + dt / 2 * k22, x3 + dt / 2 * k23, x4 + dt / 2 * k24);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k43 = function3(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);
        float k44 = function4(x1 + dt * k31, x2 + dt * k32, x3 + dt * k33, x4 + dt * k34);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);
        x3 = x3 + dt / 6.0 * (k13 + 2 * k23 + 2 * k33 + k43);
        x4 = x4 + dt / 6.0 * (k14 + 2 * k24 + 2 * k34 + k44);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        diffEq3.push_back(x3);
        diffEq4.push_back(x4);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, diffEq3, diffEq4, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::vector<float> X, std::vector<float> Y)
{

    plt::title("A Pendulum Attached to an Oscillating Mass ");
    plt::named_plot("positions", X, Y);
    plt::xlabel("mass pos");
    plt::ylabel("angle");
    plt::legend();
    plt::xlabel("mass pos");
    plt::ylabel("angle");
    plt::show();
}

//---------------------------------------------------------------
int main()
{

    auto dyn = methodRuneKuttaDynamics();

    auto pos = std::get<0>(dyn);
    auto ang = std::get<1>(dyn);

    plot2D(pos, ang);
}
