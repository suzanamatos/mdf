#include "timeSpent.h"

#include <iostream>

TimeSpent::TimeSpent(std::string _name)
    : name(_name), accumulator(0)//, start(0), end(0)
{
}

TimeSpent::~TimeSpent()
{
}

void TimeSpent::startCount()
{
    start = std::chrono::system_clock::now();
}

void TimeSpent::endCount()
{
    end = std::chrono::system_clock::now();
    accumulator += end - start;
}

void TimeSpent::print()
{
    std::cout << name << " \t=> levou um tempo total de "
              << accumulator.count() << "s" << std::endl << std::endl;
}

void TimeSpent::setZero()
{
    accumulator = end-end;
}

std::chrono::duration<double> TimeSpent::getAccumulator() const
{
    return accumulator;
}
